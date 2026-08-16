"""Preparation of a GROMACS QM/MM topology with link atoms.

The module takes a finished MM system (a Gromacs ``.top`` plus a coordinate
file), carves a QM region out of it and writes everything the patched
GROMACS/DFTB+ build needs:

* ``*.top``  -- the QM/MM topology: QM-QM bonds turned into connections, link
  atoms and charge points added as ``[ virtual_sites2 ]``, charges adjusted;
* ``*.gro``  -- coordinates in exactly the same atom order as the topology;
* ``*.ndx``  -- the ``[ QM ]`` group (QM atoms + link atoms) and two helper groups;
* ``dftb_in.hsd`` -- the DFTB+/xTB input, with the QM atoms in the same order.

Typical use::

    import QMMMtools

    qm = QMMMtools.QM('npt.gro', 'topol.top', 'qm.gro', 'qm.top', 'qm.ndx')
    qm.choose_qm_to_extend('@21313,21391,21408')          # grow to breakable bonds
    qm.choose_qm_manually(f'(({qm.qm_input_mask})<:3.0)&(:SOL)')  # + nearby water
    qm.job()                                              # charge is detected automatically
    qm.make_hsd('dftb_in.hsd', method='dftb3-d4', skpath='./3ob-3-1-ophyd/')
    qm.check_consistency()

The chemistry-specific tables (which bonds may be cut, link-atom distances,
residue names, QM methods) live in :mod:`QMMMtools.data`; switching from a protein to
a nucleic acid or a lipid system means replacing ``breakable_bonds`` and
``H_dist``, nothing else.
"""

import logging
import math
import re
from pathlib import Path

import numpy as np
import parmed as pmd
from parmed import Atom
from parmed.geometry import box_lengths_and_angles_to_vectors, reduce_box_vectors

from . import data
from .data import (AMINO_ACIDS, ELEMENTS, HUBBARD_DERIVS, H_DIST_BY_ELEMENT,  # noqa: F401
                       H_DIST_FALLBACK, LIPID_BREAKABLE_BONDS, LIPID_H_DIST,
                       MAX_ANGULAR_MOMENTUM, NUCLEIC_ACIDS, NUCLEIC_BREAKABLE_BONDS,
                       NUCLEIC_H_DIST, POLYMER_RESIDUES, PROTEIN_BREAKABLE_BONDS,
                       PROTEIN_H_DIST, QM_METHODS, DEFAULT_QM_METHOD, SOLVENT_AND_IONS,
                       TYPE2ELEMENT, WATER_RESIDUES, ION_RESIDUES)

__all__ = ['QM', 'QMMMError', 'QMGeometry', 'HsdFile', 'set_log_level',
           'read_qm_geometry', 'read_index_file', 'write_hsd', 'rewrite_hsd',
           'hamiltonian_block', 'list_methods', 'get_method']

LOGGER = logging.getLogger('QMMMtools')


class QMMMError(RuntimeError):
    """Raised for every problem this module detects in the input or the setup."""


def _init_logging():
    """Print progress to stdout unless the host application configured logging."""
    if not LOGGER.handlers:
        handler = logging.StreamHandler()
        handler.setFormatter(logging.Formatter('[QMMMtools] %(message)s'))
        LOGGER.addHandler(handler)
        LOGGER.setLevel(logging.INFO)
        LOGGER.propagate = False


def set_log_level(level):
    """``set_log_level('DEBUG')`` for a detailed trace, ``'WARNING'`` to keep quiet."""
    _init_logging()
    LOGGER.setLevel(level)


_init_logging()


# directives that belong to a [ moleculetype ] block; used to find where the
# first moleculetype of the input topology ends
_MOLECULETYPE_SECTIONS = frozenset({
    'atoms', 'bonds', 'pairs', 'pairs_nb', 'angles', 'dihedrals', 'exclusions',
    'constraints', 'settles', 'cmap', 'position_restraints', 'distance_restraints',
    'dihedral_restraints', 'orientation_restraints', 'angle_restraints',
    'angle_restraints_z', 'virtual_sites1', 'virtual_sites2', 'virtual_sites3',
    'virtual_sites4', 'virtual_sitesn', 'polarization', 'thole_polarization',
    'water_polarization',
})

_DIRECTIVE_RE = re.compile(r'^\s*\[\s*([A-Za-z_0-9]+)\s*\]')
_INCLUDE_RE = re.compile(r'^(\s*#include\s+")([^"]+)("\s*.*)$')


def _directive(line):
    """Return the name inside ``[ ... ]`` or None."""
    match = _DIRECTIVE_RE.match(line)
    return match.group(1).lower() if match else None


# How a fractional force-field charge becomes the integer the QM code needs.
# Ordinary rounding switches over at a fractional part of 0.5; here the threshold
# is a parameter, because a QM/MM cut that runs through a charge group leaves part
# of a formal charge behind on the MM side and the sum comes out too small in
# magnitude.  A lower threshold therefore rounds up in magnitude more readily.
DEFAULT_CHARGE_ROUNDING = 0.25

#: convenience names for the two extremes
CHARGE_ROUNDING_ALIASES = {'nearest': 0.5, 'away': 0.0, 'default': DEFAULT_CHARGE_ROUNDING}


def _rounding_threshold(rounding):
    """Accept a number or one of :data:`CHARGE_ROUNDING_ALIASES`."""
    if rounding is None:
        return DEFAULT_CHARGE_ROUNDING
    if isinstance(rounding, str):
        key = rounding.strip().lower()
        if key in CHARGE_ROUNDING_ALIASES:
            return CHARGE_ROUNDING_ALIASES[key]
        try:
            rounding = float(key)
        except ValueError:
            raise QMMMError(f'unknown charge rounding {rounding!r}: give a threshold between '
                            '0 and 1, or one of '
                            + ', '.join(sorted(CHARGE_ROUNDING_ALIASES))) from None
    threshold = float(rounding)
    if not 0.0 <= threshold <= 1.0:
        raise QMMMError(f'the charge rounding threshold must lie between 0 and 1, got {threshold}')
    return threshold


def _round_charge(value, rounding=None):
    """Round ``value`` to an integer, switching over at a fractional part of ``rounding``.

    ``rounding=0.5`` is ordinary rounding; the default ``0.25`` rounds up in
    magnitude as soon as more than a quarter of an electron is left over, so
    ``1.43 -> 2`` and ``-1.43 -> -2`` while ``0.04 -> 0`` and ``2.0 -> 2``.
    """
    threshold = _rounding_threshold(rounding)
    magnitude = abs(value)
    whole = math.floor(magnitude + 1e-9)          # the epsilon absorbs 1.9999999999
    if magnitude - whole > threshold + 1e-9:
        whole += 1
    return int(whole) if value >= 0 else -int(whole)


# ===========================================================================
# DFTB+ / xTB input files
# ===========================================================================
class QMGeometry:
    """The QM atoms in the form the ``.hsd`` needs them.

    Produced either by :meth:`QM.qm_geometry` (topology available, elements taken
    from the atomic numbers) or by :func:`read_qm_geometry` (only a coordinate
    file and an index file, elements guessed from the atom names or inherited
    from the ``.hsd`` that is being rewritten).

    Parameters
    ----------
    coordinates : array-like, (natoms, 3)
        Cartesian coordinates in **Angstrom**.
    elements : list of str or None
        chemical element per atom.  ``None`` means "unknown": such a geometry can
        still be used to refresh the coordinates of an existing ``.hsd``, whose
        ``TypeNames`` are then kept.
    labels : list of str, optional
        ``RES123:NAME`` per atom, only used in messages.
    """

    def __init__(self, coordinates, elements=None, labels=None):
        self.coordinates = np.asarray(coordinates, dtype=float).reshape(-1, 3)
        self.elements = list(elements) if elements is not None else None
        self.labels = list(labels) if labels is not None else None
        if self.elements is not None and len(self.elements) != len(self.coordinates):
            raise QMMMError('elements and coordinates have different lengths')

    def __len__(self):
        return len(self.coordinates)

    @property
    def type_names(self):
        """Distinct elements in first-appearance order.

        Never a ``set``: the order defines the species indices of the geometry
        block and has to be reproducible from run to run.
        """
        if self.elements is None:
            raise QMMMError('the elements of this geometry are unknown')
        names = []
        for element in self.elements:
            if element not in names:
                names.append(element)
        return names

    def block(self, type_names=None, species=None):
        """The ``Geometry = { ... }`` text.

        ``type_names`` and ``species`` allow reusing the assignment of an existing
        file (see :func:`rewrite_hsd`); otherwise both come from :attr:`elements`.
        """
        if type_names is None or species is None:
            type_names = self.type_names
            species = [type_names.index(e) + 1 for e in self.elements]
        names = ' '.join(f'"{name}"' for name in type_names)
        out = [f'Geometry = {{\n  TypeNames = {{ {names} }}\n'
               f'  TypesAndCoordinates [Angstrom] = {{\n']
        for kind, position in zip(species, self.coordinates):
            out.append(f'{int(kind):4d} {position[0]:12.6f} {position[1]:12.6f} '
                       f'{position[2]:12.6f}\n')
        out.append('  }\n}\n')
        return ''.join(out)


class HsdFile:
    """A minimal editor for DFTB+ HSD input.

    HSD is a tree of ``Name = Value`` and ``Name = Tag { ... }`` nodes.  This
    class only needs to find a named block, replace it, and patch scalar
    assignments, which is enough to change the geometry, the Hamiltonian, the
    charge or the Slater-Koster path of an existing input without touching
    anything else in the file.
    """

    def __init__(self, text=''):
        self.text = text

    @classmethod
    def read(cls, path):
        path = Path(path).expanduser()
        if not path.is_file():
            raise QMMMError(f'cannot read {path}: file not found')
        return cls(path.read_text())

    def write(self, path):
        path = Path(path).expanduser()
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(self.text)
        return path

    # -- blocks ------------------------------------------------------------
    def block_span(self, name, span=None):
        """``(start, stop, body_start, body_stop)`` of ``name = ... { ... }``.

        Only blocks at the top level of ``span`` (the whole file by default) are
        considered, so ``block_span('Charge')`` never matches a ``Charge`` that
        sits inside some other block.
        """
        start, stop = span if span else (0, len(self.text))
        pattern = re.compile(r'(?m)^[ \t]*' + re.escape(name) + r'\b[^\n{]*\{')
        position, depth = start, 0
        while position < stop:
            character = self.text[position]
            if character == '{':
                depth += 1
            elif character == '}':
                depth -= 1
            elif depth == 0:
                match = pattern.match(self.text, position)
                if match and (position == start or self.text[position - 1] == '\n'):
                    body_start = match.end()
                    end = self._matching_brace(match.end() - 1)
                    return match.start(), end + 1, body_start, end
            position += 1
        return None

    def _matching_brace(self, opening):
        depth = 0
        for position in range(opening, len(self.text)):
            if self.text[position] == '{':
                depth += 1
            elif self.text[position] == '}':
                depth -= 1
                if depth == 0:
                    return position
        raise QMMMError('unbalanced braces in the HSD input')

    def get_block(self, name, span=None):
        found = self.block_span(name, span)
        return None if found is None else self.text[found[0]:found[1]]

    def set_block(self, name, text, span=None):
        """Replace a block, append it when it is not there, drop it on ``None``."""
        found = self.block_span(name, span)
        if text is None:
            if found:
                self.text = self.text[:found[0]] + self.text[found[1]:].lstrip('\n')
                LOGGER.debug('removed the %s block', name)
            return
        text = text if text.endswith('\n') else text + '\n'
        if found:
            self.text = self.text[:found[0]] + text + self.text[found[1]:]
        else:
            self.text = self.text.rstrip('\n') + '\n' + text
            LOGGER.debug('appended a new %s block', name)

    # -- scalars -----------------------------------------------------------
    def get_value(self, key, span=None):
        start, stop = span if span else (0, len(self.text))
        match = re.search(r'(?m)^[ \t]*' + re.escape(key) + r'[ \t]*=[ \t]*([^\n{]+)$',
                          self.text[start:stop])
        return match.group(1).strip() if match else None

    def set_value(self, key, value, span=None, indent='  '):
        """Patch ``key = value``; insert it at the top of ``span`` when missing."""
        start, stop = span if span else (0, len(self.text))
        pattern = re.compile(r'(?m)^([ \t]*)' + re.escape(key) + r'([ \t]*=[ \t]*)([^\n{]*)$')
        match = pattern.search(self.text, start, stop)
        if match:
            self.text = (self.text[:match.start()]
                         + f'{match.group(1)}{key}{match.group(2)}{value}'
                         + self.text[match.end():])
            return True
        self.text = self.text[:start] + f'\n{indent}{key} = {value}' + self.text[start:]
        return False


def _scc_settings(scc_tolerance, max_scc_iterations, mixer, indent='  '):
    """The SCC convergence lines shared by the DFTB and the xTB Hamiltonian."""
    out = []
    if scc_tolerance is not None:
        out.append(f'{indent}SCCTolerance = {scc_tolerance}\n')
    if max_scc_iterations:
        out.append(f'{indent}MaxSCCIterations = {int(max_scc_iterations)}\n')
    if mixer:
        out.append(indent + _mixer_block(mixer).strip() + '\n')
    return ''.join(out)


def _mixer_block(mixer):
    """``mixer='anderson'|'broyden'`` shorthands, or a raw HSD string."""
    if isinstance(mixer, str) and mixer.lower() == 'anderson':
        # verified to converge a Mg(2+)/pyrophosphate QM region that the DFTB+
        # default Broyden mixer could not converge in 500 iterations
        return 'Mixer = Anderson {\n    MixingParameter = 0.05\n    Generations = 8\n  }'
    if isinstance(mixer, str) and mixer.lower() == 'broyden':
        return 'Mixer = Broyden {\n    MixingParameter = 0.05\n  }'
    return str(mixer)


def get_method(method):
    """Look a QM method up by name, with a helpful error for typos."""
    if isinstance(method, data.QMMethod):
        return method
    key = str(method).lower()
    if key not in QM_METHODS:
        raise QMMMError(f'unknown QM method {method!r}; available: ' + ', '.join(sorted(QM_METHODS)))
    return QM_METHODS[key]


def list_methods():
    """Print the available QM methods."""
    for name in sorted(QM_METHODS):
        entry = QM_METHODS[name]
        kind = 'needs Slater-Koster files' if entry.needs_slater_koster else 'self-contained'
        print(f'  {name:<12s} {entry.description} ({kind})')


def hamiltonian_block(method, elements, charge, skpath=None, sk_suffix='-c.spl',
                      sk_separator='', sk_lowercase=True, scc_tolerance='1e-6',
                      max_scc_iterations=250, mixer=None):
    """Build the ``Hamiltonian`` block of a ``dftb_in.hsd``."""
    method = get_method(method)
    charge = round(charge)
    scc = _scc_settings(scc_tolerance, max_scc_iterations, mixer)
    if method.kind == 'xtb':
        return (f'Hamiltonian = xTB {{\n'
                f'  Method = "{method.xtb_method}"\n'
                f'{scc}'
                f'  Charge = {charge}\n'
                f'}}\n')

    missing = [e for e in elements if e not in MAX_ANGULAR_MOMENTUM]
    if missing:
        raise QMMMError(f'no DFTB angular momentum known for {missing}; '
                        'add them to QMMMtools.data.MAX_ANGULAR_MOMENTUM')
    if skpath is None:
        raise QMMMError(f'method {method.name!r} needs Slater-Koster files: '
                        'pass skpath="/path/to/3ob-3-1/"')
    skpath = str(skpath)
    if not skpath.endswith('/'):
        skpath += '/'
    if not Path(skpath).expanduser().is_dir():
        LOGGER.warning('the Slater-Koster directory %s does not exist here; make sure it '
                       'does on the machine that runs mdrun', skpath)

    text = ['Hamiltonian = DFTB {\n', '  SCC = Yes\n', scc,
            f'  Charge = {charge}\n', '  MaxAngularMomentum {\n']
    text += [f'    {e} = "{MAX_ANGULAR_MOMENTUM[e]}"\n' for e in elements]
    text.append('  }\n')
    text += ['  SlaterKosterFiles = Type2FileNames {\n',
             f'    Prefix = {skpath}\n',
             f'    Separator = "{sk_separator}"\n',
             f'    LowerCaseTypeName = {"Yes" if sk_lowercase else "No"}\n',
             f'    Suffix = "{sk_suffix}"\n  }}\n']
    if method.third_order:
        unknown = [e for e in elements if e not in HUBBARD_DERIVS]
        if unknown:
            raise QMMMError(f'no Hubbard derivative known for {unknown}; '
                            'add them to QMMMtools.data.HUBBARD_DERIVS')
        text.append('  ThirdOrderFull = Yes\n  HubbardDerivs {\n')
        text += [f'    {e} = {HUBBARD_DERIVS[e]}\n' for e in elements]
        text.append('  }\n')
    for block in (method.hcorrection, method.dispersion):
        if block:
            text.append('  ' + block + '\n')
    text.append('}\n')
    return ''.join(text)


ANALYSIS_BLOCK = ('Analysis = {\n'
                  '  CalculateForces = Yes\n'
                  '  ProjectStates = {}\n'
                  '  WriteEigenvectors = No\n'
                  '  WriteBandOut = No\n'
                  '  MullikenAnalysis = No\n'
                  '  AtomResolvedEnergies = No\n'
                  '}\n')

OPTIONS_BLOCK = ('Options = {\n'
                 '  WriteDetailedOut = No\n'
                 '  WriteAutotestTag = No\n'
                 '  WriteDetailedXML = No\n'
                 '  WriteResultsTag = No\n'
                 '  RestartFrequency = 2000\n'
                 '  RandomSeed = 0\n'
                 '  WriteHS = No\n'
                 '  WriteRealHS = No\n'
                 '  MinimiseMemoryUsage = No\n'
                 '  ShowFoldedCoords = No\n'
                 '  TimingVerbosity = 0\n'
                 '  WriteChargesAsText = No\n'
                 '}\n')


def write_hsd(file_hsd, geometry, charge, method=DEFAULT_QM_METHOD, skpath=None,
              analysis=ANALYSIS_BLOCK, options=OPTIONS_BLOCK, extra='', **hamiltonian_kwargs):
    """Write a complete ``dftb_in.hsd`` from a :class:`QMGeometry`."""
    if geometry.elements is None:
        raise QMMMError('the elements of the QM atoms are unknown; pass elements= to '
                        'read_qm_geometry(), or use rewrite_hsd() which can inherit the '
                        'TypeNames of the file it rewrites')
    elements = geometry.type_names
    text = geometry.block()
    text += hamiltonian_block(method, elements, charge, skpath=skpath, **hamiltonian_kwargs)
    text += (analysis or '') + (options or '')
    if extra:
        text += extra if extra.endswith('\n') else extra + '\n'
    path = HsdFile(text).write(file_hsd)
    LOGGER.info('wrote %s: %s, %d atoms, elements %s, charge %+d', path,
                get_method(method).name, len(geometry), ' '.join(elements), round(charge))
    return path


def rewrite_hsd(file_hsd, geometry=None, source_hsd=None, keep_types=True, method=None,
                charge=None, skpath=None, scc_tolerance=None, max_scc_iterations=None,
                mixer=None, analysis=None, options=None, blocks=None, **hamiltonian_kwargs):
    """Update parts of an existing ``dftb_in.hsd`` in place.

    Everything that is not addressed stays byte for byte as it was, so
    hand-tuned settings survive.

    Parameters
    ----------
    file_hsd : path-like
        file to write; also the file that is read when ``source_hsd`` is None.
    geometry : QMGeometry, optional
        new coordinates.  With ``keep_types`` and an unchanged atom count the
        ``TypeNames`` and the species column of the old file are reused, which is
        what makes a rewrite from a bare ``.gro`` + ``.ndx`` possible.
    method : str, optional
        replace the whole ``Hamiltonian`` block with a freshly built one.  Needs
        the elements (from ``geometry`` or the old ``TypeNames``) and a charge
        (from ``charge`` or the old file).
    charge, skpath, scc_tolerance, max_scc_iterations, mixer
        patched into the existing ``Hamiltonian`` when ``method`` is None.
    analysis, options : str, optional
        replace the ``Analysis`` / ``Options`` blocks; ``''`` removes them.
    blocks : dict, optional
        any other top-level blocks, ``{'Driver': 'Driver = {}'}``; a value of
        ``None`` deletes the block.
    """
    hsd = HsdFile.read(source_hsd if source_hsd else file_hsd)
    old_types, old_species = _read_geometry_types(hsd)
    changed = []

    if geometry is not None:
        types, species = None, None
        if keep_types and old_species is not None and len(old_species) == len(geometry):
            old_elements = [old_types[s - 1] for s in old_species] if old_types else None
            if geometry.elements is None or old_elements == list(geometry.elements):
                types, species = old_types, old_species
            else:
                differing = [i for i, (o, n) in enumerate(zip(old_elements, geometry.elements))
                             if o != n]
                examples = ', '.join(
                    f'#{i + 1}{" " + geometry.labels[i] if geometry.labels else ""} '
                    f'{old_elements[i]}->{geometry.elements[i]}' for i in differing[:3])
                LOGGER.warning('%d of %d atoms have a different element than the file being '
                               'rewritten (%s); rebuilding TypeNames from the new geometry. '
                               'Check that the .hsd and the QM group really describe the same '
                               'selection.', len(differing), len(geometry), examples)
        elif old_species is not None and len(old_species) != len(geometry):
            LOGGER.warning('the old geometry had %d atoms, the new one has %d -- regenerate '
                           'the index file and the .tpr as well', len(old_species), len(geometry))
        hsd.set_block('Geometry', geometry.block(types, species))
        changed.append(f'geometry ({len(geometry)} atoms)')

    hamiltonian = hsd.block_span('Hamiltonian')
    if method is not None:
        elements = (geometry.type_names if geometry is not None and geometry.elements is not None
                    else old_types)
        if not elements:
            raise QMMMError('cannot rebuild the Hamiltonian: the elements are unknown')
        new_charge = charge
        if new_charge is None:
            old_charge = hsd.get_value('Charge', hamiltonian[2:4] if hamiltonian else None)
            if old_charge is None:
                raise QMMMError('cannot rebuild the Hamiltonian: no charge given and none '
                                'found in the old file')
            new_charge = float(old_charge)
        hsd.set_block('Hamiltonian',
                      hamiltonian_block(method, elements, new_charge, skpath=skpath,
                                        scc_tolerance=scc_tolerance or '1e-6',
                                        max_scc_iterations=max_scc_iterations or 250,
                                        mixer=mixer, **hamiltonian_kwargs))
        changed.append(f'Hamiltonian -> {get_method(method).name}')
    else:
        if hamiltonian is None and any(v is not None for v in
                                       (charge, skpath, scc_tolerance, max_scc_iterations, mixer)):
            raise QMMMError(f'{file_hsd}: no Hamiltonian block to patch')
        body = hamiltonian[2:4] if hamiltonian else None
        if charge is not None:
            hsd.set_value('Charge', round(charge), body)
            changed.append(f'charge {round(charge):+d}')
        if scc_tolerance is not None:
            hsd.set_value('SCCTolerance', scc_tolerance, body)
            changed.append('SCCTolerance')
        if max_scc_iterations is not None:
            hsd.set_value('MaxSCCIterations', int(max_scc_iterations), body)
            changed.append('MaxSCCIterations')
        if mixer is not None:
            hsd.set_block('Mixer', '  ' + _mixer_block(mixer) if mixer else None,
                          hsd.block_span('Hamiltonian')[2:4])
            changed.append('Mixer')
        if skpath is not None:
            path = str(skpath) if str(skpath).endswith('/') else str(skpath) + '/'
            sk = hsd.block_span('SlaterKosterFiles', hsd.block_span('Hamiltonian')[2:4])
            if sk is None:
                raise QMMMError(f'{file_hsd}: no SlaterKosterFiles block to point elsewhere')
            hsd.set_value('Prefix', path, sk[2:4], indent='    ')
            changed.append(f'skpath {path}')

    for name, text in (('Analysis', analysis), ('Options', options)):
        if text is not None:
            hsd.set_block(name, text or None)
            changed.append(name)
    for name, text in (blocks or {}).items():
        hsd.set_block(name, text)
        changed.append(name)

    path = hsd.write(file_hsd)
    LOGGER.info('rewrote %s: %s', path, ', '.join(changed) if changed else 'nothing to change')
    return path


def _read_geometry_types(hsd):
    """``TypeNames`` and the species column of the ``Geometry`` block of a file."""
    block = hsd.get_block('Geometry')
    if block is None:
        return None, None
    # both the braced form written here and the bare "TypeNames = S O N C H"
    # that older hand-made inputs use
    names = re.search(r'TypeNames\s*=\s*(?:\{([^}]*)\}|([^\n{]+))', block)
    types = ([t.strip('" \t,') for t in (names.group(1) or names.group(2)).split()]
             if names else None)
    species = [int(line.split()[0]) for line in block.splitlines()
               if re.match(r'^\s*\d+\s+-?[\d.]', line)]
    return types, (species or None)


def read_qm_geometry(file_gro, file_ndx, group='QM', elements=None, warn_unknown=True):
    """Read the QM atoms straight from a coordinate file and an index file.

    This is the "I only have ``qm.gro`` and ``qm.ndx``" entry point: no topology
    is needed, so it also works for files somebody else produced.

    Because a ``.gro`` carries no atomic numbers the elements are guessed from
    the atom and residue names (see :func:`QMMMtools.data.guess_element`).  Pass
    ``elements`` to override, or leave the guessing to :func:`rewrite_hsd`, which
    prefers the ``TypeNames`` already present in the file it rewrites.

    Any format ParmEd can read works, not just ``.gro``.
    """
    file_gro, file_ndx = Path(file_gro).expanduser(), Path(file_ndx).expanduser()
    for path in (file_gro, file_ndx):
        if not path.is_file():
            raise QMMMError(f'input file not found: {path}')
    groups = read_index_file(file_ndx)
    if group not in groups:
        raise QMMMError(f'{file_ndx} has no [ {group} ] group; it has: '
                        + ', '.join(groups) if groups else f'{file_ndx} defines no group')
    indexes = groups[group]
    if not indexes:
        raise QMMMError(f'the [ {group} ] group of {file_ndx} is empty')

    frame = pmd.load_file(str(file_gro), skip_bonds=True)
    natoms = len(frame.atoms)
    out_of_range = [i for i in indexes if i < 1 or i > natoms]
    if out_of_range:
        raise QMMMError(f'{file_ndx}: [ {group} ] refers to atom {out_of_range[0]} but '
                        f'{file_gro} has only {natoms} atoms')

    atoms = [frame.atoms[i - 1] for i in indexes]
    coordinates = np.array([[a.xx, a.xy, a.xz] for a in atoms])   # ParmEd works in Angstrom
    labels = [f'{a.residue.name}{a.residue.number}:{a.name}' for a in atoms]
    if elements is None:
        elements = [data.guess_element(a.name, a.residue.name) for a in atoms]
        unknown = [labels[i] for i, e in enumerate(elements) if e is None]
        if unknown:
            if warn_unknown:
                LOGGER.warning('could not guess the element of %d atom(s), e.g. %s; the elements '
                               'of this geometry stay unknown (fine for a coordinate-only '
                               'rewrite) -- pass elements=[...] to set them',
                               len(unknown), ', '.join(unknown[:5]))
            elements = None
        else:
            counts = {}
            for element in elements:
                counts[element] = counts.get(element, 0) + 1
            LOGGER.info('guessed the elements of the %d atoms of [ %s ]: %s', len(atoms), group,
                        ', '.join(f'{n}x{e}' for e, n in sorted(counts.items())))
    elif len(elements) != len(atoms):
        raise QMMMError(f'{len(elements)} elements given for {len(atoms)} atoms')
    LOGGER.info('read %d atoms of [ %s ] from %s', len(atoms), group, file_gro.name)
    return QMGeometry(coordinates, elements, labels)


def read_index_file(path):
    """``{group name: [1-based atom numbers]}`` of a Gromacs ``.ndx``."""
    groups, current = {}, None
    for line in Path(path).expanduser().read_text().splitlines():
        stripped = line.strip()
        if stripped.startswith('['):
            current = stripped.strip('[] \t')
            groups.setdefault(current, [])
        elif current is not None and stripped:
            try:
                groups[current].extend(int(token) for token in stripped.split())
            except ValueError as exc:
                raise QMMMError(f'{path}: cannot parse a number in group [ {current} ]: '
                                f'{stripped[:40]!r}') from exc
    return groups


class QM:
    """A QM/MM system built on top of a Gromacs topology.

    Parameters
    ----------
    file_gro, file_top : path-like
        input coordinates and topology.
    file_ogro, file_otop, file_ondx : path-like
        output coordinates, topology and index file.  Parent directories are
        created when needed.
    """

    def __init__(self, file_gro, file_top, file_ogro, file_otop, file_ondx):
        # ---------------- files -------------------------------------------------
        self.file_gro = Path(file_gro).expanduser().resolve()
        self.file_top = Path(file_top).expanduser().resolve()
        self.file_ogro = Path(file_ogro).expanduser()
        self.file_otop = Path(file_otop).expanduser()
        self.file_ondx = Path(file_ondx).expanduser()

        # ---------------- selection ---------------------------------------------
        self.qm_manual_mask = ''      # user masks, in *input* atom numbering
        self.qm_extend_mask = ''      # masks grown to the breakable bonds, input numbering
        self.qm_input_mask = ''       # the two combined -- reusable inside further masks
        self.qm_mask = ''             # the QM region in *output* numbering

        # ---------------- structures --------------------------------------------
        self.itop = None              # the whole input system
        self._velocities = None       # (natoms, 3) from the input coordinate file, A/ps
        self._keep_input_idx = None   # input indices that ended up in qm_mol
        self._rest_input_idx = None   # input indices that ended up in rest
        self.qm_mol = None            # merged moleculetype: QM + polymer + QM solvent
        self.rest = None              # everything that stays untouched (bulk solvent, ions)
        self.qm = None                # view of qm_mol holding the QM atoms

        self.qm_idx = set()           # indices in qm_mol of the real QM atoms
        self.la_idx = []              # indices in qm_mol of the link atoms
        self.cp_idx = []              # indices in qm_mol of the charge points
        self.qm_group_idx = []        # QM atoms + link atoms, sorted -> [ QM ] group
        self.n_input_atoms = 0        # size of qm_mol before link atoms were added

        # ---------------- chemistry knobs ---------------------------------------
        # order is important: (QM atom name, MM atom name)
        self.breakable_bonds = set(PROTEIN_BREAKABLE_BONDS)
        self.H_dist = dict(PROTEIN_H_DIST)          # QM atom -> link atom, Angstrom
        self.solvent_and_ions = set(SOLVENT_AND_IONS)
        self.redistr_residues = set(POLYMER_RESIDUES)  # may accept redistributed charge
        self.strict_h_dist = False    # True -> never guess a link-atom distance
        # fractional part at which the MM charge sum rounds up in magnitude
        self.charge_rounding = DEFAULT_CHARGE_ROUNDING

        # ---------------- results -----------------------------------------------
        self.qmmm_bonds = []
        self.qmqm_bonds = []
        self.mm1_atoms = []
        self.vs2 = []
        self.LA_indexes = []
        self.qm_charge = None         # sum of the force-field charges of the QM atoms
        self.total_charge = None      # of the whole input system
        self.aim_qm_charge = None     # integer charge handed to the QM code
        self.count = {}               # residue name -> molecules pulled into the QM block

        # ---------------- topology bookkeeping ----------------------------------
        self._pre_lines = []          # input topology before the first [ moleculetype ]
        self._post_lines = []         # after the first moleculetype, before [ system ]
        self._system_lines = []       # [ system ] block
        self._molecule_entries = []   # [(moleculetype name, count), ...]
        self._merged_per_entry = []   # how many instances of each entry joined the QM block
        self._merged_name = None      # name of the merged moleculetype
        self._merged_nrexcl = 3

        self.read_inputs()

    # ==================================================================== input
    def read_inputs(self):
        """Load the topology and the coordinates and pre-parse the topology text."""
        for path in (self.file_gro, self.file_top):
            if not path.is_file():
                raise QMMMError(f'input file not found: {path}')
        LOGGER.info('reading %s + %s', self.file_top.name, self.file_gro.name)
        try:
            # read the coordinate file on its own so that the velocities survive:
            # load_file(top, xyz=...) keeps only the positions
            frame = pmd.load_file(str(self.file_gro), skip_bonds=True)
            if getattr(frame, 'coordinates', None) is None:
                raise QMMMError(f'{self.file_gro} carries no coordinates')
            self._velocities = getattr(frame, 'velocities', None)
            # ParmEd resolves #include relative to the directory of the topology,
            # so the module works from any working directory.
            self.itop = pmd.load_file(str(self.file_top), xyz=frame.coordinates,
                                      box=getattr(frame, 'box', None))
        except QMMMError:
            raise
        except Exception as exc:
            raise QMMMError(f'could not read {self.file_top} / {self.file_gro}: {exc}') from exc
        if len(self.itop.atoms) != len(frame.atoms):
            raise QMMMError(f'{self.file_top} has {len(self.itop.atoms)} atoms but '
                            f'{self.file_gro} has {len(frame.atoms)}')
        if self._velocities is not None:
            self._velocities = np.asarray(self._velocities).reshape(-1, 3)
            LOGGER.debug('velocities read for %d atoms', len(self._velocities))
        LOGGER.info('input system: %d atoms, %d residues',
                    len(self.itop.atoms), len(self.itop.residues))
        self._parse_topology_text()

    def _parse_topology_text(self):
        """Split the topology file into the parts that are copied verbatim.

        The body of the first ``[ moleculetype ]`` is regenerated by ParmEd, but
        everything around it -- force-field includes, water/ion includes, the
        position-restraint blocks, ``[ system ]`` -- has to survive unchanged.
        """
        lines = self.file_top.read_text().splitlines(keepends=True)

        first_mol = next((i for i, l in enumerate(lines) if _directive(l) == 'moleculetype'), None)
        if first_mol is None:
            raise QMMMError(f'{self.file_top}: no [ moleculetype ] directive found')
        self._pre_lines = lines[:first_mol]

        atoms_at = next((i for i in range(first_mol, len(lines)) if _directive(lines[i]) == 'atoms'), None)
        if atoms_at is None:
            raise QMMMError(f'{self.file_top}: the first moleculetype has no [ atoms ] section')

        body_end = self._find_body_end(lines, atoms_at)
        # the blank/comment lines in front of the include block introduce it, so
        # they belong to the tail rather than to the moleculetype body
        while body_end > atoms_at:
            previous = lines[body_end - 1].strip()
            if previous and not previous.startswith(';'):
                break
            body_end -= 1

        system_at = next((i for i in range(body_end, len(lines)) if _directive(lines[i]) == 'system'), None)
        if system_at is None:
            raise QMMMError(f'{self.file_top}: no [ system ] directive found')
        self._post_lines = lines[body_end:system_at]
        self._system_lines = lines[system_at:]
        self._parse_molecule_list()
        LOGGER.debug('topology text: %d header lines, %d lines between moleculetype and [ system ]',
                     len(self._pre_lines), len(self._post_lines))

    @staticmethod
    def _find_body_end(lines, start):
        """Index of the first line after ``start`` that is no longer part of the
        first moleculetype: an ``#include`` (possibly wrapped in ``#ifdef``) or a
        directive that cannot appear inside a moleculetype."""
        i = start
        while i < len(lines):
            stripped = lines[i].strip()
            directive = _directive(lines[i])
            if directive is not None and directive not in _MOLECULETYPE_SECTIONS:
                return i
            if stripped.startswith('#include'):
                return i
            if stripped.startswith(('#ifdef', '#ifndef')):
                # a conditional block that pulls in a file (e.g. POSRES) ends the body;
                # one that only holds parameters (e.g. FLEXIBLE water) does not
                depth, j = 1, i + 1
                has_include = False
                while j < len(lines) and depth:
                    s = lines[j].strip()
                    if s.startswith(('#ifdef', '#ifndef')):
                        depth += 1
                    elif s.startswith('#endif'):
                        depth -= 1
                    elif s.startswith('#include'):
                        has_include = True
                    j += 1
                if has_include:
                    return i
                i = j
                continue
            i += 1
        return len(lines)

    def _parse_molecule_list(self):
        """Read ``[ molecules ]`` and check it against the loaded structure."""
        entries = []
        in_molecules = False
        for line in self._system_lines:
            directive = _directive(line)
            if directive is not None:
                in_molecules = directive == 'molecules'
                continue
            if not in_molecules:
                continue
            text = line.split(';')[0].strip()
            if not text:
                continue
            fields = text.split()
            if len(fields) != 2 or not fields[1].lstrip('-').isdigit():
                raise QMMMError(f'{self.file_top}: cannot parse "[ molecules ]" line: {line.strip()!r}')
            entries.append((fields[0], int(fields[1])))
        if not entries:
            raise QMMMError(f'{self.file_top}: the [ molecules ] section is empty')

        known = getattr(self.itop, 'molecules', None)
        if not known:
            raise QMMMError('ParmEd did not expose the moleculetype templates; '
                            'a Gromacs .top (not .itp) is required')
        total = 0
        for name, count in entries:
            if name not in known:
                raise QMMMError(f'{self.file_top}: [ molecules ] refers to the unknown '
                                f'moleculetype {name!r}')
            total += count * len(known[name][0].atoms)
        if total != len(self.itop.atoms):
            raise QMMMError(f'{self.file_top}: [ molecules ] accounts for {total} atoms but the '
                            f'structure has {len(self.itop.atoms)}')
        self._molecule_entries = entries
        LOGGER.debug('[ molecules ]: %s', ', '.join(f'{n}x{c}' for n, c in entries))

    # ================================================================ selection
    def choose_qm_manually(self, qm_manual_mask):
        """Add an Amber mask (input numbering) to the QM region as it is."""
        self.qm_manual_mask = (qm_manual_mask if not self.qm_manual_mask
                               else f'({self.qm_manual_mask})|({qm_manual_mask})')
        self._update_input_mask()
        return self.qm_input_mask

    def choose_qm_to_extend(self, qm_extend_mask):
        """Add an Amber mask and grow it along bonds until a breakable bond."""
        grown = self.extend_until_break(qm_extend_mask)
        self.qm_extend_mask = grown if not self.qm_extend_mask else f'({self.qm_extend_mask})|({grown})'
        self._update_input_mask()
        return self.qm_input_mask

    def _update_input_mask(self):
        """Keep ``qm_input_mask`` in sync so it can be reused in later masks."""
        parts = [m for m in (self.qm_manual_mask, self.qm_extend_mask) if m]
        self.qm_input_mask = parts[0] if len(parts) == 1 else '|'.join(f'({p})' for p in parts)

    def dfs_extend(self, start):
        """Atoms reachable from ``start`` without stepping over a breakable bond."""
        atoms = self.itop.atoms
        visited = set()
        stack = [start]
        while stack:
            idx = stack.pop()
            if idx in visited:
                continue
            visited.add(idx)
            atom = atoms[idx]
            for partner in atom.bond_partners:
                if partner.idx not in visited and (atom.name, partner.name) not in self.breakable_bonds:
                    stack.append(partner.idx)
        return visited

    def extend_until_break(self, amber_mask):
        """Turn a mask into an explicit ``@i,j,k`` mask grown to the breakable bonds."""
        seeds = [atom.idx for atom in self._select(amber_mask)]
        if not seeds:
            raise QMMMError(f'mask {amber_mask!r} selects no atoms')
        visited = set()
        for index in seeds:
            if index not in visited:
                visited |= self.dfs_extend(index)
        LOGGER.info('extended %d seed atom(s) of %r to %d atoms', len(seeds), amber_mask, len(visited))
        return '@' + ','.join(str(i + 1) for i in sorted(visited))

    def _select(self, amber_mask):
        """Atoms of the input system matching ``amber_mask`` (a view, so
        ``atom.idx`` is the index in the input system)."""
        try:
            view = self.itop.view[amber_mask]
        except Exception as exc:
            raise QMMMError(f'invalid selection mask {amber_mask!r}: {exc}') from exc
        return list(getattr(view, 'atoms', []))

    # ================================================================ QM region
    def determine_qm(self):
        """Split the system into the merged QM moleculetype and the untouched rest.

        A molecule of the input system joins the merged moleculetype when it
        contains a QM atom, or when it may accept redistributed charge (the
        protein / nucleic acid).  Everything else -- bulk water, ions -- stays in
        its own moleculetype and keeps its ``[ molecules ]`` entry.
        """
        self._update_input_mask()
        if not self.qm_input_mask:
            raise QMMMError('the QM region is empty: call choose_qm_manually() '
                            'and/or choose_qm_to_extend() first')

        n_atoms = len(self.itop.atoms)
        qm_flag = np.zeros(n_atoms, dtype=bool)
        qm_input_idx = [atom.idx for atom in self._select(self.qm_input_mask)]
        if not qm_input_idx:
            raise QMMMError(f'mask {self.qm_input_mask!r} selects no atoms')
        qm_flag[qm_input_idx] = True
        self._merged_name = None       # determine_qm() may be called more than once

        acceptor_flag = np.fromiter(
            (atom.residue.name.upper() in self.redistr_residues for atom in self.itop.atoms),
            dtype=bool, count=n_atoms)

        keep_flag, merged_per_entry = self._merge_molecules(qm_flag, acceptor_flag)
        self._merged_per_entry = merged_per_entry
        self._check_whole_molecules(qm_flag, keep_flag)

        LOGGER.info('QM moleculetype %r: %d atoms; untouched rest: %d atoms',
                    self._merged_name, int(keep_flag.sum()), int((~keep_flag).sum()))

        self.qm_mol = self.itop[keep_flag]
        self.rest = self.itop[~keep_flag]
        self.n_input_atoms = len(self.qm_mol.atoms)
        self._keep_input_idx = np.flatnonzero(keep_flag)
        self._rest_input_idx = np.flatnonzero(~keep_flag)

        # input index -> index inside qm_mol (ParmEd keeps the original order)
        new_index = np.cumsum(keep_flag) - 1
        self.qm_idx = {int(new_index[i]) for i in qm_input_idx}
        self.qm_group_idx = sorted(self.qm_idx)
        self.qm_mask = '@' + ','.join(str(i + 1) for i in self.qm_group_idx)
        self.qm = self.qm_mol.view[self.qm_mask]
        self.la_idx, self.cp_idx = [], []

        residues = {self.itop.atoms[i].residue.idx for i in qm_input_idx}
        # how many molecules each [ molecules ] entry lost to the merged moleculetype
        self.count = {}
        for (name, _), merged in zip(self._molecule_entries, merged_per_entry):
            if merged and name != self._merged_name:
                self.count[name] = self.count.get(name, 0) + merged
        LOGGER.info('QM region: %d atoms in %d residues%s', len(self.qm_idx), len(residues),
                    (' (+ ' + ', '.join(f'{c} {n}' for n, c in sorted(self.count.items()))
                     + ' moved out of [ molecules ])') if self.count else '')

    def _merge_molecules(self, qm_flag, acceptor_flag):
        """Decide, molecule by molecule, what goes into the merged moleculetype."""
        templates = self.itop.molecules
        keep = np.zeros(len(qm_flag), dtype=bool)
        merged_per_entry = []
        first = 0
        for name, count in self._molecule_entries:
            size = len(templates[name][0].atoms)
            merged = 0
            for _ in range(count):
                stop = first + size
                if qm_flag[first:stop].any() or acceptor_flag[first:stop].any():
                    keep[first:stop] = True
                    merged += 1
                    if self._merged_name is None:
                        self._merged_name = name
                        self._merged_nrexcl = templates[name][1]
                first = stop
            merged_per_entry.append(merged)
            if merged and merged != count:
                LOGGER.info('%d of %d %s molecules move into the QM moleculetype', merged, count, name)
        if self._merged_name is None:
            raise QMMMError('no molecule ended up in the QM moleculetype -- this should not happen')
        for (name, _), merged in zip(self._molecule_entries, merged_per_entry):
            nrexcl = templates[name][1]
            if merged and nrexcl != self._merged_nrexcl and len(templates[name][0].atoms) > 3:
                LOGGER.warning('moleculetype %s has nrexcl=%d but is merged into %s (nrexcl=%d); '
                               'check the exclusions of that molecule',
                               name, nrexcl, self._merged_name, self._merged_nrexcl)
        return keep, merged_per_entry

    def _check_whole_molecules(self, qm_flag, keep_flag):
        """A molecule must be either completely inside or completely outside the
        merged block -- otherwise the coordinate file and ``[ molecules ]`` cannot
        stay consistent.  (``_merge_molecules`` guarantees this; the check is here
        to catch a future change in the selection logic.)"""
        broken = []
        for res in self.itop.residues:
            flags = {bool(keep_flag[a.idx]) for a in res.atoms}
            if len(flags) > 1:
                broken.append(f'{res.name}{res.number}')
        if broken:
            raise QMMMError('these residues are split between the QM moleculetype and the rest: '
                            + ', '.join(broken[:10]))
        partial = [f'{self.itop.atoms[i].residue.name}{self.itop.atoms[i].residue.number}'
                   for i in np.flatnonzero(qm_flag)
                   if self.itop.atoms[i].residue.name.upper() in self.solvent_and_ions
                   and not all(qm_flag[a.idx] for a in self.itop.atoms[i].residue.atoms)]
        if partial:
            LOGGER.warning('these solvent/ion residues are only partly quantum: %s -- '
                           'usually you want to select them completely (e.g. with "<:")',
                           ', '.join(sorted(set(partial))[:10]))

    # =================================================================== charge
    def calculate_charge_qm(self):
        """Sum of the force-field charges of the QM atoms."""
        self.qm_charge = sum(atom.charge for atom in self.qm.atoms)
        self.total_charge = (sum(a.charge for a in self.qm_mol.atoms)
                             + sum(a.charge for a in self.rest.atoms))
        return self.qm_charge

    def detect_qm_charge(self, rounding=None):
        """Integer QM charge implied by the force field.

        The whole system carries an integer charge and so does the MM part, so
        ``round(total) - round(total - q_QM)`` is the integer the QM code has to
        be told about.

        ``rounding`` (default :attr:`charge_rounding`, 0.25) is the fractional
        part at which the MM sum rounds up in magnitude instead of down.  Ordinary
        rounding would use 0.5; a QM/MM cut through a charge group leaves part of
        a formal charge behind on the MM side, so the sum comes out too small in
        magnitude and a lower threshold recovers it:

        ============ ============ ============ ============
        MM sum       0.25 (default) 0.5 ('nearest') 0.0 ('away')
        ============ ============ ============ ============
        1.43         2            1            2
        -1.43        -2           -1           -2
        0.04         0            0            1
        2.00         2            2            2
        ============ ============ ============ ============

        When ordinary rounding would give a different answer both candidates are
        reported: the choice is genuinely ambiguous then, and the wrong one can
        keep the SCC from converging at all.
        """
        threshold = _rounding_threshold(self.charge_rounding if rounding is None else rounding)
        q_mm = self.total_charge - self.qm_charge
        total = _round_charge(self.total_charge, 'nearest')   # this one really is an integer
        aim = total - _round_charge(q_mm, threshold)
        other = total - _round_charge(q_mm, 'nearest')
        LOGGER.info('charges: system %.4f, QM %.4f, MM %.4f -> QM charge %+d '
                    '(rounding threshold %.2f)',
                    self.total_charge, self.qm_charge, q_mm, aim, threshold)
        if other != aim:
            LOGGER.warning('the MM part carries %.3f e, so the QM charge is ambiguous: the %.2f '
                           'threshold gives %+d, ordinary rounding would give %+d. The QM/MM '
                           'boundary probably cuts through a charge group. If the SCC will '
                           'not converge, try %+d (qm_aim_charge=%d, or '
                           '"qmmmtools rewrite-hsd --charge %d" plus QMcharge in the .mdp).',
                           q_mm, threshold, aim, other, other, other, other)
        return aim

    def _charge_acceptors(self):
        """Indices in ``qm_mol`` of MM atoms that may take up excess charge."""
        acceptors = [atom.idx for atom in self.qm_mol.atoms
                     if atom.idx not in self.qm_idx
                     and atom.residue.name.upper() in self.redistr_residues]
        if acceptors:
            return acceptors
        # lipid-only or ligand-only systems have no standard polymer residues
        fallback = [atom.idx for atom in self.qm_mol.atoms
                    if atom.idx not in self.qm_idx
                    and atom.residue.name.upper() not in self.solvent_and_ions]
        if fallback:
            LOGGER.warning('no protein/nucleic-acid atom found for the charge redistribution; '
                           'falling back to all %d non-solvent MM atoms. Set '
                           'qm.redistr_residues to control this.', len(fallback))
        return fallback

    def redistribute_charge_from_qm_to_mm(self, aim_charge=None, rounding=None):
        """Shift the QM charge to ``aim_charge`` and compensate on the MM atoms.

        The QM atoms take up the difference evenly (their MM charges are not used
        by the QM code, but keeping the book straight makes the topology readable)
        and the same amount is removed from the MM acceptors, so the total charge
        of the system does not change.
        """
        if aim_charge is None:
            aim_charge = self.detect_qm_charge(rounding)
        self.aim_qm_charge = aim_charge
        delta = aim_charge - self.qm_charge
        if abs(delta) < 1e-9:
            LOGGER.info('QM charge already equals %+d, nothing to redistribute', aim_charge)
            return
        acceptors = self._charge_acceptors()
        if not acceptors:
            raise QMMMError(f'{delta:+.4f} e has to be moved out of the QM region but there is no '
                            'MM atom to take it (water and ions are never used). '
                            'Set qm.redistr_residues to the residue names that may accept charge.')
        qm_indices = sorted(self.qm_idx)
        dq_qm = delta / len(qm_indices)
        dq_mm = delta / len(acceptors)
        atoms = self.qm_mol.atoms
        for i in qm_indices:
            atoms[i].charge += dq_qm
        for i in acceptors:
            atoms[i].charge -= dq_mm
        LOGGER.info('moved %+.4f e from the QM region onto %d MM atoms (%+.2e e each)',
                    delta, len(acceptors), -dq_mm)

    def full_charge(self, system=None):
        """Total charge of ``system`` (default: the whole output system)."""
        if system is not None:
            return sum(atom.charge for atom in system.atoms)
        return sum(a.charge for a in self.qm_mol.atoms) + sum(a.charge for a in self.rest.atoms)

    # ==================================================================== bonds
    def find_qmqm_bonds(self):
        """Bonds with both atoms in the QM region."""
        self.qmqm_bonds = [b for b in self.qm_mol.bonds
                           if b.atom1.idx in self.qm_idx and b.atom2.idx in self.qm_idx]
        return self.qmqm_bonds

    def find_qmmm_bonds(self):
        """Bonds that cross the QM/MM boundary; they get a link atom each."""
        self.qmmm_bonds = []
        mm1 = []
        for bond in self.qm_mol.bonds:
            in1 = bond.atom1.idx in self.qm_idx
            in2 = bond.atom2.idx in self.qm_idx
            if in1 == in2:
                continue
            self.qmmm_bonds.append(bond)
            mm1.append(bond.atom2.idx if in1 else bond.atom1.idx)
        self.mm1_atoms = sorted(set(mm1))
        LOGGER.info('%d QM/MM boundary bond(s) on %d MM atom(s)',
                    len(self.qmmm_bonds), len(self.mm1_atoms))
        for bond in self.qmmm_bonds:
            aqm, amm = ((bond.atom1, bond.atom2) if bond.atom1.idx in self.qm_idx
                        else (bond.atom2, bond.atom1))
            LOGGER.debug('  cut %s%d:%s -- %s%d:%s', aqm.residue.name, aqm.residue.number,
                         aqm.name, amm.residue.name, amm.residue.number, amm.name)
        return self.qmmm_bonds

    def process_bonds(self):
        """Turn QM-QM bonds into Gromacs "connections" (funct 5).

        This is not optional: a plain funct-1 bond would be converted into a
        constraint by ``constraints = h-bonds`` before grompp removes the QM
        bonded terms, and the QM hydrogens would end up rigid.
        """
        self.find_qmqm_bonds()
        for bond in self.qm_mol.bonds:
            if bond.atom1.idx in self.qm_idx and bond.atom2.idx in self.qm_idx:
                bond.funct = 5
                bond.type = None      # a funct 5 bond must not carry b0/kb
        LOGGER.info('%d QM-QM bond(s) converted to connections (funct 5)', len(self.qmqm_bonds))

    # ------------------------------------------------------- MM term retention
    def _n_qm(self, *atoms):
        return sum(1 for atom in atoms if atom.idx in self.qm_idx)

    def process_mm_terms(self, mode):
        """Remove the MM bonded terms of the QM region from the written topology.

        ``mode``
            ``'no'``       leave them in place.  grompp of the DFTB+ build removes
                           them itself and reports what it removed, so this is the
                           default.
            ``'classic'``  the scheme grompp uses: a term goes as soon as all but
                           one of its atoms are QM, and the 1-4 pairs between a QM
                           atom and a directly bonded MM atom go as well.
            ``'amber'``    only terms whose atoms are *all* QM are removed.
        """
        mode = (mode or 'no').lower()
        if mode == 'no':
            LOGGER.info('MM bonded terms of the QM region are left in the topology '
                        '(grompp removes them itself)')
            return
        if mode not in ('classic', 'amber'):
            raise QMMMError(f'unknown mm_retention {mode!r}, use "no", "classic" or "amber"')

        angle_min, dihedral_min = (2, 3) if mode == 'classic' else (3, 4)
        mm1 = set(self.mm1_atoms)
        removed = {}

        def prune(container, n_atoms, minimum):
            kept, gone = [], 0
            for term in container:
                atoms = [getattr(term, f'atom{i}') for i in range(1, n_atoms + 1)]
                if self._n_qm(*atoms) >= minimum:
                    if hasattr(term, 'delete'):
                        term.delete()
                    gone += 1
                else:
                    kept.append(term)
            return pmd.TrackedList(kept), gone

        self.qm_mol.angles, removed['angles'] = prune(self.qm_mol.angles, 3, angle_min)
        self.qm_mol.dihedrals, removed['dihedrals'] = prune(self.qm_mol.dihedrals, 4, dihedral_min)
        if len(self.qm_mol.impropers):   # Gromacs keeps funct-4 impropers in .dihedrals
            self.qm_mol.impropers, removed['impropers'] = prune(self.qm_mol.impropers, 4, dihedral_min)

        kept_pairs, gone_pairs = [], 0
        for pair in self.qm_mol.adjusts:
            n_qm = self._n_qm(pair.atom1, pair.atom2)
            boundary = mode == 'classic' and n_qm == 1 and (
                pair.atom1.idx in mm1 or pair.atom2.idx in mm1)
            if n_qm == 2 or boundary:
                gone_pairs += 1
            else:
                kept_pairs.append(pair)
        self.qm_mol.adjusts = pmd.TrackedList(kept_pairs)
        removed['1-4 pairs'] = gone_pairs

        LOGGER.info("mm_retention='%s' removed %s", mode,
                    ', '.join(f'{v} {k}' for k, v in removed.items() if v) or 'nothing')

    # =============================================================== link atoms
    def _link_atom_distance(self, aqm, amm):
        """QM atom -> link atom distance for a boundary bond, in Angstrom."""
        for key in ((aqm.name, amm.name), (amm.name, aqm.name)):
            if key in self.H_dist:
                return self.H_dist[key]
        label = (f'{aqm.residue.name}{aqm.residue.number}:{aqm.name} -- '
                 f'{amm.residue.name}{amm.residue.number}:{amm.name}')
        if self.strict_h_dist:
            raise QMMMError(f'no link-atom distance for the boundary bond {label}; '
                            'add it to qm.H_dist')
        element = self._element_of(aqm)
        distance = H_DIST_BY_ELEMENT.get(element, H_DIST_FALLBACK)
        LOGGER.warning('no entry in H_dist for %s; using the %s-H default of %.2f A. '
                       'Add ("%s", "%s") to qm.H_dist to control it.',
                       label, element, distance, aqm.name, amm.name)
        return distance

    def _add_dummy(self, name, charge, xyz, res_n):
        """Append a massless site (link atom or charge point) to the QM moleculetype."""
        atom = Atom(name=name, atomic_number=1, type=name, charge=charge, mass=0)
        atom.xx, atom.xy, atom.xz = (float(v) for v in xyz)
        # without velocities ParmEd would drop the velocities of the whole system
        atom.vx = atom.vy = atom.vz = 0.0
        self.qm_mol.add_atom(atom, 'XXX', res_n, chain='')
        return atom

    def vs2_and_LA(self, link_la_to_mm1=True):
        """Place a link atom on every QM/MM bond as a two-body virtual site.

        The site sits on the QM--MM axis at ``H_dist`` from the QM atom, so its
        Gromacs construction weight is ``d / |r_QM-MM|``.  With
        ``link_la_to_mm1`` a funct-5 bond to the MM atom is added as well, which
        is what makes grompp generate the exclusions around the link atom.
        """
        self.vs2 = []
        self.LA_indexes = []
        self.la_idx = []
        res_n = len(self.qm_mol.residues)
        atom_n = len(self.qm_mol.atoms)
        xyz = self.qm_mol.coordinates
        for bond in self.qmmm_bonds:
            aqm, amm = ((bond.atom1, bond.atom2) if bond.atom1.idx in self.qm_idx
                        else (bond.atom2, bond.atom1))
            qm_xyz, mm_xyz = xyz[aqm.idx], xyz[amm.idx]
            r = float(np.linalg.norm(mm_xyz - qm_xyz))
            if r < 1e-6:
                raise QMMMError(f'atoms {aqm.idx + 1} and {amm.idx + 1} sit on top of each other; '
                                'is the coordinate file the right one for this topology?')
            d = self._link_atom_distance(aqm, amm)
            link = self._add_dummy('LA', 0.0, qm_xyz + d * (mm_xyz - qm_xyz) / r, res_n)
            self.vs2.append([str(atom_n + 1), str(aqm.idx + 1), str(amm.idx + 1),
                             '1', f'{d / r:.3f}', '; qmmm'])
            self.LA_indexes.append(atom_n)
            self.la_idx.append(atom_n)
            if link_la_to_mm1:
                bond5 = pmd.Bond(self.qm_mol.atoms[amm.idx], link)
                bond5.funct = 5
                self.qm_mol.bonds.append(bond5)
            res_n += 1
            atom_n += 1
        self.qm_group_idx = sorted(self.qm_idx | set(self.la_idx))
        self.qm_mask = '@' + ','.join(str(i + 1) for i in self.qm_group_idx)
        self.qm = self.qm_mol.view[self.qm_mask]
        LOGGER.info('%d link atom(s) added%s', len(self.la_idx),
                    ', each bonded to its MM atom (funct 5)' if link_la_to_mm1 else '')

    # ================================================== boundary charge schemes
    def _boundary_mm1_atoms(self):
        """The MM atoms of the boundary bonds, each one only once."""
        seen = set()
        for bond in self.qmmm_bonds:
            mm1 = bond.atom2 if bond.atom1.idx in self.qm_idx else bond.atom1
            if mm1.idx not in seen:
                seen.add(mm1.idx)
                yield mm1

    def _mm2_partners(self, mm1):
        """Real MM neighbours of a boundary atom (no QM atom, no dummy site)."""
        return [a for a in mm1.bond_partners
                if a.idx not in self.qm_idx and a.idx < self.n_input_atoms
                and a.name not in ('LA', 'CP')]

    def redistribute_boundary_charge(self, scheme):
        """Deal with the charge of the MM boundary atom.

        ``'no'``     keep it (the link atom then sits next to a full MM charge);
        ``'amber'``  zero it and spread the charge over the MM acceptors;
        ``'RC'``     zero it and put ``q/n`` on the midpoint of every MM1-MM2 bond;
        ``'RCD'``    like RC but with ``2q/n`` on the midpoint and ``-q/n`` on MM2,
                     which preserves the MM1-MM2 dipole;
        ``'CS'``     charge shift: ``q/n`` onto MM2 plus a ``+q/n``/``-q/n`` pair
                     around it, which preserves both charge and dipole.
        """
        scheme = (scheme or 'no').lower()
        if scheme == 'no':
            LOGGER.info('boundary charges are kept as they are (redistr_scheme="no")')
            return
        if scheme == 'amber':
            return self.amber_redist()
        if scheme == 'rc':
            return self.RC_redist()
        if scheme == 'rcd':
            return self.RCD_redist()
        if scheme == 'cs':
            return self.CS_redist()
        raise QMMMError(f'unknown redistr_scheme {scheme!r}, use "no", "amber", "RC", "RCD" or "CS"')

    def amber_redist(self):
        """Zero the boundary charges and spread them over the MM acceptors."""
        moved = 0.0
        boundary = set()
        atoms = self.qm_mol.atoms
        for mm1 in self._boundary_mm1_atoms():
            moved += atoms[mm1.idx].charge
            atoms[mm1.idx].charge = 0.0
            boundary.add(mm1.idx)
        acceptors = [i for i in self._charge_acceptors() if i not in boundary]
        if not acceptors:
            raise QMMMError('no MM atom left to accept the charge of the boundary atoms')
        dq = moved / len(acceptors)
        for i in acceptors:
            atoms[i].charge += dq
        LOGGER.info('amber scheme: %+.4f e taken from %d boundary atom(s) and spread over %d MM atoms',
                    moved, len(boundary), len(acceptors))

    def _redist_points(self, points, comment, label):
        """Shared machinery of RC / RCD / CS.

        ``points(dq)`` returns ``(fraction along MM1->MM2, charge of the point,
        charge added to MM2)`` triples.
        """
        res_n = len(self.qm_mol.residues)
        atom_n = len(self.qm_mol.atoms)
        xyz = self.qm_mol.coordinates
        atoms = self.qm_mol.atoms
        self.cp_idx = []
        moved, n_mm1 = 0.0, 0
        for mm1 in self._boundary_mm1_atoms():
            mm2_list = self._mm2_partners(mm1)
            if not mm2_list:
                LOGGER.warning('boundary atom %s%d:%s has no MM neighbour, its charge is left alone',
                               mm1.residue.name, mm1.residue.number, mm1.name)
                continue
            charge = atoms[mm1.idx].charge
            dq = charge / len(mm2_list)
            atoms[mm1.idx].charge = 0.0
            moved += charge
            n_mm1 += 1
            for mm2 in mm2_list:
                axis = xyz[mm2.idx] - xyz[mm1.idx]
                for fraction, point_charge, mm2_dq in points(dq):
                    self._add_dummy('CP', point_charge, xyz[mm1.idx] + axis * fraction, res_n)
                    atoms[mm2.idx].charge += mm2_dq
                    self.vs2.append([str(atom_n + 1), str(mm1.idx + 1), str(mm2.idx + 1),
                                     '1', f'{fraction:.3f}', comment])
                    self.cp_idx.append(atom_n)
                    res_n += 1
                    atom_n += 1
        LOGGER.info('%s scheme: %+.4f e taken from %d boundary atom(s), %d charge point(s) created',
                    label, moved, n_mm1, len(self.cp_idx))

    def RC_redist(self):
        self._redist_points(lambda dq: [(0.500, dq, 0.0)], '; RC', 'RC')

    def RCD_redist(self):
        self._redist_points(lambda dq: [(0.500, 2 * dq, -dq)], '; RCD', 'RCD')

    def CS_redist(self):
        self._redist_points(lambda dq: [(0.940, dq, dq), (1.060, -dq, 0.0)],
                            '; charge shift', 'CS')

    # ================================================================== outputs
    def write_outputs(self):
        """Write the index file, the topology and the coordinates."""
        for path in (self.file_ondx, self.file_otop, self.file_ogro):
            path.parent.mkdir(parents=True, exist_ok=True)
        self.write_ndx()
        self.write_otop()
        self.write_ogro()

    @staticmethod
    def _print_group(handle, name, indexes, per_line=15):
        handle.write(f'[ {name} ]\n')
        for start in range(0, len(indexes), per_line):
            handle.write(' '.join(str(i) for i in indexes[start:start + per_line]) + '\n')

    def write_ndx(self):
        """``[ QM ]`` (QM atoms + link atoms), the QM moleculetype and the rest."""
        n_qm_mol = len(self.qm_mol.atoms)
        with open(self.file_ondx, 'w') as handle:
            self._print_group(handle, 'QM', [i + 1 for i in self.qm_group_idx])
            self._print_group(handle, 'freeze', list(range(1, n_qm_mol + 1)))
            self._print_group(handle, 'Water_and_ions',
                              list(range(n_qm_mol + 1, n_qm_mol + len(self.rest.atoms) + 1)))
        LOGGER.info('wrote %s ([ QM ] = %d atoms)', self.file_ondx, len(self.qm_group_idx))

    def write_otop(self):
        """Write the QM/MM topology, keeping every include of the input file."""
        tmp = self.file_otop.with_suffix(self.file_otop.suffix + '.parmed')
        try:
            self.qm_mol.save(str(tmp), overwrite=True, format='gromacs', combine='all')
            body = self._extract_moleculetype_body(tmp)
        finally:
            tmp.unlink(missing_ok=True)

        src_dir, dst_dir = self.file_top.parent, self.file_otop.resolve().parent
        with open(self.file_otop, 'w') as handle:
            for line in self._pre_lines:
                handle.write(self._fix_include(line, src_dir, dst_dir))
            handle.write('[ moleculetype ]\n; Name            nrexcl\n')
            handle.write(f'{self._merged_name:<18s}{self._merged_nrexcl}\n\n')
            handle.writelines(body)
            self._write_virtual_sites(handle)
            for line in self._post_lines:
                handle.write(self._fix_include(line, src_dir, dst_dir))
            handle.writelines(self._rebuild_system_section())
        LOGGER.info('wrote %s (%d atoms in moleculetype %s)',
                    self.file_otop, len(self.qm_mol.atoms), self._merged_name)

    @staticmethod
    def _extract_moleculetype_body(path):
        """Take everything from ``[ atoms ]`` to just before ``[ system ]`` out of
        the topology ParmEd wrote; the header and the tail come from the input."""
        body, started = [], False
        with open(path) as handle:
            for line in handle:
                directive = _directive(line)
                if directive == 'atoms':
                    started = True
                elif directive == 'system':
                    break
                if started:
                    body.append(line)
        if not body:
            raise QMMMError('ParmEd wrote a topology without an [ atoms ] section')
        return body

    @staticmethod
    def _fix_include(line, src_dir, dst_dir):
        """Make a relative ``#include`` survive writing the topology elsewhere.

        Gromacs resolves an include relative to the including file, so a path that
        only works next to the *input* topology is rewritten to an absolute one.
        Includes that resolve from the output directory, or that come from GMXLIB,
        are left untouched.
        """
        match = _INCLUDE_RE.match(line)
        if not match:
            return line
        target = Path(match.group(2))
        if target.is_absolute():
            return line
        if (dst_dir / target).exists():
            return line
        source = src_dir / target
        if not source.exists():
            return line          # comes from GMXLIB or a -I path; nothing we can fix
        LOGGER.debug('rewriting #include %s -> %s', target, source.resolve())
        return f'{match.group(1)}{source.resolve()}{match.group(3)}'

    def _write_virtual_sites(self, handle):
        """The ``[ virtual_sites2 ]`` section with the link atoms and charge points."""
        if not self.vs2:
            return
        handle.write('\n[ virtual_sites2 ]\n;  ai    aj    ak funct            c0\n')
        widths = [max(len(row[i]) for row in self.vs2) for i in range(len(self.vs2[0]))]
        for row in self.vs2:
            handle.write('    ' + '     '.join(f'{item:>{widths[i]}}'
                                               for i, item in enumerate(row)) + '\n')
        handle.write('\n')

    def _rebuild_system_section(self):
        """``[ system ]`` unchanged plus a ``[ molecules ]`` with corrected counts.

        The merged moleculetype is listed first because its atoms come first in
        the coordinate file; the remaining entries keep their original order and
        lose the molecules that moved into the QM block.
        """
        out, in_molecules, written = [], False, False
        cursor = 0                      # position in self._molecule_entries
        for line in self._system_lines:
            directive = _directive(line)
            if directive is not None:
                in_molecules = directive == 'molecules'
                out.append(line)
                if in_molecules:
                    out.append(f'; Compound        #mols\n{self._merged_name:<18s}1\n')
                    written = True
                continue
            if not in_molecules:
                out.append(line)
                continue
            if not line.split(';')[0].strip():
                continue
            name, count = self._molecule_entries[cursor]
            left = count - self._merged_per_entry[cursor]
            cursor += 1
            if left:
                out.append(f'{name:<18s}{left}\n')
            elif cursor > 1 or count > 1:
                LOGGER.info('the [ molecules ] entry "%s %d" is now covered by the merged '
                            'moleculetype %s', name, count, self._merged_name)
        if not written:
            raise QMMMError('the input topology has no [ molecules ] section to rebuild')
        return out

    def write_ogro(self):
        """Coordinates in exactly the order of the topology: QM block, then rest.

        Written directly instead of through ParmEd: ``Structure.save`` reorders the
        atoms per molecule unless ``combine='all'`` is given, which silently breaks
        the correspondence with the topology and the index file.
        """
        atoms = list(self.qm_mol.atoms) + list(self.rest.atoms)
        n_qm_mol, n_qm_res = len(self.qm_mol.atoms), len(self.qm_mol.residues)
        velocities = self._output_velocities(len(atoms))
        with open(self.file_ogro, 'w') as handle:
            handle.write(f'QM/MM system prepared by QMMMtools: '
                         f'{len(self.qm_group_idx)} QM atoms\n')
            handle.write(f'{len(atoms):5d}\n')
            for n, atom in enumerate(atoms):
                res = atom.residue
                resid = res.idx + 1 if n < n_qm_mol else res.idx + 1 + n_qm_res
                handle.write(f'{resid % 100000:5d}{res.name[:5]:<5s}{atom.name[:5]:>5s}'
                             f'{(n + 1) % 100000:5d}'
                             f'{atom.xx / 10:8.3f}{atom.xy / 10:8.3f}{atom.xz / 10:8.3f}')
                if velocities is not None:
                    v = velocities[n]
                    handle.write(f'{v[0] / 10:8.4f}{v[1] / 10:8.4f}{v[2] / 10:8.4f}')
                handle.write('\n')
            handle.write(self._box_line(atoms))
        LOGGER.info('wrote %s (%d atoms%s)', self.file_ogro, len(atoms),
                    ', velocities kept' if velocities is not None else ', no velocities')

    def _output_velocities(self, n_out):
        """Velocities in output order; the massless sites get zero."""
        if self._velocities is None:
            return None
        n_dummies = len(self.qm_mol.atoms) - self.n_input_atoms
        return np.concatenate([
            self._velocities[self._keep_input_idx],
            np.zeros((n_dummies, 3)),
            self._velocities[self._rest_input_idx],
        ])[:n_out]

    def _box_line(self, atoms):
        """The Gromacs box line, triclinic when the input box is."""
        box = self.itop.box
        if box is None:
            crd = np.array([[a.xx, a.xy, a.xz] for a in atoms])
            diff = (crd.max(axis=0) - crd.min(axis=0)) / 10 + 0.5
            LOGGER.warning('the input has no box; writing a bounding box with 0.5 nm padding')
            return f'{diff[0]:10.5f}{diff[1]:10.5f}{diff[2]:10.5f}\n'
        a, b, c = reduce_box_vectors(*box_lengths_and_angles_to_vectors(*box))
        if all(abs(x - 90) < 1e-5 for x in box[3:]):
            return f'{a[0] / 10:10.5f}{b[1] / 10:10.5f}{c[2] / 10:10.5f}\n'
        return (f'{a[0] / 10:10.5f}{b[1] / 10:10.5f}{c[2] / 10:10.5f}'
                f'{a[1] / 10:10.5f}{a[2] / 10:10.5f}{b[0] / 10:10.5f}'
                f'{b[2] / 10:10.5f}{c[0] / 10:10.5f}{c[1] / 10:10.5f}\n')

    # ============================================================ DFTB+ / xTB
    def _element_of(self, atom):
        """Chemical element of an atom, from the atomic number ParmEd assigned."""
        if atom.atomic_number and atom.atomic_number in ELEMENTS:
            return ELEMENTS[atom.atomic_number]
        # massless sites and hand-made types may have no usable atomic number
        for key in (atom.type, atom.name):
            if key in TYPE2ELEMENT:
                return TYPE2ELEMENT[key]
        raise QMMMError(
            f'cannot tell the element of atom {atom.idx + 1} '
            f'({atom.residue.name}{atom.residue.number}:{atom.name}, type {atom.type!r}): '
            'ParmEd found no atomic number. Fix the at.num column of [ atomtypes ] '
            'or add the type to QMMMtools.data.TYPE2ELEMENT.')

    def qm_elements(self):
        """Distinct elements of the QM group, in first-appearance order."""
        return self.qm_geometry().type_names

    def qm_geometry(self):
        """The QM group as a :class:`QMGeometry`, in ``[ QM ]`` order."""
        if self.qm is None:
            raise QMMMError('there is no QM region yet: run determine_qm() or job() first')
        xyz = self.qm_mol.coordinates
        atoms = list(self.qm.atoms)
        return QMGeometry(
            coordinates=[xyz[atom.idx] for atom in atoms],
            elements=[self._element_of(atom) for atom in atoms],
            labels=[f'{a.residue.name}{a.residue.number}:{a.name}' for a in atoms])

    def make_hsd(self, file_hsd, method=DEFAULT_QM_METHOD, skpath=None, charge=None, **kwargs):
        """Write ``dftb_in.hsd`` for the QM region.

        Parameters
        ----------
        method : str
            key of :data:`QMMMtools.data.QM_METHODS`; ``'dftb3-d4'`` reproduces the
            original setup, ``'dftb3-d3h5'`` and ``'gfn2-xtb'`` are the other
            common choices.  :func:`list_methods` prints them all.
        skpath : path-like
            directory with the Slater-Koster files (DFTB methods only).  Give the
            path as mdrun will see it -- a relative one is fine.
        charge : int, optional
            overrides the charge worked out by :meth:`job`.
        **kwargs
            passed on to :func:`write_hsd` / :func:`hamiltonian_block`:
            ``sk_suffix``, ``sk_separator``, ``sk_lowercase``, ``scc_tolerance``,
            ``max_scc_iterations``, ``mixer``, ``analysis``, ``options``, ``extra``.
            Charged metal sites often need ``mixer='anderson'``.
        """
        if charge is None:
            charge = self.aim_qm_charge
        if charge is None:
            raise QMMMError('the QM charge is unknown: run job() first, or pass charge=')
        return write_hsd(file_hsd, self.qm_geometry(), charge,
                         method=method, skpath=skpath, **kwargs)

    def rewrite_hsd(self, file_hsd, source_hsd=None, geometry=True, **kwargs):
        """Update an existing ``dftb_in.hsd`` from this QM region.

        Thin wrapper around :func:`rewrite_hsd`; ``geometry=False`` changes only
        the settings and leaves the coordinates alone.  See that function for the
        full list of things that can be replaced.
        """
        if kwargs.get('charge') is None and kwargs.get('method') is not None:
            kwargs['charge'] = self.aim_qm_charge
        return rewrite_hsd(file_hsd, geometry=self.qm_geometry() if geometry else None,
                           source_hsd=source_hsd, **kwargs)

    get_method = staticmethod(get_method)
    list_methods = staticmethod(list_methods)

    # ================================================================== checks
    def check_consistency(self, file_hsd=None):
        """Re-read the written files and verify that the QM atoms line up.

        Compares the ``[ QM ]`` group of the index file with the coordinate file
        and, when given, with the geometry of the ``.hsd``.  Raises on any
        mismatch, so it is worth calling before every production run.
        """
        gro_lines = Path(self.file_ogro).read_text().splitlines()
        n_atoms = int(gro_lines[1])
        if n_atoms != len(self.qm_mol.atoms) + len(self.rest.atoms):
            raise QMMMError(f'{self.file_ogro} holds {n_atoms} atoms, expected '
                            f'{len(self.qm_mol.atoms) + len(self.rest.atoms)}')

        groups = read_index_file(self.file_ondx)
        if 'QM' not in groups:
            raise QMMMError(f'{self.file_ondx} has no [ QM ] group')
        qm_group = groups['QM']
        expected = [i + 1 for i in self.qm_group_idx]
        if qm_group != expected:
            raise QMMMError(f'the [ QM ] group of {self.file_ondx} does not match the selection')

        for number, atom in zip(qm_group, self.qm.atoms):
            line = gro_lines[1 + number]
            gro_name, gro_res = line[10:15].strip(), line[5:10].strip()
            if gro_name != atom.name[:5] or gro_res != atom.residue.name[:5]:
                raise QMMMError(f'atom {number} is {gro_res}:{gro_name} in {self.file_ogro.name} '
                                f'but {atom.residue.name}:{atom.name} in the topology')

        if file_hsd is not None:
            types, species = _read_geometry_types(HsdFile.read(file_hsd))
            if not types or not species:
                raise QMMMError(f'{file_hsd}: no usable Geometry block found')
            if len(species) != len(self.qm.atoms):
                raise QMMMError(f'{file_hsd} holds {len(species)} atoms, the QM group has '
                                f'{len(self.qm.atoms)}')
            for kind, atom in zip(species, self.qm.atoms):
                if types[kind - 1] != self._element_of(atom):
                    raise QMMMError(
                        f'{file_hsd}: atom {atom.idx + 1} ({atom.residue.name}:{atom.name}) is '
                        f'{types[kind - 1]} in the .hsd but {self._element_of(atom)} in the topology')
            LOGGER.info('consistency check passed: .ndx, .gro and %s agree on all %d QM atoms',
                        Path(file_hsd).name, len(qm_group))
        else:
            LOGGER.info('consistency check passed: .ndx and .gro agree on all %d QM atoms',
                        len(qm_group))
        return True

    # ===================================================================== run
    def job(self, qm_aim_charge=None, mm_retention='no', redistr_scheme='no',
            link_la_to_mm1=True, charge_rounding=None):
        """Build the QM/MM system and write the topology, coordinates and index file.

        Parameters
        ----------
        qm_aim_charge : int or None
            total charge of the QM region.  ``None`` (the default) derives it from
            the force-field charges, see :meth:`detect_qm_charge`.
        charge_rounding : float or {'default', 'nearest', 'away'}
            fractional part at which the derived charge rounds up in magnitude;
            0.25 by default, 0.5 is ordinary rounding.  See
            :meth:`detect_qm_charge`.
        mm_retention : {'no', 'classic', 'amber'}
            whether the MM bonded terms of the QM region are already removed from
            the written topology, see :meth:`process_mm_terms`.
        redistr_scheme : {'no', 'amber', 'RC', 'RCD', 'CS'}
            what happens to the charge of the MM boundary atoms, see
            :meth:`redistribute_boundary_charge`.
        link_la_to_mm1 : bool
            add a funct-5 bond between each link atom and its MM atom so that
            grompp generates the exclusions around the link atom.
        """
        self.determine_qm()
        self.calculate_charge_qm()
        self.redistribute_charge_from_qm_to_mm(aim_charge=qm_aim_charge,
                                               rounding=charge_rounding)
        self.find_qmmm_bonds()
        self.process_mm_terms(mm_retention)
        self.process_bonds()
        self.vs2_and_LA(link_la_to_mm1=link_la_to_mm1)
        self.redistribute_boundary_charge(redistr_scheme)
        self.write_outputs()
        LOGGER.info('done: total charge of the written system %+.4f e', self.full_charge())
        return self

    # ------------------------------------------------- backwards compatibility
    @property
    def o_qm_protein(self):
        """Deprecated alias of :attr:`qm_mol`."""
        return self.qm_mol

    @property
    def i_qm_protein(self):
        """Deprecated alias of :attr:`qm_mol`."""
        return self.qm_mol

    @property
    def o_rest(self):
        """Deprecated alias of :attr:`rest`."""
        return self.rest
