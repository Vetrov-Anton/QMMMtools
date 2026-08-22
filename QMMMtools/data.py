"""Reference tables for :mod:`QMMMtools`.

Everything that has to be adapted when a new kind of system (protein, nucleic
acid, lipid, ...) or a new QM method is used lives here, so that ``QMMMtools.core``
itself stays free of chemistry-specific tables.

The three groups of tables are

* residue-name sets -- used to tell water/ions from the polymer that may accept
  redistributed charge;
* link-atom tables -- which bonds may be cut and how far from the QM atom the
  link atom is placed;
* QM-method definitions -- the ``Hamiltonian`` block that goes into
  ``dftb_in.hsd``.
"""

# ---------------------------------------------------------------------------
# residue names
# ---------------------------------------------------------------------------

#: water, under the many names the common force fields / tools use for it
WATER_RESIDUES = frozenset({
    'SOL', 'HOH', 'WAT', 'H2O', 'OH2', 'DOD', 'D2O',
    'TIP', 'TIP2', 'TIP3', 'TIP4', 'TIP5', 'TIP3P', 'TIP4P', 'TIP5P',
    'T3P', 'T4P', 'T5P', 'TP3', 'TP4', 'TP5',
    'SPC', 'SPCE', 'SPCF', 'OPC', 'OPC3', 'FBW', 'WATER',
})

#: monoatomic ions; the names cover Amber, CHARMM and Gromacs conventions
ION_RESIDUES = frozenset({
    'NA', 'NA+', 'SOD', 'CL', 'CL-', 'CLA', 'K', 'K+', 'POT',
    'LI', 'LI+', 'RB', 'RB+', 'CS', 'CS+', 'CES', 'F', 'F-', 'BR', 'BR-', 'I', 'I-',
    'MG', 'MG2', 'MG2+', 'CA', 'CA2', 'CA2+', 'CAL', 'ZN', 'ZN2', 'ZN2+',
    'FE', 'FE2', 'FE3', 'CU', 'CU1', 'CU2', 'MN', 'MN2', 'CO', 'NI', 'CD', 'HG',
    'SR', 'BA', 'AL', 'IB+', 'IOD', 'ACT',
})

SOLVENT_AND_IONS = WATER_RESIDUES | ION_RESIDUES


def _with_termini(names):
    """Amber/CHARMM also write terminal residues as N<XXX> / C<XXX>."""
    out = set(names)
    for name in names:
        out.add('N' + name)
        out.add('C' + name)
    return frozenset(out)


#: the 20 standard residues plus the usual protonation / disulfide variants
AMINO_ACIDS = _with_termini({
    'ALA', 'ARG', 'ASN', 'ASP', 'ASH', 'CYS', 'CYX', 'CYM', 'GLN', 'GLU', 'GLH',
    'GLY', 'HIS', 'HID', 'HIE', 'HIP', 'HSD', 'HSE', 'HSP', 'ILE', 'LEU', 'LYS',
    'LYN', 'MET', 'PHE', 'PRO', 'HYP', 'SER', 'THR', 'TRP', 'TYR', 'VAL',
    'ACE', 'NME', 'NHE', 'NH2', 'FOR',
})

#: DNA / RNA, in the Amber (DA/DA3/DA5/RA...) and plain (A/T/G/C/U) spellings
NUCLEIC_ACIDS = frozenset({
    'DA', 'DT', 'DG', 'DC', 'DU', 'DI',
    'DA3', 'DT3', 'DG3', 'DC3', 'DU3', 'DA5', 'DT5', 'DG5', 'DC5', 'DU5',
    'DAN', 'DTN', 'DGN', 'DCN', 'DUN',
    'RA', 'RU', 'RG', 'RC', 'RT',
    'RA3', 'RU3', 'RG3', 'RC3', 'RA5', 'RU5', 'RG5', 'RC5',
    'RAN', 'RUN', 'RGN', 'RCN',
    'A', 'T', 'G', 'C', 'U', 'I',
    'ADE', 'THY', 'GUA', 'CYT', 'URA',
})

#: default acceptors of the charge that is moved out of the QM region
POLYMER_RESIDUES = AMINO_ACIDS | NUCLEIC_ACIDS


# ---------------------------------------------------------------------------
# link atoms
# ---------------------------------------------------------------------------
# ``breakable_bonds`` is a set of *directed* (QM atom name, MM atom name) pairs.
# The depth-first expansion stops as soon as it would step from the QM name onto
# the MM name, so the direction encodes "this side stays quantum".
# ``h_dist`` gives the QM atom -> link atom distance in Angstrom for the same
# pairs; it is looked up in both directions.
#
# The link atom caps the *QM* atom of the cut bond, so the distance is the
# equilibrium X-H bond length of that atom's element -- see H_DIST_BY_ELEMENT.
# Writing anything else there shortens or stretches a real chemical bond inside
# the QM calculation and shows up directly in the QM energy and forces.

#: backbone cuts of a protein: side chain only (CB-CA) or a peptide unit
PROTEIN_BREAKABLE_BONDS = {('CB', 'CA'), ('CA', 'CB'), ('N', 'CA'),
                           ('CA', 'C'), ('C', 'CA')}
PROTEIN_H_DIST = {
    ('CB', 'CA'): 1.09,   # C(sp3)-H, caps the side-chain CB
    ('CA', 'CB'): 1.09,   # C(sp3)-H, caps the alpha carbon (backbone kept quantum)
    ('N', 'CA'): 1.01,    # N-H, caps the amide nitrogen
    ('CA', 'C'): 1.09,    # C(sp3)-H, caps the alpha carbon
    ('C', 'CA'): 1.09,    # C(sp2)-H, caps the carbonyl carbon
}

#: nucleic acids: glycosidic bond (base only) and the sugar-phosphate ester bonds
NUCLEIC_BREAKABLE_BONDS = {
    ("N9", "C1'"), ("N1", "C1'"),           # base -> sugar (purine / pyrimidine)
    ("C1'", "N9"), ("C1'", "N1"),           # sugar -> base
    ("C5'", "O5'"), ("O5'", "C5'"),         # 5' ester
    ("C3'", "O3'"), ("O3'", "C3'"),         # 3' ester
}
NUCLEIC_H_DIST = {
    ("N9", "C1'"): 1.01, ("N1", "C1'"): 1.01,
    ("C1'", "N9"): 1.09, ("C1'", "N1"): 1.09,
    ("C5'", "O5'"): 1.09, ("O5'", "C5'"): 0.97,
    ("C3'", "O3'"): 1.09, ("O3'", "C3'"): 0.97,
}

#: lipids have no agreed atom naming, so this is only a starting point for the
#: glycerol/ester region of CHARMM-style names.  Extend it for your force field.
LIPID_BREAKABLE_BONDS = {
    ('C1', 'C2'), ('C2', 'C1'), ('C2', 'C3'), ('C3', 'C2'),
    ('O21', 'C21'), ('C21', 'O21'), ('O31', 'C31'), ('C31', 'O31'),
}
LIPID_H_DIST = {pair: (0.97 if pair[0].startswith('O') else 1.09)
                for pair in LIPID_BREAKABLE_BONDS}

#: equilibrium X-H bond lengths in Angstrom; the link-atom distance of a cut is
#: the value for the element of the QM atom that is capped.  Used as the
#: fall-back (with a warning) when a boundary bond is not listed in h_dist.
H_DIST_BY_ELEMENT = {'C': 1.09, 'N': 1.01, 'O': 0.97, 'S': 1.34, 'P': 1.42, 'Si': 1.48}
H_DIST_FALLBACK = 1.09


# ---------------------------------------------------------------------------
# elements
# ---------------------------------------------------------------------------

ELEMENTS = {
    1: "H", 2: "He", 3: "Li", 4: "Be", 5: "B", 6: "C", 7: "N", 8: "O", 9: "F", 10: "Ne",
    11: "Na", 12: "Mg", 13: "Al", 14: "Si", 15: "P", 16: "S", 17: "Cl", 18: "Ar", 19: "K", 20: "Ca",
    21: "Sc", 22: "Ti", 23: "V", 24: "Cr", 25: "Mn", 26: "Fe", 27: "Co", 28: "Ni", 29: "Cu", 30: "Zn",
    31: "Ga", 32: "Ge", 33: "As", 34: "Se", 35: "Br", 36: "Kr", 37: "Rb", 38: "Sr", 39: "Y", 40: "Zr",
    41: "Nb", 42: "Mo", 43: "Tc", 44: "Ru", 45: "Rh", 46: "Pd", 47: "Ag", 48: "Cd", 49: "In", 50: "Sn",
    51: "Sb", 52: "Te", 53: "I", 54: "Xe", 55: "Cs", 56: "Ba", 57: "La", 58: "Ce", 59: "Pr", 60: "Nd",
    61: "Pm", 62: "Sm", 63: "Eu", 64: "Gd", 65: "Tb", 66: "Dy", 67: "Ho", 68: "Er", 69: "Tm", 70: "Yb",
    71: "Lu", 72: "Hf", 73: "Ta", 74: "W", 75: "Re", 76: "Os", 77: "Ir", 78: "Pt", 79: "Au", 80: "Hg",
    81: "Tl", 82: "Pb", 83: "Bi", 84: "Po", 85: "At", 86: "Rn", 87: "Fr", 88: "Ra", 89: "Ac", 90: "Th",
    91: "Pa", 92: "U", 93: "Np", 94: "Pu", 95: "Am", 96: "Cm", 97: "Bk", 98: "Cf", 99: "Es", 100: "Fm",
}

#: Fall-back force-field-type -> element map.  It is only consulted for atoms
#: whose atomic number ParmEd could not work out (atomic_number <= 0), which in
#: practice means massless sites and exotic hand-made types.  Prefer fixing the
#: ``at.num`` column of your ``[ atomtypes ]`` over extending this table.
TYPE2ELEMENT = {
    'LA': 'H', 'CP': 'H', 'MW': 'O', 'MNH3': 'N', 'MCH3': 'C', 'LP': 'O', 'EP': 'O',
    'HW': 'H', 'OW': 'O', 'DUM': 'H',
}

#: massless / virtual sites, recognised by their atom name alone
SPECIAL_SITE_ELEMENTS = {'LA': 'H', 'CP': 'H', 'MW': 'O', 'DUM': 'H', 'LP': 'O', 'EP': 'O'}

ELEMENT_SYMBOLS = frozenset(ELEMENTS.values())

#: two-letter elements that may appear as a complete atom name.  They are only
#: accepted outside the standard residues, because CA/CD/CE/HG/... are ordinary
#: protein atom names, not calcium/cadmium/cerium/mercury.
STANDALONE_TWO_LETTER = frozenset({
    'CL', 'BR', 'MG', 'NA', 'ZN', 'FE', 'MN', 'CU', 'LI', 'RB', 'CS', 'SE', 'SI',
    'AL', 'NI', 'CO', 'CD', 'HG', 'CA', 'SR', 'BA', 'PT', 'AU', 'AG',
})

_STANDARD_RESIDUES = AMINO_ACIDS | NUCLEIC_ACIDS


def guess_element(atom_name, residue_name=''):
    """Best guess of the chemical element from a PDB/GRO style atom name.

    Only needed when no topology is available (``rewrite_hsd`` from a bare
    ``.gro`` + ``.ndx``).  Whenever a topology is at hand the atomic number that
    ParmEd reads from ``[ atomtypes ]`` is used instead, which is always right.

    The tricky part is that ``CA``, ``CD``, ``CE``, ``HG``, ``NA`` and friends are
    ordinary protein atom names as well as element symbols, so a two-letter
    reading is only accepted outside the standard residues.

    Returns ``None`` when nothing sensible can be deduced.
    """
    name = str(atom_name).strip().upper()
    residue = str(residue_name).strip().upper()
    if not name:
        return None
    if name in SPECIAL_SITE_ELEMENTS:
        return SPECIAL_SITE_ELEMENTS[name]

    letters = ''.join(c for c in name if c.isalpha())
    if not letters:
        return None

    # a monoatomic ion: the residue name usually is the element
    if residue in ION_RESIDUES:
        for candidate in (letters, ''.join(c for c in residue if c.isalpha())):
            symbol = candidate[:2].capitalize()
            if candidate[:2] in STANDALONE_TWO_LETTER and symbol in ELEMENT_SYMBOLS:
                return symbol
            if candidate[:1].capitalize() in ELEMENT_SYMBOLS:
                return candidate[:1].capitalize()

    # selenomethionine and friends: SE really is selenium even in a residue that
    # otherwise follows the protein naming
    if letters == 'SE' and residue in ('MSE', 'SEC', 'CSE'):
        return 'Se'

    if residue not in _STANDARD_RESIDUES and letters[:2] in STANDALONE_TWO_LETTER \
            and letters == name:      # the whole name, no trailing digits
        return letters.capitalize()

    first = letters[0].capitalize()
    return first if first in ELEMENT_SYMBOLS else None


# ---------------------------------------------------------------------------
# QM methods
# ---------------------------------------------------------------------------

#: Maximum angular momenta of the 3ob-3-1 set, verbatim from the README that
#: ships with it.  These fifteen elements are the whole of 3ob: anything else
#: needs a different Slater-Koster set.
MAX_ANGULAR_MOMENTUM_3OB = {
    'Br': 'd', 'C': 'p', 'Ca': 'p', 'Cl': 'd', 'F': 'p', 'H': 's', 'I': 'd',
    'K': 'p', 'Mg': 'p', 'N': 'p', 'Na': 'p', 'O': 'p', 'P': 'd', 'S': 'd',
    'Zn': 'd',
}

#: Elements that 3ob does not cover, for use with other Slater-Koster sets
#: (mio, matsci, trans3d, ...).  Writing one of these with skpath pointing at
#: 3ob produces an input DFTB+ rejects, because the SK file does not exist.
MAX_ANGULAR_MOMENTUM_OTHER = {
    'Fe': 'd', 'Cu': 'd', 'Ti': 'd', 'Ni': 'd', 'Co': 'd', 'Mn': 'd',
    'Li': 'p', 'Rb': 'p', 'Cs': 'p', 'Si': 'p', 'Al': 'p', 'B': 'p',
}

MAX_ANGULAR_MOMENTUM = {**MAX_ANGULAR_MOMENTUM_3OB, **MAX_ANGULAR_MOMENTUM_OTHER}

#: DFTB3 Hubbard derivatives (atomic units) of the 3ob-3-1 set, verbatim from
#: its README.  A DFTB3 method refuses to write an element that is missing here,
#: because a wrong derivative is worse than no calculation.
HUBBARD_DERIVS = {
    'Br': -0.0573, 'C': -0.1492, 'Ca': -0.0340, 'Cl': -0.0697, 'F': -0.1623,
    'H': -0.1857, 'I': -0.0433, 'K': -0.0339, 'Mg': -0.0200, 'N': -0.1535,
    'Na': -0.0454, 'O': -0.1575, 'P': -0.1400, 'S': -0.1100, 'Zn': -0.0300,
}


class QMMethod:
    """One entry of :data:`QM_METHODS`.

    Parameters
    ----------
    kind : {'dftb', 'xtb'}
        ``dftb`` needs Slater-Koster files, angular momenta and (for DFTB3)
        Hubbard derivatives; ``xtb`` is self-contained.
    third_order : bool
        write ``ThirdOrderFull = Yes`` + ``HubbardDerivs`` (i.e. DFTB3).
    hcorrection, dispersion : str or None
        ready-made HSD sub-blocks, indented by four spaces when written.
    xtb_method : str or None
        value of ``Method`` inside ``Hamiltonian = xTB``.
    """

    def __init__(self, name, kind, description, third_order=False,
                 hcorrection=None, dispersion=None, xtb_method=None):
        self.name = name
        self.kind = kind
        self.description = description
        self.third_order = third_order
        self.hcorrection = hcorrection
        self.dispersion = dispersion
        self.xtb_method = xtb_method

    @property
    def needs_slater_koster(self):
        return self.kind == 'dftb'

    def __repr__(self):
        return f'<QMMethod {self.name}: {self.description}>'


_D4_3OB = """Dispersion = DftD4 {
      s6 = 1
      s8 = 0.4727337
      s9 = 0
      s10 = 0
      a1 = 0.5467502
      a2 = 4.4955068
    }"""

_D3_ZERO_H5 = """Dispersion = DftD3 {
      Damping = ZeroDamping {
        sr6 = 1.25
        alpha6 = 29.61
      }
      s6 = 1.0
      s8 = 0.49
      HHRepulsion = Yes
    }"""

# DFTB3-D3(BJ) parameters as published with 3ob-3-1 (README of the set,
# J. Chem. Theory Comput. 2015, 11, 332): a1 = 0.746, a2 = 4.191, s8 = 3.209
_D3_BJ = """Dispersion = DftD3 {
      Damping = BeckeJohnson {
        a1 = 0.746
        a2 = 4.191
      }
      s6 = 1.0
      s8 = 3.209
    }"""

# zeta of the gamma^h function.  The 3ob-3-1 README prescribes 4.00 for every
# DFTB3/3OB calculation; 4.05 in the DFTB+ manual is only an illustration.
_H_DAMPING = """HCorrection = Damping {
      Exponent = 4.00
    }"""

_H5 = "HCorrection = H5 {}"

#: Everything below has been checked against DFTB+ 21.2 (the version inside the
#: gmx+DFTB+ container) -- each one runs a single point without complaint.
QM_METHODS = {
    'dftb3-d4': QMMethod(
        'dftb3-d4', 'dftb', 'DFTB3/3ob with D4 dispersion and H damping',
        third_order=True, hcorrection=_H_DAMPING, dispersion=_D4_3OB),
    'dftb3-d3h5': QMMethod(
        'dftb3-d3h5', 'dftb',
        'DFTB3/3ob with D3(zero) dispersion and the H5 hydrogen-bond correction (default)',
        third_order=True, hcorrection=_H5, dispersion=_D3_ZERO_H5),
    'dftb3-d3bj': QMMethod(
        'dftb3-d3bj', 'dftb', 'DFTB3/3ob with D3(BJ) dispersion and H damping',
        third_order=True, hcorrection=_H_DAMPING, dispersion=_D3_BJ),
    'dftb3': QMMethod(
        'dftb3', 'dftb', 'DFTB3/3ob with H damping, no dispersion',
        third_order=True, hcorrection=_H_DAMPING),
    'dftb2': QMMethod(
        'dftb2', 'dftb', 'plain SCC-DFTB (DFTB2, e.g. the mio set), no third order'),
    'gfn2-xtb': QMMethod(
        'gfn2-xtb', 'xtb', 'GFN2-xTB through the tblite library, no Slater-Koster files',
        xtb_method='GFN2-xTB'),
    'gfn1-xtb': QMMethod(
        'gfn1-xtb', 'xtb', 'GFN1-xTB through the tblite library, no Slater-Koster files',
        xtb_method='GFN1-xTB'),
}

DEFAULT_QM_METHOD = 'dftb3-d3h5'
