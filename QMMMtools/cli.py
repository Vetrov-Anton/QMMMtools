"""Command line front end of QMMMtools.

    qmmmtools prepare      build a QM/MM topology, coordinates, index file and dftb_in.hsd
    qmmmtools rewrite-hsd  update an existing dftb_in.hsd (needs only a .gro and a .ndx)
    qmmmtools check        verify that a .gro, a .ndx and a .hsd describe the same QM atoms
    qmmmtools methods      list the available QM methods
    qmmmtools tables       show the link-atom and residue-name tables

Every option has a Python equivalent; see the README and ``examples.ipynb`` if you
need more control than the command line offers.
"""

import argparse
import sys
from pathlib import Path

from . import data
from .data import DEFAULT_QM_METHOD
from .core import (DEFAULT_CHARGE_ROUNDING, DEFAULT_DFTBPLUS_VERSION, LOGGER, QM, QMMMError, HsdFile,
                   _read_geometry_types, _rounding_threshold, get_method, list_methods,
                   read_index_file, read_qm_geometry, rewrite_hsd, set_log_level)


def _rounding(text):
    """argparse type: a threshold between 0 and 1, or one of the alias names."""
    try:
        return _rounding_threshold(text)
    except QMMMError as error:
        raise argparse.ArgumentTypeError(str(error)) from None

BOND_PRESETS = {
    'protein': (data.PROTEIN_BREAKABLE_BONDS, data.PROTEIN_H_DIST),
    'nucleic': (data.NUCLEIC_BREAKABLE_BONDS, data.NUCLEIC_H_DIST),
    'lipid': (data.LIPID_BREAKABLE_BONDS, data.LIPID_H_DIST),
}


# ---------------------------------------------------------------- arguments
def _add_hsd_options(parser, with_method_default=False):
    """Options shared by ``prepare`` and ``rewrite-hsd``."""
    group = parser.add_argument_group('QM method')
    group.add_argument('--method', default=DEFAULT_QM_METHOD if with_method_default else None,
                       help=f'QM method, default {DEFAULT_QM_METHOD} (see "qmmmtools methods")')
    group.add_argument('--skpath', help='directory with the Slater-Koster files, as mdrun '
                                        'will see it (DFTB methods only)')
    group.add_argument('--sk-suffix', default='-c.spl', help='Slater-Koster file suffix')
    group.add_argument('--scc-tolerance', help='SCC convergence threshold, e.g. 1e-6')
    group.add_argument('--max-scc-iterations', type=int,
                       help='maximum number of SCC iterations')
    group.add_argument('--mixer', help='"anderson" (recommended for charged metal sites), '
                                       '"broyden", or a raw HSD block')
    group.add_argument('--dftbplus-version', metavar='V',
                       help='DFTB+ release the input is written for: 24.1 and later spell '
                            'the Analysis switch "PrintForces", 21.x-23.x "CalculateForces" '
                            f'(default {DEFAULT_DFTBPLUS_VERSION})')


def build_parser():
    # -v/-q are accepted both before and after the subcommand; SUPPRESS keeps the
    # subparser from overwriting a value that was given in front of it
    common = argparse.ArgumentParser(add_help=False)
    verbosity = common.add_mutually_exclusive_group()
    verbosity.add_argument('-v', '--verbose', action='store_true', default=argparse.SUPPRESS,
                           help='report every decision')
    verbosity.add_argument('-q', '--quiet', action='store_true', default=argparse.SUPPRESS,
                           help='warnings and errors only')

    parser = argparse.ArgumentParser(
        prog='qmmmtools', description=__doc__.split('\n\n')[0], parents=[common],
        formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    parser.add_argument('-V', '--version', action='store_true', help='print the version and exit')
    sub = parser.add_subparsers(dest='command')

    # ---- prepare ---------------------------------------------------------
    prepare = sub.add_parser('prepare', help='build a complete QM/MM system', parents=[common],
                             description='Carve a QM region out of an MM system and write '
                                         'the topology, coordinates, index file and QM input.')
    prepare.add_argument('-c', '--gro', required=True, help='input coordinates (.gro, .pdb, ...)')
    prepare.add_argument('-p', '--top', required=True, help='input Gromacs topology (.top)')
    prepare.add_argument('-e', '--extend', action='append', default=[], metavar='MASK',
                         help='Amber mask grown along bonds to the breakable bonds; repeatable')
    prepare.add_argument('-s', '--select', action='append', default=[], metavar='MASK',
                         help='Amber mask taken as it is, without growing; repeatable')
    prepare.add_argument('--solvate', type=float, metavar='R',
                         help='also take every solvent residue within R Angstrom of the '
                              'selection made so far')

    out = prepare.add_argument_group('output files')
    out.add_argument('-d', '--outdir', default='.', help='directory for the default names')
    out.add_argument('--out-gro', help='output coordinates (default: <outdir>/qm.gro)')
    out.add_argument('--out-top', help='output topology (default: <outdir>/qm.top)')
    out.add_argument('--out-ndx', help='output index file (default: <outdir>/qm.ndx)')
    out.add_argument('--hsd', nargs='?', const='', metavar='FILE',
                     help='also write a DFTB+ input (default: <outdir>/dftb_in.hsd)')

    setup = prepare.add_argument_group('QM/MM setup')
    setup.add_argument('--charge', type=int,
                       help='total charge of the QM region (default: derived from the force field)')
    setup.add_argument('--charge-rounding', type=_rounding, default=DEFAULT_CHARGE_ROUNDING,
                       metavar='T',
                       help='fractional part at which the derived charge rounds up in '
                            'magnitude: 0.25 by default (1.43 -> 2, 0.04 -> 0), 0.5 or '
                            '"nearest" for ordinary rounding, 0 or "away" to always round up')
    setup.add_argument('--mm-retention', choices=('no', 'classic', 'amber'), default='no',
                       help='strip the MM bonded terms of the QM region from the topology '
                            '(default: no, grompp does it itself)')
    setup.add_argument('--redistr-scheme', default='no',
                       choices=('no', 'amber', 'RC', 'RCD', 'CS'),
                       help='what happens to the charge of the MM boundary atoms')
    setup.add_argument('--no-link-bond', action='store_true',
                       help='do not connect the link atoms to their MM1 atom (funct 5)')
    setup.add_argument('--link-la-to-mm2', action='store_true',
                       help='also connect each link atom to the MM2 atoms behind its MM1 '
                            'atom (funct 5), widening the exclusion shell')
    setup.add_argument('--preset', action='append', choices=sorted(BOND_PRESETS), default=[],
                       help='breakable-bond table to use; repeatable, default protein')
    setup.add_argument('--breakable', action='append', default=[], metavar='QM:MM[:DIST]',
                       help='extra breakable bond, e.g. CG:CB:1.09; repeatable')
    setup.add_argument('--no-breakable', action='append', default=[], metavar='QM:MM',
                       help='remove a breakable bond from the table; repeatable')
    setup.add_argument('--strict-h-dist', action='store_true',
                       help='fail instead of guessing a missing link-atom distance')
    _add_hsd_options(prepare, with_method_default=True)

    # ---- rewrite-hsd -----------------------------------------------------
    rewrite = sub.add_parser('rewrite-hsd', help='update an existing dftb_in.hsd', parents=[common],
                             description='Replace parts of a DFTB+ input in place. Everything '
                                         'that is not addressed stays byte for byte, so '
                                         'hand-tuned settings survive.')
    rewrite.add_argument('hsd', help='file to update')
    rewrite.add_argument('-c', '--gro', metavar='FILE',
                         help='coordinates to take the new geometry from; a .pdb is read '
                              'just as well as a .gro')
    rewrite.add_argument('-n', '--ndx', help='index file holding the QM group')
    rewrite.add_argument('--out-gro', metavar='FILE',
                         help='also write the whole input system as a .gro, every atom in '
                              'the input order. A .pdb input is converted under its own '
                              'name by default; --no-out-gro switches that off')
    rewrite.add_argument('--no-out-gro', action='store_true',
                         help='do not write a .gro even when the input is a .pdb')
    rewrite.add_argument('--group', default='QM', help='index group to use (default: QM)')
    rewrite.add_argument('--source', help='read from this file and write to HSD instead')
    rewrite.add_argument('--charge', type=int, help='new total charge of the QM region')
    rewrite.add_argument('--no-keep-types', action='store_true',
                         help='do not inherit TypeNames from the file being rewritten')
    _add_hsd_options(rewrite)

    # ---- check -----------------------------------------------------------
    check = sub.add_parser('check', help='verify that .gro, .ndx and .hsd agree', parents=[common],
                           description='Compare the [ QM ] group of an index file with a '
                                       'coordinate file and a DFTB+ input.')
    check.add_argument('-c', '--gro', required=True)
    check.add_argument('-n', '--ndx', required=True)
    check.add_argument('--hsd', help='DFTB+ input to compare against')
    check.add_argument('--group', default='QM')

    sub.add_parser('methods', help='list the available QM methods', parents=[common])

    tables = sub.add_parser('tables', help='show the link-atom and residue-name tables', parents=[common])
    tables.add_argument('what', nargs='?', default='all',
                        choices=('all', 'bonds', 'residues', 'elements'))
    return parser


# ------------------------------------------------------------------ helpers
def _parse_bond(text, with_distance):
    fields = text.split(':')
    if len(fields) < 2 or not all(fields[:2]):
        raise QMMMError(f'--breakable expects QM:MM[:DIST], got {text!r}')
    pair = (fields[0], fields[1])
    if not with_distance:
        return pair, None
    if len(fields) == 3:
        try:
            return pair, float(fields[2])
        except ValueError:
            raise QMMMError(f'{fields[2]!r} is not a distance in Angstrom') from None
    return pair, None


def _apply_bond_tables(qm, args):
    presets = args.preset or ['protein']
    bonds, distances = set(), {}
    for name in presets:
        preset_bonds, preset_dist = BOND_PRESETS[name]
        bonds |= set(preset_bonds)
        distances.update(preset_dist)
    for text in args.breakable:
        pair, distance = _parse_bond(text, with_distance=True)
        bonds.add(pair)
        if distance is not None:
            distances[pair] = distance
    for text in args.no_breakable:
        pair, _ = _parse_bond(text, with_distance=False)
        bonds.discard(pair)
        distances.pop(pair, None)
    qm.breakable_bonds = bonds
    qm.H_dist = distances
    qm.strict_h_dist = args.strict_h_dist
    LOGGER.info('breakable bonds: %s', ', '.join(f'{a}-{b}' for a, b in sorted(bonds)) or 'none')


def _hsd_kwargs(args):
    """The hsd options the user actually gave, so defaults are not forced on a rewrite."""
    mapping = {'skpath': args.skpath, 'sk_suffix': getattr(args, 'sk_suffix', None),
               'scc_tolerance': args.scc_tolerance, 'mixer': args.mixer,
               'max_scc_iterations': args.max_scc_iterations,
               'dftbplus_version': args.dftbplus_version}
    return {key: value for key, value in mapping.items() if value is not None}


# ----------------------------------------------------------------- commands
def cmd_prepare(args):
    outdir = Path(args.outdir)
    gro = args.out_gro or outdir / 'qm.gro'
    top = args.out_top or outdir / 'qm.top'
    ndx = args.out_ndx or outdir / 'qm.ndx'

    if not args.extend and not args.select:
        raise QMMMError('nothing selected: give at least one --extend or --select mask')

    qm = QM(args.gro, args.top, gro, top, ndx)
    _apply_bond_tables(qm, args)
    for mask in args.extend:
        qm.choose_qm_to_extend(mask)
    for mask in args.select:
        qm.choose_qm_manually(mask)
    if args.solvate:
        solvent = ':' + ','.join(sorted(qm.solvent_and_ions & data.WATER_RESIDUES))
        qm.choose_qm_manually(f'(({qm.qm_input_mask})<:{args.solvate})&({solvent})')

    qm.job(qm_aim_charge=args.charge, mm_retention=args.mm_retention,
           redistr_scheme=args.redistr_scheme, link_la_to_mm1=not args.no_link_bond,
           link_la_to_mm2=args.link_la_to_mm2, charge_rounding=args.charge_rounding)

    hsd = None
    if args.hsd is not None:
        hsd = Path(args.hsd) if args.hsd else outdir / 'dftb_in.hsd'
        kwargs = _hsd_kwargs(args)
        kwargs.setdefault('sk_suffix', args.sk_suffix)
        qm.make_hsd(hsd, method=args.method, **kwargs)
    qm.check_consistency(hsd)

    print(f'\nQM region : {len(qm.qm_idx)} atoms + {len(qm.la_idx)} link atoms'
          f'{f" + {len(qm.cp_idx)} charge points" if qm.cp_idx else ""}')
    print(f'QM charge : {qm.aim_qm_charge:+d}   '
          f'(force field: {qm.qm_charge:+.4f})')
    print(f'written   : {top}, {gro}, {ndx}' + (f', {hsd}' if hsd else ''))
    print(f'\nremember to set QMcharge = {qm.aim_qm_charge} in the .mdp')
    return 0


def cmd_rewrite_hsd(args):
    if bool(args.gro) != bool(args.ndx):
        raise QMMMError('--gro and --ndx go together: both are needed for a new geometry')
    geometry = None
    if args.gro:
        write_gro_to = False if args.no_out_gro else args.out_gro
        geometry = read_qm_geometry(args.gro, args.ndx, group=args.group,
                                    write_gro_to=write_gro_to)
    elif args.out_gro:
        raise QMMMError('--out-gro needs an input to convert: give -c/--gro as well')
    kwargs = _hsd_kwargs(args)
    kwargs.pop('sk_suffix', None) if args.method is None else None
    rewrite_hsd(args.hsd, geometry=geometry, source_hsd=args.source,
                keep_types=not args.no_keep_types, method=args.method,
                charge=args.charge, **kwargs)
    if geometry is not None and geometry.gro_path is not None:
        print(f'coordinates written to {geometry.gro_path}')
    return 0


def cmd_check(args):
    groups = read_index_file(args.ndx)
    if args.group not in groups:
        raise QMMMError(f'{args.ndx} has no [ {args.group} ] group; it has: '
                        + ', '.join(groups))
    geometry = read_qm_geometry(args.gro, args.ndx, group=args.group, warn_unknown=False)
    print(f'[ {args.group} ] holds {len(geometry)} atoms of {args.gro}')
    if geometry.elements is not None:
        counts = {}
        for element in geometry.elements:
            counts[element] = counts.get(element, 0) + 1
        print('  elements guessed from the names:',
              ', '.join(f'{n}x{e}' for e, n in sorted(counts.items())))

    if not args.hsd:
        return 0
    types, species = _read_geometry_types(HsdFile.read(args.hsd))
    if not species:
        raise QMMMError(f'{args.hsd}: no usable Geometry block found')
    if len(species) != len(geometry):
        raise QMMMError(f'{args.hsd} holds {len(species)} atoms but [ {args.group} ] has '
                        f'{len(geometry)} -- the .hsd belongs to a different selection')
    print(f'{args.hsd} holds {len(species)} atoms, TypeNames = {" ".join(types or [])}')
    if geometry.elements is None:
        print('  the elements could not be guessed from the names, '
              'so only the atom count was checked')
        return 0
    mismatch = [(i, types[s - 1], geometry.elements[i]) for i, s in enumerate(species)
                if types[s - 1] != geometry.elements[i]]
    if mismatch:
        for index, in_hsd, guessed in mismatch[:10]:
            label = geometry.labels[index] if geometry.labels else f'#{index + 1}'
            print(f'  atom {index + 1} ({label}): {in_hsd} in the .hsd, {guessed} from the name')
        raise QMMMError(f'{len(mismatch)} of {len(species)} atoms disagree')
    print('OK: the index file, the coordinates and the DFTB+ input agree on every QM atom')
    return 0


def cmd_tables(args):
    if args.what in ('all', 'bonds'):
        print('breakable bonds and link-atom distances (QM atom -> MM atom, Angstrom)')
        for name, (bonds, distances) in sorted(BOND_PRESETS.items()):
            print(f'  [{name}]')
            for pair in sorted(bonds):
                print(f'     {pair[0]:>5s} -> {pair[1]:<5s} {distances.get(pair, float("nan")):.2f}')
        print('  fall-back by element:',
              ', '.join(f'{e}-H {d}' for e, d in sorted(data.H_DIST_BY_ELEMENT.items())))
    if args.what in ('all', 'residues'):
        print('\nresidue names')
        print(f'  water      ({len(data.WATER_RESIDUES):3d}):', ' '.join(sorted(data.WATER_RESIDUES)))
        print(f'  ions       ({len(data.ION_RESIDUES):3d}):', ' '.join(sorted(data.ION_RESIDUES)))
        print(f'  amino acids({len(data.AMINO_ACIDS):3d}) and nucleotides '
              f'({len(data.NUCLEIC_ACIDS):3d}) may accept redistributed charge')
    if args.what in ('all', 'elements'):
        print('\nDFTB parameters per element')
        for element in sorted(data.MAX_ANGULAR_MOMENTUM):
            hubbard = data.HUBBARD_DERIVS.get(element)
            print(f'  {element:<3s} max. angular momentum {data.MAX_ANGULAR_MOMENTUM[element]}'
                  + (f', Hubbard derivative {hubbard}' if hubbard is not None else
                     ', no Hubbard derivative (DFTB3 unavailable)'))
    return 0


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.version:
        from . import __version__
        print(f'QMMMtools {__version__}')
        return 0
    if not args.command:
        parser.print_help()
        return 1
    verbose, quiet = getattr(args, 'verbose', False), getattr(args, 'quiet', False)
    set_log_level('DEBUG' if verbose else 'WARNING' if quiet else 'INFO')

    handlers = {'prepare': cmd_prepare, 'rewrite-hsd': cmd_rewrite_hsd, 'check': cmd_check,
                'tables': cmd_tables}
    try:
        if args.command == 'methods':
            list_methods()
            return 0
        return handlers[args.command](args)
    except QMMMError as error:
        print(f'error: {error}', file=sys.stderr)
        return 2
    except KeyboardInterrupt:
        return 130


if __name__ == '__main__':
    sys.exit(main())
