"""QMMMtools -- QM/MM topology preparation for GROMACS + DFTB+/xTB.

Carves a QM region out of a finished classical GROMACS system and writes the
topology, coordinates, index file and DFTB+ input that a link-atom QM/MM run
needs, all of them describing the same atoms in the same order.

    import QMMMtools

    qm = QMMMtools.QM('npt.gro', 'topol.top', 'qm.gro', 'qm.top', 'qm.ndx')
    qm.choose_qm_to_extend('@21313,21391,21408')
    qm.job()
    qm.make_hsd('dftb_in.hsd', method='dftb3-d4', skpath='./3ob-3-1/')
    qm.check_consistency('dftb_in.hsd')

Everything chemistry-specific -- which bonds may be cut, link-atom distances,
residue names, QM methods -- lives in :mod:`QMMMtools.data` and can be replaced
per object at run time.  There is also a command line front end, ``qmmmtools``.
"""

from . import data
from .core import (QM, QMGeometry, QMMMError, HsdFile, set_log_level, LOGGER,
                   read_qm_geometry, read_index_file, write_hsd, rewrite_hsd,
                   hamiltonian_block, list_methods, get_method,
                   write_gro, gro_box_line, convert_to_gro,
                   analysis_block, forces_keyword, parse_dftbplus_version,
                   ANALYSIS_BLOCK, OPTIONS_BLOCK, DEFAULT_DFTBPLUS_VERSION)
from .data import QMMethod, QM_METHODS, DEFAULT_QM_METHOD

__version__ = '1.2.0'

__all__ = [
    'QM', 'QMGeometry', 'QMMMError', 'HsdFile', 'QMMethod',
    'set_log_level', 'LOGGER',
    'read_qm_geometry', 'read_index_file',
    'write_hsd', 'rewrite_hsd', 'hamiltonian_block',
    'write_gro', 'gro_box_line', 'convert_to_gro',
    'list_methods', 'get_method',
    'analysis_block', 'forces_keyword', 'parse_dftbplus_version',
    'ANALYSIS_BLOCK', 'OPTIONS_BLOCK', 'DEFAULT_DFTBPLUS_VERSION',
    'QM_METHODS', 'DEFAULT_QM_METHOD',
    'data', '__version__',
]
