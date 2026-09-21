# QMMMtools

QM/MM topology preparation with link atoms for **GROMACS + DFTB+/xTB**.
Carves a QM region out of a finished classical GROMACS system and writes four files
that are guaranteed to describe the same atoms in the same order:

| file | content |
|---|---|
| `qm.top` | QM–QM bonds turned into connections (`funct 5`), link atoms and charge points as `[ virtual_sites2 ]`, adjusted charges, corrected `[ molecules ]` |
| `qm.gro` | coordinates **in the topology's atom order**, velocities preserved |
| `qm.ndx` | `[ QM ]` (QM atoms + link atoms), `[ Biomolecule ]`, `[ Water_and_ions ]`, `[ System ]` |
| `dftb_in.hsd` | DFTB+/xTB input with the QM atoms **in the same order as `[ QM ]`** |

NOTE: all ligands should be introduced directly into forcefield as .rtps (NOT .itps) !!!

We recommend to use our python module **[ligtools](https://github.com/Vetrov-Anton/ligtools)** for these purposes. 

Works with proteins, nucleic acids, lipids and arbitrary ligands — switching system
type means replacing two tables, nothing else. Usable both as a Python library and as
a command line tool.

---

## Contents

- [Installation](#installation)
- [Quick start](#quick-start) · [command line](#command-line) · [Python](#python)
- [Choosing the QM region](#choosing-the-qm-region)
- [Charge](#charge)
- [The QM/MM boundary](#the-qmmm-boundary)
- [MM bonded terms](#mm-bonded-terms-inside-the-qm-region)
- [Running GROMACS](#running-gromacs)
- [DFTB+ / xTB input](#dftb--xtb-input)
- [Adapting the tables](#adapting-qmmmtoolsdata)
- [Troubleshooting](#troubleshooting)
- [Limitations](#limitations)

---

## Installation

Everything installs straight from GitHub.

### For the command line tool — `pipx`

[pipx](https://pipx.pypa.io) puts the `qmmmtools` command on your `PATH` in its own
isolated environment, so nothing is added to your system or conda Python:

```bash
pipx install git+https://github.com/Vetrov-Anton/QMMMtools.git
```

```bash
qmmmtools --help
qmmmtools methods
```

Upgrading and removing:

```bash
pipx upgrade QMMMtools
pipx uninstall QMMMtools
```

### For `import QMMMtools` — `pip`

**pipx deliberately hides the package from your interpreter**, so a pipx install does
*not* make `import QMMMtools` work in a script or a notebook. Use `pip` for that, into the
environment you actually work in:

```bash
pip install git+https://github.com/Vetrov-Anton/QMMMtools.git
```

Installing both ways is perfectly reasonable: pipx for the command, pip for the library.

### From a clone

```bash
git clone https://github.com/Vetrov-Anton/QMMMtools.git
cd QMMMtools
pip install -e .        # editable: your edits to QMMMtools/data.py take effect at once
pipx install .          # or the command line tool from the same clone
```

### Requirements

* Python ≥ 3.8, [ParmEd](https://github.com/ParmEd/ParmEd) ≥ 4.0, NumPy — installed automatically
* a GROMACS built with the DFTB+ QM/MM interface (Kubař *et al.*), to run the result
* Slater–Koster files (e.g. `3ob-3-1`) for the DFTB methods; the xTB methods need none

---

## Quick start

### Command line

```bash
qmmmtools prepare \
    -c npt.gro -p topol.top \
    -e '@21313,21391,21408' \
    --solvate 3.0 \
    --hsd --skpath ./3ob-3-1-ophyd/ --mixer anderson
```

```
QM region : 202 atoms + 9 link atoms
QM charge : -2   (force field: -1.4453)
written   : qm.top, qm.gro, qm.ndx, dftb_in.hsd

remember to set QMcharge = -2 in the .mdp
```

| subcommand | what it does |
|---|---|
| `prepare` | build the whole QM/MM system |
| `rewrite-hsd` | update an existing `dftb_in.hsd` — needs only a `.gro` and a `.ndx` |
| `check` | verify that a `.gro`, a `.ndx` and a `.hsd` describe the same QM atoms |
| `methods` | list the available QM methods |
| `tables` | print the link-atom, residue-name, element and Slater–Koster tables |

`qmmmtools prepare --help` lists every option.

### Python

```python
import QMMMtools

qm = QMMMtools.QM('npt.gro', 'topol.top', 'qm.gro', 'qm.top', 'qm.ndx')

qm.choose_qm_to_extend('@21313,21391,21408')      # grow each seed to the breakable bonds
qm.job()                                          # charge derived from the force field
qm.make_hsd('dftb_in.hsd', method='dftb3-d3h5',
            skpath='./3ob-3-1-ophyd/', mixer='anderson')
qm.check_consistency('dftb_in.hsd')               # .ndx / .gro / .hsd must agree
```

`examples.ipynb` walks through everything in detail.

---

## Choosing the QM region

Two selectors, freely combinable and repeatable. Both take
[Amber masks](https://parmed.github.io/ParmEd/html/amber.html#amber-selection-masks)
in the numbering of the **input** system.

```python
qm.choose_qm_to_extend('@21313')          # depth-first growth, stops at breakable bonds
qm.choose_qm_manually('@1455,4052')       # exactly these atoms, no growth
```

`qm.qm_input_mask` always holds everything selected so far, so later masks can build on it:

```python
qm.choose_qm_manually(f'(({qm.qm_input_mask})<:3.0)&(:SOL)')   # + every water within 3 Å
```

On the command line the same thing is `-e`/`--extend`, `-s`/`--select` (both repeatable)
and `--solvate R`.

### Ready-made masks

`QMMMtools.data` carries the masks that come up in almost every setup, so they do not have
to be typed out:

```python
from QMMMtools import data

qm.choose_qm_manually(f'({data.PROTEIN_SIDECHAIN_MASK})&(:HIS)')       # side chains of every His
qm.choose_qm_manually(f'({data.PROTEIN_SIDECHAIN_MASK})&(:ARG,LYS)')   # of Arg and Lys
qm.choose_qm_to_extend(f'({data.PROTEIN_BACKBONE_MASK})&(:167-169)')   # a stretch of backbone
```

| mask | selects |
|---|---|
| `PROTEIN_SIDECHAIN_MASK` | side chains of the amino acids |
| `NOT_PROTEIN_SIDECHAIN_MASK` | its negation: backbone, solvent, ions, ligands |
| `PROTEIN_BACKBONE_MASK` | backbone of the amino acids |
| `PROTEIN_BACKBONE_ATOM_MASK` | just the names, `@C,CA,H,...`, unrestricted by residue |
| `NUCLEIC_BASE_MASK`, `NUCLEIC_BACKBONE_MASK`, `NOT_NUCLEIC_BASE_MASK` | the same split for nucleotides: base against sugar-phosphate |
| `AMINO_ACID_MASK`, `NUCLEIC_ACID_MASK` | `:ALA,ARG,…` and `:DA,DC,…` |

Only the **backbone** is enumerated by name — twenty names covering the Gromacs
(`HA1/HA2`, `OC1/OC2`), Amber (`HA2/HA3`, `OXT`) and CHARMM (`HN`, `OT1/OT2`) spellings —
and the side chain is its negation, so the list never has to grow with the residue types.
The side-chain masks are restricted to amino-acid residues on purpose: a bare
`!@N,CA,C,O,…` would take in every water molecule and every ligand.

Build your own the same way:

```python
data.atom_mask({'CA', 'CB'})        # '@CA,CB'
data.residue_mask({'ALA', 'GLY'})   # ':ALA,GLY'
```

Growth stops at the *directed* pairs in `qm.breakable_bonds`. `('CB', 'CA')` means
"walking from an atom named CB onto one named CA is a cut". The protein table lists the
CA–CB bond in both directions, so a side-chain selection stops at CB and a backbone
selection stops at CA; the peptide cuts `N–CA`, `CA–C` and `C–CA` let a single peptide
unit be taken on its own.

### No MM atom between two QM atoms

Whatever the masks ask for, the region is closed so that **no MM atom is bonded to more
than one QM atom**. Such an atom — the `…QM–MM–QM…` case — cannot be capped: it would
need one link atom per QM neighbour, several massless hydrogens a bond length apart on the
same centre, each pretending to terminate the region on its own, and the atom they cap
would be described by the QM calculation and by the force field at the same time.

It happens more easily than it looks. Growing from residues 496 and 498 walks along both
peptide bonds into residue 497, stops at `N–CA` from one side and at `C–CA` from the other,
and leaves exactly one MM atom — `CA` — wedged between two QM ones:

```
... C(496) ── N(497) ── CA(497) ── C(497) ── N(498) ...
      QM        QM        MM!        QM        QM
```

Such an atom joins the QM region, and the region grows on from it to the next breakable
bond exactly as it grows from a seed — above, `CA` and its `HA` come in and the cut moves
to `CA–CB`. Taking one bridge in can expose the next, so the closure repeats until none is
left. `qm.close_qm_region()` does it, and both selectors as well as `determine_qm()` call
it, so it also catches a bridge produced by a manual mask or by the combination of several
masks. A selection that has no such atom is left exactly as it was.

## Charge

`job()` with no `qm_aim_charge` derives the integer QM charge from the force field:

```
aim = round(q_total) − round(q_total − q_QM, threshold=0.25)
```

Ordinary rounding switches over at a fractional part of 0.5. Here the threshold is
**0.25**, because a cut that runs through a charge group leaves part of a formal charge
behind on the MM side, so the MM sum comes out too small in magnitude:

| MM sum | 0.25 (default) | 0.5 (`'nearest'`) | 0.0 (`'away'`) |
|---:|---:|---:|---:|
| 1.43 | **2** | 1 | 2 |
| −1.43 | **−2** | −1 | −2 |
| 0.04 | **0** | 0 | 1 |
| 2.00 | **2** | 2 | 2 |

Any threshold between 0 and 1 works — `charge_rounding=0.4`,
`--charge-rounding 0.4` — and `'nearest'` / `'away'` are names for 0.5 and 0.0. When the
chosen threshold and ordinary rounding disagree, both candidates are logged, so an
ambiguous case never passes silently.

The difference between the chosen integer and the force-field sum is moved onto the MM
atoms so the total system charge is unchanged. Acceptors are **only** protein /
nucleic-acid atoms (`qm.redistr_residues`) — never water and never ions.

## The QM/MM boundary

Each QM–MM bond gets a hydrogen link atom, written as a two-body virtual site of function
type 2 — GROMACS's "2fd", *on a line with a fixed distance*:

```
[ virtual_sites2 ]
;   LA    QM    MM  funct      d (nm)
 10707   450   448      2      0.1090     ; qmmm
```

The site is built as `r_LA = r_QM + d · r_QM→MM / |r_QM→MM|`, so it stays on the QM–MM
axis at exactly `H_dist` from the QM atom however the QM–MM bond stretches during the run.
(Function type 1 would store a fraction of the bond instead, frozen from the starting
geometry, and the QM–link-atom distance would breathe with the bond.) The charge points of
the `RC`, `RCD` and `CS` schemes stay type 1 on purpose: they are defined as fractions of
the MM1–MM2 bond.

`H_dist` is the equilibrium X–H bond length of the **capped QM atom** — C 1.09 Å,
N 1.01 Å, O 0.97 Å, S 1.34 Å. Anything else stretches or compresses a real bond inside
the QM calculation.

Two independent flags add `funct 5` bonds around the link atom — "connections" that carry
no potential and exist only so that grompp generates exclusions from them:

| flag | effect |
|---|---|
| `link_la_to_mm1=True` (default) | one bond to the MM atom the link atom caps; with `nrexcl = 3` this already excludes the link atom from everything within three bonds of MM1 |
| `link_la_to_mm2=False` | one bond to each further MM neighbour of MM1, pushing the exclusion shell one bond further out |

On the command line: `--no-link-bond` turns the first off, `--link-la-to-mm2` turns the
second on. They can be combined in any way.

What happens to the charge of the MM boundary atom is `redistr_scheme`
(`--redistr-scheme`):

| value | effect |
|---|---|
| `'no'` (default) | left as it is |
| `'amber'` | zeroed, spread over the MM acceptors |
| `'RC'` | zeroed, `q/n` on the midpoint of each MM1–MM2 bond |
| `'RCD'` | as RC but `2q/n` on the midpoint and `−q/n` on MM2 — keeps the dipole |
| `'CS'` | charge shift: `q/n` onto MM2 plus a `±q/n` pair around it |

## MM bonded terms inside the QM region

`mm_retention` (`--mm-retention`) controls whether they are already stripped from the
written topology:

| value | effect |
|---|---|
| `'no'` (default) | left in — grompp of this build removes them itself and prints a table of what it removed |
| `'classic'` | grompp's own rule: a term goes once all but one of its atoms are QM; also drops the 1-4 pairs between a QM atom and a bonded MM atom |
| `'amber'` | only terms whose atoms are *all* QM |

`'classic'`/`'amber'` also pre-empt `GMX_QMMM_BONDED_SCHEME`.

**QM–QM bonds are always converted to `funct 5`**, independent of `mm_retention`. This is
not cosmetic: left as `funct 1` they are turned into constraints by `constraints = h-bonds`
*before* grompp removes the QM bonded terms, and the QM hydrogens end up rigid (measured:
50 degrees of freedom silently lost).

## Running GROMACS

```bash
gmx grompp -f qm.mdp -c qm.gro -p qm.top -n qm.ndx -o qm.tpr

GMX_QMMM_NREXCL=3 GMX_QMMM_VARIANT=1 gmx mdrun -deffnm qm
```

with an `.mdp` containing

```
QMMM      = yes
QMMM-grps = QM
QMmethod  = RHF        ; grompp insists on a value, this build ignores it
QMbasis   = STO-3G     ; likewise
QMcharge  = -2         ; must equal Charge in dftb_in.hsd
QMmult    = 1
```

`dftb_in.hsd` must sit in the run directory. mdrun reads it once and then overwrites the
coordinates every step, so the *numbers* in `Geometry` do not matter for the run — but the
atom **count** and **order** do.

### The index groups

| group | atoms |
|---|---|
| `QM` | what the QM code sees: the QM atoms **plus the link atoms** |
| `Biomolecule` | the whole rewritten moleculetype — the protein or nucleic acid together with everything quantum: ligands, cofactors, metals, link atoms, and any solvent taken into the QM region |
| `Water_and_ions` | the bulk that was left untouched |
| `System` | everything |

`Biomolecule` and `Water_and_ions` are the natural pair for temperature coupling and for
trajectory output — the interesting part of the system against the bath around it:

```
tc-grps           = Biomolecule Water_and_ions
tau-t             = 0.1 0.1
ref-t             = 300 300
compressed-x-grps = Biomolecule
```

`Biomolecule` is **not** the QM region — it is the whole molecule the QM region was cut out
of, usually a few thousand atoms against a couple of hundred quantum ones. Ask for `QM` when
you want the quantum atoms themselves.

A group named `freeze` with the atoms of `Biomolecule` is written at the end of the file so
that an older `.mdp` keeps working; the first three groups also keep their positions, so
selecting a group by number from a pipe still means what it did.

---

## DFTB+ / xTB input

```bash
qmmmtools methods
```

| method | notes |
|---|---|
| `dftb3-d3h5` | DFTB3/3ob + D3(zero) + H5 hydrogen-bond correction — the default |
| `dftb3-d4` | DFTB3/3ob + D4 + H damping |
| `dftb3-d3bj` | DFTB3/3ob + D3(BJ) + H damping |
| `dftb3` | DFTB3/3ob, no dispersion |
| `dftb2` | plain SCC-DFTB (mio), no third order |
| `gfn2-xtb`, `gfn1-xtb` | via tblite, no Slater–Koster files |

Every DFTB3 method uses the 3ob-3-1 set as its authors specify it: the maximum angular
momenta and the Hubbard derivatives are the fifteen-element tables from the README that
ships with the set, and the γ^h damping exponent is ζ = 4.00. The D4 parameters are the
`dftb_3ob` two-body entry of `dftd4`; the D3(BJ) ones (a1 = 0.746, a2 = 4.191, s8 = 3.209)
are those published with 3ob itself. `dftb3-d3h5` replaces the γ^h damping with the H5
correction and adds D3 with zero damping and H–H repulsion, which is the recipe the DFTB+
manual gives for DFTB3-D3H5.

GFN2-xTB was verified to run through the GROMACS interface.

### Slater–Koster file names

DFTB+ is never given a list of files; it builds every name from the two element names:

```
Prefix + first + Separator + second + Suffix
```

A parameter set usually ships the same data twice — `3ob-3-1` holds both `Mg-C.skf` and
`mgc-c.spl` — so the spelling has to match the files that are really in the directory,
and the wrong one makes DFTB+ stop with *SK file not found* in the middle of an mdrun.
`sk_format` (`--sk-format`) picks it:

| format | file | separator | suffix | lower case |
|---|---|---|---|---|
| `skf` — default | `Mg-C.skf` | `-` | `.skf` | no |
| `spl` | `mgc-c.spl` | *(none)* | `-c.spl` | yes |

```python
qm.make_hsd('dftb_in.hsd', skpath=SKPATH)                     # Mg-C.skf
qm.make_hsd('dftb_in.hsd', skpath=SKPATH, sk_format='spl')    # mgc-c.spl
```

```bash
qmmmtools prepare ... --hsd --skpath ./3ob-3-1/ --sk-format spl
qmmmtools tables slater-koster
```

`sk_separator`, `sk_suffix` and `sk_lowercase` (`--sk-separator`, `--sk-suffix`,
`--sk-lowercase` / `--no-sk-lowercase`) override one field each, for a set that spells
its files in yet another way:

```python
qm.make_hsd('dftb_in.hsd', skpath=SKPATH, sk_separator='_', sk_suffix='.par')
```

When the directory is readable from here it is checked against the choice, and a
mismatch names the convention that would have worked:

```
49 of the 49 Slater-Koster files DFTB+ will ask for are not in ./3ob-3-1/,
e.g. C-C.par, C-H.par, C-O.par; the files of the "skf" convention (Mg-C.skf)
are all there -- pass sk_format="skf"
```

An input that already works only needs its `SlaterKosterFiles` block re-pointed:

```bash
qmmmtools rewrite-hsd dftb_in.hsd --skpath /new/3ob-3-1/ --sk-format skf
```

### DFTB+ release

DFTB+ 24.1 renamed the `Analysis` switch that returns the forces from `CalculateForces` to
`PrintForces`, and the old spelling is not accepted as an alias. `dftbplus_version`
(`--dftbplus-version`) picks the dialect:

```python
qm.make_hsd('dftb_in.hsd', skpath=SKPATH)                        # 25.1 -> PrintForces
qm.make_hsd('dftb_in.hsd', skpath=SKPATH, dftbplus_version=21)   # 21.x -> CalculateForces
```

The default is 25.1, a bare year means the `.1` release of that series, and anything older
than 21 is refused. Nothing else this module writes differs between 21.x and 25.x, so an
existing input moves between the two by rewriting that one keyword:

```bash
qmmmtools rewrite-hsd dftb_in.hsd --dftbplus-version 21
```

### Updating an existing `dftb_in.hsd`

`rewrite_hsd` edits an input in place — everything it is not asked to change stays byte
for byte, so hand-tuned settings survive:

```python
qm.rewrite_hsd('dftb_in.hsd')                                   # coordinates only
qm.rewrite_hsd('dftb_in.hsd', charge=-3, mixer='anderson')      # + settings
qm.rewrite_hsd('dftb_in.hsd', method='gfn2-xtb', charge=-2)     # swap the Hamiltonian
qm.rewrite_hsd('dftb_in.hsd', skpath='/new/3ob/', sk_format='spl')   # other SK files
qm.rewrite_hsd('dftb_in.hsd', blocks={'Driver': 'Driver = {}'}) # any top-level block
qm.rewrite_hsd('dftb_in.hsd', geometry=False, scc_tolerance='1e-7')  # settings only
```

### …without a topology

Only a coordinate file and an index file are needed, so a `.hsd` can be refreshed even
for files somebody else produced. Any format ParmEd reads is accepted — a `.pdb` just as
well as a `.gro`:

```bash
qmmmtools rewrite-hsd dftb_in.hsd -c qm.gro  -n qm.ndx --charge -2 --mixer anderson
qmmmtools rewrite-hsd dftb_in.hsd -c sys.pdb -n qm.ndx --charge -2
```

```python
geometry = QMMMtools.read_qm_geometry('qm.gro', 'qm.ndx', group='QM')
QMMMtools.rewrite_hsd('dftb_in.hsd', geometry=geometry, skpath='/new/path/')
```

A `.pdb` input is also written out as a `.gro` — the **whole** system, every atom in the
input's own order, so the index file and the topology stay valid across the conversion.
By default it lands next to the input under the same name (`sys.pdb` → `sys.gro`):

| | command line | Python |
|---|---|---|
| default for a `.pdb` | — | `write_gro_to=None` |
| pick the name | `--out-gro other.gro` | `write_gro_to='other.gro'` |
| do not write one | `--no-out-gro` | `write_gro_to=False` |

A `.gro` input is left alone, since there is nothing to convert. The conversion is also
available on its own:

```python
QMMMtools.convert_to_gro('sys.pdb')                 # -> sys.gro
QMMMtools.convert_to_gro('sys.pdb', 'other.gro')
```

Coordinates survive the round trip to the 0.001 nm a `.gro` stores, and the box comes from
the `CRYST1` record; a file without one gets a bounding box and a warning.

Elements are taken from the `TypeNames` of the file being rewritten whenever the atom
count matches; otherwise they are guessed from the atom and residue names (`CA` in `ALA`
is a carbon, `CA` in residue `CA` is calcium, `SE` in `MSE` is selenium). A disagreement
between the two sources is reported — that check has already caught a stale `.hsd` in
practice. Note that a `.gro` stores 0.001 nm, so a rewrite from one is accurate to 0.01 Å;
irrelevant for a QM/MM start, since mdrun overwrites the coordinates anyway.

---

## Adapting `QMMMtools.data`

Everything system-specific lives in `QMMMtools/data.py`, and every table can be replaced
per object at run time instead of editing the file.

```python
import QMMMtools
from QMMMtools import data

# --- nucleic acids -------------------------------------------------------
qm.breakable_bonds = set(data.NUCLEIC_BREAKABLE_BONDS)   # N9-C1', C5'-O5', C3'-O3'
qm.H_dist          = dict(data.NUCLEIC_H_DIST)

# --- one extra cut -------------------------------------------------------
qm.breakable_bonds.add(('CG', 'CB'))
qm.H_dist[('CG', 'CB')] = 1.09

# --- remove a cut you do not want ---------------------------------------
qm.breakable_bonds.discard(('C', 'CA'))
qm.H_dist.pop(('C', 'CA'), None)

# --- protein + DNA complex: merge both tables ---------------------------
qm.breakable_bonds = set(data.PROTEIN_BREAKABLE_BONDS) | set(data.NUCLEIC_BREAKABLE_BONDS)
qm.H_dist          = {**data.PROTEIN_H_DIST, **data.NUCLEIC_H_DIST}

# --- who may accept redistributed charge ---------------------------------
qm.redistr_residues = set(data.NUCLEIC_ACIDS)
qm.redistr_residues.discard('PRO')

# --- a water model with an unusual residue name --------------------------
qm.solvent_and_ions.add('T4P')

# --- a new element, from a set other than 3ob ----------------------------
data.MAX_ANGULAR_MOMENTUM_OTHER['Se'] = 'd'      # merged into MAX_ANGULAR_MOMENTUM
data.MAX_ANGULAR_MOMENTUM['Se'] = 'd'
data.HUBBARD_DERIVS['Se'] = -0.11                # only needed by the DFTB3 methods

# --- a new QM method -----------------------------------------------------
data.QM_METHODS['gfn0-xtb'] = data.QMMethod(
    'gfn0-xtb', 'xtb', 'GFN0-xTB via tblite', xtb_method='GFN0-xTB')
data.QM_METHODS.pop('dftb2', None)          # or remove one

# --- a set that names its files differently ------------------------------
data.SK_FORMATS['par'] = data.SKFormat('par', '_', '.par', False)
```

The same is reachable from the command line with `--preset`, `--breakable QM:MM:DIST` and
`--no-breakable QM:MM`.

Missing entries fail loudly rather than silently: an unknown link-atom distance falls back
to the element default *with a warning* (or raises, with `qm.strict_h_dist = True`), and a
missing `MaxAngularMomentum` / `HubbardDerivs` raises with the name of the element to add.

| table | purpose |
|---|---|
| `WATER_RESIDUES`, `ION_RESIDUES` | never accept redistributed charge |
| `AMINO_ACIDS`, `NUCLEIC_ACIDS`, `POLYMER_RESIDUES` | default charge acceptors |
| `PROTEIN_SIDECHAIN_MASK`, `PROTEIN_BACKBONE_MASK`, `NUCLEIC_BASE_MASK`, … | ready-made selection masks |
| `PROTEIN_BACKBONE_ATOMS`, `NUCLEIC_BACKBONE_ATOMS`, `atom_mask`, `residue_mask` | what those masks are built from |
| `*_BREAKABLE_BONDS`, `*_H_DIST` | where the QM region may be cut, and the link-atom distance |
| `H_DIST_BY_ELEMENT` | fall-back X–H bond lengths |
| `ELEMENTS`, `TYPE2ELEMENT`, `SPECIAL_SITE_ELEMENTS`, `guess_element` | element determination |
| `MAX_ANGULAR_MOMENTUM_3OB`, `HUBBARD_DERIVS` | the fifteen elements of 3ob-3-1, straight from its README |
| `MAX_ANGULAR_MOMENTUM_OTHER` | elements from other Slater–Koster sets, merged into `MAX_ANGULAR_MOMENTUM` |
| `QM_METHODS` | the `Hamiltonian` blocks |
| `SK_FORMATS`, `DEFAULT_SK_FORMAT`, `sk_file_name` | how the Slater–Koster file names are spelled |

## Consistency

```bash
qmmmtools check -c qm.gro -n qm.ndx --hsd dftb_in.hsd
```

```python
qm.check_consistency('dftb_in.hsd')
```

re-reads the written files and verifies that `[ QM ]`, the `.gro` and the `.hsd` describe
the same atoms in the same order. Cheap, and worth running before every production job —
a silent reordering between these three files is the failure mode that is hardest to
notice afterwards.

## Logging

```python
QMMMtools.set_log_level('DEBUG')     # every cut bond, every include rewritten
QMMMtools.set_log_level('WARNING')   # quiet
```

Command line: `-v` / `-q`.

---

## Troubleshooting

**`SCC is NOT converged` in `mdrun`.** Two independent causes, both seen on the same system:

1. *The mixer.* Charged metal sites (Mg²⁺ with a pyrophosphate, for instance) do not
   converge with the DFTB+ default Broyden mixer no matter how many iterations you allow.
   Use `mixer='anderson'` / `--mixer anderson`, which writes a small-step Anderson mixer —
   measured: no convergence in 500 iterations before, converged in 62 after.
2. *The charge.* If preparation printed *"the QM charge is ambiguous"*, the chosen
   threshold and ordinary rounding disagreed, and the wrong integer can keep the SCC from
   converging at all. Measured on the example system: the force-field sum was −1.4453;
   ordinary rounding gives −1, which never converged, the 0.25 threshold gives −2, which
   converged in 157 iterations. Try the other candidate — one command, no re-preparation:

   ```bash
   qmmmtools rewrite-hsd dftb_in.hsd --charge -1   # and QMcharge = -1 in the .mdp
   ```

**`No default Bond types` from grompp.** A ligand whose bonded parameters are written
inline (GAFF/ACPYPE) lost them. This module only strips parameters from the QM–QM bonds it
converts to `funct 5`; if you see this, check that nothing else edited the topology.

**`cannot tell the element of atom …`.** ParmEd found no atomic number, i.e. the `at.num`
column of `[ atomtypes ]` is missing for that type. Fix the force field, or add the type
to `data.TYPE2ELEMENT`.

**DFTB+ halts on `dftb_in.hsd` without saying much.** Usually a dialect mismatch: a 21.x–23.x
binary reading `PrintForces`, or a 24.1+ binary reading `CalculateForces`. Point
`--dftbplus-version` at the release you actually run.

**`Could not open SK-file` / *SK file not found*.** The parameter directory spells its
files the other way — `mgc-c.spl` where `Mg-C.skf` was asked for, or the reverse. Switch
with `--sk-format` (`skf` is the default), no re-preparation needed:

```bash
qmmmtools rewrite-hsd dftb_in.hsd --sk-format spl
```

Preparation warns about this by itself whenever `skpath` is readable from the machine that
writes the input.

**`no DFTB angular momentum known for [...]` or `no Hubbard derivative known for [...]`.**
The element is outside 3ob-3-1, which covers only Br, C, Ca, Cl, F, H, I, K, Mg, N, Na, O,
P, S and Zn. Use a Slater–Koster set that has it and add the element to
`data.MAX_ANGULAR_MOMENTUM` / `data.HUBBARD_DERIVS` with that set's own values.

**`MM atom(s) are bonded to more than one QM atom`.** Only reachable by driving the
low-level methods by hand: the QM region was handed to `find_qmmm_bonds()` without being
closed, and one MM atom would carry a link atom per QM neighbour. Pass the selection
through `qm.close_qm_region()`; `choose_qm_to_extend()`, `choose_qm_manually()` and
`determine_qm()` all do it already, so `job()` cannot hit this.

**`[ molecules ] accounts for N atoms but the structure has M`.** The topology and the
coordinate file do not describe the same system.

**`import QMMMtools` fails after `pipx install`.** Expected — pipx isolates the package
from your interpreter. Use `pip install git+https://github.com/Vetrov-Anton/QMMMtools.git` for library use.

---

## Limitations

* The input must be a standalone GROMACS `.top` (with `[ molecules ]`), not a bare `.itp`.
* Molecules that contain a QM atom, or that may accept charge, are merged into a single
  moleculetype; a molecule cannot be split between the QM block and the rest.
* The link-atom placement assumes a single bond is cut per QM–MM pair. One link atom per
  MM atom is guaranteed — an MM atom bonded to two QM atoms is taken into the QM region
  instead — but cutting a double bond or an aromatic ring is chemically unsound and is not
  checked for.
* The `[ molecules ]` bookkeeping assumes the coordinate file follows the topology order,
  which is the GROMACS convention.

## License

MIT.
