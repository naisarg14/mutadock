# MUTADOCK

[![License: GPL v3](https://img.shields.io/badge/license-GPL--3.0-blue.svg)](LICENSE)
[![PyPI version](https://img.shields.io/pypi/v/mutadock.svg)](https://pypi.org/project/mutadock/)
[![Python](https://img.shields.io/pypi/pyversions/mutadock.svg)](https://pypi.org/project/mutadock/)
[![Docs](https://readthedocs.org/projects/mutadock/badge/?version=latest)](https://mutadock.readthedocs.io/en/latest/)
[![Coverage](https://img.shields.io/badge/coverage-45.45%25-yellow.svg)]()

## Introduction

MUTADOCK is a comprehensive library for protein mutation studies and multi-receptor/multi-ligand docking. It integrates automated protein mutation via PyRosetta with an AutoDock Vina docking pipeline, enabling systematic exploration of how mutations affect receptor–ligand binding affinity.

### Key Features

**Automated Protein Mutation**
- Single-, double-, and triple-point mutation prediction using PyRosetta
- ΔΔG scoring to rank mutations by stability change
- Flexible substitution matrix support: built-in PAM250 and BLOSUM62, any matrix from the [NCBI BLAST FTP](https://ftp.ncbi.nih.gov/blast/matrices/) (downloaded automatically), or a custom file

**Docking Pipeline**
- Batch docking of N receptors × M ligands with AutoDock Vina
- Automatic receptor preparation via PDBFixer + meeko (`mk_receptor`)
- Ligand preparation via meeko + RDKit
- Accepts both `.pdb` and `.cif` receptor files

**Usability**
- Simple CLI for each workflow step
- Python API for scripting and integration into existing pipelines


## System Requirements

| Requirement | Minimum |
|-------------|---------|
| Python | 3.11 or 3.12 |
| RAM | 8 GB (16 GB recommended for PyRosetta) |
| Disk | ~6 GB (PyRosetta installation) |
| OS | Linux, macOS, Windows |

AutoDock Vina must be installed separately — see [Vina installation](https://autodock-vina.readthedocs.io/en/latest/installation.html).


## Installation

```bash
pip install mutadock
```

Install PyRosetta (required for ΔΔG calculations):

```bash
md_install_dependencies
```

Install receptor preparation dependencies (conda recommended for pdbfixer/openmm):

```bash
conda install -c conda-forge pdbfixer openmm
```


## Quick Start

The `data/` directory in this repository contains a sample receptor (`4QJR.cif`) and ligand (`Ligand.sdf`) you can use to try the workflow immediately.

### The 30-second demo: `md_quick`

One command takes a **PDB ID**, a **mutation**, and a **ligand code**, and gives you a ΔΔG plus a docked pose — fetching the structure from RCSB and the ligand from PubChem automatically:

```bash
md_quick --pdb-id 4QJR --mutation A:386:ASN:HIS --ligand-code imatinib
```

```
  Mutation : ASN-A386-HIS
  ΔΔG      : -27.30 REU (negative = stabilizing)
  Affinity : -5.66 kcal/mol
  Report   : mdquick_4QJR/reports/report.html
```

It fetches + cleans the structure, computes the ΔΔG for that single mutation, builds the mutant, fetches the ligand, finds the pocket with AutoSite (or pass `-c config.txt`), docks, and writes the standard `reports/` bundle. Use `-i protein.pdb` for a local structure, `--ligand-file lig.sdf` for a local ligand, or `--ligand-code 2244` / `cid:5291` / `name:aspirin` for PubChem. Requires the `autosite` binary (ADFRsuite) on PATH when no config is given.

The staged tools below give you full control over each step.

### 1. Generate mutation candidates

```bash
md_csv_generator -i data/4QJR.cif -o mutations.csv -O mutations_all.csv
```

Use a different substitution matrix (downloaded automatically if not present locally):

```bash
md_csv_generator -i data/4QJR.cif --matrix BLOSUM62
```

Use a custom matrix file:

```bash
md_csv_generator -i data/4QJR.cif --matrix-file /path/to/my_matrix.txt
```

### 2. Run the full mutation pipeline

```bash
md_mutate -i data/4QJR.cif
```

This generates all output files described in the [Mutation Output](#mutation-output) table and produces `4QJR_modified_mutants.txt` listing every mutated PDB for use with `md_dock`.

### 3. Dock receptors against a ligand

Create a ligand list:

```bash
echo "data/Ligand.sdf" > ligands.txt
```

Then dock against all mutants (or any receptor list):

```bash
md_dock -r 4QJR_modified_mutants.txt -l ligands.txt -c config.txt
```


## How-To Guide

### Mutation Studies

`md_mutate` takes a PDB or CIF file and runs the complete mutation and ΔΔG pipeline.

```bash
md_mutate -i protein.pdb
md_mutate --pdb-id 4QJR              # fetch the structure from RCSB instead of -i
md_mutate -i protein.pdb -o results/ # write all outputs to results/
md_mutate -h                         # all options
```

Provide exactly one of `-i/--input` (a local `.pdb`/`.cif` file) or `--pdb-id` (a
4-character RCSB accession, downloaded automatically). `md_mutate` also warns when
the input contains multiple models, alternate conformations (altlocs), or
insertion codes, since these affect residue numbering.

By default every output file is written next to the input structure. Pass
`-o/--output-dir DIR` to collect them in `DIR` instead (created if absent); the
paths recorded in `*_mutants.txt` point into `DIR`, so they remain valid input
for `md_dock`.

#### Mutation Output

| # | File | Description |
|---|------|-------------|
| 1 | `protein_modified_mutations_all.csv` | All possible single-residue substitutions |
| 2 | `protein_modified_mutations.csv` | Substitutions with positive matrix score |
| 3 | `protein_modified_mutations_ddG.csv` | ΔΔG for each mutation in file 2 |
| 4 | `protein_modified_mutations_ddG_sorted.csv` | File 3 sorted lowest→highest ΔΔG |
| 5 | `protein_modified_double_ddg.csv` | Double-mutation ΔΔG combinations |
| 6 | `protein_modified_double_ddg_sorted.csv` | File 5 sorted |
| 7 | `protein_modified_triple_ddg.csv` | Triple-mutation ΔΔG combinations |
| 8 | `protein_modified_triple_ddg_sorted.csv` | File 7 sorted |
| 9 | `protein_modified_mutants.txt` | List of mutated PDB paths (direct input for `md_dock`) |

#### ΔΔG rigor and units (read this before trusting the numbers)

ΔΔG is computed as `score(mutant) − score(wild-type self-mutation reference)`, where **both** sides run the identical repack(+minimization) protocol — so a null WT→WT mutation scores ≈ 0 and the values are not biased toward "stabilizing."

- **Units:** `ddG_value` is in **Rosetta Energy Units (REU), not kcal/mol.** A `ddG_kcal` column is also written, and reports show both. REU ≈ but ≠ kcal/mol.
- **Scaling factor:** the REU→kcal/mol factor is **`REU_TO_KCAL_SCALE` in `src/mutadock/mutation/predict_ddG.py`** (default `0.34`, ≈ 1/2.94; the ref2015 `cartesian_ddg` convention, Park et al. 2016). **To change it, edit that constant** or pass `--reu-to-kcal FACTOR` per run.
- **Protocol (`--protocol`, default `min`):** the recorded protocol is written to a `ddG_protocol` column, and reports flag screening-only runs.

  ```bash
  md_mutate -i protein.pdb                       # default: 'min' = repack + backbone/side-chain minimization (reliable)
  md_mutate -i protein.pdb --protocol cartesian  # ref2015_cart + Cartesian minimization (most accurate, cartesian_ddg style; slowest)
  md_mutate -i protein.pdb --protocol fast        # single repack, NO minimization — SCREENING ONLY (see warning below)
  md_mutate -i protein.pdb --replicates 3         # mean ± SD over 3 stochastic repacks (adds ddG_sd, n_replicates)
  md_mutate -i protein.pdb --pack-radius 10 --reu-to-kcal 0.29
  ```

  > **Why `min` is the default, not `fast`:** a single repack with **no minimization** can leave clashes it can't relieve — and because ΔΔG uses a wild-type self-reference, a failed WT repack can make a *destabilizing* mutation rank as the top *stabilizer* (observed: a site flipping from +69 REU to −375 REU). Minimization relieves those clashes and fixes it. `--protocol fast` (~1 min/proteome vs ~12 min for `min`, ~30 min for `cartesian`) is retained for **screening only** and its reports carry a "SCREENING ONLY — values unreliable" banner.
  >
  > For research-grade numbers prefer `--protocol cartesian` and `--replicates 3`+. The cartesian_ddg protocol is cited in Park, H. *et al.* (2016) *J. Chem. Theory Comput.* **12**(12):6201–6212, doi:10.1021/acs.jctc.6b00819.

### Generating Mutant PDB Files

Use `md_generate_pdb` to apply one or more specific mutations to a PDB and write the mutant structure(s) directly — no ΔΔG scoring required.

**Single mutation on the command line:**

```bash
md_generate_pdb -i protein.pdb -m A:386:ASN:ALA
# omit the wild-type AA if unknown:
md_generate_pdb -i protein.pdb -m A:386:ALA
```

**From a CSV file (one independent PDB per row):**

```bash
md_generate_pdb -i protein.pdb -c mutations.csv
```

CSV format:

```
chain,position,wtAA,prAA
A,386,ASN,ALA
B,45,GLY,VAL
```

`wtAA` is optional (used for output file naming). Column names are case-insensitive and common aliases (`old_aa`, `new_aa`, `from`, `to`) are accepted.

**Compound mutant** (all CSV mutations applied to a single pose):

```bash
md_generate_pdb -i protein.pdb -c mutations.csv --compound
```

**Custom output folder:**

```bash
md_generate_pdb -i protein.pdb -c mutations.csv -o ./mutants/
```

Output files are named `{stem}_{WTAA}-{CHAIN}{POS}-{NEWAA}.pdb` (e.g. `protein_ASN-A386-A.pdb`).

### Docking Studies

```bash
md_dock -r receptors.txt -l ligands.txt -c config.txt
md_dock -r receptors.txt -l ligands.txt -c config.txt -o results/  # collect outputs in results/
md_dock -h   # all options
```

Every receptor in `receptors.txt` is docked against every ligand in `ligands.txt`. Receptors can be `.pdb` or `.cif` — they are fixed and converted to PDBQT automatically.

By default docking outputs go to an `out/` folder next to each receptor. Pass
`-o/--output-dir DIR` to collect poses, logs, `docking_results.csv`, and the
`*_completed.txt` resume file in a single `DIR` instead. Prepared PDBQT files and
AutoSite caches still live next to their inputs so they can be reused across runs.

#### Docking Output

| # | Output | Description |
|---|--------|-------------|
| 1 | PDBQT files | Prepared receptor and ligand files |
| 2 | Log file | Vina output with binding scores per combination |
| 3 | Output PDB | Top 5 docking poses per combination |
| 4 | Output PDBQT | Best pose (Vina split) per combination |
| 5 | Output SDF | Best pose as SDF for visualization |
| 6 | `docking_results.csv` | All affinities tabulated for easy analysis |

### Reports

Every `md_mutate` and `md_dock` run automatically produces a **self-contained
report** so you don't have to open the raw CSVs by hand:

- `report.html` — a single file (opens offline / emails cleanly) with a ΔΔG
  distribution, a top-stabilizing-mutations table + chart, a docking-affinity
  chart + table, and an **interactive 3D viewer of the top mutant pose** (NGL,
  inlined).
- `report.pptx` — the same content as presentation slides.
- `figures/` — each chart also saved as a standalone PNG, ready to drop into
  your own slides or manuscript.

Everything is collected in a **`reports/` folder** inside the run's output
directory:

```
results/reports/
├── figures/
│   ├── ddg_distribution.png
│   ├── top_mutations.png
│   ├── docking_affinity.png
│   └── top_mutant_structure.png
├── report.html
└── report.pptx
```

Add `--no-report` to skip report generation. Rebuild (or build from an older
run) at any time with `md_report`:

```bash
md_report -d results/                 # writes results/reports/{report.html,report.pptx,figures/}
md_report -d results/ --html-only     # HTML only
md_report -d results/ -o custom_dir/ --top 20   # write reports into custom_dir/ instead
md_report -h
```

Using the **same** `-o/--output-dir` for `md_mutate` and then `md_dock` yields a
single combined report showing both ΔΔG and docking affinities.

Reporting needs `matplotlib`, `python-pptx`, and `jinja2` (installed
automatically with `mutadock`).


## CLI Reference

| Command | Description |
|---------|-------------|
| `md_quick` | One-shot demo: PDB ID + mutation + ligand code → ΔΔG + docked pose + report |
| `md_mutate` | Full mutation + ΔΔG pipeline from a PDB/CIF file |
| `md_dock` | Batch receptor–ligand docking |
| `md_vina_dock` | Direct AutoDock Vina CLI wrapper |
| `md_csv_generator` | Generate all possible substitutions with matrix scoring |
| `md_csv_sort` | Sort any CSV by column name or number |
| `md_ddg_single` | Calculate single-mutation ΔΔG from a mutations CSV |
| `md_ddg_double` | Calculate double-mutation ΔΔG combinations |
| `md_ddg_triple` | Calculate triple-mutation ΔΔG combinations |
| `md_generate_pdb` | Generate mutant PDB file(s) from a single mutation or a CSV list |
| `md_report` | Build a self-contained HTML + PPTX report from a run directory |


## Python API

```python
# --- Mutation ---
from mutadock.mutation.csv_generator import generate_csv
from mutadock.mutation.helpers import resolve_matrix, load_matrix

# Generate mutation CSV with default matrix (PAM250)
generate_csv("data/4QJR.cif")

# Use BLOSUM62 (downloaded automatically if absent)
generate_csv("data/4QJR.cif", matrix="BLOSUM62")

# Load a matrix directly
score_dict = load_matrix("data/PAM250")   # dict[str, dict[str, int]], 3-letter keys
score_dict = resolve_matrix("PAM30")      # downloads PAM30 from NCBI if needed

# --- Generate mutant PDB ---
from mutadock.mutation.generate_mutant_pdb import generate_pdb

# Single mutation
outputs = generate_pdb("protein.pdb", [{"chain": "A", "position": 386, "wtAA": "ASN", "prAA": "ALA"}])

# Multiple independent mutants from a list
mutations = [
    {"chain": "A", "position": 386, "wtAA": "ASN", "prAA": "ALA"},
    {"chain": "B", "position": 45,  "wtAA": "GLY", "prAA": "VAL"},
]
outputs = generate_pdb("protein.pdb", mutations, output_folder="./mutants/")

# One compound mutant (all mutations on the same pose)
outputs = generate_pdb("protein.pdb", mutations, output_folder="./mutants/", compound=True)

# --- Docking ---
from mutadock.docking.vina_helper import prepare_receptor, prepare_ligand, dock_vina

prepare_receptor("data/4QJR.cif", "receptor.pdbqt")        # PDB or CIF
prepare_ligand("data/Ligand.sdf", "ligand.pdbqt")

dock_vina(
    receptor="receptor.pdbqt",
    ligand="ligand.pdbqt",
    output="output.pdbqt",
    log_file="vina.log",
    center=[10.0, 5.0, 20.0],
    box_size=[20.0, 20.0, 20.0],
)
```


## Troubleshooting

**`mk_receptor: command not found`**
meeko is not installed or not on PATH.
```bash
pip install meeko
```

**`pdbfixer` / `openmm` import error during receptor preparation**
```bash
conda install -c conda-forge pdbfixer openmm
# or
pip install pdbfixer openmm
```

**PyRosetta fails to install**
Run the bundled installer which handles license and platform detection:
```bash
md_install_dependencies
```

**`vina: command not found`**
Vina is not bundled with mutadock. Install it from the [official guide](https://autodock-vina.readthedocs.io/en/latest/installation.html).

**Matrix download fails**
If NCBI FTP is unreachable, download the matrix manually and use `--matrix-file`:
```bash
md_csv_generator -i protein.pdb --matrix-file /path/to/PAM30
```

**CIF file not recognized**
PDBFixer and BioPython both support `.cif` natively. Make sure the file extension is `.cif` or `.pdb` — other extensions are not accepted.


## Applications

- **Protein Engineering:** Identify stabilising mutations for therapeutic proteins
- **Drug Discovery:** Screen mutant variants for changes in binding affinity
- **Biochemical Research:** Study how point mutations alter protein–ligand interactions


## Contributing

Bug reports and pull requests are welcome at [github.com/naisarg14/mutadock](https://github.com/naisarg14/mutadock/issues). Please open an issue before submitting a large PR so we can discuss the approach.


## Documentation

Full API reference: [mutadock.readthedocs.io](https://mutadock.readthedocs.io/en/latest/)


## License

[GNU General Public License v3.0](LICENSE) — © 2026 Naisarg Patel


## Contact

[naisarg.patel14@hotmail.com](mailto:naisarg.patel14@hotmail.com)
