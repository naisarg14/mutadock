---
name: mutadock-python-api
description: Use MUTADOCK from Python instead of the CLI — generate_csv, resolve_matrix/load_matrix, clean_pdb, pose_residue_map, structure_warnings, generate_pdb, calc_ddg, prepare_receptor/prepare_ligand/dock_vina and read_config. Use when embedding MUTADOCK in a script, notebook or pipeline, when custom filtering or a non-standard loop is needed, or when deciding between importing the library and shelling out to the md_* commands.
---

# MUTADOCK Python API

Shell out to the `md_*` commands for standard workflows — they handle resume,
reporting and output naming for you. Import the library when you need custom
filtering, a non-standard loop, or integration into an existing pipeline.

## Mutation enumeration and matrices

`generate_csv` takes `.pdb` or `.cif`. A CIF is converted to PDB first (the
underlying parser is PDB-only), and the converted path is what the default
output names and the CSV's `pdb` column refer to — so the result feeds the ΔΔG
calculators unchanged. Converting yourself with `convert_cif_pdb` is still fine;
it is just no longer required.

```python
from mutadock.mutation.csv_generator import generate_csv
from mutadock.mutation.helpers import convert_cif_pdb, resolve_matrix, load_matrix

# Returns (positive_csv, all_csv); defaults are derived from the input name
out_op, out_all = generate_csv("data/4QJR.cif")       # writes data/4QJR.pdb too
generate_csv("data/4QJR.pdb", matrix="BLOSUM62")      # downloaded if absent
generate_csv("data/4QJR.pdb", out_all="all.csv", out_op="positive.csv")

pdb = convert_cif_pdb("data/4QJR.cif")                # explicit, if you want the path

scores = resolve_matrix("PAM250")      # bundled; works installed or from source
scores = resolve_matrix("PAM30")       # fetched from NCBI, cached
scores = load_matrix("/path/to/custom_matrix")
```

CSV columns: `sr, pdb, chain, position, wtAA, prAA, wtProb, prProb`.

## Preparing and inspecting a structure

```python
from mutadock.mutation.helpers import (
    clean_pdb, convert_cif_pdb, fetch_pdb, pose_residue_map, structure_warnings,
)

path = fetch_pdb("4QJR")                      # download from RCSB
pdb = convert_cif_pdb(path) if path.endswith(".cif") else path

for w in structure_warnings(pdb):             # multi-model, chains, altlocs, …
    print("WARNING:", w)

clean_pdb(pdb, "protein_clean.pdb")           # strips HETATM, one model/conformer
```

`pose_residue_map` takes a **PyRosetta pose**, not a path, and returns
`{chain: [residue numbers]}` in PDB author numbering — the numbering a
`CHAIN:POSITION` mutation spec uses:

```python
from pyrosetta import init, pose_from_pdb

init("-mute all")
numbering = pose_residue_map(pose_from_pdb("protein_clean.pdb"))
print(sorted(numbering), min(numbering["A"]), max(numbering["A"]))
```

`clean_pdb` removes HETATM records, so cofactors, metals, waters and ligands are
absent from everything scored afterwards. Every ΔΔG call expects the **cleaned**
structure.

## ΔΔG scoring

```python
from mutadock.mutation.ddg_calc import calc_ddg

out = calc_ddg(
    pdb_file="protein_clean.pdb",
    in_file="mutations.csv",
    out_file="ddG.csv",
    replicates=3,
    pack_radius=8.0,
    backbone_minimization=True,   # 'min'
    cartesian=False,              # True + backbone_minimization → 'cartesian'
)
```

Protocol mapping (`resolve_ddg_params` in `mutation/helpers.py`):

| Protocol | `cartesian` | `backbone_minimization` | Use |
|---|---|---|---|
| `fast` | `False` | `False` | screening only — never report |
| `min` *(default)* | `False` | `True` | reliable general use |
| `cartesian` | `True` | `True` | research-grade; slowest |

`fast` runs a single repack with no minimization. Because ΔΔG uses a wild-type
self-reference, a failed WT repack can rank a destabilizing mutation as the top
stabilizer — never quote those numbers.

Output columns: `sr, pdb, chain, position, wtAA, prAA, wtProb, prProb,
ddG_value, ddG_sd, ddG_kcal, n_replicates, ddG_protocol`. **`ddG_value` is in
Rosetta Energy Units, not kcal/mol**; `ddG_kcal` applies a `0.34` factor
(ref2015 `cartesian_ddg` convention, Park et al. 2016). Lower is more
stabilizing; a null WT→WT mutation scores ≈ 0.

Rosetta packing is stochastic — use `replicates=3` or more so `ddG_sd` exists,
or the ranking has no error bar.

For a site-focused study, filter the enumeration CSV and pass it straight to
`calc_ddg`; there is no distance or binding-site parameter.

## Mutant structures

```python
from mutadock.mutation.generate_mutant_pdb import generate_pdb

muts = [
    {"chain": "A", "position": 386, "wtAA": "ASN", "prAA": "ALA"},
    {"chain": "B", "position": 45,  "wtAA": "GLY", "prAA": "VAL"},
]
generate_pdb("protein.pdb", muts[:1])                              # one mutant
generate_pdb("protein.pdb", muts, output_folder="./mutants/")      # independent mutants
generate_pdb("protein.pdb", muts, output_folder="./mutants/", compound=True)  # one pose
```

`wtAA` is optional and used only for naming; 1- or 3-letter codes are accepted.
Files are written as `{stem}_{WT3}-{CHAIN}{POS}-{NEW1}.pdb`, with the new residue
as a **one-letter** code — e.g. `4QJR_clean_ASN-A386-A.pdb`. `md_mutate` names its
mutants differently (`{sr}_{WT}-{CHAIN}{POS}-{NEW}.pdb`), so match on chain and
position rather than on an assumed filename.

`generate_pdb` returns the list of paths it wrote — use that instead of
reconstructing names.

## Docking

```python
from mutadock.docking.vina_helper import (
    prepare_receptor, prepare_ligand, dock_vina, read_config, fetch_ligand,
)

prepare_receptor("data/4QJR.cif", "receptor.pdbqt")   # PDB or CIF; PDBFixer + meeko
prepare_ligand("data/Ligand.sdf", "ligand.pdbqt")     # needs 3D coords + explicit H

center, box, exhaustiveness, n_poses, n_poses_write, overwrite, seed = read_config("config.txt")

dock_vina(
    receptor="receptor.pdbqt",
    ligand="ligand.pdbqt",
    output="output.pdbqt",
    log_file="vina.log",
    center=center,
    box_size=box,
    seed=seed,
)
```

`read_config` returns the tuple above and falls back to safe defaults for
missing keys — `center` defaults to the **origin**, so always set `center_x/y/z`.
The default seed is 19.

An affinity is only comparable to another computed in the same box with the same
seed. When recording results, use `add_score_to_csv(...)` with `center`,
`box_size`, `exhaustiveness` and `seed` so the provenance columns are populated.

## Reports

```python
from mutadock.report.data import discover_run, summary_stats

run = discover_run("results/")      # classifies the CSVs and mutant PDBs
print(summary_stats(run))
```

Or shell out: `md_report -d results/ --top 20`.

## When to shell out instead

The CLI adds resume checkpoints, output naming, structure warnings, report
generation and combination sizing. Reimplementing those in a script is the usual
source of inconsistent results — call `md_mutate` / `md_dock` via `subprocess`
and import the library only for the parts you are genuinely customizing.
