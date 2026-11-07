---
name: mutadock-docking
description: Run MUTADOCK batch docking with AutoDock Vina — every receptor in a list against every ligand, with receptor/ligand PDBQT preparation handled automatically. Use for "dock these mutants against this ligand", writing or debugging a Vina config.txt, defining or auto-detecting the search box with AutoSite, making a docking run reproducible with a seed, resuming an interrupted run, or fixing Vina/meeko timeouts.
---

# Batch docking with `md_dock`

Every receptor in `receptors.txt` is docked against every ligand in
`ligands.txt`. Receptors may be `.pdb` or `.cif`; they are repaired with
PDBFixer and converted to PDBQT automatically.

```bash
md_dock -r receptors.txt -l ligands.txt -c config.txt -o results/
```

## Inputs

**`receptors.txt`** — one structure path per line. `md_mutate` writes exactly
this file as `*_mutants.txt`:

```
results/protein_ASN-A386-HIS.pdb
results/protein_ASN-A386-GLN.pdb
```

**`ligands.txt`** — one ligand path per line (`.sdf` recommended):

```
data/Ligand.sdf
```

The ligand must already have **3D coordinates and explicit hydrogens** — this is
validated up front and is the most common input rejection. See
`mutadock-debugging`.

**`config.txt`** — Vina search parameters as `key = value` lines; `#` starts a
comment. Missing keys fall back to the defaults shown:

```
center_x = 10.0
center_y = 5.0
center_z = 20.0
size_x = 30.0
size_y = 30.0
size_z = 30.0
exhaustiveness = 32
n_poses = 20
n_poses_write = 5
overwrite = True
seed = 19
```

| Key | Default |
|---|---|
| `center_x/y/z` | `0.0` — **always set these**; the default box is at the origin |
| `size_x/y/z` | `30.0` |
| `exhaustiveness` | `32` |
| `n_poses` | `20` |
| `n_poses_write` | `5` |
| `overwrite` | `True` |
| `seed` | `19` |

## Options

| Flag | Meaning |
|---|---|
| `-r/--receptor_txt FILE` | Text file of receptor paths |
| `-l/--ligand_txt FILE` | Text file of ligand paths |
| `-c/--config FILE` | Vina config file |
| `-a/--autosite PDB` | A pocket PDB produced by AutoSite, used to define the box |
| `-o/--output-dir DIR` | Collect poses, logs, `docking_results.csv` and the resume file in DIR |
| `-q/--quiet` | Quiet mode |
| `-i/--ignore_existing` | Re-dock pairs already marked complete |
| `--seed N` | Vina RNG seed for **every** docking in the run |
| `--no-report` | Skip the HTML + PPTX report |

Without `-o`, outputs go to an `out/` folder next to each receptor. Prepared
PDBQT files and AutoSite caches always stay next to their inputs so they can be
reused across runs.

## Defining the box

Three ways, in order of preference:

1. **Explicit `-c config.txt`** — reproducible, comparable, no extra dependency.
2. **`-a pocket.pdb`** — a pocket PDB from a previous AutoSite run.
3. **Automatic AutoSite** — needs the `autosite` binary (ADFRsuite) on PATH.
   It creates `{stem}_autosite_out/` next to the receptor and **refuses to
   reuse a directory that already exists**, so stale output must be removed
   before re-running rather than being silently reused with old grids.

A Vina affinity is only interpretable relative to its box. **Two receptors
scored in different boxes are not comparable, and the difference between them is
not a ΔΔG_binding.** When comparing mutants to wild type, dock them all with the
same `config.txt` and the same seed.

## Reproducibility

Seed precedence: **CLI `--seed` > `seed =` in the config > the built-in default
of 19.** At fixed settings a seed change can move a single affinity by about
2 kcal/mol by finding a different binding mode, so a result without its seed
cannot be replicated. `docking_results.csv` records the box, exhaustiveness and
seed for every row for exactly this reason.

Change the seed deliberately to produce an independent replicate of a whole
pipeline; keep it fixed when comparing variants.

## Resume

Completed receptor–ligand pairs are listed in `*_completed.txt` and skipped
automatically, so an interrupted run resumes when re-issued. Pass
`-i/--ignore_existing` to force everything to dock again.

## Outputs

| Output | Description |
|---|---|
| `*.pdbqt` | Prepared receptor and ligand files |
| `*_log.txt` | Vina output with per-combination scores |
| `*_out.pdbqt` | Raw multi-pose Vina output |
| `*_out.sdf` | Best pose, extracted for visualization |
| `docking_results.csv` | Affinities plus provenance |
| `*_completed.txt` | Resume checkpoint |

`docking_results.csv` columns: `sr, name, affinity, center_x, center_y,
center_z, size_x, size_y, size_z, exhaustiveness, seed`. The column order is
append-only, so a CSV written by an older version stays readable by name; if an
existing file has fewer columns the new row is truncated to match and a warning
is logged rather than producing a ragged file. Start a fresh CSV when you want
full provenance.

Interpretation, and how to turn these into ΔΔG_binding, is in
`mutadock-results`.

## Combined report

Using the **same** `-o/--output-dir` for `md_mutate` and then `md_dock` yields a
single report showing both ΔΔG and docking affinities.

## Troubleshooting

**`mk_receptor: command not found`** — meeko is missing: `pip install meeko`.

**`vina: command not found`** — Vina is not bundled; install it separately
(`mutadock-setup`).

**A command times out.** Defaults: 900 s receptor preparation, 1800 s AutoSite,
3600 s Vina. Override with `MUTADOCK_RECEPTOR_PREP_TIMEOUT`,
`MUTADOCK_AUTOSITE_TIMEOUT`, `MUTADOCK_VINA_TIMEOUT`; `0` disables a timeout.
Large receptors and high `exhaustiveness` are the usual causes.

**`FileExistsError` on `{stem}_autosite_out`** — remove the stale directory, or
pass `-c`/`-a` to skip AutoSite.

**Ligand rejected for missing 3D coordinates or hydrogens** — see
`mutadock-debugging`.

**All affinities look identical or implausible** — check `center_x/y/z` in the
config. An unset centre puts a 30 Å box at the origin, far from the protein.
