---
name: mutadock-results
description: Interpret MUTADOCK output — the ΔΔG and docking CSV schemas, REU vs kcal/mol, the sign convention, which mutations are genuinely stabilizing, computing ΔΔG_binding from docking affinities, rebuilding HTML/PPTX reports with md_report, and the caveats (HETATM stripping, single conformer, screening-only runs) that limit what the numbers support.
---

# Reading MUTADOCK results

## ΔΔG CSVs

Columns: `sr, pdb, chain, position, wtAA, prAA, wtProb, prProb, ddG_value,
ddG_sd, ddG_kcal, n_replicates, ddG_protocol`.

| Column | Meaning |
|---|---|
| `wtProb` / `prProb` | Substitution-matrix scores for the wild-type and proposed residue |
| `ddG_value` | **ΔΔG in Rosetta Energy Units (REU), not kcal/mol** |
| `ddG_sd` | Standard deviation across replicates (blank when `n_replicates` is 1) |
| `ddG_kcal` | `ddG_value × 0.34` unless `--reu-to-kcal` overrode the factor |
| `n_replicates` | Independent estimates averaged |
| `ddG_protocol` | `fast`, `min` or `cartesian` — **check this before quoting anything** |

**Sign convention: lower is more stabilizing.** ΔΔG is
`score(mutant) − score(wild type run through the identical protocol)`, so a null
WT→WT mutation scores ≈ 0 and the values are not biased toward "stabilizing".
`*_sorted.csv` files are sorted lowest → highest, so the most stabilizing
candidates are at the top.

**Units.** REU ≈ but ≠ kcal/mol. The `0.34` factor (≈1/2.94) is the ref2015
`cartesian_ddg` convention from Park et al. 2016
(doi:10.1021/acs.jctc.6b00819). Quote REU, or quote kcal/mol *and* the factor.

**Before trusting a ranking, check `ddG_protocol`:**

- `fast` — **screening only.** A single repack with no minimization can leave
  clashes it cannot relieve, and because of the wild-type self-reference a
  failed WT repack can put a destabilizing mutation at the top of the
  stabilizer list. Reports from these runs carry a "SCREENING ONLY — values
  unreliable" banner. Re-run the shortlist with `min` or `cartesian` before
  drawing any conclusion.
- `min` — the reliable default.
- `cartesian` — research-grade; pair with `--replicates 3`+ so `ddG_sd` exists.

**Without `ddG_sd`, a ranking has no error bar.** Rosetta packing is stochastic;
differences smaller than the replicate spread are not real. Re-run with
`--replicates` before claiming one mutation beats another.

Sort any CSV by column name (`-n`) or column number (`-N`), ascending (`-a`,
the useful direction for ΔΔG) or descending (`-d`):

```bash
md_csv_sort -i ddG.csv -n ddG_value -a          # most stabilizing first
md_csv_sort -i ddG.csv -n affinity -a -o by_affinity.csv
```

## Docking CSV

`docking_results.csv` columns: `sr, name, affinity, center_x, center_y,
center_z, size_x, size_y, size_z, exhaustiveness, seed`. `affinity` is the best
pose's Vina score in kcal/mol (more negative = tighter predicted binding).
`name` is derived from the pose file with `_out` and the extension stripped,
so `1_ASN-A386-HIS_Ligand_out.sdf` → `1_ASN-A386-HIS_Ligand`. Mutants produced by
`md_mutate` carry a serial-number prefix and those from `md_generate_pdb` carry
the structure stem, so parse the chain and position out of the name rather than
assuming one layout.

**The provenance columns are not decoration.** An affinity is only interpretable
relative to the box it was computed in, and at fixed settings a different seed
can move a single affinity by ~2 kcal/mol by finding a different binding mode.
Two rows with different `center_*`/`size_*`/`seed` values **are not comparable**,
and their difference is not a meaningful Δ.

## Computing ΔΔG_binding

MUTADOCK does not emit this; compute it yourself. A difference is only
meaningful between two rows that docked **the same ligand** into **the same box
with the same seed and exhaustiveness**:

> ΔΔG_binding = affinity(mutant) − affinity(wild type)

Negative ⇒ the mutant is predicted to bind more tightly. Dock the wild-type
receptor in the same `md_dock` run (add it to `receptors.txt`) so the comparison
is matched by construction.

**Identify the wild type by its exact `name`, and say which ligand you are
comparing.** Substring tests do not work here: mutants are usually derived from
the cleaned structure, so the wild-type stem `4QJR_clean` is a *prefix of every
mutant row* (`4QJR_clean_ASN-A386-A_Ligand`), and `"clean" in name` matches them
all. `name` is `{receptor}_{ligand}`, so give the baseline row's full name and
the ligand suffix, and run once per ligand.

```bash
python compare_affinity.py results/docking_results.csv --wt 4QJR_clean_Ligand --ligand Ligand
```

```python
# compare_affinity.py
import argparse
import csv
import sys

PROVENANCE = ("center_x", "center_y", "center_z",
              "size_x", "size_y", "size_z", "exhaustiveness", "seed")

ap = argparse.ArgumentParser()
ap.add_argument("csv_path")
ap.add_argument("--wt", required=True, help="exact 'name' of the wild-type row")
ap.add_argument("--ligand", required=True, help="ligand suffix shared by comparable rows")
args = ap.parse_args()

rows = list(csv.DictReader(open(args.csv_path)))


def settings(row):
    """Box/seed/exhaustiveness the row was scored with, or None if unrecorded."""
    values = [row.get(col) for col in PROVENANCE]
    # A legacy 3-column CSV, or one truncated to match an older file, cannot be
    # compared at all: missing provenance is not a match.
    return None if any(v in (None, "") for v in values) else tuple(values)


baselines = [r for r in rows if r["name"] == args.wt]
if len(baselines) != 1:
    sys.exit(f"expected exactly one row named '{args.wt}', found {len(baselines)}")
wt = baselines[0]
if settings(wt) is None:
    sys.exit(f"wild-type row '{args.wt}' has no recorded box/seed; cannot compare")

suffix = "_" + args.ligand
compared = 0
for row in rows:
    name = row["name"]
    if name == args.wt or not name.endswith(suffix):
        continue          # different ligand, or the baseline itself
    if settings(row) is None:
        print(f"{name}: SKIPPED - missing box/seed provenance")
    elif settings(row) != settings(wt):
        print(f"{name}: SKIPPED - different box/seed/exhaustiveness from the wild type")
    else:
        print(f"{name}: ddG_binding = {float(row['affinity']) - float(wt['affinity']):+.2f} kcal/mol")
        compared += 1

if not compared:
    print(f"no comparable rows for ligand '{args.ligand}'", file=sys.stderr)
```

The script refuses rather than guesses: an unrecorded box or seed, a ligand it
cannot match, or a duplicated baseline all produce a skip, because a number
computed across mismatched settings is wrong rather than approximate.

A Vina score is a docking score, not a measured affinity. Treat the ranking as a
hypothesis to test, and do not quote absolute values as binding free energies.

## Reports

`md_mutate` and `md_dock` build a `reports/` bundle automatically unless
`--no-report` is passed:

```
results/reports/
├── figures/{ddg_distribution,top_mutations,docking_affinity,top_mutant_structure}.png
├── report.html     # self-contained; includes an interactive NGL 3D viewer
└── report.pptx
```

Rebuild from any finished run directory:

```bash
md_report -d results/
md_report -d results/ --html-only
md_report -d results/ --pptx-only
md_report -d results/ -o custom_dir/ --top 20 --title "4QJR imatinib scan"
```

| Flag | Meaning |
|---|---|
| `-d DIR` | Run directory to read |
| `-o DIR` | Write reports elsewhere (default `<run>/reports/`) |
| `--top N` | Number of top mutations to tabulate/chart |
| `--title` | Report title |
| `--html-only` / `--pptx-only` | Build just one format |
| `-q` | Quiet |

Using the **same** `-o/--output-dir` for `md_mutate` and then `md_dock` makes one
combined report with both ΔΔG and affinities.

## Caveats that limit the conclusions

- **HETATM records are stripped** by `clean_pdb`. ΔΔG is computed on the apo
  protein — cofactors, metals, waters and the ligand itself are absent. A
  mutation that coordinates a metal is scored without it.
- **One conformer, one model.** Alternate locations and NMR ensembles are
  reduced before scoring; a structure that warned about altlocs was scored on a
  single arbitrary conformer.
- **Every chain is mutated and scored.** On a homodimer, a position appears once
  per chain and the hit list is duplicated.
- **ΔΔG is folding stability, not binding.** It says nothing about the ligand.
  Only the docking step addresses binding, and only as a docking score.
- **The two halves use different structures.** ΔΔG is scored on the cleaned apo
  PDB; docking uses a PDBFixer-repaired receptor converted to PDBQT. They are
  not the same coordinates.
- **Matrix scores are a prior, not evidence.** `prProb` only decided which
  substitutions got enumerated.

## Which mutant to take forward

Require all of: `ddG_protocol` is `min` or `cartesian`; `ddG_value` is
meaningfully below zero relative to `ddG_sd`; the position is one you can
rationalize structurally; and — if binding matters — a ΔΔG_binding computed
against a wild type docked in the same box with the same seed.
