---
name: mutadock-mutation
description: Plan and run MUTADOCK mutation studies — enumerate substitutions with a scoring matrix, score ΔΔG with PyRosetta for single, double and triple mutants, and write mutant PDB files. Use for "which mutations stabilize this protein?", "scan residues near a site", "compute ΔΔG for these substitutions", "make a mutant PDB", choosing between the fast/min/cartesian protocols, or sizing a run that would otherwise take hours.
---

# Mutation studies and ΔΔG scanning

Pick the route before running anything — the wrong one costs hours.

| Question | Route |
|---|---|
| ΔΔG for a handful of named substitutions | `md_csv_generator` → filter → `md_ddg_single` |
| Stabilizing mutations near one site/chain | same filtered route, filtered by position |
| Unbiased scan of the whole protein | `md_mutate` — **pilot first** |
| Just the mutant structure, no ΔΔG | `md_generate_pdb` |
| One mutation, end to end with docking | `md_quick` (see `mutadock-quickstart`) |

## The cost model — read before choosing `md_mutate`

A default `md_mutate` run takes **hours**. The flags do not all control runtime:

- The **single-mutation pass** scores *every* positive-score substitution from
  `md_csv_generator`. **No flag bounds it.** On a few-hundred-residue protein
  that is hundreds to thousands of PyRosetta jobs.
- `-s/--double N` (default 100) takes the top-N singles and scores **C(N,2)**
  pairs — 100 gives **4950** double ΔΔG calculations.
- `-t/--triple N` (default 40) takes the top-N double rows, extracts their
  unique single mutations, and scores **C(k,3)** triples.
- `-S/-D/-T` (default 15 each) only cap how many mutant **PDB files** get
  written. They do not reduce ΔΔG work at all.

So "make it smaller" means lowering `-s` and `-t`, or skipping `md_mutate`
entirely for a site-focused question.

## Route A — site-focused (preferred when the site is known)

There is no distance or binding-site flag. Filter the enumeration CSV instead:

```bash
# 1. Enumerate every substitution (fast, no PyRosetta)
md_csv_generator -i protein.pdb -o mutations.csv -O mutations_all.csv

# 2. Keep only the positions of interest
python - <<'PY'
import csv
keep = {380, 381, 382, 383, 384, 385, 386, 387}
rows = [r for r in csv.DictReader(open("mutations.csv"))
        if r["chain"] == "A" and int(r["position"]) in keep]
w = csv.DictWriter(open("site.csv", "w", newline=""), fieldnames=rows[0].keys())
w.writeheader(); w.writerows(rows)
print(len(rows), "substitutions kept")
PY

# 3. Score only those, at research-grade settings
md_ddg_single -p protein_clean.pdb -i site.csv -o site_ddg.csv \
    --protocol cartesian --replicates 3

md_csv_sort -i site_ddg.csv -n ddG_value -a   # lowest (most stabilizing) first
```

**`-p` must be the cleaned PDB.** `md_mutate` produces `*_clean.pdb` via
`clean_pdb` (HETATM stripped, one model, one conformer); PyRosetta scoring
expects that structure. Produce it first if you are skipping `md_mutate`:

```bash
python -c "from mutadock.mutation.helpers import clean_pdb; clean_pdb('protein.pdb','protein_clean.pdb')"
```

## Route B — full scan with `md_mutate`

**`md_mutate` has no pilot-sized mode.** `-s`/`-t` shrink only the combination
steps; the single-mutation pass still scores every positive-score substitution
(hundreds to thousands of PyRosetta jobs), so even `-s 10 -t 0` runs for hours.

Pilot with Route A instead — a handful of substitutions through `md_ddg_single`,
which finishes in seconds and exercises the same cleaning, numbering and scoring
code:

```bash
md_csv_generator -i protein.pdb -o mutations.csv -O mutations_all.csv
head -4 mutations.csv > pilot.csv          # header + 3 substitutions
python -c "from mutadock.mutation.helpers import clean_pdb; clean_pdb('protein.pdb','protein_clean.pdb')"
md_ddg_single -p protein_clean.pdb -i pilot.csv -o pilot_ddg.csv --protocol min
```

Check the row count `md_csv_generator` reported — that is how many single-point
ΔΔG jobs a full `md_mutate` will run — then confirm the chain and numbering in
`pilot_ddg.csv` are what you expect before committing to:

```bash
md_mutate -i protein.pdb -s 60 -t 20 -o full/ --protocol cartesian --replicates 3
```

Run it detached (or under `nohup`/`tmux`) and let its resume behaviour handle
interruptions — do not wait on it interactively.

Flags:

| Flag | Meaning |
|---|---|
| `-i FILE` / `--pdb-id ID` | Local `.pdb`/`.cif`, or an RCSB accession to fetch (exactly one) |
| `-s/--double N` | Top-N singles fed to double ΔΔG (default 100 → C(N,2) jobs) |
| `-t/--triple N` | Top-N double rows fed to triple ΔΔG (default 40) |
| `-S/-D/-T N` | Mutant PDBs written for single/double/triple (default 15 each) |
| `-o/--output-dir DIR` | Collect all outputs in DIR (created if absent) |
| `--no-append` | Fresh run — do **not** resume from existing CSVs |
| `--no-report` | Skip the HTML + PPTX report |
| `--quiet` | Suppress informational logging |

Plus the shared protocol flags below. By default outputs are written next to the
input structure; `-o` keeps them together and the paths recorded in
`*_mutants.txt` stay valid for `md_dock`.

`md_mutate` **resumes by default**: it appends to its checkpoint CSVs and skips
completed items, so an interrupted run can simply be re-issued. Use
`--no-append` when a fresh run is intended.

It also warns about multiple models, multiple chains, similar chain sequences,
altlocs and insertion codes. These warnings matter — every chain is mutated and
scored, so a homodimer doubles the work and duplicates the hits. See
`mutadock-debugging` before ignoring one.

## Protocol and units — the part that is easy to get wrong

`--protocol`, `--replicates`, `--pack-radius` and `--reu-to-kcal` are shared by
`md_quick`, `md_mutate` and all three `md_ddg_*` commands
(`add_ddg_protocol_args`, `mutation/helpers.py`).

| `--protocol` | What it does | Use for |
|---|---|---|
| `fast` | single repack, **no minimization** | screening only — never report |
| `min` *(default)* | repack + backbone/side-chain minimization | reliable general use |
| `cartesian` | ref2015_cart + Cartesian minimization | research-grade; slowest |

> **Why `min` and not `fast`:** a repack with no minimization can leave clashes
> it cannot relieve. ΔΔG is computed against a **wild-type self-reference**
> (`score(mutant) − score(WT run through the same protocol)`), so a failed WT
> repack can make a destabilizing mutation rank as the top stabilizer — observed
> flipping a site from +69 REU to −375 REU. Reports from `fast` runs carry a
> "SCREENING ONLY — values unreliable" banner.

- `--replicates N` — N independent estimates per mutation; adds `ddG_sd` and
  `n_replicates`. Rosetta packing is stochastic, so a single estimate has no
  error bar. Use ≥3 for anything reported.
- `--pack-radius Å` — neighbourhood repacked/minimized around the mutation
  (default 8.0).
- **`ddG_value` is in Rosetta Energy Units (REU), not kcal/mol.** `ddG_kcal`
  applies a `0.34` factor (≈1/2.94, the ref2015 `cartesian_ddg` convention,
  Park et al. 2016, doi:10.1021/acs.jctc.6b00819). Override per run with
  `--reu-to-kcal FACTOR`; never edit installed package code to change it.
- **Lower ΔΔG = more stabilizing.** A null WT→WT mutation scores ≈ 0.

## Enumerating substitutions — `md_csv_generator`

```bash
md_csv_generator -i protein.pdb                       # PAM250 (default)
md_csv_generator -i protein.pdb --matrix BLOSUM62     # auto-downloaded if absent
md_csv_generator -i protein.pdb --matrix-file ./my_matrix.txt
md_csv_generator -i protein.pdb -o positive.csv -O all.csv
```

`-o/--positive` keeps substitutions with a positive matrix score; `-O/--all`
writes every substitution. Columns: `sr, pdb, chain, position, wtAA, prAA,
wtProb, prProb`. PAM250 and BLOSUM62 ship with the package; any other NCBI BLAST
matrix is downloaded to `~/.cache/mutadock/matrices` (relocate with
`MUTADOCK_DATA_DIR`).

## Scoring ΔΔG directly — `md_ddg_single|double|triple`

```bash
md_ddg_single -p protein_clean.pdb -i mutations.csv -o ddG.csv --protocol min
md_ddg_double -p protein_clean.pdb -i ddG_sorted.csv --resume
md_ddg_triple -p protein_clean.pdb -i double_ddg_sorted.csv --resume
```

`-p` is the cleaned PDB, `-i` the input CSV, `-o` the output CSV (defaults to
`<input>_ddG.csv`). `--resume` keeps rows already in the output CSV and scores
only the rest — it matches on the recorded `ddG_protocol` and `n_replicates`, so
resuming with different protocol flags correctly recomputes instead of mixing
incompatible numbers.

Output columns: `sr, pdb, chain, position, wtAA, prAA, wtProb, prProb,
ddG_value, ddG_sd, ddG_kcal, n_replicates, ddG_protocol`.

## Writing mutant PDBs — `md_generate_pdb`

No ΔΔG scoring; just applies mutations and writes structures.

```bash
md_generate_pdb -i protein.pdb -m A:386:ASN:ALA          # one mutation
md_generate_pdb -i protein.pdb -m A:386:ALA              # WT AA optional
md_generate_pdb -i protein.pdb -c mutations.csv          # one PDB per CSV row
md_generate_pdb -i protein.pdb -c mutations.csv --compound   # all on one pose
md_generate_pdb -i protein.pdb -c mutations.csv -o ./mutants/
```

`-m` and `-c` are mutually exclusive. CSV columns are `chain,position,wtAA,prAA`
(case-insensitive; `old_aa`/`new_aa`/`from`/`to` also accepted; `wtAA` optional
and used only for naming). Either form accepts 1- or 3-letter amino acid codes.

**Output naming differs between the two generators** — match on position, not on
a filename you assumed:

| Written by | Pattern | Example |
|---|---|---|
| `md_generate_pdb` | `{stem}_{WT3}-{CHAIN}{POS}-{NEW1}.pdb` — new residue as a **one-letter** code | `4QJR_clean_ASN-A386-A.pdb` |
| `md_mutate` (singles) | `{sr}_{WT}-{CHAIN}{POS}-{NEW}.pdb` — serial-number prefix, no stem | `1_ASN-A386-HIS.pdb` |
| `md_mutate` (double/triple) | `{sr}_{combination}.pdb` | `10_80_N_386_Q+H_387_E.pdb` |

`--compound` applies every CSV mutation to a **single** pose — use it for a
designed multi-mutant, not to generate a panel of independent variants.

## Feeding docking

`md_mutate` writes `*_mutants.txt` listing every mutant PDB it generated; that
file is direct input for `md_dock -r`. Using the same `-o` directory for both
commands produces one combined report. See `mutadock-docking`.

## Reference

- `references/outputs.md` — every file `md_mutate` writes, in order.
