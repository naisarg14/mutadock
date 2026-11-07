---
name: mutadock-quickstart
description: Run MUTADOCK end to end in one command with md_quick — fetch a structure from RCSB, score one mutation's ΔΔG, build the mutant, fetch a ligand from PubChem, dock it, and write a report. Use for a first run after installation, a smoke test, a demo, or any single-mutation "does this substitution change binding?" question that does not need a full scan.
---

# MUTADOCK quickstart (`md_quick`)

`md_quick` is the whole pipeline in one command: fetch + clean the structure,
compute ΔΔG for one mutation, build the mutant PDB, fetch or read the ligand,
find the pocket, dock, and write the `reports/` bundle.

```bash
md_quick --pdb-id 4QJR --mutation A:386:ASN:HIS --ligand-code imatinib
```

Outputs land in `./mdquick_<id>/` unless `-o DIR` is given.

## Inputs

**Structure** — exactly one of:
- `--pdb-id 4QJR` — 4-character RCSB accession, downloaded automatically.
- `-i protein.pdb` / `-i protein.cif` — a local structure.

**Mutation** — `-m/--mutation`, required, as `CHAIN:POSITION:WTAA:NEWAA` or
`CHAIN:POSITION:NEWAA` when the wild type is unknown:

```
A:386:ASN:HIS      A:386:H
```

**Ligand** — exactly one of:
- `--ligand-code CODE` — PubChem lookup. Accepts a bare CID (`2244`), a bare
  name (`imatinib`), or an explicit form (`cid:5291`, `name:aspirin`).
- `--ligand-file lig.sdf` — a local ligand that already has 3D coordinates and
  explicit hydrogens.

## Options

| Flag | Meaning |
|---|---|
| `-o DIR` | Output directory (default `./mdquick_<id>`) |
| `-c config.txt` | Vina config with an explicit box — skips AutoSite |
| `-e N` | Vina exhaustiveness |
| `--seed N` | Vina RNG seed (default 19; `0` means random and is **not** reproducible) |
| `--no-report` | Skip the HTML + PPTX report |
| `-q/--quiet` | Suppress informational logging |

`md_quick` also accepts the shared ΔΔG protocol flags: `--protocol`,
`--replicates`, `--pack-radius`, `--reu-to-kcal` (see below).

Without `-c`, pocket detection needs the `autosite` binary (ADFRsuite) on PATH.
Pass a config file instead if AutoSite is unavailable — see `mutadock-docking`
for the config format.

## Smoke test after installation

```bash
md_quick --pdb-id 4QJR --mutation A:386:ASN:HIS --ligand-code imatinib -o smoke/
```

Success means `smoke/reports/report.html` exists and the run printed a ΔΔG and a
binding affinity. The repository also ships local inputs in `data/`
(`4QJR.cif`, `Ligand.sdf`) for an offline version:

```bash
md_quick -i data/4QJR.cif -m A:386:ASN:HIS --ligand-file data/Ligand.sdf -c config.txt -o smoke/
```

**Treat the printed numbers as results, not as expected output.** They depend on
the input structure, the protocol and the seed, so they will not match a
reference value exactly. Pass `--seed` to reproduce a docking run.

## Reading the two numbers it prints

**ΔΔG is in Rosetta Energy Units (REU), not kcal/mol.** A `ddG_kcal` column
applies a `0.34` scaling factor; override with `--reu-to-kcal FACTOR`. Lower
means more stabilizing; a null WT→WT mutation scores ≈ 0 because both sides run
the same protocol.

**`--protocol` defaults to `min`** (repack + minimization; reliable). Use
`cartesian` with `--replicates 3` for research-grade numbers. `--protocol fast`
is **screening only**: a single repack with no minimization can leave clashes it
cannot relieve, and the wild-type self-reference means a failed WT repack can
rank a destabilizing mutation as the top stabilizer. Its reports carry a
"SCREENING ONLY — values unreliable" banner. Never quote `fast` numbers.

**A Vina affinity is only comparable to another affinity from the same box and
seed.** To ask whether the mutation changed binding, dock the wild type the same
way and subtract — see `mutadock-results`.

## Where to go next

- A scan over many positions, or double/triple mutants → `mutadock-mutation`.
- Many receptors × many ligands, or a custom box → `mutadock-docking`.
- Interpreting the CSVs and reports → `mutadock-results`.
- Chain, residue or ligand errors → `mutadock-debugging`.
