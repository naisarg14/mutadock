# MUTADOCK agent guide

MUTADOCK runs protein mutation studies and receptor–ligand docking as a set of
composable CLI steps: substitution enumeration → PyRosetta ΔΔG scoring → mutant
PDB generation → AutoDock Vina batch docking → HTML/PPTX reports. Act as a
computational-biology assistant: make the question concrete, run the smallest
informative job first, read the numbers with their caveats attached, and spend
CPU hours only once a pilot has shown the setup is correct.

## Route to one skill

Open `.agents/skills/INDEX.md`, then load only the skill the current action
needs. Do not read several skills to find the right one.

- `mutadock-setup` — install or verify PyRosetta, Vina, meeko, pdbfixer, AutoSite.
- `mutadock-quickstart` — one-shot `md_quick` demo or a first end-to-end check.
- `mutadock-mutation` — plan or run a ΔΔG scan; build specific mutant PDBs.
- `mutadock-docking` — batch docking, search box, reproducibility.
- `mutadock-results` — read the CSVs and reports; compute ΔΔG_binding.
- `mutadock-debugging` — chain, residue, ligand, or structure input failures.
- `mutadock-python-api` — script MUTADOCK from Python instead of the CLI.

## Rules that hold regardless of which skill is loaded

**ΔΔG units.** `ddG_value` is in **Rosetta Energy Units (REU), not kcal/mol**.
The `ddG_kcal` column applies a scaling factor of `0.34` (ref2015 `cartesian_ddg`
convention, Park et al. 2016). Override with `--reu-to-kcal FACTOR`. Lower ΔΔG =
more stabilizing; a null WT→WT mutation scores ≈ 0 because both sides run the
identical protocol.

**Protocol.** `--protocol` defaults to `min` (repack + minimization; reliable).
Use `cartesian` for research-grade numbers, ideally with `--replicates 3`.
`fast` is **screening only** — a single repack with no minimization can leave
clashes it cannot relieve, and because ΔΔG uses a wild-type self-reference a
failed WT repack can make a destabilizing mutation rank as the top stabilizer.
Never report `fast` numbers as results.

**Cost.** A default `md_mutate` run takes hours. The runtime drivers are
`-s/--double` (top-N singles → C(N,2) pairs; default 100 → 4950) and
`-t/--triple` (top-N double rows → C(k,3)); `-S/-D/-T` only cap how many mutant
PDBs are written. The single-mutation ΔΔG pass covers **every** positive-score
substitution and is bounded by no flag. Always pilot before a full run, and for
a site-focused question filter the `md_csv_generator` CSV and call
`md_ddg_single` directly rather than running `md_mutate`.

**Resume.** `md_mutate` appends to its checkpoint CSVs and skips completed
items; `md_ddg_single|double|triple` take `--resume`; `md_dock` skips pairs
listed in `*_completed.txt`. Use `--no-append` / `--ignore_existing` only when a
fresh run is intended.

**Don't re-implement, and don't patch the package.** Reuse the helpers in
`src/mutadock/mutation/helpers.py` (matrices, cleaning, structure warnings,
residue maps) and `src/mutadock/docking/vina_helper.py` (preparation, config
parsing, AutoSite, Vina). Every parameter is a CLI flag or an environment
variable — `MUTADOCK_DATA_DIR`, `MUTADOCK_RECEPTOR_PREP_TIMEOUT`,
`MUTADOCK_AUTOSITE_TIMEOUT`, `MUTADOCK_VINA_TIMEOUT` — so editing installed code
is never the right fix.

**Interpretation.** `clean_pdb` strips HETATM records and keeps one conformer,
so ΔΔG is computed on the apo, single-conformer protein. A Vina affinity is only
comparable against another affinity scored in the **same box with the same
seed**; `docking_results.csv` records both for exactly that reason.

## Repository conventions

Source lives in `src/mutadock/{mutation,docking,report}`; tests mirror that
layout under `tests/`. Run `pytest -m "not slow"` for the fast suite. Formatting
is black + ruff + isort (line length 88) via pre-commit. The agent skills are
canonical in `.agents/skills/`; `.claude/skills/` is a generated mirror — edit
the former and run `python scripts/sync_agent_skills.py`.
