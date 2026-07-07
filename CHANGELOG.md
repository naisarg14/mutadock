# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Per-item checkpoint/resume for ΔΔG.** `calc_ddg` / `calc_double_ddg` / `calc_triple_ddg` (and `md_ddg_single` / `md_ddg_double` / `md_ddg_triple` via a new `--resume` flag) can now continue an interrupted run: the output CSV is its own checkpoint, so an interruption partway through keeps the rows already computed and only the missing mutations/combinations are recomputed — instead of discarding everything and starting over. Each completed row is flushed immediately, a torn trailing line from a mid-write kill is dropped and recomputed cleanly, and resume refuses to mix rows written under a different `ddG_protocol` or `n_replicates` (it starts fresh in that case). `md_mutate` uses this automatically (driven by the existing append default), upgrading its previous coarse "skip only if the whole file is complete" check to true per-item resume — the same granularity the docking side already had via `*_completed.txt`.
- **Subprocess timeouts for external tools.** `mk_prepare_receptor`, AutoSite, and Vina are now each run with a wall-clock timeout so a single hung job can no longer stall an entire batch — a timeout raises the tool's domain error and the batch loop skips that receptor-ligand combination. Defaults are 15 min (receptor prep), 30 min (AutoSite), and 60 min (Vina), each overridable via `MUTADOCK_RECEPTOR_PREP_TIMEOUT` / `MUTADOCK_AUTOSITE_TIMEOUT` / `MUTADOCK_VINA_TIMEOUT` (seconds; set to `0` to disable that timeout).
- `read_partial_ddg()` helper in `mutation/helpers.py` — reads the complete rows of a partial ΔΔG CSV for resume, dropping any torn trailing line and enforcing a settings match.
- Regression tests: ΔΔG resume identity for single/double/triple (resumed run reproduces the exact row set, `sr` numbering, and `combination` names of a from-scratch run), `read_partial_ddg` guards, and TimeoutExpired→domain-error routing for all three subprocesses (including a real sleeping-subprocess timeout).

## [2.1.0] - 2026-07-03

### Added
- **Research-grade ΔΔG protocols.** `md_mutate` / `md_quick` / `md_ddg_*` now accept `--protocol {fast,min,cartesian}` (**default `min`**), `--replicates N` (independent estimates → **mean ± SD**, since Rosetta packing is stochastic), `--pack-radius Å`, and `--reu-to-kcal F`. `min` = repack + backbone/side-chain minimization (reliable, ~12 min/proteome); `cartesian` = `ref2015_cart` + Cartesian minimization (most accurate, `cartesian_ddg` style — Park et al. 2016, *J. Chem. Theory Comput.* 12(12):6201, doi:10.1021/acs.jctc.6b00819; ~30 min); `fast` = single repack, no minimization (~1 min, **screening only** — see the Changed note on the default protocol). The chosen protocol is recorded in a `ddG_protocol` column and surfaced in the report (fast-protocol reports carry a "SCREENING ONLY" banner).
- **ΔΔG reported in both REU and kcal/mol.** ΔΔG CSVs gain `ddG_sd`, `ddG_kcal`, `n_replicates`, and `ddG_protocol` columns (and the `double_`/`triple_` equivalents); `ddG_value` remains Rosetta Energy Units (REU). The REU→kcal/mol factor lives in **`src/mutadock/mutation/predict_ddG.py` as `REU_TO_KCAL_SCALE` (default 0.34 ≈ 1/2.94)** and is overridable per run with `--reu-to-kcal`. Reports now show REU + kcal/mol + ± SD and state the units, protocol, and scaling factor explicitly.
- **`md_quick` — the one-command demo.** `md_quick --pdb-id 4QJR --mutation A:386:ASN:HIS --ligand-code imatinib` fetches the structure from RCSB, computes the ΔΔG for that single mutation, builds the mutant, fetches the ligand from PubChem, docks it against the mutant (AutoSite auto-box, or `-c` config), and writes the standard `reports/` bundle — end-to-end in well under a minute. Accepts `-i` local PDB/CIF and `--ligand-file` as alternatives; docking failures still yield the ΔΔG + report.
- **`fetch_ligand(code)`** in `docking/vina_helper.py` — downloads a 3D ligand SDF from PubChem by CID or name (auto-detected; explicit `cid:` / `name:` prefixes accepted), URL-encoding names and validating that the download is a usable SDF (guards against PubChem error pages).
- **Auto-generated run reports** — every `md_mutate` and `md_dock` run now writes a self-contained `report.html` and a `report.pptx` deck (opt out with `--no-report`). The HTML bundles a ΔΔG distribution, a top-stabilizing-mutations table + chart, a docking-affinity chart + table, and an **interactive 3D view of the top mutant pose** (NGL, inlined so the single file opens offline / when emailed). The PPTX mirrors it as presentation slides. Everything is collected in a **`reports/` folder** inside the run's output directory, with each chart also saved as a standalone PNG under `reports/figures/`.
- New `md_report -d <run-dir>` command to (re)build reports from an existing output directory (`--html-only` / `--pptx-only` / `--top N` / `-o/--output-dir`). Using the same `-o` directory for `md_mutate` then `md_dock` yields one combined report.
- New `mutadock.report` subpackage with a reusable `generate_report(output_dir, ...)` API and a documented `RunData` data contract (`discover_run`) that auto-discovers and re-sorts a run's CSVs (robust to the pre-existing triple-sort column quirk) and resolves the top mutant PDB.
- New runtime dependencies for reporting: `matplotlib`, `python-pptx`, `jinja2`. Vendored `ngl.min.js` (NGL 2.3.1) ships in the wheel for the offline 3D viewer.
- `md_mutate -o/--output-dir DIR` writes every output (converted PDB, cleaned PDB, mutation/ddG CSVs, `*_mutants.txt`, and the `mutation_*` folder of mutant PDBs) to `DIR` instead of alongside the input file; the directory is created if absent. Paths recorded in `*_mutants.txt` point into `DIR`, so they remain valid input for `md_dock`
- `md_dock -o/--output-dir DIR` writes docking outputs (poses, logs, `docking_results.csv`, split SDFs, and the `*_completed.txt` resume file) to `DIR` instead of a per-receptor `out/` folder. Prepared PDBQT and AutoSite caches still live next to their inputs so they can be reused across runs
- `md_mutate --pdb-id <ID>` fetches a structure directly from RCSB (e.g. `--pdb-id 4QJR`) instead of requiring a local file with `-i`; the two are mutually exclusive and one is required
- `fetch_pdb()` helper in `mutation/helpers.py` (downloads `https://files.rcsb.org/download/<ID>.pdb`, validates the 4-character code, reuses an already-downloaded file)
- Actionable "residue not found" errors: `md_mutate` now reports the available chains (when the chain is absent) or the chain's residue-numbering span (when only the position is absent) instead of a bare `Residue A386 not found`
- `md_mutate` now warns when the input structure has multiple models (e.g. an NMR ensemble), alternate location indicators, or residue insertion codes — features that silently affect residue numbering
- Test suite with 131 tests across mutation and docking modules (45.45% coverage)
- Regression tests for every correctness fix below (32 new tests): `sort_csv` direction/column handling, first-residue enumeration, non-standard-residue skipping, `clean_pdb` multi-model handling, triple-ddG append predicate, `read_pdb_file` ATOM/HETATM parsing, zero-atom guards, `read_config` values containing `=`, docking-CSV header, `overwrite` flag round-trip, and `.cif`/auto-box docking naming
- `load_matrix()`, `download_matrix()`, `resolve_matrix()` helpers in `mutation/helpers.py` for reading any NCBI-format substitution matrix file
- `--matrix NAME` CLI flag for `md_csv_generator` — auto-downloads from NCBI BLAST FTP if not present locally
- `--matrix-file FILE` CLI flag for `md_csv_generator` — load a custom substitution matrix
- `pyrightconfig.json` for Pylance/pyright import resolution in IDEs
- Coverage badge (45.45%) in README

### Changed
- **Corrected ΔΔG reference (results change).** ΔΔG is now `score(mutant) − score(wild-type *self-mutation* reference)`, where the wild-type side runs the identical repack (+minimization) protocol as the mutant. Previously the wild type was scored raw (unrepacked) while the mutant was repacked, so repacking alone lowered the mutant's energy and biased **every** ΔΔG toward "stabilizing" — a null WT→WT mutation scored ≈ −27 REU instead of ~0. **All ΔΔG values change with this release** (they are lower-magnitude and no longer systematically negative; prior "stabilizing" counts were inflated). Verified: a null self-mutation now yields ≈ 0.00 REU.
- **Default ΔΔG protocol is now `min` (repack + minimization), not a bare repack.** With the symmetric reference, a no-minimization repack can fail to relieve clashes on the wild-type side and *invert* a site's ranking (e.g. a destabilizing mutation at +69 REU appearing as the top stabilizer at −375 REU; minimization corrects it to ≈ −5.5 REU). `min` is reliable and costs ~12 min for a whole-proteome single scan; `--protocol fast` keeps the old no-minimization behavior for quick screening only and its outputs/reports are labelled unreliable.
- `md_csv_generator` defaults to PAM250 matrix (matrix name is case-insensitive)
- PAM250 scoring data now read from `data/PAM250` file instead of a hardcoded function
- Receptor preparation (`prepare_receptor`) now uses `pdbfixer` + `openmm` for structure fixing and `mk_receptor` (meeko CLI) for PDBQT conversion, replacing AutoDockTools_py3
- Both `.pdb` and `.cif` receptor formats are accepted by `md_dock` and `prepare_receptor`
- README rewritten with badges, system requirements table, Quick Start examples, Python API section, CLI reference, and troubleshooting guide
- `read_config` returns typed values `(center: list[float], box_size: list[float], exhaustiveness: int, n_poses: int, n_poses_write: int, overwrite: bool)` instead of a 7-element string tuple
- **Breaking:** `md_mutate` regeneration flag renamed `--noappend` → `--no-append` with corrected semantics — output files are reused by default, and `--no-append` forces regeneration (the old flag's name, help, and default were mutually contradictory)
- **Breaking:** `md_vina_dock` / `vina_dock.py` now expose `--overwrite` / `--no-overwrite` (`argparse.BooleanOptionalAction`, default `overwrite=True`) matching the `vina_dock()` default; `dock_vina` emits the matching flag so the value reaches `write_poses()` correctly
- AutoSite auto-box sizing (`np_docking`) now tracks the cluster radius plus a small margin (`radius * 2 + 8 Å`) instead of a ≥100 Å floor that turned targeted docking into a whole-protein search
- `clean_pdb` gained an optional `keep_hetatm` parameter (default `False`, preserving the existing HETATM-stripping behavior)
- `docking_results.csv` is now written with a `sr,name,affinity` header row

### Removed
- Hardcoded `get_scfn_250()` and `get_scfn_1()` functions (~900 lines) from `mutation/Amino.py`
- Dependency on `AutoDockTools_py3` (git-only package from Valdes-Tresanco-MS)
- `requirements.txt` (consolidated into `pyproject.toml`)

### Fixed
- **`prepare_receptor` failed on every real (multi-chain) receptor** with `RDKit molecule has implicit Hs / has N fragments` — it mistakenly used meeko's *ligand* API (`MoleculePreparation`), which only handles single-fragment molecules. It now invokes meeko's proper receptor CLI (`mk_prepare_receptor.py`, the tool it was documented to use) after PDBFixer, correctly preparing multi-chain receptors. This is what unblocks real docking in `md_dock` / `md_quick`.
- `md_dock` wrote the docked pose to a file named `<rec>_<lig>_out.pdb` that actually contained SDF/MOL text (the raw multi-pose Vina output was overwritten in place by `vina_split`). The two artifacts are now separate files with honest extensions: `<rec>_<lig>_out.pdbqt` (raw multi-pose Vina output) and `<rec>_<lig>_out.sdf` (extracted best pose). `add_score_to_csv` derives the run name extension-agnostically (`<stem>` minus `_out`), so `docking_results.csv` is unchanged
- `md_mutate` on a `.cif` input: the pipeline now converts CIF → PDB (via `convert_cif_pdb`) before cleaning, instead of text-cleaning the CIF as if it were PDB — which produced a malformed/near-empty structure and broke every downstream step
- Windows file-descriptor leak in `prepare_receptor`: `tempfile.mkstemp` fd is now closed immediately with `os.close()` before the temp file is written, preventing `PermissionError` on cleanup
- `sort_csv` ignored the sort direction — `ascending` was hard-coded `True`, so `md_csv_sort -d` silently produced ascending output (now honors the `order`/`-a`/`-d` flag)
- `sort_csv` could not sort by column 0 (treated as "not provided" and dropped into an interactive `input()` prompt that hung in non-interactive use), and raised on string column indices from the CLI — now coerces the index to `int`, accepts column 0, and raises `MutationError` instead of prompting when no column is given
- `csv_generator` silently dropped the first residue (N-terminal of the first chain) from every enumerated CSV — an off-by-one that omitted a residue from all downstream ΔΔG and docking sets
- `csv_generator` raised `KeyError` and aborted the whole run on non-standard/modified residues (e.g. `MSE`, `HOH`); such residues are now skipped with a warning without corrupting serial numbering
- `clean_pdb` truncated multi-model / NMR structures at the first `END` marker; it now retains all models and emits a single terminal `END`
- Double- and triple-ddG `--append` "skip if exists" checks never fired, so re-runs always recomputed the (combinatorially expensive) double/triple ΔΔG steps even when complete, correct outputs already existed. The checks compared the CSV's line count against a *permutation* count (`permutations(num, 2/3)`), but the files hold *combination* rows plus a header — counts that can only coincide at n=2. The triple check was additionally broken by a filename mismatch (it looked for `<stem>_triple_ddg.csv` while the calculation writes `<stem>_clean_triple_ddg.csv`). Both now compare against the true expected line count (`combinations(...) + 1`) computed by reusing the calculation modules' own row-selection logic, and the triple output filename is passed through explicitly so the check and the write agree. Single-mutation steps were unaffected. Output filenames are unchanged.
- `read_pdb_file` referenced an undefined/stale `match` for non-`ATOM` lines (`NameError` / mis-parsing) and its regex never matched `HETATM` records despite allowing them; it now resets per line and parses both
- `calculate_geometric_center` / `calculate_radius` raised a bare `ZeroDivisionError` on structures with zero atoms; they now raise a clear `PDBFileError`
- `dock_vina` emitted a `--nooverwrite` flag that `vina_dock.py` did not define, breaking the `overwrite=False` path and inverting the `overwrite=True` path (see Changed)
- `vina_split` referenced an undefined `score` (`NameError`) when a Vina result line was absent, masking the real "no result" condition; it now raises a clear `DockingError` and parses affinities robustly (including scientific notation)
- `read_config` split on every `=`, raising `ValueError` for values containing `=`; now splits on the first `=` only
- `md_dock` produced an invalid `<name>.cifqt` PDBQT filename for `.cif` receptors (only `.pdb` inputs worked); the name is now derived with `Path(receptor).with_suffix(".pdbqt")`
- Hard-coded `python3` interpreter in the docking subprocess call replaced with `sys.executable`, fixing failures on Windows and in virtualenvs
- Broken error messages in `np_docking` that chained `.format(...)` onto an f-string, dropping the error detail; removed the dead `naisarg = np_docking` alias and fixed user-facing typos ("Easaies"→"Easiest", "ingoring"→"ignoring", "EOFE Error")

## [2.0.0] - 2026-02-10

### Added
- Multi-receptor × multi-ligand batch docking via `md_dock` / `np_docking.py`
- `vina_split` to extract best docking pose as SDF with embedded affinity score
- `add_score_to_csv` to aggregate docking affinities into a single CSV
- `read_config` to parse AutoDock Vina configuration files
- `dock_vina` subprocess wrapper for AutoDock Vina CLI
- `calculate_geometric_center` and `calculate_radius` for binding-site estimation
- `md_vina_dock` CLI entry point using the AutoDock Vina Python bindings
- Logger support across all modules
- `md_install_dependencies` script for guided PyRosetta installation
- Double- and triple-mutation ΔΔG combination commands (`md_ddg_double`, `md_ddg_triple`)
- `md_csv_sort` CLI for sorting any output CSV by column
- `openmm` and `pdbfixer` as explicit dependencies in `pyproject.toml`

### Changed
- Migrated packaging to `pyproject.toml` / hatchling build backend
- Receptor preparation pipeline updated to `pdbfixer` + `openmm` + `mk_receptor`

## [1.0.0] - 2024-09-08

### Added
- Initial release
- Single-, double-, and triple-point mutation prediction using PyRosetta
- ΔΔG scoring with PAM250 substitution matrix
- `md_mutate` pipeline: CIF → PDB → mutation CSV → ΔΔG → sorted output
- `md_csv_generator` to enumerate all possible single-residue substitutions
- `md_ddg_single` for per-mutation ΔΔG calculation
- `clean_pdb` utility to strip non-ATOM records from PDB files
- `convert_cif_pdb` utility using BioPython's MMCIFParser
- `Mutation` class for representing amino-acid substitutions
- `in_directory` context manager (workaround for Rosetta whitespace path bug)
- PAM250 matrix bundled in `data/PAM250`
- Sample data: `data/4QJR.cif` and `data/Ligand.sdf`

[Unreleased]: https://github.com/naisarg14/mutadock/compare/v2.0.0...HEAD
[2.0.0]: https://github.com/naisarg14/mutadock/compare/v1.0.0...v2.0.0
[1.0.0]: https://github.com/naisarg14/mutadock/releases/tag/v1.0.0
