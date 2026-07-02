# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Test suite with 131 tests across mutation and docking modules (45.45% coverage)
- Regression tests for every correctness fix below (32 new tests): `sort_csv` direction/column handling, first-residue enumeration, non-standard-residue skipping, `clean_pdb` multi-model handling, triple-ddG append predicate, `read_pdb_file` ATOM/HETATM parsing, zero-atom guards, `read_config` values containing `=`, docking-CSV header, `overwrite` flag round-trip, and `.cif`/auto-box docking naming
- `load_matrix()`, `download_matrix()`, `resolve_matrix()` helpers in `mutation/helpers.py` for reading any NCBI-format substitution matrix file
- `--matrix NAME` CLI flag for `md_csv_generator` — auto-downloads from NCBI BLAST FTP if not present locally
- `--matrix-file FILE` CLI flag for `md_csv_generator` — load a custom substitution matrix
- `pyrightconfig.json` for Pylance/pyright import resolution in IDEs
- Coverage badge (45.45%) in README

### Changed
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
- Windows file-descriptor leak in `prepare_receptor`: `tempfile.mkstemp` fd is now closed immediately with `os.close()` before the temp file is written, preventing `PermissionError` on cleanup
- `sort_csv` ignored the sort direction — `ascending` was hard-coded `True`, so `md_csv_sort -d` silently produced ascending output (now honors the `order`/`-a`/`-d` flag)
- `sort_csv` could not sort by column 0 (treated as "not provided" and dropped into an interactive `input()` prompt that hung in non-interactive use), and raised on string column indices from the CLI — now coerces the index to `int`, accepts column 0, and raises `MutationError` instead of prompting when no column is given
- `csv_generator` silently dropped the first residue (N-terminal of the first chain) from every enumerated CSV — an off-by-one that omitted a residue from all downstream ΔΔG and docking sets
- `csv_generator` raised `KeyError` and aborted the whole run on non-standard/modified residues (e.g. `MSE`, `HOH`); such residues are now skipped with a warning without corrupting serial numbering
- `clean_pdb` truncated multi-model / NMR structures at the first `END` marker; it now retains all models and emits a single terminal `END`
- Triple-ddG "skip if exists" check compared against the double-mutation count (`permutations(num_double_ddg, 2)`) instead of `permutations(num_triple_ddg, 3)`, so `--append` re-runs could skip or recompute triple results incorrectly
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
