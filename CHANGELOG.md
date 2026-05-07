# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Test suite with 131 tests across mutation and docking modules (45.45% coverage)
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

### Removed
- Hardcoded `get_scfn_250()` and `get_scfn_1()` functions (~900 lines) from `mutation/Amino.py`
- Dependency on `AutoDockTools_py3` (git-only package from Valdes-Tresanco-MS)
- `requirements.txt` (consolidated into `pyproject.toml`)

### Fixed
- Windows file-descriptor leak in `prepare_receptor`: `tempfile.mkstemp` fd is now closed immediately with `os.close()` before the temp file is written, preventing `PermissionError` on cleanup

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
