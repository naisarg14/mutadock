---
name: mutadock-setup
description: Install or verify a working MUTADOCK environment — PyRosetta, AutoDock Vina, meeko, pdbfixer/openmm, RDKit, and the optional AutoSite binary. Use when MUTADOCK is not installed, when an md_* command fails with ImportError or "command not found" (mk_receptor, vina, autosite), when PyRosetta will not install, or when a collaborator needs the environment reproduced via Docker or conda.
---

# Setting up MUTADOCK

MUTADOCK needs Python 3.11 or 3.12 plus several heavy native dependencies.
Install them in this order and verify before running any real workload — a
half-installed environment fails late, in the middle of an hours-long job.

| Component | Needed for | Comes from |
|---|---|---|
| `mutadock` | everything | PyPI |
| PyRosetta | all ΔΔG scoring and mutant generation | `pyrosetta-installer` |
| AutoDock Vina (`import vina`) | all docking | conda-forge `vina` package |
| meeko (`mk_receptor`) | receptor/ligand PDBQT preparation | pip / conda |
| pdbfixer + openmm | receptor repair before PDBQT conversion | conda-forge |
| RDKit | ligand handling | pip / conda |
| AutoSite (ADFRsuite) | automatic box detection | ADFRsuite, **optional** |

AutoSite is only required when `md_quick` or `md_dock` must find the pocket
itself. Supplying `-c config.txt` with an explicit box avoids it entirely.

## Route 1 — Docker (most reliable)

The bundled `Dockerfile` builds the whole conda-forge stack plus PyRosetta.

```bash
docker build -t mutadock .
```

**The image deliberately does not install AutoSite** (it has no conda/PyPI
package), so `md_quick` inside the container must be given an explicit search
box with `-c config.txt` instead of auto-detecting the pocket:

```bash
docker run --rm -it -v "$PWD/work:/work" -w /work mutadock \
    md_quick --pdb-id 4QJR --mutation A:386:ASN:HIS --ligand-code imatinib -c /work/config.txt
```

See `mutadock-docking` for the config format. For automatic pocket detection,
install ADFRsuite into the image yourself or run outside Docker.

Workflows that need no search box work unchanged:

```bash
docker run --rm -it -v "$PWD/work:/work" -w /work mutadock md_mutate -i /work/protein.pdb
```

## Route 2 — conda

`environment.yml` pins the full stack including the Vina Python bindings.

```bash
conda env create -f environment.yml
conda activate mutadock
pip install .      # REQUIRED - environment.yml pins the dependencies, not MUTADOCK
python -c "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()"
```

`environment.yml` installs the dependency stack — including the conda-forge
`vina` package, which supplies the Python bindings MUTADOCK actually uses — but
**not MUTADOCK itself**. Skip `pip install .` (or `pip install mutadock`) and
none of the `md_*` commands will exist. PyRosetta is installed separately
because its licence keeps it off conda-forge.

## Route 3 — pip

```bash
pip install mutadock
python -m pip install pyrosetta_installer
python -c "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()"
conda install -c conda-forge pdbfixer openmm     # pip wheels for these are fragile
```

Install the AutoDock Vina **Python bindings** separately — they are not a pip
dependency, and the docking module imports `vina` rather than calling a binary:

```bash
conda install -c conda-forge vina
```

See also the [official installation guide](https://autodock-vina.readthedocs.io/en/latest/installation.html).

## Verify before doing real work

Run this check and fix whatever it reports before starting a job:

```bash
python - <<'PY'
import shutil, importlib
for mod in ("mutadock", "pyrosetta", "rdkit", "meeko", "pdbfixer", "openmm", "vina"):
    try:
        importlib.import_module(mod)
        print(f"OK      import {mod}")
    except Exception as e:
        print(f"MISSING import {mod}: {type(e).__name__}: {e}")
for exe in ("mk_receptor", "autosite"):
    path = shutil.which(exe)
    print(f"{'OK     ' if path else 'MISSING'} {exe} {path or ''}")
PY
```

**`vina` is checked as an import, not an executable, and that is deliberate.**
MUTADOCK docks through `vina_dock.py`, which does `from vina import Vina` — the
Python bindings. A `vina` binary on PATH without the module does not work, and
the module without a binary works fine, so `shutil.which("vina")` answers the
wrong question.

`autosite` showing MISSING is acceptable if every docking run passes an explicit
`-c config.txt`. Everything else must be OK.

Then run the smoke test from `mutadock-quickstart`, which exercises fetch →
clean → ΔΔG → mutant → dock → report in one command.

## Install the agent skills for other projects

MUTADOCK ships these skills inside the wheel, so they can be installed into a
coding agent's own skills directory without cloning the repository:

```bash
md_skills list                                # show bundled skills
md_skills install --claude --codex --user     # ~/.claude/skills + ~/.agents/skills
md_skills install --project /path/to/work     # that project only
```

## Common installation failures

**`pyrosetta` import fails.** Use `pyrosetta-installer` rather than a manual
wheel; it handles licence acceptance and platform detection. PyRosetta needs
~6 GB of disk and prefers 16 GB RAM.

**`mk_receptor: command not found`.** meeko is missing or its console scripts
are not on PATH: `pip install meeko`, then re-run the verification block.

**`pdbfixer` / `openmm` import error.** Install both from conda-forge; the pip
route is the usual cause of a broken OpenMM.

**`ModuleNotFoundError: No module named 'vina'`.** Docking needs the Vina
Python bindings, which are intentionally not a pip dependency:
`conda install -c conda-forge vina`. A `vina` executable on PATH is not a
substitute.

**Matrix download fails.** PAM250 and BLOSUM62 ship with the package; others are
fetched from the NCBI BLAST FTP into `~/.cache/mutadock/matrices`. If NCBI is
unreachable, download the matrix by hand and pass `--matrix-file /path/to/PAM30`.
Set `MUTADOCK_DATA_DIR` to relocate the cache.

**External command times out.** Defaults are 900 s for receptor preparation,
1800 s for AutoSite and 3600 s for Vina. Override with
`MUTADOCK_RECEPTOR_PREP_TIMEOUT`, `MUTADOCK_AUTOSITE_TIMEOUT`,
`MUTADOCK_VINA_TIMEOUT`; `0` disables a timeout.

MUTADOCK has no environment-doctor subcommand — use the verification block
above rather than looking for one.
