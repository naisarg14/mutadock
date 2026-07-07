# syntax=docker/dockerfile:1
#
# MUTADOCK image: conda-forge stack (OpenMM, PDBFixer, meeko, RDKit, AutoDock
# Vina bindings) + PyRosetta, pre-wired so `md_mutate` / `md_dock` / `md_quick`
# work out of the box. This is the recommended way to run mutadock — it
# sidesteps building PyRosetta/OpenMM from source and matching system
# toolchains by hand.
#
# Build:
#   docker build -t mutadock .
# Run the bundled demo:
#   docker run --rm -it mutadock md_quick --pdb-id 4QJR --mutation A:386:ASN:HIS --ligand-code imatinib
# Mount your own data and collect outputs on the host:
#   docker run --rm -it -v "$PWD/work:/work" -w /work mutadock md_mutate -i /work/protein.pdb

FROM condaforge/miniforge3:26.3.2-3

LABEL org.opencontainers.image.title="mutadock" \
      org.opencontainers.image.description="Protein mutation studies (PyRosetta) and multi-receptor/multi-ligand docking (AutoDock Vina)" \
      org.opencontainers.image.source="https://github.com/naisarg14/mutadock" \
      org.opencontainers.image.licenses="GPL-3.0"

WORKDIR /app

# Install the conda-forge stack first (see environment.yml for what's pinned
# and why). Kept as its own layer, ahead of the source copy, so it's only
# rebuilt when the environment spec changes, not on every source edit.
COPY environment.yml .
RUN mamba env create -f environment.yml -n mutadock && \
    mamba clean -afy

# PyRosetta has its own academic/commercial license and isn't distributed on
# conda-forge/PyPI directly — pyrosetta-installer (pinned in environment.yml)
# fetches the matching prebuilt wheel from RosettaCommons at build time.
# This step downloads several GB and needs network access. It doesn't depend
# on mutadock's own source, so it stays cached across source-only rebuilds.
RUN mamba run -n mutadock python -c \
    "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()"

# Now bring in the rest of the source and install mutadock itself. This is
# the layer that actually changes on every source edit.
COPY . .
RUN mamba run -n mutadock python -m pip install --no-deps .

# Make the mutadock env the default without requiring `conda activate` or
# `mamba run` in every command.
ENV PATH=/opt/conda/envs/mutadock/bin:$PATH
SHELL ["/bin/bash", "-c"]

# AutoSite (ADFRsuite), used by `md_quick` for automatic pocket detection
# when no config file is given, has no conda/PyPI package and isn't installed
# here — see README.md. `md_quick`/`md_dock` still work fully with an
# explicit `-c config.txt` / center+box.

ENTRYPOINT []
CMD ["bash"]
