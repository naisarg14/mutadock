#!/usr/bin/env bash
set -euo pipefail

VENV_DIR="${1:-.venv}"

echo "==> Creating virtual environment in '$VENV_DIR'..."
python -m venv "$VENV_DIR"

echo "==> Activating virtual environment..."
# shellcheck disable=SC1091
source "$VENV_DIR/bin/activate"

echo "==> Installing mutadock and its dependencies..."
python -m pip install --upgrade pip
python -m pip install .

echo "==> Installing PyRosetta..."
python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'

echo ""
echo "Installation complete."
echo "Activate the environment with: source $VENV_DIR/bin/activate"
