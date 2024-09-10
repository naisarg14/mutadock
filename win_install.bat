echo Installing required Python packages...
pip install -r requirements.txt

echo Installing PyRosetta...
python3 -c "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()"

echo Done.
pause
