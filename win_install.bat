@echo off
setlocal EnableDelayedExpansion

set VENV_DIR=.venv
if not "%~1"=="" set VENV_DIR=%~1

echo =^> Creating virtual environment in '%VENV_DIR%'...
python -m venv "%VENV_DIR%"
if errorlevel 1 (
    echo ERROR: Failed to create virtual environment. Is Python installed?
    exit /b 1
)

echo =^> Activating virtual environment...
call "%VENV_DIR%\Scripts\activate.bat"
if errorlevel 1 (
    echo ERROR: Failed to activate virtual environment.
    exit /b 1
)

echo =^> Installing mutadock and its dependencies...
python -m pip install --upgrade pip
python -m pip install .
if errorlevel 1 (
    echo ERROR: Installation failed.
    exit /b 1
)

echo =^> Installing PyRosetta...
python -c "import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()"
if errorlevel 1 (
    echo ERROR: PyRosetta installation failed.
    exit /b 1
)

echo.
echo Installation complete.
echo Activate the environment with: %VENV_DIR%\Scripts\activate.bat
pause
