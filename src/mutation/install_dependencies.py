def main():
    packages = [
        "git+https://github.com/Valdes-Tresanco-MS/AutoDockTools_py3.git",
        "biopython==1.84",
        "meeko==0.5.0",
        "numpy>=2.1.1",
        "pyrosetta-installer==0.1.1",
        "tqdm>=4.66.5",
        "pyarrow==17.0.0",
        "pandas==2.2.2",
    ]
    import subprocess
    for package in packages:
        subprocess.check_call(["python3", "-m", "pip", "install", package])
    try:
        import pyrosetta_installer
        pyrosetta_installer.install_pyrosetta()
    except ImportError:
        msg = "Error when installing pyrosetta module.\n"
        msg += "Easiest way to fix this is to install pyrosetta_installer using the following command and run again:\n\n"
        msg += "python -m pip install pyrosetta-installer\n"
        msg += "If you already have pyrosetta-installer installed, please check the installation.\n"
        msg += "If the problem persists, please create a github issue or contact developer at naisarg.patel14@hotmail.com"
        print(msg)
        import sys; sys.exit(2)


if __name__ == '__main__':
    main()