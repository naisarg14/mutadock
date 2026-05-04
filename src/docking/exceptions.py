class DockingError(Exception):
    """Base exception for the docking package."""


class PDBFileError(DockingError):
    """Raised when a PDB file cannot be read or parsed."""


class ConfigError(DockingError):
    """Raised when a Vina config file cannot be read or parsed."""


class LigandPreparationError(DockingError):
    """Raised when ligand PDBQT preparation fails."""


class ReceptorPreparationError(DockingError):
    """Raised when receptor PDBQT preparation fails."""


class DockingRunError(DockingError):
    """Raised when an AutoDock Vina run fails."""
