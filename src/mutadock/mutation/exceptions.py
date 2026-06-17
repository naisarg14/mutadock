class MutationError(Exception):
    """Base exception for the mutation package."""


class PDBFileError(MutationError):
    """Raised when a PDB file cannot be read or parsed."""


class CSVGenerationError(MutationError):
    """Raised when CSV generation fails."""
