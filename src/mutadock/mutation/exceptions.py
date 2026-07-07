class MutationError(Exception):
    """Base exception for the mutation package."""


class PDBFileError(MutationError):
    """Raised when a PDB file cannot be read or parsed."""


class CSVGenerationError(MutationError):
    """Raised when CSV generation fails."""


class ResidueMismatchError(MutationError):
    """Raised when the expected wild-type residue doesn't match the structure.

    Typically means the caller's (chain, position) numbering convention
    doesn't match the PDB file's own residue numbering.
    """
