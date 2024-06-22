"""Scripts to handle the various chemical identifier inputs provided to the ADME class."""

from rdkit import Chem
from rdkit.Chem import Mol


def mol_from_identifier(identifier: str) -> Mol:
    """Create RDKit Mol object from identifier."""
    if not isinstance(identifier, str):
        raise ValueError("Input identifier must be a string.")

    if _is_smiles(identifier):
        mol = Chem.MolFromSmiles(identifier)
    elif _is_smarts(identifier):
        mol = Chem.MolFromSmarts(identifier)
    elif _is_inchi(identifier):
        mol = Chem.MolFromInchi(identifier)
    else:
        raise ValueError("Unsupported string identifier format.")

    return mol


def _is_smiles(identifier: str) -> bool:
    """Check if the identifier is a valid SMILES string."""
    try:
        return Chem.MolFromSmiles(identifier) is not None
    except Exception:
        return False


def _is_smarts(identifier: str) -> bool:
    """Check if the identifier is a valid SMARTS string."""
    try:
        return Chem.MolFromSmarts(identifier) is not None
    except Exception:
        return False


def _is_inchi(identifier: str) -> bool:
    """Check if the identifier is a valid Inchi string."""
    try:
        return Chem.MolFromInchi(identifier) is not None
    except Exception:
        return False
