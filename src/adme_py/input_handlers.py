"""Scripts to handle the various chemical identifier inputs provided to the ADME class."""

from rdkit import Chem
from rdkit.Chem import Mol


def mol_from_identifier(identifier: str) -> Mol:
    """Create an RDKit Mol object from a chemical identifier.

    Parameters
    ----------
    identifier : str
        The chemical identifier string. Supported formats are SMILES, SMARTS, and InChI.

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        The RDKit Mol object corresponding to the provided identifier.

    Raises
    ------
    ValueError
        If the identifier is not a string or if its format is not supported.
    """
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
    """Check if the identifier is a valid SMILES string.

    Parameters
    ----------
    identifier : str
        The identifier string to check.

    Returns
    -------
    bool
        True if the identifier is a valid SMILES, False otherwise.
    """
    try:
        return Chem.MolFromSmiles(identifier) is not None
    except Exception:
        return False


def _is_smarts(identifier: str) -> bool:
    """Check if the identifier is a valid SMARTS string.

    Parameters
    ----------
    identifier : str
        The identifier string to check.

    Returns
    -------
    bool
        True if the identifier is a valid SMARTS, False otherwise.
    """
    try:
        return Chem.MolFromSmarts(identifier) is not None
    except Exception:
        return False


def _is_inchi(identifier: str) -> bool:
    """Check if the identifier is a valid InChI string.

    Parameters
    ----------
    identifier : str
        The identifier string to check.

    Returns
    -------
    bool
        True if the identifier is a valid InChI, False otherwise.
    """
    try:
        return Chem.MolFromInchi(identifier) is not None
    except Exception:
        return False
