"""Methods to calculate solubility parameters."""

from typing import Union

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors


def calculate_all_solubility(mol: Chem.Mol) -> dict[str, Union[str, float]]:
    """Calculate solubility-related properties of a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    dict[str, Union[str, float]]
        A dictionary containing the calculated solubility properties:
        - "log_s_esol": The estimated ESOL value (log mol/L) using the Delaney model (float).
        - "solubility_esol": The estimated solubility in mol/L (float).
        - "class_esol": The solubility class based on the ESOL value (str).
    """
    esol: float = calculate_esol(mol)
    solubility_class: str = _calculate_solubility_class(esol)

    properties: dict[str, Union[str, float]] = {
        "log_s_esol": esol,
        "solubility_esol": 10**esol,
        "class_esol": solubility_class,
    }

    return properties


def calculate_esol(mol: Chem.Mol) -> float:
    """Calculate the aqueous solubility (ESOL) using the Delaney model.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    float
        The estimated ESOL value (log mol/L).

    Notes
    -----
    - This function utilizes the method and coefficients calculated by Pat Walters:
      https://github.com/PatWalters/solubility/blob/master/esol.py
    - The ESOL estimation method is provided by Delaney 2004:
      https://doi.org/10.1021/ci034243x
    """
    logp: float = Crippen.MolLogP(mol)
    mw: float = Descriptors.MolWt(mol)
    rotors: int = Descriptors.NumRotatableBonds(mol)
    aromatic_prop: float = _calculate_aromatic_proportion(mol)

    intercept: float = 0.26121066137801696
    coef: dict[str, float] = {
        "logp": -0.7416739523408995,
        "mw": -0.0066138847738667125,
        "rotors": 0.003451545565957996,
        "ap": -0.42624840441316975,
    }

    esol: float = (
        intercept
        + coef["logp"] * logp
        + coef["mw"] * mw
        + coef["rotors"] * rotors
        + coef["ap"] * aromatic_prop
    )

    return esol


def _calculate_aromatic_proportion(mol: Chem.Mol) -> float:
    """Calculate the proportion of aromatic atoms in a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    float
        The aromatic proportion of the molecule.
    """
    return len(mol.GetAromaticAtoms()) / mol.GetNumAtoms()


def _calculate_solubility_class(log_s: float) -> str:
    """Determine the solubility class based on the given LogS value.

    Parameters
    ----------
    log_s : float
        The LogS value.

    Returns
    -------
    str
        The solubility class.
    """
    if log_s < -10:
        return "Insoluble"
    elif -10 <= log_s < -6:
        return "Poorly Soluble"
    elif -6 <= log_s < -4:
        return "Moderately Soluble"
    elif -4 <= log_s < -2:
        return "Soluble"
    elif -2 <= log_s < 0:
        return "Very Soluble"
    else:
        return "Highly Soluble"
