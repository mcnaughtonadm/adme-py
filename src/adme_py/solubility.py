"""Methods to calculate solubility parameters."""

from typing import Union

from rdkit import Chem
from rdkit.Chem import Crippen, Descriptors


def calculate_all_solubility(mol: Chem.Mol) -> dict[str, Union[str, float, int]]:
    """Calculate all physiochemical properties of a given molecule.

    Parameters
    ----------
        mol : Chem.Mol
             The input rdkit Mol object
    """
    esol = calculate_esol(mol)

    properties = {
        "log_s_esol": esol,
        "solubility_esol": 10**esol,
        "class_esol": _calculate_solubility_class(esol),
    }

    return properties


def calculate_esol(mol: Chem.Mol) -> float:
    """Calculate the aqueous solubility (ESOL) using the Delaney model.

    Method and coefficients calculated by PatWalters https://github.com/PatWalters/solubility/blob/master/esol.py

    Estimation method provided by Delaney 2004 (https://doi.org/10.1021/ci034243x)

    Parameters
    ----------
        mol: Chem.Mol
             The input RDKit Mol object.

    Returns
    -------
        The estimated ESOL value (log mol/L).
    """
    logp = Crippen.MolLogP(mol)
    mw = Descriptors.MolWt(mol)
    rotors = Descriptors.NumRotatableBonds(mol)
    aromatic_prop = _calculate_aromatic_proportion(mol)

    intercept = 0.26121066137801696
    coef = {
        "logp": -0.7416739523408995,
        "mw": -0.0066138847738667125,
        "rotors": 0.003451545565957996,
        "ap": -0.42624840441316975,
    }

    esol = (
        intercept
        + coef["logp"] * logp
        + coef["mw"] * mw
        + coef["rotors"] * rotors
        + coef["ap"] * aromatic_prop
    )

    return esol


def _calculate_aromatic_proportion(mol: Chem.Mol) -> float:
    """Calculate the aromatic proportion of the molecule.

    Paramters
    ---------
        mol: Chem.Mol
             The input RDKit Mol object.

    Returns
    -------
        result: float
             The aromatic proportion of the molecule
    """
    return len(mol.GetAromaticAtoms()) / mol.GetNumAtoms()


def _calculate_solubility_class(log_s: float) -> str:
    """Determine the solubility class of the given LogS value.

    Parameters
    ----------
        log_s: float
            The LogS value calculated.

    Returns
    -------
        result: str
            The solubility class given the LogS value.
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
