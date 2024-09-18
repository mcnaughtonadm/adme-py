"""Script containing tools to calculate Medicinal Chemistry properties."""

import os
import sys
from typing import Union

from rdkit import Chem
from rdkit.Chem import RDConfig
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

from adme_py.lipophilicity import calculate_logp_crippen
from adme_py.physiochemical import calculate_molecular_weight, calculate_number_rotatable_bonds

sys.path.append(os.path.join(RDConfig.RDContribDir, "SA_Score"))
import sascorer


def calculate_all_medicinal(mol: Chem.Mol) -> dict[str, Union[bool, str, dict[str, str], float]]:
    """Calculate all properties in the medicinal chemistry script.

    Parameters
    ----------
        mol : Chem.Mol
             The input rdkit Mol object

    Returns
    -------
        properties : dict[str, Union[str,dict[str,str]]]
            Dictionary of properties calculated
    """
    pains: bool = calculate_pains(mol)
    brenk: bool = calculate_brenk(mol)
    zinc: bool = calculate_zinc(mol)
    leadlikeness: Union[str, dict[str, str]] = calculate_leadlikeness(mol)
    synthetic_accessibility: float = calculate_synthetic_accessiblity(mol)

    properties: dict[str, Union[bool, str, dict[str, str], float]] = {
        "pains": pains,
        "brenk": brenk,
        "zinc": zinc,
        "leadlikeness": leadlikeness,
        "synthetic_accessibility": synthetic_accessibility,
    }

    return properties


def calculate_pains(mol: Chem.Mol) -> bool:
    """Calculate whether the molecule triggers the PAINS filter.

    Parameters
    ----------
        mol : Chem.Mol
             The input rdkit Mol object

    Returns
    -------
        bool
             Whether the molecule triggers the PAINS filter.
    """
    params_pains = FilterCatalogParams()
    params_pains.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params_pains.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params_pains.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)

    catalog_pains = FilterCatalog(params_pains)
    return catalog_pains.HasMatch(mol)


def calculate_brenk(mol: Chem.Mol) -> bool:
    """Calculate whether the molecule triggers the Brenk filter.

    Parameters
    ----------
        mol : Chem.Mol
             The input rdkit Mol object

    Returns
    -------
        bool
             Whether the molecule triggers the Brenk filter.
    """
    params_brenk = FilterCatalogParams()
    params_brenk.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)

    catalog_brenk = FilterCatalog(params_brenk)
    return catalog_brenk.HasMatch(mol)


def calculate_zinc(mol: Chem.Mol) -> bool:
    """Calculate whether the molecule triggers the Zinc filter.

    Parameters
    ----------
        mol : Chem.Mol
             The input rdkit Mol object

    Returns
    -------
        bool
             Whether the molecule triggers the ZINC filter.
    """
    params_zinc = FilterCatalogParams()
    params_zinc.AddCatalog(FilterCatalogParams.FilterCatalogs.ZINC)

    catalog_zinc = FilterCatalog(params_zinc)
    return catalog_zinc.HasMatch(mol)


def calculate_leadlikeness(mol: Chem.Mol) -> Union[str, dict[str, str]]:
    """Calculate if the molecule exhibits leadlikeness.

    Parameters
    ----------
        mol : Chem.Mol
            The input rdkit Mol object

    Returns
    -------
        Union[str,dict[str,str]]
            The violations to the leadlikeness filter
    """
    violation = {}

    logp = calculate_logp_crippen(mol)
    if logp > 3.5:
        violation["LogP"] = f"LogP: {logp} > 3.5"

    molecular_weight = calculate_molecular_weight(mol)
    if molecular_weight < 250 or molecular_weight > 350:
        violation["MW"] = f"MW: {molecular_weight} is outside the acceptable range (250-350)"

    number_rotatable_bonds = calculate_number_rotatable_bonds(mol)
    if number_rotatable_bonds > 7:
        violation["num_rot_bonds"] = f"Number of Rotatable Bonds: {number_rotatable_bonds} > 7"

    if violation:
        return violation
    else:
        return "Yes"


def calculate_synthetic_accessiblity(mol: Chem.Mol) -> float:
    """Calculate the synthetic accessibility of the molecule.

    Parameters
    ----------
        mol : Chem.Mol
            The input rdkit Mol object

    Returns
    -------
        float
            The Synthetic Accessibility score of the molecule
    """
    return sascorer.calculateScore(mol)
