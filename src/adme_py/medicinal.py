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
    """Assess a molecule against various medicinal chemistry filters and properties.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    properties : dict[str, Union[bool, str, Dict[str, str], float]]
        A dictionary containing the calculated medicinal chemistry properties:
        - "pains": True if the molecule triggers the PAINS filter, False otherwise.
        - "brenk": True if the molecule triggers the Brenk filter, False otherwise.
        - "zinc": True if the molecule triggers the ZINC filter, False otherwise.
        - "leadlikeness": "Yes" if the molecule exhibits lead-likeness,
                          or a dictionary of violations if it doesn't.
        - "synthetic_accessibility": The synthetic accessibility score of the molecule.
    """
    pains: bool = calculate_pains(mol)
    brenk: bool = calculate_brenk(mol)
    zinc: bool = calculate_zinc(mol)
    leadlikeness: Union[str, dict[str, str]] = calculate_leadlikeness(mol)
    synthetic_accessibility: float = calculate_synthetic_accessibility(mol)

    properties: dict[str, Union[bool, str, dict[str, str], float]] = {
        "pains": pains,
        "brenk": brenk,
        "zinc": zinc,
        "leadlikeness": leadlikeness,
        "synthetic_accessibility": synthetic_accessibility,
    }

    return properties


def calculate_pains(mol: Chem.Mol) -> bool:
    """Check if a molecule triggers the PAINS (Pan Assay Interference Compounds) filter.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    bool
        True if the molecule triggers the PAINS filter, False otherwise.

    Notes
    -----
    The PAINS filter is described in:
    J.B. Baell and G.A. Holloway, J. Med. Chem. 2010, 53, 7, 2719–2740
    https://doi.org/10.1021/jm901137j
    """
    params_pains = FilterCatalogParams()
    params_pains.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_A)
    params_pains.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_B)
    params_pains.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS_C)

    catalog_pains = FilterCatalog(params_pains)
    return catalog_pains.HasMatch(mol)


def calculate_brenk(mol: Chem.Mol) -> bool:
    """Check if a molecule triggers the Brenk filter.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    bool
        True if the molecule triggers the Brenk filter, False otherwise.

    Notes
    -----
    The Brenk filter is described in:
    R. Brenk, et al., ChemMedChem, 2008, 3: 435-444.
    https://doi.org/10.1002/cmdc.200700139
    """
    params_brenk = FilterCatalogParams()
    params_brenk.AddCatalog(FilterCatalogParams.FilterCatalogs.BRENK)

    catalog_brenk = FilterCatalog(params_brenk)
    return catalog_brenk.HasMatch(mol)


def calculate_zinc(mol: Chem.Mol) -> bool:
    """
    Check if a molecule triggers the ZINC filter.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    bool
        True if the molecule triggers the ZINC filter, False otherwise.

    Notes
    -----
    The ZINC filter is described in:
    J.J. Irwin, et al. JCIM, 2012, 52 (7), 1757-1768
    https://doi.org/10.1021/ci3001277
    """
    params_zinc = FilterCatalogParams()
    params_zinc.AddCatalog(FilterCatalogParams.FilterCatalogs.ZINC)

    catalog_zinc = FilterCatalog(params_zinc)
    return catalog_zinc.HasMatch(mol)


def calculate_leadlikeness(mol: Chem.Mol) -> Union[str, dict[str, str]]:
    """Assess if a molecule exhibits lead-like properties.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    Union[str, dict[str, str]]
        - "Yes" if the molecule meets the lead-likeness criteria.
        - A dictionary with the violated criteria as keys and descriptive messages as values.

    Notes
    -----
    The method to calculate leadlikeness is described in:
    S.J. Teague, et al. Angew Chem Int Ed Engl. 1999 Dec 16;38(24):3743-3748.
    https://doi.org/10.1002/(SICI)1521-3773(19991216)38:24<3743::AID-ANIE3743>3.0.CO;2-U
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


def calculate_synthetic_accessibility(mol: Chem.Mol) -> float:
    """Calculate the synthetic accessibility score of a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    float
        The synthetic accessibility score (SA_Score) of the molecule.

    Notes
    -----
    The method for calculating synthetic accessibility is detailed in:
    P. Ertl and A. Schuffenhauer, J Cheminform 1, 8 (2009).
    https://doi.org/10.1186/1758-2946-1-8

    and adapted for RDKit by P. Ertl and G. Landrum
    """
    return sascorer.calculateScore(mol)
