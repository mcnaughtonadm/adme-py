"""Tools to calculate druglikeness filters."""

from typing import Union

from rdkit import Chem

from adme_py.lipophilicity import calculate_logp_crippen
from adme_py.physiochemical import (
    calculate_molar_refractivity,
    calculate_molecular_weight,
    calculate_number_hbond_acceptors,
    calculate_number_hbond_donors,
    calculate_number_of_atoms,
    calculate_number_rotatable_bonds,
    calculate_tpsa,
)


def calculate_all_druglikeness(mol: Chem.Mol) -> dict[str, Union[str, dict[str, str]]]:
    """Calculate the all the lipophilicity properties of a given molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input rdkit Mol object

    Returns
    -------
    properties : dict[str, Union[str,dict[str,str]]]
        A dictionary containing the results of the druglikeness filters:
        - "lipinski": "Pass" or a dictionary of violations for Lipinski's Rule of 5
        - "ghose": "Pass" or a dictionary of violations for the Ghose filter
        - "veber": "Pass" or a dictionary of violations for Veber's Rule
    """
    lipinksi: Union[str, dict] = calculate_lipinksi(mol)
    ghose: Union[str, dict] = calculate_ghose(mol)
    veber: Union[str, dict] = calculate_veber(mol)

    properties: dict[str, Union[str, dict]] = {
        "lipinski": lipinksi,
        "ghose": ghose,
        "veber": veber,
    }

    return properties


def calculate_lipinksi(mol: Chem.Mol) -> Union[str, dict[str, str]]:
    """Assess if a molecule violates Lipinski's Rule of 5.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    Union[str, dict[str, str]]
        - "Pass" if the molecule adheres to all rules.
        - A dictionary with the violated rules as keys and descriptive messages as values.

    Notes
    -----
    The Lipinksi Rule of 5 is described in:
    C.A Lipinksi, et al. Adv. Drug Delivery Rev. 2001, 46, 1-3, 3-26
    https://doi.org/10.1016/S0169-409X(00)00129-0
    """
    violation = {}

    # Add explicit hydrogens for filter
    mol = Chem.AddHs(mol)

    molecular_weight = calculate_molecular_weight(mol)
    if molecular_weight > 500:
        violation["MW"] = f"MW: {molecular_weight} > 500 Dalton"

    logp = calculate_logp_crippen(mol)
    if logp > 5:
        violation["LogP"] = f"LogP: {logp} >5"

    hbond_donors = calculate_number_hbond_donors(mol)
    if hbond_donors > 5:
        violation["hbond_donors"] = f"hbond donors: {hbond_donors} >5"

    hbond_acceptors = calculate_number_hbond_acceptors(mol)
    if hbond_acceptors > 10:
        violation["hbond_acceptors"] = f"hbond acceptors: {hbond_acceptors} >10"

    if violation:
        return violation
    else:
        return "Pass"


def calculate_ghose(mol: Chem.Mol) -> Union[str, dict[str, str]]:
    """Evaluate a molecule against the Ghose filter.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    properties : Union[str, dict[str, str]]
        - "Pass" if the molecule adheres to all criteria.
        - A dictionary with the violated criteria as keys and descriptive messages as values.

    Notes
    -----
    The Ghose Filter is described in:
    A.K Ghose, et al., J. Comb. Chem., 1999, 1, 1, 55-68
    https://doi.org/10.1021/cc9800071
    """
    violation = {}

    # Add explicit hydrogens for filter
    mol = Chem.AddHs(mol)

    logp = calculate_logp_crippen(mol)
    if logp < -0.4 or logp > 5.6:
        violation["LogP"] = f"LogP ({logp}) is outside the acceptable range (-0.4 to +5.6)"

    molecular_weight = calculate_molecular_weight(mol)
    if molecular_weight < 180 or molecular_weight > 480:
        violation["MW"] = f"MW: {molecular_weight} is outside the acceptable range (180-480)"

    molar_refractivity = calculate_molar_refractivity(mol)
    if molar_refractivity < 40 or molar_refractivity > 480:
        violation["MR"] = f"MR: {molar_refractivity} is outside the acceptable range (40-480)"

    number_of_atoms = calculate_number_of_atoms(mol)
    if number_of_atoms < 20 or number_of_atoms > 70:
        violation[
            "num_atoms"
        ] = f"Number of atoms: {number_of_atoms} is outside the acceptable range (20-70)"

    if violation:
        return violation
    else:
        return "Pass"


def calculate_veber(mol: Chem.Mol) -> Union[str, dict[str, str]]:
    """Check if a molecule complies with Veber's Rule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    properties : Union[str, dict[str, str]]
        - "Pass" if the molecule adheres to both criteria.
        - A dictionary with the violated criteria as keys and descriptive messages as values.

    Notes
    -----
    The Veber Filter is described in:
    D.F. Veber, et al. J. Med. Chem 2002, 45, 12, 2615-2623
    https://doi.org/10.1021/jm020017n
    """
    violation = {}

    # Add explicit hydrogens for filter
    mol = Chem.AddHs(mol)

    num_rotatable_bonds = calculate_number_rotatable_bonds(mol)
    if num_rotatable_bonds > 10:
        violation["rotatable_bonds"] = f"Number of Rotatable Bonds: {num_rotatable_bonds} > 10"

    tpsa = calculate_tpsa(mol)
    if tpsa > 140:
        violation["TPSA"] = f"TPSA: {tpsa} > 140 angstrom-squared"

    if violation:
        return violation
    else:
        return "Pass"
