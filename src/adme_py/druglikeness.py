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


def calculate_all_druglikeness(mol: Chem.Mol):
    """Calculate the all the lipophilicity properties of a given molecule.

    Parameters
    ----------
        mol : Chem.Mol
             The input rdkit Mol object
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


def calculate_lipinksi(mol: Chem.Mol) -> Union[str, dict]:
    """Calculate whether the molecule passes the Lipinski Rule of 5."""
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


def calculate_ghose(mol: Chem.Mol) -> Union[str, dict]:
    """Calculate whether the molecule passes the Ghose Filter."""
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


def calculate_veber(mol: Chem.Mol) -> Union[str, dict]:
    """Calculate whether the molecule passes the Veber's Rule."""
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
