"""Module file for calculating physiochemical properties."""

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


def calculate_all_physiochemical(mol: Chem.Mol) -> dict[str, str | float | int]:
    """Calculate all physiochemical properties of a given molecule.

    Parameters
    ----------
         mol : Chem.Mol
              The input rdkit Mol object
    """
    properties: dict[str, str | float | int] = {
        "formula": calculate_formula(mol),
        "molecular_weight": calculate_molecular_weight(mol),
        "num_heavy_atoms": calculate_number_heavy_atoms(mol),
        "num_aromatic_atoms": calculate_number_aromatic_atoms(mol),
        "sp3_carbon_ratio": calculate_sp3_carbon_ratio(mol),
        "num_rotatable_bonds": calculate_number_rotatable_bonds(mol),
        "num_h_donors": calculate_number_hbond_donors(mol),
        "num_h_acceptors": calculate_number_hbond_acceptors(mol),
        "molecular_refractivity": calculate_molecular_refractivity(mol),
        "tpsa": calculate_tpsa(mol),
    }

    return properties


def calculate_formula(mol: Chem.Mol) -> str:
    """Calculate the formula of a given molecule.

    Parameters
    ----------
         mol : Chem.Mol
              The input rdkit Mol object

    Returns
    -------
         formula : str
              The chemical formula of the molecule
    """
    return rdMolDescriptors.CalcMolFormula(mol)


def calculate_molecular_weight(mol: Chem.Mol) -> float:
    """Calculate the molecular weight of a given molecule.

    Parameters
    ----------
         mol : Chem.Mol
              The input rdkit Mol object.

    Returns
    -------
         result : float
              The molecular weight of the molecule.
    """
    return Descriptors.MolWt(mol)


def calculate_number_heavy_atoms(mol: Chem.Mol) -> int:
    """Calculate the number of heavy atoms in a given molecule.

    Parameters
    ----------
         mol : Chem.Mol
              The input rdkit Mol object.

    Returns
    -------
         result : int
              The number of heavy atoms in the molecule.
    """
    return rdMolDescriptors.CalcNumHeavyAtoms(mol)


def calculate_number_aromatic_atoms(mol: Chem.Mol) -> int:
    """Calculate the number of aromatic atoms in a given molecule.

    Parameters
    ----------
         mol : Chem.Mol
              The input rdkit Mol object.

    Returns
    -------
         result : int
              The number of aromatic atoms in the molecule.
    """
    return len(mol.GetAromaticAtoms())


def calculate_sp3_carbon_ratio(mol: Chem.Mol) -> float:
    """Calculate the ratio of sp3-hybridized carbons to total carbons.

    Parameters
    ----------
         mol : Chem.Mol
              The input rdkit Mol object.

    Returns
    -------
         result : float
              The calculated ratio of sp3-hydrizized carbons to total carbons.
    """
    return Chem.Lipinski.FractionCSP3(mol)


def calculate_number_rotatable_bonds(mol: Chem.Mol) -> int:
    """Calculate the number of rotatable bonds.

    Parameters
    ----------
         mol : Chem.Mol
              The input rdkit Mol object.

    Returns
    -------
         result : float
              The number of rotatable bonds in the molecule.
    """
    return rdMolDescriptors.CalcNumRotatableBonds(mol)


def calculate_number_hbond_donors(mol: Chem.Mol) -> int:
    """Calculate the number of H-bond donors.

    Parameters
    ----------
         mol : Chem.Mol
              The input rdkit Mol object.

    Returns
    -------
         result : int
              The number of h-bond donors in the molecule.
    """
    return Chem.Lipinski.NumHDonors(mol)


def calculate_number_hbond_acceptors(mol: Chem.Mol) -> int:
    """Calculate the number of H-bond acceptors.

    Parameters
    ----------
         mol : Chem.Mol
              The input rdkit Mol object.

    Returns
    -------
         result : int
              The number of h-bond acceptors in the molecule.
    """
    return Chem.Lipinski.NumHAcceptors(mol)


def calculate_molecular_refractivity(mol: Chem.Mol) -> float:
    """Calculate the Molecular Refractivity using the Crippen algorithm.

    Parameters
    ----------
         mol : Chem.Mol
              The input rdkit Mol object.

    Returns
    -------
         result : float
              The molecules calculated Molecular Refractivity.
    """
    return Chem.Crippen.MolMR(mol)


def calculate_tpsa(mol: Chem.Mol) -> float:
    """Calculate the Topological Polar Surface Area.

    Parameters
    ----------
         mol : Chem.Mol
              The input rdkit Mol object.

    Returns
    -------
         result : float
              The molecules calculated Topological Polar Surface Area.
    """
    return rdMolDescriptors.CalcTPSA(mol)


def calculate_molar_refractivity(mol: Chem.Mol):
    """Calculate the Molar Refractivity using the Wildman-Crippen algorithm.

    Parameters
    ----------
         mol : Chem.Mol
              The input rdkit Mol object.

    Returns
    -------
         result : float
              The molecules calculated Molar Refractivity.
    """
    return Chem.Crippen.MolMR(mol)


def calculate_number_of_atoms(mol: Chem.Mol):
    """Calculate the total number of atoms.

    Parameters
    ----------
         mol : Chem.Mol
              The input rdkit Mol object.

    Returns
    -------
         result : int
              The number of atoms in the molecule.
    """
    return mol.GetNumAtoms()
