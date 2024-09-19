"""Module file for calculating physiochemical properties."""

from typing import Union

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


def calculate_all_physiochemical(mol: Chem.Mol) -> dict[str, Union[str, float, int]]:
    """Calculate various physiochemical properties of a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    dict[str, Union[str, float, int]]
        A dictionary containing the calculated physiochemical properties:
        - "formula": The chemical formula of the molecule (str).
        - "molecular_weight": The molecular weight of the molecule (float).
        - "num_heavy_atoms": The number of heavy atoms in the molecule (int).
        - "num_aromatic_atoms": The number of aromatic atoms in the molecule (int).
        - "sp3_carbon_ratio": The ratio of sp3-hybridized carbons to total carbons (float).
        - "num_rotatable_bonds": The number of rotatable bonds in the molecule (int).
        - "num_h_donors": The number of hydrogen bond donors in the molecule (int).
        - "num_h_acceptors": The number of hydrogen bond acceptors in the molecule (int).
        - "molar_refractivity": The molar refractivity calculated using the Crippen algorithm (float).
        - "tpsa": The topological polar surface area (TPSA) of the molecule (float).
    """
    properties: dict[str, Union[str, float, int]] = {
        "formula": calculate_formula(mol),
        "molecular_weight": calculate_molecular_weight(mol),
        "num_heavy_atoms": calculate_number_heavy_atoms(mol),
        "num_aromatic_atoms": calculate_number_aromatic_atoms(mol),
        "sp3_carbon_ratio": calculate_sp3_carbon_ratio(mol),
        "num_rotatable_bonds": calculate_number_rotatable_bonds(mol),
        "num_h_donors": calculate_number_hbond_donors(mol),
        "num_h_acceptors": calculate_number_hbond_acceptors(mol),
        "molar_refractivity": calculate_molar_refractivity(mol),
        "tpsa": calculate_tpsa(mol),
        "num_atoms": calculate_number_of_atoms(mol),
    }

    return properties


def calculate_formula(mol: Chem.Mol) -> str:
    """Calculate the chemical formula of a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
              The input RDKit molecule object.

    Returns
    -------
    str
        The chemical formula of the molecule.
    """
    return rdMolDescriptors.CalcMolFormula(mol)


def calculate_molecular_weight(mol: Chem.Mol) -> float:
    """Calculate the molecular weight of a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    float
        The molecular weight of the molecule.
    """
    return Descriptors.MolWt(mol)


def calculate_number_heavy_atoms(mol: Chem.Mol) -> int:
    """Calculate the number of heavy atoms (non-hydrogen atoms) in a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    int
        The number of heavy atoms in the molecule.
    """
    return rdMolDescriptors.CalcNumHeavyAtoms(mol)


def calculate_number_aromatic_atoms(mol: Chem.Mol) -> int:
    """Calculate the number of aromatic atoms in a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    int
        The number of aromatic atoms in the molecule.
    """
    return len(mol.GetAromaticAtoms())


def calculate_sp3_carbon_ratio(mol: Chem.Mol) -> float:
    """Calculate the ratio of sp3-hybridized carbons to total carbons in a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    float
        The calculated ratio of sp3-hybridized carbons to total carbons.
    """
    return Chem.Lipinski.FractionCSP3(mol)


def calculate_number_rotatable_bonds(mol: Chem.Mol) -> int:
    """Calculate the number of rotatable bonds in a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    int
        The number of rotatable bonds in the molecule.
    """
    return rdMolDescriptors.CalcNumRotatableBonds(mol)


def calculate_number_hbond_donors(mol: Chem.Mol) -> int:
    """Calculate the number of hydrogen bond donors in a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    int
        The number of hydrogen bond donors in the molecule.
    """
    return Chem.Lipinski.NumHDonors(mol)


def calculate_number_hbond_acceptors(mol: Chem.Mol) -> int:
    """Calculate the number of hydrogen bond acceptors in a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    int
        The number of hydrogen bond acceptors in the molecule.
    """
    return Chem.Lipinski.NumHAcceptors(mol)


def calculate_tpsa(mol: Chem.Mol) -> float:
    """Calculate the topological polar surface area (TPSA) of a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    float
        The calculated TPSA of the molecule.

    Notes
    -----
    The method for calculating the topological polar surface area is described in:
    P. Ertl, B. Rohde, and P. Selzer J. Med. Chem. 43:3714-7, (2000)
    https://doi.org/10.1021/jm000942e
    """
    return rdMolDescriptors.CalcTPSA(mol)


def calculate_molar_refractivity(mol: Chem.Mol) -> float:
    """Calculate the molar refractivity using the Wildman-Crippen algorithm.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    float
        The calculated molar refractivity of the molecule.

    Notes
    -----
    The method for calculating Molar Refractivity is described in:
    S.A. Wildman and G.M. Crippen, JCICS _39_ 868-873 (1999)
    https://doi.org/10.1021/ci990307l
    """
    return Chem.Crippen.MolMR(mol)


def calculate_number_of_atoms(mol: Chem.Mol) -> int:
    """Calculate the total number of atoms in a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    int
        The number of atoms in the molecule.
    """
    return mol.GetNumAtoms()
