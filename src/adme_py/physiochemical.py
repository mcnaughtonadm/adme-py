"""Module file for calculating physiochemical properties."""

from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors


def calculate_formula(mol: Chem.Mol) -> str:
    """Calculate the formula of a given molecule.

    Arguments
    ---------
         mol (Chem.Mol): The input rdkit Mol object

    Returns
    -------
         formula (str):  The chemical formula of the molecule

    """
    return rdMolDescriptors.CalcMolFormula(mol)


def calculate_molecular_weight(mol: Chem.Mol) -> float:
    """Calculate the molecular weight of a given molecule."""
    return Descriptors.MolWt(mol)


def calculate_number_heavy_atoms(mol: Chem.Mol) -> int:
    """Calculate the number of heavy atoms in a given molecule."""
    return rdMolDescriptors.CalcNumHeavyAtoms(mol)


def calculate_number_aromatic_atoms(mol: Chem.Mol) -> int:
    """Calculate the number of aromatic atoms in a given molecule."""
    return len(mol.GetAromaticAtoms())
