"""Fixture file for pytest."""

import pytest
from rdkit import Chem


@pytest.fixture
def rdkit_mol():
    """Fixture to create an RDKit molecule from a SMILES string."""
    smiles = "Cc1ccccc1"
    return Chem.MolFromSmiles(smiles)
