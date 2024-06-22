"""Fixture file for pytest."""

import pytest
from adme_py import ADME


@pytest.fixture
def rdkit_mol():
    """Fixture to create an RDKit molecule from a SMILES string."""
    return ADME("Cc1ccccc1").mol


@pytest.fixture
def test_adme():
    """Fixture to create an ADME class for testing."""
    return ADME("Cc1ccccc1")
