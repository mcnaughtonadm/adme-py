"""Tests for input_handlers.py scripts."""

import pytest
from adme_py.input_handlers import mol_from_identifier
from rdkit import Chem


# Valid Identifier Test-Case
@pytest.mark.parametrize(
    "identifier",
    [
        "Cc1ccccc1",
        "InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3",
        "[#8H]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1",
    ],
)
def test_valid_identifiers(identifier):
    """Tests valid identifiers on the mol_from_identifier method."""
    mol = mol_from_identifier(identifier)

    assert isinstance(mol, Chem.Mol)


# Invalid Identifier Test-Case
@pytest.mark.parametrize(
    "identifier, expected_error",
    [
        ("invalid_smiles", "Unsupported string identifier format."),
        (0, "Input identifier must be a string."),
        ({}, "Input identifier must be a string."),
    ],
)
def test_invalid_identifiers(identifier, expected_error):
    """Tests invalid identifiers on the mol_from_identifier method."""
    with pytest.raises(ValueError, match=expected_error):
        mol_from_identifier(identifier)
