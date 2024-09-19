"""Tests for methods in medicinal.py file."""

from adme_py.medicinal import (
    calculate_brenk,
    calculate_leadlikeness,
    calculate_pains,
    calculate_synthetic_accessibility,
    calculate_zinc,
)
from rdkit import Chem


def test_calculate_pains_pass(rdkit_mol):
    """Test a passing case of calculate_pains()."""
    assert not calculate_pains(rdkit_mol)


def test_calculate_pains_fail():
    """Test a failing case of calculate_pains()."""
    pains_mol = Chem.MolFromSmiles("CC(=O)NCC1CN(C(=O)O1)C2=CC(=C(C=C2)N3CCOCC3)F")
    assert calculate_pains(pains_mol)


def test_calculate_brenk_pass(rdkit_mol):
    """Test a passing case of calculate_brenk()."""
    assert not calculate_brenk(rdkit_mol)


def test_calculate_brenk_fail():
    """Test a failing case of calculate_brenk()."""
    brenk_mol = Chem.MolFromSmiles("OCOCC(COc1ccc(c(c1C(=O)Oc1ccc(cc1N)N1CCNCC1)N)OC)O")
    assert calculate_brenk(brenk_mol)


def test_calculate_zinc_pass(rdkit_mol):
    """Test a passing case of calculate_zinc()."""
    assert not calculate_zinc(rdkit_mol)


def test_calculate_zinc_fail():
    """Test a failing case of calculate_zinc()."""
    zinc_mol = Chem.MolFromSmiles(
        "CC1=C2C(C(=O)C3(C(CC4C(C3C(C(C2(C)C)(CC1OC(=O)C(C(C5=CC=CC=C5)NC(=O)C6=CC=CC=C6)O)O)OC(=O)C7=CC=CC=C7)(CO4)OC(=O)C)O)C)OC(=O)C"
    )
    assert calculate_zinc(zinc_mol)


def test_calculate_leadlikeness_pass(rdkit_mol):
    """Test a passing case of calculate_leadlikeness()."""
    mol = Chem.MolFromSmiles("COc1ccc2c(c1)nc([nH]2)S(=O)Cc1ncc(c(c1C)OC)C")
    result = calculate_leadlikeness(mol)
    expected_result = "Yes"
    assert result == expected_result


def test_calculate_leadlikeness_fail(rdkit_mol):
    """Test a failing case of calculate_leadlikeness()."""
    result = calculate_leadlikeness(rdkit_mol)
    expected_violations = {
        "MW": "MW: 92.14099999999999 is outside the acceptable range (250-350)",
    }

    assert isinstance(result, dict)
    assert result == expected_violations


def test_calculate_synthetic_accessibility(rdkit_mol):
    """Test the result of calculate_synthetic_accessibility()."""
    result = calculate_synthetic_accessibility(rdkit_mol)
    expected_result = 1.0

    assert result == expected_result
