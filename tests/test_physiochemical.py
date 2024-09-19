"""Test the physiochemical properties."""

from adme_py.physiochemical import (
    calculate_formula,
    calculate_molar_refractivity,
    calculate_molecular_weight,
    calculate_number_aromatic_atoms,
    calculate_number_hbond_acceptors,
    calculate_number_hbond_donors,
    calculate_number_heavy_atoms,
    calculate_number_rotatable_bonds,
    calculate_sp3_carbon_ratio,
    calculate_tpsa,
)


def test_calculate_formula(rdkit_mol):
    """Test the chemical formula calculator."""
    formula = calculate_formula(rdkit_mol)
    expected_formula = "C7H8"
    assert formula == expected_formula


def test_molecular_weight(rdkit_mol):
    """Test the chemical weight calculator."""
    weight = calculate_molecular_weight(rdkit_mol)
    expected_weight = 92.14099999999999
    assert weight == expected_weight


def test_calculate_number_heavy_atoms(rdkit_mol):
    """Test the number of heavy atoms calculator."""
    heavy_atoms = calculate_number_heavy_atoms(rdkit_mol)
    expected_number = 7
    assert expected_number == heavy_atoms


def test_calculate_number_aromatic_atoms(rdkit_mol):
    """Test the number of aromatic atoms calculator."""
    aromatic_atoms = calculate_number_aromatic_atoms(rdkit_mol)
    expected_number = 6
    assert expected_number == aromatic_atoms


def test_calculate_sp3_carbon_ratio(rdkit_mol):
    """Test the calculation of the csp3 ratio."""
    csp3_ratio = calculate_sp3_carbon_ratio(rdkit_mol)
    expected_ratio = 0.14
    assert expected_ratio == round(csp3_ratio, 2)


def test_calculate_number_rotatable_bonds(rdkit_mol):
    """Test the number of rotatble bonds calculator."""
    rotatable_bonds = calculate_number_rotatable_bonds(rdkit_mol)
    expected_number = 0
    assert expected_number == rotatable_bonds


def test_calculate_number_hbond_donors(rdkit_mol):
    """Test the H-bond donor calculator."""
    hbond_donors = calculate_number_hbond_donors(rdkit_mol)
    expected_number = 0
    assert expected_number == hbond_donors


def test_calculate_number_hbond_acceptors(rdkit_mol):
    """Test the H-bond donor calculator."""
    hbond_acceptors = calculate_number_hbond_acceptors(rdkit_mol)
    expected_number = 0
    assert expected_number == hbond_acceptors


def test_calculate_molar_refractivity(rdkit_mol):
    """Test the molar refractivity calculator."""
    molar_refractivity = calculate_molar_refractivity(rdkit_mol)
    expected_value = 31.41
    assert round(expected_value, 0) == round(molar_refractivity, 0)


def test_calculate_tpsa(rdkit_mol):
    """Test the TPSA calculator."""
    calculated_tpsa = calculate_tpsa(rdkit_mol)
    expected_value = 0.0
    assert calculated_tpsa == expected_value
