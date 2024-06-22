"""Tests for ADME class methods."""

import pytest
from adme_py import ADME


# Test-Case Initialization
@pytest.mark.parametrize(
    "identifier",
    [
        "Cc1ccccc1",
        "InChI=1S/C7H8/c1-7-5-3-2-4-6-7/h2-6H,1H3",
        "[#8H]-[#6]1:[#6]:[#6]:[#6]:[#6]:[#6]:1",
    ],
)
def test_adme_initialization(identifier):
    """Test ADME class initialization."""
    adme = ADME(identifier)
    assert isinstance(adme, ADME)


# Test-Case Physiochemical
def test_adme_physiochemical_formula(test_adme):
    """Test physiochemical property formula prediction from ADME class."""
    formula = test_adme.properties["physiochemical"]["formula"]
    expected_formula = "C7H8"

    assert expected_formula == formula


def test_adme_physiochemical_molecular_weight(test_adme):
    """Test physiochemical property MolWt prediction from ADME class."""
    weight = test_adme.properties["physiochemical"]["molecular_weight"]
    expected_weight = 92.14099999999999

    assert expected_weight == weight


def test_adme_physiochemical_number_heavy_atoms(test_adme):
    """Test physiochemical property number of heavy atoms from ADME class."""
    heavy_atoms = test_adme.properties["physiochemical"]["num_heavy_atoms"]
    expected_number = 7
    assert expected_number == heavy_atoms


def test_adme_number_aromatic_atoms(test_adme):
    """Test physiochemical property number of aromatic atoms from ADME class."""
    aromatic_atoms = test_adme.properties["physiochemical"]["num_aromatic_atoms"]
    expected_number = 6
    assert expected_number == aromatic_atoms


def test_adme_sp3_carbon_ratio(test_adme):
    """Test physiochemical property calculation of the csp3 ratio from ADME class."""
    csp3_ratio = test_adme.properties["physiochemical"]["sp3_carbon_ratio"]
    expected_ratio = 0.14
    assert expected_ratio == round(csp3_ratio, 2)


def test_adme_number_rotatable_bonds(test_adme):
    """Test physiochemical property number of rotatble bonds from ADME class."""
    rotatable_bonds = test_adme.properties["physiochemical"]["num_rotatable_bonds"]
    expected_number = 0
    assert expected_number == rotatable_bonds


def test_adme_number_hbond_donors(test_adme):
    """Test physiochemical property number of H-bond donor from ADME class."""
    hbond_donors = test_adme.properties["physiochemical"]["num_h_donors"]
    expected_number = 0
    assert expected_number == hbond_donors


def test_adme_number_hbond_acceptors(test_adme):
    """Test physiochemical property number of H-bond donor from ADME class."""
    hbond_acceptors = test_adme.properties["physiochemical"]["num_h_acceptors"]
    expected_number = 0
    assert expected_number == hbond_acceptors


def test_adme_molecular_refractivity(test_adme):
    """Test physiochemical property molecular refractivity from ADME class."""
    molecular_refractivity = test_adme.properties["physiochemical"]["molecular_refractivity"]
    expected_value = 31.41
    assert round(expected_value, 0) == round(molecular_refractivity, 0)


def test_adme_tpsa(test_adme):
    """Test physiochemical property TPSA from ADME class."""
    calculated_tpsa = test_adme.properties["physiochemical"]["tpsa"]
    expected_value = 0.0
    assert calculated_tpsa == expected_value
