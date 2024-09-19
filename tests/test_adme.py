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


def test_adme_molar_refractivity(test_adme):
    """Test physiochemical property molecular refractivity from ADME class."""
    molar_refractivity = test_adme.properties["physiochemical"]["molar_refractivity"]
    expected_value = 31.41
    assert round(expected_value, 0) == round(molar_refractivity, 0)


def test_adme_tpsa(test_adme):
    """Test physiochemical property TPSA from ADME class."""
    calculated_tpsa = test_adme.properties["physiochemical"]["tpsa"]
    expected_value = 0.0
    assert calculated_tpsa == expected_value


def test_adme_log_s_esol(test_adme):
    """Test solubility property LogS (ESOL) from ADME class."""
    calculated_log_s_esol = test_adme.properties["solubility"]["log_s_esol"]
    expected_value = -2.1932094391812655
    assert calculated_log_s_esol == expected_value


def test_adme_solubility_esol(test_adme):
    """Test solubility property solubility (ESOL) from ADME class."""
    calculated_solubility_esol = test_adme.properties["solubility"]["solubility_esol"]
    expected_value = 10**-2.1932094391812655
    assert calculated_solubility_esol == expected_value


def test_adme_esol_class(test_adme):
    """Test the solubility property solubility class (ESOL) from ADME class."""
    calculated_esol_class = test_adme.properties["solubility"]["class_esol"]
    expected_class = "Soluble"
    assert calculated_esol_class == expected_class


def test_adme_lipophilicity_wlopg(test_adme):
    """Test the lipophilicity property WLOGP from ADME class."""
    calculated_wlogp = test_adme.properties["lipophilicity"]["wlogp"]
    expected_value = 1.99502
    assert calculated_wlogp == expected_value


def test_bbb_absorption(test_adme):
    """Test the pharmacokinetic property Blood Brain Barrier Permeance from ADME class."""
    calculated_bbb = test_adme.properties["pharmacokinetics"]["blood_brain_barrier_permeant"]
    expected_value = False
    assert calculated_bbb == expected_value


def test_gastrointestinal_absorption(test_adme):
    """Test the pharmacokinetic property gastrointestinal absorption from ADME class."""
    calculated_hia = test_adme.properties["pharmacokinetics"]["gastrointestinal_absorption"]
    expected_value = "Low"
    assert calculated_hia == expected_value


def test_skin_permeability_logkp(test_adme):
    """Test the pharmacokinetic property Skin Permeability (LogKp) from ADME class."""
    calculated_logkp = test_adme.properties["pharmacokinetics"]["skin_permeability_logkp"]
    expected_value = -5.4455959
    assert calculated_logkp == expected_value


def test_druglikeness_lipinski(test_adme):
    """Test the druglikeness filter for lipinski rule of 5."""
    calculated_lipinski = test_adme.properties["druglikeness"]["lipinski"]
    expected_value = "Pass"

    assert calculated_lipinski == expected_value


def test_druglikeness_ghose(test_adme):
    """Test the druglikeness filter for Ghose."""
    calculated_ghose = test_adme.properties["druglikeness"]["ghose"]
    expected_value = {
        "MW": "MW: 92.14099999999995 is outside the acceptable range (180-480)",
        "MR": "MR: 31.17899999999999 is outside the acceptable range (40-480)",
        "num_atoms": "Number of atoms: 15 is outside the acceptable range (20-70)",
    }

    assert calculated_ghose == expected_value


def test_druglikeness_veber(test_adme):
    """Test the druglikeness filter for Veber."""
    calculated_veber = test_adme.properties["druglikeness"]["veber"]
    expected_value = "Pass"

    assert calculated_veber == expected_value


def test_medicinal_pains(test_adme):
    """Test the medicinal PAINS filter."""
    assert not test_adme.properties["medicinal"]["pains"]


def test_medicinal_brenk(test_adme):
    """Test the medicinal Brenk filter."""
    assert not test_adme.properties["medicinal"]["brenk"]


def test_medicinal_zinc(test_adme):
    """Test the medicinal Zinc filter."""
    assert not test_adme.properties["medicinal"]["zinc"]


def test_medicinal_leadlikeness(test_adme):
    """Test the medicinal leadlikeness filter."""
    result = test_adme.properties["medicinal"]["leadlikeness"]
    expected_result = {"MW": "MW: 92.14099999999999 is outside the acceptable range (250-350)"}

    assert expected_result == result


def test_medicinal_synthetic_accessibility(test_adme):
    """Test the medicinal synthetic accesssibility calculator."""
    result = test_adme.properties["medicinal"]["synthetic_accessibility"]
    expected_result = 1.0

    assert expected_result == result
