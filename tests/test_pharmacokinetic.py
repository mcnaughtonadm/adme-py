"""Test the pharmacokinetic property calculators."""

from adme_py.pharmacokinetics import calculate_skin_permeability, predict_bbb_hia


def test_calculate_skin_permeability(rdkit_mol):
    """Test the calculation of LogKp."""
    calculated_value = calculate_skin_permeability(rdkit_mol)
    expected_value = -5.4455959
    assert calculated_value == expected_value


def test_bbb_hia(rdkit_mol):
    """Test the Blood Brain Barrier and Gastrointestinal Calculators."""
    in_hia, in_bbb = predict_bbb_hia(rdkit_mol)

    assert not in_hia
    assert not in_bbb
