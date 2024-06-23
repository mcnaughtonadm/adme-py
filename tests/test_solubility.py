"""Test the water solubility properties."""

from adme_py.solubility import calculate_esol


def test_calculate_esol(rdkit_mol):
    """Test the calculation of esol."""
    esol = calculate_esol(rdkit_mol)
    actual_esol = -2.1932094391812655

    assert esol == actual_esol
    assert 10**esol == 10**actual_esol
