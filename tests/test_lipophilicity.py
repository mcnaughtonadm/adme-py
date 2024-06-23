"""Test lipophilicity calculators."""

from adme_py.lipophilicity import (
    calculate_logp_crippen,
)


def test_calculate_logp_crippen(rdkit_mol):
    """Test the calculation of the WLogP."""
    mlogp_value = calculate_logp_crippen(rdkit_mol)
    expected_value = 1.99502

    assert expected_value == mlogp_value
