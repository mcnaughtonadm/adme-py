"""Test the pharmacokinetic property calculators."""

from adme_py.druglikeness import calculate_ghose, calculate_lipinksi, calculate_veber
from rdkit import Chem


def test_calculate_lipinski_pass(rdkit_mol):
    """Test that a molecule passing the Lipinksi filter, returns Pass."""
    result = calculate_lipinksi(rdkit_mol)
    assert result == "Pass", f"Expected 'Pass' but got {result}"


def test_calculate_lipinski_fail():
    """Test that a molecule violating the Lipinksi filter, returns the correct violations."""
    smiles = "CC1(CCC(=C(C1)C2=CC=C(C=C2)Cl)CN3CCN(CC3)C4=CC(=C(C=C4)C(=O)NS(=O)(=O)C5=CC(=C(C=C5)NCC6CCOCC6)[N+](=O)[O-])OC7=CN=C8C(=C7)C=CN8)C"
    mol = Chem.MolFromSmiles(smiles)
    result = calculate_lipinksi(mol)

    expected_violations = {
        "MW": "MW: 868.4570000000017 > 500 Dalton",
        "LogP": "LogP: 8.659900000000006 >5",
        "hbond_acceptors": "hbond acceptors: 11 >10",
    }

    assert isinstance(result, dict)
    assert result == expected_violations


def test_calculate_ghose_pass():
    """Test that a molecule passing the Ghose filter, returns Pass."""
    smiles = "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"
    mol = Chem.MolFromSmiles(smiles)
    result = calculate_ghose(mol)

    assert result == "Pass", f"Expected 'Pass' but got {result}"


def test_calculate_ghose_fail():
    """Test that a molecule violating the Ghose filter, returns the correct violations."""
    smiles = "CC1(CCC(=C(C1)C2=CC=C(C=C2)Cl)CN3CCN(CC3)C4=CC(=C(C=C4)C(=O)NS(=O)(=O)C5=CC(=C(C=C5)NCC6CCOCC6)[N+](=O)[O-])OC7=CN=C8C(=C7)C=CN8)C"
    mol = Chem.MolFromSmiles(smiles)
    result = calculate_ghose(mol)

    expected_violations = {
        "LogP": "LogP (8.659900000000006) is outside the acceptable range (-0.4 to +5.6)",
        "MW": "MW: 868.4570000000017 is outside the acceptable range (180-480)",
        "num_atoms": "Number of atoms: 111 is outside the acceptable range (20-70)",
    }

    assert isinstance(expected_violations, dict)
    assert result == expected_violations


def test_calculate_veber_pass(rdkit_mol):
    """Test that a molecule passing the Veber filter, returns Pass."""
    result = calculate_veber(rdkit_mol)

    assert result == "Pass", f"Expected 'Pass' but got {result}"


def test_calculate_veber_fail():
    """Test that a molecule violating the Veber filter, returns the correct violations."""
    smiles = "CC1(CCC(=C(C1)C2=CC=C(C=C2)Cl)CN3CCN(CC3)C4=CC(=C(C=C4)C(=O)NS(=O)(=O)C5=CC(=C(C=C5)NCC6CCOCC6)[N+](=O)[O-])OC7=CN=C8C(=C7)C=CN8)C"
    mol = Chem.MolFromSmiles(smiles)
    result = calculate_ghose(mol)

    expected_violations = {
        "LogP": "LogP (8.659900000000006) is outside the acceptable range (-0.4 to +5.6)",
        "MW": "MW: 868.4570000000017 is outside the acceptable range (180-480)",
        "num_atoms": "Number of atoms: 111 is outside the acceptable range (20-70)",
    }

    assert isinstance(result, dict)
    assert result == expected_violations
