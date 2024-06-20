"""Calculate Lipophilicity parameters."""

from rdkit import Chem


def calculate_logp_crippen(mol):
    """Calculate the LogP of a molecule using Wildman-Crippen Method.

    Parameters
    ----------
         mol : Chem.Mol
              The input rdkit Mol object.

    Returns
    -------
         result : float
              The calculated WLogP of the molecule.
    """
    return Chem.Crippen.MolLogP(mol)
