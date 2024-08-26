"""Calculate Lipophilicity parameters."""

from rdkit import Chem


def calculate_all_lipophilicity(mol: Chem.Mol) -> dict[str, float]:
    """Calculate the all the lipophilicity properties of a given molecule.

    Parameters
    ----------
        mol : Chem.Mol
             The input rdkit Mol object
    """
    wlogp: float = calculate_logp_crippen(mol)

    properties: dict[str, float] = {
        "wlogp": wlogp,
    }

    return properties


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
