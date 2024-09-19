"""Calculate Lipophilicity parameters."""

from rdkit import Chem


def calculate_all_lipophilicity(mol: Chem.Mol) -> dict[str, float]:
    """Calculate lipophilicity properties of a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    dict[str, float]
        A dictionary containing the calculated lipophilicity properties:
        - "wlogp": Wildman-Crippen LogP value.
    """
    wlogp: float = calculate_logp_crippen(mol)

    properties: dict[str, float] = {
        "wlogp": wlogp,
    }

    return properties


def calculate_logp_crippen(mol) -> float:
    """Calculate the LogP of a molecule using Wildman-Crippen Method.

    Parameters
    ----------
         mol : rdkit.Chem.rdchem.Mol
              The input RDKit molecule object.

    Returns
    -------
         result : float
              The calculated Wildman-Crippen LogP (WLogP) of the molecule.

    Notes
    -----
    The Wildman-Crippen method for LogP calculation is described in:
    S.A. Wildman and G.M. Crippen, JCICS _39_ 868-873 (1999)
    https://doi.org/10.1021/ci990307l
    """
    return Chem.Crippen.MolLogP(mol)
