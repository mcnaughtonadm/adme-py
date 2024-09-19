"""Module for calculating Pharmacokinetic properties."""

from pathlib import Path
from typing import Union

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Polygon
from rdkit import Chem

from adme_py.lipophilicity import calculate_logp_crippen
from adme_py.physiochemical import calculate_molecular_weight, calculate_tpsa

HERE = Path(__file__).parent.resolve()
HIA_COORDS = HERE.joinpath("hia_coords.tsv")
BBB_COORDS = HERE.joinpath("bbb_coords.tsv")


def calculate_all_pharmacokinetics(mol: Chem.Mol) -> dict[str, Union[str, float]]:
    """Calculate various pharmacokinetic properties of a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    dict[str, Union[str, float]]
        A dictionary containing the calculated pharmacokinetic properties:
        - "gastrointestinal_absorption": "High" or "Low" based on the BOILED-Egg model.
        - "blood_brain_barrier_permeant": True if the molecule is predicted to permeate the blood-brain barrier, False otherwise.
        - "skin_permeability_logkp": The predicted skin permeability (LogKp) of the molecule.
    """
    in_hia, in_bbb = predict_bbb_hia(mol)  # type: bool, bool
    hia: str = "High" if in_hia else "Low"
    logkp: float = calculate_skin_permeability(mol)

    properties: dict[str, float] = {
        "gastrointestinal_absorption": hia,
        "blood_brain_barrier_permeant": in_bbb,
        "skin_permeability_logkp": logkp,
    }

    return properties


def predict_bbb_hia(mol: Chem.Mol) -> tuple[bool, bool]:
    """Predict the human intestinal absorption (HIA) and blood-brain barrier (BBB) permeability of a molecule using the BOILED-Egg model.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    tuple[bool, bool]
        A tuple containing two boolean values:
        - The first value indicates whether the molecule is predicted to have high HIA (True) or low HIA (False).
        - The second value indicates whether the molecule is predicted to permeate the BBB (True) or not (False).

    Notes
    -----
    The BOILED-Egg model is described in:
    A. Daina and V. Zoete, ChemMedChem 2016, 11, 1117.
    https:/doi.org/10.1002/cmdc.201600182
    """
    point = _get_tpsa_logp(mol)

    hia_polygon, bbb_polygon = _get_polygons()

    in_hia: bool = hia_polygon.contains_point(point)
    in_bbb: bool = bbb_polygon.contains_point(point)

    return in_hia, in_bbb


def plot_boiled_egg(mols: list[Chem.Mol]) -> None:
    """Plot the BOILED-Egg model for a list of molecules.

    Parameters
    ----------
    mols : List[rdkit.Chem.rdchem.Mol]
        A list of RDKit molecule objects.
    """
    _, ax = plt.subplots()
    hia_polygon, bbb_polygon = _get_polygons()

    ax.set_facecolor("gainsboro")

    ax.add_patch(hia_polygon)
    hia_polygon.set_facecolor("white")
    hia_polygon.set_edgecolor("white")

    ax.add_patch(bbb_polygon)
    bbb_polygon.set_facecolor("yellow")
    bbb_polygon.set_edgecolor("yellow")

    for mol in mols:
        in_hia, in_bbb = predict_bbb_hia(mol)
        prediction = _check_prediction(in_hia, in_bbb)

        color = (
            "black"
            if prediction == "Outside BOILED-Egg"
            else ("blue" if prediction == "HIA" else "red")
        )
        tpsa, logp = _get_tpsa_logp(mol)
        ax.plot(tpsa, logp, "o", color=color)

    ax.set_xlim(-20, 200)
    ax.set_ylim(-3, 8)

    ax.set_xlabel("TPSA")
    ax.set_ylabel("WLOGP")

    # Show the plot
    plt.show()


def calculate_skin_permeability(mol: Chem.Mol) -> float:
    """Calculate the skin permeability (LogKp) of a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    float
        The predicted skin permeability (LogKp) of the molecule.

    Notes
    -----
    This function uses the QSPR model described in:
    R.O Potts and R.H Guy, Pharm Res. 1992 May;9(5):663-9.
    https://doi.org/10.1023/a:1015810312465
    """
    logp = calculate_logp_crippen(mol)
    mw = calculate_molecular_weight(mol)
    return 0.71 * logp - 0.0061 * mw - 6.3


def _get_tpsa_logp(mol: Chem.Mol) -> tuple[float, float]:
    """Calculate the TPSA and LogP of a molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input RDKit molecule object.

    Returns
    -------
    tuple[float, float]
        A tuple containing the TPSA and LogP of the molecule.
    """
    tpsa = calculate_tpsa(mol)
    logp = calculate_logp_crippen(mol)
    return tpsa, logp


def _get_polygons() -> tuple[Polygon, Polygon]:
    """Get the polygons representing the HIA and BBB regions in the BOILED-Egg model.

    Returns
    -------
    tuple[matplotlib.patches.Polygon, matplotlib.patches.Polygon]
        A tuple containing two polygons:
        - The first polygon represents the HIA region.
        - The second polygon represents the BBB region.
    """
    hia_coords = pd.read_csv(HIA_COORDS, sep="\t")
    bbb_coords = pd.read_csv(BBB_COORDS, sep="\t")

    hia_polygon = Polygon(hia_coords, closed=True, label="HIA (Polygon)")
    bbb_polygon = Polygon(bbb_coords, closed=True, label="BBB (Polygon)")

    return hia_polygon, bbb_polygon


def _check_prediction(in_hia: bool, in_bbb: bool) -> str:
    """Check the prediction based on the HIA and BBB status.

    Parameters
    ----------
    in_hia : bool
        True if the molecule is in the HIA region, False otherwise
    in_bbb : bool
        True if the molecule is in the BBB region, False otherwise

    Returns
    -------
    str
        A string representing the prediction:
        - "HIA and BBB" if the molecule is in both regions
        - "HIA" if the molecule is only in the HIA region
        - "BBB" if the molecule is only in the BBB region
        - "Outside BOILED-Egg" if the molecule is in neither region
    """
    if in_hia and in_bbb:
        return "HIA and BBB"
    elif in_hia and not in_bbb:
        return "HIA"
    elif in_bbb and not in_hia:
        return "BBB"
    else:
        return "Outside BOILED-Egg"
