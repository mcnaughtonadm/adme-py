"""Module for calculating Pharmacokinetic properties."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Polygon
from rdkit import Chem

from adme_py.lipophilicity import calculate_logp_crippen
from adme_py.physiochemical import calculate_molecular_weight, calculate_tpsa

HERE = Path(__file__).parent.resolve()
HIA_COORDS = HERE.joinpath("hia_coords.tsv")
BBB_COORDS = HERE.joinpath("bbb_coords.tsv")


def calculate_all_pharmacokinetics(mol: Chem.Mol):
    """Calculate all the pharmacokinetic properties for a given molecule.

    Parameters
    ----------
        mol : Chem.Mol
             The input rdkit Mol object
    """
    in_hia, in_bbb = predict_bbb_hia(mol)
    hia = "High" if in_hia else "Low"
    logkp = calculate_skin_permeability(mol)

    properties: dict[str, float] = {
        "gastrointestinal_absorption": hia,
        "blood_brain_barrier_permeant": in_bbb,
        "skin_permeability_logkp": logkp,
    }

    return properties


def predict_bbb_hia(mol: Chem.Mol) -> tuple[bool, bool]:
    """Predict whether a molecule lands in the BOILED-Egg Model."""
    point = _get_tpsa_logp(mol)

    hia_polygon, bbb_polygon = _get_polygons()

    in_hia: bool = hia_polygon.contains_point(point)
    in_bbb: bool = bbb_polygon.contains_point(point)

    return in_hia, in_bbb


def plot_boiled_egg(mols: list[Chem.Mol]):
    """Plot the results of the BOILED-Egg model for list of molecules."""
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


def calculate_skin_permeability(mol: Chem.Mol):
    """Calculate the skin permeability (LogKp)of a molecule.

    Leverages the QSPR model present in Potts and Guy 1992.
    """
    logp = calculate_logp_crippen(mol)
    mw = calculate_molecular_weight(mol)
    return 0.71 * logp - 0.0061 * mw - 6.3


def _get_tpsa_logp(mol: Chem.Mol):
    tpsa = calculate_tpsa(mol)
    logp = calculate_logp_crippen(mol)
    return tpsa, logp


def _get_polygons():
    hia_coords = pd.read_csv(HIA_COORDS, sep="\t")
    bbb_coords = pd.read_csv(BBB_COORDS, sep="\t")

    hia_polygon = Polygon(hia_coords, closed=True, label="HIA (Polygon)")
    bbb_polygon = Polygon(bbb_coords, closed=True, label="BBB (Polygon)")

    return hia_polygon, bbb_polygon


def _check_prediction(in_hia: bool, in_bbb: bool):
    if in_hia and in_bbb:
        return "HIA and BBB"
    elif in_hia and not in_bbb:
        return "HIA"
    elif in_bbb and not in_hia:
        return "BBB"
    else:
        return "Outside BOILED-Egg"
