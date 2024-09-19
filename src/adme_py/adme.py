"""Establishing ADME Class to wrap calculators."""

from typing import Union

from adme_py.druglikeness import calculate_all_druglikeness
from adme_py.input_handlers import mol_from_identifier
from adme_py.lipophilicity import calculate_all_lipophilicity
from adme_py.medicinal import calculate_all_medicinal
from adme_py.pharmacokinetics import calculate_all_pharmacokinetics
from adme_py.physiochemical import calculate_all_physiochemical
from adme_py.solubility import calculate_all_solubility


class ADME:
    """ADME Class object for aggregating ADME Properties.

    Parameters
    ----------
    identifier : str
        The identifier of the molecule you'd like to analyze.

    Attributes
    ----------
    mol : rdkit.Chem.rdchem.Mol
        The input molecule as an RDKit molecule object.
    properties : dict[str, dict[str, Union[str, float, int, bool, dict[str, str]]]]
        A dictionary containing calculated ADME properties, where the keys are property categories
        (e.g., "physiochemical", "solubility") and the values are dictionaries of property names
        and their values.
    """

    def __init__(self, identifier) -> None:
        self.mol = mol_from_identifier(identifier)
        self.properties = self.calculate()

    def calculate(self) -> dict[str, dict[str, Union[str, float, int, bool, dict[str, str]]]]:
        """Run calculators to make property predictions.

        Returns
        -------
        dict[str, dict[str, Union[str, float, int, bool, dict[str, str]]]]
            A dictionary containing calculated ADME properties.
        """
        properties = {
            "physiochemical": calculate_all_physiochemical(self.mol),
            "solubility": calculate_all_solubility(self.mol),
            "lipophilicity": calculate_all_lipophilicity(self.mol),
            "pharmacokinetics": calculate_all_pharmacokinetics(self.mol),
            "druglikeness": calculate_all_druglikeness(self.mol),
            "medicinal": calculate_all_medicinal(self.mol),
        }

        return properties
