"""Establishing ADME Class to wrap calculators."""

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
    """

    def __init__(self, identifier):
        self.mol = mol_from_identifier(identifier)
        self.properties = self.calculate()

    def calculate(self):
        """Run calculators to make property predictions."""
        properties = {
            "physiochemical": calculate_all_physiochemical(self.mol),
            "solubility": calculate_all_solubility(self.mol),
            "lipophilicity": calculate_all_lipophilicity(self.mol),
            "pharmacokinetics": calculate_all_pharmacokinetics(self.mol),
            "druglikeness": calculate_all_druglikeness(self.mol),
            "medicinal": calculate_all_medicinal(self.mol),
        }

        return properties
