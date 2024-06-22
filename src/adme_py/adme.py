"""Establishing ADME Class to wrap calculators."""

from adme_py.input_handlers import mol_from_identifier
from adme_py.physiochemical import calculate_all_physiochemical


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
        properties = {}

        properties["physiochemical"] = calculate_all_physiochemical(self.mol)

        return properties
