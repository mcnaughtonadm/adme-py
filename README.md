
![ADME-py-logo_large](https://github.com/mcnaughtonadm/adme-py/assets/47539433/a2b53c05-8d36-4798-ad51-111823527217)

# adme-py 💊
[![Rye](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/rye/main/artwork/badge.json)](https://rye.astral.sh/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

`adme-py` is a Python package that leverages open-source tools such as `rdkit` to generate *Absorption*, *Distribution*, *Metabolism*, and *Excretion* (and sometimes including *Toxicity*) summaries for molecules. 

It seeks to mimic the outputs provided by commonly-used tools such as [SwissADME](http://www.swissadme.ch/index.php) and [others](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9864198/) but is not restricted to a web application and can be imported or used through the command-line. This allows for nicer integration between other Python scripts.

# Installation 🪛 

To install `adme-py` run the following:

```bash
pip install adme-py
```

or install from github directly with:

```bash
pip install git+https://github.com/mcnaughtonadm/adme-py.git
```

# How to Run 🏃

To run `adme-py` in it's most basic form, you can do the following:

```python
from adme_py import ADME

summary = ADME("Cc1ccccc1").calculate()

summary.properties["physiochemical"]["tpsa"]
```

# License 📜

The code in this package is licensed under the [MIT](LICENSE) License.

