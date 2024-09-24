adme-py Documentation
=====================

Welcome to the documentation for **adme-py**, a Python package that leverages open-source tools such as `rdkit` to generate *Absorption*, *Distribution*, *Metabolism*, and *Excretion* summaries for molecules.

Project Description
-------------------

**adme-py** seeks to mimic the outputs provided by commonly-used tools such as `SwissADME` and others but is not restricted to a web application and can be imported or used through the command-line. This allows for nicer integration between other Python scripts.

Installation
------------
To install **adme-py**, simply run:

.. code-block:: bash

   pip install adme-py

or install it directly from GitHub:

.. code-block:: bash

   pip install git+https://github.com/mcnaughtonadm/adme-py.git

Usage
-----
To run `adme-py` in its most basic form, you can do the following:

.. code-block:: python

   from adme_py import ADME

   summary = ADME("Cc1ccccc1").calculate()

   summary.properties

License
-------
The code in this package is licensed under the MIT License.

Contents
--------
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules

Indices and tables
------------------
* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

FAQ
---
* **How do I install adme-py?**
  See the `Installation`_ section.

* **Where can I find examples?**
  Check the `Usage`_ section.

* **How can I contribute?**
  Check out the `Contribution`_ section.
