# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- General information -----------------------------------------------------
import os
import sys

sys.path.insert(0, os.path.abspath("../../src"))


# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = "adme-py"
copyright = "2024, Andrew McNaughton"
author = "Andrew McNaughton"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
]

templates_path = ["_templates"]
exclude_patterns = []


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "alabaster"

html_theme_options = {
    "logo": "logo.png",
    "logo_name": True,  # Display the project name next to the logo
    "description": "Tool to calculate ADME properties in python.",  # Optional: Add a short description
    "github_user": "mcnaughtonadm",
    "github_repo": "adme-py",
}
html_static_path = ["_static"]
