# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

# type: ignore
import os
import sys
project = 'Trojan Detection'
copyright = '2026, Souren Ishkhanian'
author = 'Souren Ishkhanian'
release = '2026-02-11'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

sys.path.insert(0, os.path.abspath('../..'))
extensions = [
    'sphinx.ext.autodoc',   # reads docstrings from your modules
    'sphinx.ext.napoleon',  # supports Google/NumPy style docstrings
    'sphinx.ext.viewcode',  # adds links to source code
]

templates_path = ['_templates']
exclude_patterns = []

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']

html_theme_options = {
    "light_logo": "RUB_Logo_Full_Black.jpg",
    "dark_logo": "RUB_Logo_Full_Inverted.jpg"
}
html_favicon = '_static/RUB_Emblem.svg'
