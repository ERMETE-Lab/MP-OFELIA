# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import importlib
from pathlib import Path

import os
import sys

sys.path.insert(0, os.path.abspath('..'))
sys.path.insert(0, os.path.abspath('../models'))
sys.path.insert(0, os.path.abspath('../Tutorials'))

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information
author = "Lorenzo Loi & Stefano Riva, Carolina Introini, Antonio Cammi"
project = 'OFELIA'
copyright = f"2024, {author}"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_parser",
    'sphinx.ext.mathjax',
    "nbsphinx",
    # "sphinx.ext.autodoc",
    "sphinx.ext.todo",
    "sphinx.ext.viewcode",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.intersphinx",
    "sphinxcontrib.bibtex",
    "sphinx_proof",
    "IPython.sphinxext.ipython_console_highlighting"
]

myst_enable_extensions = ["dollarmath", "amsmath"]

mathjax3_config = {
    'tex': {'tags': 'ams', 'useLabelIds': True},
}

source_suffix = {
    '.rst': 'restructuredtext',
    '.txt': 'markdown',
    '.md': 'markdown',
}

autosummary_generate = True  # Turn on sphinx.ext.autosummary

bibtex_bibfiles = ['references.bib']

nb_execution_mode = "off"

apidoc_module_dir = f"../models/"
apidoc_excluded_paths = ["tests", "build", "pyforce.egg-info"]
apidoc_toc_file = False

autodoc_default_options = {"members": True}

# language = "Python"
suppress_warnings = [
    'nbsphinx',
]
nbsphinx_execute = 'never'

here = Path(__file__).parent.resolve()

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']
autodoc_mock_imports = ["dolfinx", "ufl", "petsc4py", "pyvista", "h5py"]

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

# Options for HTML output
html_theme = 'sphinx_rtd_theme'
# html_theme = 'pydata_sphinx_theme'
# html_theme = 'sphinx_book_theme'