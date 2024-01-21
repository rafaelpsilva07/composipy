# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import os
import sys
import sphinx_readable_theme
import pydata_sphinx_theme

sys.path.insert(0, "../composipy/")

from version import __version__


project = 'composipy'
copyright = '2023, Rafael Pereira da Silva'
author = 'Rafael Pereira da Silva'
release = __version__

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "nbsphinx",
    "m2r2"
]

# Don't run notebooks
nbsphinx_execute = "never"

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']



# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

#html_theme = 'alabaster'
html_theme = 'pydata_sphinx_theme'
#html_theme = 'renku'
#html_theme_path = [pydata_sphinx_theme.get_html_theme_path()]
#html_theme = 'yummy_sphinx_theme'
#html_theme = 'sphinx-readable-theme'
html_static_path = ['_static']
#html_theme_path = [sphinx_readable_theme.get_html_theme_path()]
#html_theme = 'readable'
