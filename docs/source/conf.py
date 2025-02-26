# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'Broadband_albedo_calculator'
copyright = '2025, Vidya Venkatesan'
author = 'Vidya Venkatesan'
release = '2025'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

# Required Extensions
extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.napoleon',           # For Google-style docstrings
    'sphinx.ext.viewcode',           # Link to source code
    'sphinx.ext.autodoc.typehints',  # Type hints in docstrings
    'sphinx.ext.todo',               # Support for TODO directives
    'sphinx.ext.githubpages'         # For GitHub Pages compatibility
]
import os
import sys
sys.path.insert(0, os.path.abspath('/Users/astrovidee/Dropbox/Broadband_albedo_Calculator'))


templates_path = ['_templates']
exclude_patterns = []

language = 'Python'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_static_path = ['_static']
html_theme = 'sphinx_rtd_theme'
html_theme_options = {
    'collapse_navigation': False,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False
}
autodoc_typehints = 'description'
autodoc_member_order = 'bysource'

html_show_sourcelink = True
suppress_warnings = ['autodoc.mock']
