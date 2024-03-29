# Configuration file for the Sphinx documentation builder.
#
# This file only contains a selection of the most common options. For a full
# list see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Path setup --------------------------------------------------------------

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#
#import os
#import sys
#sys.path.insert(0, os.path.abspath('.'))
import os
import sys
from pathlib import Path

localpath = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(localpath)
import ucdeconvolve

# -- Project information -----------------------------------------------------

project = 'UniCell Deconvolve'
copyright = '2023, Daniel Charytonowicz'
author = 'Daniel Charytonowicz'
repository_url = "https://github.com/dchary/ucdeconvolve"
github_url = "https://github.com/dchary/ucdeconvolve"
# The full version, including alpha/beta/rc tags
release = '0.1.0'


# -- General configuration ---------------------------------------------------


# Add any Sphinx extension module names here, as strings. They can be
# extensions coming with Sphinx (named 'sphinx.ext.*') or your custom
# ones.

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.intersphinx",
    "sphinx.ext.doctest",
    "sphinx.ext.coverage",
    "sphinx.ext.mathjax",
    "sphinx.ext.napoleon",
    "sphinx.ext.autosummary",
    "nbsphinx",
    "myst_parser"
]

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = "bysource"
issues_github_path = "dchary/ucdeconvolve"
#napoleon_google_docstring = False
napoleon_numpy_docstring = True
#napoleon_include_init_with_doc = False
#napoleon_use_rtype = True  # having a separate entry generally helps readability
#napoleon_use_param = True
napoleon_custom_sections = [("Params", "Parameters")]
#typehints_defaults = "braces"
#todo_include_todos = False
#suppress_warnings = ["ref.citation"]


# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
# This pattern also affects html_static_path and html_extra_path.
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']


# -- Options for HTML output -------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.
#
html_theme = 'sphinx_rtd_theme'

html_theme_options = {
    'navigation_depth': 9,
    "logo_only": True,
    "collapse_navigation" : False
}

add_function_parentheses = False
nbsphinx_allow_errors = True
nbsphinx_execute = 'never'
html_css_files = [
    'css/override.css',
]

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']
html_show_sphinx = False
html_logo = '_static/img/logo_unicell.png'
html_title = "ucdeconvolve"

html_sidebars = {
    "**": ["logo-text.html", "globaltoc.html", "localtoc.html", "searchbox.html"]
}
