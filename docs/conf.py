# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import os
import sys

# Make both source packages importable during the build.
sys.path.insert(0, os.path.abspath("../src"))

# -- Project information -----------------------------------------------------

project = "MutaDock"
copyright = "2026, Naisarg Patel"
author = "Naisarg Patel"

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",  # Pull docstrings from source
    "sphinx.ext.napoleon",  # Support Google-style (Args/Returns/Raises) docstrings
    "sphinx.ext.viewcode",  # Add [source] links to API pages
    "sphinx.ext.intersphinx",  # Cross-reference to Python stdlib docs
]

intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
}

# autodoc settings
autodoc_default_options = {
    "members": True,
    "undoc-members": False,
    "show-inheritance": True,
    "member-order": "bysource",
}
autodoc_typehints = "description"  # Render type hints in the description section
autodoc_typehints_description_target = "documented"

# napoleon settings (match the Google-style used in the codebase)
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = True
napoleon_attr_annotations = True

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

# -- Options for HTML output -------------------------------------------------

html_theme = "alabaster"
html_static_path = ["_static"]

html_theme_options = {
    "description": "Mutation analysis and multi-receptor docking toolkit",
    "github_user": "naisarg14",
    "github_repo": "mutadock",
    "fixed_sidebar": True,
}
