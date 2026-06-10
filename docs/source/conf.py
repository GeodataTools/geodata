# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

import re
from pathlib import Path

# Read version without importing geodata (RTD installs with --no-deps).
_version_file = Path(__file__).resolve().parents[2] / "src" / "geodata" / "_version.py"
_release = re.search(
    r'^__version__\s*=\s*["\']([^"\']+)["\']', _version_file.read_text(), re.M
)
if _release is None:
    raise RuntimeError(f"Could not parse __version__ from {_version_file}")
release = _release.group(1)

project = "Geodata"
copyright = "2025, Geodata Contributors"
author = "Geodata Contributors"

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
    "myst_nb",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "sphinx.ext.intersphinx",
    "sphinx.ext.mathjax",
    "autoapi.extension",
]

templates_path = ["_templates"]
exclude_patterns = []

myst_heading_anchors = 3
nb_execution_mode = "off"
myst_enable_extensions = [
    "amsmath",
    "colon_fence",
    "deflist",
    "dollarmath",
    "html_image",
]


intersphinx_mapping = {
    "python": ("https://docs.python.org/3", None),
    "numpy": ("https://numpy.org/doc/stable/", None),
    "xarray": ("https://docs.xarray.dev/en/stable/", None),
    "shapely": ("https://shapely.readthedocs.io/en/stable/", None),
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = "sphinx_book_theme"
html_static_path = ["_static"]


autoapi_dirs = ["../../src/geodata"]
autoapi_options = [
    "members",
    "undoc-members",
    "inherited-members",
    "show-module-summary",
    "imported-members",
]
autoapi_own_page_level = "method"
