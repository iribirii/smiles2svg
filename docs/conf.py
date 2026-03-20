# Configuration file for the Sphinx documentation builder.
# https://www.sphinx-doc.org/en/master/usage/configuration.html

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "src"))

project = "smiles2svg"
copyright = "2023, Inigo Iribarren Aguirre"
author = "Inigo Iribarren Aguirre"

release = "0.2.0"
version = "0.2.0"

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
]

templates_path = ["_templates"]
exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]

html_theme = "furo"
html_title = f"{project} {version}"
html_static_path = ["_static"]

html_theme_options = {
    "source_repository": "https://github.com/iribirii/smiles2svg",
    "source_branch": "master",
    "source_directory": "docs/",
    "light_css_variables": {
        "color-brand-primary": "#156090",
        "color-brand-content": "#156090",
    },
    "dark_css_variables": {
        "color-brand-primary": "#5eb3e0",
        "color-brand-content": "#5eb3e0",
    },
}

autodoc_member_order = "bysource"
napoleon_google_docstring = True
