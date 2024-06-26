# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

def skip(app, what, name, obj, would_skip, options):
    if name == "__init__":
        return False
    return would_skip

def setup(app):
    app.connect("autodoc-skip-member", skip)

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'PySEP'
copyright = '2024, adjTomo Dev Team'
author = 'adjTomo Dev Team'
release = ''
# Grab version number from 'pyproject.toml'
with open("../pyproject.toml", "r") as f:
    _lines = f.readlines()
for _line in _lines:
    if _line.startswith("version"):
        version = _line.split('"')[1].strip()

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [
        'sphinx.ext.doctest',
        'sphinx.ext.intersphinx',
        'sphinx.ext.todo',
        'sphinx.ext.coverage',
        'sphinx.ext.viewcode',
        'sphinx.ext.napoleon',
        'sphinx.ext.autosummary',
        "autoapi.extension",
        "myst_parser",
        "IPython.sphinxext.ipython_console_highlighting",
        ]

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

# Need to tell the autoapi that our source code is one level up
autoapi_type = "python"
autoapi_dirs = [
        "../pysep",
        "../pysep/utils",
        "../pysep/tests"
        ]
autoapi_add_toctree_entry = True
autoapi_python_class_content = 'both'
# autoclass_content = 'both'  # show init

source_suffix = {'.rst': 'restructuredtext',
                 '.md': 'markdown'}


# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'furo'
html_static_path = ['_static']
