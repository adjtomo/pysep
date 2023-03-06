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
copyright = '2023, adjTomo Dev Team'
author = 'adjTomo Dev Team'
release = ''
version = '0.3.2'

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
        "myst_parser"
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
