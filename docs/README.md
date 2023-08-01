# PySEP Documentation

The PySEP documentation is built with Sphinx and ReadTheDocs. The official 
documentation can be found at: https://pysep.readthedocs.io

Docs building is automatically triggered when updates are pushed to PySEP. 
Each branch of the code may have a different documentation. 
The 'latest' version of the docs points to the 'devel' branch of the code.

## Building Docs Locally

In order to build the docs locally, you will first need to create a separate 
Conda environment with a few packages, you can do this by running:

``` bash
conda env create --file environment.yaml
conda activate pysep-docs
```

You can then run the make command to generate the .html files. You can find your 
local docs in the *_build/html* directory

```bash
make html
```

See your locally built documentation by opening *_build/html/index.html*

## Publishing Package on PyPi 
*Useful link: https://realpython.com/pypi-publish-python-package/*

1. Ensure your `pyproject.toml` file is set up properly; required fields are name and version
2. Set dependencies, do **not(( pin exact versions but allow for upper and lower bounds; only list direct dependencies
3. Include `tests/`, `docs/`, license, and MANIFEST files (MANIFIST used for including non-source code material
4. Ensure you have an account on PyPi and TestPyPi (for testing publishing)
5. Install `twine` and `build` which are used to build and push packages to PyPi
6. Build your packages locally, which creates the `.tar.gz` and `.whl` dist files
```bash
python -m build
```
6. Check that files in your .whl (zip file) are as expected (including everything in 3)
7. Check dist files with:
```bash
twine check dist/*
```
8. Upload test package (note: requires TestPyPi account)
```bash
twine upload -r testpypi dist/*
```
9. Upload real package (note: requires PyPi account)
```bash
twine upload dist/*
```
