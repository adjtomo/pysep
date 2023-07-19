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

