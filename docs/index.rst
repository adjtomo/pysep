.. pysep documentation master file, created by
   sphinx-quickstart on Wed Feb 15 13:04:35 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Python Seismogram Extraction and Processing
===========================================

`PySEP` is a package revolving around the collection, visualization and curation
of seismic data and metadata. It is built around ObsPy and uses a number of 
core classes to achieve its goals:

- ``PySep``: data download and processing tool
- ``RecSec``: record section plotter for data visualization
- ``Declust``: event and station declustering and weighting

As part of the `SCOPED toolkit <https://github.com/SeisSCOPED>`_, `PySEP` 
has been
`containerized <https://github.com/SeisSCOPED/pysep/pkgs/container/pysep>`_
using Docker.

`PySEP` is hosted on `GitHub <https://github.com/adjtomo/pysep>`_ as part of the  
`adjTomo organization <https://github.com/adjtomo>`_.

---------------------------------

Installation
------------

We recommend installing PySEP into a Conda environment. Most dependencies are 
installed via Conda while PySEP itself is installed in editable mode via Pip.

This installation also creates two command line entry points for 
`pysep` and `recsec`. See the `pysep <pysep.html>`_ and `recsec <recsec.html>`_
documentation pages for how to use these command line tools.

Installing to a new Conda environment (recommended)
```````````````````````````````````````````````````
To create a new Conda environment for PySEP and its dependencies:

.. code:: bash

    git clone --branch devel https://github.com/uafgeotools/pysep.git
    cd pysep
    conda env create -f environment.yml
    conda activate pysep

Updating an existing Conda environment
``````````````````````````````````````

If you have an existing Conda environment that you want to install PySEP and its
dependencies into, you can do so with the following commands:

.. code:: bash

    git clone -- branch devel https://github.com/uafgeotools/pysep.git
    cd pysep
    conda activate <your environment>
    conda env update -f environment.yml

.. note:: 
   
   We are currently working to get `PySEP hosted on PyPi and Conda 
   <https://github.com/adjtomo/pysep/issues/55>`_


---------------------------------

Running Tests
-------------

PySEP comes with unit testing which should be run before and after making any
changes to see if your changes have altered the expected code behavior.

.. code:: bash

    cd pysep/tests
    pytest

---------------------------------

.. toctree::
   :maxdepth: 2
   :hidden:

   pysep
   recsec
   declust
   cookbook
   usage


