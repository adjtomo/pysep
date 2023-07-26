.. pysep documentation master file, created by
   sphinx-quickstart on Wed Feb 15 13:04:35 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

PySEP
===========================================

``PySEP`` (`Python Seismogram Extraction and Processing`) is a package revolving 
around the collection, visualization and curation
of seismic data and metadata. It is built around ObsPy and uses a number of 
core classes to achieve its goals:

- `PySep <pysep.html>`__: data download and processing tool
- `RecSec <recsec.html>`__: record section plotter for data visualization
- `Declust <declust.html>`__: event and station declustering and weighting

`PySEP` is hosted on `GitHub <https://github.com/adjtomo/pysep>`_ as part of the  
`adjTomo organization <https://github.com/adjtomo>`_.

As part of the `SCOPED toolkit <https://github.com/SeisSCOPED>`_, `PySEP` 
has been
`containerized <https://github.com/SeisSCOPED/pysep/pkgs/container/pysep>`_
using Docker.

---------------------------------

Installation
------------

We recommend installing PySEP into a `Conda <https://conda.io>`__ environment 
to avoid package or dependency conflicts with other Python packages.

Installation via Pip is recommended for the latest, stable version of 
PySEP. Note that the package name `PySEP` was already take on PyPi, so our 
package must be installed via the name: ``pysep-adjtomo``.

.. note::

    Cartopy must be installed separately via Conda otherwise you 
    will encounter Pip dependency errors


.. code:: bash

    conda create -n pysep 
    conda activate pysep
    conda install cartopy
    pip install pysep-adjtomo

------------------------------------

Installing Development Version
``````````````````````````````

PySEP is under an active state of development, so for the latest version of the
codebase, installation must take place directly from the ``devel`` branch of the 
code.

.. warning::

    API and code stability is subject to change without warning when using the
    ``devel`` branch

.. code:: bash

    git clone --branch devel https://github.com/adjtomo/pysep.git
    cd pysep/
    conda env create -f environment.yml
    conda activate pysep

The above code installs PySEP in *editable* mode via Pip. This means that any 
source code changes that occur in the *pysep/* directory are directly 
accesible in your development version. 

To update the development branch with the most up-to-date changes, run:

.. code:: bash

    git pull origin devel


Running Tests
`````````````

PySEP comes with unit testing which should be run before and after making any
changes to see if your changes have altered the expected code behavior. 

If you have installed the Development version of PySEP, you can run tests with:

.. code:: bash

    conda install pytest
    cd pysep/tests
    pytest

---------------------------------

.. toctree::
   :maxdepth: 2
   :hidden:

   pysep
   recsec
   declust
   tutorials
   lab_record_section
   cookbook
   usage


