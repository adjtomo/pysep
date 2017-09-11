.. pyseis documentation master file, created by
   sphinx-quickstart on Sat Sep  9 20:05:32 2017.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.


``pyseis`` documentation
==================================
.. contents::

Most users will want to use the latest version of obspy and the latest version of python. Here is how to install the latest version of obspy with python 3.6.


1. Install Miniconda and create conda environment
2. Install obpsy inside the environment
3. Clone ``pyseis`` and run the example


* Before starting update the software on your computer ::
    
    $ sudo yum update


Check out here if you already have conda installed and want to update the enviornment

Install miniconda
-----------------

* At the bottom of your ~/.bashrc, add these lines ::

    # to over-ride the PYTHONPATH settings from antelope startup
    unset PYTHONPATH
    set --

* For a clean installation, we suggest creating a python environment and install python, obspy (and dependcies) inside it. To create environments you will need to install miniconda first ::

    $ wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh
    $ bash Miniconda2-latest-Linux-x86_64.sh
    $ source ~/.bashrc

Answer the prompts. This will append your .bashrc file with the path to miniconda.
This will install several packages, including a version of python.

* Install the conda client and check the conda version ::

    $ conda install anaconda-client
    $ conda --version

This will install the directory ~/miniconda2. This will take several minutes.


Create conda environment
++++++++++++++++++++++++

* Make sure you’re in your home dir and that you are not inside the python shell.
  Create the sln environment with python ::

   $ conda create --name sln python=3.6
   $ source activate sln
   (sln)$ python --version


Install obspy
-------------

* Install obspy and its dependencies ::

   (sln)$ conda config --add channels conda-forge
   (sln)$ conda install obspy jupyter basemap

Answer the prompts. This will take several minutes.

.. note::
   
     If you want to always have python 3.6 and the latest obspy in your path, then add this line at the end of your .bashrc file ::
   
       # for obspy and python
       source activate sln

Run ``pyseis`` example
----------------------

* Get the pyseis repository ::

    git clone https://GITHUBUSERNAME@github.com/uafseismo/pyseis.git pyseis

* Run ``pyseis`` test ::
   
     (sln)$ cd pyseis
     (sln)$ python run_getwaveform.py


Test history
------------

+------------+-----------+------------+--------+--------+-------+--------+--------+
| date       | user      | machine    | conda  | python | obspy | pyseis | wfdiff |
+============+===========+============+========+========+=======+========+========+
| 2017-09-07 | crichards | otter      | 4.3.25 | 3.6.2  | 1.0.3 | works  | works  |
+------------+-----------+------------+--------+--------+-------+--------+--------+
| 2017-08-23 | carltape  | grizzly    | 4.3.25 | 3.6.2  | 1.0.3 | works  | works  |
+------------+-----------+------------+--------+--------+-------+--------+--------+ 
| 2017-08-31 | ksmith    | ptarmigan  | 4.3.25 | 3.6.2  | 1.0.3 | works  | works  | 
+------------+-----------+------------+--------+--------+-------+--------+--------+


..
   :maxdepth: 2
   :caption: Contents:

..
   Indices and tables
   ==================

   * :ref:`genindex`
   * :ref:`modindex`
   * :ref:`search`
