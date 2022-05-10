import os
from setuptools import setup, find_packages


setup(name="pysep",
      version='0.1.0',
      description="Python Seismogram Extraction and Processing",
      url="http://github.com/uafgeotools/pysep",
      author='UAFGeotools',
      license='GPL-3.0',
      python_requires=">=3.7",
      packages=find_packages(),
      install_requires=["obspy", "cartopy", "pyyaml"],
      zip_safe=False
      )
