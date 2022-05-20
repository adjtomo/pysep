import os
from setuptools import setup, find_packages


install_requires = [
        "obspy>=1.2.2", 
        "cartopy>=0.20.2", 
        "pyyaml>=5.4",
        "llnl_db_client @ git+https://github.com/krischer/llnl_db_client.git",
        ]

setup(name="pysep",
      version='0.1.0',
      description="Python Seismogram Extraction and Processing",
      url="http://github.com/uafgeotools/pysep",
      author='UAFGeotools',
      license='GPL-3.0',
      python_requires=">=3.7",
      packages=find_packages(),
      install_requires=[
      zip_safe=False
      )
