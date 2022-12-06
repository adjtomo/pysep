import os
from setuptools import setup, find_packages


install_requires = [
        "obspy>=1.2.2", 
        "cartopy>=0.20.2", 
        "pyyaml>=5.4",
        "llnl_db_client @ git+https://github.com/krischer/llnl_db_client.git",
        ]

setup(name="pysep",
      version='0.2.0',
      description="Python Seismogram Extraction and Processing",
      url="http://github.com/adjtomo/pysep",
      author='adjTomo Dev Team',
      license='GPL-3.0',
      python_requires=">=3.7",
      packages=find_packages(),
      install_requires=install_requires,
      entry_points={"console_scripts": ["pysep=pysep.pysep:main",
                                        "recsec=pysep.recsec:main"]
                                        },
      zip_safe=False
      )
