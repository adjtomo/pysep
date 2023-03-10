"""
Quick script to update example Config files if API changes occur
"""
import os
from glob import glob
from pysep import Pysep

for dir_ in glob("*"):
    if not os.path.isdir(dir_):
        continue
    for fid in glob(os.path.join(dir_, "*.yaml")):
        sep = Pysep(config_file=fid)
        sep.load()
        sep.write_config(f"{fid}")


