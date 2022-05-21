import logging

logger = logging.getLogger("pysep")
logger.setLevel("INFO")
logger.propagate = 0
ch = logging.StreamHandler()
FORMAT = "[%(asctime)s] - %(name)s - %(levelname)s: %(message)s"
formatter = logging.Formatter(FORMAT, datefmt="%Y-%m-%d %H:%M:%S")
ch.setFormatter(formatter)
logger.addHandler(ch)

from pysep.pysep import Pysep
from pysep.recsec import RecordSection
