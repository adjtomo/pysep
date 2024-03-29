import logging
from importlib.metadata import version

# Defines version number from 'pyproject.toml'
__version__ = version("pysep-adjtomo")

logger = logging.getLogger("pysep")
logger.setLevel("INFO")
logger.propagate = 0
ch = logging.StreamHandler()
FORMAT = "[%(asctime)s] - %(name)s - %(levelname)s: %(message)s"
formatter = logging.Formatter(FORMAT, datefmt="%Y-%m-%d %H:%M:%S")
ch.setFormatter(formatter)
logger.addHandler(ch)

# from pysep.pysep import Pysep, get_data  # NOQA
from pysep.pysep import Pysep, get_data  # NOQA
from pysep.recsec import RecordSection, plotw_rs  # NOQA
from pysep.declust import Declust  # NOQA
from pysep.utils.io import read_sem, read_stations, read_events_plus  # NOQA
