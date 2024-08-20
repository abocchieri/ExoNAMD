import importlib.metadata as metadata
from datetime import date

from .__version__ import __version__

# load package info
__pkg_name__ = metadata.metadata("exonamd")["Name"]
__author__ = metadata.metadata("exonamd")["Author"]
__license__ = metadata.metadata("exonamd")["license"]
__copyright__ = "2024-{:d}, {}".format(date.today().year, __author__)
__summary__ = metadata.metadata("exonamd")["Summary"]

from loguru import logger

# instantiate logger
logger.level("Announce", no=100, color="<magenta>")
