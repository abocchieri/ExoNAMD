"""
This module was directly inspired by the log module 
in the `PAOS <https://github.com/arielmission-space/PAOS>`_ package.
"""

import sys
from exonamd import logger


def addLogFile(fname="exonamd.log"):
    """
    Adds a new log file for logging.

    Parameters
    ----------
    fname : str, optional
        The filename for the log file. Defaults to 'exonamd.log'.

    Returns
    -------
    None
    """
    logger.add(fname)


def setLogLevel(level="INFO"):
    """
    Configures the logging level for the logger.

    Parameters
    ----------
    level : str, optional
        The logging level to set. Defaults to 'INFO'. Possible values include
        'DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL', etc.

    Returns
    -------
    None
    """
    logger.configure(handlers=[{"sink": sys.stderr, "level": level}])


def disableLogging(name="exonamd"):
    """
    Disables the logger for the given module name.

    Parameters
    ----------
    name : str, optional
        The module name to disable logging for. Defaults to 'exonamd'.

    Returns
    -------
    None
    """
    logger.disable(name)


def enableLogging(name="exonamd"):
    """
    Enables the logger for the given module name.

    Parameters
    ----------
    name : str, optional
        The module name to enable logging for. Defaults to 'exonamd'.

    Returns
    -------
    None
    """
    logger.enable(name)
