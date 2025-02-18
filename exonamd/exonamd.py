import sys

from exonamd import logger
from exonamd.log import addLogFile, setLogLevel

from exonamd import __pkg_name__
from exonamd import __version__
from exonamd.run import run


def main():
    """
    Main entry point for the ExoNAMD command line interface.

    This function parses arguments using argparse and runs the `run` function
    with the appropriate arguments. The arguments are:

    * `-u` or `--update`: Update the exo.csv database.
    * `-d` or `--debug`: Enable debug mode.
    * `-l` or `--log`: Enable logging to file.

    If debug mode is not enabled, the logging level is set to INFO and
    logging to file is disabled. If logging to file is enabled, a log file
    is created with the current date and time in the format
    `{__pkg_name__}_{time.strftime('%Y%m%d_%H:%M:%S')}.log`.

    The function logs a message at the `Announce` level when it starts and
    finishes, and logs the elapsed time at the `INFO` level.
    """
    import argparse
    import time

    start = time.time()

    parser = argparse.ArgumentParser(description="Run ExoNAMD from the command line.")

    parser.add_argument(
        "-u",
        "--update",
        dest="update",
        default=False,
        required=False,
        action="store_true",
        help="Update the exo.csv database",
    )

    parser.add_argument(
        "-d",
        "--debug",
        dest="debug",
        default=False,
        required=False,
        action="store_true",
        help="Enable debug mode",
    )

    parser.add_argument(
        "-l",
        "--log",
        dest="log",
        default=False,
        required=False,
        action="store_true",
        help="Enable logging to file",
    )

    args = parser.parse_args()

    if not args.debug:
        logger.remove(0)
        setLogLevel("INFO")

    if args.log:
        addLogFile(f"{__pkg_name__}_{time.strftime('%Y%m%d_%H:%M:%S')}.log")
        logger.info("Logging to file enabled")

    logger.log("Announce", f"Starting {__pkg_name__} v{__version__}")
    run(from_scratch=args.update)

    end = time.time()
    logger.info(f"Finished in {end - start:.2f} seconds")
    logger.log("Announce", f"Exiting {__pkg_name__}")


if __name__ == "__main__":
    main()
