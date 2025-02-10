Usage Guide
===========

Basic Usage
-----------

The main workflow consists of three steps:

1. Create/update the database:

.. code-block:: python

    from exonamd.run import create_db
    
    # Download full database
    df = create_db(from_scratch=True)
    
    # Or update existing database
    df = create_db(from_scratch=False)

2. Interpolate missing values:

.. code-block:: python

    from exonamd.run import interp_db
    
    df = interp_db(df)

3. Calculate NAMD:

.. code-block:: python

    from exonamd.run import calc_namd
    
    df = calc_namd(df, core=True)

Command Line Interface
----------------------

ExoNAMD can also be run from the command line:

.. code-block:: bash

    # Update database and run calculations
    exonamd -u

    # Run with existing database
    exonamd

    # Enable debug mode
    exonamd -d

    # Enable logging to file
    exonamd -l

Core Functions
--------------

NAMD Calculation
^^^^^^^^^^^^^^^^

The package calculates both relative and absolute NAMD:

- Relative NAMD: Uses relative inclination with respect to the most massive planet
- Absolute NAMD: Uses true obliquity values

Monte Carlo Analysis
^^^^^^^^^^^^^^^^^^^^

Uncertainty estimation is performed using Monte Carlo sampling:

.. code-block:: python

    df = calc_namd(df, core=True)  # Includes Monte Carlo analysis
