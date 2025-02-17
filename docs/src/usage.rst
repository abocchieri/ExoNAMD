Usage Guide
===========

Basic Usage
-----------

The main workflow consists of three steps:

1. Create/update the database of multiplanetary systems:

.. code-block:: python

    from exonamd.run import create_db
    
    # Download full database
    df = create_db(from_scratch=True)
    
    # Or update existing database
    df = create_db(from_scratch=False)

The data is downloaded from the `NASA Exoplanet Archive <https://exoplanetarchive.ipac.caltech.edu/>`_.

The database will be populated with the following columns:

- `hostname`: Host star name
- `pl_name`: Planet name
- `default_flag`: Default flag, marking whether it is the default value in the NASA database
- `rowupdate`: Last update
- `sy_pnum`: Number of planets in the system
- `st_rad`: Stellar radius
- `st_mass`: Stellar mass
- `pl_bmasse`: Planet mass
- `pl_rade`: Planet radius
- `pl_orbper`: Orbital period
- `pl_orbsmax`: Semi-major axis
- `pl_orbeccen`: Eccentricity
- `pl_orbincl`: Inclination
- `pl_trueobliq`: True obliquity
- `pl_ratdor`: Ratio of the semi-major axis to the stellar radius
- `pl_ratror`: Ratio of the planet radius to the stellar radius

The `create_db` function also deals with aliases, computes missing values from simple relations, and stores the curated database.


2. Interpolate missing values:

.. code-block:: python

    from exonamd.run import interp_db
    
    df = interp_db(df)

This function reloads the database, thins it down by using the median of the values for each parameter of each planet, and then interpolates missing values of:

- eccentricities;
- planetary masses;
- inclinations;
- semi-major axis uncertainties;
- true obliquity.

The function then stores the curated+interpolated database in a new file.

3. Calculate NAMD:

.. code-block:: python

    from exonamd.run import calc_namd
    
    df = calc_namd(df, core=True)

This function calculates the NAMD for each system in the database. The `core` parameter is used to restrict the calculation to the "core" sample, defined as the systems with all planets having a flag of either 0, 05+, 05-, or 05+-, i.e. nothing or only the obliquity has been interpolated.

Finally, the function stores the database with the NAMD values in a new file.

.. note::

    The `calc_namd` function computes both the relative and absolute NAMD values. The relative NAMD is calculated using the relative inclination with respect to the most massive planet, while the absolute NAMD uses the true obliquity values.

.. tip::

    The relative NAMD is defined in `Turrini et al. (2021) <https://doi.org/10.1051/0004-6361/201936301>`_ and the equation is given also in our paper `Bocchieri et al. (2025) <https://doi.org/TODO>`_, which contains the definition of the absolute NAMD.

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
