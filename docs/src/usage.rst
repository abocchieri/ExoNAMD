Usage Guide
===========

Basic Usage
-----------

The main workflow consists of three steps:

1. Create/update the database of multiplanetary systems:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from exonamd.run import create_db

    # Download full database
    df = create_db(from_scratch=True, out_path="path/to/output")

    # Or update existing database
    df = create_db(from_scratch=False, out_path="path/to/output")

The data is downloaded from the `NASA Exoplanet Archive <https://exoplanetarchive.ipac.caltech.edu/>`_.

The query will retrieve the following values:

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

Of these, `st_rad`, `st_mass`, `pl_ratdor`, and `pl_ratror` are used to solve for e.g. `pl_rade` if this is not available but `st_rad` and `pl_ratror` are.

Then, they are discarded from the produced database as they are not needed for the NAMD calculation. If the user supplies their own custom database instead of relying on the NASA Exoplanet Archive data, these columns are not needed and can be omitted.

For `pl_bmasse`, `pl_rade`, `pl_orbsmax`, `pl_orbeccen`, `pl_orbincl`, and `pl_trueobliq`, the associated error bars are also retrieved.

.. warning:: 

    When providing custom databases, make sure that the lower error bars are provided as negative values (this follows the NASA Exoplanet Archive convention). Otherwise, the code **will** break later on.

.. note::

    Throughout this documentation, we use the terms true obliquity or stellar obliquity to indicate the spin-orbit angle. As we comment in the paper, there is a subtle distinction between stellar obliquity and spin-orbit angle for multi-planet systems.

The `create_db` function also deals with aliases (inspired from the `gen_tso code <https://pcubillos.github.io/gen_tso/>`_).

The function then stores the curated database in a new file.


2. Interpolate missing values:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from exonamd.run import interp_db
    
    df = interp_db(df, out_path="path/to/output")

This function reloads the database, thins it down by using the median of the values for each parameter of each planet, and then interpolates missing values of:

- eccentricities;
- planetary masses;
- inclinations;
- true obliquity.

For any parameter, including the semi-major axis, if the uncertainties are NaNs we put them to zero to increase the sample.

The flagging system, to keep track of interpolated values, is described later in :ref:`Flags`.

The function then stores the curated+interpolated database in a new file.

3. Calculate NAMD:
^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from exonamd.run import calc_namd
    
    df = calc_namd(df, core=True, out_path="path/to/output")

This function calculates the NAMD for each system in the database. The `core` parameter is used to restrict the calculation to the "core" sample, defined by default as the systems with all planets having a flag of either 0, 05+, 05-, or 05+-, i.e. nothing or only the obliquity has been interpolated.

.. tip::

    The "core" sample can be customised by providing a custom filter to the `filt` parameter:

    .. code-block:: python
    
        filt = lambda x: all(x["flag"].isin(["0"]))

Finally, the function stores the database with the NAMD values in a new file.

.. note::

    The `calc_namd` function computes both the relative and absolute NAMD values by default. It can be customised by the `which` parameter. 

    The relative NAMD is calculated using the relative inclination with respect to the most massive planet, while the absolute NAMD uses the true obliquity values.

.. tip::

    The relative NAMD is defined in `Turrini et al. (2021) <https://doi.org/10.1051/0004-6361/201936301>`_ and the equation is given also in our paper `Bocchieri et al. (2025) <https://doi.org/TODO>`_, which contains the definition of the absolute NAMD.

.. note::

    The code divides the absolute NAMD value by 2 to ensure it is normalized to 1.

4. Plot the results:
^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from exonamd.run import plot_sample_namd

    plot_sample_namd(df, title="NAMD vs. Multiplicity", kind="rel", out_path="path/to/output")

This function plots the NAMD values for the systems in the database. It produces a scatter plot similar to the one shown in `Turrini et al. (2021) <https://doi.org/10.1051/0004-6361/201936301>`_, their Figure 2.

.. tip:: 

    Set the `kind` parameter to "abs" to plot the absolute NAMD vs multiplicity instead.

.. _Flags:

Flags
^^^^^

Flags are used to keep track of the interpolated values. The flags are stored in the database produced by ``ExoNAMD`` and are used to interpret the results. The flags are as follows:

- [1]: Eccentricity
- [2]: Mass
- [3]: Inclination
- [4]: Semi-major axis
- [5]: Stellar obliquity
- [-]: Associated lower errorbar
- [+]: Associated upper errorbar
- [d]: Do not use

"0" is set at the beginning of the process, and the flags are updated as the values are interpolated.

For example, if we interpolated the eccentricity and stellar obliquity, together with their uncertainties, the resulting flag would be [01+-5+-]. 

.. note::

    Missing error bars are interpolated by setting them to zero by default to keep more targets in the sample. As a consequence, the resulting NAMD values from our Monte Carlo procedure provide a lower limit by definition. This bias is most prominent when the error bars of the eccentricity and stellar obliquity are missing.

.. warning::

    If "d" is present in the flag, the parameter is absent for all planets in the system and we mark it as a "do not use" system.

Monte Carlo Analysis
^^^^^^^^^^^^^^^^^^^^

Uncertainty estimation is performed using Monte Carlo sampling following the methodology in `Turrini et al. (2021) <https://doi.org/10.1051/0004-6361/201936301>`_:

- We draw from a truncated normal: for each parameter, we draw 100,000 random samples from a truncated normal distribution centered at the expected value, with an uncertainty obtained from the arithmetic mean of the upper and lower error bars, and within the physical bounds of each parameter (e.g. eccentricity between 0 and 1);
- We perform this step on all parameters needed for the relative and absolute NAMD calculations.

.. note::

    It is possible to use sampling from a normal distribution with rejection for samples outside of physical bounds, by setting the parameter `use_trunc_normal` to False in the function :func:`exonamd.run.calc_namd`. However, this is not recommended as the remaining number of samples after rejection may be too low. In any case, this approach yields the same results as the sampling from a truncated normal in the limit of large sample size.

Command Line Interface
----------------------

``ExoNAMD`` can be run from the command line:

.. code-block:: bash

    # Update existing database and run calculations
    exonamd -u

    # Create database from scratch and run calculations
    exonamd

    # Enable debug mode
    exonamd -d

    # Enable logging to file
    exonamd -l

These options are shown when running ``exonamd -h``, i.e. the help command.

.. tip::

    The optional arguments can be used simultaneoulsy.
