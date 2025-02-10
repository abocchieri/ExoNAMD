Configuration Reference
=====================

Config File Format
----------------

ExoNAMD uses YAML for configuration. Below is a complete reference of all available options.

System Configuration
------------------

.. code-block:: yaml

    system:
      temperature: 300  # Kelvin
      pressure: 1.0    # atm
      timestep: 2.0    # fs
      periodic_boundary: true

Atmospheric Parameters
-------------------

.. code-block:: yaml

    atmosphere:
      composition:
        H2: 0.85
        He: 0.15
      altitude_range: [0, 1000]  # km
      pressure_levels: [1e-6, 1.0]  # bar

Force Field Options
----------------

.. code-block:: yaml

    forcefield:
      type: "custom"
      parameters: "params.inp"
      cutoff: 12.0
