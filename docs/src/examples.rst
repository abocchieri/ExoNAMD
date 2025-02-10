Examples and Tutorials
====================

Basic Examples
------------

Calculate NAMD for All Systems
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from exonamd.run import run
    
    # Full workflow from scratch
    run(from_scratch=True)

Plot NAMD for a Specific System
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from exonamd.plot import plot_host_namd
    import pandas as pd
    
    # Load data
    df = pd.read_csv("exo_namd.csv")
    
    # Plot relative NAMD for WASP-47
    plot_host_namd(df, "WASP-47", kind="rel")
    
    # Plot absolute NAMD
    plot_host_namd(df, "WASP-47", kind="abs")

Advanced Analysis
--------------

Working with Core Sample
^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: python

    from exonamd.run import calc_namd
    
    # Select only well-constrained systems
    df = calc_namd(df, core=True)
    
    # Visualize results
    plot_sample_namd(df, title="Core Sample")
