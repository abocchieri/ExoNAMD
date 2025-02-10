Installation Guide
==================

Prerequisites
-------------

* Python 3.8+
* numpy
* pandas
* matplotlib
* swifter
* requests
* loguru
* astroquery

Installing ExoNAMD
------------------

From PyPI
^^^^^^^^^

.. code-block:: bash

    pip install exonamd

From Source
^^^^^^^^^^^

.. code-block:: bash

    git clone https://github.com/abocchieri/ExoNAMD.git
    cd ExoNAMD
    pip install -e .

Verifying Installation
----------------------

You can verify the installation by running a simple example:

.. code-block:: python

    from exonamd.run import run
    run(from_scratch=True)
