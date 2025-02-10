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

Installing ``ExoNAMD``
----------------------

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

To test for correct setup you can do

.. code-block:: console

    python -c "import exonamd"

If no errors appeared then it was successfully installed.

Additionally the ``ExoNAMD`` program should now be available in the command line

.. code-block:: console

    exonamd


Uninstall ``ExoNAMD``
---------------------

``ExoNAMD`` is installed in your system as a standard python package:
you can uninstall it from your Environment as

.. code-block:: console

    pip uninstall exonamd


Update ``ExoNAMD``
------------------

If you have installed ``ExoNAMD`` using Pip, now you can update the package simply as

.. code-block:: console

    pip install exonamd --upgrade

If you have installed ``ExoNAMD`` from GitHub, you can download or pull a newer version of ``ExoNAMD`` over the old one.

Then you have to place yourself inside the installation directory with the console

.. code-block:: console

    cd /your_path/ExoNAMD

Now you can update ``ExoNAMD`` simply as

.. code-block:: console

    pip install . --upgrade

or simply

.. code-block:: console

    pip install .

Modify ``ExoNAMD``
------------------

You can modify ``ExoNAMD`` main code, editing as you prefer, but in order to make the changes effective

.. code-block:: console

    pip install . --upgrade

or simply

.. code-block:: console

    pip install .

To produce new ``ExoNAMD`` functionalities and contribute to the code, please see :ref:`Developer Guide`.