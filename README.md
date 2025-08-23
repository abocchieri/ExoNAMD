# ``ExoNAMD``

[![PyPI version](https://badge.fury.io/py/exonamd.svg?icon=si%3Apython)](https://badge.fury.io/py/exonamd)
[![GitHub version](https://badge.fury.io/gh/abocchieri%2Fexonamd.svg)](https://badge.fury.io/gh/abocchieri%2Fexonamd)
[![Downloads](https://static.pepy.tech/badge/exonamd)](https://pepy.tech/project/exonamd)
[![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Documentation Status](https://readthedocs.org/projects/exonamd/badge/?version=stable)](https://exonamd.readthedocs.io/en/stable/?badge=stable)

## Introduction

``ExoNAMD`` is a Python package to compute the Normalized Angular Momentum Deficit (NAMD) of exoplanetary systems. The NAMD is a measure of the dynamical excitation of a planetary system, and it can be defined as the difference between the total angular momentum of the system and the angular momentum it would have if all planets were on circular, coplanar orbits.

``ExoNAMD`` is designed to be fast, modern, and reliable. It is built using the latest Python features and libraries, and it is tested extensively to ensure its accuracy and reliability.

## Table of contents

- [``ExoNAMD``](#exonamd)
  - [Introduction](#introduction)
  - [Table of contents](#table-of-contents)
  - [How to install](#how-to-install)
    - [Install from PyPI](#install-from-pypi)
    - [Install from source code](#install-from-source-code)
      - [Test your installation](#test-your-installation)
  - [Documentation](#documentation)
    - [Build the html documentation](#build-the-html-documentation)
    - [Build the pdf documentation](#build-the-pdf-documentation)
  - [How to contribute](#how-to-contribute)
  - [How to cite](#how-to-cite)

## How to install

Instructions on how to install ``ExoNAMD``.

### Install from PyPI

``ExoNAMD`` is available on PyPI and can be installed via pip as

    pip install exonamd

### Install from source code

``ExoNAMD`` is compatible (tested) with Python 3.8+

To install from source, clone the [repository](https://github.com/abocchieri/ExoNAMD) and move inside the directory.

Then use `pip` as

    pip install .

#### Test your installation

Try importing ``ExoNAMD`` as

    python -c "import exonamd; print(exonamd.__version__)"

Or running ``ExoNAMD`` itself with the `help` flag as

    exonamd -h

If there are no errors then the installation was successful!

## Documentation

``ExoNAMD`` comes with an extensive documentation, which can be built using Sphinx.
The documentation includes a tutorial, a user guide and a reference guide.

To build the documentation, install the needed packages first via `poetry`:

    pip install poetry
    poetry install --with docs

### Build the html documentation

To build the html documentation, move into the `docs` directory and run

    make html

The documentation will be produced into the `build/html` directory inside `docs`.
Open `index.html` to read the documentation.

### Build the pdf documentation

To build the pdf, move into the `docs` directory and run

    make latexpdf

The documentation will be produced into the `build/latex` directory inside `docs`.
Open `exonamd.pdf` to read the documentation.

The developers use `pdflatex`; if you have another compiler for LaTex, please refer to [sphinx documentation](https://www.sphinx-doc.org/en/master/usage/configuration.html#latex-options).

## How to contribute

You can contribute to ``ExoNAMD`` by reporting bugs, suggesting new features, or contributing to the code itself.
If you wish to contribute to the code, please follow the steps described in the documentation under `Developer Guide`.

## How to cite

```bibtex
@ARTICLE{Bocchieri2025,
       author = {{Bocchieri}, Andrea and {Zak}, Jiri and {Turrini}, Diego},
        title = "{ExoNAMD: Leveraging the spin-orbit angle to constrain the dynamics of multiplanetary systems}",
      journal = {in preparation},
         year = 2025,
}
```

