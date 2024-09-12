.. MutaDock documentation master file, created by
   sphinx-quickstart on Thu Sep 12 15:21:49 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MutaDock documentation
======================
.. toctree::
   :maxdepth: 2
   :caption: Contents:

Overview
--------

MUTADOCK is a comprehensive library designed for mutation studies and
multiple receptor-ligand docking. It provides tools and methods to
analyze and predict the effects of mutations on receptor-ligand
interactions, enabling researchers to study protein function and drug
binding affinity in a detailed manner.

Features
--------

-  Mutation Studies: Analyze the impact of various mutations on protein
   structure and function.
-  Multiple Receptor-Ligand Docking: Perform docking simulations
   involving multiple receptor and ligand combinations.
-  User-Friendly Interface: Easy-to-use commands and functions for quick
   analysis and simulations.

Installation
------------

To install MUTADOCK, follow these steps:

1. Install using pip

::

   pip install mutadock

OR
~~

1. Clone the repository:

::

   git clone https://github.com/naisarg14/mutadock.git

Navigate to the project directory:

::

   cd MUTADOCK

3. Install dependencies including Pyrosetta:

-  For Linux:

::

   ./install.sh

-  For Windows:

::

   win_install.bat

3. Install dependencies without PyRosetta (PyRosetta is required and
   should be installed seperately):

::

   pip install -r requirements.txt

Usage
-----

Here’s a brief guide to using MUTADOCK:

1. Common Commands:

::

   md_mutate --input <input_file> --output <output_file>

::

   md_dock --input <receptor_txt> --output <ligand_txt> --config <config_file>

2. For more detailed usage options, see the
   `documentation <https://mutadock.readthedocs.io/en/latest/>`__.

License
-------

This project is licensed under the GNU General Public License v3.0
(GPL-3.0). See the `LICENSE <LICENSE>`__ file for details.

Contact
-------

For any questions or feedback, please contact
naisarg.patel14@hotmail.com.

Acknowledgments
---------------

-  Thanks to the contributors and libraries that made this project
   possible.
