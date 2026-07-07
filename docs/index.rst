MutaDock documentation
======================
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   api/index

Introduction
------------

MUTADOCK is a comprehensive library designed for mutation studies and multiple receptor-ligand docking. It provides tools and methods to analyze and predict the effects of mutations on receptor-ligand interactions, enabling researchers to study protein function and drug binding affinity in a detailed manner.

Description
^^^^^^^^^^^

Our software is designed to facilitate protein mutation analysis and molecular docking. It integrates automated protein mutation using PyRosetta and a docking library capable of docking multiple proteins with multiple ligands.

Key Features
^^^^^^^^^^^^

Automated Protein Mutation:
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Utilizes PyRosetta for systematic protein mutations.
* Supports various mutation strategies (e.g., single-point mutations, double-point mutations and triple-point mutations).
* Allows customization of mutation and docking parameters.

Docking Library:
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* Capable of docking a list of proteins against a list of ligands.
* Employs AutoDock Vina to predict binding affinities and best poses.
* Provides detailed output files with docking scores and poses.

User Interface:
~~~~~~~~~~~~~~~~~~~~~~~~~~~

* **Command-Line Interface:** Simple CLI for both beginner and expert users.
* **Python Bindings:** The library can be imported in other codes for increased customizability by expert users

Try it Now
----------

The basic codes to perform mutation studies as used by us for our project can be found in a Jupyter Notebook here. The same notebook can be found on collab here.

How To Guide
------------

Installation
^^^^^^^^^^^^

MutaDock has been deployed on PyPi, making installation quick and simple

.. code-block::

   pip install mutadock

The Pyrosetta Installer will be automatically installed but Pyrosetta should be installed using

.. code-block::

   md_install_dependencies

this will install all dependencies including Pyrosetta.


* Currently there is a problem with the vina on PyPi, so vina needs to be installed separately, the installation guide can be found at https://autodock-vina.readthedocs.io/en/latest/installation.html

Quick demo (md_quick)
^^^^^^^^^^^^^^^^^^^^^^

One command takes a PDB ID, a mutation, and a ligand code and produces a ΔΔG and
a docked pose, fetching the structure from RCSB and the ligand from PubChem:

.. code-block::

   md_quick --pdb-id 4QJR --mutation A:386:ASN:HIS --ligand-code imatinib

It fetches and cleans the structure, computes the ΔΔG for the single mutation,
builds the mutant, fetches the ligand (``--ligand-code`` accepts a PubChem CID or
name; ``--ligand-file`` takes a local file), locates the pocket with AutoSite (or
pass ``-c config.txt``), docks, and writes the standard ``reports/`` bundle into
``mdquick_<id>/`` (or ``-o``). ``-i protein.pdb`` uses a local structure instead
of ``--pdb-id``. Requires the ``autosite`` binary on PATH when no config is given.

Mutation Studies
^^^^^^^^^^^^^^^^

Mutation Studies for a protein is a very fast process with just a PDB file of the protein as the input. (We assume for the tutorial that the name of the PDB file is “protein.pdb”)

.. code-block::

   md_mutate -i protein.pdb

By default all outputs are written next to the input structure. Pass
``-o/--output-dir DIR`` to collect them in ``DIR`` instead (created if absent):

.. code-block::

   md_mutate -i protein.pdb -o results/

Other optional arguments can be changed as required, to check the usage run

.. code-block::

   md_mutate -h

**ΔΔG rigor and units.** ΔΔG is ``score(mutant) − score(wild-type self-mutation
reference)`` with the identical protocol on both sides (so a null WT→WT mutation
scores ≈ 0). ``ddG_value`` is in **Rosetta Energy Units (REU), not kcal/mol**; a
``ddG_kcal`` column is also written. The REU→kcal/mol factor is
``REU_TO_KCAL_SCALE`` in ``src/mutadock/mutation/predict_ddG.py`` (default 0.34 ≈
1/2.94, Park et al. 2016), overridable with ``--reu-to-kcal``. The default
``--protocol min`` (repack + minimization) is reliable; ``--protocol cartesian``
is most accurate; ``--protocol fast`` skips minimization and is **screening
only** (absolute values unreliable — its reports carry a warning banner). Use
``--replicates N`` for mean ± SD and ``--pack-radius`` to size the neighbourhood.

.. code-block::

   md_mutate -i protein.pdb                       # default: min (reliable)
   md_mutate -i protein.pdb --protocol cartesian --replicates 3

The md_mutate will output CSV files and one text file, their description is in the table below:

.. list-table::
   :header-rows: 1

   * - **No.**
     - **File Name**
     - **Description**
   * - 1.
     - protein_modified_mutations_all.csv
     - Contains all possible mutations for the given protein
   * - 2.
     - protein_modified_mutations.csv
     - Contains mutations that are possible according to the PAM matrix for the given protein
   * - 3.
     - protein_modified_mutations_ddG.csv
     - The single mutation ddG values for the mutation in the File-2
   * - 4.
     - protein_modified_mutations_ddG_sorted.csv
     - Sorted File-3 from lowest to highest ddG values
   * - 5.
     - protein_modified_double_ddg.csv
     - The ddG values of the double mutation for all the combinations of the most negative single ddG compounds
   * - 6.
     - protein_modified_double_ddg_sorted.csv
     - Sorted File-5 from lowest to highest ddG values
   * - 7.
     - protein_modified_triple_ddg.csv
     - The ddG values of the triple mutation for all the combinations of the most negative double ddG compounds
   * - 8.
     - protein_modified_triple_ddg_sorted.csv
     - Sorted File-7 from lowest to highest ddG values
   * - 9.
     - protein_modified_mutants.txt
     - Generates a list of all the mutated PDB files created. Can be directly used as input for the md_dock command in our mutadock library

Choosing a substitution matrix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The substitution matrix is a **cheap pre-filter that decides which mutations are
worth the expensive ΔΔG step**. For every residue, MUTADOCK divides the matrix
score for each candidate amino acid by 100 to get ``prProb``; substitutions with
``prProb > 0`` (i.e. a *positive* matrix score) are kept in
``protein_modified_mutations.csv`` and go on to the ΔΔG pipeline, while the full
set is preserved in ``protein_modified_mutations_all.csv``. A more permissive
matrix therefore lets more candidates through — a broader, more thorough search,
but a slower one — and a stricter matrix keeps only the most conservative
substitutions, giving fewer candidates and a faster run.

**Where the choice takes effect.** Only ``md_csv_generator`` exposes
``--matrix``/``--matrix-file``. The one-shot ``md_mutate`` pipeline is
**PAM250-only** (the matrix is not configurable there), and ``md_quick`` docks a
single mutation you name explicitly, so no matrix is used. To run the full
pipeline with a non-default matrix, generate the candidate CSV yourself and feed
it to the ΔΔG stages:

.. code-block::

   md_csv_generator -i protein.pdb --matrix BLOSUM62      # -> protein_mutations.csv
   md_ddg_single -i protein_mutations.csv                 # then double/triple as needed

**Available matrices.** ``PAM250`` (default) and ``BLOSUM62`` ship with MUTADOCK.
Any other matrix name from the `NCBI BLAST FTP
<https://ftp.ncbi.nih.gov/blast/matrices/>`_ (e.g. ``PAM30``, ``PAM70``,
``BLOSUM45``, ``BLOSUM80``) is downloaded automatically into ``data/`` on first
use. A fully custom matrix in NCBI format can be supplied with
``--matrix-file /path/to/matrix`` (this overrides ``--matrix``).

**When to use which.** PAM matrices model accepted point mutations over
evolutionary distance: high-numbered ones (PAM250) are permissive, low-numbered
ones (PAM30) are stringent. BLOSUM matrices are derived from conserved sequence
blocks and run the other way — low-numbered (BLOSUM45) is permissive,
high-numbered (BLOSUM80) is stringent. For reference, on a standard 20×20 alphabet
PAM250 keeps ~18% of off-diagonal substitutions while BLOSUM62 keeps ~11%.

.. list-table::
   :header-rows: 1

   * - **Goal**
     - **Suggested matrix**
   * - Broad exploration; don't miss unusual but tolerated substitutions (default)
     - ``PAM250`` (or ``BLOSUM45``)
   * - General-purpose, balanced screen
     - ``BLOSUM62``
   * - Only closely related / highly conservative substitutions; fastest run
     - ``PAM30`` / ``PAM70`` (or ``BLOSUM80``)
   * - Domain-specific scoring (e.g. membrane proteins)
     - custom file via ``--matrix-file``

Docking Studies
^^^^^^^^^^^^^^^

Docking for multiple receptors and ligands is made simple and efficient by mutadock. The text files containing the names of the receptors and ligands need to be given as input, after that everything is automated. (If md_mutate is used, the text file for receptor is generated automatically)
Every receptor in the receptor file will be docked with every ligand in the ligand file. A standard Vina configuration file or an AutoSIte prediction output is required.
Example:

.. code-block::

   md_dock -r receptors.txt -l ligands.txt -c config.txt

By default docking outputs go to an ``out/`` folder next to each receptor. Pass
``-o/--output-dir DIR`` to collect poses, logs, ``docking_results.csv``, and the
resume file in a single ``DIR`` instead. Prepared PDBQT files and AutoSite caches
still live next to their inputs so they can be reused across runs.

.. code-block::

   md_dock -r receptors.txt -l ligands.txt -c config.txt -o results/

Other optional arguments can be changed as required, to check the usage run

.. code-block::

   md_dock -h

The output of md_dock with their description is in the table below:

.. list-table::
   :header-rows: 1

   * - **No.**
     - **Output**
     - **Description**
   * - 1.
     - PDBQT files
     - The receptors and ligands will be converted to PDBQT files for AutoDock Vina.
   * - 2.
     - Output Log
     - The output of AutoDock Vina with the docking scores will be stored in a log file for each combination.
   * - 3.
     - Output PDB
     - The output of AutoDock Vina with the 5 best docking poses will be stored in a PDB file for each combination.
   * - 4.
     - Output PDBQT
     - The output of AutoDock Vina Split with the best pose will be stored in a PDBQT file for each combination.
   * - 5.
     - Output SDF
     - The best pose after docking will be stored in a SDF file for visualization and better usability.
   * - 6.
     - Docking Results CSV
     - All the docking affinities are tabulated in a CSV to make analysis trivial.


Choosing the docking box
^^^^^^^^^^^^^^^^^^^^^^^^^^

Vina only searches inside a rectangular **box**, so every docking run needs a
center and a size. MUTADOCK offers three ways to supply them, and understanding
how they interact avoids the most common surprises.

* **``-c config.txt`` — an explicit Vina config.** You provide the box and the
  search parameters directly. Recognised ``key = value`` keys are ``center_x``,
  ``center_y``, ``center_z`` (box center), ``size_x``, ``size_y``, ``size_z``
  (box dimensions in Å — these are independent, so the box may be
  **non-cubic**), plus ``exhaustiveness`` (default 32), ``n_poses`` (20),
  ``n_poses_write`` (5) and ``overwrite``. Lines beginning with ``#`` are
  ignored. Use this when you already know the pocket, want a reproducible box, or
  need to tune the search. **Watch out:** any omitted ``center_*`` key silently
  defaults to ``0.0`` — a config without an explicit center places an empty box
  at the origin, and docking will find nothing.

  .. code-block::

     center_x = 12.4
     center_y = -3.1
     center_z = 25.8
     size_x = 24
     size_y = 24
     size_z = 24
     exhaustiveness = 32

* **``-a autosite.pdb`` — a fixed AutoSite cluster.** Pass a cluster PDB that
  AutoSite already produced. MUTADOCK sets the center to the cluster's geometric
  center and the box to a **cube** of side ``2 × radius + 8 Å`` (the 8 Å is
  ``DEFAULT_BOX_MARGIN``), where *radius* is the farthest atom from that center.
  The same box is reused for every receptor.

* **Neither flag — automatic per-receptor pocket detection.** If you pass no
  config and no AutoSite PDB, MUTADOCK runs the ``autosite`` binary on each
  prepared receptor and derives the center/box the same way as ``-a`` (cached per
  receptor). This is the default for ``md_quick`` and requires the ``autosite``
  binary (from ADFRsuite) on your ``PATH``; without it, ``md_dock`` exits with an
  error asking for one of the options above. Prefer this when docking mutants
  whose pockets may shift, since each receptor gets its own box.

**How ``-c`` and ``-a`` interact.** If you pass both, AutoSite wins for the
*geometry* — the cluster's center and cube override whatever ``center_*`` /
``size_*`` were in the config — but the config's **search parameters**
(``exhaustiveness``, ``n_poses``, ``n_poses_write``, ``overwrite``) still apply.
This lets you take the box from AutoSite while keeping a tuned search from your
config.

Reports
^^^^^^^

Every ``md_mutate`` and ``md_dock`` run automatically writes a ``reports/``
folder into the run's output directory (pass ``--no-report`` to skip). It
contains a self-contained ``report.html``, a ``report.pptx``, and a
``figures/`` sub-folder with each chart as a standalone PNG. The HTML bundles a
ΔΔG distribution, a top-mutations table and chart, a docking-affinity chart and
table, and an interactive 3D view of the top mutant pose (NGL, inlined so the
file opens offline). Rebuild from an existing run directory at any time:

.. code-block::

   md_report -d results/
   md_report -d results/ --html-only --top 20
   md_report -h

Using the same ``-o/--output-dir`` for ``md_mutate`` then ``md_dock`` produces a
single combined report with both ΔΔG and docking affinities.


All CLI Scripts
^^^^^^^^^^^^^^^

.. list-table::
   :header-rows: 1

   * - **No.**
     - **Command**
     - **Description**
   * - 1.
     - md_mutate
     - Predicts the best mutation of the given protein
   * - 2.
     - md_dock
     - Docked all combinations from a list of receptors and ligands
   * - 3.
     - md_vina_dock
     - CLI for AutoDock Vina
   * - 4.
     - md_csv_generator
     - Generates all possible mutations for a protein and also the mutations possible according to PAM Matrix
   * - 5.
     - md_csv_sort
     - Can sort any CSV file according to the column name or number chosen
   * - 6.
     - md_ddg_single
     - Calculates single ddG values for a given CSV of mutations
   * - 7.
     - md_ddg_double
     - Calculates double ddG values for all combinations using a given CSV of mutations
   * - 8.
     - md_ddg_triple
     - Calculates triple ddG values for all combinations using a given CSV of mutations


Applications
------------


* **Protein Engineering:** Designing mutated proteins with enhanced stability or new functionalities.
* **Drug Discovery:** Screening potential drug candidates by predicting binding affinities.
* **Biochemical Research:** Studying protein-ligand interactions to understand biological processes.

Documentation
-------------


* README is included in the repository to serve as a comprehensive guide
* ReadtheDocs Page for updated documentation can be found `here <https://mutadock.readthedocs.io/en/latest/#>`_

Future Developments
-------------------


* **Developing a Graphical User Interface (GUI):** Enhancing user experience by providing a user-friendly interface for easier interaction with the software.
* **Creating a Web Server:** Allowing remote access and usage of the software through a web-based platform, making it accessible from anywhere.
* **Increasing Parameter Customizability:** Offering more options for users to fine-tune mutation and docking parameters to suit specific research needs and conditions.

Acknowledgements
----------------


* Open-source tools and libraries used in the development.

Contact
-------


* For questions, suggestions, or collaboration, please contact `Naisarg Patel <mailto:naisarg.patel14@hotmail.com>`_.
