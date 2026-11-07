Agent skills
============

MUTADOCK ships **agent skills**: focused Markdown guides that let a coding agent
drive the library without re-deriving the CLI from this documentation. Each one
is a ``SKILL.md`` file following the `Agent Skills specification
<https://agentskills.io/specification>`_.

.. list-table::
   :header-rows: 1
   :widths: 30 70

   * - Skill
     - Use it for
   * - ``mutadock-setup``
     - Installing or verifying PyRosetta, AutoDock Vina, meeko, pdbfixer/openmm, AutoSite
   * - ``mutadock-quickstart``
     - A one-command end-to-end run with ``md_quick``
   * - ``mutadock-mutation``
     - Planning and running |ddG| scans; writing mutant PDB files
   * - ``mutadock-docking``
     - Batch docking, search boxes, reproducible seeds
   * - ``mutadock-results``
     - Reading the output CSVs and reports; computing |ddG|\ :sub:`binding`
   * - ``mutadock-debugging``
     - Chain, residue, ligand and structure input failures
   * - ``mutadock-python-api``
     - Scripting MUTADOCK from Python

.. |ddG| replace:: ΔΔG

Working in a clone
------------------

Start the agent at the repository root. Claude Code reads ``CLAUDE.md`` and
``.claude/skills/``; Codex reads ``AGENTS.md`` and ``.agents/skills/``. Both
trees are committed, so discovery is automatic:

.. code-block:: bash

   git clone https://github.com/naisarg14/mutadock && cd mutadock
   claude    # or: codex

Ask in plain language ("plan a ΔΔG scan around residue A386 of 4QJR"), or invoke
a skill explicitly with ``/mutadock-mutation`` (Claude Code) or
``$mutadock-mutation`` (Codex).

Working from a pip install
--------------------------

The skills ship inside the wheel, so they can be installed without cloning:

.. code-block:: bash

   md_skills list                             # show the bundled skills
   md_skills install --claude --codex --user  # ~/.claude/skills and ~/.agents/skills
   md_skills install --project ~/my-study     # that project only

Only directories named ``mutadock-*`` are written, so a shared skills folder
containing other skills is left untouched. Pass ``--force`` to replace an
existing installation.

Skills provide guidance, not tooling: running a workflow still requires a
working MUTADOCK environment.

Contributing a change
---------------------

``.agents/skills/`` is the source of truth; ``.claude/skills/`` is a generated
mirror (plain copies, not symlinks, so Windows checkouts work). After editing a
skill:

.. code-block:: bash

   python scripts/sync_agent_skills.py
   pytest tests/test_skills

The test suite checks the frontmatter, that every skill is routable from
``INDEX.md``, that the two trees match, and that no skill names an ``md_*``
command absent from ``[project.scripts]``.
