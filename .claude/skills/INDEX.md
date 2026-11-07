# MUTADOCK skills index

Route here first, then open the single `<name>/SKILL.md` that matches the work in
hand. Do not survey this directory or read several skills to decide which one
applies — each restates the rules it needs, so the matching one is sufficient.

| Skill | Open it when you need to |
|---|---|
| `mutadock-setup` | Install or verify PyRosetta, AutoDock Vina, meeko, pdbfixer/openmm, AutoSite. |
| `mutadock-quickstart` | Run the one-shot `md_quick` demo, or confirm an installation end to end. |
| `mutadock-mutation` | Plan or run a ΔΔG mutation scan, or build specific mutant PDB files. |
| `mutadock-docking` | Batch-dock receptors × ligands, define the search box, reproduce a run. |
| `mutadock-results` | Read the output CSVs and reports, compute ΔΔG_binding, apply the caveats. |
| `mutadock-debugging` | Fix input failures: chain IDs, residue numbering, ligands, CIF/PDB parsing. |
| `mutadock-python-api` | Script MUTADOCK from Python instead of the CLI. |

Project-wide rules that every skill assumes live in `AGENTS.md` at the repository
root. Skills installed standalone (via `md_skills install`) restate the rules they
depend on, so they work without it.
