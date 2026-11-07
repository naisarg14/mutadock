# `md_mutate` output files

Written next to the input structure, or into `-o/--output-dir` when given.
`protein` below is the input stem; `_modified` appears when the input was
converted or cleaned.

| # | File | Description |
|---|---|---|
| 1 | `protein_modified_mutations_all.csv` | Every possible single-residue substitution with its matrix score |
| 2 | `protein_modified_mutations.csv` | Substitutions with a positive matrix score (the ΔΔG input set) |
| 3 | `protein_modified_mutations_ddG.csv` | ΔΔG for each row of file 2 |
| 4 | `protein_modified_mutations_ddG_sorted.csv` | File 3, lowest → highest ΔΔG (most stabilizing first) |
| 5 | `protein_modified_double_ddg.csv` | Double-mutation ΔΔG combinations |
| 6 | `protein_modified_double_ddg_sorted.csv` | File 5, sorted |
| 7 | `protein_modified_triple_ddg.csv` | Triple-mutation ΔΔG combinations |
| 8 | `protein_modified_triple_ddg_sorted.csv` | File 7, sorted |
| 9 | `protein_modified_mutants.txt` | Paths of the generated mutant PDBs — direct input for `md_dock -r` |

Plus `protein_clean.pdb` (the cleaned structure every ΔΔG run scores against),
the mutant PDB files themselves, and a `reports/` bundle unless `--no-report`.

Mutant PDBs from `md_mutate` are named `{sr}_{WT}-{CHAIN}{POS}-{NEW}.pdb` for
singles (e.g. `1_ASN-A386-HIS.pdb`) and `{sr}_{combination}.pdb` for doubles and
triples (e.g. `10_80_N_386_Q+H_387_E.pdb`). `md_generate_pdb` uses a different
pattern — see the main skill.

## Intermediate file reuse and resume

`md_mutate` checks each output before recomputing it, and the double/triple
steps compare the file's **line count** against the number of combinations the
current `-s`/`-t` imply. A file produced with different `-s`/`-t` is therefore
recomputed rather than silently reused at the wrong size.

Re-running the same command after an interruption resumes. `--no-append` forces
a fresh computation of everything.
