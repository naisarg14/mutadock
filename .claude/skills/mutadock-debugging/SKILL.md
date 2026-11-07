---
name: mutadock-debugging
description: Diagnose MUTADOCK input failures — structure warnings about multiple models, multiple chains, similar chains, altlocs and insertion codes; residue numbering and chain-ID mismatches between a PDB and the PyRosetta pose; ligand SDF rejections for missing 3D coordinates or hydrogens; CIF/PDB parsing errors; mutations that silently apply to the wrong residue. Use when an md_* command errors, warns, or produces results that do not match the structure you expected.
---

# Debugging MUTADOCK inputs

Most MUTADOCK failures are input problems that surface late. Read the warnings
the tools print at the start of a run — they are the cheap diagnosis.

## Structure warnings from `structure_warnings`

`md_mutate` and `md_quick` call this on the input and print any of the
following. None aborts the run; all of them change what the results mean.

**"Structure contains N models (e.g. an NMR ensemble)"**
Only one model is used after cleaning. The rest of the pipeline silently scores
a single arbitrary conformer. Extract the model you want first:

```bash
python - <<'PY'
from Bio.PDB import PDBParser, PDBIO, Select
s = PDBParser(QUIET=True).get_structure("x", "protein.pdb")
io = PDBIO(); io.set_structure(s[0])          # model 0
io.save("protein_model1.pdb")
PY
```

**"Structure contains alternate location indicators (A, B, …)"**
Alternate side-chain conformations change the packing environment, so ΔΔG
depends on which one survived cleaning. Keep altloc A (or the highest-occupancy
one) explicitly before scoring if the site of interest has altlocs.

**"Structure contains N residue(s) with insertion codes"**
Insertion codes (e.g. `100A`, `100B`) make `position` ambiguous: a mutation
spec of `A:100:ASN:HIS` cannot distinguish them. Renumber the structure, or
work on a chain without insertion codes.

**"Structure contains N chains (…); every chain is mutated/scored"**
This is the most expensive warning. A homodimer doubles the ΔΔG workload and
produces two rows for every real position. Isolate the chain you care about:

```bash
python - <<'PY'
from Bio.PDB import PDBParser, PDBIO, Select
class OneChain(Select):
    def accept_chain(self, c): return c.id == "A"
s = PDBParser(QUIET=True).get_structure("x", "protein.pdb")
io = PDBIO(); io.set_structure(s); io.save("protein_A.pdb", OneChain())
PY
```

**"Chains X, Y are similar in the input (e.g. a homodimer)"**
Confirms the duplication above. Deduplicate before interpreting the hit list,
or the same substitution will appear to be found twice.

## Chain and residue numbering

**A mutation applies to the wrong residue, or errors as not found.** PDB author
numbering is not PyRosetta pose numbering. `pose_residue_map` builds the
mapping; use it to check what a `CHAIN:POSITION` spec actually resolves to
before scoring a long run:

```python
from pyrosetta import init, pose_from_pdb
from mutadock.mutation.helpers import pose_residue_map

init("-mute all")
pose = pose_from_pdb("protein_clean.pdb")   # takes a POSE, not a path
numbering = pose_residue_map(pose)          # {chain: [residue numbers]}
print("chains:", sorted(numbering))
print("chain A spans:", min(numbering["A"]), "-", max(numbering["A"]))
print("386 present:", 386 in numbering["A"])
```

The numbers come from the pose's PDB info, so they are the author numbering a
`CHAIN:POSITION` mutation spec uses. When a residue really is absent,
`format_missing_residue` turns that into an actionable message naming the
available chains or the chain's numbering span.

Check against the **cleaned** PDB, not the original: `clean_pdb` removes HETATM
records and non-standard residues, which shifts pose indices.

**ΔΔG scoring expects `*_clean.pdb`.** `md_ddg_single|double|triple -p` must
receive the cleaned structure `md_mutate` produces. Build it yourself when
skipping `md_mutate`:

```bash
python -c "from mutadock.mutation.helpers import clean_pdb; clean_pdb('protein.pdb','protein_clean.pdb')"
```

Remember that cleaning **strips HETATM**: cofactors, metals, waters and ligands
are not present during ΔΔG scoring.

## Ligand failures

**"ligand has no 3D coordinates" / missing hydrogens.** `md_dock` and `md_quick`
validate ligands up front (`_validate_sdf`, `_require_3d_ligand`). A 2D SDF
drawn in a sketcher, or one exported without explicit hydrogens, is rejected.
Generate a 3D conformer first:

```bash
python - <<'PY'
from rdkit import Chem
from rdkit.Chem import AllChem
m = Chem.AddHs(Chem.MolFromMolFile("ligand_2d.sdf"))
AllChem.EmbedMolecule(m, randomSeed=19)
AllChem.MMFFOptimizeMolecule(m)
Chem.MolToMolFile(m, "ligand_3d.sdf")
PY
```

**A PubChem `--ligand-code` lookup fails.** Disambiguate the form explicitly:
`cid:5291` or `name:aspirin` instead of a bare token. Fall back to
`--ligand-file` with a locally prepared SDF.

**`mk_receptor: command not found`.** meeko is missing from the active
environment: `pip install meeko`. See `mutadock-setup`.

## Structure file parsing

**"CIF file not recognized".** Only `.pdb` and `.cif` extensions are accepted —
a correct file with a different extension is rejected. CIF inputs are converted
up front (`convert_cif_pdb`), by both the CLIs and `generate_csv`, so a `.cif`
also leaves a `.pdb` next to it. Downstream files carry the converted stem,
which is why output names can differ from the name you passed in.

**`ValueError: invalid literal for int() with base 10` from a structure call.**
A CIF reached a PDB-only parser. `generate_csv` converts for you; other helpers
that take a structure path may not, so pass `convert_cif_pdb(path)` when a
function's docstring says PDB.

**PDBFixer errors during receptor preparation.** Missing residues, nonstandard
residues or broken chain breaks. Check the structure opens in BioPython first,
and consider isolating the chain of interest.

## Resume behaving unexpectedly

**A rerun recomputes everything.** `--no-append` (`md_mutate`) and
`-i/--ignore_existing` (`md_dock`) force fresh work — drop them to resume.
`md_ddg_*` needs `--resume` explicitly.

**A rerun skips work you wanted redone.** `md_dock` skips pairs in
`*_completed.txt`; delete that file or pass `-i`. `md_ddg_*` `--resume` matches
on the recorded `ddG_protocol` and `n_replicates`, so changing protocol flags
correctly forces recomputation rather than mixing incompatible numbers.

**The double/triple step recomputes after you changed `-s`/`-t`.** Expected:
the step compares the existing file's row count against the number of
combinations the current flags imply, so a file built at a different size is
rebuilt instead of being reused at the wrong size.

## Timeouts

Defaults are 900 s for receptor preparation, 1800 s for AutoSite and 3600 s for
Vina. Raise them with `MUTADOCK_RECEPTOR_PREP_TIMEOUT`,
`MUTADOCK_AUTOSITE_TIMEOUT`, `MUTADOCK_VINA_TIMEOUT`; `0` disables a timeout.

## AutoSite

`FileExistsError` on `{stem}_autosite_out` is deliberate — the directory is
created with `exist_ok=False` so a stale run's grids are never silently reused.
Delete it, or skip AutoSite with `-c config.txt` / `-a pocket.pdb`.

## Matrix download

PAM250 and BLOSUM62 ship with the package. Others are fetched from the NCBI
BLAST FTP into `~/.cache/mutadock/matrices` (`MUTADOCK_DATA_DIR` relocates it).
If NCBI is unreachable, download the file by hand and pass
`--matrix-file /path/to/PAM30`.
