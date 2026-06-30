# imgt_renumber

Renumber **TCR** (and antibody) PDB/mmCIF files with **IMGT** numbering, without the
usual ANARCI install pain. The classic ANARCI needs a system HMMER binary and a
network build step; this tool stands on a fully **pip-installable** stack instead:

- `pyhmmer` — bundles HMMER3, no separate binary
- `biopython` — structure parsing/writing
- ANARCI (the `prihoda` pyhmmer fork) + its prebuilt IMGT HMMs

ANARCI does the actual numbering, so results match ANARCI/IMGT exactly. This script
adds the part ANARCI doesn't give you: extracting chain sequences from the structure,
running the numbering, and mapping IMGT numbers + insertion codes back onto the
residues of a renumbered output file.

## Install (one time)

```bash
python imgt_renumber.py setup
```

This installs the backend, copies the prebuilt IMGT HMMs into place, and applies a
tiny str/bytes compatibility shim so it works with current `pyhmmer`. It's idempotent
— safe to re-run. (Needs `pip` and `git`, and network access to PyPI + GitHub.)

## Use

```bash
# basic: writes <input>_imgt.pdb
python imgt_renumber.py number 1ao7.pdb

# choose the output name
python imgt_renumber.py number 1ao7.pdb -o 1ao7_imgt.pdb

# mmCIF in/out works too (by file extension)
python imgt_renumber.py number model.cif -o model_imgt.cif
```

Options for the `number` command:

| flag | default | meaning |
|------|---------|---------|
| `-o, --output` | `<input>_imgt.<ext>` | output path |
| `--scheme` | `imgt` | `imgt` or `aho` (the two TCR-valid schemes) |
| `--species` | all | restrict, e.g. `--species human,mouse` |
| `--chains` | all | restrict, e.g. `--chains D,E` |
| `--constant` | `renumber` | how to treat non-variable-domain residues (see below) |
| `--bit-score-threshold` | `80` | ANARCI acceptance threshold |

### What gets renumbered

Only the **variable domain** carries true IMGT numbering (positions 1–128, with CDR3
insertion codes placed symmetrically about 111/112, e.g. `111A 111B … 112B 112A`).
IMGT has no V-domain-style scheme for the constant domain, so `--constant` controls
the rest of each chain:

- `renumber` (default): non-V residues are numbered sequentially **after** the V
  domain (129, 130, …), guaranteeing unique, ordered IDs.
- `original`: keep their original author numbering (errors out if that collides with
  the new V-domain IDs).
- `drop`: remove non-V residues entirely.

Non-amino-acid records (waters, ligands) are left untouched.

> Note: with `renumber`, any residues *N-terminal* to the V domain (e.g. an expression
> tag still present in the coordinates) also land after 128. For the usual V-then-C
> chain layout this is exactly what you want.

## Use as a library

```python
from Bio.PDB import PDBParser, PDBIO
from imgt_renumber import renumber_structure

s = PDBParser(QUIET=True).get_structure("x", "1ao7.pdb")
report = renumber_structure(s, scheme="imgt", constant="renumber")
# report -> [{'chain': 'D', 'chain_type': 'A', 'label': 'TCR alpha (Valpha)',
#             'species': 'human', 'residues_numbered': 114, 'bitscore': 123.2, ...}, ...]
io = PDBIO(); io.set_structure(s); io.save("1ao7_imgt.pdb")
```

## Validation

Tested on αβ TCR chains: conserved IMGT anchors land exactly (1st-CYS at 23,
CONSERVED-TRP at 41, 2nd-CYS at 104, J-PHE/TRP at 118), and long CDR3s receive the
correct symmetric insertion codes around 111/112.
