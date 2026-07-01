# imgt_renumber

Renumber **TCR** (and antibody) PDB/mmCIF files with **IMGT** numbering, without the
usual ANARCI install pain. The classic ANARCI needs a system HMMER binary and a
network build step; this tool stands on a fully **pip-installable** stack instead:

- `pyhmmer` — bundles HMMER3, no separate binary
- `biopython` — structure parsing/writing
- ANARCI (a pinned [pyhmmer fork](https://github.com/LilySnow/ANARCI-pyhmmer-fork)) + its prebuilt IMGT HMMs

ANARCI does the actual numbering, so results match ANARCI/IMGT exactly. This tool
adds the part ANARCI doesn't give you: extracting chain sequences from the structure,
running the numbering, and mapping IMGT numbers + insertion codes back onto the
residues of a renumbered output file.

## Install

Install the package from GitHub, then run the one-time backend setup:

```bash
pip install "git+https://github.com/X-lab-3D/imgt_renumber.git"
imgt-renumber setup
imgt-renumber number 1ao7.pdb -o 1ao7_imgt.pdb
```

Using a virtual environment (or conda env) is recommended, since `setup` installs and
patches the ANARCI backend inside whatever Python is active:

```bash
python -m venv tcr-imgt && source tcr-imgt/bin/activate   # Windows: tcr-imgt\Scripts\activate
pip install "git+https://github.com/LilySnow/imgt-renumber.git"
imgt-renumber setup
```

`setup` installs the backend, copies the prebuilt IMGT HMMs into place (from the pinned
fork), and applies a tiny str/bytes compatibility shim so it works with current
`pyhmmer`. It is idempotent — safe to re-run. It needs `pip` and `git`, plus network
access to PyPI and GitHub.

> Not on PyPI: because it depends on a not-on-PyPI ANARCI fork (a direct URL
> dependency), this package installs from GitHub rather than via `pip install
> imgt-renumber`.

### Without installing the package

You can also just run the single script directly — no `pip install .` needed:

```bash
python imgt_renumber.py setup
python imgt_renumber.py number 1ao7.pdb -o 1ao7_imgt.pdb
```

## Use

```bash
# basic: writes <input>_imgt.pdb
imgt-renumber number 1ao7.pdb

# choose the output name
imgt-renumber number 1ao7.pdb -o 1ao7_imgt.pdb

# mmCIF in/out works too (by file extension)
imgt-renumber number model.cif -o model_imgt.cif
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

## Troubleshooting

**`error: externally-managed-environment` (pip refuses to install).**
Modern Debian/Ubuntu Python blocks installs into the system interpreter. Use a
virtual environment (see Install above) — that's the clean fix. `setup` also retries
with `--break-system-packages` automatically as a fallback.

**Network / GitHub errors during `setup`.**
`setup` needs to reach PyPI and GitHub to fetch the backend and the IMGT HMMs. Behind
a proxy or firewall, make sure `pip` and `git` can reach `github.com` and
`pypi.org`. Fully offline installs aren't supported out of the box (you'd need to
pre-stage the ANARCI wheel and the `dat/` HMM folder).

**`ANARCI backend not installed` when running `number`.**
Run `imgt-renumber setup` first (once per environment).

**No variable domains were numbered.**
The chain(s) may be below the detection threshold, or not antigen-receptor variable
domains. Try lowering `--bit-score-threshold`, or check you're pointing at the right
chains with `--chains`.

## License

The code in this repository (`imgt_renumber.py`) is released under the Apache
License 2.0 — see `LICENSE` and `NOTICE`.

This repository does **not** bundle any third-party source or data; the dependencies
below are fetched by pip / the `setup` step at install time and remain under their own
licenses. All of them are permissive (no copyleft), so this tool can be shared and
relicensed freely.

## Third-party components

| Component | Role | License |
|-----------|------|---------|
| [ANARCI](https://github.com/oxpig/ANARCI), used via the pinned fork [`LilySnow/ANARCI-pyhmmer-fork`](https://github.com/LilySnow/ANARCI-pyhmmer-fork) (forked from [`prihoda/ANARCI`](https://github.com/prihoda/ANARCI)) | performs the actual IMGT numbering; its prebuilt HMMs are copied in by `setup` | BSD 3-Clause |
| MUSCLE | used only in ANARCI's HMM build pipeline (not at runtime) | Public domain |
| [pyhmmer](https://github.com/althonos/pyhmmer) (bundles HMMER3 + Easel) | sequence/HMM alignment engine | MIT / BSD-3-Clause / BSD-2-Clause |
| [Biopython](https://biopython.org/) | structure parsing and writing | Biopython License Agreement (BSD-style) |

The ANARCI dependency is pinned to commit `4995e15` of the fork for reproducibility.
If you ever redistribute ANARCI's source or its HMM data files alongside this tool
(rather than fetching them at install time), retain ANARCI's BSD-3-Clause copyright
notice and disclaimer.

## Citing

If you use this in published work, please cite the tools that do the heavy lifting:

- Dunbar J, Deane CM. *ANARCI: antigen receptor numbering and receptor
  classification.* Bioinformatics, 2016.
- Larralde M, Zeller G. *PyHMMER: a Python library binding to HMMER.*
  Bioinformatics, 2023.
- The IMGT unique numbering: Lefranc M-P et al., *Dev. Comp. Immunol.*, 2003.

