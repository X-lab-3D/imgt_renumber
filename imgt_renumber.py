#!/usr/bin/env python3
# SPDX-License-Identifier: Apache-2.0
# Copyright (c) 2026 <Li Xue>
"""
imgt_renumber.py - Renumber TCR (and antibody) PDB/mmCIF files with IMGT numbering.

It uses ANARCI for the actual numbering, but unlike the classic ANARCI install
(which needs a system HMMER binary and a network build step), this script relies
on a fully pip-installable stack:

    pyhmmer            (bundles HMMER3, no system binary)
    biopython          (structure parsing / writing)
    ANARCI (prihoda fork, the pyhmmer branch)  + its bundled IMGT HMMs

Run the one-time setup, then renumber:

    python imgt_renumber.py setup
    python imgt_renumber.py number 1ao7.pdb -o 1ao7_imgt.pdb

`setup` installs the backend, copies the prebuilt IMGT HMMs into place, and
applies a tiny str/bytes compatibility patch so it works with current pyhmmer.
It is idempotent - safe to run repeatedly.

You can also import and call renumber_structure(...) from your own code.
"""

import argparse
import os
import subprocess
import sys

# ANARCI chain_type -> human-readable label (TCR + antibody)
CHAIN_TYPE_LABELS = {
    "A": "TCR alpha (Valpha)",
    "B": "TCR beta (Vbeta)",
    "G": "TCR gamma (Vgamma)",
    "D": "TCR delta (Vdelta)",
    "H": "antibody heavy (VH)",
    "K": "antibody kappa (VL-kappa)",
    "L": "antibody lambda (VL-lambda)",
}

ANARCI_PIN = "4995e15925f6584115b76ee790c19a7468fbcc8a"  # pinned commit for reproducibility
ANARCI_REPO = "https://github.com/LilySnow/ANARCI-pyhmmer-fork.git"
ANARCI_GIT = f"git+{ANARCI_REPO}@{ANARCI_PIN}"


# --------------------------------------------------------------------------- #
#  SETUP / BACKEND BOOTSTRAP
# --------------------------------------------------------------------------- #
def _pip(*pkgs):
    cmd = [sys.executable, "-m", "pip", "install", "--break-system-packages", *pkgs]
    # --break-system-packages is harmless on normal/venv installs; drop it if unsupported
    try:
        subprocess.check_call(cmd)
    except subprocess.CalledProcessError:
        subprocess.check_call([sys.executable, "-m", "pip", "install", *pkgs])


def _patch_decode(anarci_py: str) -> bool:
    """Make ANARCI's .decode() sites tolerant of pyhmmer returning str.
    Idempotent: a sentinel marker prevents re-patching. Returns True if changed."""
    SENTINEL = "# imgt_renumber: str/bytes decode patch applied\n"
    with open(anarci_py) as fh:
        src = fh.read()
    if SENTINEL in src:
        return False  # already patched
    reps = [
        ("hsp.alignment.hmm_name.decode()",
         "(hsp.alignment.hmm_name.decode() if isinstance(hsp.alignment.hmm_name, bytes) "
         "else hsp.alignment.hmm_name)"),
        ("ali.hmm_name.decode(), ali.hmm_accession.decode()",
         "(ali.hmm_name.decode() if isinstance(ali.hmm_name, bytes) else ali.hmm_name), "
         "(ali.hmm_accession.decode() if isinstance(ali.hmm_accession, bytes) "
         "else ali.hmm_accession)"),
        ("ali.hmm_name.decode().split('_')",
         "(ali.hmm_name.decode() if isinstance(ali.hmm_name, bytes) "
         "else ali.hmm_name).split('_')"),
        ("int(hit.accession.decode())",
         "int(hit.accession.decode() if isinstance(hit.accession, bytes) else hit.accession)"),
    ]
    changed = False
    for old, new in reps:
        if old in src:
            src = src.replace(old, new)
            changed = True
    if changed:
        src = SENTINEL + src
        with open(anarci_py, "w") as fh:
            fh.write(src)
    return changed


def _ensure_hmms(anarci_pkg_dir: str):
    """Copy the prebuilt IMGT HMMs from the repo into the installed package if missing."""
    import shutil
    import tempfile

    target = os.path.join(anarci_pkg_dir, "dat", "HMMs", "ALL.hmm")
    if os.path.exists(target):
        return
    print("  IMGT HMM database missing - fetching prebuilt HMMs from the pinned fork...")
    tmp = tempfile.mkdtemp(prefix="anarci_hmm_")
    src = os.path.join(tmp, "src")
    # Shallow-fetch exactly the pinned commit so HMMs match the installed code.
    subprocess.check_call(["git", "init", "-q", src])
    subprocess.check_call(["git", "-C", src, "remote", "add", "origin", ANARCI_REPO])
    subprocess.check_call(["git", "-C", src, "fetch", "-q", "--depth", "1", "origin", ANARCI_PIN])
    subprocess.check_call(["git", "-C", src, "checkout", "-q", "FETCH_HEAD"])
    src_dat = os.path.join(src, "lib", "python", "anarci", "dat")
    if not os.path.exists(os.path.join(src_dat, "HMMs", "ALL.hmm")):
        raise RuntimeError("Prebuilt HMMs not found in the ANARCI repo - layout may have changed.")
    dst_dat = os.path.join(anarci_pkg_dir, "dat")
    if os.path.exists(dst_dat):
        shutil.rmtree(dst_dat)
    shutil.copytree(src_dat, dst_dat)
    shutil.rmtree(tmp, ignore_errors=True)


def _anarci_pkg_dir():
    """Return the installed anarci package dir, or None if not properly installed.
    Uses a fresh interpreter so in-process import caching can't interfere."""
    code = ("import anarci, os, sys;"
            "f=getattr(anarci,'__file__',None);"
            "sys.stdout.write(os.path.dirname(f) if f else '')")
    try:
        out = subprocess.check_output([sys.executable, "-c", code],
                                      stderr=subprocess.DEVNULL).decode().strip()
    except subprocess.CalledProcessError:
        return None
    # A bare namespace-package leftover (e.g. a stray dat/ dir) has no anarci.py.
    if out and os.path.exists(os.path.join(out, "anarci.py")):
        return out
    return None


def run_setup():
    print("[1/4] Installing pyhmmer + biopython ...")
    _pip("pyhmmer", "biopython")

    print("[2/4] Installing ANARCI (pyhmmer fork) ...")
    if _anarci_pkg_dir():
        print("  anarci already installed.")
    else:
        _pip(ANARCI_GIT)

    pkg_dir = _anarci_pkg_dir()
    if not pkg_dir:
        raise RuntimeError(
            "ANARCI did not install correctly. Try:\n"
            f"  pip install '{ANARCI_GIT}'"
        )

    print("[3/4] Ensuring IMGT HMM database is present ...")
    _ensure_hmms(pkg_dir)

    print("[4/4] Applying pyhmmer str/bytes compatibility patch ...")
    changed = _patch_decode(os.path.join(pkg_dir, "anarci.py"))
    print("  patched." if changed else "  already patched.")

    # Verify end to end
    print("\nVerifying backend on a test TCR beta chain ...")
    # reimport in a clean subprocess so the patched module is loaded fresh
    test = (
        "from anarci import anarci;"
        "seq=[('t','DAGVIQSPRHEVTEMGQEVTLRCKPISGHNSLFWYRQTMMRGLELLIYFNNNVPIDDSGMPEDRFSAKMPNASFSTLKIQPSEPRDSAVYFCASSYVGNTGELFFGEGSRLTVL')];"
        "n,d,h=anarci(seq,scheme='imgt',output=False);"
        "print('  OK -> chain_type=%s species=%s score=%.1f' % "
        "(d[0][0]['chain_type'], d[0][0]['species'], d[0][0]['bitscore']))"
    )
    subprocess.check_call([sys.executable, "-c", test])
    print("\nSetup complete. You can now run:  python imgt_renumber.py number IN.pdb -o OUT.pdb")


# --------------------------------------------------------------------------- #
#  STRUCTURE I/O HELPERS
# --------------------------------------------------------------------------- #
def _load_structure(path):
    from Bio.PDB import PDBParser, MMCIFParser

    ext = os.path.splitext(path)[1].lower()
    if ext in (".cif", ".mmcif"):
        parser = MMCIFParser(QUIET=True)
        fmt = "cif"
    else:
        parser = PDBParser(QUIET=True)
        fmt = "pdb"
    structure = parser.get_structure("input", path)
    return structure, fmt


def _save_structure(structure, path, fmt):
    from Bio.PDB import PDBIO, MMCIFIO

    ext = os.path.splitext(path)[1].lower()
    use_cif = ext in (".cif", ".mmcif") or (fmt == "cif" and ext not in (".pdb", ".ent"))
    io = MMCIFIO() if use_cif else PDBIO()
    io.set_structure(structure)
    io.save(path)


def _three_to_one():
    from Bio.Data.PDBData import protein_letters_3to1
    return protein_letters_3to1


def _chain_sequence(chain, aa_map):
    """Return (seq_string, [residue,...]) for standard amino-acid residues in order."""
    from Bio.PDB.Polypeptide import is_aa

    seq, residues = [], []
    for res in chain:
        # standard=False so modified residues are still considered amino acids
        if not is_aa(res, standard=False):
            continue
        one = aa_map.get(res.get_resname(), "X")
        seq.append(one)
        residues.append(res)
    return "".join(seq), residues


# --------------------------------------------------------------------------- #
#  CORE RENUMBERING
# --------------------------------------------------------------------------- #
def renumber_structure(structure, scheme="imgt", allowed_species=None,
                       bit_score_threshold=80, constant="renumber",
                       restrict_chains=None, verbose=True):
    """Renumber variable domains in `structure` in place. Returns a report list.

    constant: what to do with residues that are NOT part of a numbered V domain
        'renumber' -> continue sequentially after the last IMGT number (unique, ordered)
        'original' -> keep their original author numbering (may collide; a warning is given)
        'drop'     -> remove them from the structure
    """
    from anarci import anarci

    aa_map = _three_to_one()
    report = []

    # Number using model 0 sequences, then apply identical mapping to every model.
    model0 = next(structure.get_models())

    for chain in model0:
        cid = chain.id
        if restrict_chains and cid not in restrict_chains:
            continue
        seq, residues = _chain_sequence(chain, aa_map)
        if len(seq) < 30:
            continue  # too short to be a variable domain

        results = anarci([(cid, seq)], scheme=scheme,
                         allowed_species=allowed_species,
                         bit_score_threshold=bit_score_threshold, output=False)
        numbering, details, _ = results
        if numbering[0] is None:
            if verbose:
                print(f"  chain {cid}: no Ig/TCR variable domain detected - left unchanged")
            continue

        # Build the new id for every amino-acid residue index in this chain.
        # Default: keep original; we overwrite the V-domain region below.
        new_ids = {i: None for i in range(len(residues))}
        domain_info = []

        for dom_num, dom_det in zip(numbering[0], details[0]):
            dom_states = dom_num[0]            # list of ((num, icode), aa)
            qstart = dom_det["query_start"]    # 0-based index of first V residue
            qend = dom_det["query_end"]        # exclusive end (one past last V residue)
            ctype = dom_det["chain_type"]
            species = dom_det["species"]

            seq_idx = qstart
            assigned = 0
            for (num, icode), aa in dom_states:
                if aa == "-":
                    continue
                ic = icode if icode and icode != " " else " "
                new_ids[seq_idx] = (" ", num, ic)
                seq_idx += 1
                assigned += 1
            if seq_idx != qend:
                raise RuntimeError(
                    f"chain {cid}: alignment/sequence mismatch "
                    f"(consumed up to {seq_idx}, expected {qend})"
                )
            domain_info.append((ctype, species, qstart, qend - 1, assigned,
                                dom_det.get("bitscore")))

        if not domain_info:
            continue

        # Determine the highest IMGT number used, for sequential constant numbering.
        max_imgt = max((v[1] for v in new_ids.values() if v is not None), default=0)
        running = max_imgt + 1

        # Fill in residues not covered by any V domain.
        for i in range(len(residues)):
            if new_ids[i] is not None:
                continue
            if constant == "renumber":
                new_ids[i] = (" ", running, " ")
                running += 1
            elif constant == "original":
                new_ids[i] = residues[i].id  # keep as-is
            elif constant == "drop":
                new_ids[i] = "DROP"
            else:
                raise ValueError(f"Unknown constant mode: {constant}")

        _apply_chain_numbering(structure, cid, new_ids, constant)

        # Report
        for ctype, species, qs, qe, assigned, score in domain_info:
            label = CHAIN_TYPE_LABELS.get(ctype, f"type {ctype}")
            report.append({
                "chain": cid, "chain_type": ctype, "label": label,
                "species": species, "v_start_resindex": qs, "v_end_resindex": qe,
                "residues_numbered": assigned,
                "bitscore": round(score, 1) if score is not None else None,
            })
            if verbose:
                sc = f"{score:.1f}" if score is not None else "?"
                print(f"  chain {cid}: {label}, species={species}, "
                      f"score={sc}, V-domain residues={assigned}")

    return report


def _apply_chain_numbering(structure, chain_id, new_ids, constant):
    """Apply new_ids (index -> new residue id or 'DROP') to `chain_id` in all models."""
    from Bio.PDB.Polypeptide import is_aa

    for model in structure:
        if chain_id not in model:
            continue
        chain = model[chain_id]
        aa_residues = [r for r in chain if is_aa(r, standard=False)]
        if len(aa_residues) != len(new_ids):
            raise RuntimeError(
                f"chain {chain_id}: residue count differs across models "
                f"({len(aa_residues)} vs {len(new_ids)}); cannot map safely."
            )

        # Drops first.
        to_drop = [aa_residues[i] for i, nid in new_ids.items() if nid == "DROP"]
        for res in to_drop:
            chain.detach_child(res.id)
        keep = [(aa_residues[i], new_ids[i]) for i in range(len(aa_residues))
                if new_ids[i] != "DROP"]

        # Two-pass reassignment via temporary ids to avoid transient collisions.
        for j, (res, _nid) in enumerate(keep):
            res.id = (" ", 900000 + j, " ")
        for res, nid in keep:
            res.id = nid

        # Final uniqueness check & child_dict rebuild.
        seen = {}
        for res in chain:
            if res.id in seen:
                raise RuntimeError(
                    f"chain {chain_id}: duplicate residue id {res.id} after renumbering "
                    f"(try --constant renumber)."
                )
            seen[res.id] = res
        chain.child_dict = {r.id: r for r in chain}


# --------------------------------------------------------------------------- #
#  CLI
# --------------------------------------------------------------------------- #
def run_number(args):
    try:
        import anarci  # noqa: F401
    except ImportError:
        print("ANARCI backend not installed. Run:  python imgt_renumber.py setup",
              file=sys.stderr)
        sys.exit(1)

    structure, fmt = _load_structure(args.input)
    out = args.output or _default_out(args.input)
    species = [s.strip() for s in args.species.split(",")] if args.species else None
    chains = set(c.strip() for c in args.chains.split(",")) if args.chains else None

    print(f"Renumbering {args.input}  (scheme={args.scheme}, constant={args.constant})")
    report = renumber_structure(
        structure, scheme=args.scheme, allowed_species=species,
        bit_score_threshold=args.bit_score_threshold, constant=args.constant,
        restrict_chains=chains, verbose=True,
    )
    if not report:
        print("No variable domains were numbered - output not written.", file=sys.stderr)
        sys.exit(2)

    _save_structure(structure, out, fmt)
    print(f"\nWrote {out}  ({len(report)} domain(s) numbered)")


def _default_out(inp):
    base, ext = os.path.splitext(inp)
    return f"{base}_imgt{ext or '.pdb'}"


def build_parser():
    p = argparse.ArgumentParser(
        description="Renumber TCR/antibody PDB or mmCIF files with IMGT numbering (via ANARCI/pyhmmer).")
    sub = p.add_subparsers(dest="cmd", required=True)

    sub.add_parser("setup", help="One-time install/patch of the numbering backend.")

    n = sub.add_parser("number", help="Renumber a structure file.")
    n.add_argument("input", help="Input .pdb / .cif file")
    n.add_argument("-o", "--output", help="Output file (default: <input>_imgt.<ext>)")
    n.add_argument("--scheme", default="imgt", choices=["imgt", "aho"],
                   help="Numbering scheme (TCR supports imgt or aho). Default: imgt")
    n.add_argument("--species", default=None,
                   help="Comma-separated allowed species (e.g. human,mouse). Default: all")
    n.add_argument("--chains", default=None,
                   help="Comma-separated chain IDs to restrict to (default: all chains)")
    n.add_argument("--constant", default="renumber",
                   choices=["renumber", "original", "drop"],
                   help="Handling of non-variable-domain residues (e.g. constant domain). "
                        "Default: renumber sequentially after the V domain.")
    n.add_argument("--bit-score-threshold", type=float, default=80,
                   help="ANARCI bit-score threshold for accepting a domain. Default: 80")
    return p


def main(argv=None):
    args = build_parser().parse_args(argv)
    if args.cmd == "setup":
        run_setup()
    elif args.cmd == "number":
        run_number(args)


if __name__ == "__main__":
    main()
