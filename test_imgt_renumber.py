"""
Tests for imgt_renumber.

These exercise the full sequence -> ANARCI -> remap path on small synthetic
structures (geometry is irrelevant to renumbering; only residue identity and
order matter). They are skipped automatically if the ANARCI backend / IMGT HMMs
aren't installed yet -- run `imgt-renumber setup` first to enable them.
"""
import io

import pytest

# Skip the whole module unless the backend is importable and the HMM DB is present.
anarci = pytest.importorskip("anarci")


def _backend_ready():
    try:
        from anarci import anarci as _run
        _run([("t", "DAGVIQSPRHEVTEMGQEVTLRCKPISGHNSLFWYRQTMMRGLELLIY"
                    "FNNNVPIDDSGMPEDRFSAKMPNASFSTLKIQPSEPRDSAVYFCASSY"
                    "VGNTGELFFGEGSRLTVL")],
             scheme="imgt", output=False)
        return True
    except Exception:
        return False


pytestmark = pytest.mark.skipif(
    not _backend_ready(),
    reason="ANARCI backend/HMMs not set up (run `imgt-renumber setup`)",
)

AA3 = {'A': 'ALA', 'R': 'ARG', 'N': 'ASN', 'D': 'ASP', 'C': 'CYS', 'Q': 'GLN',
       'E': 'GLU', 'G': 'GLY', 'H': 'HIS', 'I': 'ILE', 'L': 'LEU', 'K': 'LYS',
       'M': 'MET', 'F': 'PHE', 'P': 'PRO', 'S': 'SER', 'T': 'THR', 'W': 'TRP',
       'Y': 'TYR', 'V': 'VAL'}

VALPHA = ("KQEVTQIPAALSVPEGENLVLNCSFTDSAIYNLQWFRQDPGKGLTSLLLIQSSQREQTSGRL"
          "NASLDKSSGRSTLYIAASQPGDSATYLCAVRPLLDGTYIPTFGRGTSLIVHP")
VBETA = ("DAGVIQSPRHEVTEMGQEVTLRCKPISGHNSLFWYRQTMMRGLELLIYFNNNVPIDDSGMPE"
         "DRFSAKMPNASFSTLKIQPSEPRDSAVYFCASSYVGNTGELFFGEGSRLTVL")
# Beta with an artificially long CDR3 to force IMGT insertion codes.
VBETA_LONG = ("DAGVIQSPRHEVTEMGQEVTLRCKPISGHNSLFWYRQTMMRGLELLIYFNNNVPIDDSGMPE"
              "DRFSAKMPNASFSTLKIQPSEPRDSAVYFCASSYRRDDGGTTPPNNTGELFFGEGSRLTVL")


def _make_pdb(chains):
    """chains: list of (chain_id, one_letter_seq). Returns a Structure object."""
    from Bio.PDB import PDBParser

    lines, serial = [], 1
    for cid, seq in chains:
        x = 0.0
        for i, aa in enumerate(seq, start=1):
            res = AA3[aa]
            lines.append(
                "ATOM  " + f"{serial:>5}" + "  CA  " + f"{res:>3}" + " " + cid +
                f"{i:>4}" + "    " + f"{x:8.3f}{0.0:8.3f}{0.0:8.3f}"
                f"{1.0:6.2f}{0.0:6.2f}" + " " * 10 + " C")
            serial += 1
            x += 3.8
        lines.append("TER")
    lines.append("END")
    handle = io.StringIO("\n".join(lines) + "\n")
    return PDBParser(QUIET=True).get_structure("t", handle)


def _numbering(chain):
    from Bio.PDB.Polypeptide import is_aa
    return {res.id[1]: res.get_resname() for res in chain if is_aa(res, standard=False)}


def test_alpha_beta_detected_and_anchored():
    from imgt_renumber import renumber_structure

    s = _make_pdb([("D", VALPHA), ("E", VBETA)])
    report = renumber_structure(s, verbose=False)

    types = {r["chain"]: r["chain_type"] for r in report}
    assert types == {"D": "A", "E": "B"}

    model = next(s.get_models())
    for cid in ("D", "E"):
        num = _numbering(model[cid])
        # IMGT conserved anchors
        assert num[23] == "CYS"
        assert num[41] == "TRP"
        assert num[104] == "CYS"
        assert num[118] == "PHE"


def test_long_cdr3_gets_symmetric_insertion_codes():
    from Bio.PDB.Polypeptide import is_aa
    from imgt_renumber import renumber_structure

    s = _make_pdb([("E", VBETA_LONG)])
    renumber_structure(s, verbose=False)

    model = next(s.get_models())
    icodes = [(r.id[1], r.id[2].strip())
              for r in model["E"] if is_aa(r, standard=False) and r.id[2].strip()]
    # IMGT places CDR3 insertions symmetrically about 111 and 112.
    assert any(n == 111 for n, ic in icodes)
    assert any(n == 112 for n, ic in icodes)
    assert all(ic.isalpha() for _, ic in icodes)


def test_constant_drop_removes_tail():
    from Bio.PDB.Polypeptide import is_aa
    from imgt_renumber import renumber_structure

    tail = "PNIQNPDPAVYQLRDSKSSDKSVCLFTDFDSQTNVSQSKDSDVYITDKTVLDMRSMDFK"
    s = _make_pdb([("D", VALPHA + tail)])
    renumber_structure(s, constant="drop", verbose=False)

    model = next(s.get_models())
    nums = [r.id[1] for r in model["D"] if is_aa(r, standard=False)]
    # V domain only; nothing numbered beyond the IMGT range.
    assert max(nums) <= 128
    assert len(nums) == len(VALPHA)
