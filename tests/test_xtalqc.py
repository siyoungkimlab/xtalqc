import os

import gemmi
import pytest

import xtalqc as xq

ACETATE = """\
HETATM    1  C1  ACT A   1       0.000   0.000   0.000  1.00 10.00           C
HETATM    2  C2  ACT A   1       1.520   0.000   0.000  1.00 10.00           C
HETATM    3  O1  ACT A   1       2.140   1.070   0.000  1.00 10.00           O
HETATM    4  O2  ACT A   1       2.140  -1.070   0.000  1.00 10.00           O
END
"""


def acetate():
    st = gemmi.read_pdb_string(ACETATE)
    st.setup_entities()
    return st


def test_parse_system_id():
    s = xq.parse_system_id("7x11__1__1.D__1.P")
    assert (s.pdb_id, s.receptors, s.ligands) == ("7x11", ("1.D",), ("1.P",))
    assert xq.asym_id("1.P") == "P"


def test_symmetry_rmsd():
    ref, mob = acetate(), acetate()
    o1, o2 = xq.ligand_atoms(mob, "ACT")[2:]
    o1.pos, o2.pos = gemmi.Position(o2.pos), gemmi.Position(o1.pos)
    assert xq.ligand_rmsd(mob, ref, "ACT", "ACT", align=False) == pytest.approx(0)
    for k, a in enumerate(xq.ligand_atoms(mob, "ACT")):  # no shared names: graph matching
        a.name = f"X{k}"
    assert xq.ligand_rmsd(mob, ref, "ACT", "ACT", align=False) == pytest.approx(0)


def test_mae_roundtrip(tmp_path):
    st = acetate()
    xq.write_mae(tmp_path / "a.mae", st)
    (st2, bonds), = xq.read_mae(tmp_path / "a.mae")
    assert sorted((i, j) for i, j, _ in bonds) == [(0, 1), (1, 2), (1, 3)]
    atoms = [a for ch in st2[0] for r in ch for a in r]
    assert [a.name for a in atoms] == ["C1", "C2", "O1", "O2"]
    assert atoms[3].pos.dist(gemmi.Position(2.14, -1.07, 0)) < 1e-3


@pytest.mark.skipif(not os.environ.get("XTALQC_NETWORK"), reason="needs RCSB")
def test_check_7x11():
    rep = xq.check("7x11__1__1.D__1.P", "1.P")
    assert (rep.resname, rep.chain, rep.resid, rep.asym) == ("86I", "D", "603", "P")
    assert rep.contact_chains == ["D"]
    assert rep.all_ions and not rep.ions  # the Mn ions are 13+ A away
    assert rep.all_other_ligands and not rep.other_ligands and rep.clean
    assert rep.resolution == pytest.approx(2.07) and rep.rscc == pytest.approx(0.855)


@pytest.mark.skipif(not os.environ.get("XTALQC_NETWORK"), reason="needs RCSB")
def test_align_ignores_chain_names():
    ref, mob = xq.load("7x11"), xq.load("7x11")
    for ch, name in zip(mob[0], "ZYXW"):
        ch.name = name
    assert xq.match_chains(mob, ref) == {"Z": "A", "Y": "B", "X": "C", "W": "D"}
    one = xq.load("7x11")  # one copy of a homotetramer, renamed: the ligand picks the partner
    for name in "ACD":
        one[0].remove_chain(name)
    one[0]["B"].name = "A"
    assert xq.ligand_rmsd(one, ref, "J", "J") == pytest.approx(0, abs=1e-6)
