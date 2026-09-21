"""Maestro (.mae/.maegz) structures: atoms and bonds only, no force field blocks."""

from __future__ import annotations

import gzip
import re

import gemmi

from .rmsd import infer_bonds

_TOKEN = re.compile(r'#[^#\n]*#?|"(?:[^"\\]|\\.)*"|[{}]|[^\s{}"]+')
_ARRAY = re.compile(r"(.+)\[(\d+)\]")
_ESC = re.compile(r"\\(.)")
# atomic number -> (maestro color, mmod type)
_MMOD = {1: (21, 48), 3: (4, 11), 6: (2, 14), 7: (43, 40), 8: (70, 23), 9: (8, 56),
         11: (4, 66), 12: (4, 72), 14: (14, 60), 15: (15, 53), 16: (13, 52),
         17: (13, 102), 19: (4, 67), 20: (4, 70)}


def _open(path, mode):
    path = str(path)
    return gzip.open(path, mode) if path.endswith("gz") else open(path, mode)


def _unquote(t: str) -> str:
    return _ESC.sub(r"\1", t[1:-1]) if t[:1] == '"' else t


def _value(t: str, kind: str):
    if t == "<>":
        return None
    t = _unquote(t)
    return float(t) if kind == "r" else int(float(t)) if kind in "ib" else t


def _keys(toks, k):
    keys = []
    while toks[k] != ":::":
        keys.append(toks[k])
        k += 1
    return keys, k + 1


def _children(toks, k, out):
    """Sub-blocks up to the closing '}'; returns the index after it."""
    while toks[k] != "}":
        name, k = toks[k], k + 1
        arr = _ARRAY.fullmatch(name)
        if arr:  # array block: name[n] { keys ::: rows ::: sub-blocks }
            name, n = arr.group(1), int(arr.group(2))
            cols, k = _keys(toks, k + 1)
            rows = []
            for _ in range(n):
                rows.append({c: _value(t, c[0]) for c, t in zip(cols, toks[k + 1:k + 1 + len(cols)])})
                k += 1 + len(cols)
            out[name] = rows
            k = _children(toks, k + 1, {})
        else:
            out[name], k = _block(toks, k)
    return k + 1


def _block(toks, k):
    """Parse ``{ keys ::: values sub-blocks }`` at toks[k] == '{'; returns (dict, next k)."""
    keys, k = _keys(toks, k + 1)
    out = {key: _value(t, key[0]) for key, t in zip(keys, toks[k:k + len(keys)])}
    return out, _children(toks, k + len(keys), out)


def read_mae(path) -> list[tuple[gemmi.Structure, list[tuple[int, int, int]]]]:
    """One ``(structure, bonds)`` per ct; bonds are 0-based ``(i, j, order)``
    over atoms in file order (the order ``structure[0].all()`` walks them)."""
    with _open(path, "rt") as f:
        toks = [t for t in _TOKEN.findall(f.read()) if t[0] != "#"]
    out, k = [], 0
    while k < len(toks):
        name = None
        if toks[k] != "{":
            name, k = toks[k], k + 1
        blk, k = _block(toks, k)
        if name in ("f_m_ct", "p_m_ct"):
            out.append(_to_structure(blk))
    return out


def _to_structure(ct):
    st = gemmi.Structure()
    st.name = ct.get("s_m_title") or ""
    model = gemmi.Model(1)
    rows = ct.get("m_atom", [])
    for row in rows:
        cname = (row.get("s_m_chain_name") or "").strip() or "A"
        if not len(model) or model[len(model) - 1].name != cname:
            model.add_chain(gemmi.Chain(cname))
        chain = model[len(model) - 1]
        num = row.get("i_m_residue_number") or 0
        icode = (row.get("s_m_insertion_code") or " ").strip() or " "
        rname = (row.get("s_m_pdb_residue_name") or "UNK").strip()
        if not len(chain) or (chain[len(chain) - 1].seqid.num, chain[len(chain) - 1].seqid.icode,
                              chain[len(chain) - 1].name) != (num, icode, rname):
            res = gemmi.Residue()
            res.name, res.seqid = rname, gemmi.SeqId(num, icode)
            chain.add_residue(res)
        a = gemmi.Atom()
        a.element = gemmi.Element(row.get("i_m_atomic_number") or 0)
        a.name = (row.get("s_m_pdb_atom_name") or "").strip() or a.element.name
        a.pos = gemmi.Position(row["r_m_x_coord"], row["r_m_y_coord"], row["r_m_z_coord"])
        a.charge = row.get("i_m_formal_charge") or 0
        a.occ = row.get("r_m_pdb_occupancy") or 1.0
        a.b_iso = row.get("r_m_pdb_tfactor") or 0.0
        chain[len(chain) - 1].add_atom(a)
    st.add_model(model)
    box = [ct.get(f"r_chorus_box_{v}{x}") for v in "abc" for x in "xyz"]
    if all(b is not None for b in box) and box[1] == box[2] == box[3] == box[5] == box[6] == box[7] == 0:
        st.cell = gemmi.UnitCell(box[0], box[4], box[8], 90, 90, 90)
    st.setup_entities()
    bonds = {}
    for b in ct.get("m_bond", []):
        i, j = sorted((b["i_m_from"] - 1, b["i_m_to"] - 1))
        bonds[i, j] = b.get("i_m_order") or 1
    return st, [(i, j, o) for (i, j), o in bonds.items()]


def _q(s: str) -> str:
    s = str(s)
    return '"' + s.replace("\\", "\\\\").replace('"', '\\"') + '"' if not s or re.search(r'[\s"]', s) else s


def write_mae(path, st: gemmi.Structure, bonds=None) -> None:
    """Write the first model; ``bonds`` as from :func:`read_mae`, else inferred (order 1)."""
    atoms = [(ch, r, a) for ch in st[0] for r in ch for a in r]
    if bonds is None:
        bonds = [(i, j, 1) for i, j in infer_bonds([a for _, _, a in atoms])]
    cols = ["r_m_x_coord", "r_m_y_coord", "r_m_z_coord", "i_m_residue_number",
            "s_m_insertion_code", "s_m_chain_name", "s_m_pdb_residue_name", "s_m_pdb_atom_name",
            "i_m_atomic_number", "i_m_formal_charge", "r_m_pdb_occupancy", "r_m_pdb_tfactor",
            "i_m_color", "i_m_mmod_type"]
    lines = ["{", "  s_m_m2io_version", "  :::", "  2.0.0", "}", "", "f_m_ct {", "  s_m_title",
             "  :::", f"  {_q(st.name)}", f"  m_atom[{len(atoms)}] {{"]
    lines += [f"    {c}" for c in cols] + ["    :::"]
    for k, (ch, r, a) in enumerate(atoms, 1):
        p = a.pos
        vals = [f"{p.x:.4f}", f"{p.y:.4f}", f"{p.z:.4f}", r.seqid.num, _q(r.seqid.icode),
                _q(ch.name), _q(f"{r.name:<4}"), _q(f" {a.name:<3}" if len(a.name) < 4 else a.name),
                a.element.atomic_number, a.charge, f"{a.occ:.2f}", f"{a.b_iso:.2f}",
                *_MMOD.get(a.element.atomic_number, (2, 64))]
        lines.append(f"    {k} " + " ".join(map(str, vals)))
    lines += ["    :::", "  }", f"  m_bond[{len(bonds)}] {{", "    i_m_from", "    i_m_to",
              "    i_m_order", "    :::"]
    lines += [f"    {k} {i + 1} {j + 1} {o}" for k, (i, j, o) in enumerate(bonds, 1)]
    lines += ["    :::", "  }", "}", ""]
    with _open(path, "wt") as f:
        f.write("\n".join(lines))
