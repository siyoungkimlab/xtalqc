"""Protein superposition, symmetry-corrected ligand RMSD and pocket dRMSD.

    pairs = residue_pairs(model, crystal[0], chains={"A": "D"})
    ligand_rmsd(model, crystal, "B", "P")   # chains paired by match_chains()
    drmsd(model, crystal, "B", "P", chains={"A": "D"})

Ligand symmetry: every mapping of mobile onto reference atoms that keeps
elements and bonds is tried, and the smallest value wins.  Bonds are
inferred from covalent radii (mmCIF entries carry none for ligands).
"""

from __future__ import annotations

import itertools
import re

import gemmi
import networkx as nx
import numpy as np

MAX_MAPPINGS = 10000


def _model(st):
    return st[0] if isinstance(st, gemmi.Structure) else st


def ligand_atoms(st, key: str) -> list[gemmi.Atom]:
    """Heavy atoms of a ligand, one conformer.

    ``key`` is a label_asym_id (subchain), else an auth chain name (non-polymer
    residues only), else a residue name.
    """
    model = _model(st)
    for match in (lambda ch, r: r.subchain == key,
                  lambda ch, r: ch.name == key and r.entity_type != gemmi.EntityType.Polymer,
                  lambda ch, r: r.name == key):
        res = [r for ch in model for r in ch if match(ch, r)]
        if res:
            break
    else:
        raise ValueError(f"no ligand {key!r}")
    if res[0].subchain != key:  # a name may occur many times: keep the first instance
        res = [r for r in res if r.subchain == res[0].subchain]
    atoms, alt = [], None
    for r in res:
        for a in r:
            if a.is_hydrogen():
                continue
            if a.altloc != "\0":
                alt = alt or a.altloc
                if a.altloc != alt:
                    continue
            atoms.append(a)
    return atoms


def infer_bonds(atoms, tol: float = 0.45) -> list[tuple[int, int]]:
    """Atom index pairs closer than the sum of covalent radii plus ``tol``.

    Metals and different alternate conformers are never bonded.
    """
    xyz = np.array([a.pos.tolist() for a in atoms]).reshape(-1, 3)
    rad = np.array([a.element.covalent_r for a in atoms])
    alt = np.array([a.altloc for a in atoms])
    metal = np.array([a.element.is_metal for a in atoms], bool)
    cells: dict = {}
    for k, c in enumerate(map(tuple, np.floor(xyz / 3.0).astype(int))):
        cells.setdefault(c, []).append(k)
    out = []
    for (x, y, z), mine in cells.items():
        near = [k for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)
                for k in cells.get((x + dx, y + dy, z + dz), ())]
        i, j = np.array(mine)[:, None], np.array(near)[None]
        ok = (i < j) & (np.linalg.norm(xyz[i] - xyz[j], axis=-1) < rad[i] + rad[j] + tol)
        ok &= (alt[i] == alt[j]) | (alt[i] == "\0") | (alt[j] == "\0")
        ok &= ~metal[i] & ~metal[j]
        a, b = np.nonzero(ok)
        out += zip(i[a, 0].tolist(), j[0, b].tolist())
    return sorted(out)


def _graph(atoms) -> nx.Graph:
    g = nx.Graph()
    g.add_nodes_from((k, {"el": a.element.name}) for k, a in enumerate(atoms))
    g.add_edges_from(infer_bonds(atoms))
    return g


def mappings(mobile, reference, limit: int = MAX_MAPPINGS) -> np.ndarray:
    """(n, natoms): each row lists the mobile atom that plays each reference atom.

    With matching atom names the name pairing is combined with the reference
    graph's automorphisms, so a distorted mobile pose cannot spoil its bonds;
    otherwise the mobile graph is inferred as well and matched to the reference.
    """
    if len(mobile) != len(reference):
        raise ValueError(f"ligands differ: {len(mobile)} vs {len(reference)} heavy atoms")
    same = lambda a, b: a["el"] == b["el"]
    gref = _graph(reference)
    names = {a.name: k for k, a in enumerate(mobile)}
    if len(names) == len(mobile) and all(a.name in names for a in reference):
        base = np.array([names[a.name] for a in reference])
        gm = nx.algorithms.isomorphism.GraphMatcher(gref, gref, node_match=same)
        maps = [base[[m[i] for i in range(len(reference))]]
                for m in itertools.islice(gm.isomorphisms_iter(), limit)]
    else:
        gm = nx.algorithms.isomorphism.GraphMatcher(gref, _graph(mobile), node_match=same)
        maps = [[m[i] for i in range(len(reference))]
                for m in itertools.islice(gm.isomorphisms_iter(), limit)]
    if not maps:
        raise ValueError("ligand bond graphs differ")
    return np.array(maps)


def _align_chains(mp, rp) -> tuple[list, float]:
    """Aligned residue pairs of two polymers and their identity (over the shorter)."""
    aln = gemmi.align_string_sequences([r.name for r in mp], [r.name for r in rp], [])
    pairs, i, j = [], 0, 0
    for n, op in re.findall(r"(\d+)([MID])", aln.cigar_str()):
        n = int(n)
        if op == "M":
            pairs += [(mp[i + k], rp[j + k]) for k in range(n) if mp[i + k].name == rp[j + k].name]
        i += n if op in "MI" else 0
        j += n if op in "MD" else 0
    return pairs, len(pairs) / max(1, min(len(mp), len(rp)))


def _fit(pairs, atom: str = "CA") -> gemmi.SupResult:
    """Superposition of the mobile onto the reference residues of ``pairs``."""
    fixed, moving = [], []
    for m, r in pairs:
        am, ar = m.find_atom(atom, "*"), r.find_atom(atom, "*")
        if am and ar:
            moving.append(am.pos)
            fixed.append(ar.pos)
    if len(fixed) < 3:
        raise ValueError(f"only {len(fixed)} aligned {atom} atoms")
    return gemmi.superpose_positions(fixed, moving)


def _nearest_chain(polys, atoms) -> int:
    lig = _xyz(atoms)
    return int(np.argmin([np.linalg.norm(_xyz([a for r in p for a in r])[:, None] - lig[None],
                                         axis=-1).min() for p in polys]))


def match_chains(mobile, reference, min_identity: float = 0.3,
                 anchor=None) -> dict[str, str]:
    """Pair polymer chains by sequence and position, ignoring chain names.

    As in ChimeraX matchmaker, the pair with the most aligned residues fixes
    the frame; each other mobile chain then takes the nearest (by CA centroid,
    in that frame) unused reference chain it aligns to with ``min_identity``,
    so copies in a homo-oligomer are told apart by where they are.
    ``anchor``: (mobile ligand atoms, reference ligand atoms); the chains
    nearest them fix the frame instead, when they align.
    """
    mnames = [c.name for c in _model(mobile) if c.get_polymer()]
    rnames = [c.name for c in _model(reference) if c.get_polymer()]
    mob = [_model(mobile)[n].get_polymer() for n in mnames]
    ref = [_model(reference)[n].get_polymer() for n in rnames]
    aln = {(i, j): _align_chains(mp, rp) for i, mp in enumerate(mob) for j, rp in enumerate(ref)}
    ok = {k: v for k, v in aln.items() if v[1] >= min_identity}
    if not ok:
        raise ValueError("no mobile chain aligns to a reference chain")
    bi, bj = max(ok, key=lambda k: len(ok[k][0]))
    if anchor is not None:
        pair = _nearest_chain(mob, anchor[0]), _nearest_chain(ref, anchor[1])
        bi, bj = pair if pair in ok else (bi, bj)
    tr = _fit(ok[bi, bj][0]).transform

    def centroid(poly, move=False):
        xyz = [r.get_ca().pos for r in poly if r.get_ca()] or [r[0].pos for r in poly]
        c = gemmi.Position(*np.mean([p.tolist() for p in xyz], axis=0))
        return gemmi.Position(tr.apply(c)) if move else c

    out, used = {bi: bj}, {bj}
    best = lambda i: max(len(v[0]) for k, v in ok.items() if k[0] == i)
    for i in sorted({i for i, _ in ok} - {bi}, key=best, reverse=True):
        cands = [j for (a, j) in ok if a == i and j not in used]
        if cands:
            c = centroid(mob[i], move=True)
            j = min(cands, key=lambda j: c.dist(centroid(ref[j])))
            out[i] = j
            used.add(j)
    return {mnames[i]: rnames[j] for i, j in sorted(out.items())}


def residue_pairs(mobile, reference, chains: dict[str, str] | None = None, anchor=None):
    """Aligned (mobile, reference) residue pairs of polymer chains.

    ``chains`` maps mobile to reference chain names; by default :func:`match_chains`.
    """
    mob, ref = _model(mobile), _model(reference)
    chains = chains or match_chains(mob, ref, anchor=anchor)
    return [p for m, r in chains.items()
            for p in _align_chains(mob[m].get_polymer(), ref[r].get_polymer())[0]]


def superpose(mobile, reference, chains=None, atom: str = "CA", anchor=None) -> gemmi.SupResult:
    """Superpose ``mobile`` onto ``reference`` by aligned ``atom`` positions, in place."""
    sup = _fit(residue_pairs(mobile, reference, chains, anchor), atom)
    _model(mobile).transform_pos_and_adp(sup.transform)
    return sup


def _xyz(atoms) -> np.ndarray:
    return np.array([a.pos.tolist() for a in atoms])


def ligand_rmsd(mobile, reference, mobile_ligand: str, reference_ligand: str,
                chains=None, align: bool = True) -> float:
    """Symmetry-corrected ligand RMSD (A) after superposing the proteins.

    ``mobile`` is moved in place when ``align`` is True.
    """
    mob, ref = ligand_atoms(mobile, mobile_ligand), ligand_atoms(reference, reference_ligand)
    if align:
        superpose(mobile, reference, chains, anchor=(mob, ref))
    maps = mappings(mob, ref)
    d2 = ((_xyz(mob)[maps] - _xyz(ref)[None]) ** 2).sum(-1).mean(-1)
    return float(np.sqrt(d2.min()))


def drmsd(mobile, reference, mobile_ligand: str, reference_ligand: str,
          chains=None, cutoff: float = 6.0, atom: str = "CA") -> float:
    """Pocket-ligand distance RMSD (A); no superposition needed.

    The pocket is the reference ``atom`` atoms (C-alpha) within ``cutoff`` of
    any reference ligand heavy atom, paired with the aligned mobile residue's.
    """
    mob, ref = ligand_atoms(mobile, mobile_ligand), ligand_atoms(reference, reference_ligand)
    lref = _xyz(ref)
    pm, pr = [], []
    for m, r in residue_pairs(mobile, reference, chains, anchor=(mob, ref)):
        am, ar = m.find_atom(atom, "*"), r.find_atom(atom, "*")
        if am and ar and np.linalg.norm(lref - ar.pos.tolist(), axis=1).min() <= cutoff:
            pm.append(am)
            pr.append(ar)
    if not pr:
        raise ValueError(f"no pocket atoms within {cutoff} A of the reference ligand")
    dist = lambda p, lig: np.linalg.norm(p[:, None] - lig[None], axis=-1)
    dref = dist(_xyz(pr), lref)
    dmob = dist(_xyz(pm), _xyz(mob))
    maps = mappings(mob, ref)
    d2 = ((dmob[:, maps] - dref[:, None]) ** 2).mean(axis=(0, 2))
    return float(np.sqrt(d2.min()))
