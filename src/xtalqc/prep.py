"""Cut a prepared (e.g. Protein Preparation Wizard) MAE down to one Runs N' Poses system.

    extract("7x11.prepped.mae", "7x11__1__1.D__1.P", "1.P", "7x11_D_P.mae")

Keeps the system's receptor chains, the ligand, and protein-bound ions; drops
waters, other ligands and other chains.  Residues are matched to the deposited
entry by position, so renamed chains and residues are fine; caps and other
prep-built residues come along, and every atom property is preserved.
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass

import gemmi
import networkx as nx
import numpy as np

from . import rcsb
from .mae import read_cts, write_cts
from .rnp import asym_id, parse_system_id

DONORS = {7, 8, 16}  # N, O, S


@dataclass
class Extracted:
    receptors: list[str]  # MAE chain names kept
    ligand: str  # "<resname> <chain>/<resid>"
    ions: list[str]  # protein-bound ions kept
    cut_bonds: list[str]  # bonds from kept to dropped residues (open valences!)
    natoms: int


def _key(row) -> tuple:
    return ((row.get("s_m_chain_name") or "").strip(), row.get("i_m_residue_number") or 0,
            (row.get("s_m_insertion_code") or "").strip(), row.get("s_m_pdb_residue_name") or "")


def _label(row) -> str:
    chain, num, icode, name = _key(row)
    return f"{name.strip()} {chain}/{num}{icode}"


def extract(mae_in, system_id: str, ligand_instance: str, mae_out=None,
            ion_cutoff: float = 3.0, min_coordination: int = 3,
            tolerance: float = 1.5) -> Extracted:
    """Write the receptor chains, ligand and protein-bound ions of a system to ``mae_out``.

    Residues of the prepared file are matched to the deposited entry by
    position, not by name: each heavy atom takes the nearest same-element
    crystal atom within ``tolerance`` A, and a residue belongs to the chain
    (label_asym_id) most of its heavy atoms match.  Residues with no match
    (caps, rebuilt loops) are kept when bonded to kept ones.  Chain names,
    numbering and residue names may differ, but the prepared file must stay in
    the crystal frame.  An ion is protein-bound when at least
    ``min_coordination`` receptor N/O/S atoms lie within ``ion_cutoff`` A.
    """
    sid = parse_system_id(system_id)
    if len({c.split(".")[0] for c in sid.receptors + (ligand_instance,)}) > 1:
        raise ValueError(f"{system_id}: chains from different assembly copies cannot be cut "
                         "from an asymmetric-unit file")
    want, lig_asym = {asym_id(c) for c in sid.receptors}, asym_id(ligand_instance)
    crystal = rcsb.load(sid.pdb_id)
    crystal.remove_hydrogens()
    model = crystal[0]
    ns = gemmi.NeighborSearch(model, gemmi.UnitCell(), max(5.0, tolerance)).populate()

    ct = read_cts(mae_in)[0]
    atoms, bonds = ct["m_atom"], ct.get("m_bond", [])
    xyz = np.array([[a["r_m_x_coord"], a["r_m_y_coord"], a["r_m_z_coord"]] for a in atoms])
    anum = np.array([a.get("i_m_atomic_number") or 0 for a in atoms])

    def asym_of(k):
        p, best, where = gemmi.Position(*xyz[k]), tolerance, None
        for mk in ns.find_atoms(p, "\0", radius=tolerance):
            cra = mk.to_cra(model)
            if cra.atom.element.atomic_number == anum[k] and cra.atom.pos.dist(p) <= best:
                best, where = cra.atom.pos.dist(p), cra.residue.subchain
        return where

    # residues: runs of rows with the same chain/number/name (names alone may repeat)
    units, unit_of = [], np.empty(len(atoms), int)
    for k, row in enumerate(atoms):
        if not units or _key(atoms[units[-1][0]]) != _key(row):
            units.append([])
        units[-1].append(k)
        unit_of[k] = len(units) - 1
    owner = []
    for u in units:
        votes = Counter(asym_of(k) for k in u if anum[k] > 1)
        top, n = votes.most_common(1)[0] if votes else (None, 0)
        owner.append(top if n * 2 >= sum(votes.values()) else None)

    g = nx.Graph()
    g.add_nodes_from(range(len(units)))
    g.add_edges_from((unit_of[b["i_m_from"] - 1], unit_of[b["i_m_to"] - 1]) for b in bonds)
    keep = {u for u, o in enumerate(owner) if o in want or o == lig_asym}
    todo = list(keep)
    while todo:  # unmatched residues bonded to kept ones: caps, rebuilt parts
        for v in g[todo.pop()]:
            if v not in keep and owner[v] is None:
                keep.add(v)
                todo.append(v)
    rec_units = [u for u in keep if owner[u] in want]
    lig_units = [u for u in keep if owner[u] == lig_asym]
    if not rec_units or not lig_units:
        raise ValueError(f"{mae_in}: receptor or ligand of {system_id} not found within "
                         f"{tolerance} A of the deposited coordinates")

    rec = np.array([k for u in rec_units for k in units[u] if anum[k] in DONORS])
    ions = []
    for u, idx in enumerate(units):
        heavy = [k for k in idx if anum[k] > 1]
        if u in keep or len(heavy) != 1 or anum[heavy[0]] in (6, 7, 8) or g.degree(u):
            continue  # ions: one unbonded heavy atom that is not C/N/O (so not water)
        if (np.linalg.norm(xyz[rec] - xyz[heavy[0]], axis=1) <= ion_cutoff).sum() >= min_coordination:
            keep.add(u)
            ions.append(_label(atoms[heavy[0]]))

    order = [k for u in sorted(keep) for k in units[u]]
    new = {old: i + 1 for i, old in enumerate(order)}
    out = {k: v for k, v in ct.items() if not isinstance(v, (list, dict))}
    out["m_atom"] = [atoms[k] for k in order]
    out["m_bond"], cut = [], []
    for b in bonds:
        i, j = b["i_m_from"] - 1, b["i_m_to"] - 1
        if i in new and j in new:
            out["m_bond"].append(b | {"i_m_from": new[i], "i_m_to": new[j]})
        elif (i in new) != (j in new) and i < j:
            cut.append(f"{_label(atoms[i])}:{atoms[i]['s_m_pdb_atom_name'].strip()} - "
                       f"{_label(atoms[j])}:{atoms[j]['s_m_pdb_atom_name'].strip()}")
    if mae_out is not None:
        write_cts(mae_out, [out])
    receptors = sorted({_key(atoms[units[u][0]])[0] for u in rec_units})
    return Extracted(receptors, _label(atoms[units[lig_units[0]][0]]), ions, cut, len(order))
