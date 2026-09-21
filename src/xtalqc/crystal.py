"""Crystal mates and a lattice-aware check that a ligand site is clean.

    st = rcsb.load("7x11")
    mates = symmates(st, cutoff=5.0)            # like PyMOL symexp / ChimeraX crystalcontacts
    mates.make_mmcif_document().write_file("mates.cif")
    check("7x11__1__1.D__1.P", "1.P")
"""

from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache

import gemmi

from . import rcsb
from .rmsd import ligand_atoms
from .rnp import asym_id, parse_system_id

# waters are always removed; these are the usual crystallization additives
SOLVENTS = {
    "GOL", "EDO", "PEG", "PGE", "PG4", "PG0", "1PE", "2PE", "P6G", "P33", "MPD", "MRD", "PGO",
    "PGR", "BU3", "IPA", "EOH", "MOH", "DMS", "DMF", "ACN", "ACE", "ACT", "ACY", "FMT", "CIT",
    "FLC", "TAR", "TLA", "MLI", "MLA", "SIN", "MAL", "TRS", "MES", "EPE", "HEZ", "BME", "DTT",
    "IMD", "PO4", "PI", "2HP", "SO4", "SUL", "NO3", "SCN", "AZI", "CO3", "BCT", "NH4", "IOD",
    "BR", "CXS", "B3P", "BTB", "MPO", "CAQ", "OXL", "DIO", "ETF", "TMA", "SGM", "PE4", "PE8",
    "12P", "15P", "7PE", "XPE", "CPS", "LDA", "BOG", "HTG",
}


_load = lru_cache(maxsize=16)(rcsb.load)  # rows of one entry share a parse; check() copies it


def _is_ion(res: gemmi.Residue) -> bool:
    return not res.is_water() and len({a.name for a in res if not a.is_hydrogen()}) == 1


def _label(ch: gemmi.Chain, r: gemmi.Residue) -> str:
    return f"{r.name} {ch.name}/{str(r.seqid).strip()} [{r.subchain}]"


def clean(st: gemmi.Structure, keep: str | None = None) -> gemmi.Structure:
    """Copy without hydrogens, waters and ``SOLVENTS`` (subchain ``keep`` is kept)."""
    st = st.clone()
    st.remove_hydrogens()
    for ch in st[0]:
        for k in reversed(range(len(ch))):
            r = ch[k]
            if r.subchain != keep and r.entity_type != gemmi.EntityType.Polymer and \
                    (r.is_water() or r.name in SOLVENTS):
                del ch[k]
    st.remove_empty_chains()
    return st


def symmates(st: gemmi.Structure, cutoff: float = 5.0, around=None) -> gemmi.Structure:
    """The asymmetric unit plus every symmetry copy of a chain that comes within
    ``cutoff`` A of ``around`` (atoms; default the whole asymmetric unit).

    Copies are whole chains named ``<chain>_<op>_<a>_<b>_<c>``, e.g. ``B_3_0_0_-1``
    for operator 3 translated by (0, 0, -1) cells.
    """
    out = st.clone()
    for k in reversed(range(1, len(out))):
        del out[k]
    if not st.cell.is_crystal():
        return out
    model = st[0]
    ns = gemmi.NeighborSearch(model, st.cell, max(cutoff, 5.0)).populate()
    around = around if around is not None else [a for ch in model for r in ch for a in r]
    found = {}
    for a in around:
        for mk in ns.find_atoms(a.pos, "\0", radius=cutoff):
            ni = st.cell.find_nearest_pbc_image(a.pos, mk.to_cra(model).atom.pos, mk.image_idx)
            key = (mk.chain_idx, mk.image_idx, tuple(ni.pbc_shift))
            if key[1:] != (0, (0, 0, 0)):
                found[key] = True
    for ci, image, shift in sorted(found):
        tag = f"_{image}_" + "_".join(map(str, shift))
        ch = model[ci].clone()
        ch.name += tag
        ft = ns.get_image_transformation(image)
        for r in ch:
            r.subchain += tag
            for a in r:
                f = ft.apply(st.cell.fractionalize(a.pos))
                a.pos = st.cell.orthogonalize(gemmi.Fractional(f.x + shift[0], f.y + shift[1],
                                                               f.z + shift[2]))
        out[0].add_chain(ch)
    return out


@dataclass
class Report:
    pdb_id: str
    resname: str  # the ligand
    chain: str  # auth_asym_id
    resid: str  # auth_seq_id (+ insertion code)
    asym: str  # label_asym_id
    contact_chains: list[str] = field(default_factory=list)  # polymer chain copies within cutoff
    other_ligands: list[str] = field(default_factory=list)  # "<resname> <chain>/<resid> [<asym>]"
    ions: list[str] = field(default_factory=list)
    resolution: float | None = None
    rscc: float | None = None

    @property
    def sandwiched(self) -> bool:
        return len(self.contact_chains) > 1

    @property
    def clean(self) -> bool:
        return not (self.sandwiched or self.other_ligands or self.ions)


def check(system_id: str, ligand_instance: str, cutoff: float = 4.0,
          validation: bool = True) -> Report:
    """Crystal-environment QC of a Runs N' Poses ligand.

    Fails (``clean`` False) if more than one polymer chain copy -- in the
    asymmetric unit or a neighbouring lattice copy -- lies within ``cutoff`` A
    of the ligand, if the entry has other ligands, or if it has ions.
    Waters and ``SOLVENTS`` are ignored.  ``validation`` also asks RCSB for the
    resolution and the ligand RSCC.
    """
    pdb_id, asym = parse_system_id(system_id).pdb_id, asym_id(ligand_instance)
    st = clean(_load(pdb_id), keep=asym)
    lig = ligand_atoms(st, asym)
    lig_ch, lig_res = next((ch, r) for ch in st[0] for r in ch if r.subchain == asym)
    mates = symmates(st, cutoff + 1.0, around=lig)

    ns = gemmi.NeighborSearch(mates[0], gemmi.UnitCell(), cutoff + 1.0).populate()
    chains = set()
    for a in lig:
        for mk in ns.find_atoms(a.pos, "\0", radius=cutoff):
            cra = mk.to_cra(mates[0])
            if cra.residue.entity_type == gemmi.EntityType.Polymer:
                chains.add(cra.chain.name)

    others, ions, seen = [], [], {asym}
    for ch in st[0]:
        for r in ch:
            if r.entity_type == gemmi.EntityType.Polymer or r.subchain in seen:
                continue
            seen.add(r.subchain)  # one entry per molecule (first residue of a glycan)
            (ions if _is_ion(r) else others).append(_label(ch, r))

    rep = Report(pdb_id, lig_res.name, lig_ch.name, str(lig_res.seqid).strip(), asym,
                 sorted(chains), others, ions)
    if validation:
        rep.resolution, rep.rscc = rcsb.resolution(pdb_id), rcsb.rscc(pdb_id, asym)
    return rep
