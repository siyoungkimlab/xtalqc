"""Runs N' Poses / PLINDER identifiers.

A system id is ``<pdb>__<assembly>__<receptor chains>__<ligand chains>``, e.g.
``7x11__1__1.D__1.P``; chains are ``<instance>.<label_asym_id>`` joined by ``_``.
"""

from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class SystemID:
    pdb_id: str
    assembly: str
    receptors: tuple[str, ...]
    ligands: tuple[str, ...]


def parse_system_id(system_id: str) -> SystemID:
    pdb_id, assembly, receptors, ligands = system_id.split("__")
    return SystemID(pdb_id.lower(), assembly, tuple(receptors.split("_")),
                    tuple(ligands.split("_")))


def asym_id(chain: str) -> str:
    """``"1.P"`` -> ``"P"``: the label_asym_id in the deposited entry.

    The instance number only says which assembly copy it is; every copy has the
    same crystal environment, so the asymmetric-unit copy stands for all of them.
    """
    return chain.split(".", 1)[-1]
