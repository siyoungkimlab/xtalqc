"""Download entries and validation data from RCSB."""

from __future__ import annotations

import json
import os
import urllib.request
from pathlib import Path

import gemmi

CACHE = Path(os.environ.get("XTALQC_CACHE", Path.home() / ".cache" / "xtalqc"))


def _get(url: str) -> bytes:
    with urllib.request.urlopen(url, timeout=60) as r:
        return r.read()


def fetch(pdb_id: str) -> Path:
    """Path to the cached ``<pdb_id>.cif.gz``, downloading it if needed."""
    path = CACHE / f"{pdb_id.lower()}.cif.gz"
    if not path.exists():
        CACHE.mkdir(parents=True, exist_ok=True)
        path.write_bytes(_get(f"https://files.rcsb.org/download/{pdb_id.upper()}.cif.gz"))
    return path


def load(pdb_id: str) -> gemmi.Structure:
    """The deposited entry (first model, entities set up)."""
    st = gemmi.read_structure(str(fetch(pdb_id)))
    st.setup_entities()
    return st


def resolution(pdb_id: str) -> float | None:
    """Best reported resolution (A), None for methods without one."""
    info = json.loads(_get(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id.upper()}"))
    res = info.get("rcsb_entry_info", {}).get("resolution_combined")
    return min(res) if res else None


def rscc(pdb_id: str, asym: str) -> float | None:
    """Real-space correlation coefficient of a ligand instance (label_asym_id)."""
    url = f"https://data.rcsb.org/rest/v1/core/nonpolymer_entity_instance/{pdb_id.upper()}/{asym}"
    info = json.loads(_get(url))
    scores = info.get("rcsb_nonpolymer_instance_validation_score") or []
    return scores[0].get("RSCC") if scores else None
