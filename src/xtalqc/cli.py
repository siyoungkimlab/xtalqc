"""``xtalqc <system_id> <ligand_instance> [-o mates.cif]``
``xtalqc csv in.csv out.csv [-j 8]``"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from concurrent.futures import ThreadPoolExecutor
from dataclasses import asdict

from .crystal import Report, check, clean, symmates
from .rcsb import load
from .rmsd import ligand_atoms
from .rnp import asym_id, parse_system_id


def _row(rep) -> dict:
    d = asdict(rep) | {"sandwiched": rep.sandwiched, "clean": rep.clean, "error": ""}
    del d["pdb_id"]
    return {f"qc_{k}": ";".join(v) if isinstance(v, list) else v for k, v in d.items()}


def batch(argv) -> None:
    p = argparse.ArgumentParser(prog="xtalqc csv",
                                description="Add qc_* columns to a CSV of Runs N' Poses ligands")
    p.add_argument("input")
    p.add_argument("output")
    p.add_argument("--system-col", default="system_id")
    p.add_argument("--ligand-col", help="default: ligand_instance or ligand_instance_chain")
    p.add_argument("--cutoff", type=float, default=4.0)
    p.add_argument("--no-validation", action="store_true", help="skip RCSB resolution/RSCC")
    p.add_argument("-j", "--jobs", type=int, default=4, help="parallel downloads/checks")
    a = p.parse_args(argv)
    with open(a.input, newline="") as f:
        reader = csv.DictReader(f)
        rows, fields = list(reader), list(reader.fieldnames or [])
    lcol = a.ligand_col or next((c for c in ("ligand_instance", "ligand_instance_chain")
                                 if c in fields), "ligand_instance")
    for c in (a.system_col, lcol):
        if c not in fields:
            sys.exit(f"{a.input}: no column {c!r}")

    def run(row):
        try:
            return _row(check(row[a.system_col], row[lcol], a.cutoff,
                              validation=not a.no_validation))
        except Exception as e:  # noqa: BLE001 -- one bad entry should not stop the table
            return {"qc_error": f"{type(e).__name__}: {e}"}

    with ThreadPoolExecutor(a.jobs) as pool:
        results = list(pool.map(run, rows))
    qc = [k for k in _row(Report("", "", "", "", "")) if k not in fields]
    with open(a.output, "w", newline="") as f:
        w = csv.DictWriter(f, fields + qc)
        w.writeheader()
        for row, res in zip(rows, results):
            w.writerow(row | res)
    bad = sum(bool(r.get("qc_error")) for r in results)
    print(f"{len(rows)} rows, {sum(r.get('qc_clean') is True for r in results)} clean, "
          f"{bad} errors -> {a.output}", file=sys.stderr)


def main(argv=None) -> None:
    argv = sys.argv[1:] if argv is None else argv
    if argv[:1] == ["csv"]:
        return batch(argv[1:])
    p = argparse.ArgumentParser(description="Crystal-environment QC of a Runs N' Poses ligand "
                                "(or: xtalqc csv in.csv out.csv)")
    p.add_argument("system_id")
    p.add_argument("ligand_instance")
    p.add_argument("--cutoff", type=float, default=4.0)
    p.add_argument("--no-validation", action="store_true", help="skip RCSB resolution/RSCC")
    p.add_argument("-o", "--mates", help="write the cleaned ligand environment (mmCIF)")
    a = p.parse_args(argv)
    rep = check(a.system_id, a.ligand_instance, a.cutoff, validation=not a.no_validation)
    print(json.dumps(asdict(rep) | {"sandwiched": rep.sandwiched, "clean": rep.clean}, indent=2))
    if a.mates:
        asym = asym_id(a.ligand_instance)
        st = clean(load(parse_system_id(a.system_id).pdb_id), keep=asym)
        symmates(st, 12.0, around=ligand_atoms(st, asym)).make_mmcif_document().write_file(a.mates)
