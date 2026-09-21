# xtalqc

Crystal-structure QC for protein-ligand complexes, built on [gemmi](https://gemmi.readthedocs.io).

- Checks the crystal environment of a [Runs N' Poses](https://github.com/plinder-org/runs-n-poses) ligand: builds crystal mates, drops waters and crystallization additives, and flags ligands that are
  - sandwiched (more than one protein chain copy within 4 Å, including lattice neighbours),
  - accompanied by other ligands,
  - in an entry with ions.
- Resolution and ligand RSCC from RCSB
- Symmetry-corrected ligand RMSD and pocket dRMSD, after aligning chains by sequence and position (chain names are ignored)
- Crystal mates (like PyMOL `symexp` / ChimeraX `crystalcontacts`)
- MAE read/write (atoms and bonds)

## Installation

```bash
git clone https://github.com/siyoungkimlab/xtalqc.git
cd xtalqc
pip install -e .
```

## Command line

```bash
# one ligand: system_id + ligand instance; -o writes the cleaned ligand environment
xtalqc 7x11__1__1.D__1.P 1.P -o mates.cif

# a CSV with system_id and ligand_instance (or ligand_instance_chain) -> adds qc_* columns
xtalqc csv annotations.csv annotations_qc.csv -j 8
```

A ligand is `qc_clean` when it is not sandwiched, has no other ligands, and there are no ions.

## Python

```python
import gemmi
import xtalqc as xq

rep = xq.check("7x11__1__1.D__1.P", "1.P")   # .clean, .resname, .chain, .resid, .resolution, .rscc
st = xq.load("7x11")                          # RCSB mmCIF, cached in ~/.cache/xtalqc
mates = xq.symmates(st, cutoff=5.0)

pred = gemmi.read_structure("pred.cif")       # e.g. a Boltz / AF3 model
xq.ligand_rmsd(pred, st, "B", "P")            # ligand: label_asym_id, chain, or residue name
xq.drmsd(pred, st, "B", "P")                  # CA atoms within 6 Å of the ligand

(mae, bonds), = xq.read_mae("in.mae")
xq.write_mae("out.mae", mae, bonds)
```

## License

MIT
