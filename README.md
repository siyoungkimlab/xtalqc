# xtalqc

Crystal-structure QC for protein-ligand complexes, built on [gemmi](https://gemmi.readthedocs.io).

- Checks the crystal environment of a [Runs N' Poses](https://github.com/plinder-org/runs-n-poses) ligand: builds crystal mates, drops waters and crystallization additives, and flags a ligand that has, within 4 Å (lattice neighbours included),
  - more than one protein chain copy (sandwiched),
  - another ligand,
  - an ion.
- Resolution and ligand RSCC from RCSB
- Symmetry-corrected ligand RMSD and pocket dRMSD, after aligning chains by sequence and position (chain names are ignored)
- Crystal mates (like PyMOL `symexp` / ChimeraX `crystalcontacts`)
- MAE read/write (atoms and bonds)
- Cut a prepared MAE down to one system: receptor chains, the ligand and protein-bound ions

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

# prepared MAE (e.g. Protein Preparation Wizard, crystal frame) -> receptor + ligand + bound ions
xtalqc extract 7x11.prepped.mae 7x11__1__1.D__1.P 1.P 7x11_D_P.mae
```

`extract` matches residues to the PDB entry by position, so renamed chains and residues are
fine. An ion is kept when ≥3 receptor N/O/S atoms are within 3 Å. Bonds cut to dropped
residues (e.g. glycans, inter-chain disulfides) are reported as `cut_bonds`.

A ligand is `qc_clean` when nothing above is found. Ligands and ions anywhere in the entry are
listed in `qc_all_other_ligands` and `qc_all_ions` for reference.

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
