#!/usr/bin/env python3
"""Cut a small region around the impact point out of a finished run, for viewing in OVITO / VMD / ase gui.

    python 06_export_defect.py Ar_random_E90_run11 Ar_random_E100_run9 --radius 8

For each run it writes  defect_views/<tag>_final.extxyz  (the relaxed structure after the impact) and
defect_views/<tag>_initial.extxyz (the same region before the impact), keeping all atoms within --radius (A) of the
impact point, in the periodic images closest to it. Each atom carries extra columns:
    dz     height above the sheet plane (A)
    coord  number of C neighbours (cut-off 1.85 A)
    index  original atom index (matches the numbers printed by 05_inspect_run.py)
In OVITO colour the atoms by 'dz' or 'coord'; in ase gui use  View > Colors.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from ase.geometry import find_mic
from ase.io import read, write
from ase.neighborlist import neighbor_list

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("tags", nargs="+")
p.add_argument("--indir", default="results_v2")
p.add_argument("--workdir", default="bombardment_runs_v2")
p.add_argument("--radius", type=float, default=8.0, help="keep atoms within this distance (A) of the impact point")
p.add_argument("--cc-cut", type=float, default=1.85)
p.add_argument("--out", default="defect_views")
args = p.parse_args()
Path(args.out).mkdir(exist_ok=True)


def crop(atoms, centre_xy, ref_z, with_props=True):
    atoms = atoms.copy()
    n = len(atoms)
    c_idx = np.array([k for k, s in enumerate(atoms.get_chemical_symbols()) if s == "C"])
    # coordination on the full periodic structure
    carbon = atoms[c_idx]
    i, _ = neighbor_list("ij", carbon, args.cc_cut)
    coord = np.zeros(n, dtype=int)
    coord[c_idx] = np.bincount(i, minlength=len(c_idx))
    # shift each atom to the periodic image closest to the impact point (x, y only)
    vec = atoms.positions.copy()
    vec[:, :2] -= centre_xy
    vec[:, 2] = 0.0
    vec, _ = find_mic(vec, atoms.cell, pbc=[True, True, False])
    r = np.linalg.norm(vec[:, :2], axis=1)
    keep = np.where(r <= args.radius)[0]
    new = atoms[keep]
    new.positions[:, :2] = vec[keep, :2]            # centred on the impact point (0, 0)
    new.set_pbc(False)
    new.set_cell(None)
    new.arrays["dz"] = atoms.positions[keep, 2] - ref_z
    new.arrays["coord"] = coord[keep]
    new.arrays["index"] = keep
    return new


for tag in args.tags:
    rundir = Path(args.workdir) / tag
    final = read(Path(args.indir) / f"graphene_defective_{tag}.xyz")
    pristine = read(rundir / "initial.xyz")
    xy0 = np.array(json.loads((rundir / "meta.json").read_text())["impact_xy"])
    is_c = np.array(final.get_chemical_symbols()) == "C"
    z_final = float(np.median(final.positions[is_c, 2]))
    z_init = float(np.median(pristine.positions[np.array(pristine.get_chemical_symbols()) == "C", 2]))
    write(Path(args.out) / f"{tag}_final.extxyz", crop(final, xy0, z_final))
    write(Path(args.out) / f"{tag}_initial.extxyz", crop(pristine, xy0, z_init))
    print(f"{tag}: wrote {args.out}/{tag}_final.extxyz and {tag}_initial.extxyz")