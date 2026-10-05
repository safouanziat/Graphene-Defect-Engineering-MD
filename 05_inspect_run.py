#!/usr/bin/env python3
"""Look inside one run: which C atoms ended up far from the sheet plane, and why.

    python 05_inspect_run.py Ar_random_E90_run11
    python 05_inspect_run.py Ar_random_E60_run6 Ar_random_E70_run23 --site-z 1.5

For every atom further than --site-z from the sheet plane it prints the height, the in-plane distance to the
impact point, the number of C neighbours (cutoff --cc-cut), the displacement from its original lattice site and
its neighbours. It also prints how rippled the rest of the sheet is, to tell a few displaced atoms from a bulge.
"""
import argparse
import json
from pathlib import Path

import numpy as np
from ase.geometry import find_mic
from ase.io import read
from ase.neighborlist import neighbor_list

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("tags", nargs="+", help="run tags such as Ar_random_E90_run11")
p.add_argument("--indir", default="results_v2")
p.add_argument("--workdir", default="bombardment_runs_v2")
p.add_argument("--site-z", type=float, default=1.5)
p.add_argument("--eject", type=float, default=3.0)
p.add_argument("--cc-cut", type=float, default=1.85)
args = p.parse_args()

for tag in args.tags:
    final = read(Path(args.indir) / f"graphene_defective_{tag}.xyz")
    rundir = Path(args.workdir) / tag
    pristine = read(rundir / "initial.xyz")
    xy0 = np.array(json.loads((rundir / "meta.json").read_text())["impact_xy"])

    sym = np.array(final.get_chemical_symbols())
    is_c = sym == "C"
    c_idx = np.where(is_c)[0]
    pos = final.positions[is_c]
    z_sheet = float(np.median(pos[:, 2]))
    dz = pos[:, 2] - z_sheet

    disp, _ = find_mic(pos - pristine.positions[: len(pos)], pristine.cell, pbc=[True, True, False])
    disp = disp - np.median(disp, axis=0)

    carbon = final[is_c]
    i, j = neighbor_list("ij", carbon, args.cc_cut)
    coord = np.bincount(i, minlength=len(carbon))
    nbrs = {k: sorted(j[i == k]) for k in range(len(carbon))}

    # in-plane distance from the impact point (use the ORIGINAL site of each atom, so a displaced atom is measured
    # from where it started)
    d0 = pristine.positions[: len(pos), :2] - xy0
    d0 = np.column_stack([d0, np.zeros(len(d0))])
    d0, _ = find_mic(d0, pristine.cell, pbc=[True, True, False])
    r_imp = np.linalg.norm(d0[:, :2], axis=1)

    far = np.where(np.abs(dz) > args.site_z)[0]
    ripple = np.abs(dz[np.abs(dz) <= args.site_z])
    print(f"\n=== {tag} ===")
    print(f"sheet plane z = {z_sheet:.2f} A; rest of the sheet: median |dz| = {np.median(ripple):.2f} A, "
          f"99th percentile = {np.percentile(ripple, 99):.2f} A, max = {ripple.max():.2f} A")
    print(f"atoms with |dz| > {args.site_z} A: {len(far)}  "
          f"({int((np.abs(dz) > args.eject).sum())} of them beyond {args.eject} A = ejected)")
    print("  idx |   dz(A) | r_from_impact(A) | coord | |disp|(A) | neighbours (idx:dz)")
    for k in far[np.argsort(-np.abs(dz[far]))]:
        nb = ", ".join(f"{n}:{dz[n]:+.2f}" for n in nbrs[k])
        print(f"  {k:4d} | {dz[k]:+7.2f} | {r_imp[k]:16.2f} | {coord[k]:5d} | {np.linalg.norm(disp[k]):9.2f} | {nb}")
    # the under-coordinated atoms that stay in the sheet (the neighbours of a vacancy)
    under = [k for k in range(len(pos)) if coord[k] < 3 and abs(dz[k]) <= args.eject]
    if under:
        print("under-coordinated atoms still in the sheet (idx:dz, coord): "
              + ", ".join(f"{k}:{dz[k]:+.2f},c{coord[k]}" for k in under))