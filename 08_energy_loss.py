#!/usr/bin/env python3
"""How much of the Ar energy does the sheet absorb?  (the quantity of Fig. 4 in Bellido & Seminario, JPCC 2012)

    python 08_energy_loss.py            # reads impact_analysis.csv and bombardment_runs_v2/<tag>/traj.dump

The Ar velocity after the collision is taken from the last two frames of the impact stage (Ar is in free flight there:
the Ar-C force is zero beyond ~5 A). KE_after = 1/2 m v^2, E_deposited = E_0 - KE_after. Writes energy_loss.csv and
energy_loss.png.  Assumes the dump interval of 02_defect_engineering_stats.py (--dump-fs, default 1 fs).
"""
import argparse
import json
from collections import deque
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--csv", default="impact_analysis.csv")
p.add_argument("--workdir", default="bombardment_runs_v2")
p.add_argument("--dump-fs", type=float, default=1.0)
p.add_argument("--mass", type=float, default=39.948, help="projectile mass (u)")
p.add_argument("--out", default=".")
args = p.parse_args()
EV_PER_AMU_A2_FS2 = 103.6427          # 1 u * (A/fs)^2 in eV


def last_two_projectile_positions(path, n_atoms_hint=None):
    """Return the positions of the projectile (type 2) in the last two frames of a LAMMPS custom dump."""
    with open(path) as f:
        head = [next(f) for _ in range(9)]
        n = int(head[3])
    block = n + 9
    with open(path) as f:
        tail = deque(f, maxlen=2 * block)
    lines = list(tail)
    out = []
    for fr in (lines[:block], lines[block:]):
        box = [list(map(float, fr[5 + k].split()[:2])) for k in range(3)]
        for ln in fr[9:]:
            w = ln.split()
            if w[1] == "2":
                out.append([float(w[2]), float(w[3]), float(w[4])])
                break
    L = np.array([b[1] - b[0] for b in box])
    return np.array(out), L


rows = []
df = pd.read_csv(args.csv)
df = df[df.target == "random"]
for _, r in df.iterrows():
    tag = f"{r.projectile}_{r.target}_E{int(r.energy_eV)}_run{int(r.run)}"
    d = Path(args.workdir) / tag
    try:
        meta = json.loads((d / "meta.json").read_text())
        (x1, x2), L = last_two_projectile_positions(d / "traj.dump")
    except Exception as e:                                  # missing or unreadable run
        print("skip", tag, e)
        continue
    dx = x2 - x1
    dx[:2] -= L[:2] * np.round(dx[:2] / L[:2])              # periodic in x, y
    v = np.linalg.norm(dx) / args.dump_fs                   # A/fs
    ke = 0.5 * args.mass * v * v * EV_PER_AMU_A2_FS2
    rows.append(dict(tag=tag, E0=r.energy_eV, run=r.run, d_impact=r.d_impact_C, fate=r.fate, n_vacant=r.n_vacant,
                     ke_after=ke, ratio=ke / r.energy_eV, deposited=r.energy_eV - ke))
res = pd.DataFrame(rows)
res.to_csv(f"{args.out}/energy_loss.csv", index=False)
print(res.groupby("E0")[["ratio", "deposited"]].agg(["mean", "std"]).round(2).to_string())

BLUE, ORANGE, GREEN, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#9aa0a6"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(1, 2, figsize=(11, 4.4))
for fate, c in (("reflected", ORANGE), ("transmitted", BLUE)):
    g = res[res.fate == fate]
    ax[0].scatter(g.E0 * (1 + 0.012 * np.random.default_rng(0).normal(size=len(g))), g.ratio, c=c, s=22, alpha=.7, label=fate)
ax[0].set_xscale("log")
ax[0].set_xlabel("Ar kinetic energy (eV)")
ax[0].set_ylabel("KE after / KE before")
ax[0].set_title("Ar keeps most of its energy once it crosses", loc="left", fontsize=12)
ax[0].legend(frameon=False)
vac = res.n_vacant >= 1
ax[1].scatter(res[~vac].d_impact, res[~vac].deposited, c=GREY, s=22, label="no vacancy")
ax[1].scatter(res[vac].d_impact, res[vac].deposited, c=BLUE, s=26, label="vacancy")
ax[1].set_xlabel("distance to nearest C atom (Å)")
ax[1].set_ylabel("energy deposited in the sheet (eV)")
ax[1].set_title("Energy deposited decides the defect", loc="left", fontsize=12)
ax[1].legend(frameon=False)
fig.tight_layout()
fig.savefig(f"{args.out}/energy_loss.png", dpi=200)
print("wrote energy_loss.csv, energy_loss.png")