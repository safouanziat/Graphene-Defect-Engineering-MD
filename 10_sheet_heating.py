#!/usr/bin/env python3
"""How much does the sheet heat up during the impact stage?  (finite-size check, no new MD needed)

    python 10_sheet_heating.py

The impact stage is microcanonical, so the energy the Ar leaves in the sheet raises its temperature. In a small periodic
cell this can be large. The script reads the thermo output of each run's log.lammps (sheet temperature every 1000 steps),
and reports T at the start and at the end of the impact stage, dT = T_end - T_start, and the mean dT per energy.
If 08_energy_loss.py has been run (energy_loss.csv) it also compares dT with the equipartition estimate
dT ~ E_dep / (3 N kB). This estimate holds once the deposited energy is shared equally between kinetic and
potential energy (harmonic lattice). It is read at the end of the impact stage, before that has happened, so the
measured kinetic temperature can lie above it (up to about 2x if all the energy were still kinetic): it is an
estimate of the scale, not an upper bound.
Writes heating.csv and heating.png.
"""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--csv", default="impact_analysis.csv")
p.add_argument("--workdir", default="bombardment_runs_v2")
p.add_argument("--energy-loss", default="energy_loss.csv")
p.add_argument("--n-sheet", type=int, default=2048, help="number of C atoms in the sheet")
p.add_argument("--out", default=".")
args = p.parse_args()
KB = 8.617333e-5


def thermo_blocks(path):
    """Return a list of blocks, each an array with columns (step, time, c_tsheet, ...)."""
    blocks, cur, inside = [], [], False
    for line in Path(path).read_text().splitlines():
        w = line.split()
        if w[:1] == ["Step"] and "c_tsheet" in w:
            inside, cur = True, []
            continue
        if inside:
            if line.startswith("Loop time"):
                if cur:
                    blocks.append(np.array(cur))
                inside = False
                continue
            try:
                cur.append([float(x) for x in w])
            except ValueError:
                inside = False
                if cur:
                    blocks.append(np.array(cur))
    return blocks


rows = []
df = pd.read_csv(args.csv)
df = df[df.target == "random"]
for _, r in df.iterrows():
    tag = f"{r.projectile}_{r.target}_E{int(r.energy_eV)}_run{int(r.run)}"
    d = Path(args.workdir) / tag
    try:
        meta = json.loads((d / "meta.json").read_text())
        blocks = thermo_blocks(d / "log.lammps")
        n_imp = int(meta["n_impact"])
        imp = [b for b in blocks if abs((b[-1, 0] - b[0, 0]) - n_imp) <= 1 and len(b) > 2]
        b = imp[0]
    except Exception as e:
        print("skip", tag, e)
        continue
    rows.append(dict(tag=tag, E0=r.energy_eV, run=r.run, T_start=b[0, 2], T_end=b[-1, 2], dT=b[-1, 2] - b[0, 2]))
res = pd.DataFrame(rows)
if Path(args.energy_loss).exists() and len(res):
    el = pd.read_csv(args.energy_loss)[["tag", "deposited"]]
    res = res.merge(el, on="tag", how="left")
    res["dT_equipartition"] = res.deposited / (3 * args.n_sheet * KB)
res.to_csv(f"{args.out}/heating.csv", index=False)
print(res.groupby("E0")[["T_start", "T_end", "dT"] + (["dT_equipartition"] if "dT_equipartition" in res else [])]
      .mean().round(1).to_string())

plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})
fig, ax = plt.subplots(figsize=(6.4, 4.4))
ax.scatter(res.E0, res.dT, s=22, c="#2a78d6", alpha=0.7, label="measured")
if "dT_equipartition" in res:
    ax.scatter(res.E0, res.dT_equipartition, s=14, c="#9aa0a6", marker="x", label="equipartition estimate $E_{dep}/3Nk_B$")
ax.set_xlabel("Ar kinetic energy (eV)")
ax.set_ylabel("sheet temperature rise, end of impact stage (K)")
ax.set_title(f"Heating of a {args.n_sheet}-atom periodic sheet", loc="left", fontsize=12)
ax.legend(frameon=False)
fig.tight_layout()
fig.savefig(f"{args.out}/heating.png", dpi=200)
print("wrote heating.csv, heating.png")