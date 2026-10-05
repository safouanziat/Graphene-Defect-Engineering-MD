#!/usr/bin/env python3
"""Plot the yield curves from yield_summary.csv (written by 03_analyze_impacts.py).

    python 04_plot_yield.py                       # reads yield_summary.csv -> yield_curve.png
    python 04_plot_yield.py --csv yield_summary.csv --projectile Ar --target random --out figures/yield_curve.png
"""
import argparse
import csv
from pathlib import Path

import matplotlib
matplotlib.use("Agg")           # write the PNG without needing a display
import matplotlib.pyplot as plt
import numpy as np

p = argparse.ArgumentParser()
p.add_argument("--csv", default="yield_summary.csv")
p.add_argument("--projectile", default="Ar")
p.add_argument("--target", default="random", help="impact mode: random (statistical) / top / bridge")
p.add_argument("--out", default="yield_curve.png")
p.add_argument("--vacancy", choices=["site", "ejected"], default="site",
               help="site: vacancy = empty lattice site (includes atoms pushed through the sheet); "
                    "ejected: only atoms that left the sheet by more than 3 A")
p.add_argument("--title", default=None)
p.add_argument("--show", action="store_true", help="also open an interactive window")
args = p.parse_args()

# ---- read the table (the CSV stores everything as text) ----
rows = [r for r in csv.DictReader(open(args.csv))
        if r["projectile"] == args.projectile and r["target"] == args.target]
rows.sort(key=lambda r: float(r["energy_eV"]))
if not rows:
    raise SystemExit(f"no rows for projectile={args.projectile}, target={args.target} in {args.csv}")
col = lambda name: np.array([float(r[name]) for r in rows])
E, N = col("energy_eV"), col("n_impacts").astype(int)
if args.vacancy == "site" and "P_vac_site" in rows[0]:
    P, LO, HI = col("P_vac_site"), col("P_vac_site_lo95"), col("P_vac_site_hi95")
    vac_label = "Any vacancy (empty lattice site, 95% Wilson)"
else:
    P, LO, HI = col("P_vacancy"), col("P_vac_lo95"), col("P_vac_hi95")
    vac_label = "Vacancy, atom ejected (95% Wilson)"
di, heal, trans = col("f_divacancy"), col("f_self_healed"), col("f_transmitted")

# ---- style ----
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"      # colour-blind-safe trio
INK, SEC, GRID, BG = "#0b0b0b", "#52514e", "#e6e5e1", "#fcfcfb"
plt.rcParams.update({"text.color": INK, "axes.labelcolor": SEC, "xtick.color": SEC,
                     "ytick.color": SEC, "axes.edgecolor": GRID})

fig, (a, b) = plt.subplots(1, 2, figsize=(11, 4.4), facecolor=BG)
for ax in (a, b):
    ax.set_facecolor(BG)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.set_xticks(E)
    ax.set_xticklabels([f"{e:g}" if (e % 10 == 0 and e not in (110, 130, 140, 160, 170, 180, 190)) or e in (150,) else "" for e in E])
    ax.set_ylim(0, 1.05)
    ticks = np.arange(0, 1.01, 0.25)
    ax.set_yticks(ticks)
    ax.set_yticklabels([f"{int(t * 100)}%" for t in ticks])

# left: vacancy probability with its 95% interval (the bars show LO..HI, not +/- a symmetric error)
a.errorbar(E, P, yerr=[P - LO, HI - P], fmt="o-", color=BLUE, lw=2, ms=8, mec=BG, mew=2,
           capsize=4, elinewidth=1.5, label=vac_label)
a.plot(E, di, "s--", color=ORANGE, lw=2, ms=7, mec=BG, mew=2, label="Di-vacancy")
a.set_xlabel(f"{args.projectile} kinetic energy (eV)")
a.set_ylabel("Fraction of impacts")
a.set_title("Vacancy probability per impact", loc="left", fontsize=12, color=INK)
a.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.27), ncol=2, fontsize=9)

# right: what happens to the projectile and to the sheet
b.plot(E, trans, "o-", color=BLUE, lw=2, ms=8, mec=BG, mew=2, label=f"{args.projectile} transmitted")
b.plot(E, heal, "s-", color=AQUA, lw=2, ms=7, mec=BG, mew=2, label="Self-healed bond breaking")
b.set_xlabel(f"{args.projectile} kinetic energy (eV)")
b.set_title("What happens to projectile and sheet", loc="left", fontsize=12, color=INK)
b.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.27), ncol=2, fontsize=9)

fig.suptitle(args.title or f"{args.projectile} on free-standing graphene, AIREBO+ZBL, 300 K, {args.target} impact points",
             x=0.01, ha="left", fontsize=11, color=SEC, y=0.99)
groups, start = [], 0
for i in range(1, len(E) + 1):
    if i == len(E) or N[i] != N[start]:
        lo_e, hi_e = E[start], E[i - 1]
        groups.append(f"{N[start]} at {lo_e:g} eV" if lo_e == hi_e else f"{N[start]} at {lo_e:g}-{hi_e:g} eV")
        start = i
fig.text(0.01, 0.01, "Impacts per energy: " + ", ".join(groups) + ".  Error bars: 95% Wilson interval.",
         fontsize=9, color=SEC, ha="left")
fig.tight_layout(rect=(0, 0.04, 1, 0.95))
Path(args.out).parent.mkdir(parents=True, exist_ok=True)
fig.savefig(args.out, dpi=200, facecolor=BG)
print("saved", args.out)
if args.show:
    matplotlib.use("TkAgg", force=True); plt.show()