#!/usr/bin/env python3
"""Two physics figures from impact_analysis.csv.

    python 07_plot_physics.py

1. outcome_map.png    every impact as a point: x = distance from the impact point to the nearest C atom,
                      y = Ar energy, colour = what happened. Shows WHERE you have to hit to make a vacancy.
2. cross_section.png  effective vacancy radius r_eff(E) = sqrt(P_vac * A_atom / pi), where A_atom is the area per
                      carbon atom. P_vac is the chance that a random impact lands inside a disc of radius r_eff
                      around some atom, so r_eff is the size of the 'target' each atom presents. 95% Wilson intervals.
"""
import argparse

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("--csv", default="impact_analysis.csv")
p.add_argument("--a0", type=float, default=2.419, help="lattice constant used for the sheet (A)")
p.add_argument("--target", default="random")
p.add_argument("--out", default=".")
args = p.parse_args()

BLUE, ORANGE, GREEN, PURPLE, GREY = "#2a78d6", "#eb6834", "#1baf7a", "#8a5cd0", "#9aa0a6"
plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})

df = pd.read_csv(args.csv)
df = df[df.target == args.target].copy()


def category(r):
    o = str(r["outcome"])
    if r["n_vacant"] >= 2:
        return "di-vacancy"
    if r["n_vacant"] == 1:
        return "single vacancy"
    if o.startswith(("bond defect", "topological", "extended")):
        return "other defect (no atom lost)"
    if o.startswith("transient"):
        return "self-healed"
    return "elastic"


df["cat"] = df.apply(category, axis=1)
style = {"single vacancy": (BLUE, "o"), "di-vacancy": (ORANGE, "s"), "other defect (no atom lost)": (PURPLE, "^"),
         "self-healed": (GREEN, "D"), "elastic": (GREY, "x")}


def wilson(k, n, z=1.96):
    p_ = k / n
    d = 1 + z * z / n
    c = (p_ + z * z / (2 * n)) / d
    h = z * np.sqrt(p_ * (1 - p_) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


A = np.sqrt(3) / 4 * args.a0 ** 2            # area per carbon atom (A^2)
r_full = np.sqrt(A / np.pi)                   # radius at which the discs tile the sheet
bond = args.a0 / np.sqrt(3)                   # C-C bond length

E = sorted(df.energy_eV.unique())
rows = []
for e in E:
    g = df[df.energy_eV == e]
    k = int((g.n_vacant >= 1).sum())
    n = len(g)
    lo, hi = wilson(k, n)
    f = lambda q: np.sqrt(q * A / np.pi)
    rows.append((e, n, k / n, f(k / n), f(lo), f(hi)))
res = pd.DataFrame(rows, columns=["E", "N", "P", "r", "r_lo", "r_hi"])

# ---- figure 1: outcome map
fig, ax = plt.subplots(figsize=(8.5, 5.2))
rng = np.random.default_rng(1)
pos = {e: i for i, e in enumerate(E)}
for cat, (c, m) in style.items():
    g = df[df.cat == cat]
    if g.empty:
        continue
    y = g.energy_eV.map(pos) + rng.uniform(-0.22, 0.22, len(g))
    ax.scatter(g.d_impact_C, y, c=c, marker=m, s=46, label=f"{cat} ({len(g)})", linewidths=1.6,
               edgecolors="white" if m != "x" else None)
ax.plot(res.r, [pos[e] for e in res.E], color="black", lw=1.4, ls="--", label="r_eff from P_vac")
ax.axvline(bond / 2, color="black", lw=0.8, ls=":")
ax.set_ylim(-0.6, len(E) - 0.1)
ax.text(bond / 2 + 0.01, len(E) - 0.2, "half a C-C bond", va="top", ha="left", fontsize=9)
ax.set_yticks(range(len(E)))
ax.set_yticklabels([f"{e:.0f}" for e in E])
ax.set_xlabel("distance from impact point to the nearest C atom (Å)")
ax.set_ylabel("Ar kinetic energy (eV)")
ax.set_title("Where the ion lands decides the outcome", loc="left", fontsize=13)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.16), ncol=3, frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(f"{args.out}/outcome_map.png", dpi=200)

# ---- figure 2: effective radius
fig, ax = plt.subplots(figsize=(7.2, 4.6))
ax.errorbar(res.E, res.r, yerr=[res.r - res.r_lo, res.r_hi - res.r], fmt="o-", color=BLUE, capsize=4, lw=2,
            label="effective vacancy radius")
ax.axhline(r_full, color="black", lw=1, ls="--")
ax.text(res.E.min(), r_full + 0.02, "discs tile the sheet: every impact makes a vacancy", fontsize=9, va="bottom")
ax.axhline(bond / 2, color=GREY, lw=1, ls=":")
ax.text(res.E.min(), bond / 2 + 0.02, "half a C-C bond", fontsize=9, color="#555", va="bottom")
last = res.iloc[-1]
if last.P >= 1.0:
    ax.annotate("lower bound\n(all impacts hit)", (last.E, last.r), xytext=(last.E - 38, last.r - 0.33), fontsize=9,
                arrowprops=dict(arrowstyle="->", lw=0.8))
ax.set_ylim(0, 1.1)
ax.set_xlabel("Ar kinetic energy (eV)")
ax.set_ylabel("r_eff (Å)")
ax.set_title("How big a target is each carbon atom?", loc="left", fontsize=13)
fig.tight_layout()
fig.savefig(f"{args.out}/cross_section.png", dpi=200)

# ---- figure 3: outcome composition per energy (analogue of the probability-vs-energy plot of Bellido & Seminario)
order = ["elastic", "self-healed", "other defect (no atom lost)", "single vacancy", "di-vacancy"]
cols = {"elastic": GREY, "self-healed": GREEN, "other defect (no atom lost)": PURPLE, "single vacancy": BLUE, "di-vacancy": ORANGE}
comp = df.groupby(["energy_eV", "cat"]).size().unstack(fill_value=0).reindex(columns=order, fill_value=0)
frac = comp.div(comp.sum(axis=1), axis=0)
fig, ax = plt.subplots(figsize=(8, 4.6))
x = np.arange(len(frac))
bottom = np.zeros(len(frac))
for c in order:
    ax.bar(x, frac[c], bottom=bottom, color=cols[c], label=c, width=0.72, edgecolor="white", linewidth=0.8)
    bottom += frac[c].values
ax.set_xticks(x)
ax.set_xticklabels([f"{e:.0f}\n(N={int(n)})" for e, n in zip(frac.index, comp.sum(axis=1))], fontsize=9)
ax.set_xlabel("Ar kinetic energy (eV)")
ax.set_ylabel("fraction of impacts")
ax.set_title("What a random impact does, vs energy", loc="left", fontsize=13)
ax.legend(loc="upper center", bbox_to_anchor=(0.5, -0.2), ncol=3, frameon=False, fontsize=9)
fig.tight_layout()
fig.savefig(f"{args.out}/outcome_composition.png", dpi=200)

print(res.round(3).to_string(index=False))
print(f"\nA_atom = {A:.3f} A^2, full-coverage radius = {r_full:.3f} A, bond = {bond:.3f} A")
print("outcome_map.png, cross_section.png, outcome_composition.png written")