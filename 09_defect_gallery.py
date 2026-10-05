#!/usr/bin/env python3
"""Gallery of defect structures: top view + side view of the region around the impact point.

    python 09_defect_gallery.py                                   # the six reference cases
    python 09_defect_gallery.py Ar_random_E90_run11 Ar_random_E60_run6 --radius 6 --out gallery.png

Top row: bonds (cutoff 1.85 A) between the C atoms near the impact point. Rings that are not hexagons are shaded
(pentagon orange, heptagon blue, octagon and larger purple). Atoms with fewer than 3 neighbours have a red outline;
atoms more than 1.5 A out of the sheet plane are drawn in blue. Bottom row: side view (x against height).
Titles come from the 'outcome' column of impact_analysis.csv when present.
"""
import argparse
import json
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from ase.geometry import find_mic
from ase.io import read
from ase.neighborlist import neighbor_list
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon

DEFAULT = ["Ar_random_E60_run6", "Ar_random_E90_run11", "Ar_random_E70_run20",
           "Ar_random_E70_run15", "Ar_random_E100_run9", "Ar_random_E80_run1"]

REMARKS = {
    "Ar_random_E60_run6": "One atom knocked out of its site; its three neighbours are left under-coordinated. Ar is reflected.",
    "Ar_random_E90_run11": "One vacancy (open 9-ring with a pentagon) and the missing atom still hangs from the sheet by one bond.",
    "Ar_random_E70_run20": "A hole with a chain of atoms hanging below the sheet: the side view shows it dropping about 3 Å.",
    "Ar_random_E70_run15": "No atom lost: a 5-7 ring pair and one atom 1.6 Å below the plane. A Stone-Wales-like topological defect.",
    "Ar_random_E100_run9": "Open hole with a pentagon and an under-coordinated rim; two atoms sit about 2 Å out of plane. Not a clean vacancy.",
    "Ar_random_E80_run1": "An atom was pushed about 1.2 Å and came back: two under-coordinated atoms, but nothing is missing.",
}

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("tags", nargs="*", default=DEFAULT)
p.add_argument("--indir", default="results_v2")
p.add_argument("--workdir", default="bombardment_runs_v2")
p.add_argument("--csv", default="impact_analysis.csv")
p.add_argument("--radius", type=float, default=6.0)
p.add_argument("--cc-cut", type=float, default=1.85)
p.add_argument("--cols", type=int, default=3)
p.add_argument("--out", default="defect_gallery.png")
p.add_argument("--no-remarks", action="store_true", help="do not print the one-line remark under each panel")
args = p.parse_args()

GREY, RED, BLUE = "#8d949c", "#d6342c", "#2a78d6"
RING_COL = {5: "#eb6834", 7: "#2a78d6", 8: "#8a5cd0", 9: "#8a5cd0"}
plt.rcParams.update({"font.size": 10})

labels = {}
if Path(args.csv).exists():
    d = pd.read_csv(args.csv)
    for _, r in d.iterrows():
        labels[f"{r.projectile}_{r.target}_E{int(r.energy_eV)}_run{int(r.run)}"] = str(r.outcome)


def faces(adj, max_size=9):
    """Chordless rings (faces of the bond graph) up to max_size, each returned once as an ordered node list."""
    out = []
    for s in adj:
        stack = [(s, [s])]
        while stack:
            node, path = stack.pop()
            for nb in adj[node]:
                if nb == s and len(path) >= 3 and path[1] < path[-1]:
                    pset = set(path)
                    if all(len(adj[v] & pset) == 2 for v in path):
                        out.append(list(path))
                elif nb > s and nb not in path and len(path) < max_size:
                    stack.append((nb, path + [nb]))
    return out


def prepare(tag):
    final = read(Path(args.indir) / f"graphene_defective_{tag}.xyz")
    xy0 = np.array(json.loads((Path(args.workdir) / tag / "meta.json").read_text())["impact_xy"])
    is_c = np.array(final.get_chemical_symbols()) == "C"
    c_idx = np.where(is_c)[0]
    carbon = final[is_c]
    i, j = neighbor_list("ij", carbon, args.cc_cut)
    coord = np.bincount(i, minlength=len(carbon))
    z0 = float(np.median(carbon.positions[:, 2]))
    vec = carbon.positions.copy()
    vec[:, :2] -= xy0
    vec[:, 2] = 0.0
    vec, _ = find_mic(vec, carbon.cell, pbc=[True, True, False])
    r = np.linalg.norm(vec[:, :2], axis=1)
    keep = np.where(r <= args.radius)[0]
    pos = np.column_stack([vec[:, :2], carbon.positions[:, 2] - z0])
    return pos, coord, keep, i, j


def draw(axt, axs, axr, tag):
    pos, coord, keep, i, j = prepare(tag)
    ks = set(keep.tolist())
    adj = {int(k): set() for k in keep}
    segs, sides = [], []
    for a, b in zip(i, j):
        if a in ks and b in ks and a < b and np.linalg.norm(pos[a, :2] - pos[b, :2]) < 3.0:
            adj[int(a)].add(int(b)); adj[int(b)].add(int(a))
            segs.append([pos[a, :2], pos[b, :2]])
            sides.append([(pos[a, 0], pos[a, 2]), (pos[b, 0], pos[b, 2])])
    for ring in faces(adj):
        if len(ring) != 6 and len(ring) in RING_COL:
            axt.add_patch(Polygon(pos[ring, :2], closed=True, fc=RING_COL[len(ring)], ec="none", alpha=0.35, zorder=0))
    axt.add_collection(LineCollection(segs, colors="#444", linewidths=1.6, zorder=1))
    axs.add_collection(LineCollection(sides, colors="#444", linewidths=1.2, zorder=1))
    out = np.abs(pos[keep, 2]) > 1.5
    col = np.where(out, BLUE, GREY)
    edge = np.where(coord[keep] < 3, RED, "white")
    axt.scatter(pos[keep, 0], pos[keep, 1], s=150, c=col, edgecolors=edge, linewidths=2.0, zorder=2)
    axs.scatter(pos[keep, 0], pos[keep, 2], s=70, c=col, edgecolors=edge, linewidths=1.6, zorder=2)
    axt.plot(0, 0, marker="+", color="black", ms=9, mew=1.2, zorder=3)
    R = args.radius
    axt.set_xlim(-R - 0.4, R + 0.4); axt.set_ylim(-R - 0.4, R + 0.4); axt.set_aspect("equal"); axt.axis("off")
    zmin = min(-3.0, pos[keep, 2].min() - 0.5)
    zmax = max(1.5, pos[keep, 2].max() + 0.5)
    axs.set_xlim(-R - 0.4, R + 0.4); axs.set_ylim(zmin, zmax)
    axs.axhline(0, color="#bbb", lw=0.8, zorder=0)
    axs.set_xlabel("x (Å)", fontsize=8); axs.set_ylabel("height (Å)", fontsize=8)
    axs.tick_params(labelsize=8)
    for s in ("top", "right"):
        axs.spines[s].set_visible(False)
    e = int(tag.split("_E")[1].split("_")[0]); run = tag.split("_run")[1]
    outcome = labels.get(tag, "")
    short = outcome.split(" (")[0].split(" /")[0].strip()
    axt.set_title(f"{e} eV  ·  run {run}", fontsize=11, fontweight="bold", loc="left", pad=16)
    axt.text(0.0, 1.015, textwrap.fill(short, 34), transform=axt.transAxes, fontsize=9, color="#555", va="bottom", ha="left")
    if not args.no_remarks and tag in REMARKS:
        axr.text(0.0, 1.0, textwrap.fill(REMARKS[tag], 50), transform=axr.transAxes, fontsize=9, color="#222",
                 va="top", ha="left", style="italic")
    axr.axis("off")


n = len(args.tags)
cols = min(args.cols, n)
rows = int(np.ceil(n / cols))
H = 6.4 * rows + 1.6
fig = plt.figure(figsize=(4.8 * cols, H))
outer = fig.add_gridspec(rows, cols, hspace=0.12, wspace=0.10, left=0.04, right=0.99,
                         top=1 - 0.9 / H, bottom=1.1 / H)
for k, tag in enumerate(args.tags):
    r_, c_ = divmod(k, cols)
    inner = outer[r_, c_].subgridspec(3, 1, height_ratios=[3, 1.1, 0.7], hspace=0.32)
    axt = fig.add_subplot(inner[0]); axs = fig.add_subplot(inner[1]); axr = fig.add_subplot(inner[2])
    try:
        draw(axt, axs, axr, tag)
    except Exception as ex:
        axt.text(0.5, 0.5, f"{tag}\n{ex}", ha="center", va="center", fontsize=8, transform=axt.transAxes); axt.axis("off")
        print("failed:", tag, ex)

from matplotlib.lines import Line2D
from matplotlib.patches import Patch
handles = [
    Line2D([], [], marker="o", ls="", ms=10, mfc=GREY, mec="white", label="carbon atom in the sheet"),
    Line2D([], [], marker="o", ls="", ms=10, mfc=GREY, mec=RED, mew=2, label="fewer than 3 neighbours (dangling bonds)"),
    Line2D([], [], marker="o", ls="", ms=10, mfc=BLUE, mec="white", label="more than 1.5 Å out of the sheet plane"),
    Patch(fc=RING_COL[5], alpha=0.35, label="pentagon"),
    Patch(fc=RING_COL[7], alpha=0.35, label="heptagon"),
    Patch(fc=RING_COL[8], alpha=0.35, label="octagon or larger ring"),
    Line2D([], [], marker="+", ls="", ms=10, mec="black", mew=1.4, label="impact point"),
]
fig.legend(handles=handles, loc="lower center", ncol=4, frameon=False, fontsize=10, bbox_to_anchor=(0.5, 0.0))
fig.suptitle("Defects created by Ar impacts on graphene (top view and side view of the region around the impact)",
             x=0.02, ha="left", fontsize=13, y=0.99)
fig.savefig(args.out, dpi=200, bbox_inches="tight")
print("wrote", args.out)