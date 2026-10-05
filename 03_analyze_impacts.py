#!/usr/bin/env python3
"""Classify the outcome of each bombardment run and summarize defect yields per energy.

For every graphene_defective_<proj>_<target>_E<energy>_run<k>.xyz (the relaxed final structure) it reports:
  * projectile fate: reflected, transmitted (through the sheet), or trapped on/in the sheet,
    plus the closest approach of the projectile to the sheet plane (penetration depth),
  * carbon atoms ejected from the sheet (|z - z_sheet| > --eject Angstrom) = sputtered atoms = vacancies:
    1 -> single vacancy, 2 -> di-vacancy, >=3 -> multi-vacancy,
  * under-/over-coordinated carbon atoms (coordination != 3, cutoff --cc-cut) in the final frame,
  * transient bond breaking along the trajectory (bonds that break and heal during the impact),
  * carbon atoms permanently displaced in-plane by more than --disp Angstrom from the initial lattice,
  * distance from the impact point to the nearest lattice atom (when meta.json is available).

Then, per (projectile, target, energy), it writes a yield table: number of impacts, probability of
creating a vacancy (with a 95% Wilson interval), fractions of single/di/multi-vacancies, self-healed
pass-throughs, and projectile fates.

Limits: the structure is the relaxed end state of the run (impact + quench + minimization), not a long
anneal; ring statistics (Stone-Wales, 5-8-5) are not analysed; coordination counts are a diagnostic,
not a proof.

Usage (folders of 02_defect_engineering_stats.py):
  python 03_analyze_impacts.py --indir results_v2 --workdir bombardment_runs_v2
Older runs (no meta.json / initial.xyz):
  python 03_analyze_impacts.py --indir . --workdir none --structure graphene_32x32_relaxed.xyz
"""
import argparse
import csv
import json
import os
import re
from concurrent.futures import ProcessPoolExecutor
from collections import defaultdict
from pathlib import Path

import numpy as np
from ase.geometry import find_mic
from ase.io import iread, read
from ase.neighborlist import neighbor_list

PATTERN = re.compile(r"graphene_defective_(?P<proj>[A-Za-z]+)_(?P<target>[a-z]+)_E(?P<e>\d+(?:\.\d+)?)_run(?P<run>\d+)\.xyz$")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--structure", default="graphene_32x32_relaxed.xyz",
                   help="reference sheet, used when a run has no initial.xyz")
    p.add_argument("--indir", default="results_v2", help="folder with graphene_defective_*.xyz / bombardment_*.traj")
    p.add_argument("--workdir", default="bombardment_runs_v2",
                   help="folder with one subfolder per run (initial.xyz, meta.json); 'none' to disable")
    p.add_argument("--out", default="impact_analysis.csv")
    p.add_argument("--summary", default="yield_summary.csv")
    p.add_argument("--disp", type=float, default=0.8, help="in-plane displacement threshold (A) for 'displaced' C atoms")
    p.add_argument("--eject", type=float, default=3.0, help="distance from the sheet plane (A) for 'ejected' C atoms")
    p.add_argument("--cc-cut", type=float, default=1.85, help="C-C bonding cutoff (A)")
    p.add_argument("--transient-min", type=int, default=0,
                   help="min. number of transiently under-coordinated C atoms to report 'self-healed'")
    p.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) // 2),
                   help="analyse this many runs in parallel (default: half of the CPU threads)")
    p.add_argument("--reanalyse", action="store_true",
                   help="ignore the cached results in --out and analyse every run again")
    p.add_argument("--no-traj", action="store_true", help="skip the trajectory (no closest-approach/transient columns)")
    return p.parse_args()


def scan_trajectory(traj_path, args, stride=5):
    """Stream the trajectory frame by frame (low memory) and return two numbers:

    * lowest height of the projectile above the sheet plane (A); the plane is the median z of the
      carbon atoms in each frame, so a drift of the whole sheet is ignored;
    * max number of under-coordinated C atoms (coordination < 3) in any sampled frame. Atoms that
      left the sheet plane (|z - median z| > --eject) are ignored. A positive value with a clean final
      frame means bonds broke and healed during the impact.
    """
    lowest, transient = np.inf, 0
    for k, f in enumerate(iread(traj_path, index=":")):
        sym = np.array(f.get_chemical_symbols())
        is_c = sym == "C"
        z_plane = float(np.median(f.positions[is_c, 2]))
        lowest = min(lowest, float(f.positions[-1, 2]) - z_plane)
        if k % stride == 0:
            carbon = f[is_c]
            i, j_nb = neighbor_list("ij", carbon, args.cc_cut)
            coord = np.bincount(i, minlength=len(carbon))
            stay = np.abs(carbon.positions[:, 2] - z_plane) <= args.eject
            transient = max(transient, int(((coord < 3) & stay).sum()))
    return lowest, transient


def analyse(final, pristine, args):
    sym = np.array(final.get_chemical_symbols())
    is_c = sym == "C"
    # sheet plane of the final structure (robust to a few ejected atoms and to a drift of the whole sheet)
    z_sheet = float(np.median(final.positions[is_c, 2]))
    proj_idx = int(np.where(~is_c)[0][0])
    n_c = int(is_c.sum())

    # --- carbon displacements relative to the initial lattice (minimum image in x, y)
    disp_vec = final.positions[is_c] - pristine.positions[:n_c]
    disp_vec, _ = find_mic(disp_vec, pristine.cell, pbc=[True, True, False])
    # remove the rigid drift of the free sheet (it absorbs the projectile momentum); median is robust to ejected atoms
    drift = np.median(disp_vec, axis=0)
    disp_vec = disp_vec - drift
    disp = np.linalg.norm(disp_vec, axis=1)
    disp_xy = np.linalg.norm(disp_vec[:, :2], axis=1)  # in-plane only: an out-of-plane bump is not a defect

    # --- ejected carbon atoms (left the sheet plane)
    dz = final.positions[is_c, 2] - z_sheet
    ejected = np.abs(dz) > args.eject
    n_ejected = int(ejected.sum())

    # --- coordination of the carbon atoms that stay in the sheet
    carbon = final[is_c]
    i, j_nb = neighbor_list("ij", carbon, args.cc_cut)
    coord = np.bincount(i, minlength=n_c)
    # A C atom is MISSING from the sheet if it has fewer than 2 C neighbours (free: 0 neighbours; bonded adatom: 1
    # neighbour, e.g. an atom hanging by one bond) or if it left the sheet plane (|dz| > --eject).
    # Height alone is not used to decide: a local dimple or bulge moves whole groups of bonded atoms by > 1 A.
    missing = (coord <= 1) | ejected
    # Peel dangling chains: an atom whose only remaining bonds lead to missing atoms is itself hanging off the sheet
    # (e.g. a 3-atom chain hanging below a hole), so repeat until nothing changes.
    while True:
        keep_edge = ~missing[i] & ~missing[j_nb]
        c_kept = np.bincount(i[keep_edge], minlength=n_c)
        newly = (~missing) & (c_kept <= 1)
        if not newly.any():
            break
        missing = missing | newly
    n_vacant = int(missing.sum())                      # empty lattice sites
    n_free = int((missing & (coord == 0)).sum())       # atom has left, no bond to the sheet
    n_adatom = int((missing & (coord >= 1)).sum())     # atom still bonded to the sheet by one bond
    kept = ~missing
    # --- ring statistics of the carbon network that stays in the sheet (non-hexagonal rings = topological defects)
    adj = {int(k): set() for k in np.where(kept)[0]}
    for a, b in zip(i, j_nb):
        if kept[a] and kept[b]:
            adj[int(a)].add(int(b))
    rings = ring_signature(ring_counts(adj))
    n_under = int(((coord < 3) & kept).sum())
    n_over = int(((coord > 3) & kept).sum())

    # --- projectile fate
    z_p = float(final.positions[proj_idx, 2]) - z_sheet
    d_min = float(final.get_distances(proj_idx, np.where(is_c)[0], mic=True).min())
    fate = "reflected" if z_p > 3.0 else ("transmitted" if z_p < -2.0 else "trapped")

    # --- outcome label
    if n_vacant == 0:
        if n_under >= 4:
            outcome = (f"extended defect (open hole, {n_under} under-coordinated atoms"
                       + (f", rings {rings})" if rings else ")"))
        elif rings:
            outcome = f"topological defect ({rings} rings, no atom missing" + (f", {n_under} under-coordinated" if n_under else "") + ")"
        elif n_under or n_over:
            outcome = "bond defect / buckled region (no atom missing from the sheet" + (f", rings {rings})" if rings else ")")
        elif rings:
            outcome = f"topological defect ({rings} rings, no atom missing)"
        elif (disp_xy > args.disp).any():
            outcome = "in-plane rearrangement, no coordination or ring defect"
        else:
            outcome = "elastic / no permanent defect"
    else:
        if n_vacant == 1:
            outcome = "single vacancy" if n_free == 1 else "Frenkel pair (vacancy + bonded adatom)"
        elif n_vacant == 2:
            outcome = "di-vacancy"
        else:
            outcome = f"multi-vacancy ({n_vacant} C)"
        if n_vacant >= 2 and n_adatom:
            outcome += f" ({n_adatom} atom{'s' if n_adatom > 1 else ''} still bonded as adatom)"

    return dict(z_proj_final=round(z_p, 2), closest_C=round(d_min, 2), fate=fate,
                n_ejected=n_ejected, n_vacant=n_vacant, n_adatom=n_adatom, n_undercoord=n_under, n_overcoord=n_over,
                rings=rings, n_displaced=int((disp_xy > args.disp).sum()), max_disp=round(float(disp.max()), 2),
                outcome=outcome), z_sheet


def ring_counts(adj, max_size=9):
    """Count chordless rings of size 3..max_size in a bond graph (dict: node -> set of neighbours).

    Each ring is found once: the start node is the smallest index in the ring and the traversal
    direction is fixed. Chordless rings of a planar sp2 network are its faces, so pristine graphene
    gives only 6-rings, a relaxed single vacancy gives 5+9, a Stone-Wales defect 5-7-7-5, etc.
    """
    counts = defaultdict(int)
    for s in adj:
        stack = [(s, [s])]
        while stack:
            node, path = stack.pop()
            for nb in adj[node]:
                if nb == s and len(path) >= 3 and path[1] < path[-1]:
                    pset = set(path)
                    if all(len(adj[v] & pset) == 2 for v in path):   # chordless
                        counts[len(path)] += 1
                elif nb > s and nb not in path and len(path) < max_size:
                    stack.append((nb, path + [nb]))
    return dict(counts)


def ring_signature(counts):
    """e.g. {5: 2, 7: 2} -> '5-5-7-7' ; empty if only hexagons."""
    odd = [k for k in sorted(counts) if k != 6 for _ in range(counts[k])]
    return "-".join(map(str, odd))


def wilson(k, n, z=1.96):
    """95% Wilson score interval for a binomial proportion."""
    if n == 0:
        return 0.0, 0.0
    p = k / n
    denom = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / denom
    half = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / denom
    return max(0.0, centre - half), min(1.0, centre + half)


def summarize(rows):
    groups = defaultdict(list)
    for r in rows:
        groups[(r["projectile"], r["target"], r["energy_eV"])].append(r)
    out = []
    for (proj, target, e), rs in sorted(groups.items()):
        n = len(rs)
        vac = sum(r["n_ejected"] >= 1 for r in rs)
        lo, hi = wilson(vac, n)
        vsite = sum(r["n_vacant"] >= 1 for r in rs)
        slo, shi = wilson(vsite, n)
        n_defect = sum(not (r["outcome"].startswith("elastic") or r["outcome"].startswith("transient")) for r in rs)
        out.append(dict(
            projectile=proj, target=target, energy_eV=e, n_impacts=n,
            P_vacancy=round(vac / n, 3), P_vac_lo95=round(lo, 3), P_vac_hi95=round(hi, 3),
            P_vac_site=round(vsite / n, 3), P_vac_site_lo95=round(slo, 3), P_vac_site_hi95=round(shi, 3),
            f_frenkel=round(sum(r["outcome"].startswith("Frenkel") for r in rs) / n, 3),
            P_any_defect=round(n_defect / n, 3),
            f_single=round(sum(r["n_vacant"] == 1 for r in rs) / n, 3),
            f_divacancy=round(sum(r["n_vacant"] == 2 for r in rs) / n, 3),
            f_multi=round(sum(r["n_vacant"] > 2 for r in rs) / n, 3),
            f_topological=round(sum(r["outcome"].startswith("topological") for r in rs) / n, 3),
            f_bond_defect=round(sum(r["outcome"].startswith("bond defect") for r in rs) / n, 3),
            f_extended=round(sum(r["outcome"].startswith("extended defect") for r in rs) / n, 3),
            f_self_healed=round(sum(r["outcome"].startswith("transient") for r in rs) / n, 3),
            f_reflected=round(sum(r["fate"] == "reflected" for r in rs) / n, 3),
            f_transmitted=round(sum(r["fate"] == "transmitted" for r in rs) / n, 3),
            f_trapped=round(sum(r["fate"] == "trapped" for r in rs) / n, 3),
            mean_C_ejected=round(float(np.mean([r["n_ejected"] for r in rs])), 2)))
    return out


ANALYSIS_VERSION = "5"   # bump whenever the classification logic changes: cached rows of other versions are redone

NUMERIC = {"energy_eV": float, "run": int, "d_impact_C": float, "min_height": float, "z_proj_final": float,
           "closest_C": float, "n_ejected": int, "n_undercoord": int, "max_transient_undercoord": int,
           "n_overcoord": int, "n_vacant": int, "n_adatom": int, "n_displaced": int, "max_disp": float}


def load_cache(path):
    """Rows of a previous analysis, keyed by run tag (empty if the file is missing or from an older version)."""
    if not Path(path).is_file():
        return {}
    cache = {}
    with open(path, newline="") as fh:
        rd = csv.DictReader(fh)
        if rd.fieldnames is None or "analysis_version" not in rd.fieldnames:
            return {}                      # written by an older version: redo everything
        for r in rd:
            if r.get("analysis_version") != ANALYSIS_VERSION:
                continue                   # produced by a different version of the analysis: redo it
            for k, conv in NUMERIC.items():
                v = r.get(k, "")
                r[k] = conv(v) if v not in ("", None) else ""
            cache[f"{r['projectile']}_{r['target']}_E{r['energy_eV']:g}_run{r['run']}"] = r
    return cache


def process_one(task):
    """Analyse one run (executed in a worker process); returns the result row."""
    xyz, args = task
    m = PATTERN.search(xyz.name)
    indir = Path(args.indir)
    workdir = None if args.workdir.lower() == "none" else Path(args.workdir)
    default_ref = read(args.structure) if Path(args.structure).is_file() else None
    tag = f"{m['proj']}_{m['target']}_E{m['e']}_run{m['run']}"
    rundir = workdir / tag if workdir else None
    ref_path = rundir / "initial.xyz" if rundir else None
    if ref_path is not None and ref_path.is_file():
        pristine = read(ref_path)
    elif default_ref is not None:
        pristine = default_ref
    else:
        raise SystemExit(f"No reference lattice for {tag}: pass --structure")
    res, _ = analyse(read(xyz), pristine, args)

    res["d_impact_C"] = ""
    if rundir is not None and (rundir / "meta.json").is_file():
        xy = np.array(json.loads((rundir / "meta.json").read_text())["impact_xy"])
        dvec = pristine.positions[:, :2] - xy
        dvec3 = np.column_stack([dvec, np.zeros(len(dvec))])
        dvec3, _ = find_mic(dvec3, pristine.cell, pbc=[True, True, False])
        res["d_impact_C"] = round(float(np.linalg.norm(dvec3[:, :2], axis=1).min()), 2)

    traj = indir / f"bombardment_{tag}.traj"
    res["min_height"], res["max_transient_undercoord"] = "", ""
    if not args.no_traj and traj.exists():
        try:
            lowest, transient = scan_trajectory(traj, args)
            res["min_height"] = round(lowest, 2)
            res["max_transient_undercoord"] = transient
            if res["outcome"].startswith("elastic") and transient > args.transient_min:
                res["outcome"] = "transient bond breaking, self-healed"
        except Exception as exc:  # damaged or truncated trajectory: keep the final-structure result
            print(f"WARNING: could not read {traj.name} ({type(exc).__name__}); trajectory columns left empty. "
                  f"Rebuild it from {workdir}/{tag}/traj.dump if needed.", flush=True)
            res["outcome"] += " [trajectory unreadable]"
    return dict(projectile=m["proj"], energy_eV=float(m["e"]), target=m["target"], run=int(m["run"]),
                analysis_version=ANALYSIS_VERSION, **res)


def main():
    args = parse_args()
    indir = Path(args.indir)
    cache = {} if (args.reanalyse or args.no_traj) else load_cache(args.out)
    rows, todo = [], []
    for xyz in sorted(indir.glob("graphene_defective_*.xyz")):
        m = PATTERN.search(xyz.name)
        if not m:
            continue
        tag = f"{m['proj']}_{m['target']}_E{m['e']}_run{m['run']}"
        if tag in cache and "unreadable" not in cache[tag]["outcome"] and cache[tag]["min_height"] != "":
            rows.append(cache[tag])
        else:
            todo.append(xyz)
    print(f"{len(rows)} runs reused from {args.out}, {len(todo)} to analyse with {args.jobs} worker(s)", flush=True)
    if todo:
        with ProcessPoolExecutor(max_workers=args.jobs) as ex:
            for n, row in enumerate(ex.map(process_one, [(x, args) for x in todo]), 1):
                rows.append(row)
                print(f"  analysed {n}/{len(todo)}: {row['projectile']} E{row['energy_eV']:g} run{row['run']}", flush=True)

    rows.sort(key=lambda r: (r["projectile"], r["target"], r["energy_eV"], r["run"]))
    cols = ["projectile", "energy_eV", "target", "run", "d_impact_C", "fate", "min_height", "z_proj_final",
            "closest_C", "n_ejected", "n_vacant", "n_adatom", "n_undercoord", "max_transient_undercoord", "n_overcoord", "rings", "n_displaced",
            "max_disp", "outcome", "analysis_version"]
    header = ["proj", "E(eV)", "target", "run", "d_imp(A)", "fate", "min_h(A)", "z_end(A)", "d_C(A)", "ejec", "vac", "adat",
              "under", "transient", "over", "rings", "displ", "max_d(A)", "outcome"]
    print(" | ".join(header))
    for r in rows:
        print(" | ".join(f"{r[c]:g}" if isinstance(r[c], float) else str(r[c]) for c in cols if c != "analysis_version"))
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows([{c: r[c] for c in cols} for r in rows])
    print(f"\n{len(rows)} runs analysed -> {args.out}")

    summary = summarize(rows)
    if summary:
        print("\nYield summary (P_vac = probability that >= 1 C atom is ejected, 95% Wilson interval):")
        print("proj | target | E(eV) | N | P_vac(ejected) [95% CI] | P_vac(empty site) [95% CI] | frenkel | any defect"
              " | single | di | multi | topo | healed | refl | trans | trap")
        for s in summary:
            print(f"{s['projectile']} | {s['target']} | {s['energy_eV']:g} | {s['n_impacts']} | "
                  f"{s['P_vacancy']:.2f} [{s['P_vac_lo95']:.2f}-{s['P_vac_hi95']:.2f}] | "
                  f"{s['P_vac_site']:.2f} [{s['P_vac_site_lo95']:.2f}-{s['P_vac_site_hi95']:.2f}] | "
                  f"{s['f_frenkel']:.2f} | {s['P_any_defect']:.2f} | {s['f_single']:.2f} | "
                  f"{s['f_divacancy']:.2f} | {s['f_multi']:.2f} | {s['f_topological']:.2f} | {s['f_self_healed']:.2f} | "
                  f"{s['f_reflected']:.2f} | {s['f_transmitted']:.2f} | {s['f_trapped']:.2f}")
        with open(args.summary, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(summary[0]))
            w.writeheader()
            w.writerows(summary)
        print(f"\nYield table -> {args.summary}")


if __name__ == "__main__":
    main()