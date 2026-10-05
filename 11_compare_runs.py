#!/usr/bin/env python3
"""Compare two analysis tables: baseline vs a test run (bigger sheet, longer quench, longer window, other cutoffs...).

    python 11_compare_runs.py impact_analysis.csv conv_big/impact_analysis.csv --label "4608 atoms"
    python 11_compare_runs.py impact_analysis.csv conv_quench5/impact_analysis.csv --label "quench 5 ps" --paired

Per energy it prints P(vacancy) and the self-healed fraction with 95% Wilson intervals for both tables, and says
whether the intervals overlap. With --paired (same seeds, so the same impact point and thermal state, as in the
quench and time-window tests) it also counts how many individual runs keep the same outcome class and lists the runs
that changed, which is a much sharper test than comparing the two proportions.
"""
import argparse

import numpy as np
import pandas as pd

p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
p.add_argument("baseline")
p.add_argument("test")
p.add_argument("--label", default="test")
p.add_argument("--paired", action="store_true")
args = p.parse_args()


def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0
    p_ = k / n
    d = 1 + z * z / n
    c = (p_ + z * z / (2 * n)) / d
    h = z * np.sqrt(p_ * (1 - p_) / n + z * z / (4 * n * n)) / d
    return max(0.0, c - h), min(1.0, c + h)


def prep(path):
    d = pd.read_csv(path)
    d = d[d.target == "random"].copy()
    d["vac"] = d.n_vacant >= 1
    d["healed"] = d.outcome.astype(str).str.startswith("transient")
    d["cls"] = np.where(d.vac, np.where(d.n_vacant >= 2, "di-vacancy", "single vacancy"),
                        np.where(d.healed, "self-healed", np.where(d.outcome.astype(str).str.startswith("elastic"),
                                                                   "elastic", "other defect")))
    return d


a, b = prep(args.baseline), prep(args.test)
energies = sorted(set(a.energy_eV) & set(b.energy_eV))
print(f"baseline: {args.baseline}\n{args.label}: {args.test}\n")
print(f"{'E(eV)':>6} | {'N':>3} {'P_vac base':>22} | {'N':>3} {'P_vac ' + args.label:>22} | overlap | {'healed base':>11} {'healed test':>11}")
bad = 0
for e in energies:
    ga, gb = a[a.energy_eV == e], b[b.energy_eV == e]
    if args.paired:
        common = sorted(set(ga.run) & set(gb.run))
        ga, gb = ga[ga.run.isin(common)], gb[gb.run.isin(common)]
    ka, na, kb, nb = int(ga.vac.sum()), len(ga), int(gb.vac.sum()), len(gb)
    la, ha = wilson(ka, na)
    lb, hb = wilson(kb, nb)
    ok = max(la, lb) <= min(ha, hb)
    bad += (not ok)
    print(f"{e:6.0f} | {na:3d} {ka/na:6.2f} [{la:.2f}-{ha:.2f}]      | {nb:3d} {kb/nb:6.2f} [{lb:.2f}-{hb:.2f}]      | "
          f"{'yes' if ok else 'NO ':>7} | {ga.healed.mean():11.2f} {gb.healed.mean():11.2f}")
print("\nVerdict:", "all vacancy probabilities agree within 95% intervals" if bad == 0
      else f"{bad} energy(ies) differ beyond the intervals: investigate")
if args.paired:
    m = a.merge(b, on=["energy_eV", "run"], suffixes=("_a", "_b"))
    same = (m.cls_a == m.cls_b)
    print(f"\npaired runs: {len(m)}; same outcome class: {int(same.sum())} ({same.mean():.0%}); "
          f"same vacancy yes/no: {int((m.vac_a == m.vac_b).sum())} ({(m.vac_a == m.vac_b).mean():.0%})")
    diff = m[~same][["energy_eV", "run", "cls_a", "cls_b"]]
    if len(diff):
        print("runs whose class changed:\n" + diff.to_string(index=False))