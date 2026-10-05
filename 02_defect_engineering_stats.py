#!/usr/bin/env python3
"""Statistical defect engineering of graphene by projectile bombardment (native LAMMPS dynamics).

Each impact is a complete, reproducible mini-experiment run inside one LAMMPS process:

  0. (optional) rescale the sheet to the equilibrium lattice constant of the potential
     and minimize it, so the support is not strained;
  1. thermalize the sheet at T (Nose-Hoover NVT, projectile held still above it);
  2. impact in NVE: the projectile gets its kinetic energy, the collision runs with a short step;
  3. quench: Langevin thermostat on the sheet at T (removes the impact heat), then a
     minimization, so metastable defects (e.g. a divacancy) can reconstruct;
  4. the relaxed final structure is written for the analysis script.

The impact point is deterministic ('top', 'bridge') or random ('random', uniform over the
periodic sheet) with a reproducible seed per (projectile, energy, run). Run many impacts per
energy (--impacts N) to get vacancy probabilities instead of single yes/no outcomes.

Physical models:
  * Ar : hybrid AIREBO + ZBL (short-range repulsion), units metal
  * N  : ReaxFF with charge equilibration (qeq/reaxff), units real. The ReaxFF parameter set
         MUST be validated for graphene first (python check_potentials.py); ffield.reax.CHN
         from the LAMMPS examples describes nitramines (RDX) and is NOT suitable.

Examples:
  python 02_defect_engineering_stats.py --parallel 14 --energies 40 44 48 52 56 60 --targets random --impacts 20
  python 02_defect_engineering_stats.py --parallel 14 --energies 50 --targets top bridge --impacts 1 --timestep 0.01
  python 02_defect_engineering_stats.py --dry-run --energies 50 --impacts 2
"""
import argparse
import csv
import json
import os
import shlex
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from ase import Atom, units
from ase.data import atomic_masses, atomic_numbers
from ase.io import read, write

TARGET_MODES = ("top", "bridge", "random")
A0_AIREBO = 2.419  # graphene lattice constant of AIREBO (Angstrom), from check_potentials.py


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--structure", default="graphene_32x32_relaxed.xyz", help="graphene sheet")
    p.add_argument("--projectiles", nargs="+", default=["Ar"], choices=["N", "Ar"])
    p.add_argument("--energies", nargs="+", type=float, default=[20, 30, 40, 45, 50, 55, 60, 70, 80],
                   help="projectile kinetic energies in eV")
    p.add_argument("--targets", nargs="+", default=["random"], choices=TARGET_MODES)
    p.add_argument("--impacts", type=int, default=10, help="impacts per (projectile, energy, target)")
    p.add_argument("--seed", type=int, default=2026, help="base random seed (impact points and velocities)")
    p.add_argument("--temperature", type=float, default=300.0, help="sheet temperature (K)")
    p.add_argument("--timestep", type=float, default=0.05, help="impact time step in fs")
    p.add_argument("--time-fs", type=float, default=750.0, help="simulated time of the impact stage (fs)")
    p.add_argument("--therm-ps", type=float, default=2.0, help="thermalization time (ps)")
    p.add_argument("--quench-ps", type=float, default=2.0, help="post-impact Langevin quench time (ps)")
    p.add_argument("--therm-dt", type=float, default=None, help="time step (fs) for thermalization/quench "
                   "(default 0.5 for Ar, 0.25 for N)")
    p.add_argument("--dump-fs", type=float, default=1.0, help="trajectory output interval (fs)")
    p.add_argument("--height", type=float, default=6.0, help="initial projectile height above the sheet (A)")
    p.add_argument("--a0", type=float, default=None,
                   help="rescale the sheet to this lattice constant (A). Default: 2.419 (AIREBO) for Ar, none for N")
    p.add_argument("--no-min", action="store_true", help="skip the initial and final minimizations")
    p.add_argument("--airebo", default="CH.airebo", help="path to the AIREBO potential file")
    p.add_argument("--reaxff", default="ffield.reax.CHN", help="path to the ReaxFF force field")
    p.add_argument("--parallel", type=int, default=1, help="run N impacts at once, 1 core each")
    p.add_argument("--cores", type=int, default=int(os.environ.get("LMP_CORES", 1)),
                   help="MPI ranks per impact (ignored if --parallel > 1)")
    p.add_argument("--mpi-flags", default=os.environ.get("LMP_MPI_FLAGS", "--map-by :OVERSUBSCRIBE"))
    p.add_argument("--lmp", default=os.environ.get("LMP_CMD", "lmp"), help="LAMMPS executable")
    p.add_argument("--workdir", default="bombardment_runs_v2", help="folder for LAMMPS input/log/dump files")
    p.add_argument("--outdir", default="results_v2", help="folder for the final .traj / .xyz files")
    p.add_argument("--force", action="store_true", help="re-run impacts whose final .xyz already exists")
    p.add_argument("--dry-run", action="store_true", help="write the LAMMPS files but do not run LAMMPS")
    return p.parse_args()


def impact_point(atoms, mode, rng):
    """xy position of the impact on the sheet."""
    cell = np.array(atoms.get_cell())
    if mode == "random":
        u, v = rng.random(2)
        return u * cell[0, :2] + v * cell[1, :2]
    center = (cell[0, :2] + cell[1, :2]) / 2.0
    d = np.linalg.norm(atoms.positions[:, :2] - center, axis=1)
    idx = int(np.argmin(d))
    atom_xy = atoms.positions[idx, :2]
    if mode == "top":
        return atom_xy
    d[idx] = np.inf
    return (atom_xy + atoms.positions[int(np.argmin(d)), :2]) / 2.0  # midpoint of the nearest C-C bond


def model_settings(proj, airebo, reaxff):
    """Return (LAMMPS units, pair-style lines, extra fix lines)."""
    if proj == "Ar":
        z_proj = atomic_numbers["Ar"]
        pair = [
            "pair_style hybrid airebo 3.0 zbl 4.0 5.0",
            f"pair_coeff * * airebo {airebo} C NULL",
            f"pair_coeff 2 2 zbl {z_proj} {z_proj}",
            f"pair_coeff 1 2 zbl {atomic_numbers['C']}.0 {z_proj}",
        ]
        return "metal", pair, []
    pair = ["pair_style reaxff NULL", f"pair_coeff * * {reaxff} C N"]
    return "real", pair, ["fix qeq all qeq/reaxff 1 0.0 10.0 1.0e-6 reaxff"]


def projectile_velocity(proj, e_kin, lmp_units):
    """Initial z-velocity (downwards) in LAMMPS units (A/ps for metal, A/fs for real)."""
    mass = atomic_masses[atomic_numbers[proj]]
    v_per_fs = np.sqrt(2.0 * e_kin / mass) * units.fs
    return -v_per_fs * (1000.0 if lmp_units == "metal" else 1.0)


def write_lammps_input(path, proj, e_kin, projectile_id, args, seeds):
    lmp_units, pair, fixes = model_settings(proj, os.path.abspath(args.airebo), os.path.abspath(args.reaxff))
    metal = lmp_units == "metal"
    fs = 1.0e-3 if metal else 1.0                       # fs expressed in LAMMPS time units
    dt = args.timestep * fs
    dt_th_fs = args.therm_dt if args.therm_dt else (0.5 if proj == "Ar" else 0.25)
    dt_th = dt_th_fs * fs
    n_therm = int(round(args.therm_ps * 1000.0 / dt_th_fs))
    n_quench = int(round(args.quench_ps * 1000.0 / dt_th_fs))
    n_impact = int(round(args.time_fs / args.timestep))
    dump_every = max(1, int(round(args.dump_fs / args.timestep)))
    tdamp = 100.0 * fs                                  # 100 fs thermostat damping
    vz = projectile_velocity(proj, e_kin, lmp_units)
    T = args.temperature
    lines = [
        f"# {proj} projectile, {e_kin:g} eV, T = {T:g} K; thermalize -> impact (NVE) -> quench -> minimize",
        f"units {lmp_units}",
        "boundary p p m",
        "atom_style charge",
        "read_data system.data",
        f"mass 1 {atomic_masses[atomic_numbers['C']]:.4f}",
        f"mass 2 {atomic_masses[atomic_numbers[proj]]:.4f}",
        *pair,
        "neighbor 2.0 bin",
        "neigh_modify every 1 delay 0 check yes",
        *fixes,
        f"group proj id {projectile_id}",
        "group sheet subtract all proj",
        "compute tsheet sheet temp",
        "thermo_style custom step time c_tsheet pe ke etotal",
        "thermo 1000",
    ]
    if not args.no_min:
        lines += ["# 0. relax the support", "min_style cg", "minimize 1.0e-10 1.0e-12 2000 20000"]
    lines += [
        "# 1. thermalization (projectile is not integrated: it stays at its starting height)",
        "reset_timestep 0",
        f"timestep {dt_th:.8g}",
        f"velocity sheet create {T} {seeds['velocity']} mom yes rot no dist gaussian",
        f"fix therm sheet nvt temp {T} {T} {tdamp:.8g}",
        f"run {n_therm}",
        "unfix therm",
        "# 2. impact, microcanonical",
        f"velocity proj set 0.0 0.0 {vz:.8f} units box",
        f"timestep {dt:.8g}",
        "fix integ all nve",
        f"dump traj all custom {dump_every} traj.dump id type x y z",
        "dump_modify traj sort id",
        f"run {n_impact}",
        "undump traj",
        "unfix integ",
        "# 3. quench the sheet and relax",
        f"timestep {dt_th:.8g}",
        f"fix quench sheet langevin {T} {T} {tdamp:.8g} {seeds['langevin']}",
        "fix qnve sheet nve",
        f"run {n_quench}",
        "unfix quench",
        "unfix qnve",
    ]
    if not args.no_min:
        lines += ["min_style cg", "minimize 1.0e-8 1.0e-10 3000 30000"]
    lines += ["write_dump all custom final.dump id type x y z modify sort id"]
    Path(path).write_text("\n".join(lines) + "\n")
    return vz, lmp_units, dict(n_therm=n_therm, n_impact=n_impact, n_quench=n_quench, dt_th_fs=dt_th_fs)


def lammps_command(args, input_name, log_name):
    cmd = []
    if args.cores > 1 and args.parallel <= 1:
        cmd += ["mpirun", *shlex.split(args.mpi_flags), "-np", str(args.cores)]
    cmd += shlex.split(args.lmp) + ["-in", input_name, "-log", log_name, "-screen", "none"]
    return cmd


def check_files(args):
    needed = [args.structure]
    if "Ar" in args.projectiles:
        needed.append(args.airebo)
    if "N" in args.projectiles:
        needed.append(args.reaxff)
    missing = [f for f in needed if not Path(f).is_file()]
    if missing:
        sys.exit("Missing required file(s): " + ", ".join(missing))


def rescaled_sheet(base, a0_target, a0_current):
    """Isotropic in-plane rescaling to the target lattice constant (z untouched)."""
    atoms = base.copy()
    s = a0_target / a0_current
    cell = np.array(atoms.get_cell())
    cell[:2] *= s
    atoms.set_cell(cell, scale_atoms=True)
    return atoms


def run_job(job, args, base, a0_current, workdir, outdir):
    """Build, run and convert one impact. Returns (job, status, wall_time_s, message)."""
    proj, e_kin, mode, run_id, tag = job
    final_xyz = outdir / f"graphene_defective_{tag}.xyz"
    traj_out = outdir / f"bombardment_{tag}.traj"

    # reproducible seeds per (projectile, energy, target, run)
    ss = np.random.SeedSequence([args.seed, atomic_numbers[proj], int(round(e_kin * 10)),
                                 TARGET_MODES.index(mode), run_id])
    rng = np.random.default_rng(ss)
    seeds = dict(velocity=int(rng.integers(1, 800_000_000)), langevin=int(rng.integers(1, 800_000_000)))

    a0_target = args.a0 if args.a0 else (A0_AIREBO if proj == "Ar" else None)
    sheet = rescaled_sheet(base, a0_target, a0_current) if a0_target else base.copy()
    z_sheet = float(np.mean(sheet.positions[:, 2]))
    xy = impact_point(sheet, mode, rng)
    atoms = sheet.copy()
    atoms.append(Atom(proj, position=(xy[0], xy[1], z_sheet + args.height)))
    atoms.set_initial_charges(np.zeros(len(atoms)))
    projectile_id = len(atoms)

    rundir = workdir / tag
    rundir.mkdir(parents=True, exist_ok=True)
    write(rundir / "system.data", atoms, format="lammps-data", specorder=["C", proj], atom_style="charge")
    write(rundir / "initial.xyz", sheet)  # reference lattice for the analysis
    vz, lmp_units, steps = write_lammps_input(rundir / "in.bombard", proj, e_kin, projectile_id, args, seeds)
    d_c = np.linalg.norm(sheet.positions[:, :2] - xy, axis=1)
    meta = dict(projectile=proj, energy_eV=e_kin, target=mode, run=run_id, impact_xy=[float(xy[0]), float(xy[1])],
                a0_used=a0_target or a0_current, temperature_K=args.temperature, timestep_fs=args.timestep,
                seeds=seeds, **steps)
    (rundir / "meta.json").write_text(json.dumps(meta, indent=2))
    cmd = lammps_command(args, "in.bombard", "log.lammps")
    if args.dry_run:
        return job, "dry-run", 0, f"impact at ({xy[0]:.2f}, {xy[1]:.2f}), v_z = {vz:.5f}; files in {rundir}"

    t0 = time.time()
    proc = subprocess.run(cmd, cwd=rundir, capture_output=True, text=True)
    wall = time.time() - t0
    traj_dump, final_dump = rundir / "traj.dump", rundir / "final.dump"
    if proc.returncode != 0 or not traj_dump.exists() or not final_dump.exists():
        log = rundir / "log.lammps"
        tail = log.read_text().splitlines()[-15:] if log.exists() else []
        return job, "failed", wall, f"LAMMPS failed (exit {proc.returncode}):\n" + ("\n".join(tail) or proc.stderr[-1500:])

    frames = read(traj_dump, format="lammps-dump-text", index=":", specorder=["C", proj])
    write(traj_out, frames)
    final = read(final_dump, format="lammps-dump-text", index=-1, specorder=["C", proj])
    write(final_xyz, final)  # written last: its existence marks a finished run
    return job, "ok", wall, f"impact ({xy[0]:.1f},{xy[1]:.1f}) -> {final_xyz.name} [{wall / 60:.1f} min]"


def main():
    args = parse_args()
    if not args.dry_run:
        check_files(args)
    elif not Path(args.structure).is_file():
        sys.exit(f"Missing structure file: {args.structure}")
    if "N" in args.projectiles:
        print("WARNING: N impacts need a ReaxFF set validated for graphene (see check_potentials.py).", flush=True)

    workdir, outdir = Path(args.workdir), Path(args.outdir)
    workdir.mkdir(parents=True, exist_ok=True)
    outdir.mkdir(parents=True, exist_ok=True)

    base = read(args.structure)
    a0_current = float(np.linalg.norm(np.array(base.get_cell())[0, :2])) / int(round(
        np.linalg.norm(np.array(base.get_cell())[0, :2]) / 2.46))  # lattice constant of the input sheet
    jobs = [(proj, e_kin, mode, run_id, f"{proj}_{mode}_E{e_kin:g}_run{run_id}")
            for proj in args.projectiles for e_kin in args.energies
            for mode in args.targets for run_id in range(1, (args.impacts if mode == "random" else 1) + 1)]
    print(f"--- {len(jobs)} impacts, T = {args.temperature:g} K, dt = {args.timestep} fs, impact stage "
          f"{args.time_fs:g} fs; {max(args.parallel, 1)} at a time ---", flush=True)

    summary_path = outdir / "runs_summary.csv"
    new_summary = (not summary_path.exists()) or summary_path.stat().st_size == 0
    failures, todo = [], []
    with open(summary_path, "a", newline="") as fh:
        writer = csv.writer(fh)
        if new_summary:
            writer.writerow(["projectile", "energy_eV", "target", "run", "status", "wall_time_s"])
        for job in jobs:
            proj, e_kin, mode, run_id, tag = job
            if (outdir / f"graphene_defective_{tag}.xyz").exists() and not args.force and not args.dry_run:
                writer.writerow([proj, e_kin, mode, run_id, "skipped", 0])
            else:
                todo.append(job)
        fh.flush()
        print(f"{len(jobs) - len(todo)} already done (skipped), {len(todo)} to run", flush=True)

        done = 0
        with ThreadPoolExecutor(max_workers=max(args.parallel, 1)) as pool:
            futures = {pool.submit(run_job, job, args, base, a0_current, workdir, outdir): job for job in todo}
            for fut in as_completed(futures):
                proj, e_kin, mode, run_id, tag = futures[fut]
                done += 1
                try:
                    _, status, wall, msg = fut.result()
                except Exception as exc:
                    status, wall, msg = "failed", 0, f"{type(exc).__name__}: {exc}"
                if status == "failed":
                    failures.append(tag)
                print(f"[{done}/{len(todo)}] {tag}: {status} -> {msg}", flush=True)
                writer.writerow([proj, e_kin, mode, run_id, status, f"{wall:.0f}"])
                fh.flush()

    if failures:
        sys.exit(f"\n{len(failures)} impact(s) failed: " + ", ".join(failures))
    print("\nAll impacts finished." if not args.dry_run else "\nDry run finished.")


if __name__ == "__main__":
    main()