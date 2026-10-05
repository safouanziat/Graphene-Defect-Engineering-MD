#!/usr/bin/env bash
# Overnight convergence + statistics campaign. Run from the Molecular_Dynamics folder:
#
#   nohup bash run_overnight.sh > overnight.log 2>&1 &          # survives closing the terminal
#   DEADLINE=07:30 PAR=10 nohup bash run_overnight.sh > overnight.log 2>&1 &
#
# Steps run in order of scientific value, so whatever is finished by the morning is useful. Every MD step is
# resumable (finished runs are skipped): if the machine reboots, start the script again.
#   PAR           single-core LAMMPS jobs at once (default 10)
#   JOBS          analysis workers (default 7)
#   DEADLINE      HH:MM; no new step starts after this time (a running step is NOT killed)
#   MIN_FREE_GB   a step is skipped if less disk is free (default 20; each run writes ~60 MB of trajectory,
#                 ~140 MB on the 4608-atom sheet)
# Stop cleanly at any time with:  touch STOP_OVERNIGHT   (the current step finishes its running jobs first)
set -uo pipefail
export PAR="${PAR:-10}" JOBS="${JOBS:-7}"
MIN_FREE_GB="${MIN_FREE_GB:-20}"
mkdir -p overnight_logs
MD="python 02_defect_engineering_stats.py --projectiles Ar --parallel $PAR"
AN="python 03_analyze_impacts.py --jobs $JOBS"
say() { echo "[$(date '+%a %H:%M:%S')] $*"; }

too_late() {
  [[ -f STOP_OVERNIGHT ]] && return 0
  if [[ -n "${DEADLINE:-}" ]]; then
    local now dl; now=$(date +%s); dl=$(date -d "$DEADLINE" +%s 2>/dev/null || echo 0)
    (( dl < now - 43200 )) && dl=$(( dl + 86400 ))      # deadline given as tomorrow's clock time
    (( now >= dl )) && return 0
  fi
  return 1
}
enough_disk() {
  local free; free=$(df -Pk . | awk 'NR==2 {printf "%d", $4/1024/1024}')
  (( free >= MIN_FREE_GB )) || { say "only ${free} GB free (< $MIN_FREE_GB): skipping"; return 1; }
}
step() {                    # step <name> <command...>
  local name="$1"; shift
  if too_late; then say "SKIP $name (deadline or STOP file)"; return 0; fi
  enough_disk || return 0
  say "START $name"; local t0=$SECONDS
  "$@" > "overnight_logs/$name.log" 2>&1; local rc=$?
  say "END   $name  rc=$rc  $(( (SECONDS - t0) / 60 )) min"
}
analyse() { $AN --indir "$1/results" --workdir "$1/runs" ${3:+--structure $3} --out "$1/impact_analysis.csv" --summary "$1/yield_summary.csv"; }

q_quench()  { $MD --energies 70 100 --impacts 12 --quench-ps 5 --workdir conv_quench5/runs --outdir conv_quench5/results \
              && analyse conv_quench5 && python 11_compare_runs.py impact_analysis.csv conv_quench5/impact_analysis.csv --label "quench 5 ps" --paired; }
q_window()  { $MD --energies 70 100 --impacts 12 --time-fs 1500 --workdir conv_window1500/runs --outdir conv_window1500/results \
              && analyse conv_window1500 && python 11_compare_runs.py impact_analysis.csv conv_window1500/impact_analysis.csv --label "window 1500 fs" --paired; }
q_size()    { python -c "from ase.build import graphene; from ase.io import write; write('graphene_48x48_relaxed.xyz', graphene(a=2.46, size=(48,48,1), vacuum=15.0))" \
              && $MD --structure graphene_48x48_relaxed.xyz --energies 70 100 --impacts 24 --workdir conv_big/runs --outdir conv_big/results \
              && analyse conv_big x graphene_48x48_relaxed.xyz \
              && python 11_compare_runs.py impact_analysis.csv conv_big/impact_analysis.csv --label "4608 atoms" \
              && python 10_sheet_heating.py --csv conv_big/impact_analysis.csv --workdir conv_big/runs --n-sheet 4608 --out conv_big; }
q_stats()   { $MD --energies 70 100 --impacts "$1" && $AN; }
q_onset()   { $MD --energies 62 64 66 68 --impacts "$1" && $AN; }
q_report()  { $AN && python 04_plot_yield.py && python 07_plot_physics.py && python 08_energy_loss.py && python 10_sheet_heating.py; }
q_cutoff()  { mkdir -p conv_cutoff
              for ej in 2.5 3.5; do $AN --eject $ej --out conv_cutoff/eject_$ej.csv --summary conv_cutoff/s_eject_$ej.csv \
                  && python 11_compare_runs.py impact_analysis.csv conv_cutoff/eject_$ej.csv --label "eject $ej A"; done
              for cc in 1.75 1.95; do $AN --cc-cut $cc --out conv_cutoff/cc_$cc.csv --summary conv_cutoff/s_cc_$cc.csv \
                  && python 11_compare_runs.py impact_analysis.csv conv_cutoff/cc_$cc.csv --label "bond cutoff $cc A"; done; }

say "campaign start (PAR=$PAR, JOBS=$JOBS, DEADLINE=${DEADLINE:-none})"
q_heat() { python 08_energy_loss.py && python 10_sheet_heating.py; }
step heating_existing   q_heat
step 1_quench_5ps       q_quench
step 2_window_1500fs    q_window
step 3_size_4608        q_size
step 4_stats_N48        q_stats 48
step 5_onset_N24        q_onset 24
step 6_stats_N72        q_stats 72
step 7_onset_N48        q_onset 48
step 8_stats_N100       q_stats 100
step 9_cutoffs          q_cutoff
step 10_final_report    q_report
say "campaign end"
echo; echo "=== verdicts ==="; grep -H "Verdict\|paired runs" overnight_logs/*.log 2>/dev/null