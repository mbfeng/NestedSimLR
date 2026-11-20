#!/bin/sh
# Full Phase-2 reruns (November 2026). Can be started from any folder. Parts
# whose output already exists are skipped, so the script can be restarted
# after an interruption. Progress: logs/run_all_reruns.log
#
# Environment variables: RSCRIPT (path to Rscript), CORES (parallel workers,
# default 12), MEM_BUDGET_GB (Figure 3(b) memory cap, default 24).
# Example on a Linux server, detached from the terminal:
#   RSCRIPT=Rscript CORES=24 MEM_BUDGET_GB=150 setsid -f nice -n 10 \
#     sh tools/run_all_reruns.sh >> logs/run_all_reruns.log 2>&1 < /dev/null
cd "$(dirname "$0")/.." || exit 1
# One BLAS/OpenMP thread per R process: parallelism comes from the workers.
# (A multithreaded OpenBLAS, as on Linux servers, otherwise starts one thread
# per core in every worker and oversubscribes the machine.)
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
R="${RSCRIPT:-/c/Program Files/R/R-4.3.3/bin/Rscript.exe}"
CORES="${CORES:-12}"
OUT=results/rerun_full
BIG=results/rerun_10000_largeM
TRUTH="$OUT/sec71_table1_truth.csv"
mkdir -p "$OUT" "$BIG" logs
run() { log=$1; shift; echo "== $(date) start $*"
        "$R" "$@" > "logs/$log" 2>&1 && echo "== $(date) done" ||
        echo "== $(date) FAILED (see logs/$log)"; }

[ -f "$OUT/sec73_table4_per_macro.csv.gz" ] ||
  run sec73_table4.log scripts/sec73_simulate.R table4 "$CORES" "$OUT"
[ -f "$OUT/sec73_fig3a_var_ratio.csv" ] ||
  run sec73_fig3a.log scripts/sec73_simulate.R fig3a "$CORES" "$OUT"
[ -f "$OUT/sec71_fig2_bands.csv" ] ||
  run sec71_fig2.log scripts/sec71_simulate.R fig2 10000 "$CORES" "$OUT"
[ -f "$TRUTH" ] ||
  run sec71_truth.log scripts/sec71_simulate.R truth 1 1 "$OUT"
[ -f "$OUT/sec72_per_macro.csv" ] ||
  run sec72_full.log scripts/sec72_simulate.R 1000 "$CORES" "$OUT" full
[ -f "$BIG/sec71_table1_per_macro.csv.gz" ] ||
  run sec71_table1_largeM.log scripts/sec71_simulate.R table1 10000 "$CORES" "$BIG" 2048,4096
[ -f "$TRUTH" ] && cp "$TRUTH" "$BIG/"
# Figure 3(b) resumes by itself (skips M values already saved)
run sec73_fig3b.log scripts/sec73_simulate.R fig3b "$CORES" "$OUT"
echo "== $(date) all done"
