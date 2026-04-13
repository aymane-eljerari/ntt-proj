#!/bin/bash
set -uo pipefail

# ================= CONFIGURATION =================
N_VALUES=(16 18)
L_VALUES=(24 32)
BARRETT_FLAGS=("" "--barrett")
ALGORITHMS=("--gpu-radix2" "--gpu-radix4")

TARGET_APP="./build/bench_ntt"

OUT_DIR="ncu_reps"
mkdir -p "$OUT_DIR"

# ================= NCU SETTINGS =================
NCU_SECTIONS="--section SpeedOfLight \
              --section MemoryWorkloadAnalysis \
              --section ComputeWorkloadAnalysis \
              --section Occupancy \
              --section SchedulerStats"

# ─── Preflight check ────────────────────────────────────────────────────────
if [[ ! -x "$TARGET_APP" ]]; then
    echo "ERROR: $TARGET_APP not found or not executable." >&2
    exit 1
fi

if ! command -v ncu &> /dev/null; then
    echo "ERROR: ncu (Nsight Compute) not found in PATH." >&2
    exit 1
fi

# ================= EXECUTION LOOP =================
total=$(( ${#ALGORITHMS[@]} * ${#N_VALUES[@]} * ${#L_VALUES[@]} * ${#BARRETT_FLAGS[@]} ))
current=0

echo "Starting NTT GPU Parameter Sweep Profiling..."
echo "Configurations: $total"
echo "Reports will be saved to ./${OUT_DIR}/"
echo "---------------------------------------------------"

for algo in "${ALGORITHMS[@]}"; do
    algo_name=$(echo "$algo" | sed 's/--//')
    for n in "${N_VALUES[@]}"; do
        for l in "${L_VALUES[@]}"; do
            for b_flag in "${BARRETT_FLAGS[@]}"; do
                current=$((current + 1))

                b_name="NoBarrett"
                if [ -n "$b_flag" ]; then
                    b_name="Barrett"
                fi

                FILE_BASE="${OUT_DIR}/${algo_name}_N${n}_L${l}_${b_name}"

                echo "[$current/$total] ${algo_name}  N=2^${n}  L=${l}  ${b_name}"

                CMD="$TARGET_APP -N $n -L $l -T 1 $algo $b_flag"

                # 1. Generate .ncu-rep for GUI inspection
                if ! ncu -o "$FILE_BASE" -f $NCU_SECTIONS $CMD > /dev/null 2>&1; then
                    echo "  WARNING: ncu-rep generation failed" >&2
                fi

                # 2. Extract flat CSV for aggregation
                if ! ncu --csv --page details $NCU_SECTIONS $CMD > "${FILE_BASE}.csv" 2>/dev/null; then
                    echo "  WARNING: CSV export failed" >&2
                fi

            done
        done
    done
done

echo "---------------------------------------------------"
echo "Sweep complete. $total configurations profiled."
echo "  .ncu-rep files: open in Nsight Compute UI"
echo "  .csv files:     aggregate with Python"