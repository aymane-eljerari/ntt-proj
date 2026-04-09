#!/usr/bin/env bash
set -euo pipefail

# ─── Configuration ───────────────────────────────────────────────────────────
BENCH_BIN="./build/bench_ntt"
OUTPUT_CSV="bench_results3.csv"
THREADS=32

N_EXPONENTS=(10 11 12 13 14 15 16 17 18 19 20)
LIMBS=(24 32 48 64 128 256)

declare -A IMPL_FLAGS=(
    ["cpu_fast"]="--cpu-fast"
    ["cpu_prod"]="--cpu-prod"
    ["omp_openfhe"]="--openfhe"
    ["gpu_r2"]="--gpu-radix2"
    ["gpu_r4"]="--gpu-radix4"
)

# Sweep definitions: "suffix|extra_flags"
SWEEPS=(
    "|"                  # Sweep 1: no extra flags
    "_barrett|--barrett" # Sweep 2: Barrett reduction
)

# ─── Preflight check ────────────────────────────────────────────────────────
if [[ ! -x "$BENCH_BIN" ]]; then
    echo "ERROR: $BENCH_BIN not found or not executable." >&2
    exit 1
fi

# ─── Main ────────────────────────────────────────────────────────────────────
echo "Implementation,N,L,Result" > "$OUTPUT_CSV"

configs_per_sweep=$(( ${#IMPL_FLAGS[@]} * ${#N_EXPONENTS[@]} * ${#LIMBS[@]} ))
total_configs=$(( configs_per_sweep * ${#SWEEPS[@]} ))
current=0

for sweep in "${SWEEPS[@]}"; do
    IFS='|' read -r suffix extra_flags <<< "$sweep"
    echo ""
    echo "════════════════════════════════════════════════════════"
    echo "  Sweep: ${suffix:-default} ${extra_flags:+(flags: $extra_flags)}"
    echo "════════════════════════════════════════════════════════"

    for impl in "${!IMPL_FLAGS[@]}"; do
        flag="${IMPL_FLAGS[$impl]}"
        for n_exp in "${N_EXPONENTS[@]}"; do
            N=$((1 << n_exp))
            for limbs in "${LIMBS[@]}"; do
                current=$((current + 1))
                echo "[$current/$total_configs] ${impl}${suffix}  N=2^${n_exp}=$N  L=$limbs"

                runtime=$("$BENCH_BIN" $flag $extra_flags -N "$n_exp" -L "$limbs" -T "$THREADS" \
                          | grep -oE '[0-9]+(\.[0-9]+)?' | tail -1)

                if [[ -z "$runtime" ]]; then
                    echo "  WARNING: no numeric output, writing N/A" >&2
                    runtime="N/A"
                fi

                echo "${impl}${suffix},$N,$limbs,$runtime" >> "$OUTPUT_CSV"
            done
        done
    done
done

echo ""
echo "Done. Results written to $OUTPUT_CSV"
echo "Total configurations: $total_configs"