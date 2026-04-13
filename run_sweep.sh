#!/usr/bin/env bash
set -uo pipefail

# ─── Configuration ───────────────────────────────────────────────────────────
BENCH_BIN="./build/bench_ntt"
OUTPUT_CSV="bench_results.csv"
THREADS=32

N_EXPONENTS=(10 11 12 13 14 15 16 17 18 19 20)
LIMBS=(24 32 48 64 128 256)

# Ordered: GPU first, then CPU/OpenFHE
IMPL_ORDER=(gpu_r2 gpu_r4 cpu_fast cpu_prod omp_openfhe)
declare -A IMPL_FLAGS=(
    ["gpu_r2"]="--gpu-radix2"
    ["gpu_r4"]="--gpu-radix4"
    ["cpu_fast"]="--cpu-fast"
    ["cpu_prod"]="--cpu-prod"
    ["omp_openfhe"]="--openfhe"
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

configs_per_sweep=$(( ${#IMPL_ORDER[@]} * ${#N_EXPONENTS[@]} * ${#LIMBS[@]} ))
total_configs=$(( configs_per_sweep * ${#SWEEPS[@]} ))
current=0
failures=0

for sweep in "${SWEEPS[@]}"; do
    IFS='|' read -r suffix extra_flags <<< "$sweep"
    echo ""
    echo "════════════════════════════════════════════════════════"
    echo "  Sweep: ${suffix:-default} ${extra_flags:+(flags: $extra_flags)}"
    echo "════════════════════════════════════════════════════════"

    for impl in "${IMPL_ORDER[@]}"; do
        flag="${IMPL_FLAGS[$impl]}"
        for n_exp in "${N_EXPONENTS[@]}"; do
            N=$((1 << n_exp))
            for limbs in "${LIMBS[@]}"; do
                current=$((current + 1))
                echo "[$current/$total_configs] ${impl}${suffix}  N=2^${n_exp}=$N  L=$limbs"

                # Run benchmark, capture output; || true prevents exit on failure
                raw_output=$("$BENCH_BIN" $flag $extra_flags -N "$n_exp" -L "$limbs" -T "$THREADS" 2>&1) || true

                # Extract last number from output
                runtime=$(echo "$raw_output" | grep -oE '[0-9]+(\.[0-9]+)?' | tail -1 || true)

                if [[ -z "$runtime" ]]; then
                    echo "  WARNING: no numeric output, writing N/A" >&2
                    echo "  Output was: $raw_output" >&2
                    runtime="N/A"
                    failures=$((failures + 1))
                fi

                echo "${impl}${suffix},$N,$limbs,$runtime" >> "$OUTPUT_CSV"
            done
        done
    done
done

echo ""
echo "Done. Results written to $OUTPUT_CSV"
echo "Total configurations: $total_configs"
if [[ $failures -gt 0 ]]; then
    echo "WARNING: $failures configuration(s) produced no output (recorded as N/A)"
fi