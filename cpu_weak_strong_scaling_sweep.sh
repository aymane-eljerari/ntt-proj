#!/bin/bash

# Configuration
EXEC="./build/bench_ntt"
CSV_STRONG="strong_scaling_results.csv"
CSV_WEAK1="weak_scaling_limb_sweep_results.csv"
CSV_WEAK2="weak_scaling_poly_sweep_results.csv"

# Thread sweep parameters
THREADS=(1 2 4 8 16 32 64)

# Verify executable exists
if [ ! -f "$EXEC" ]; then
    echo "Error: Executable not found at $EXEC. Please build the project first."
    exit 1
fi

echo "======================================================="
echo "Starting NTT Benchmark Parameter Sweeps"
echo "======================================================="

# ---------------------------------------------------------
# 1. Strong Scaling: Fixed N=2^18, Fixed L=64
# ---------------------------------------------------------
echo ""
echo ">>> Running Strong Scaling (Fixed N=18, Fixed L=64) <<<"
echo "Results will be appended to: $CSV_STRONG"
N=18
L=64

for T in "${THREADS[@]}"; do
    echo "Executing -> N=$N, L=$L, Threads=$T"
    $EXEC -N $N -L $L -T $T --cpu-prod --openfhe --csv-out "$CSV_STRONG" --barrett
done

# ---------------------------------------------------------
# 2. Weak Scaling 1: Fixed N=2^18, Varying L
# Sweeping through L values 32, 48, 64, 72, 96, 128
# ---------------------------------------------------------
echo ""
echo ">>> Running Weak Scaling 1 (Fixed N=18, Varying L) <<<"
echo "Results will be appended to: $CSV_WEAK1"
LIMB_SWEEP=(32 48 64 72 96 128)
N=18

for L in "${LIMB_SWEEP[@]}"; do
    for T in "${THREADS[@]}"; do
        echo "Executing -> N=$N, L=$L, Threads=$T"
        $EXEC -N $N -L $L -T $T --cpu-prod --openfhe --csv-out "$CSV_WEAK1" --barrett
    done
done

# ---------------------------------------------------------
# 3. Weak Scaling 2: Varying N (2^12 to 2^18), Fixed L=64
# ---------------------------------------------------------
echo ""
echo ">>> Running Weak Scaling 2 (Varying N, Fixed L=64) <<<"
echo "Results will be appended to: $CSV_WEAK2"
N_SWEEP=(12 13 14 15 16 17 18)
L=64

for N in "${N_SWEEP[@]}"; do
    for T in "${THREADS[@]}"; do
        echo "Executing -> N=$N, L=$L, Threads=$T"
        $EXEC -N $N -L $L -T $T --cpu-prod --openfhe --csv-out "$CSV_WEAK2" --barrett
    done
done

echo ""
echo "======================================================="
echo "All benchmarks completed successfully."
echo "======================================================="