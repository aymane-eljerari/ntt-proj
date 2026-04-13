#!/usr/bin/env python3
"""
Aggregate per-stage NCU kernel metrics into an end-to-end (I)NTT view.

Reads all CSVs from the ncu_reps/ directory, groups kernel invocations
by configuration (algo, N, L, barrett), and produces a single summary CSV
where each row represents the full NTT pipeline for that configuration.

Aggregation rules:
  - Durations, cycles, instruction counts  → SUM across stages
  - Throughput %, hit rates, occupancy, IPC → weighted average by Duration
  - Frequencies                             → plain average (should be constant)
"""

import os
import re
import sys
import glob
import pandas as pd
import numpy as np

# ─── Config ──────────────────────────────────────────────────────────────────
NCU_DIR = "ncu_reps"
OUTPUT_CSV = "ncu_aggregated.csv"

# ─── Filters (set to None to include all) ────────────────────────────────────
FILTER_N = 18
FILTER_L = 32
# FILTER_ALGO = None      # e.g. "gpu-radix2" or None for all
# FILTER_BARRETT = None   # e.g. "Barrett" or None for all

# The 5 kernels to track
KERNEL_PATTERNS = {
    "dif_gs_ntt_stage":      "ntt",
    "dit_ct_intt_stage":     "intt",
    "normalize_intt":        "intt",
    "dif_radix4_ntt_stage":  "ntt",
    "dit_radix4_ntt_stage":  "intt",
}

# Aggregation rules by metric name substring
# "sum"     → add across all kernel stages
# "wavg"    → duration-weighted average
# "avg"     → simple average
SUM_KEYWORDS = [
    "duration", "elapsed cycles", "sm active cycles",
    "memory throughput",  # the byte/s one gets summed then we derive effective BW
]

AVG_KEYWORDS = [
    "frequency",
]

# Everything else (percentages, rates, IPC) defaults to weighted average


# ─── Helpers ─────────────────────────────────────────────────────────────────
def parse_config_from_filename(filename: str) -> dict:
    """Extract algo, N, L, barrett from filename like gpu-radix2_N14_L24_Barrett.csv"""
    base = os.path.basename(filename).replace(".csv", "")
    parts = base.split("_")

    algo = parts[0]  # e.g. "gpu-radix2"
    n_val = None
    l_val = None
    barrett = "NoBarrett"

    for p in parts[1:]:
        if p.lower() in ("barrett",):
            barrett = "Barrett"
        elif p.lower() in ("nobarrett",):
            barrett = "NoBarrett"
        elif p.startswith("N") and p[1:].isdigit():
            n_val = int(p[1:])
        elif p.startswith("L") and p[1:].isdigit():
            l_val = int(p[1:])

    return {"algo": algo, "N": n_val, "L": l_val, "barrett": barrett}


def classify_agg_rule(metric_name: str, metric_unit: str) -> str:
    """Decide how to aggregate a metric: sum, wavg, or avg."""
    name_lower = metric_name.lower()
    unit_lower = metric_unit.lower()

    for kw in SUM_KEYWORDS:
        if kw in name_lower:
            return "sum"

    for kw in AVG_KEYWORDS:
        if kw in name_lower:
            return "avg"

    # Percentages and rates → weighted average
    if unit_lower == "%" or "rate" in name_lower or "ipc" in name_lower:
        return "wavg"

    # Counts (cycles, instructions, bytes) → sum
    if unit_lower in ("cycle", "inst", "byte", "sector", "request", "warp"):
        return "sum"

    # Default to weighted average
    return "wavg"


def is_ntt_kernel(kernel_name: str) -> bool:
    """Check if this kernel is one of our 5 NTT kernels."""
    return any(pat in kernel_name for pat in KERNEL_PATTERNS)


def load_ncu_csv(filepath: str) -> pd.DataFrame:
    """Load an NCU CSV, skipping non-data preamble lines and rule rows."""
    # NCU CSVs may have preamble text before the header row
    with open(filepath, "r") as f:
        lines = f.readlines()

    # Find the header line (starts with "ID")
    header_idx = None
    for i, line in enumerate(lines):
        if line.strip().startswith('"ID"'):
            header_idx = i
            break

    if header_idx is None:
        return pd.DataFrame()

    from io import StringIO
    csv_text = "".join(lines[header_idx:])
    df = pd.read_csv(StringIO(csv_text), quotechar='"')

    # Drop rule/recommendation rows (they have empty Metric Name)
    df = df[df["Metric Name"].notna() & (df["Metric Name"].str.strip() != "")]

    # Keep only our NTT kernels
    df = df[df["Kernel Name"].apply(is_ntt_kernel)].copy()

    # Clean up metric values
    df["Metric Value"] = pd.to_numeric(
        df["Metric Value"].astype(str).str.replace(",", ""), errors="coerce"
    )

    return df


def aggregate_config(df: pd.DataFrame) -> pd.DataFrame:
    """
    Aggregate per-kernel-invocation metrics into a single end-to-end row.

    Each unique (Section Name, Metric Name, Metric Unit) becomes one output column.
    """
    if df.empty:
        return pd.DataFrame()

    # Get per-invocation durations for weighting
    durations = (
        df[df["Metric Name"] == "Duration"]
        .groupby("Kernel Name")["Metric Value"]
        .first()
    )

    # Build per-invocation duration lookup: map each row's kernel to its duration
    df = df.merge(
        durations.rename("_kernel_duration"),
        left_on="Kernel Name",
        right_index=True,
        how="left",
    )
    df["_kernel_duration"] = df["_kernel_duration"].fillna(1)

    results = {}
    total_duration_ns = durations.sum()
    results[("Summary", "Total Duration", "ns")] = total_duration_ns
    results[("Summary", "Num Kernel Launches", "count")] = len(
        df[["Kernel Name", "ID"]].drop_duplicates()
    )

    # Group by metric identity
    for (section, metric, unit), group in df.groupby(
        ["Section Name", "Metric Name", "Metric Unit"]
    ):
        values = group["Metric Value"].dropna()
        weights = group.loc[values.index, "_kernel_duration"]

        if values.empty:
            continue

        rule = classify_agg_rule(metric, unit)

        if rule == "sum":
            results[(section, metric, unit)] = values.sum()
        elif rule == "avg":
            results[(section, metric, unit)] = values.mean()
        else:  # wavg
            if weights.sum() > 0:
                results[(section, metric, unit)] = np.average(values, weights=weights)
            else:
                results[(section, metric, unit)] = values.mean()

    return results


# ─── Main ────────────────────────────────────────────────────────────────────
def main():
    csv_files = sorted(glob.glob(os.path.join(NCU_DIR, "*.csv")))

    if not csv_files:
        print(f"No CSV files found in {NCU_DIR}/")
        sys.exit(1)

    # Apply filters
    def matches_filter(filepath):
        cfg = parse_config_from_filename(filepath)
        if FILTER_N is not None and cfg["N"] != FILTER_N:
            return False
        if FILTER_L is not None and cfg["L"] != FILTER_L:
            return False
        # if FILTER_ALGO is not None and cfg["algo"] != FILTER_ALGO:
        #     return False
        # if FILTER_BARRETT is not None and cfg["barrett"] != FILTER_BARRETT:
        #     return False
        return True

    csv_files = [f for f in csv_files if matches_filter(f)]

    print(f"Found {len(csv_files)} NCU CSV files (after filtering: N={FILTER_N}, L={FILTER_L})")

    all_rows = []

    for filepath in csv_files:
        config = parse_config_from_filename(filepath)
        print(f"  Processing: {os.path.basename(filepath)} -> {config}")

        df = load_ncu_csv(filepath)
        if df.empty:
            print(f"    WARNING: no kernel data found, skipping")
            continue

        agg = aggregate_config(df)
        if not agg:
            continue

        # Flatten: config columns + one column per (section, metric, unit)
        row = dict(config)
        for (section, metric, unit), value in agg.items():
            col_name = f"{metric} ({unit})"
            row[col_name] = value

        all_rows.append(row)

    if not all_rows:
        print("No data aggregated.")
        sys.exit(1)

    result_df = pd.DataFrame(all_rows)

    # Sort for readability
    result_df = result_df.sort_values(["algo", "barrett", "N", "L"]).reset_index(
        drop=True
    )

    # Move config columns to front
    config_cols = ["algo", "N", "L", "barrett"]
    other_cols = [c for c in result_df.columns if c not in config_cols]
    result_df = result_df[config_cols + sorted(other_cols)]

    result_df.to_csv(OUTPUT_CSV, index=False, float_format="%.4f")
    print(f"\nAggregated results written to {OUTPUT_CSV}")
    print(f"  Configurations: {len(result_df)}")
    print(f"  Metrics:        {len(other_cols)}")

    # Print a quick summary
    print(f"\n{'algo':<16} {'N':>4} {'L':>4} {'barrett':<10} {'Total Duration (ns)':>20}")
    print("-" * 60)
    for _, r in result_df.iterrows():
        dur = r.get("Total Duration (ns)", "N/A")
        if isinstance(dur, float):
            dur = f"{dur:>20,.0f}"
        print(f"{r['algo']:<16} {r['N']:>4} {r['L']:>4} {r['barrett']:<10} {dur}")


if __name__ == "__main__":
    main()