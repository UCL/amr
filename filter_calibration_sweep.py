"""
filter_calibration_sweep.py

Batch filter for the multi-stage calibration parameter sweep.

For each simulation_summary_{run_id}.csv in the output directory, computes the
infection resistance simulation mean (%) over the 2022-2025 calibration window
and applies the stage threshold.  If a matching sampled_parameters_{run_id}.csv
exists it is joined so the surviving multiplier values are preserved for the next
stage.

Usage
-----
# Stage 1 filter (pop ~10k runs)
python filter_calibration_sweep.py --stage 1 --input-dir amr_simulation_output_analysis_outputs

# Stage 2 filter (pop ~30k runs, after re-running survivors from stage 1)
python filter_calibration_sweep.py --stage 2 --input-dir amr_simulation_output_analysis_outputs

# Stage 3 filter (pop ~100k)
python filter_calibration_sweep.py --stage 3 --input-dir amr_simulation_output_analysis_outputs

# Stage 4 / final filter (pop ~1M)
python filter_calibration_sweep.py --stage 4 --input-dir amr_simulation_output_analysis_outputs

Optionally restrict to a specific set of run IDs (e.g. only run stage-2 on stage-1 survivors):
python filter_calibration_sweep.py --stage 2 --run-ids stage1_survivors.csv

Output
------
Writes a CSV named  sweep_stage{N}_survivors.csv  containing one row per surviving
run with columns:
  run_id, inf_res_mean_pct, [sampled_quantity_1, sampled_quantity_2, ...]

The stage 4 output also includes the overall calibration score if a
calibration_summary_{run_id}.txt file exists alongside the CSV.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Stage thresholds
# ---------------------------------------------------------------------------
# Infection resistance sim mean (%) must fall within [low, high].
# Upper caps are generous to account for small-population extinction bias
# (which suppresses the simulation mean below the true value).
# As population grows across stages the variance shrinks and the window
# tightens symmetrically.  Target mean from a large reference run: ~24%.
#
# Stage 4 uses an overall calibration score threshold instead of the
# resistance mean (the score requires a full calibration_summary.txt).

STAGE_THRESHOLDS: dict[int, dict] = {
    1: {"inf_res_min": 5.0,  "inf_res_max": 32.0},  # pop ~10k
    2: {"inf_res_min": 10.0, "inf_res_max": 30.0},  # pop ~30k
    3: {"inf_res_min": 15.0, "inf_res_max": 27.0},  # pop ~100k
    4: {"score_pool": 1.4,   "score_accept": 1.2},  # pop ~1M — needs score
}

# Simulation start year (must match src/config.rs SIMULATION_START_YEAR = 1930).
SIMULATION_START_YEAR = 1930

# Path to the resistance prevalence targets file.  Used to restrict the
# inf_res_mean computation to only the (bacteria, drug) pairs that have a
# calibration target, matching the set of combinations used by the
# calibration_summary.py "Overall Infection Resistance" block.
RESISTANCE_TARGETS_PATH = Path("data") / "resistance_prevalence_values.csv"

# Calibration window: rows with calendar_year in [2022, 2026).
CALIBRATION_WINDOW_START = 2022
CALIBRATION_WINDOW_END = 2026  # exclusive

# ---------------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------------

def _load_valid_bd_pairs(targets_path: Path) -> frozenset:
    """
    Parse resistance_prevalence_values.csv and return a frozenset of
    (b_slug, d_slug) pairs that have a non-missing resistance target.

    The bacteria slugs use the same column-name format as the simulation
    output CSV (e.g. "citrobacter_spp.", "p_stuartii").  Providencia stuartii
    is stored in the simulation as "p_stuartii" so that name is overridden.
    """
    _B_SLUG_OVERRIDES = {"providencia_stuartii": "p_stuartii"}
    try:
        df = pd.read_csv(targets_path)
    except Exception as exc:
        print(f"WARNING: could not load resistance targets {targets_path}: {exc}", file=sys.stderr)
        return frozenset()
    drug_cols = [c for c in df.columns if c not in ("Bacteria", "notes")]
    pairs: set = set()
    for _, row in df.iterrows():
        b_name = str(row.get("Bacteria", "")).strip()
        if not b_name:
            continue
        b_slug = b_name.strip().lower().replace(" ", "_")
        b_slug = _B_SLUG_OVERRIDES.get(b_slug, b_slug)
        for d_col in drug_cols:
            val = str(row.get(d_col, ".")).strip()
            if val and val not in (".", "nan"):
                pairs.add((b_slug, d_col))
    return frozenset(pairs)


def _compute_inf_res_mean(csv_path: Path, valid_bd_pairs: frozenset) -> Optional[float]:
    """
    Load a simulation_summary CSV, filter to the calibration window, and
    return the mean infection resistance (%) across the target-defined
    (bacteria, drug) pairs.

    The calculation mirrors the logic in calibration_summary.py:
      prevalence_i = sum_window(infected_with_any_r_positive_{b}_{d})
                   / sum_window({b}_currently_infected)
    averaged over the pairs in valid_bd_pairs where the denominator > 0.

    MDR-TB (rifampicin on tuberculosis) and Listeria monocytogenes are
    excluded to match the exclusions applied in calibration_summary.py.
    Only pairs present in the resistance targets file are included, so
    spurious near-zero contributions from unmodelled combinations and
    the hospital/community location-split columns are not counted.
    """
    try:
        df = pd.read_csv(csv_path, low_memory=False)
        # Defragment after read_csv on a very wide CSV to avoid PerformanceWarning
        # when adding new columns below.
        df = df.copy()
    except Exception as exc:
        print(f"  WARNING: could not read {csv_path.name}: {exc}", file=sys.stderr)
        return None

    # Build calendar year column.
    if "time_in_years" not in df.columns:
        if "time_step" in df.columns:
            df["time_in_years"] = pd.to_numeric(df["time_step"], errors="coerce") / 365.0
        else:
            print(f"  WARNING: {csv_path.name} has no time_in_years or time_step column", file=sys.stderr)
            return None

    df["calendar_year"] = SIMULATION_START_YEAR + pd.to_numeric(df["time_in_years"], errors="coerce")

    # Filter to calibration window; fall back to closest available year if window is absent.
    mask = (df["calendar_year"] >= CALIBRATION_WINDOW_START) & (df["calendar_year"] < CALIBRATION_WINDOW_END)
    year_df = df.loc[mask]
    if year_df.empty:
        available = df["calendar_year"].dropna().unique()
        if available.size == 0:
            print(f"  WARNING: {csv_path.name} has no valid calendar year data", file=sys.stderr)
            return None
        nearest = float(min(available, key=lambda v: abs(v - (CALIBRATION_WINDOW_START + CALIBRATION_WINDOW_END) / 2)))
        mask = (df["calendar_year"] >= np.floor(nearest)) & (df["calendar_year"] < np.floor(nearest) + 1.0)
        year_df = df.loc[mask]

    if year_df.empty:
        print(f"  WARNING: {csv_path.name} — no rows in calibration window", file=sys.stderr)
        return None

    # Iterate over the target-defined (bacteria, drug) pairs only.
    # This matches calibration_summary.py's "Overall Infection Resistance" block,
    # which averages over the same set of combinations (1266 after exclusions).
    # Discovered-from-columns approach picked up hospital/community split columns
    # and unmodelled drug/bacteria pairs, all near-zero, diluting the mean to ~8%.
    prevalences: list[float] = []

    for b_slug, d_slug in valid_bd_pairs:
        # Skip Listeria (very low incidence → unstable percentages at small scale).
        if "listeria" in b_slug:
            continue

        # Exclude rifampicin on tuberculosis (MDR-TB hardcoded resistant).
        if "tuberculosis" in b_slug and d_slug == "rifampicin":
            continue

        infected_col = f"{b_slug}_currently_infected"
        if infected_col not in year_df.columns:
            continue

        infected_series = pd.to_numeric(year_df[infected_col], errors="coerce").fillna(0.0)
        total_infected = float(infected_series.sum())
        if total_infected <= 0.0:
            continue

        pos_col = f"{b_slug}_infected_with_any_r_positive_{d_slug}"
        if pos_col not in year_df.columns:
            continue

        positive_series = pd.to_numeric(year_df[pos_col], errors="coerce").fillna(0.0)
        total_positive = float(positive_series.sum())

        prevalence = total_positive / total_infected  # fraction in [0, 1]
        prevalences.append(float(np.clip(prevalence, 0.0, 1.0) * 100.0))

    if not prevalences:
        print(f"  WARNING: {csv_path.name} — no valid bacteria/drug combinations found", file=sys.stderr)
        return None

    return float(np.mean(prevalences))


def _load_sampled_params(params_path: Path) -> Optional[pd.DataFrame]:
    """Load sampled_parameters_{run_id}.csv and pivot to wide format (one row, one col per axis)."""
    if not params_path.exists():
        return None
    try:
        df = pd.read_csv(params_path)
    except Exception:
        return None
    if df.empty or "sampled_quantity" not in df.columns or "sampled_value" not in df.columns:
        return None
    # Pivot: sampled_quantity values become column names.
    row: dict = {}
    if "run_id" in df.columns:
        row["run_id"] = df["run_id"].iloc[0]
    for _, rec in df.iterrows():
        col_name = str(rec["sampled_quantity"])
        row[col_name] = float(rec["sampled_value"])
    return pd.DataFrame([row])


def _load_score_from_summary_txt(txt_path: Path) -> Optional[float]:
    """Extract the overall calibration score from a calibration_summary_*.txt file."""
    if not txt_path.exists():
        return None
    try:
        text = txt_path.read_text(encoding="utf-8", errors="replace")
    except Exception:
        return None
    import re
    match = re.search(r"Overall score:\s*([\d.]+)", text)
    if match:
        return float(match.group(1))
    return None


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Filter calibration sweep runs by stage thresholds.")
    parser.add_argument("--stage", type=int, required=True, choices=[1, 2, 3, 4],
                        help="Funnel stage (1=10k, 2=30k, 3=100k, 4=1M).")
    parser.add_argument("--input-dir", type=Path,
                        default=Path("amr_simulation_output_analysis_outputs"),
                        help="Directory containing simulation_summary_*.csv files.")
    parser.add_argument("--run-ids", type=Path, default=None,
                        help="Optional CSV file with a 'run_id' column listing which run IDs to consider. "
                             "If omitted, all simulation_summary_*.csv files in --input-dir are processed.")
    parser.add_argument("--output-dir", type=Path, default=Path("."),
                        help="Directory to write the survivors CSV (default: current directory).")
    parser.add_argument("--summary-txt-dir", type=Path, default=Path("output_graphs"),
                        help="Directory to search for calibration_summary_*.txt files (stage 4 only).")
    args = parser.parse_args()

    thresholds = STAGE_THRESHOLDS[args.stage]
    input_dir: Path = args.input_dir
    output_dir: Path = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    # Load the set of (b_slug, d_slug) pairs with calibration targets.
    # This restricts inf_res_mean to the same combinations used by
    # calibration_summary.py's Overall Infection Resistance block.
    valid_bd_pairs = _load_valid_bd_pairs(RESISTANCE_TARGETS_PATH)
    if not valid_bd_pairs:
        print("ERROR: resistance targets file not found or empty — cannot compute inf_res_mean.", file=sys.stderr)
        sys.exit(1)

    # Collect candidate CSV paths.
    all_summary_csvs = sorted(input_dir.glob("simulation_summary_*.csv"))
    if not all_summary_csvs:
        print(f"No simulation_summary_*.csv files found in {input_dir}", file=sys.stderr)
        sys.exit(1)

    # Optionally restrict to a specified list of run IDs.
    allowed_ids: Optional[set[str]] = None
    if args.run_ids is not None:
        id_df = pd.read_csv(args.run_ids)
        if "run_id" not in id_df.columns:
            print("--run-ids file must have a 'run_id' column", file=sys.stderr)
            sys.exit(1)
        allowed_ids = set(id_df["run_id"].astype(str))

    candidate_csvs = [
        p for p in all_summary_csvs
        if allowed_ids is None or p.stem.replace("simulation_summary_", "") in allowed_ids
    ]

    print(f"Stage {args.stage}: processing {len(candidate_csvs)} run(s)...")

    rows: list[dict] = []

    for csv_path in candidate_csvs:
        run_id_str = csv_path.stem.replace("simulation_summary_", "")

        if args.stage in (1, 2, 3):
            # Resistance-mean filter.
            inf_res_mean = _compute_inf_res_mean(csv_path, valid_bd_pairs)
            if inf_res_mean is None:
                continue

            lo = thresholds["inf_res_min"]
            hi = thresholds["inf_res_max"]
            passes = lo <= inf_res_mean <= hi

            status = "PASS" if passes else "fail"
            print(f"  {run_id_str}: inf_res_mean={inf_res_mean:.2f}%  [{lo}–{hi}%]  {status}")

            if not passes:
                continue

            row: dict = {"run_id": run_id_str, "inf_res_mean_pct": round(inf_res_mean, 3)}

        else:
            # Stage 4: use overall calibration score from calibration_summary txt.
            txt_path = args.summary_txt_dir / f"calibration_summary_{run_id_str}.txt"
            score = _load_score_from_summary_txt(txt_path)
            if score is None:
                print(f"  {run_id_str}: no calibration_summary txt found — skipping")
                continue

            pool_thresh = thresholds["score_pool"]
            accept_thresh = thresholds["score_accept"]
            in_pool = score <= pool_thresh
            accepted = score <= accept_thresh
            status = "ACCEPTED" if accepted else ("pool" if in_pool else "fail")
            print(f"  {run_id_str}: overall_score={score:.3f}  (pool<{pool_thresh}, accept<{accept_thresh})  {status}")

            if not in_pool:
                continue

            # Also compute inf_res_mean for diversity selection reference.
            inf_res_mean = _compute_inf_res_mean(csv_path, valid_bd_pairs)
            row = {
                "run_id": run_id_str,
                "overall_score": score,
                "accepted": accepted,
                "inf_res_mean_pct": round(inf_res_mean, 3) if inf_res_mean is not None else None,
            }

        # Join sampled parameters if available.
        params_path = input_dir / f"sampled_parameters_{run_id_str}.csv"
        params_row = _load_sampled_params(params_path)
        if params_row is not None:
            for col in params_row.columns:
                if col != "run_id":
                    row[col] = params_row[col].iloc[0]

        rows.append(row)

    if not rows:
        print(f"\nNo survivors at stage {args.stage}.")
        sys.exit(0)

    survivors = pd.DataFrame(rows)
    out_path = output_dir / f"sweep_stage{args.stage}_survivors.csv"
    survivors.to_csv(out_path, index=False)

    print(f"\nStage {args.stage}: {len(survivors)} survivor(s) → {out_path}")

    if args.stage == 4:
        accepted = survivors[survivors.get("accepted", True) == True]  # noqa: E712
        print(f"  Of which {len(accepted)} meet the accept threshold (<{thresholds['score_accept']}).")
        if len(accepted) >= 5:
            print("  5 diverse sets can be selected from this pool.")
        elif len(accepted) > 0:
            print(f"  Only {len(accepted)} accepted — consider relaxing stage 3 threshold and re-running.")
        else:
            print("  No sets meet accept threshold — consider relaxing stage 4 threshold.")


if __name__ == "__main__":
    main()
