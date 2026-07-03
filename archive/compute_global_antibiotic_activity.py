#!/usr/bin/env python3
"""Compute Global Antibiotic Activity (GAA) from simulation CSV output.

For each bacterium b on each simulated day t:

    GAA(b, t) = Σ_{drugs d: intro_date(d) ≤ t}  potency(b, d) × (1 − mean_any_r(b, d, t))

where
    mean_any_r(b, d, t) = {b}_sum_any_r_{d}(t)  /  {b}_currently_infected(t)

    (treated as 0 when currently_infected = 0, i.e. full baseline activity assumed
     when no infections are present to carry resistance)

All static parameters are parsed live from src/config.rs:
  • POTENCY_EMBEDDED_DATA       →  potency(b, d)
  • DRUG_INTRODUCTION_DATES     →  intro_date(d)  [days since 1930-01-01]

Only CalibrationMode::None (and Partial) runs write the sum_any_r columns; Full
calibration runs skip the expensive B×D loops (need_full_summary = false for
the pre-window period) so those CSVs lack the necessary columns.

Usage
-----
    python compute_global_antibiotic_activity.py \\
        amr_simulation_output_analysis_outputs/simulation_summary_123456.csv \\
        [--output gaa_123456.csv]

Or import and call compute_gaa() programmatically.
"""

import re
import sys
import argparse
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Drug column order (positions 0–51) as they appear in POTENCY_EMBEDDED_DATA
# and in the CSV column names.  Must match DRUG_SHORT_NAMES in config.rs.
# ---------------------------------------------------------------------------
DRUG_NAMES = [
    "sulfanilamide",            # 0
    "penicillin_g",             # 1
    "ampicillin",               # 2
    "amoxicillin",              # 3
    "piperacillin",             # 4
    "ticarcillin",              # 5
    "cephalexin",               # 6
    "cefazolin",                # 7
    "cefuroxime",               # 8
    "ceftriaxone",              # 9
    "ceftazidime",              # 10
    "cefepime",                 # 11
    "ceftaroline",              # 12
    "meropenem",                # 13
    "imipenem_c",               # 14
    "ertapenem",                # 15
    "aztreonam",                # 16
    "erythromycin",             # 17
    "azithromycin",             # 18
    "clarithromycin",           # 19
    "clindamycin",              # 20
    "gentamicin",               # 21
    "tobramycin",               # 22
    "amikacin",                 # 23
    "ciprofloxacin",            # 24
    "levofloxacin",             # 25
    "moxifloxacin",             # 26
    "ofloxacin",                # 27
    "tetracycline",             # 28
    "doxycycline",              # 29
    "minocycline",              # 30
    "vancomycin",               # 31
    "teicoplanin",              # 32
    "dalbavancin",              # 33
    "linezolid",                # 34
    "tedizolid",                # 35
    "quinu_dalfo",              # 36
    "trim_sulf",                # 37
    "chloramphenicol",          # 38
    "nitrofurantoin",           # 39
    "retapamulin",              # 40
    "fusidic_a",                # 41
    "metronidazole",            # 42
    "furazolidone",             # 43
    "rifampicin",               # 44
    "amoxicillin_clavulanate",  # 45
    "piperacillin_tazobactam",  # 46
    "ampicillin_sulbactam",     # 47
    "ticarcillin_clavulanate",  # 48
    "ceftazidime_avibactam",    # 49
    "meropenem_vaborbactam",    # 50
    "colistin",                 # 51
    "nalidixic_acid",           # 52  (first-generation quinolone; historical UTI/GI use 1963–~1990)
]

N_DRUGS = len(DRUG_NAMES)


# ---------------------------------------------------------------------------
# Parsing helpers
# ---------------------------------------------------------------------------

def _extract_block(text: str, start_marker: str, end_pattern: str) -> str:
    """Return the substring of *text* starting at *start_marker* up to the
    first match of *end_pattern* (inclusive)."""
    idx = text.index(start_marker)
    m = re.search(end_pattern, text[idx:])
    if m is None:
        raise ValueError(f"Could not find end of block starting with {start_marker!r}")
    return text[idx: idx + m.end()]


def parse_potency_matrix(config_rs_path: str | Path) -> dict[str, np.ndarray]:
    """Return {bacteria_name: float64 array of length 52} parsed from
    POTENCY_EMBEDDED_DATA in config.rs.  None entries become 0.0."""
    text = Path(config_rs_path).read_text(encoding="utf-8")

    # Isolate the POTENCY_EMBEDDED_DATA const (ends with the line "];")
    start_marker = "const POTENCY_EMBEDDED_DATA: &[(&str, [Option<f64>; 52])] = &["
    block = _extract_block(text, start_marker, r"\n\];\n")

    # Each entry: ("bacteria_name", [ Some(v), ..., None, ... ])
    # Split on the outer tuple boundaries
    entry_pattern = re.compile(
        r'"([^"]+)",\s*\[([^\]]+)\]',
        re.DOTALL,
    )
    value_pattern = re.compile(r'Some\(([\d.eE+\-]+)\)|None')

    result: dict[str, np.ndarray] = {}
    for m in entry_pattern.finditer(block):
        bact_name = m.group(1)
        values_str = m.group(2)
        vals: list[float] = []
        for vm in value_pattern.finditer(values_str):
            vals.append(float(vm.group(1)) if vm.group(1) is not None else 0.0)
        if len(vals) != N_DRUGS:
            raise ValueError(
                f"Bacteria {bact_name!r}: expected {N_DRUGS} potency values, "
                f"got {len(vals)}"
            )
        result[bact_name] = np.array(vals, dtype=np.float64)

    return result


def parse_drug_intro_dates(config_rs_path: str | Path) -> dict[str, int]:
    """Return {drug_name: time_step} parsed from DRUG_INTRODUCTION_DATES in
    config.rs.  Time steps are days since 1930-01-01."""
    text = Path(config_rs_path).read_text(encoding="utf-8")
    start_marker = "pub static ref DRUG_INTRODUCTION_DATES"
    block = _extract_block(text, start_marker, r"\n    \};\n\}")

    pattern = re.compile(r'map\.insert\("([^"]+)",\s*(\d+)\)')
    return {m.group(1): int(m.group(2)) for m in pattern.finditer(block)}


# ---------------------------------------------------------------------------
# Main computation
# ---------------------------------------------------------------------------

def compute_gaa(
    csv_path: str | Path,
    config_rs_path: str | Path = None,
    output_path: str | Path = None,
    *,
    potency_matrix: dict[str, np.ndarray] = None,
    drug_intro_dates: dict[str, int] = None,
) -> pd.DataFrame:
    """Compute Global Antibiotic Activity for every row in *csv_path*.

    Parameters
    ----------
    csv_path:
        Path to a simulation_summary CSV (CalibrationMode::None or Partial run).
    config_rs_path:
        Path to ``src/config.rs``.  Defaults to ``src/config.rs`` relative to
        the directory containing this script.
    output_path:
        If given, write the result CSV here.
    potency_matrix, drug_intro_dates:
        Pre-parsed static tables; if supplied, *config_rs_path* is ignored.

    Returns
    -------
    pd.DataFrame with columns:
        time_step, policy_option, run_id, time_in_years,
        {bacteria}_gaa  for every bacterium present in the CSV.
    """
    csv_path = Path(csv_path)
    if potency_matrix is None or drug_intro_dates is None:
        if config_rs_path is None:
            config_rs_path = Path(__file__).parent / "src" / "config.rs"
        config_rs_path = Path(config_rs_path)
        potency_matrix = parse_potency_matrix(config_rs_path)
        drug_intro_dates = parse_drug_intro_dates(config_rs_path)

    # ------------------------------------------------------------------
    # 1.  Selective column loading (only what we need)
    # ------------------------------------------------------------------
    def _keep_col(name: str) -> bool:
        return (
            name in {"time_step", "policy_option", "run_id", "time_in_years"}
            or name.endswith("_currently_infected")
            or "_sum_any_r_" in name
        )

    print(f"Loading columns from {csv_path.name} …", flush=True)
    df = pd.read_csv(csv_path, usecols=_keep_col, low_memory=False)
    print(f"  {len(df):,} rows loaded.", flush=True)

    # ------------------------------------------------------------------
    # 2.  Identify bacteria present in this CSV
    # ------------------------------------------------------------------
    infected_cols = [
        c for c in df.columns
        if c.endswith("_currently_infected") and c != "total_currently_infected"
    ]
    bacteria_names = [c[: -len("_currently_infected")] for c in infected_cols]

    # Pre-build drug intro array aligned to DRUG_NAMES
    intro_ts = np.array(
        [drug_intro_dates.get(d, 10**9) for d in DRUG_NAMES],
        dtype=np.int64,
    )  # shape (n_drugs,)

    time_steps = df["time_step"].to_numpy(dtype=np.int64)  # shape (n_rows,)

    # ------------------------------------------------------------------
    # 3.  Compute GAA per bacterium (vectorised over rows and drugs)
    # ------------------------------------------------------------------
    # drug_available[row, drug] = 1 if intro_ts[drug] <= time_steps[row], else 0
    drug_available = (time_steps[:, None] >= intro_ts[None, :]).astype(
        np.float64
    )  # shape (n_rows, n_drugs)

    gaa_dict: dict[str, np.ndarray] = {}

    for bact in bacteria_names:
        if bact not in potency_matrix:
            # Bacteria present in CSV but absent from potency table — skip
            continue

        potency = potency_matrix[bact]  # shape (n_drugs,)
        count = df[f"{bact}_currently_infected"].to_numpy(
            dtype=np.float64
        )  # shape (n_rows,)
        safe_count = np.where(count > 0, count, np.nan)

        # Accumulate sum of  potency[d] * (1 - mean_any_r[d]) * drug_available[d]
        # across all drugs.  We iterate over drugs rather than loading a huge
        # 2-D array all at once, to keep peak memory low.
        gaa = np.zeros(len(df), dtype=np.float64)

        for d_idx, drug in enumerate(DRUG_NAMES):
            pot = potency[d_idx]
            if pot == 0.0:
                continue  # Drug has no activity against this bacterium

            col = f"{bact}_sum_any_r_{drug}"
            if col not in df.columns:
                # Column absent (e.g. filtered CSV) – treat as zero resistance
                contribution = pot * drug_available[:, d_idx]
            else:
                sum_any_r = df[col].to_numpy(dtype=np.float64)
                mean_any_r = np.where(count > 0, sum_any_r / safe_count, 0.0)
                mean_any_r = mean_any_r.clip(0.0, 1.0)  # guard numerical noise
                contribution = pot * (1.0 - mean_any_r) * drug_available[:, d_idx]

            gaa += contribution

        gaa_dict[f"{bact}_gaa"] = gaa

    # ------------------------------------------------------------------
    # 4.  Assemble output DataFrame
    # ------------------------------------------------------------------
    id_cols = ["time_step", "policy_option", "run_id", "time_in_years"]
    id_cols_present = [c for c in id_cols if c in df.columns]
    out = df[id_cols_present].copy()
    for col, arr in gaa_dict.items():
        out[col] = arr

    if output_path is not None:
        output_path = Path(output_path)
        out.to_csv(output_path, index=False)
        print(f"GAA written to {output_path}", flush=True)

    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main(argv=None):
    parser = argparse.ArgumentParser(
        description="Compute Global Antibiotic Activity from a simulation summary CSV."
    )
    parser.add_argument("csv", help="Path to simulation_summary CSV")
    parser.add_argument(
        "--config",
        default=None,
        help="Path to src/config.rs (default: src/config.rs alongside this script)",
    )
    parser.add_argument(
        "--output",
        "-o",
        default=None,
        help="Output CSV path (default: <csv_stem>_gaa.csv in same directory)",
    )
    args = parser.parse_args(argv)

    csv_path = Path(args.csv)
    output_path = (
        Path(args.output)
        if args.output
        else csv_path.with_name(csv_path.stem + "_gaa.csv")
    )
    config_path = args.config  # None → default inside compute_gaa

    compute_gaa(csv_path, config_path, output_path)


if __name__ == "__main__":
    main()
