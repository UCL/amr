"""
parse_calibration.py

Parse one or more calibration_summary_*.txt files from the Rust AMR simulation
and return structured pandas DataFrames for each section.

Also provides an aggregation function for combining multiple accepted
calibration runs into median (5th–95th percentile) estimates.

Usage
-----
Single run:
    from amr_simulation_output_analysis.parse_calibration import parse_file
    data = parse_file("output_graphs/calibration_summary_958282.txt")

Multiple runs:
    from amr_simulation_output_analysis.parse_calibration import parse_files, aggregate
    runs = parse_files(["output_graphs/calibration_summary_958282.txt", ...])
    agg  = aggregate(runs)   # numeric cells become "median (p5–p95)" strings
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Union

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Encoding correction
# ---------------------------------------------------------------------------

# These substitutions handle any residual encoding artefacts.
# Most characters are read correctly by decode('utf-8'); the entries below
# cover double-encoding edge cases (Windows-1252 bytes written into a UTF-8
# stream) and ASCII control characters.
_ENC_MAP: list[tuple[bytes, str]] = [
    # Pattern                  Replacement   Origin
    (b"\xe2\x80\x94",          "\u2014"),   # â€" → em-dash
    (b"\xe2\x80\x99",          "\u2019"),   # â€™ → right single quote
    (b"\xe2\x80\x98",          "\u2018"),   # â€˜ → left single quote
    (b"\xe2\x80\xa2",          "\u2022"),   # â€¢ → bullet
    (b"\xc3\x97",              "\u00d7"),   # Ã× → ×
    (b"\xc2\xb1",              "\u00b1"),   # Â± → ±
    (b"\x08",                  ""),         # backspace control character
]


def _fix(raw: bytes) -> str:
    """Decode bytes as UTF-8 (falling back to latin-1) and fix known artefacts."""
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")
    # Apply byte-pattern substitutions on the raw bytes first where possible,
    # then do any string-level fixes.
    for pattern, replacement in _ENC_MAP:
        if isinstance(pattern, bytes):
            raw = raw.replace(pattern, replacement.encode("utf-8"))
        else:
            text = text.replace(pattern, replacement)
    # Re-decode after byte-level fixes
    try:
        text = raw.decode("utf-8")
    except UnicodeDecodeError:
        text = raw.decode("latin-1")
    return text


def _read(path: Union[str, Path]) -> list[str]:
    raw = Path(path).read_bytes()
    return _fix(raw).splitlines()


# ---------------------------------------------------------------------------
# Section splitting
# ---------------------------------------------------------------------------

# Ordered longest-first to prevent shorter prefixes shadowing longer ones.
_SECTION_PATTERNS: list[tuple[str, str]] = [
    ("Bacteria Burden Benchmarks — Infections",   "bacteria_infections"),
    ("Bacteria Burden Benchmarks — Mortality",    "bacteria_mortality"),
    ("Serious Resistance Locus",                 "serious_resistance_locus"),
    ("Resistance Incidence Locus",                "resistance_incidence_locus"),
    ("Overall Resistance Fit",                    "overall_resistance_fit"),
    ("Per-Bacteria Mean",                         "resistance_per_bacteria"),
    ("Per-Drug Mean",                             "resistance_per_drug"),
    ("Resistance Benchmark Provenance",            "resistance_provenance"),
    ("Resistance Benchmarks",                     "resistance_benchmarks"),
    ("Headline Metrics",                          "headline_metrics"),
    ("Testing Summary",                           "testing_summary"),
    ("Syndrome Incidence Breakdown",              "syndrome_incidence"),
    ("Infection Death Rates by Age Group and Region", "age_region_death_rates"),
    ("Infection Incidence Fit Summary",           "fit_infection_incidence"),
    ("Microbiome Carriage Fit Summary",           "fit_carriage"),
    ("Infection Deaths Fit Summary",              "fit_infection_deaths"),
    ("Calibration Score",                         "calibration_score"),
    ("Block Scores",                              "block_scores"),
    ("Largest Contributors",                      "largest_contributors"),
    ("Drug Class Share (",                        "drug_class_share"),
    ("Drug Class Share History",                  "drug_class_share_history"),
    ("Overall Infection Resistance",              "overall_resistance_header"),
    ("Microbiome Resistance",                     "microbiome_resistance"),
    ("Footnotes",                                 "footnotes"),
]


def _split_sections(lines: list[str]) -> dict[str, list[str]]:
    result: dict[str, list[str]] = {"_header": []}
    current = "_header"
    for line in lines:
        s = line.strip()
        matched = None
        for prefix, key in _SECTION_PATTERNS:
            if s.startswith(prefix):
                matched = key
                break
        if matched:
            result.setdefault(matched, [])
            current = matched
        result.setdefault(current, []).append(line)
    return result


# ---------------------------------------------------------------------------
# Generic table parser
# ---------------------------------------------------------------------------

def _split_row(line: str) -> list[str]:
    """Split a fixed-width row on two or more consecutive spaces."""
    return re.split(r"\s{2,}", line.strip())


def _table_from_section(section_lines: list[str]) -> pd.DataFrame:
    """
    Find the column-header row (first non-empty line after the section title
    that contains at least one 2+-space gap), then parse all following data rows.
    """
    header_idx = None
    for idx, line in enumerate(section_lines[1:], start=1):
        s = line.strip()
        if not s:
            continue
        # Must contain a 2+-space gap AND produce multiple columns when split;
        # this skips multi-line section-title continuation lines that happen to
        # start with indentation (e.g. the "Overall Resistance Fit" preamble).
        if "  " in line and len(_split_row(line)) > 1:
            header_idx = idx
            break
    if header_idx is None:
        return pd.DataFrame()

    headers = _split_row(section_lines[header_idx])
    rows: list[list[str]] = []
    for line in section_lines[header_idx + 1:]:
        s = line.strip()
        if not s:
            continue
        if s.startswith("Note:") or s.startswith("Observation") or s.startswith("*"):
            continue
        parts = _split_row(line)
        if not any(parts):
            continue
        # Pad or truncate to match header width
        while len(parts) < len(headers):
            parts.append("")
        rows.append(parts[:len(headers)])

    return pd.DataFrame(rows, columns=headers) if rows else pd.DataFrame(columns=headers)


def _coerce_numeric(df: pd.DataFrame, skip: list[str]) -> pd.DataFrame:
    """
    Attempt to convert each non-skip column to float.
    Null sentinels ("---", "-", "—", "") become NaN.
    Non-numeric strings are kept as strings (not silently dropped).
    """
    skip_set = set(skip)
    df = df.copy()
    for col in df.columns:
        if col in skip_set:
            continue

        def _convert(v: object) -> object:
            if isinstance(v, float):
                return v
            s = str(v).strip()
            if s in ("---", "-", "", "—", "N/A"):
                return np.nan
            try:
                return float(s.replace(",", ""))
            except ValueError:
                return v  # keep as string

        df[col] = df[col].apply(_convert)
    return df


# ---------------------------------------------------------------------------
# Resistance benchmarks - wide table with provenance columns
# ---------------------------------------------------------------------------

_BENCH_COLS = [
    "Bacteria", "Drug", "Class",
    "Inf sim (%)", "Inf target (%)",
    "Inf provenance", "Inf source", "Inf rationale",
    "Avg sim (%)", "Avg target (%)",
    "Avg provenance", "Avg source", "Avg rationale",
    "Micro sim (%)",
    "Inf days", "Res days", "Carrier days", "Flags",
]


def _parse_resistance_benchmarks(section_lines: list[str]) -> pd.DataFrame:
    dynamic = _table_from_section(section_lines)
    if not dynamic.empty:
        return dynamic

    rows: list[list[str]] = []
    for line in section_lines[1:]:
        s = line.strip()
        if not s:
            continue
        if s.startswith("Bacteria") or s.startswith("Note:") or s.startswith("Observation"):
            continue
        parts = _split_row(s)
        if len(parts) < 3:
            continue
        while len(parts) < len(_BENCH_COLS):
            parts.append("")
        rows.append(parts[:len(_BENCH_COLS)])
    if not rows:
        return pd.DataFrame(columns=_BENCH_COLS)
    return pd.DataFrame(rows, columns=_BENCH_COLS)


# ---------------------------------------------------------------------------
# Bullet-list parsers (fit summaries, calibration score)
# ---------------------------------------------------------------------------

def _parse_bullets(section_lines: list[str]) -> list[str]:
    return [
        line.strip()[2:]
        for line in section_lines[1:]
        if line.strip().startswith("- ")
    ]


def _parse_cal_score(section_lines: list[str]) -> dict[str, str]:
    result: dict[str, str] = {}
    in_failed = False
    failed_gates: list[str] = []
    for line in section_lines[1:]:
        s = line.strip()
        if s.startswith("- Failed gates:"):
            in_failed = True
            continue
        if in_failed:
            if s.startswith("-"):
                failed_gates.append(s[1:].strip())
            elif s.startswith("  "):
                failed_gates.append(s.strip())
            else:
                in_failed = False
        if s.startswith("- ") and ":" in s and not in_failed:
            kv = s[2:].split(":", 1)
            result[kv[0].strip()] = kv[1].strip()
    if failed_gates:
        result["Failed gates"] = "; ".join(failed_gates)
    return result


# ---------------------------------------------------------------------------
# Main parse_file
# ---------------------------------------------------------------------------

def parse_file(path: Union[str, Path]) -> dict:
    """
    Parse a calibration_summary_*.txt file into structured DataFrames.

    Returns a dict with keys:
        meta, headline_metrics, testing_summary,
        bacteria_infections, bacteria_mortality, resistance_incidence_locus,
        syndrome_incidence, fit_stats, calibration_score, block_scores,
        largest_contributors, drug_class_share, overall_resistance_fit,
        resistance_per_bacteria, resistance_per_drug,
        microbiome_resistance (float), resistance_benchmarks
    """
    run_id = re.sub(r"^calibration_summary_", "", Path(path).stem)
    lines = _read(path)
    sec = _split_sections(lines)

    # --- meta ----------------------------------------------------------------
    meta: dict[str, str] = {"run_id": run_id, "source_file": str(path)}
    for line in sec.get("_header", []):
        s = line.strip()
        for raw_key, label in [
            ("Target year:",                "target_year"),
            ("Calibration window duration:", "window_duration"),
            ("Mean simulated population",   "mean_pop"),
            ("Final simulated population",  "final_pop"),
            ("Population scale factor",     "scale_factor"),
        ]:
            if s.startswith(raw_key):
                meta[label] = s.split(":", 1)[1].strip()

    # --- generic table sections ----------------------------------------------
    _id_cols: dict[str, list[str]] = {
        "headline_metrics":           ["Metric"],
        "testing_summary":            ["Metric"],
        "bacteria_infections":        ["Bacteria"],
        "bacteria_mortality":         ["Bacteria"],
        "resistance_incidence_locus": ["Bacteria"],
        "serious_resistance_locus":   ["Bacteria"],
        "syndrome_incidence":         ["Syndrome"],
        "block_scores":               ["Block"],
        "largest_contributors":       ["Block"],
        "drug_class_share":           ["Class"],
        "drug_class_share_history":   ["Class"],
        "overall_resistance_fit":     ["Component"],
        "resistance_per_bacteria":    ["Bacteria"],
        "resistance_per_drug":        ["Drug"],
        "resistance_provenance":      ["Component", "Provenance class"],
    }

    parsed: dict = {"meta": meta}
    for section_key, id_cols in _id_cols.items():
        df = _table_from_section(sec.get(section_key, []))
        if not df.empty:
            df = _coerce_numeric(df, skip=id_cols)
        parsed[section_key] = df

    # Backward compatibility for summaries written before the dedicated
    # calibration-window drug-class section was introduced.
    if parsed["drug_class_share"].empty and not parsed["drug_class_share_history"].empty:
        parsed["drug_class_share"] = parsed["drug_class_share_history"].copy()

    # --- resistance benchmarks (wide table) ----------------------------------
    bench = _parse_resistance_benchmarks(sec.get("resistance_benchmarks", []))
    bench = _coerce_numeric(bench, skip=["Bacteria", "Drug", "Class", "Flags"])
    parsed["resistance_benchmarks"] = bench

    # --- fit stats (bullet lists) --------------------------------------------
    parsed["fit_stats"] = {
        "infection_incidence": _parse_bullets(sec.get("fit_infection_incidence", [])),
        "carriage":            _parse_bullets(sec.get("fit_carriage", [])),
        "infection_deaths":    _parse_bullets(sec.get("fit_infection_deaths", [])),
    }

    # --- calibration score ---------------------------------------------------
    parsed["calibration_score"] = _parse_cal_score(sec.get("calibration_score", []))

    # --- microbiome resistance (scalar) --------------------------------------
    micro_val: float = np.nan
    for line in sec.get("microbiome_resistance", []):
        m = re.search(r"(\d+\.\d+)", line)
        if m:
            micro_val = float(m.group(1))
            break
    parsed["microbiome_resistance"] = micro_val

    return parsed


def parse_files(paths: list[Union[str, Path]]) -> list[dict]:
    return [parse_file(p) for p in paths]


# ---------------------------------------------------------------------------
# Aggregation
# ---------------------------------------------------------------------------

def _is_numeric(v: object) -> bool:
    return isinstance(v, (int, float)) and not (isinstance(v, float) and np.isnan(v))


def _fmt_val(v: float, mag: float) -> str:
    """Format a single float value with magnitude-appropriate precision."""
    if mag >= 10_000:
        return f"{v:,.0f}"
    if mag >= 1:
        return f"{v:.1f}"
    return f"{v:.2f}"


def _fmt_agg(values: list[float]) -> str:
    """Format a list of numeric values as 'median (p5–p95)' or plain value."""
    arr = np.array([v for v in values if _is_numeric(v)], dtype=float)
    if len(arr) == 0:
        return "—"
    median = float(np.median(arr))
    mag = abs(median) if abs(median) > 0 else 1.0
    if len(arr) == 1:
        return _fmt_val(median, mag)
    p5  = float(np.percentile(arr, 5))
    p95 = float(np.percentile(arr, 95))
    return f"{_fmt_val(median, mag)} ({_fmt_val(p5, mag)}–{_fmt_val(p95, mag)})"


def _agg_dataframes(
    dfs: list[pd.DataFrame],
    key_cols: list[str],
    passthrough_cols: list[str] | None = None,
) -> pd.DataFrame:
    """
    Aggregate N DataFrames (same structure) aligned on key_cols.
    Numeric columns → "median (p5–p95)" strings.
    Non-numeric columns → taken from the first run.
    passthrough_cols: columns to copy verbatim from the first run (no aggregation).
    """
    dfs = [d for d in dfs if d is not None and not d.empty]
    if not dfs:
        return pd.DataFrame()
    ref = dfs[0]
    n_rows = len(ref)
    key_set = set(key_cols)
    passthrough_set = set(passthrough_cols or [])
    result = ref[key_cols].copy().reset_index(drop=True)

    for col in ref.columns:
        if col in key_set:
            continue
        if col in passthrough_set:
            result[col] = ref[col].reset_index(drop=True)
            continue
        # Gather this column from every run, aligned by row index
        col_per_run: list[list] = []
        for df in dfs:
            if col in df.columns:
                col_per_run.append(df[col].tolist()[:n_rows])
            else:
                col_per_run.append([np.nan] * n_rows)

        aggregated: list[str] = []
        for i in range(n_rows):
            numeric_vals = [
                float(col_per_run[r][i])
                for r in range(len(col_per_run))
                if i < len(col_per_run[r]) and _is_numeric(col_per_run[r][i])
            ]
            if numeric_vals:
                aggregated.append(_fmt_agg(numeric_vals))
            else:
                # Fallback: use first run's string value
                fallback = col_per_run[0][i] if col_per_run and i < len(col_per_run[0]) else ""
                aggregated.append(str(fallback) if fallback not in (None, np.nan) else "—")
        result[col] = aggregated

    return result


def aggregate(parsed_list: list[dict]) -> dict:
    """
    Aggregate N parsed run dicts into a single dict.
    For N=1 runs, values are formatted without confidence intervals.
    For N>1 runs, numeric cells become "median (p5–p95)" strings.
    """
    if not parsed_list:
        return {}

    _df_key_cols: dict[str, list[str]] = {
        "headline_metrics":           ["Metric"],
        "testing_summary":            ["Metric"],
        "bacteria_infections":        ["Bacteria"],
        "bacteria_mortality":         ["Bacteria"],
        "resistance_incidence_locus": ["Bacteria"],
        "serious_resistance_locus":   ["Bacteria"],
        "syndrome_incidence":         ["Syndrome"],
        "block_scores":               ["Block"],
        "largest_contributors":       ["Block"],
        "drug_class_share":           ["Class"],
        "overall_resistance_fit":     ["Component"],
        "resistance_per_bacteria":    ["Bacteria"],
        "resistance_per_drug":        ["Drug"],
        "resistance_benchmarks":      ["Bacteria", "Drug"],
    }

    # Columns that are fixed calibration targets (not simulation outputs);
    # copy verbatim from the first run rather than aggregating.
    _passthrough_cols: dict[str, list[str]] = {
        "headline_metrics":    ["Target", "Unit"],
        "bacteria_infections": ["Infection target (%)", "Carriage target (%)"],
        "bacteria_mortality":  ["Deaths target (millions)"],
        "drug_class_share":    ["Target 2025 (%)", "Target 2000 (%)", "Target 1975 (%)", "Target 1950 (%)"],
    }

    agg: dict = {
        "meta":   parsed_list[0]["meta"],
        "n_runs": len(parsed_list),
    }

    for section, key_cols in _df_key_cols.items():
        dfs = [p.get(section, pd.DataFrame()) for p in parsed_list]
        pt = _passthrough_cols.get(section, [])
        agg[section] = _agg_dataframes(dfs, key_cols, pt)

    # Microbiome resistance scalar
    vals = [p.get("microbiome_resistance", np.nan) for p in parsed_list]
    numeric = [v for v in vals if _is_numeric(v)]
    agg["microbiome_resistance"] = _fmt_agg(numeric) if numeric else "—"

    # Pass-through fields (use first run)
    agg["fit_stats"]            = parsed_list[0].get("fit_stats", {})
    agg["calibration_score"]    = parsed_list[0].get("calibration_score", {})
    agg["calibration_scores_all"] = [p.get("calibration_score", {}) for p in parsed_list]

    return agg
