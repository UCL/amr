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
    ("Deaths among people actively infected, by bacterium and home region", "deaths_among_infected_by_bacterium_region"),
    ("Bacteria Burden Benchmarks — Infections",   "bacteria_infections"),
    ("Bacteria Burden Benchmarks — Mortality",    "bacteria_mortality"),
    ("Serious Resistance Locus",                 "serious_resistance_locus"),
    ("Resistance Incidence Locus",                "resistance_incidence_locus"),
    ("Overall Resistance Fit",                    "overall_resistance_fit"),
    ("Regional Resistance Summary",               "regional_resistance"),
    ("Overall Infection Incidence by Region",      "infection_incidence_by_region"),
    ("Per-Bacteria Mean",                         "resistance_per_bacteria"),
    ("Per-Drug Mean",                             "resistance_per_drug"),
    ("Resistance Benchmark Provenance",            "resistance_provenance"),
    ("Resistance Benchmarks",                     "resistance_benchmarks"),
    ("Headline Metrics",                          "headline_metrics"),
    ("Testing Summary",                           "testing_summary"),
    ("Syndrome Incidence Breakdown",              "syndrome_incidence"),
    ("Infection Death Rates by Age Group and Region", "age_region_death_rates"),
    ("Infection Death Counts by Age Group",        "infection_deaths_by_age"),
    ("Infection Death Counts by Region",           "infection_deaths_by_region"),
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
    """Split a compact row whose nonempty cells have whitespace delimiters."""
    return re.split(r"\s{2,}", line.strip())


def _fixed_width_boundaries(header: str, rows: list[str]) -> list[int] | None:
    """Locate column separators without discarding blank cells or header padding.

    Headers and values can have different alignments, so header token starts are
    not cell boundaries. A separator must remain blank in every row. Require one
    unambiguous separator in each header gap; compact legacy tables need not have
    consistent character positions and are handled separately.
    """
    labels = list(re.finditer(r"\S.*?(?=\s{2,}|$)", header))
    if len(labels) < 2 or not rows:
        return None

    occupied = bytearray(len(header))
    for row in rows:
        for value in re.finditer(r"\S+", row):
            start, end = value.span()
            end = min(end, len(header))
            if start < end:
                occupied[start:end] = b"\x01" * (end - start)

    boundaries = [0]
    for left, right in zip(labels, labels[1:]):
        gap = occupied[left.end():right.start()]
        separators = list(re.finditer(b"\x00{2,}", gap))
        if len(separators) != 1:
            return None
        boundaries.append(left.end() + separators[0].start())

    if any(len(_split_row(row)) == len(labels) and len(row) < boundaries[-1] for row in rows):
        return None

    # Whitespace in compact rows can coincide with every header gap. Genuine
    # fixed-width columns also align their values with one edge of the header;
    # check that before interpreting wide delimiters as empty cells.
    ends = [*boundaries[1:], None]
    for label, start, end in zip(labels, boundaries, ends):
        left_aligned = right_aligned = True
        header_end = label.start() + len(label.group().rstrip())
        for row in rows:
            cell = row[start:end]
            if not cell.strip():
                continue
            left_aligned &= start + len(cell) - len(cell.lstrip()) == label.start()
            right_aligned &= start + len(cell.rstrip()) == header_end
            if not left_aligned and not right_aligned:
                return None
    return boundaries


def _table_from_section(
    section_lines: list[str], *, require_complete_rows: bool = False,
) -> pd.DataFrame:
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

    header = section_lines[header_idx]
    headers = _split_row(header)
    data_lines: list[str] = []
    for line in section_lines[header_idx + 1:]:
        s = line.strip()
        if not s:
            continue
        if s.startswith("Note:") or s.startswith("Observation") or s.startswith("*"):
            continue
        if re.fullmatch(r"[-=]{3,}", s):
            continue
        data_lines.append(line)

    boundaries = _fixed_width_boundaries(header, data_lines)
    if boundaries is not None and len(boundaries) != len(headers):
        boundaries = None
    rows: list[list[str]] = []
    for line in data_lines:
        if boundaries is not None:
            ends = [*boundaries[1:], None]
            rows.append([line[start:end].strip() for start, end in zip(boundaries, ends)])
            continue

        parts = _split_row(line)
        if not any(parts):
            continue
        if require_complete_rows and len(parts) != len(headers):
            raise ValueError(
                "Cannot determine calibration table column boundaries: "
                f"expected {len(headers)} cells, found {len(parts)} in {line.strip()!r}"
            )
        # Preserve compact legacy/mixed-section handling when no fixed layout exists.
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
    dynamic = _table_from_section(section_lines, require_complete_rows=True)
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
        largest_contributors, drug_class_share, overall_resistance_fit, regional_resistance,
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
            ("Simulation source CSV:",       "simulation_source_csv"),
            ("Simulation summary schema:",   "simulation_summary_schema"),
            ("Legacy compatibility:",        "legacy_compatibility"),
            ("Target year:",                "target_year"),
            ("Comparison resistance target set:", "comparison_resistance_target_set"),
            ("Comparison resistance target file:", "comparison_resistance_target_file"),
            ("Comparison resistance manifest SHA-256:", "comparison_resistance_manifest_sha256"),
            ("Calibration window duration:", "window_duration"),
            ("Mean simulated population",   "mean_pop"),
            ("Final simulated population",  "final_pop"),
            ("Population scale factor",     "scale_factor"),
        ]:
            if s.startswith(raw_key):
                meta[label] = s.split(":", 1)[1].strip()

    # --- generic table sections ----------------------------------------------
    _id_cols: dict[str, list[str]] = {
        "deaths_among_infected_by_bacterium_region": ["Bacterium"],
        "headline_metrics":           ["Metric"],
        "testing_summary":            ["Metric"],
        "infection_deaths_by_age":    ["Age Group"],
        "infection_deaths_by_region": ["Region"],
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
        "regional_resistance":        ["Region"],
        "infection_incidence_by_region": ["Region"],
        "resistance_per_bacteria":    ["Bacteria"],
        "resistance_per_drug":        ["Drug"],
        "resistance_provenance":      ["Component", "Provenance class"],
    }

    parsed: dict = {"meta": meta}
    for section_key, id_cols in _id_cols.items():
        df = _table_from_section(sec.get(section_key, []))
        if not df.empty:
            df = _coerce_numeric(df, skip=id_cols)
        if section_key in ("infection_deaths_by_age", "infection_deaths_by_region"):
            # Preserve the actual window and policy: target year plus duration
            # alone cannot distinguish e.g. 2022-2025 from 2023-2026.
            observations = [
                line.strip().split(":", 1)[1].strip()
                for line in sec.get(section_key, [])
                if line.strip().startswith("Observation window:")
            ]
            if len(observations) == 1:
                df.attrs["observation_window"] = observations[0]
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
    Aggregate matching keys across runs, retaining their first-seen order.

    Rows and columns present in only some runs are retained. Missing observations
    do not contribute numeric values. Keys must be present and unique within each
    run. The run-specific trailing " *" on Bacteria labels is ignored for identity;
    the first observed display label is retained.

    Numeric columns → "median (p5–p95)" strings.
    Non-numeric columns → the first run containing that key and column.
    passthrough_cols: copy that first cell verbatim, including explicit missingness.
    """
    if not key_cols or len(set(key_cols)) != len(key_cols):
        raise ValueError("Aggregation requires nonempty, distinct key columns")

    key_set = set(key_cols)
    passthrough_set = set(passthrough_cols or [])
    value_columns: dict[str, None] = {}
    rows_by_key: dict[tuple, list[dict]] = {}
    for run_number, df in enumerate(dfs, start=1):
        if df is None or df.empty:
            continue
        if not df.columns.is_unique:
            raise ValueError(f"Run {run_number} contains duplicate column names")
        missing_keys = [col for col in key_cols if col not in df.columns]
        if missing_keys:
            raise ValueError(f"Run {run_number} is missing key columns: {missing_keys}")
        if df[key_cols].isna().any().any():
            raise ValueError(f"Run {run_number} contains null key values in {key_cols}")

        for col in df.columns:
            if col not in key_set:
                value_columns.setdefault(col, None)

        seen: set[tuple] = set()
        for row in df.to_dict(orient="records"):
            key_values = []
            for col in key_cols:
                value = row[col]
                if isinstance(value, str):
                    value = value.strip()
                    if col == "Bacteria":
                        value = value.removesuffix(" *").strip()
                    if not value:
                        raise ValueError(f"Run {run_number} contains a blank key in {col!r}")
                key_values.append(value)
            key = tuple(key_values)
            try:
                duplicate = key in seen
                seen.add(key)
            except TypeError as error:
                raise ValueError(f"Run {run_number} contains an unhashable key: {key!r}") from error
            if duplicate:
                raise ValueError(f"Run {run_number} contains duplicate keys for {key_cols}: {key!r}")
            rows_by_key.setdefault(key, []).append(row)

    if not rows_by_key:
        return pd.DataFrame()

    records = []
    for matching_rows in rows_by_key.values():
        record = {col: matching_rows[0][col] for col in key_cols}
        for col in value_columns:
            values = [row[col] for row in matching_rows if col in row]
            first_value = values[0] if values else np.nan
            if col in passthrough_set:
                record[col] = first_value
                continue
            numeric_vals = [float(value) for value in values if _is_numeric(value)]
            if numeric_vals:
                record[col] = _fmt_agg(numeric_vals)
            else:
                record[col] = "—" if pd.isna(first_value) else str(first_value)
        records.append(record)

    return pd.DataFrame(records, columns=[*key_cols, *value_columns])


def aggregate(parsed_list: list[dict]) -> dict:
    """
    Aggregate N parsed run dicts into a single dict.
    Each table matches rows by its identifying columns, independently of sorting.
    A cell with one observation is formatted without an interval; multiple
    observations become "median (p5–p95)" strings.
    """
    if not parsed_list:
        return {}

    _df_key_cols: dict[str, list[str]] = {
        "deaths_among_infected_by_bacterium_region": ["Bacterium"],
        "infection_deaths_by_age":    ["Age Group"],
        "infection_deaths_by_region": ["Region"],
        "headline_metrics":           ["Metric"],
        "testing_summary":            ["Metric"],
        "bacteria_infections":        ["Bacteria"],
        "bacteria_mortality":         ["Bacteria"],
        "resistance_incidence_locus": ["Bacteria"],
        "serious_resistance_locus":   ["Bacteria"],
        "syndrome_incidence":         ["Syndrome"],
        "block_scores":               ["Block"],
        "largest_contributors":       ["Block", "Target"],
        "drug_class_share":           ["Class"],
        "overall_resistance_fit":     ["Component"],
        "resistance_per_bacteria":    ["Bacteria"],
        "resistance_per_drug":        ["Drug"],
        "resistance_benchmarks":      ["Bacteria", "Drug"],
    }

    # Columns that are fixed calibration targets (not simulation outputs);
    # copy the first matching row/column cell verbatim rather than aggregating.
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
        try:
            agg[section] = _agg_dataframes(dfs, key_cols, pt)
        except ValueError as error:
            raise ValueError(f"Cannot aggregate {section}: {error}") from error

    # Microbiome resistance scalar
    vals = [p.get("microbiome_resistance", np.nan) for p in parsed_list]
    numeric = [v for v in vals if _is_numeric(v)]
    agg["microbiome_resistance"] = _fmt_agg(numeric) if numeric else "—"

    # Pass-through fields (use first run)
    agg["fit_stats"]            = parsed_list[0].get("fit_stats", {})
    agg["calibration_score"]    = parsed_list[0].get("calibration_score", {})
    agg["calibration_scores_all"] = [p.get("calibration_score", {}) for p in parsed_list]

    return agg
