"""Build the versioned long-form resistance target baseline.

The cleaned wide matrices remain the editable value sources. This builder
materializes their values, missing cells, provenance class, and score
eligibility in the long-form dataset consumed by calibration. Potency-informed
exclusions come from a Rust-generated projection; run-data filters remain
outside the schema.
"""

from __future__ import annotations

import csv
import re
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


TARGET_SET_VERSION = "resistance_targets_v1"
TARGET_YEAR = "2025"
POTENCY_CUTOFF = Decimal("0.15")

PREVALENCE_COMPONENT = "resistance_prevalence_any_r_positive"
SEVERITY_COMPONENT = "resistance_severity_conditional_mean_any_r"

TARGET_COLUMNS = [
    "target_set_version",
    "component",
    "bacteria",
    "drug",
    "value",
    "cell_status",
    "include_in_score",
    "score_exclusion_reason",
    "target_type",
    "source_id",
    "reference_year",
    "geography",
    "infection_or_specimen",
    "care_setting",
    "target_denominator",
    "uncertainty_lower",
    "uncertainty_upper",
    "transformation",
    "rationale",
    "evidence_weight",
    "score_row_weight",
]

SOURCE_COLUMNS = ["source_id", "source_type", "description", "url"]


def _source_slug(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "_", value.strip().lower()).strip("_")


def _target_bacteria_slug(value: str) -> str:
    return value.strip().lower().replace(" ", "_")


def _potency_bacteria_slug(value: str) -> str:
    slug = value.strip().lower()
    return "providencia_stuartii" if slug == "p_stuartii" else slug


def _read_wide_matrix(
    path: Path,
    *,
    expects_notes: bool,
) -> Tuple[List[str], List[str], Dict[str, Dict[str, str]], Dict[str, str]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None or not reader.fieldnames:
            raise ValueError(f"{path} has no header")
        if reader.fieldnames[0] != "Bacteria":
            raise ValueError(f"{path} must begin with Bacteria")

        metadata = {"notes"} if expects_notes else set()
        drug_names = [
            field for field in reader.fieldnames[1:] if field not in metadata
        ]
        expected_fields = ["Bacteria", *drug_names]
        if expects_notes:
            expected_fields.append("notes")
        if reader.fieldnames != expected_fields:
            raise ValueError(f"{path} has unexpected metadata columns or column order")

        bacteria_order: List[str] = []
        values: Dict[str, Dict[str, str]] = {}
        notes: Dict[str, str] = {}
        for row_number, row in enumerate(reader, start=2):
            bacterium = (row.get("Bacteria") or "").strip()
            if not bacterium:
                raise ValueError(f"{path}:{row_number} has an empty bacterium")
            if bacterium in values:
                raise ValueError(f"{path}:{row_number} duplicates {bacterium}")

            bacteria_order.append(bacterium)
            values[bacterium] = {}
            for drug in drug_names:
                token = (row.get(drug) or "").strip()
                if token == ".":
                    values[bacterium][drug] = token
                    continue
                try:
                    parsed = Decimal(token)
                except InvalidOperation as error:
                    raise ValueError(
                        f"{path}:{row_number} {bacterium}/{drug} has invalid value {token!r}"
                    ) from error
                if not parsed.is_finite() or not Decimal("0") <= parsed <= Decimal("1"):
                    raise ValueError(
                        f"{path}:{row_number} {bacterium}/{drug} is outside [0, 1]"
                    )
                values[bacterium][drug] = token

            if expects_notes:
                note = (row.get("notes") or "").strip()
                if not note:
                    raise ValueError(f"{path}:{row_number} has no provenance note")
                notes[bacterium] = note

    return bacteria_order, drug_names, values, notes


def _read_rust_potencies(path: Path) -> Dict[Tuple[str, str], Decimal]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"bacteria", "drug", "potency_when_no_r"}
        if reader.fieldnames is None or set(reader.fieldnames) != required:
            raise ValueError(f"{path} must contain exactly {sorted(required)}")

        result: Dict[Tuple[str, str], Decimal] = {}
        for row_number, row in enumerate(reader, start=2):
            bacteria = _potency_bacteria_slug(row["bacteria"])
            drug = row["drug"].strip().lower()
            key = (bacteria, drug)
            if key in result:
                raise ValueError(f"{path}:{row_number} duplicates {bacteria}/{drug}")
            result[key] = Decimal(row["potency_when_no_r"])
    return result


def _static_score_exclusions(
    bacterium: str,
    drug: str,
    prevalence_token: str,
    potencies: Mapping[Tuple[str, str], Decimal],
) -> List[str]:
    reasons: List[str] = []
    bacterium_slug = _target_bacteria_slug(bacterium)
    if prevalence_token == ".":
        reasons.append("legacy_prevalence_target_missing")
    if drug == "rifampicin":
        reasons.append("hard_exclusion_rifampicin")
    if "tuberculosis" in bacterium_slug:
        reasons.append("hard_exclusion_mdr_tb")
    if "listeria" in bacterium_slug:
        reasons.append("hard_exclusion_listeria")
    potency = potencies.get((bacterium_slug, drug))
    if potency is None:
        raise ValueError(f"Rust potency projection lacks {bacterium_slug}/{drug}")
    if potency < POTENCY_CUTOFF:
        reasons.append("model_baseline_potency_below_0.15")
    return reasons


def _base_row(
    *,
    component: str,
    bacterium: str,
    drug: str,
    value: str,
    status: str,
    exclusions: Sequence[str],
    target_type: str,
    source_id: str,
    denominator: str,
    transformation: str,
    rationale: str,
) -> Dict[str, str]:
    included = value != "." and not exclusions
    return {
        "target_set_version": TARGET_SET_VERSION,
        "component": component,
        "bacteria": bacterium,
        "drug": drug,
        "value": "" if value == "." else value,
        "cell_status": status,
        "include_in_score": str(included).lower(),
        "score_exclusion_reason": ";".join(exclusions),
        "target_type": target_type,
        "source_id": source_id,
        "reference_year": TARGET_YEAR if value != "." else "",
        "geography": "global_model_scope" if value != "." else "",
        "infection_or_specimen": "unspecified" if value != "." else "",
        "care_setting": "all_settings" if value != "." else "",
        "target_denominator": denominator if value != "." else "",
        "uncertainty_lower": "",
        "uncertainty_upper": "",
        "transformation": transformation if value != "." else "",
        "rationale": rationale if value != "." else "",
        "evidence_weight": "",
        "score_row_weight": "1.0" if included else "0.0",
    }


def _write_csv(path: Path, columns: Iterable[str], rows: Iterable[Mapping[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(columns), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def build_resistance_targets_v1(
    root: Optional[Path] = None,
    *,
    target_output: Optional[Path] = None,
    source_output: Optional[Path] = None,
) -> Tuple[Path, Path]:
    project_root = root or Path(__file__).resolve().parents[1]
    data_dir = project_root / "data"
    prevalence_path = data_dir / "resistance_prevalence_values.csv"
    severity_path = data_dir / "resistance_average_resistant_values.csv"
    potency_path = data_dir / "model_potency_matrix.csv"

    bacteria_order, drugs, prevalence, notes = _read_wide_matrix(
        prevalence_path, expects_notes=True
    )
    severity_order, severity_drugs, severity, _ = _read_wide_matrix(
        severity_path, expects_notes=False
    )
    if set(bacteria_order) != set(severity_order):
        raise ValueError("prevalence and severity matrices cover different bacteria")
    if drugs != severity_drugs:
        raise ValueError("prevalence and severity matrices cover different drugs or order")

    potencies = _read_rust_potencies(potency_path)
    rows: List[Dict[str, str]] = []

    for bacterium in bacteria_order:
        source_id = f"legacy_prevalence_note__{_source_slug(bacterium)}"
        for drug in drugs:
            token = prevalence[bacterium][drug]
            exclusions = _static_score_exclusions(bacterium, drug, token, potencies)
            rows.append(
                _base_row(
                    component=PREVALENCE_COMPONENT,
                    bacterium=bacterium,
                    drug=drug,
                    value=token,
                    status="active_target" if token != "." else "legacy_unclassified_missing",
                    exclusions=exclusions,
                    target_type=(
                        "evidence_informed_calibration_benchmark"
                        if token != "."
                        else "not_assigned"
                    ),
                    source_id=source_id,
                    denominator="source_definition_unrecovered",
                    transformation="legacy_cell_transformation_unrecovered",
                    rationale="legacy_bacterium_level_note",
                )
            )

    severity_source = "expert_model_scale_any_r_v1"
    for bacterium in bacteria_order:
        for drug in drugs:
            token = severity[bacterium][drug]
            prevalence_token = prevalence[bacterium][drug]
            exclusions = _static_score_exclusions(
                bacterium, drug, prevalence_token, potencies
            )
            if token == ".":
                status = "legacy_unclassified_missing"
                if "legacy_prevalence_target_missing" not in exclusions:
                    exclusions = [*exclusions, "severity_target_missing"]
            elif prevalence_token == ".":
                status = "inactive_legacy_prevalence_gate"
            else:
                status = "active_target"
            rows.append(
                _base_row(
                    component=SEVERITY_COMPONENT,
                    bacterium=bacterium,
                    drug=drug,
                    value=token,
                    status=status,
                    exclusions=exclusions,
                    target_type=(
                        "expert_assigned_model_benchmark"
                        if token != "."
                        else "not_assigned"
                    ),
                    source_id=severity_source if token != "." else "",
                    denominator="model_active_infection_person_days_with_any_r_positive",
                    transformation="expert_assignment_on_unitless_model_any_r_scale",
                    rationale="model_scale_resistance_severity_constraint",
                )
            )

    source_rows = [
        {
            "source_id": f"legacy_prevalence_note__{_source_slug(bacterium)}",
            "source_type": "legacy_bacterium_level_note",
            "description": notes[bacterium],
            "url": "",
        }
        for bacterium in bacteria_order
    ]
    source_rows.append(
        {
            "source_id": severity_source,
            "source_type": "expert_model_design",
            "description": (
                "Expert-assigned model benchmarks for mean any_r "
                "conditional on any_r > 0; these are not direct clinical surveillance estimates."
            ),
            "url": "",
        }
    )

    resolved_target_output = target_output or data_dir / "resistance_targets_v1.csv"
    resolved_source_output = source_output or data_dir / "resistance_target_sources_v1.csv"
    _write_csv(resolved_target_output, TARGET_COLUMNS, rows)
    _write_csv(resolved_source_output, SOURCE_COLUMNS, source_rows)
    return resolved_target_output, resolved_source_output


if __name__ == "__main__":
    targets_path, sources_path = build_resistance_targets_v1()
    print(f"Wrote {targets_path}")
    print(f"Wrote {sources_path}")
