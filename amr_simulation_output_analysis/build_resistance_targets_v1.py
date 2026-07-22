"""Build the versioned long-form resistance target baseline.

The cleaned wide matrices remain the editable value sources. This builder
materializes their values, missing cells, provenance class, and score
eligibility in the long-form dataset consumed by calibration. Potency-informed
exclusions come from a Rust-generated projection; run-data filters remain
outside the schema.
"""

from __future__ import annotations

import csv
import hashlib
import json
import re
from decimal import Decimal, InvalidOperation
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple


TARGET_SET_VERSION = "resistance_targets_v1"
TARGET_YEAR = "2025"
POTENCY_CUTOFF = Decimal("0.15")

PREVALENCE_COMPONENT = "resistance_prevalence_any_r_positive"
SEVERITY_COMPONENT = "resistance_severity_conditional_mean_any_r"

DEFAULT_SEVERITY_SOURCE = "expert_model_scale_any_r_v1"
RESERVE_DRUG_SEVERITY_SOURCE = "expert_reserve_drug_any_r_placeholders_v1"
RARE_POSITIVE_SEVERITY_SOURCE = "expert_rare_positive_any_r_prior_v1"

PROVENANCE_EMPIRICAL = "empirical_estimate_with_cell_level_source"
PROVENANCE_EVIDENCE_UNRESOLVED = (
    "evidence_informed_benchmark_cell_provenance_unrecovered"
)
PROVENANCE_EXPERT_PLACEHOLDER = "expert_informed_placeholder"
PROVENANCE_STRUCTURAL_PRIOR = "structural_prior"
PROVENANCE_NOT_ASSIGNED = "not_assigned"

MANIFEST_FILENAME = "resistance_targets_v1.manifest.json"

_SHARED_RESERVE_DRUG_PLACEHOLDER_BACTERIA = frozenset(
    {
        "Citrobacter spp.",
        "Enterobacter spp.",
        "Escherichia coli",
        "Klebsiella pneumoniae",
        "Morganella spp.",
        "Proteus spp.",
        "Serratia spp.",
        "Pseudomonas aeruginosa",
        "Salmonella enterica serovar typhi",
        "Salmonella enterica serovar paratyphi a",
        "Invasive non-typhoidal Salmonella spp.",
        "Shigella spp.",
        "Neisseria gonorrhoeae",
        "Haemophilus influenzae",
        "Vibrio cholerae",
        "Neisseria meningitidis",
        "Campylobacter jejuni",
        "Enterobacter cloacae",
        "Yersinia enterocolitica",
        "Moraxella catarrhalis",
        "Bordetella pertussis",
        "Bacteroides fragilis",
        "Providencia stuartii",
    }
)
RESERVE_DRUG_SEVERITY_PLACEHOLDER_PAIRS = frozenset(
    {
        *(
            (bacterium, "cefiderocol")
            for bacterium in _SHARED_RESERVE_DRUG_PLACEHOLDER_BACTERIA
        ),
        ("Acinetobacter baumannii", "cefiderocol"),
        *(
            (bacterium, "ceftolozane_tazobactam")
            for bacterium in _SHARED_RESERVE_DRUG_PLACEHOLDER_BACTERIA
        ),
        ("Legionella pneumophila", "ceftolozane_tazobactam"),
    }
)
RARE_POSITIVE_SEVERITY_PRIOR_PAIRS = frozenset(
    {
        ("Staphylococcus aureus", "vancomycin"),
        ("Staphylococcus epidermidis", "vancomycin"),
        ("Streptococcus pneumoniae", "vancomycin"),
        ("Streptococcus pneumoniae", "linezolid"),
        ("Streptococcus pneumoniae", "daptomycin"),
    }
)

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
    "provenance_class",
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

SOURCE_COLUMNS = [
    "source_id",
    "provenance_class",
    "source_type",
    "description",
    "url",
]


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


def _read_rust_resistance_reachability(
    path: Path,
) -> Dict[Tuple[str, str], Tuple[bool, Decimal]]:
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {
            "bacteria",
            "drug",
            "resistance_representable",
            "maximum_any_r",
            "positive_effect_mechanisms",
        }
        if reader.fieldnames is None or set(reader.fieldnames) != required:
            raise ValueError(f"{path} must contain exactly {sorted(required)}")

        result: Dict[Tuple[str, str], Tuple[bool, Decimal]] = {}
        for row_number, row in enumerate(reader, start=2):
            bacteria = _potency_bacteria_slug(row["bacteria"])
            drug = row["drug"].strip().lower()
            key = (bacteria, drug)
            if key in result:
                raise ValueError(f"{path}:{row_number} duplicates {bacteria}/{drug}")
            token = row["resistance_representable"].strip().lower()
            if token not in {"true", "false"}:
                raise ValueError(
                    f"{path}:{row_number} has invalid representability {token!r}"
                )
            representable = token == "true"
            try:
                maximum_any_r = Decimal(row["maximum_any_r"])
            except InvalidOperation as error:
                raise ValueError(
                    f"{path}:{row_number} has invalid maximum_any_r"
                ) from error
            if not maximum_any_r.is_finite() or not Decimal("0") <= maximum_any_r <= Decimal("1"):
                raise ValueError(f"{path}:{row_number} maximum_any_r is outside [0, 1]")
            has_mechanisms = bool(row["positive_effect_mechanisms"].strip())
            if representable != has_mechanisms:
                raise ValueError(
                    f"{path}:{row_number} representability disagrees with mechanism list"
                )
            if representable != (maximum_any_r > 0):
                raise ValueError(
                    f"{path}:{row_number} representability disagrees with maximum_any_r"
                )
            result[key] = (representable, maximum_any_r)
    return result


def _static_score_exclusions(
    bacterium: str,
    drug: str,
    prevalence_token: str,
    potencies: Mapping[Tuple[str, str], Decimal],
    resistance_reachability: Mapping[Tuple[str, str], Tuple[bool, Decimal]],
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
    reachability = resistance_reachability.get((bacterium_slug, drug))
    if reachability is None:
        raise ValueError(
            f"Rust resistance reachability projection lacks {bacterium_slug}/{drug}"
        )
    if prevalence_token != "." and not reachability[0]:
        reasons.append("model_resistance_phenotype_not_representable")
    return reasons


def _severity_provenance(bacterium: str, drug: str) -> Tuple[str, str, str]:
    pair = (bacterium, drug)
    if pair in RESERVE_DRUG_SEVERITY_PLACEHOLDER_PAIRS:
        return (
            RESERVE_DRUG_SEVERITY_SOURCE,
            "expert_best_guess_reserve_drug_any_r_placeholder",
            PROVENANCE_EXPERT_PLACEHOLDER,
        )
    if pair in RARE_POSITIVE_SEVERITY_PRIOR_PAIRS:
        return (
            RARE_POSITIVE_SEVERITY_SOURCE,
            "expert_rare_positive_any_r_structural_prior",
            PROVENANCE_STRUCTURAL_PRIOR,
        )
    return (
        DEFAULT_SEVERITY_SOURCE,
        "model_scale_resistance_severity_constraint",
        PROVENANCE_EXPERT_PLACEHOLDER,
    )


def _base_row(
    *,
    component: str,
    bacterium: str,
    drug: str,
    value: str,
    status: str,
    exclusions: Sequence[str],
    target_type: str,
    provenance_class: str,
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
        "provenance_class": provenance_class,
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


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_resistance_target_manifest(
    data_dir: Path,
    target_path: Path,
    source_path: Path,
    manifest_path: Path,
) -> None:
    artifact_paths = {
        "resistance_targets_v1.csv": target_path,
        "resistance_target_sources_v1.csv": source_path,
        "resistance_targets_v1.schema.json": data_dir
        / "resistance_targets_v1.schema.json",
        "resistance_prevalence_values.csv": data_dir
        / "resistance_prevalence_values.csv",
        "resistance_average_resistant_values.csv": data_dir
        / "resistance_average_resistant_values.csv",
        "model_potency_matrix.csv": data_dir / "model_potency_matrix.csv",
        "model_resistance_reachability_matrix.csv": data_dir
        / "model_resistance_reachability_matrix.csv",
    }
    missing = [name for name, path in artifact_paths.items() if not path.exists()]
    if missing:
        raise FileNotFoundError(
            "Cannot build resistance-target manifest; missing artifacts: "
            + ", ".join(missing)
        )

    payload = {
        "target_set_version": TARGET_SET_VERSION,
        "hash_algorithm": "sha256",
        "artifacts": {
            name: {
                "sha256": _sha256(path),
                "bytes": path.stat().st_size,
            }
            for name, path in artifact_paths.items()
        },
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def build_resistance_targets_v1(
    root: Optional[Path] = None,
    *,
    target_output: Optional[Path] = None,
    source_output: Optional[Path] = None,
    manifest_output: Optional[Path] = None,
) -> Tuple[Path, Path]:
    project_root = root or Path(__file__).resolve().parents[1]
    data_dir = project_root / "data"
    prevalence_path = data_dir / "resistance_prevalence_values.csv"
    severity_path = data_dir / "resistance_average_resistant_values.csv"
    potency_path = data_dir / "model_potency_matrix.csv"
    reachability_path = data_dir / "model_resistance_reachability_matrix.csv"

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
    resistance_reachability = _read_rust_resistance_reachability(reachability_path)
    rows: List[Dict[str, str]] = []

    for bacterium in bacteria_order:
        source_id = f"legacy_prevalence_note__{_source_slug(bacterium)}"
        for drug in drugs:
            token = prevalence[bacterium][drug]
            exclusions = _static_score_exclusions(
                bacterium,
                drug,
                token,
                potencies,
                resistance_reachability,
            )
            status = "active_target"
            if token == ".":
                status = "legacy_unclassified_missing"
            elif "model_resistance_phenotype_not_representable" in exclusions:
                status = "inactive_model_unrepresentable"
            rows.append(
                _base_row(
                    component=PREVALENCE_COMPONENT,
                    bacterium=bacterium,
                    drug=drug,
                    value=token,
                    status=status,
                    exclusions=exclusions,
                    target_type=(
                        "evidence_informed_calibration_benchmark"
                        if token != "."
                        else "not_assigned"
                    ),
                    provenance_class=(
                        PROVENANCE_EVIDENCE_UNRESOLVED
                        if token != "."
                        else PROVENANCE_NOT_ASSIGNED
                    ),
                    source_id=source_id,
                    denominator="source_definition_unrecovered",
                    transformation="legacy_cell_transformation_unrecovered",
                    rationale="legacy_bacterium_level_note",
                )
            )

    for bacterium in bacteria_order:
        for drug in drugs:
            token = severity[bacterium][drug]
            prevalence_token = prevalence[bacterium][drug]
            (
                severity_source,
                severity_rationale,
                severity_provenance_class,
            ) = _severity_provenance(bacterium, drug)
            exclusions = _static_score_exclusions(
                bacterium,
                drug,
                prevalence_token,
                potencies,
                resistance_reachability,
            )
            if token != ".":
                maximum_any_r = resistance_reachability[
                    (_target_bacteria_slug(bacterium), drug)
                ][1]
                if Decimal(token) > maximum_any_r + Decimal("1e-12"):
                    exclusions = [
                        *exclusions,
                        "severity_benchmark_above_model_representable_maximum",
                    ]
            if token == ".":
                status = "legacy_unclassified_missing"
                if "legacy_prevalence_target_missing" not in exclusions:
                    exclusions = [*exclusions, "severity_target_missing"]
            elif prevalence_token == ".":
                status = "inactive_unpaired_legacy_benchmark"
            elif "model_resistance_phenotype_not_representable" in exclusions:
                status = "inactive_model_unrepresentable"
            elif "severity_benchmark_above_model_representable_maximum" in exclusions:
                status = "inactive_above_model_representable_maximum"
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
                    provenance_class=(
                        severity_provenance_class
                        if token != "."
                        else PROVENANCE_NOT_ASSIGNED
                    ),
                    source_id=severity_source if token != "." else "",
                    denominator="model_active_infection_person_days_with_any_r_positive",
                    transformation="expert_assignment_on_unitless_model_any_r_scale",
                    rationale=severity_rationale,
                )
            )

    source_rows = [
        {
            "source_id": f"legacy_prevalence_note__{_source_slug(bacterium)}",
            "provenance_class": PROVENANCE_EVIDENCE_UNRESOLVED,
            "source_type": "legacy_bacterium_level_note",
            "description": notes[bacterium],
            "url": "",
        }
        for bacterium in bacteria_order
    ]
    source_rows.append(
        {
            "source_id": DEFAULT_SEVERITY_SOURCE,
            "provenance_class": PROVENANCE_EXPERT_PLACEHOLDER,
            "source_type": "expert_placeholder",
            "description": (
                "Expert-assigned model benchmarks for mean any_r "
                "conditional on any_r > 0; these are not direct clinical surveillance estimates."
            ),
            "url": "",
        }
    )
    source_rows.append(
        {
            "source_id": RESERVE_DRUG_SEVERITY_SOURCE,
            "provenance_class": PROVENANCE_EXPERT_PLACEHOLDER,
            "source_type": "expert_placeholder",
            "description": (
                "Coarse expert best-guess placeholders for mean any_r conditional on "
                "any_r > 0: 0.60 for cefiderocol and 0.70 for "
                "ceftolozane/tazobactam. They replace legacy cells that duplicated "
                "resistance-prevalence benchmarks and are not empirical estimates."
            ),
            "url": "",
        }
    )
    source_rows.append(
        {
            "source_id": RARE_POSITIVE_SEVERITY_SOURCE,
            "provenance_class": PROVENANCE_STRUCTURAL_PRIOR,
            "source_type": "structural_prior",
            "description": (
                "Expert best-guess structural priors for mean any_r conditional on "
                "any_r > 0 in five bacterium-drug pairs with a zero prevalence "
                "benchmark. They define severity if rare simulated positives occur "
                "and are not empirical prevalence estimates."
            ),
            "url": "",
        }
    )

    resolved_target_output = target_output or data_dir / "resistance_targets_v1.csv"
    resolved_source_output = source_output or data_dir / "resistance_target_sources_v1.csv"
    resolved_manifest_output = (
        manifest_output or resolved_target_output.parent / MANIFEST_FILENAME
    )
    _write_csv(resolved_target_output, TARGET_COLUMNS, rows)
    _write_csv(resolved_source_output, SOURCE_COLUMNS, source_rows)
    write_resistance_target_manifest(
        data_dir,
        resolved_target_output,
        resolved_source_output,
        resolved_manifest_output,
    )
    return resolved_target_output, resolved_source_output


if __name__ == "__main__":
    targets_path, sources_path = build_resistance_targets_v1()
    print(f"Wrote {targets_path}")
    print(f"Wrote {sources_path}")
    print(f"Wrote {targets_path.parent / MANIFEST_FILENAME}")
