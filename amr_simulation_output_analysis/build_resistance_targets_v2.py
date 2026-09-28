"""Build reviewed v2 references from an immutable v1 snapshot and explicit decisions.

The legacy live wide matrices are historical inputs, never a fallback for v2.
Only the current Rust projections are refreshed; their eligibility rules are
unchanged. The review file pins the original snapshot manifest and guards every
corrected value/source so newer evidence cannot be silently overwritten.
"""

from __future__ import annotations

import csv
import json
from decimal import Decimal
from pathlib import Path
from typing import Optional, Tuple

from .build_resistance_targets_v1 import (
    PREVALENCE_COMPONENT,
    SEVERITY_COMPONENT,
    SOURCE_COLUMNS,
    TARGET_COLUMNS,
    _read_rust_potencies,
    _read_rust_resistance_reachability,
    _sha256,
    _static_score_exclusions,
    _target_bacteria_slug,
    _write_csv,
)


TARGET_SET_VERSION = "resistance_targets_v2"
MANIFEST_FILENAME = f"{TARGET_SET_VERSION}.manifest.json"
UPDATE_LOCK_FILENAME = f".{TARGET_SET_VERSION}.update.lock"
SNAPSHOT_DIRECTORY = "resistance_targets_v1_snapshot"
REVIEW_FILENAME = f"{TARGET_SET_VERSION}.review.json"
REVISIONS_FILENAME = f"{TARGET_SET_VERSION}.revisions.csv"
SCHEMA_FILENAME = f"{TARGET_SET_VERSION}.schema.json"
SNAPSHOT_FILENAMES = tuple(
    f"{SNAPSHOT_DIRECTORY}/{name}"
    for name in (
        "resistance_targets_v1.csv",
        "resistance_target_sources_v1.csv",
        "resistance_targets_v1.schema.json",
        "resistance_targets_v1.manifest.json",
        "resistance_prevalence_values.csv",
        "resistance_average_resistant_values.csv",
        "model_potency_matrix.csv",
        "model_resistance_reachability_matrix.csv",
    )
)
INPUT_FILENAMES = (*SNAPSHOT_FILENAMES, REVIEW_FILENAME, SCHEMA_FILENAME)
GENERATED_FILENAMES = (
    "model_potency_matrix.csv",
    "model_resistance_reachability_matrix.csv",
    f"{TARGET_SET_VERSION}.csv",
    "resistance_target_sources_v2.csv",
    REVISIONS_FILENAME,
    MANIFEST_FILENAME,
)
MANIFEST_ARTIFACTS = (*INPUT_FILENAMES, *GENERATED_FILENAMES[:-1])
AUDIT_COLUMNS = [
    "old_target_set_version", "new_target_set_version", "component", "bacteria",
    "drug", "revision_kind", "old_value", "new_value", "old_status", "new_status",
    "old_provenance_class", "new_provenance_class", "old_source_id", "new_source_id",
    "old_include_in_score", "new_include_in_score", "old_score_exclusion_reason",
    "new_score_exclusion_reason", "rationale", "review_source_ids",
]


def _read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _verify_snapshot(data_dir: Path, review: dict) -> None:
    snapshot = data_dir / SNAPSHOT_DIRECTORY
    manifest_path = snapshot / "resistance_targets_v1.manifest.json"
    if _sha256(manifest_path) != review["baseline_manifest_sha256"]:
        raise ValueError("The v1 snapshot manifest differs from the reviewed baseline")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected = {Path(name).name for name in SNAPSHOT_FILENAMES}
    expected.remove(manifest_path.name)
    if set(manifest["artifacts"]) != expected:
        raise ValueError("The v1 snapshot manifest has an unexpected artifact set")
    for name, recorded in manifest["artifacts"].items():
        path = snapshot / name
        if path.stat().st_size != recorded["bytes"] or _sha256(path) != recorded["sha256"]:
            raise ValueError(f"The v1 snapshot is not intact: {name}")


def _refresh_eligibility(rows: list[dict[str, str]], data_dir: Path) -> None:
    potencies = _read_rust_potencies(data_dir / "model_potency_matrix.csv")
    reachability = _read_rust_resistance_reachability(
        data_dir / "model_resistance_reachability_matrix.csv"
    )
    prevalence = {
        (row["bacteria"], row["drug"]): row["value"] or "."
        for row in rows if row["component"] == PREVALENCE_COMPONENT
    }
    for row in rows:
        pair = (row["bacteria"], row["drug"])
        token = prevalence[pair]
        severity = row["component"] == SEVERITY_COMPONENT
        reasons = _static_score_exclusions(
            *pair, token, potencies, reachability, exclude_unrepresentable=severity
        )
        representable, maximum = reachability[(_target_bacteria_slug(pair[0]), pair[1])]
        if severity and row["value"] and Decimal(row["value"]) > maximum + Decimal("1e-12"):
            reasons.append("severity_benchmark_above_model_representable_maximum")
        if not row["value"]:
            status = "legacy_unclassified_missing"
            if severity and "legacy_prevalence_target_missing" not in reasons:
                reasons.append("severity_target_missing")
        elif severity and token == ".":
            status = "inactive_unpaired_legacy_benchmark"
        elif severity and "model_resistance_phenotype_not_representable" in reasons:
            status = "inactive_model_unrepresentable"
        elif severity and "severity_benchmark_above_model_representable_maximum" in reasons:
            status = "inactive_above_model_representable_maximum"
        elif not severity and not representable:
            status = "active_target_model_unrepresentable"
        else:
            status = "active_target"
        included = bool(row["value"]) and not reasons
        row.update(
            cell_status=status, include_in_score=str(included).lower(),
            score_exclusion_reason=";".join(reasons),
            score_row_weight="1.0" if included else "0.0",
        )


def build_resistance_targets_v2(
    root: Optional[Path] = None,
    *,
    target_output: Optional[Path] = None,
    source_output: Optional[Path] = None,
    manifest_output: Optional[Path] = None,
    revisions_output: Optional[Path] = None,
) -> Tuple[Path, Path]:
    project_root = root or Path(__file__).resolve().parents[1]
    data_dir = project_root / "data"
    review = json.loads((data_dir / REVIEW_FILENAME).read_text(encoding="utf-8"))
    if review["target_set_version"] != TARGET_SET_VERSION:
        raise ValueError("Review targets an unexpected version")
    _verify_snapshot(data_dir, review)
    baseline = _read_csv(data_dir / SNAPSHOT_DIRECTORY / "resistance_targets_v1.csv")
    rows = [dict(row, target_set_version=TARGET_SET_VERSION) for row in baseline]
    by_key = {(r["component"], r["bacteria"], r["drug"]): r for r in rows}
    if len(by_key) != len(rows):
        raise ValueError("Duplicate baseline target keys")
    decisions = {}
    for decision in review["decisions"]:
        key = (decision["component"], decision["bacteria"], decision["drug"])
        if key in decisions:
            raise ValueError(f"Duplicate review decision: {key}")
        row = by_key[key]
        for field, expected in decision["expected"].items():
            if row[field] != expected:
                raise ValueError(f"Review baseline differs for {key}/{field}")
        updates = decision["updates"]
        if not set(updates).issubset(TARGET_COLUMNS) or "target_set_version" in updates:
            raise ValueError(f"Invalid review fields for {key}")
        row.update(updates)
        decisions[key] = decision

    _refresh_eligibility(rows, data_dir)
    sources = _read_csv(data_dir / SNAPSHOT_DIRECTORY / "resistance_target_sources_v1.csv")
    sources.extend(review["source_registry_additions"])
    source_ids = [row["source_id"] for row in sources]
    if len(source_ids) != len(set(source_ids)):
        raise ValueError("Duplicate source IDs")
    audits = []
    for old, new in zip(baseline, rows):
        if all(old[field] == new[field] for field in TARGET_COLUMNS if field != "target_set_version"):
            continue
        key = (new["component"], new["bacteria"], new["drug"])
        decision = decisions.get(key)
        paired = decisions.get((PREVALENCE_COMPONENT, new["bacteria"], new["drug"]))
        if decision is not None:
            kind = decision["revision_kind"]
            rationale = decision["rationale"]
            citations = decision["review_source_ids"]
        elif paired and new["cell_status"] == "inactive_unpaired_legacy_benchmark":
            kind = "severity_pairing_exclusion_value_preserved"
            rationale = "Numerical severity retained; paired prevalence was withdrawn as unassigned."
            citations = paired["review_source_ids"]
        else:
            kind = "eligibility_projection_refresh"
            rationale = "Existing eligibility rules applied to the current Rust projections."
            citations = []
        audits.append({
            "old_target_set_version": old["target_set_version"],
            "new_target_set_version": TARGET_SET_VERSION,
            "component": new["component"], "bacteria": new["bacteria"], "drug": new["drug"],
            "revision_kind": kind, "old_value": old["value"], "new_value": new["value"],
            "old_status": old["cell_status"], "new_status": new["cell_status"],
            "old_provenance_class": old["provenance_class"],
            "new_provenance_class": new["provenance_class"],
            "old_source_id": old["source_id"], "new_source_id": new["source_id"],
            "old_include_in_score": old["include_in_score"],
            "new_include_in_score": new["include_in_score"],
            "old_score_exclusion_reason": old["score_exclusion_reason"],
            "new_score_exclusion_reason": new["score_exclusion_reason"],
            "rationale": rationale, "review_source_ids": ";".join(citations),
        })

    targets = target_output or data_dir / f"{TARGET_SET_VERSION}.csv"
    source_path = source_output or data_dir / "resistance_target_sources_v2.csv"
    audit_path = revisions_output or targets.parent / REVISIONS_FILENAME
    manifest_path = manifest_output or targets.parent / MANIFEST_FILENAME
    _write_csv(targets, TARGET_COLUMNS, rows)
    _write_csv(source_path, SOURCE_COLUMNS, sources)
    _write_csv(audit_path, AUDIT_COLUMNS, audits)
    artifact_paths = {name: data_dir / name for name in MANIFEST_ARTIFACTS}
    artifact_paths.update({
        f"{TARGET_SET_VERSION}.csv": targets,
        "resistance_target_sources_v2.csv": source_path,
        REVISIONS_FILENAME: audit_path,
    })
    manifest = {
        "target_set_version": TARGET_SET_VERSION,
        "hash_algorithm": "sha256",
        "baseline_target_set_version": "resistance_targets_v1",
        "review_date": review["review_date"],
        "artifacts": {
            name: {"sha256": _sha256(path), "bytes": path.stat().st_size}
            for name, path in artifact_paths.items()
        },
    }
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return targets, source_path


if __name__ == "__main__":
    targets, sources = build_resistance_targets_v2()
    print(f"Wrote {targets}")
    print(f"Wrote {sources}")
    print(f"Wrote {targets.parent / REVISIONS_FILENAME}")
    print(f"Wrote {targets.parent / MANIFEST_FILENAME}")
