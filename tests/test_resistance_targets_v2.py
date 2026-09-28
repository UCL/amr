"""Regression contract for the narrowly reviewed resistance reference revision."""

import csv
import hashlib
import json
import shutil
import unittest
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from amr_simulation_output_analysis.build_resistance_targets_v1 import (
    PREVALENCE_COMPONENT as PREVALENCE,
    SEVERITY_COMPONENT as SEVERITY,
    build_resistance_targets_v1,
)
from amr_simulation_output_analysis.build_resistance_targets_v2 import (
    GENERATED_FILENAMES,
    INPUT_FILENAMES,
    MANIFEST_ARTIFACTS,
    MANIFEST_FILENAME,
    SNAPSHOT_DIRECTORY,
    build_resistance_targets_v2,
)
from amr_simulation_output_analysis.calibration_summary import (
    _load_resistance_target_set,
    _verify_resistance_target_manifest,
)


DATA = Path(__file__).resolve().parents[1] / "data"
TP_DRUGS = {"tetracycline": "0.08", "doxycycline": "0.06", "minocycline": "0.05"}
CT_DRUGS = {
    "erythromycin": "0.1", "azithromycin": "0.05", "clarithromycin": "0.1",
    "clindamycin": "0.2", "ciprofloxacin": "0.1", "levofloxacin": "0.05",
    "moxifloxacin": "0.05", "ofloxacin": "0.1", "tetracycline": "0.05",
    "doxycycline": "0.05", "minocycline": "0.05", "tigecycline": "0.05",
    "chloramphenicol": "0.1",
}
WITHDRAWN = {
    *(('Chlamydia trachomatis', drug) for drug in CT_DRUGS),
    ("Enterococcus faecium", "tigecycline"),
    ("Enterococcus faecium", "linezolid"),
}


def read_rows(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def keyed(rows):
    return {(r["component"], r["bacteria"], r["drug"]): r for r in rows}


def copy_inputs(temp_root):
    data = temp_root / "data"
    for name in (*INPUT_FILENAMES, *GENERATED_FILENAMES[:2]):
        destination = data / name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(DATA / name, destination)
    return data


class ResistanceTargetsV2Tests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.rows = read_rows(DATA / "resistance_targets_v2.csv")
        cls.targets = keyed(cls.rows)
        cls.baseline_rows = read_rows(DATA / SNAPSHOT_DIRECTORY / "resistance_targets_v1.csv")
        cls.baseline = keyed(cls.baseline_rows)
        cls.audit = read_rows(DATA / "resistance_targets_v2.revisions.csv")
        cls.review = json.loads((DATA / "resistance_targets_v2.review.json").read_text())

    def test_three_zero_prevalence_priors_are_assigned_and_severity_is_unchanged(self):
        for drug, original in TP_DRUGS.items():
            with self.subTest(drug=drug):
                key = (PREVALENCE, "Treponema pallidum", drug)
                row = self.targets[key]
                self.assertEqual(self.baseline[key]["value"], original)
                self.assertEqual(row["value"], "0.0")
                self.assertEqual(row["provenance_class"], "structural_prior")
                self.assertEqual(row["target_type"], "expert_assigned_model_benchmark")
                self.assertEqual(row["include_in_score"], self.baseline[key]["include_in_score"])
                self.assertEqual(row["uncertainty_lower"], "")
                self.assertEqual(row["uncertainty_upper"], "")
                severity_key = (SEVERITY, "Treponema pallidum", drug)
                severity = self.targets[severity_key]
                self.assertEqual(severity["value"], self.baseline[severity_key]["value"])
                self.assertEqual(severity["provenance_class"], "structural_prior")
                self.assertIn("rare-positive", severity["rationale"])
                self.assertEqual(severity["include_in_score"], self.baseline[severity_key]["include_in_score"])
        self.assertIn("Class-based", self.targets[(PREVALENCE, "Treponema pallidum", "minocycline")]["rationale"])

    def test_exact_withdrawal_set_and_severity_pairing(self):
        audited = {
            (r["bacteria"], r["drug"]): r for r in self.audit
            if r["revision_kind"] == "withdrawal_as_unassigned"
        }
        self.assertEqual(set(audited), WITHDRAWN)
        for bacterium, drug in WITHDRAWN:
            with self.subTest(bacteria=bacterium, drug=drug):
                row = self.targets[(PREVALENCE, bacterium, drug)]
                self.assertEqual(row["value"], "")
                self.assertEqual(row["target_type"], "not_assigned")
                self.assertEqual(row["provenance_class"], "not_assigned")
                self.assertEqual(row["include_in_score"], "false")
                self.assertEqual(row["score_row_weight"], "0.0")
                key = (SEVERITY, bacterium, drug)
                severity = self.targets[key]
                self.assertEqual(severity["value"], self.baseline[key]["value"])
                self.assertEqual(severity["cell_status"], "inactive_unpaired_legacy_benchmark")
                self.assertEqual(severity["include_in_score"], "false")
                self.assertIn("legacy_prevalence_target_missing", severity["score_exclusion_reason"])
                self.assertEqual(audited[(bacterium, drug)]["old_value"], self.baseline[(PREVALENCE, bacterium, drug)]["value"])
        for drug, value in CT_DRUGS.items():
            self.assertEqual(audited[("Chlamydia trachomatis", drug)]["old_value"], value)

    def test_only_authorized_rows_change_and_every_change_has_history(self):
        allowed = {
            (component, bacterium, drug)
            for bacterium, drug in WITHDRAWN | {("Treponema pallidum", d) for d in TP_DRUGS}
            for component in (PREVALENCE, SEVERITY)
        }
        actual = set()
        for key, old in self.baseline.items():
            new = self.targets[key]
            changed = {field for field in old if field != "target_set_version" and old[field] != new[field]}
            if changed:
                actual.add(key)
            if key[0] == SEVERITY:
                self.assertEqual(old["value"], new["value"])
        self.assertEqual(actual, allowed)
        self.assertEqual(set(keyed(self.audit)), actual)
        self.assertEqual(len(self.audit), 36)
        for drug, expected in {"tetracycline": "0.5", "tedizolid": "0.1", "ampicillin": "", "amoxicillin": ""}.items():
            self.assertEqual(self.targets[(PREVALENCE, "Enterococcus faecium", drug)]["value"], expected)
        self.assertIn("tedizolid", self.review["unresolved"])

    def test_schema_units_unique_keys_and_source_consistency(self):
        schema = json.loads((DATA / "resistance_targets_v2.schema.json").read_text())
        sources = {r["source_id"]: r for r in read_rows(DATA / "resistance_target_sources_v2.csv")}
        self.assertEqual(len(self.targets), len(self.rows))
        self.assertEqual(set(self.targets), set(self.baseline))
        for row in self.rows:
            self.assertEqual(set(row), set(schema["required"]))
            for field, rule in schema["properties"].items():
                self.assertIsInstance(row[field], str)
                if "enum" in rule:
                    self.assertIn(row[field], rule["enum"])
                if "const" in rule:
                    self.assertEqual(row[field], rule["const"])
                if "minLength" in rule:
                    self.assertGreaterEqual(len(row[field]), rule["minLength"])
            if row["value"]:
                self.assertTrue(0 <= float(row["value"]) <= 1)
                self.assertEqual(row["provenance_class"], sources[row["source_id"]]["provenance_class"])
            self.assertEqual(row["evidence_weight"], "")
        for pair in WITHDRAWN:
            row = self.targets[(PREVALENCE, *pair)]
            self.assertEqual(sources[row["source_id"]]["provenance_class"], "not_assigned")
        old_tigecycline = self.baseline[(PREVALENCE, "Enterococcus faecium", "tigecycline")]
        self.assertEqual(float(old_tigecycline["value"]) * 100, 40.0)
        self.assertEqual(Counter(r["component"] for r in self.rows if r["include_in_score"] == "true"),
                         {PREVALENCE: 1161, SEVERITY: 1043})

    def test_loader_preserves_explicit_zero_and_unassigned(self):
        prevalence, severity = _load_resistance_target_set(DATA / "resistance_targets_v2.csv")
        prevalence = prevalence.set_index(["Bacteria", "drug"])
        severity = severity.set_index(["Bacteria", "drug"])
        for drug in TP_DRUGS:
            self.assertEqual(prevalence.loc[("Treponema pallidum", drug), "target"], 0.0)
        for pair in WITHDRAWN | {("Enterococcus faecium", "ampicillin"), ("Enterococcus faecium", "amoxicillin")}:
            self.assertTrue(pd.isna(prevalence.loc[pair, "target"]))
            self.assertFalse(prevalence.loc[pair, "include_in_score"])
        for pair in WITHDRAWN:
            self.assertFalse(severity.loc[pair, "include_in_score"])
            self.assertEqual(severity.loc[pair, "target"], float(self.baseline[(SEVERITY, *pair)]["value"]))

    def test_manifest_covers_all_version_bound_inputs_and_detects_corruption(self):
        manifest = json.loads((DATA / MANIFEST_FILENAME).read_text())
        self.assertEqual(set(manifest["artifacts"]), set(MANIFEST_ARTIFACTS))
        _verify_resistance_target_manifest(DATA / "resistance_targets_v2.csv")
        for name, recorded in manifest["artifacts"].items():
            payload = (DATA / name).read_bytes()
            self.assertEqual(recorded["bytes"], len(payload))
            self.assertEqual(recorded["sha256"], hashlib.sha256(payload).hexdigest())
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            data = copy_inputs(root)
            build_resistance_targets_v2(root)
            with (data / "resistance_targets_v2.review.json").open("a") as handle:
                handle.write("\n")
            with self.assertRaisesRegex(ValueError, "mismatch"):
                _verify_resistance_target_manifest(data / "resistance_targets_v2.csv")

    def test_generation_is_reproducible_and_live_legacy_matrices_cannot_restore_targets(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            data = copy_inputs(root)
            for name in ("resistance_prevalence_values.csv", "resistance_average_resistant_values.csv"):
                (data / name).write_text("INVALID: live legacy inputs must not be read")
            build_resistance_targets_v2(root)
            for name in GENERATED_FILENAMES[2:]:
                self.assertEqual((data / name).read_bytes(), (DATA / name).read_bytes(), name)
            baseline_path = data / SNAPSHOT_DIRECTORY / "resistance_targets_v1.csv"
            with baseline_path.open("a") as handle:
                handle.write("\n")
            with self.assertRaisesRegex(ValueError, "snapshot is not intact"):
                build_resistance_targets_v2(root)

    def test_preserved_v1_snapshot_remains_independently_reproducible(self):
        with TemporaryDirectory() as temporary:
            root = Path(temporary)
            shutil.copytree(DATA / SNAPSHOT_DIRECTORY, root / "data")
            build_resistance_targets_v1(root)
            for name in ("resistance_targets_v1.csv", "resistance_target_sources_v1.csv", "resistance_targets_v1.manifest.json"):
                self.assertEqual((root / "data" / name).read_bytes(), (DATA / SNAPSHOT_DIRECTORY / name).read_bytes(), name)


if __name__ == "__main__":
    unittest.main()
