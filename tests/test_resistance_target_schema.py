import csv
import json
import unittest
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from amr_simulation_output_analysis.build_resistance_targets_v1 import (
    PREVALENCE_COMPONENT,
    SEVERITY_COMPONENT,
    TARGET_COLUMNS,
    TARGET_SET_VERSION,
    build_resistance_targets_v1,
)
from amr_simulation_output_analysis.calibration_summary import (
    _load_resistance_target_set,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
TARGET_PATH = DATA / "resistance_targets_v1.csv"
SOURCE_PATH = DATA / "resistance_target_sources_v1.csv"


def _read_csv(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def _read_wide_values(path: Path):
    with path.open("r", encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle)
        drugs = [name for name in reader.fieldnames[1:] if name != "notes"]
        values = {}
        for row in reader:
            for drug in drugs:
                values[(row["Bacteria"], drug)] = "" if row[drug] == "." else row[drug]
    return drugs, values


def _row_key(row):
    return row["component"], row["bacteria"], row["drug"]


class ResistanceTargetSchemaTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls) -> None:
        cls.rows = _read_csv(TARGET_PATH)
        cls.sources = _read_csv(SOURCE_PATH)

    def test_generated_files_are_reproducible(self) -> None:
        with TemporaryDirectory() as temp_dir:
            target_output = Path(temp_dir) / TARGET_PATH.name
            source_output = Path(temp_dir) / SOURCE_PATH.name
            build_resistance_targets_v1(
                ROOT,
                target_output=target_output,
                source_output=source_output,
            )

            self.assertEqual(target_output.read_bytes(), TARGET_PATH.read_bytes())
            self.assertEqual(source_output.read_bytes(), SOURCE_PATH.read_bytes())

    def test_long_form_values_exactly_match_wide_matrices(self) -> None:
        prevalence_drugs, prevalence = _read_wide_values(
            DATA / "resistance_prevalence_values.csv"
        )
        severity_drugs, severity = _read_wide_values(
            DATA / "resistance_average_resistant_values.csv"
        )
        self.assertEqual(prevalence_drugs, severity_drugs)

        long_values = {_row_key(row): row["value"] for row in self.rows}
        self.assertEqual(len(long_values), len(self.rows))
        self.assertEqual(len(self.rows), 42 * 61 * 2)

        for (bacterium, drug), value in prevalence.items():
            self.assertEqual(
                long_values[(PREVALENCE_COMPONENT, bacterium, drug)], value
            )
        for (bacterium, drug), value in severity.items():
            self.assertEqual(
                long_values[(SEVERITY_COMPONENT, bacterium, drug)], value
            )

    def test_schema_and_provenance_contract_is_complete(self) -> None:
        with TARGET_PATH.open("r", encoding="utf-8-sig", newline="") as handle:
            self.assertEqual(csv.DictReader(handle).fieldnames, TARGET_COLUMNS)

        schema = json.loads(
            (DATA / "resistance_targets_v1.schema.json").read_text(encoding="utf-8")
        )
        self.assertEqual(schema["required"], TARGET_COLUMNS)
        self.assertEqual(set(schema["properties"]), set(TARGET_COLUMNS))

        source_ids = [row["source_id"] for row in self.sources]
        self.assertEqual(len(source_ids), 43)
        self.assertEqual(len(source_ids), len(set(source_ids)))
        known_sources = set(source_ids)

        for row in self.rows:
            self.assertEqual(row["target_set_version"], TARGET_SET_VERSION)
            if row["source_id"]:
                self.assertIn(row["source_id"], known_sources)
            self.assertEqual(row["evidence_weight"], "")
            if row["include_in_score"] == "true":
                self.assertNotEqual(row["value"], "")
                self.assertEqual(row["score_exclusion_reason"], "")
                self.assertEqual(row["score_row_weight"], "1.0")
            else:
                self.assertNotEqual(row["score_exclusion_reason"], "")
                self.assertEqual(row["score_row_weight"], "0.0")

        numeric_prevalence_types = {
            row["target_type"]
            for row in self.rows
            if row["component"] == PREVALENCE_COMPONENT and row["value"]
        }
        numeric_severity_types = {
            row["target_type"]
            for row in self.rows
            if row["component"] == SEVERITY_COMPONENT and row["value"]
        }
        self.assertEqual(
            numeric_prevalence_types,
            {"evidence_informed_calibration_benchmark"},
        )
        self.assertEqual(
            numeric_severity_types,
            {"expert_assigned_model_benchmark"},
        )

    def test_cell_statuses_make_legacy_missingness_explicit(self) -> None:
        counts = Counter((row["component"], row["cell_status"]) for row in self.rows)
        self.assertEqual(counts[(PREVALENCE_COMPONENT, "active_target")], 1230)
        self.assertEqual(
            counts[(PREVALENCE_COMPONENT, "inactive_model_unrepresentable")], 64
        )
        self.assertEqual(
            counts[(PREVALENCE_COMPONENT, "legacy_unclassified_missing")], 1268
        )
        self.assertEqual(counts[(SEVERITY_COMPONENT, "active_target")], 1229)
        self.assertEqual(
            counts[(SEVERITY_COMPONENT, "inactive_model_unrepresentable")], 63
        )
        self.assertEqual(
            counts[(SEVERITY_COMPONENT, "inactive_legacy_prevalence_gate")], 118
        )
        self.assertEqual(
            counts[(SEVERITY_COMPONENT, "legacy_unclassified_missing")], 1152
        )

    def test_include_flags_capture_static_score_eligibility(self) -> None:
        included = Counter(
            row["component"]
            for row in self.rows
            if row["include_in_score"] == "true"
        )
        self.assertEqual(included[PREVALENCE_COMPONENT], 1178)
        self.assertEqual(included[SEVERITY_COMPONENT], 1177)

        allowed_reasons = {
            "legacy_prevalence_target_missing",
            "severity_target_missing",
            "hard_exclusion_rifampicin",
            "hard_exclusion_mdr_tb",
            "hard_exclusion_listeria",
            "model_baseline_potency_below_0.15",
            "model_resistance_phenotype_not_representable",
        }
        for row in self.rows:
            reasons = {
                reason
                for reason in row["score_exclusion_reason"].split(";")
                if reason
            }
            self.assertTrue(reasons.issubset(allowed_reasons))

    def test_production_loader_preserves_component_values_and_inclusion(self) -> None:
        prevalence, severity = _load_resistance_target_set(TARGET_PATH)

        self.assertEqual(len(prevalence), 42 * 61)
        self.assertEqual(len(severity), 42 * 61)
        self.assertEqual(int(prevalence["include_in_score"].sum()), 1178)
        self.assertEqual(int(severity["include_in_score"].sum()), 1177)

        excluded = prevalence.loc[
            prevalence["Bacteria"].eq("Providencia stuartii")
            & prevalence["drug"].eq("ampicillin")
        ].iloc[0]
        self.assertEqual(excluded["target"], 0.65)
        self.assertFalse(excluded["include_in_score"])
        self.assertIn("negligible potency", excluded["reason"])

        missing_with_potency = prevalence.loc[
            prevalence["Bacteria"].eq("Acinetobacter baumannii")
            & prevalence["drug"].eq("fosfomycin")
        ].iloc[0]
        self.assertTrue(pd.isna(missing_with_potency["target"]))
        self.assertFalse(missing_with_potency["include_in_score"])
        self.assertIn("benchmark not assigned", missing_with_potency["reason"])
        self.assertNotIn("negligible potency", missing_with_potency["reason"])

        missing_with_low_potency = prevalence.loc[
            prevalence["Bacteria"].eq("Acinetobacter baumannii")
            & prevalence["drug"].eq("ampicillin")
        ].iloc[0]
        self.assertIn("benchmark not assigned", missing_with_low_potency["reason"])
        self.assertIn("negligible potency", missing_with_low_potency["reason"])

        included = prevalence.loc[
            prevalence["Bacteria"].eq("Escherichia coli")
            & prevalence["drug"].eq("ampicillin")
        ].iloc[0]
        self.assertEqual(included["target"], 0.58)
        self.assertTrue(included["include_in_score"])

        unrepresentable = prevalence.loc[
            prevalence["Bacteria"].eq("Acinetobacter baumannii")
            & prevalence["drug"].eq("cefiderocol")
        ].iloc[0]
        self.assertEqual(unrepresentable["target"], 0.05)
        self.assertFalse(unrepresentable["include_in_score"])
        self.assertIn("phenotype not represented", unrepresentable["reason"])


if __name__ == "__main__":
    unittest.main()
