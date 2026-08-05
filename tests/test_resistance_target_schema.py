import csv
import json
import shutil
import unittest
from collections import Counter
from pathlib import Path
from tempfile import TemporaryDirectory

import pandas as pd

from amr_simulation_output_analysis.build_resistance_targets_v1 import (
    PREVALENCE_COMPONENT,
    PROVENANCE_EVIDENCE_UNRESOLVED,
    PROVENANCE_EXPERT_PLACEHOLDER,
    PROVENANCE_STRUCTURAL_PRIOR,
    SEVERITY_COMPONENT,
    SOURCE_COLUMNS,
    TARGET_COLUMNS,
    TARGET_SET_VERSION,
    MANIFEST_FILENAME,
    build_resistance_targets_v1,
)
from amr_simulation_output_analysis.calibration_summary import (
    _load_resistance_target_set,
    _verify_resistance_target_manifest,
)


ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "data"
TARGET_PATH = DATA / "resistance_targets_v1.csv"
SOURCE_PATH = DATA / "resistance_target_sources_v1.csv"
MANIFEST_PATH = DATA / MANIFEST_FILENAME


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
            manifest_output = Path(temp_dir) / MANIFEST_PATH.name
            build_resistance_targets_v1(
                ROOT,
                target_output=target_output,
                source_output=source_output,
                manifest_output=manifest_output,
            )

            self.assertEqual(target_output.read_bytes(), TARGET_PATH.read_bytes())
            self.assertEqual(source_output.read_bytes(), SOURCE_PATH.read_bytes())
            self.assertEqual(manifest_output.read_bytes(), MANIFEST_PATH.read_bytes())

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
        with SOURCE_PATH.open("r", encoding="utf-8-sig", newline="") as handle:
            self.assertEqual(csv.DictReader(handle).fieldnames, SOURCE_COLUMNS)

        schema = json.loads(
            (DATA / "resistance_targets_v1.schema.json").read_text(encoding="utf-8")
        )
        self.assertEqual(schema["required"], TARGET_COLUMNS)
        self.assertEqual(set(schema["properties"]), set(TARGET_COLUMNS))
        self.assertTrue(
            {row["cell_status"] for row in self.rows}.issubset(
                set(schema["properties"]["cell_status"]["enum"])
            )
        )
        self.assertTrue(
            {row["provenance_class"] for row in self.rows}.issubset(
                set(schema["properties"]["provenance_class"]["enum"])
            )
        )

        source_ids = [row["source_id"] for row in self.sources]
        self.assertEqual(len(source_ids), 45)
        self.assertEqual(len(source_ids), len(set(source_ids)))
        known_sources = set(source_ids)
        source_provenance = {
            row["source_id"]: row["provenance_class"] for row in self.sources
        }

        for row in self.rows:
            self.assertEqual(row["target_set_version"], TARGET_SET_VERSION)
            if row["source_id"]:
                self.assertIn(row["source_id"], known_sources)
            if row["value"]:
                self.assertEqual(
                    row["provenance_class"],
                    source_provenance[row["source_id"]],
                )
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
        self.assertEqual(
            {
                row["provenance_class"]
                for row in self.rows
                if row["component"] == PREVALENCE_COMPONENT and row["value"]
            },
            {PROVENANCE_EVIDENCE_UNRESOLVED},
        )
        self.assertEqual(
            {
                row["provenance_class"]
                for row in self.rows
                if row["component"] == SEVERITY_COMPONENT and row["value"]
            },
            {PROVENANCE_EXPERT_PLACEHOLDER, PROVENANCE_STRUCTURAL_PRIOR},
        )
        self.assertFalse(
            any(
                row["provenance_class"]
                == "empirical_estimate_with_cell_level_source"
                for row in self.rows
            )
        )

    def test_hash_manifest_covers_and_verifies_every_dependency(self) -> None:
        _verify_resistance_target_manifest(TARGET_PATH)
        manifest = json.loads(MANIFEST_PATH.read_text(encoding="utf-8"))
        self.assertEqual(manifest["target_set_version"], TARGET_SET_VERSION)
        self.assertEqual(manifest["hash_algorithm"], "sha256")
        self.assertEqual(
            set(manifest["artifacts"]),
            {
                "resistance_targets_v1.csv",
                "resistance_target_sources_v1.csv",
                "resistance_targets_v1.schema.json",
                "resistance_prevalence_values.csv",
                "resistance_average_resistant_values.csv",
                "model_potency_matrix.csv",
                "model_resistance_reachability_matrix.csv",
            },
        )

        with TemporaryDirectory() as temp_dir:
            temp_data = Path(temp_dir)
            shutil.copy2(MANIFEST_PATH, temp_data / MANIFEST_PATH.name)
            for name in manifest["artifacts"]:
                shutil.copy2(DATA / name, temp_data / name)
            with (temp_data / TARGET_PATH.name).open("ab") as handle:
                handle.write(b"\n")
            with self.assertRaisesRegex(ValueError, "size mismatch"):
                _verify_resistance_target_manifest(temp_data / TARGET_PATH.name)

    def test_loader_rejects_malformed_provenance_rows(self) -> None:
        cases = {
            "unknown provenance class": lambda rows: rows[0].__setitem__(
                "provenance_class", "unsupported_claim"
            ),
            "numeric row without source": lambda rows: next(
                row for row in rows if row["value"]
            ).__setitem__("source_id", ""),
            "invented evidence weight": lambda rows: next(
                row for row in rows if row["value"]
            ).__setitem__("evidence_weight", "1.0"),
            "source provenance mismatch": lambda rows: next(
                row for row in rows if row["value"]
            ).__setitem__("provenance_class", PROVENANCE_STRUCTURAL_PRIOR),
        }
        for label, mutate in cases.items():
            with self.subTest(label=label), TemporaryDirectory() as temp_dir:
                temp_data = Path(temp_dir)
                rows = [dict(row) for row in self.rows]
                mutate(rows)
                malformed_path = temp_data / TARGET_PATH.name
                with malformed_path.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.DictWriter(handle, fieldnames=TARGET_COLUMNS)
                    writer.writeheader()
                    writer.writerows(rows)
                shutil.copy2(SOURCE_PATH, temp_data / SOURCE_PATH.name)
                with self.assertRaises(ValueError):
                    _load_resistance_target_set(
                        malformed_path,
                        verify_manifest=False,
                    )

    def test_cell_statuses_make_legacy_missingness_explicit(self) -> None:
        counts = Counter((row["component"], row["cell_status"]) for row in self.rows)
        self.assertEqual(counts[(PREVALENCE_COMPONENT, "active_target")], 1238)
        self.assertEqual(
            counts[
                (
                    PREVALENCE_COMPONENT,
                    "active_target_model_unrepresentable",
                )
            ],
            56,
        )
        self.assertEqual(
            counts[(PREVALENCE_COMPONENT, "legacy_unclassified_missing")], 1268
        )
        self.assertEqual(counts[(SEVERITY_COMPONENT, "active_target")], 1159)
        self.assertEqual(
            counts[
                (
                    SEVERITY_COMPONENT,
                    "inactive_above_model_representable_maximum",
                )
            ],
            78,
        )
        self.assertEqual(
            counts[(SEVERITY_COMPONENT, "inactive_model_unrepresentable")], 55
        )
        self.assertEqual(
            counts[(SEVERITY_COMPONENT, "inactive_unpaired_legacy_benchmark")],
            118,
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
        self.assertEqual(included[PREVALENCE_COMPONENT], 1228)
        self.assertEqual(included[SEVERITY_COMPONENT], 1108)

        allowed_reasons = {
            "legacy_prevalence_target_missing",
            "severity_target_missing",
            "hard_exclusion_rifampicin",
            "hard_exclusion_mdr_tb",
            "hard_exclusion_listeria",
            "model_baseline_potency_below_0.15",
            "model_resistance_phenotype_not_representable",
            "severity_benchmark_above_model_representable_maximum",
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
        self.assertEqual(int(prevalence["include_in_score"].sum()), 1228)
        self.assertEqual(int(severity["include_in_score"].sum()), 1108)

        structural_gap = prevalence.loc[
            prevalence["Bacteria"].eq("Enterococcus faecium")
            & prevalence["drug"].eq("quinu_dalfo")
        ].iloc[0]
        self.assertEqual(structural_gap["target"], 0.5)
        self.assertTrue(structural_gap["include_in_score"])
        self.assertIn(
            "resistance phenotype not represented by model mechanisms",
            structural_gap["reason"],
        )

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
        self.assertEqual(
            included["provenance_class"], PROVENANCE_EVIDENCE_UNRESOLVED
        )
        self.assertEqual(
            included["source_id"], "legacy_prevalence_note__escherichia_coli"
        )
        self.assertEqual(included["rationale"], "legacy_bacterium_level_note")

        represented = prevalence.loc[
            prevalence["Bacteria"].eq("Acinetobacter baumannii")
            & prevalence["drug"].eq("cefiderocol")
        ].iloc[0]
        self.assertEqual(represented["target"], 0.05)
        self.assertTrue(represented["include_in_score"])

        unrepresentable = prevalence.loc[
            prevalence["Bacteria"].eq("Neisseria gonorrhoeae")
            & prevalence["drug"].eq("cefiderocol")
        ].iloc[0]
        self.assertEqual(unrepresentable["target"], 0.01)
        self.assertTrue(unrepresentable["include_in_score"])
        self.assertIn("phenotype not represented", unrepresentable["reason"])

        unpaired = severity.loc[
            severity["Bacteria"].eq("Helicobacter pylori")
            & severity["drug"].eq("cefiderocol")
        ].iloc[0]
        self.assertFalse(unpaired["include_in_score"])
        self.assertIn("paired prevalence benchmark not assigned", unpaired["reason"])

    def test_expert_severity_placeholder_provenance_is_explicit(self) -> None:
        severity = [
            row for row in self.rows if row["component"] == SEVERITY_COMPONENT
        ]
        reserve = [
            row
            for row in severity
            if row["source_id"] == "expert_reserve_drug_any_r_placeholders_v1"
        ]
        self.assertEqual(len(reserve), 48)
        self.assertEqual(
            Counter((row["drug"], row["value"]) for row in reserve),
            Counter(
                {
                    ("cefiderocol", "0.60"): 24,
                    ("ceftolozane_tazobactam", "0.70"): 24,
                }
            ),
        )
        self.assertEqual(
            {row["rationale"] for row in reserve},
            {"expert_best_guess_reserve_drug_any_r_placeholder"},
        )
        self.assertEqual(
            {row["provenance_class"] for row in reserve},
            {PROVENANCE_EXPERT_PLACEHOLDER},
        )

        rare_positive = [
            row
            for row in severity
            if row["source_id"] == "expert_rare_positive_any_r_prior_v1"
        ]
        self.assertEqual(
            {(row["bacteria"], row["drug"]) for row in rare_positive},
            {
                ("Staphylococcus aureus", "vancomycin"),
                ("Staphylococcus epidermidis", "vancomycin"),
                ("Streptococcus pneumoniae", "vancomycin"),
                ("Streptococcus pneumoniae", "linezolid"),
                ("Streptococcus pneumoniae", "daptomycin"),
            },
        )
        self.assertEqual({row["value"] for row in rare_positive}, {"0.60"})
        self.assertEqual(
            {row["rationale"] for row in rare_positive},
            {"expert_rare_positive_any_r_structural_prior"},
        )
        self.assertEqual(
            {row["provenance_class"] for row in rare_positive},
            {PROVENANCE_STRUCTURAL_PRIOR},
        )


if __name__ == "__main__":
    unittest.main()
