import unittest
from pathlib import Path

import pandas as pd

from amr_simulation_output_analysis.empirical.provenance import (
    GENERATED_BEST_GUESS_PLACEHOLDER,
    OBSERVED_COMPARISON,
    SOURCE_INFORMED_BEST_GUESS_PLACEHOLDER,
    OverlayProvenanceError,
    annotate_overlay_provenance,
    filter_overlay_rows,
    load_overlay_provenance_registry,
)
from amr_simulation_output_analysis.empirical.acquire_empirical_data import (
    EmpiricalDataAcquirer,
)
from amr_simulation_output_analysis.empirical.enhanced_empirical_loader import (
    IntegratedEmpiricalLoader,
)


PROJECT_ROOT = Path(__file__).resolve().parents[1]


def _row(provenance_class: str, generated):
    return {
        "year": 2023,
        "drug": "ciprofloxacin",
        "bacteria": "escherichia_coli",
        "mean": 0.2,
        "source_quality": "WHO GLASS",
        "overlay_provenance_class": provenance_class,
        "generated": generated,
        "generation_method": (
            "random_source_pattern_generator" if generated else "not_generated"
        ),
        "source_id": "source-id",
        "source_url_or_doi": "https://example.test/source" if generated is False else "",
        "reference_year": "2023" if generated is False else "",
        "uncertainty": "reported interval" if generated is False else "not_empirical_uncertainty",
        "rationale": "Test provenance row.",
        "last_reviewed": "2026-07-22",
    }


class EmpiricalOverlayProvenanceTests(unittest.TestCase):
    def test_authoritative_source_label_alone_is_not_observed(self) -> None:
        frame = pd.DataFrame(
            {
                "year": [2023],
                "mean": [0.2],
                "source_quality": ["high_quality_surveillance"],
            }
        )

        default = filter_overlay_rows(frame)
        opted_in = filter_overlay_rows(
            frame,
            include_best_guess_placeholders=True,
        )

        self.assertTrue(default.empty)
        self.assertEqual(len(opted_in), 1)
        self.assertEqual(
            opted_in.iloc[0]["overlay_provenance_class"],
            SOURCE_INFORMED_BEST_GUESS_PLACEHOLDER,
        )
        self.assertFalse(bool(opted_in.iloc[0]["eligible_as_observed_comparison"]))

    def test_generated_placeholder_requires_explicit_opt_in(self) -> None:
        frame = pd.DataFrame([_row(GENERATED_BEST_GUESS_PLACEHOLDER, True)])

        self.assertTrue(filter_overlay_rows(frame).empty)
        self.assertEqual(
            len(
                filter_overlay_rows(
                    frame,
                    include_best_guess_placeholders=True,
                )
            ),
            1,
        )

    def test_complete_observed_row_is_visible_by_default(self) -> None:
        frame = pd.DataFrame([_row(OBSERVED_COMPARISON, False)])

        result = filter_overlay_rows(frame)

        self.assertEqual(len(result), 1)
        self.assertTrue(bool(result.iloc[0]["eligible_as_observed_comparison"]))

    def test_observed_row_without_source_url_is_rejected(self) -> None:
        row = _row(OBSERVED_COMPARISON, False)
        row["source_url_or_doi"] = ""

        with self.assertRaisesRegex(OverlayProvenanceError, "source_url_or_doi"):
            annotate_overlay_provenance(pd.DataFrame([row]))

    def test_registry_covers_every_tracked_legacy_source_label(self) -> None:
        registry = load_overlay_provenance_registry()
        tracked_files = sorted((PROJECT_ROOT / "data" / "empirical").glob("calibration_*_empirical.csv"))
        self.assertEqual(len(tracked_files), 7)

        for path in tracked_files:
            labels = set(
                pd.read_csv(
                    path,
                    usecols=["source_quality"],
                    keep_default_na=False,
                    na_values=[""],
                )["source_quality"]
                .fillna("")
                .astype(str)
            )
            relative = path.relative_to(PROJECT_ROOT).as_posix()
            rules = registry.loc[registry["relative_path"].eq(relative)]
            self.assertFalse(rules.empty, relative)
            selectors = set(rules["source_quality"])
            self.assertTrue("*" in selectors or labels <= selectors, relative)

            annotated = annotate_overlay_provenance(
                pd.DataFrame({"source_quality": sorted(labels)}),
                source_path=path,
            )
            self.assertFalse(annotated["eligible_as_observed_comparison"].any(), relative)

    def test_generated_source_files_cannot_be_upgraded_by_source_name(self) -> None:
        registry = load_overlay_provenance_registry()
        generated_paths = {
            "data/who/glass_amr_surveillance.csv",
            "data/ecdc/ears_net_surveillance.csv",
            "data/australia/nndss_surveillance.csv",
            "data/cddep/resistancemap_surveillance.csv",
        }
        rows = registry.loc[registry["relative_path"].isin(generated_paths)]

        self.assertEqual(set(rows["relative_path"]), generated_paths)
        self.assertTrue(
            rows["overlay_provenance_class"]
            .eq(GENERATED_BEST_GUESS_PLACEHOLDER)
            .all()
        )
        self.assertTrue(rows["generated"].str.lower().eq("true").all())

    def test_generation_helper_emits_explicit_placeholder_metadata(self) -> None:
        generated = EmpiricalDataAcquirer._mark_generated_placeholder(
            pd.DataFrame({"year": [2023], "resistance_percentage": [20.0]}),
            source_id="who_glass_pattern_template",
            generation_method="test_generator",
            rationale="Generated test pattern.",
        )

        annotated = annotate_overlay_provenance(generated)
        self.assertEqual(
            annotated.iloc[0]["overlay_provenance_class"],
            GENERATED_BEST_GUESS_PLACEHOLDER,
        )
        self.assertTrue(bool(annotated.iloc[0]["generated"]))
        self.assertFalse(bool(annotated.iloc[0]["eligible_as_observed_comparison"]))

    def test_integrated_loader_preserves_generated_classification(self) -> None:
        source = EmpiricalDataAcquirer._mark_generated_placeholder(
            pd.DataFrame(
                {
                    "year": [2023],
                    "pathogen": ["Escherichia coli"],
                    "antibiotic": ["Ciprofloxacin"],
                    "resistance_percentage": [20.0],
                }
            ),
            source_id="who_glass_pattern_template",
            generation_method="test_generator",
            rationale="Generated test pattern.",
        )
        loader = IntegratedEmpiricalLoader(include_best_guess_placeholders=True)

        processed = loader._process_resistance_surveillance_data(source, "who_glass")

        self.assertEqual(len(processed), 1)
        self.assertEqual(processed.iloc[0]["source_quality"], "who_glass")
        self.assertEqual(
            processed.iloc[0]["overlay_provenance_class"],
            GENERATED_BEST_GUESS_PLACEHOLDER,
        )
        self.assertNotIn("high_quality_surveillance", set(processed["source_quality"]))


if __name__ == "__main__":
    unittest.main()
