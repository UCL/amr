import unittest
from contextlib import ExitStack, redirect_stdout
from io import StringIO
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

import pandas as pd

import amr_simulation_output_analysis.make_paper_tables as paper_tables
from amr_simulation_output_analysis.make_paper_tables import (
    _CSV_HEADER_CACHE,
    _legacy_non_sf5_schema_validation,
    _paper_build_provenance_html,
    _paper_build_provenance_text,
    _paper_schema_contract_note,
    _preflight_simulation_csv_schemas,
    _read_csv_selected,
    _validate_reported_calibration_schemas,
)
from amr_simulation_output_analysis.summary_schema import (
    SUMMARY_SCHEMA_VERSION_COLUMN,
    SimulationSummarySchemaError,
)


def _write_summary(path: Path, versions: list[int]) -> None:
    pd.DataFrame(
        {
            SUMMARY_SCHEMA_VERSION_COLUMN: versions,
            "time_in_years": [92.0 + index / 365.0 for index in range(len(versions))],
            "policy_option": [0] * len(versions),
            "example_count": list(range(len(versions))),
        }
    ).to_csv(path, index=False)


_PAPER_OUTPUT_BUILDERS = (
    "make_t1",
    "make_supplementary_table_s2_resistance_benchmarks",
    "make_figure_1_calibration_headline_metrics",
    "make_figure_2_paper_parts",
    "make_figure_2a_hospital_resistance_fit",
    "make_figure_2b_community_resistance_fit",
    "make_figure_3_calibration_drug_class_share",
    "make_figure_4_calibration_infection_deaths",
    "make_figure_5_calibration_carriage_prevalence",
    "make_supplementary_figure_s1_potential_activity_retained",
    "make_supplementary_figure_s2_microbiome_resistance_reservoir",
    "make_supplementary_figure_s3_carrier_vs_non_carrier_incidence",
    "make_supplementary_figure_s5_diagnostic_testing_targeted_treatment_cascade",
    "make_supplementary_figure_s6_new_active_infection_denominators",
    "make_supplementary_figure_s8_infection_outcome_pathway",
    "make_figure_6_resistance_trend",
    "make_figure_6b_resistance_trend_by_bacterium",
    "make_figure_6c_serious_r_trend_by_bacterium",
    "make_figure_20_serious_r_by_hospital_community",
    "make_figure_7_infection_death_rate_by_region",
    "make_figure_8_antibiotic_use_by_context",
    "make_figure_11_sepsis_context_effective_therapy",
    "make_figure_15_mean_activity_by_bacteria",
    "make_figure_12_resistance_mechanisms_by_bacterium",
    "make_figure_13_active_infection_incidence",
)


def _run_mocked_paper_build(
    calibration_path: Path,
    csv_path: Path,
    output_dir: Path,
    *,
    legacy_without_sf5: bool,
) -> tuple[str, dict[str, object]]:
    schema_version = int(
        pd.read_csv(csv_path, usecols=[SUMMARY_SCHEMA_VERSION_COLUMN]).iloc[0, 0]
    )
    run = {
        "meta": {
            "source_file": str(calibration_path),
            "simulation_source_csv": str(csv_path),
            "simulation_summary_schema": f"{schema_version} "
            + ("(current)" if schema_version == 3 else "(legacy)"),
        }
    }
    mocks: dict[str, object] = {}
    args = [str(calibration_path)]
    if legacy_without_sf5:
        args.insert(0, paper_tables.LEGACY_WITHOUT_SF5_FLAG)

    stdout = StringIO()
    with ExitStack() as stack:
        stack.enter_context(patch.object(paper_tables, "OUT_DIR", output_dir))
        stack.enter_context(patch.object(paper_tables, "parse_files", return_value=[run]))
        stack.enter_context(
            patch.object(paper_tables, "aggregate", return_value={"n_runs": 1})
        )
        stack.enter_context(
            patch.object(
                paper_tables,
                "_discover_f1_simulation_csvs",
                return_value=[csv_path],
            )
        )
        stack.enter_context(
            patch.object(
                paper_tables,
                "_discover_simulation_csvs_with_scale",
                return_value=[(csv_path, 1.0)],
            )
        )
        for name in _PAPER_OUTPUT_BUILDERS:
            mocks[name] = stack.enter_context(patch.object(paper_tables, name))
        with redirect_stdout(stdout):
            paper_tables.main(args)

    return stdout.getvalue(), mocks


class LegacyPaperSchemaPreflightTests(unittest.TestCase):
    def setUp(self) -> None:
        _CSV_HEADER_CACHE.clear()

    def tearDown(self) -> None:
        _CSV_HEADER_CACHE.clear()

    def test_v1_and_v2_require_explicit_legacy_preflight(self) -> None:
        for schema_version in (1, 2):
            with self.subTest(schema_version=schema_version), TemporaryDirectory() as tmp:
                path = Path(tmp) / f"simulation_summary_v{schema_version}.csv"
                _write_summary(path, [schema_version, schema_version])

                with self.assertRaisesRegex(SimulationSummarySchemaError, "unsupported"):
                    _preflight_simulation_csv_schemas([path], allow_legacy=False)

                versions = _preflight_simulation_csv_schemas(
                    [path],
                    allow_legacy=True,
                )

                self.assertEqual(versions, {path: schema_version})

    def test_v3_uses_the_normal_strict_preflight(self) -> None:
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "simulation_summary_v3.csv"
            _write_summary(path, [3, 3])

            versions = _preflight_simulation_csv_schemas(
                [path],
                allow_legacy=False,
            )

        self.assertEqual(versions, {path: 3})

    def test_future_unversioned_and_mixed_inputs_remain_rejected(self) -> None:
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            future = root / "future.csv"
            unversioned = root / "unversioned.csv"
            mixed = root / "mixed.csv"
            _write_summary(future, [4])
            pd.DataFrame({"time_in_years": [92.0]}).to_csv(unversioned, index=False)
            _write_summary(mixed, [1, 3])

            for path in (future, unversioned, mixed):
                with self.subTest(path=path.name), self.assertRaises(
                    SimulationSummarySchemaError
                ):
                    _preflight_simulation_csv_schemas([path], allow_legacy=True)


class LegacyPaperEntrypointTests(unittest.TestCase):
    def setUp(self) -> None:
        _CSV_HEADER_CACHE.clear()

    def tearDown(self) -> None:
        _CSV_HEADER_CACHE.clear()

    def test_default_entrypoint_rejects_v1_before_output_cleanup(self) -> None:
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            calibration_path = root / "calibration_summary_078562.txt"
            calibration_path.write_text("Calibration Snapshot\n", encoding="utf-8")
            csv_path = root / "simulation_summary_078562.csv"
            _write_summary(csv_path, [1])
            output_dir = root / "paper_tables"

            with patch.object(paper_tables, "OUT_DIR", output_dir), patch.object(
                paper_tables,
                "parse_files",
                return_value=[{"meta": {"source_file": str(calibration_path)}}],
            ), patch.object(
                paper_tables,
                "aggregate",
                return_value={"n_runs": 1},
            ), patch.object(
                paper_tables,
                "_discover_f1_simulation_csvs",
                return_value=[csv_path],
            ), patch.object(
                paper_tables,
                "_discover_simulation_csvs_with_scale",
                return_value=[(csv_path, 1.0)],
            ), patch.object(paper_tables, "_prepare_output_dirs") as prepare:
                with redirect_stdout(StringIO()), self.assertRaises(SystemExit) as raised:
                    paper_tables.main([str(calibration_path)])

            self.assertEqual(raised.exception.code, 2)
            prepare.assert_not_called()

    def test_legacy_flag_omits_sf5_removes_stale_output_and_marks_index(self) -> None:
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            calibration_path = root / "calibration_summary_078562.txt"
            calibration_path.write_text("Calibration Snapshot\n", encoding="utf-8")
            csv_path = root / "simulation_summary_078562.csv"
            _write_summary(csv_path, [1])
            output_dir = root / "paper_tables"
            stale_sf5 = (
                output_dir
                / "Figures"
                / "Supplementary_Figure_S5__diagnostic_testing_targeted_treatment_cascade.html"
            )
            stale_sf5.parent.mkdir(parents=True)
            stale_sf5.write_text("stale", encoding="utf-8")

            stdout, mocks = _run_mocked_paper_build(
                calibration_path,
                csv_path,
                output_dir,
                legacy_without_sf5=True,
            )

            sf5_mock = mocks[
                "make_supplementary_figure_s5_diagnostic_testing_targeted_treatment_cascade"
            ]
            sf5_mock.assert_not_called()
            self.assertFalse(stale_sf5.exists())
            self.assertIn("Legacy paper compatibility is enabled", stdout)
            self.assertIn("will be omitted from the entire build", stdout)

            provenance = (output_dir / "build_provenance.txt").read_text(
                encoding="utf-8"
            )
            index = (output_dir / "index.html").read_text(encoding="utf-8")
            self.assertIn("schema 1 (legacy)", provenance)
            self.assertIn("Supplementary Figure S5: omitted", provenance)
            self.assertIn("Legacy compatibility build", index)
            self.assertNotIn(
                "href='Figures/Supplementary_Figure_S5__diagnostic_testing_targeted_treatment_cascade.html'",
                index,
            )

    def test_current_v3_entrypoint_invokes_sf5_without_legacy_warning(self) -> None:
        with TemporaryDirectory() as tmp:
            root = Path(tmp)
            calibration_path = root / "calibration_summary_123456.txt"
            calibration_path.write_text("Calibration Snapshot\n", encoding="utf-8")
            csv_path = root / "simulation_summary_123456.csv"
            _write_summary(csv_path, [3])
            output_dir = root / "paper_tables"

            stdout, mocks = _run_mocked_paper_build(
                calibration_path,
                csv_path,
                output_dir,
                legacy_without_sf5=False,
            )

            sf5_mock = mocks[
                "make_supplementary_figure_s5_diagnostic_testing_targeted_treatment_cascade"
            ]
            sf5_mock.assert_called_once()
            self.assertNotIn("Legacy paper compatibility is enabled", stdout)
            provenance = (output_dir / "build_provenance.txt").read_text(
                encoding="utf-8"
            )
            self.assertIn("Validation mode: current schema only", provenance)
            self.assertIn("Supplementary Figure S5: enabled", provenance)


class ScopedLegacyPaperValidationTests(unittest.TestCase):
    def setUp(self) -> None:
        _CSV_HEADER_CACHE.clear()

    def tearDown(self) -> None:
        _CSV_HEADER_CACHE.clear()

    def test_direct_reader_stays_strict_outside_legacy_scope(self) -> None:
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "simulation_summary_v1.csv"
            _write_summary(path, [1])

            with self.assertRaisesRegex(SimulationSummarySchemaError, "unsupported"):
                _read_csv_selected(path, {"example_count"})

            with _legacy_non_sf5_schema_validation(True):
                frame = _read_csv_selected(path, {"example_count"})

            self.assertEqual(frame[SUMMARY_SCHEMA_VERSION_COLUMN].tolist(), [1])
            with self.assertRaisesRegex(SimulationSummarySchemaError, "unsupported"):
                _read_csv_selected(path, {"example_count"})

    def test_legacy_scope_is_restored_after_an_exception(self) -> None:
        with TemporaryDirectory() as tmp:
            path = Path(tmp) / "simulation_summary_v2.csv"
            _write_summary(path, [2])

            with self.assertRaisesRegex(RuntimeError, "deliberate"):
                with _legacy_non_sf5_schema_validation(True):
                    self.assertIn("Supplementary Figure S5 is omitted", _paper_schema_contract_note())
                    raise RuntimeError("deliberate")

            self.assertIn("only the current", _paper_schema_contract_note())
            with self.assertRaisesRegex(SimulationSummarySchemaError, "unsupported"):
                _read_csv_selected(path, {"example_count"})


class LegacyPaperProvenanceTests(unittest.TestCase):
    def test_reported_legacy_schema_requires_explicit_mode(self) -> None:
        runs = [
            {
                "meta": {
                    "source_file": "calibration_summary_078562.txt",
                    "simulation_summary_schema": "1 (legacy)",
                }
            }
        ]

        with self.assertRaisesRegex(SimulationSummarySchemaError, "requires schema 3"):
            _validate_reported_calibration_schemas(runs, allow_legacy=False)

        _validate_reported_calibration_schemas(runs, allow_legacy=True)

    def test_legacy_provenance_records_omitted_sf5_and_each_csv_schema(self) -> None:
        input_paths = [Path("calibration_summary_078562.txt")]
        csv_v1 = Path("simulation_summary_078562.csv")
        csv_v3 = Path("simulation_summary_123456.csv")
        runs = [
            {
                "meta": {
                    "source_file": str(input_paths[0]),
                    "simulation_source_csv": str(csv_v1),
                    "simulation_summary_schema": "1 (legacy)",
                }
            }
        ]

        provenance = _paper_build_provenance_text(
            input_paths,
            runs,
            {csv_v1: 1, csv_v3: 3},
            legacy_without_sf5=True,
        )
        rendered = _paper_build_provenance_html(
            provenance,
            legacy_without_sf5=True,
        )

        self.assertIn("legacy compatibility (--legacy-without-sf5)", provenance)
        self.assertIn("Supplementary Figure S5: omitted", provenance)
        self.assertIn("simulation_summary_078562.csv", provenance)
        self.assertIn("schema 1 (legacy)", provenance)
        self.assertIn("schema 3 (current)", provenance)
        self.assertIn("not generated or represented by a placeholder", provenance)
        self.assertIn("Legacy compatibility build", rendered)
        self.assertIn("Supplementary Figure S5 was intentionally omitted", rendered)

    def test_current_provenance_does_not_claim_legacy_compatibility(self) -> None:
        csv_path = Path("simulation_summary_123456.csv")
        provenance = _paper_build_provenance_text(
            [Path("calibration_summary_123456.txt")],
            [],
            {csv_path: 3},
            legacy_without_sf5=False,
        )
        rendered = _paper_build_provenance_html(
            provenance,
            legacy_without_sf5=False,
        )

        self.assertIn("Validation mode: current schema only", provenance)
        self.assertIn("Supplementary Figure S5: enabled", provenance)
        self.assertNotIn("Legacy compatibility build", rendered)


if __name__ == "__main__":
    unittest.main()
