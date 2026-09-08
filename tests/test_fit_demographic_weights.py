"""Check demographic age geometry and normalization against independent invariants."""

import tempfile
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from amr_simulation_output_analysis.fit_demographic_weights import (
    AGE_EDGES, AGE_GROUPS, INITIAL_EDGES, REGIONS, fit_weights, overlap_kernel,
    parameter_name, read_baseline, read_age_targets,
)


class DemographicFitTests(unittest.TestCase):
    def test_overlap_matches_enumerated_integer_ages(self):
        days = np.array([0, 2189, 33580, 35039])
        expected = np.zeros((5, 18))
        for day in days:
            for band, (low, high) in enumerate(zip(INITIAL_EDGES[:-1], INITIAL_EDGES[1:])):
                ages = np.arange(low, high) + day + 1
                for group, (young, old) in enumerate(zip(AGE_EDGES[:-1], AGE_EDGES[1:])):
                    expected[group, band] += np.mean((ages >= young) & (ages < old)) / len(days)
        np.testing.assert_allclose(overlap_kernel(days), expected, atol=1e-15)

    def test_legacy_sampler_masks_only_the_zero_upper_band(self):
        records = [{"parameter": parameter_name(region, low, high), "configured_weight": 1.0}
                   for region in REGIONS for low, high in zip(INITIAL_EDGES[:-1], INITIAL_EDGES[1:])]
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "weights.csv"
            pd.DataFrame(records).to_csv(path, index=False)
            configured, legacy = read_baseline(path, "legacy-negzero")
            _, corrected = read_baseline(path, "corrected")
        np.testing.assert_equal(configured, corrected)
        np.testing.assert_equal(np.flatnonzero(np.any(configured != legacy, axis=0)), [9])
        np.testing.assert_equal(legacy[:, 9], 0)
        self.assertEqual(parameter_name("asia", -4000, 0), "demo_asia_age_neg4000_0")

    def test_feasible_synthetic_targets_global_normalization_and_unborn_prior(self):
        # Build feasible changed targets independently from a known reweighting.
        configured = np.tile(np.linspace(0.02, 0.005, 18), (6, 1))
        effective = configured.copy()
        effective[:, 9] = 0
        kernel = overlap_kernel(np.arange(33580, 35040))
        baseline = effective / effective.sum()
        age_survival = np.array([1, 0.95, 0.85, 0.6, 0.15])
        population = (baseline @ kernel.T) * age_survival
        observed = population / population.sum()
        known_weights = configured.copy()
        known_weights[:, 1:4] *= np.linspace(0.8, 1.2, 6)[:, None]
        known_population = (known_weights @ kernel.T) * age_survival
        joint_targets = known_population / known_population.sum()
        region_targets = joint_targets.sum(axis=1)
        age_targets = joint_targets / region_targets[:, None]
        result = fit_weights(configured, effective, kernel, observed, age_targets, region_targets, 0.25)
        np.testing.assert_allclose(result["prediction"], joint_targets, rtol=1e-8)
        self.assertTrue(np.all(result["weights"] > 0))
        self.assertAlmostEqual(result["weights"].sum(), effective.sum())
        self.assertEqual(result["unobserved_bands"], [0])
        expected_unborn = (result["prior_probability_before_normalization"][:, 0]
                           * result["common_normalization_multiplier"])
        np.testing.assert_allclose(result["probability"][:, 0], expected_unborn)
        # A sampler uses global normalization, so globally scaling the config
        # cannot change either predicted shares or the living-population ratio.
        scaled = fit_weights(configured * 7, effective * 7, kernel, observed,
                             age_targets, region_targets, 0.25)
        np.testing.assert_allclose(scaled["probability"], result["probability"])
        self.assertAlmostEqual(scaled["global_living_population_ratio"],
                               result["global_living_population_ratio"])

    def test_duplicate_or_missing_age_target_is_rejected(self):
        rows = [{"region": region, "age_group": age, "age_share": 0.2}
                for region in REGIONS for age in AGE_GROUPS]
        rows[-1] = rows[0].copy()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "targets.csv"
            pd.DataFrame(rows).to_csv(path, index=False)
            with self.assertRaisesRegex(ValueError, "exactly one row"):
                read_age_targets(path)


if __name__ == "__main__":
    unittest.main()
