#!/usr/bin/env python3
"""Fit the 108 initial human age/region weights to demographic targets.

This is a frozen-response, first-pass calibration, not a demographic simulation.
For each observed time step, a geometric overlap kernel maps uniform initial-age
bands to the five reported current-age groups. Let K be its time average, p the
baseline sampling probabilities (including the legacy missing-zero-band bug),
and o the mean observed daily joint region/age shares. The response coefficient
c[r,g] = o[r,g] / (K @ p[r])[g] makes B[r] = diag(c[r]) @ K reproduce o exactly.
These coefficients combine survival, migration, population normalization, and
sampling effects; they are NOT independently estimated survival probabilities.

For each region, solve B[r] @ q[r] = target_joint_share[r]. Among nonnegative
solutions, minimize squared relative departures from the region-scaled prior
plus a penalty on adjacent differences of those relative multipliers. Cohorts
not yet born anywhere in the window retain their region-scaled prior. Candidate
weights use the corrected sampler, including initial age [-4000, 0). Finally
normalize all 108 weights together to the baseline effective weight sum.

The normalized predicted living shares are B @ q / sum(B @ q). Global sampling
normalization also changes the implied living fraction for a fixed initial
number of agents; its ratio is recorded explicitly. Changing mortality or other
model dynamics invalidates the frozen-response assumption and requires a new
simulation. Broad age cells cannot identify cohort-specific survival, especially
within 80+, or validate the newly restored cohort's survival response.

Inputs and intermediate quantities are exported for audit. The exported daily
summary CSV and baseline_weights.csv can replace the original large parquet and
Rust config on subsequent runs; use --baseline-weights for the latter.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import platform
import re
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
from scipy.optimize import Bounds, LinearConstraint, minimize


REGIONS = ("asia", "africa", "europe", "north_america", "south_america", "oceania")
AGE_GROUPS = ("0_5", "6_14", "15_49", "50_79", "80plus")
INITIAL_EDGES = np.arange(-40000, 32001, 4000, dtype=int)
AGE_EDGES = np.array([0, 6 * 365, 15 * 365, 50 * 365, 80 * 365, np.inf])
ZERO_UPPER_BAND = 9


def parameter_name(region: str, low: int, high: int) -> str:
    def token(value: int) -> str:
        return f"neg{abs(value)}" if value < 0 else str(value)
    return f"demo_{region}_age_{token(low)}_{token(high)}"


def read_baseline(path: Path, sampler: str) -> tuple[np.ndarray, np.ndarray]:
    if path.suffix == ".csv":
        records = pd.read_csv(path, float_precision="round_trip")
        if records["parameter"].duplicated().any():
            raise ValueError("Duplicate parameters in baseline CSV")
        values = dict(zip(records["parameter"], records["configured_weight"]))
    else:
        text = path.read_text(encoding="utf-8")
        matches = re.findall(
            r'map\.insert\("(demo_[^"]+)"\.to_string\(\),\s*([0-9.eE+-]+)\)', text
        )
        if len(matches) != 108 or len(dict(matches)) != 108:
            raise ValueError("Expected exactly 108 unique demo_* entries in baseline config")
        values = {key: float(value) for key, value in matches}
    configured = np.array([
        [values[parameter_name(region, low, high)]
         for low, high in zip(INITIAL_EDGES[:-1], INITIAL_EDGES[1:])]
        for region in REGIONS
    ], dtype=float)
    if not np.isfinite(configured).all() or np.any(configured <= 0):
        raise ValueError("Baseline configured weights must be finite and positive")
    effective = configured.copy()
    if sampler == "legacy-negzero":
        effective[:, ZERO_UPPER_BAND] = 0
    return configured, effective


def read_summary(path: Path, start: int, end: int) -> pd.DataFrame:
    columns = ["time_step", "total_population"]
    for region in REGIONS:
        columns += [f"{region}_population"]
        columns += [f"{region}_prop_age_{age}" for age in AGE_GROUPS]
    frame = (pd.read_csv(path, usecols=columns, float_precision="round_trip") if path.suffix == ".csv"
             else pd.read_parquet(path, columns=columns))
    frame = frame.loc[frame.time_step.between(start, end), columns].sort_values("time_step")
    if frame.empty or frame.time_step.duplicated().any():
        raise ValueError("Calibration window must contain unique time steps")
    if not np.isfinite(frame.to_numpy(dtype=float)).all():
        raise ValueError("Summary contains missing or nonfinite demographic values")
    if not np.equal(frame.time_step, np.floor(frame.time_step)).all():
        raise ValueError("Summary time steps must be integer days")
    # Parquet preprocessing may downcast to float32. Promote before arithmetic
    # so the portable CSV reproduces the same calculations exactly.
    frame = frame.astype({column: (np.int64 if column == "time_step" else np.float64)
                          for column in columns})
    if np.any(frame.total_population <= 0):
        raise ValueError("Total living population must be positive")
    populations = frame[[f"{region}_population" for region in REGIONS]].to_numpy()
    if np.any(populations <= 0):
        raise ValueError("Every region must have a positive living population")
    if not np.allclose(populations.sum(axis=1), frame.total_population, rtol=1e-8):
        raise ValueError("Regional living populations do not sum to total population")
    for region in REGIONS:
        fractions = frame[[f"{region}_prop_age_{age}" for age in AGE_GROUPS]].to_numpy()
        if np.any(fractions < 0) or not np.allclose(fractions.sum(axis=1), 1, atol=1e-5):
            raise ValueError(f"Invalid reported age fractions for {region}")
    return frame.reset_index(drop=True)


def overlap_kernel(time_steps: np.ndarray, age_offset_days: int = 1) -> np.ndarray:
    """Exact integer-day band overlap, then equal weighting of reported days."""
    days = np.asarray(time_steps, dtype=float) + age_offset_days
    lower = np.maximum(INITIAL_EDGES[:-1][None, None, :],
                       AGE_EDGES[:-1][None, :, None] - days[:, None, None])
    upper = np.minimum(INITIAL_EDGES[1:][None, None, :],
                       AGE_EDGES[1:][None, :, None] - days[:, None, None])
    return (np.maximum(upper - lower, 0) / np.diff(INITIAL_EDGES)).mean(axis=0)


def observed_shares(frame: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    joint = []
    for region in REGIONS:
        population_share = (frame[f"{region}_population"] / frame.total_population).to_numpy()
        ages = frame[[f"{region}_prop_age_{age}" for age in AGE_GROUPS]].to_numpy(copy=True)
        # Printed fractions may have been rounded; normalize before averaging.
        ages /= ages.sum(axis=1, keepdims=True)
        joint.append((population_share[:, None] * ages).mean(axis=0))
    joint = np.array(joint)
    return joint, joint.sum(axis=1)


def read_age_targets(path: Path) -> tuple[np.ndarray, pd.DataFrame]:
    frame = pd.read_csv(path, float_precision="round_trip")
    expected = {(region, age) for region in REGIONS for age in AGE_GROUPS}
    if len(frame) != 30 or set(zip(frame.region, frame.age_group)) != expected:
        raise ValueError("Targets must contain exactly one row per region/age cell (30 rows)")
    indexed = frame.set_index(["region", "age_group"])
    values = np.array([[indexed.loc[(region, age), "age_share"]
                        for age in AGE_GROUPS] for region in REGIONS], dtype=float)
    if not np.isfinite(values).all() or np.any(values <= 0):
        raise ValueError("Target age shares must be positive finite fractions")
    if not np.allclose(values.sum(axis=1), 1, atol=1e-6, rtol=0):
        raise ValueError("Target age shares must sum to one within every region")
    return values / values.sum(axis=1, keepdims=True), frame


def fit_weights(configured: np.ndarray, effective: np.ndarray, kernel: np.ndarray,
                observed: np.ndarray, age_targets: np.ndarray, region_targets: np.ndarray,
                smoothness: float) -> dict:
    baseline = effective / effective.sum()
    baseline_geometry = baseline @ kernel.T
    if np.any(baseline_geometry <= 0) or np.any(observed <= 0):
        raise ValueError("Every age/region cell needs positive baseline and observed support")
    response = observed / baseline_geometry
    operator = response[:, :, None] * kernel[None, :, :]
    target = region_targets[:, None] * age_targets
    prior = configured / effective.sum() * (region_targets / observed.sum(axis=1))[:, None]
    fitted = np.zeros_like(prior)
    fits = []
    unobserved = np.sum(kernel, axis=0) == 0
    differences = np.diff(np.eye(18), axis=0)
    penalty = np.eye(18) + smoothness * differences.T @ differences
    for index, region in enumerate(REGIONS):
        scaled_operator = operator[index] * prior[index][None, :] / target[index][:, None]
        lower = np.full(18, 1e-10)
        upper = np.full(18, np.inf)
        lower[unobserved] = upper[unobserved] = 1

        def objective(multiplier):
            delta = multiplier - 1
            return float(delta @ penalty @ delta)

        def gradient(multiplier):
            return 2 * penalty @ (multiplier - 1)

        result = minimize(objective, np.ones(18), jac=gradient, method="SLSQP",
                          bounds=Bounds(lower, upper),
                          constraints=LinearConstraint(scaled_operator, 1, 1),
                          options={"ftol": 1e-12, "maxiter": 2000})
        residual = np.max(np.abs(scaled_operator @ result.x - 1))
        if not result.success or residual > 1e-7:
            raise RuntimeError(f"{region} fit failed: {result.message}; residual={residual}")
        fitted[index] = prior[index] * result.x
        fits.append({"region": region, "success": bool(result.success),
                     "iterations": int(result.nit), "objective": float(result.fun),
                     "max_relative_cell_residual": float(residual),
                     "minimum_prior_multiplier": float(result.x.min()),
                     "maximum_prior_multiplier": float(result.x.max())})
    probability = fitted / fitted.sum()
    raw_prediction = np.einsum("rgb,rb->rg", operator, probability)
    prediction = raw_prediction / raw_prediction.sum()
    return {"weights": probability * effective.sum(), "probability": probability,
            "prior_probability_before_normalization": prior, "response": response,
            "target": target, "prediction": prediction, "fits": fits,
            "global_living_population_ratio": float(raw_prediction.sum()),
            "common_normalization_multiplier": float(1 / fitted.sum()),
            "unobserved_bands": np.flatnonzero(unobserved).tolist()}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    baseline_input = parser.add_mutually_exclusive_group(required=True)
    baseline_input.add_argument("--baseline-config", type=Path)
    baseline_input.add_argument("--baseline-weights", type=Path)
    parser.add_argument("--baseline-sampler", choices=("legacy-negzero", "corrected"),
                        default="legacy-negzero", help="Sampler used by the recorded baseline run")
    parser.add_argument("--summary", type=Path, required=True, help="Parquet or exported daily CSV")
    parser.add_argument("--targets", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--start-day", type=int, default=33580)
    parser.add_argument("--end-day", type=int, default=35039)
    parser.add_argument("--age-offset-days", type=int, default=1,
                        help="Post-rule age increment at day zero; default 1")
    parser.add_argument("--africa-share", type=float, default=0.184351)
    parser.add_argument("--asia-share", type=float, default=0.589659)
    parser.add_argument("--smoothness", type=float, default=0.25)
    args = parser.parse_args()
    if min(args.africa_share, args.asia_share) <= 0 or args.africa_share + args.asia_share >= 1:
        parser.error("Africa/Asia shares must be positive and leave a positive remaining share")
    if args.smoothness < 0:
        parser.error("--smoothness must be nonnegative")
    baseline_path = args.baseline_config or args.baseline_weights
    configured, effective = read_baseline(baseline_path, args.baseline_sampler)
    frame = read_summary(args.summary, args.start_day, args.end_day)
    observed, regional = observed_shares(frame)
    age_targets, target_frame = read_age_targets(args.targets)
    kernel = overlap_kernel(frame.time_step.to_numpy(), args.age_offset_days)
    region_targets = regional.copy()
    remaining = np.array([region not in ("asia", "africa") for region in REGIONS])
    region_targets[remaining] *= (1 - args.africa_share - args.asia_share) / regional[remaining].sum()
    region_targets[REGIONS.index("asia")] = args.asia_share
    region_targets[REGIONS.index("africa")] = args.africa_share
    result = fit_weights(configured, effective, kernel, observed, age_targets,
                         region_targets, args.smoothness)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    weights = []
    rust_lines = ["// Frozen-response first-pass demographic fit; validate with a new simulation.",
                  "// Requires corrected sampler lookup for the age_neg4000_0 cohort."]
    for r, region in enumerate(REGIONS):
        rust_lines.append(f"// {region}")
        for band, (low, high) in enumerate(zip(INITIAL_EDGES[:-1], INITIAL_EDGES[1:])):
            name = parameter_name(region, low, high)
            value = result["weights"][r, band]
            weights.append({"parameter": name, "region": region,
                            "initial_age_min_days": int(low), "initial_age_max_days": int(high),
                            "configured_weight": float(configured[r, band]),
                            "baseline_effective_weight": float(effective[r, band]),
                            "fitted_weight": float(value),
                            "fitted_sampling_probability": float(result["probability"][r, band]),
                            "fitted_to_configured_ratio": float(value / configured[r, band]),
                            "unobserved_in_calibration_window": band in result["unobserved_bands"]})
            rust_lines.append(f'        map.insert("{name}".to_string(), {value:.12f});')
    weights_frame = pd.DataFrame(weights)
    weights_frame[["parameter", "region", "initial_age_min_days", "initial_age_max_days",
                   "configured_weight", "baseline_effective_weight"]].to_csv(
                       args.output_dir / "baseline_weights.csv", index=False)
    weights_frame.to_csv(args.output_dir / "fitted_weights.csv", index=False)
    (args.output_dir / "fitted_demographic_config.rs").write_text("\n".join(rust_lines) + "\n", encoding="utf-8")
    summary_csv = frame.to_csv(index=False, float_format="%.17g")
    (args.output_dir / "observed_demographic_summary.csv").write_text(summary_csv, encoding="utf-8")
    if args.targets.resolve() != (args.output_dir / "age_targets.csv").resolve():
        target_frame.to_csv(args.output_dir / "age_targets.csv", index=False)
    kernel_frame = pd.DataFrame(kernel, index=AGE_GROUPS,
                               columns=[f"{low}_{high}" for low, high in
                                        zip(INITIAL_EDGES[:-1], INITIAL_EDGES[1:])])
    kernel_frame.to_csv(args.output_dir / "age_overlap_kernel.csv", index_label="age_group")
    cells = []
    for r, region in enumerate(REGIONS):
        for g, age in enumerate(AGE_GROUPS):
            cells.append({"region": region, "age_group": age,
                          "baseline_region_share": float(regional[r]),
                          "target_region_share": float(region_targets[r]),
                          "baseline_age_share": float(observed[r, g] / regional[r]),
                          "target_age_share": float(age_targets[r, g]),
                          "predicted_age_share": float(result["prediction"][r, g] /
                                                       result["prediction"][r].sum()),
                          "baseline_joint_share": float(observed[r, g]),
                          "target_joint_share": float(result["target"][r, g]),
                          "predicted_joint_share": float(result["prediction"][r, g]),
                          "frozen_response_coefficient": float(result["response"][r, g])})
    pd.DataFrame(cells).to_csv(args.output_dir / "fit_cell_diagnostics.csv", index=False)
    diagnostics = {
        "status": "first_pass_frozen_response_fit_not_simulation_validated",
        "baseline_input": str(baseline_path),
        "baseline_input_sha256": hashlib.sha256(baseline_path.read_bytes()).hexdigest(),
        "summary_input": str(args.summary), "targets_input": str(args.targets),
        "targets_sha256": hashlib.sha256(args.targets.read_bytes()).hexdigest(),
        "selected_summary_sha256": hashlib.sha256(summary_csv.encode()).hexdigest(),
        "baseline_sampler": args.baseline_sampler, "candidate_sampler": "corrected",
        "calibration_day_min": int(frame.time_step.min()),
        "calibration_day_max": int(frame.time_step.max()), "calibration_rows": len(frame),
        "age_offset_days": args.age_offset_days, "model_days_per_year": 365,
        "smoothness": args.smoothness, "mean_baseline_population": float(frame.total_population.mean()),
        "software_versions": {"python": platform.python_version(), "numpy": np.__version__,
                              "pandas": pd.__version__, "scipy": scipy.__version__},
        "baseline_configured_weight_sum": float(configured.sum()),
        "baseline_effective_weight_sum": float(effective.sum()),
        "fitted_weight_sum": float(result["weights"].sum()),
        "global_living_population_ratio_frozen_response": result["global_living_population_ratio"],
        "common_normalization_multiplier": result["common_normalization_multiplier"],
        "unobserved_initial_age_bands": [[int(INITIAL_EDGES[b]), int(INITIAL_EDGES[b + 1])]
                                        for b in result["unobserved_bands"]],
        "max_absolute_joint_share_residual": float(np.max(np.abs(result["prediction"] - result["target"]))),
        "region_targets": dict(zip(REGIONS, region_targets.tolist())),
        "region_fits": result["fits"],
        "limitations": [
            "Response is inferred per broad region-age cell, not from individual survival histories.",
            "Oldest-age response and the restored [-4000, 0) cohort are weakly identified.",
            "Equal time weighting uses mean daily joint shares; kernel-response factorization is approximate.",
            "Unborn cohorts preserve region-scaled prior before the unavoidable common normalization.",
            "Age targets describe their source year, while the fit averages the stated calibration window.",
            "Changed mortality, migration, and other model dynamics require a new full-model validation run.",
        ],
    }
    (args.output_dir / "fit_diagnostics.json").write_text(json.dumps(diagnostics, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({key: diagnostics[key] for key in
                      ("status", "calibration_rows", "fitted_weight_sum",
                       "global_living_population_ratio_frozen_response",
                       "max_absolute_joint_share_residual", "region_targets")}, indent=2))


if __name__ == "__main__":
    main()
