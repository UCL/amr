"""Reproduce a conservative, parameter-only age-profile proposal.

Run from the repository root after extract_targets.py. Production model code is
not changed by this script. It writes an audit and a proposed values CSV.

Only the eight initial birth cohorts contributing to ages below 80 during
2022-2025 are fitted. Each region retains its current total sampling weight.
The not-yet-observed cohort and all pure-80+ overlap columns remain fixed:
the pre-fix baseline does not identify the restored cohort's survival.
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(ROOT))

import numpy as np
import pandas as pd
from scipy.optimize import Bounds, LinearConstraint, minimize

from amr_simulation_output_analysis.fit_demographic_weights import (
    AGE_GROUPS, INITIAL_EDGES, REGIONS, observed_shares, overlap_kernel,
    parameter_name, read_age_targets, read_baseline, read_summary,
)

HERE = Path(__file__).resolve().parent
BASELINE = ROOT / "archive/calibration_snapshots/regional_first_pass_2026-09-08/baseline_weights.csv"
SOURCE = ROOT / "amr_simulation_output_analysis_outputs/simulation_summary_103646.csv"
CURRENT = HERE / "before_weights.csv"
FREE = np.arange(1, 9)
SMOOTHNESS = 0.25
BOUNDS = (0.25, 2.0)


def main() -> None:
    # before_weights.csv is a frozen snapshot of the approved regional pass,
    # making reproduction independent of subsequent changes to src/config.rs.
    prior, _ = read_baseline(CURRENT, "corrected")
    _, legacy_effective = read_baseline(BASELINE, "legacy-negzero")
    frame = read_summary(SOURCE, 33580, 35039)
    observed, regional = observed_shares(frame)
    targets, _ = read_age_targets(HERE / "age_targets.csv")
    kernel = overlap_kernel(frame.time_step.to_numpy())
    baseline_probability = legacy_effective / legacy_effective.sum()
    response = observed / (baseline_probability @ kernel.T)
    fixed = np.setdiff1d(np.arange(18), FREE)
    assert np.all(kernel[:4, fixed] == 0)

    proposed = prior.copy()
    checks = []
    fit_rows = []
    differences = np.diff(np.eye(len(FREE)), axis=0)
    penalty = np.eye(len(FREE)) + SMOOTHNESS * differences.T @ differences
    for r, region in enumerate(REGIONS):
        target = targets[r, :4] / targets[r, :4].sum()
        current = prior[r, FREE]
        # Common regional/global normalization factors cancel in conditional
        # under-80 shares. Infer response from the OLD baseline, then apply it
        # to the CURRENT weights once; do not reapply regional multipliers.
        operator = response[r, :4, None] * kernel[:4, FREE] * current[None, :]
        equations = operator[:3] - target[:3, None] * operator.sum(axis=0)[None, :]
        equations /= np.max(np.abs(equations), axis=1)[:, None]
        constraints = np.vstack([equations, current / current.sum()])
        rhs = np.array([0.0, 0.0, 0.0, 1.0])
        assert np.linalg.matrix_rank(constraints) == 4

        result = minimize(
            lambda x: float((x - 1) @ penalty @ (x - 1)),
            np.ones(len(FREE)),
            jac=lambda x: 2 * penalty @ (x - 1),
            method="SLSQP",
            bounds=Bounds(np.full(len(FREE), BOUNDS[0]), np.full(len(FREE), BOUNDS[1])),
            constraints=LinearConstraint(constraints, rhs, rhs),
            options={"ftol": 1e-12, "maxiter": 2000},
        )
        if not result.success or np.max(np.abs(constraints @ result.x - rhs)) > 1e-8:
            raise RuntimeError(f"{region} fit failed: {result.message}")
        proposed[r, FREE] = current * result.x
        prediction = operator @ result.x
        prediction /= prediction.sum()
        assert np.allclose(prediction, target, atol=1e-10, rtol=0)
        assert abs(proposed[r].sum() - prior[r].sum()) < 1e-12
        assert np.array_equal(proposed[r, fixed], prior[r, fixed])
        checks.append({
            "region": region, "success": bool(result.success),
            "iterations": int(result.nit), "objective": float(result.fun),
            "minimum_multiplier": float(result.x.min()),
            "maximum_multiplier": float(result.x.max()),
            "current_region_weight": float(prior[r].sum()),
            "proposed_region_weight": float(proposed[r].sum()),
            "max_under80_share_error": float(np.max(np.abs(prediction - target))),
        })
        observed_conditional = observed[r, :4] / observed[r, :4].sum()
        for g, age in enumerate(AGE_GROUPS[:4]):
            fit_rows.append({
                "region": region, "age_group": age,
                "observed_within_under80_share": float(observed_conditional[g]),
                "target_within_under80_share": float(target[g]),
                "predicted_within_under80_share": float(prediction[g]),
                "observed_within_region_share": float(observed[r, g] / regional[r]),
                "un_within_region_share": float(targets[r, g]),
            })

    # Preserve frozen literals and ensure rounded fitted values retain each
    # region's total exactly. Put the <= a few 1e-12 rounding correction in the
    # largest fitted cohort, where it has negligible relative effect.
    with CURRENT.open(encoding="utf-8", newline="") as handle:
        before = {row["parameter"]: Decimal(row["configured_weight"]) for row in csv.DictReader(handle)}
    after = dict(before)
    for r, region in enumerate(REGIONS):
        names = [parameter_name(region, int(INITIAL_EDGES[b]), int(INITIAL_EDGES[b + 1])) for b in FREE]
        for b, name in zip(FREE, names):
            after[name] = Decimal(f"{proposed[r, b]:.12f}")
        residual = sum(before[name] for name in names) - sum(after[name] for name in names)
        largest = max(names, key=lambda name: after[name])
        after[largest] += residual
        assert sum(after[name] for name in names) == sum(before[name] for name in names)
    records = []
    for r, region in enumerate(REGIONS):
        for b, (low, high) in enumerate(zip(INITIAL_EDGES[:-1], INITIAL_EDGES[1:])):
            name = parameter_name(region, int(low), int(high))
            records.append({
                "parameter": name, "region": region,
                "before": str(before[name]), "after": str(after[name]),
                "multiplier": str(after[name] / before[name]),
                "fitted": b in FREE,
            })
    with (HERE / "proposed_weights.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)
    pd.DataFrame(fit_rows).to_csv(HERE / "fit_age_diagnostics.csv", index=False)
    frame.to_csv(HERE / "observed_demographic_summary.csv", index=False, float_format="%.17g")
    diagnostics = {
        "status": "provisional_under80_age_fit_not_simulation_validated",
        "baseline_run": "103646",
        "baseline_sampler": "legacy-negzero",
        "baseline_weights": str(BASELINE.relative_to(ROOT)),
        "candidate_sampler": "corrected",
        "current_weights_sha256": hashlib.sha256(CURRENT.read_bytes()).hexdigest(),
        "target_file_sha256": hashlib.sha256((HERE / "age_targets.csv").read_bytes()).hexdigest(),
        "window": {"start_day": 33580, "end_day": 35039, "rows": len(frame)},
        "fitted_cohorts_per_region": len(FREE),
        "fixed_cohorts_per_region": len(fixed),
        "multiplier_bounds": list(BOUNDS),
        "smoothness": SMOOTHNESS,
        "region_checks": checks,
        "limitations": [
            "Only the conditional age distribution below age 80 is fitted.",
            "Regional total sampling weights are preserved; living regional shares may still change with survival.",
            "Pure-80+ and not-yet-observed cohorts are fixed because their response is not identified by this baseline.",
            "The 70-85 overlapping cohort can still affect the 80+ population; its size is not independently calibrated here.",
            "The baseline predates the sampling-key fix and regional pass; its original weights and legacy mask are used explicitly.",
            "The response approximation holds mortality, migration and other population responses fixed.",
            "A new simulation is required to measure achieved age profiles and effects on other outputs.",
        ],
    }
    (HERE / "fit_diagnostics.json").write_text(json.dumps(diagnostics, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"status": diagnostics["status"], "checks": checks}, indent=2))


if __name__ == "__main__":
    main()
