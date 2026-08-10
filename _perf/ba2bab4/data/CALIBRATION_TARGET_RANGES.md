# Calibration target plausible ranges

`calibration_target_ranges_v1.csv` is the display-only uncertainty companion to the
existing 2025 calibration targets. It covers:

- four headline targets;
- 28 antibiotic class-share targets;
- 42 bacterium-specific infection-incidence targets;
- 42 bacterium-specific carriage targets; and
- 42 bacterium-specific infection-death targets.

The canonical central values remain in `calibration_targets.json`,
`drug_class_share_history_targets.csv`, `infection_incidence_by_bacteria.csv`,
`microbiome_carriage_by_bacteria.csv`, and `deaths_by_bacteria.csv`. The registry
copies each central value so automated validation can detect drift.

## Interval meanings

- `published_uncertainty_range`: bounds reported by the named source. This label
  does not imply a 95% interval unless the source says so.
- `derived_plausible_range`: a published range has been transformed to the model's
  target definition, or a source-informed range has been rounded to contain the
  existing target.
- `expert_plausible_range`: an explicitly judgement-based range. It is not a
  confidence interval.
- `design_constraint`: a fixed target imposed by model scope or ecology rather
  than estimated with sampling uncertainty.

The broad default tiers are deliberately transparent:

- source-informed central estimates without recoverable intervals use a
  0.5-fold to 2-fold range;
- explicit placeholders use a 0.25-fold to 4-fold range;
- source-unresolved values use a one-third-fold to 3-fold range.

Explicit overrides replace these tiers where the canonical note contains a
quantitative range or where a narrower source-informed range is defensible.
Proportions are bounded to the interval from zero to one.

## Interpretation

These ranges communicate uncertainty in the target itself. They are separate
from:

- stochastic uncertainty in simulation means;
- calibration score tolerances;
- uncertainty across alternative parameterisations; and
- policy-effect uncertainty.

The ranges do not enter the calibration score. Drug-class ranges are marginal
ranges and are not a joint compositional confidence region; their lower and
upper bounds are not expected to sum to 100%.

The resistance target set is excluded from this first registry. Its versioned
file already has `uncertainty_lower` and `uncertainty_upper` fields, but the
cell-level source intervals have not been recovered. Those fields should remain
blank until ranges can be assigned without presenting generic generated bounds
as empirical uncertainty.

Rebuild the registry after an intentional target change with:

```powershell
python amr_simulation_output_analysis\build_calibration_target_ranges_v1.py
```
