# AMR Simulation Framework

An individual-based simulation framework for studying antimicrobial use,
antimicrobial resistance (AMR), infection outcomes, and potential policy
interventions across a broad bacterial and antibacterial ecosystem.

The current executable configuration represents:

- 42 bacteria
- 62 individual antibacterial drugs
- 39 internal drug classes
- 46 resistance mechanisms
- 6 world regions
- daily simulation from 1930 to 2025 for calibration, or to 2035 for policy
  branches

The framework is under active development and calibration. It is intended for
research and policy analysis, not for diagnosis, prescribing, or individual
clinical decision-making. A formal open-source licence has not yet been
selected, so this repository is not yet a released software version.

## Scientific Scope

The model follows simulated people through demographic change, bacterial
infection and carriage, symptoms and sepsis, diagnostic testing, antibacterial
treatment, resistance dynamics, and mortality.

Resistance can be introduced or altered through:

- resistance profiles sampled at infection or carriage acquisition
- de novo mechanism emergence under relevant drug pressure
- horizontal gene transfer for eligible mobile mechanisms
- mechanism reversion and fitness costs
- treatment selection and differential bacterial clearance
- regional and care-setting resistance history
- a bounded local persistence archive representing resistant strains circulating
  outside the finite simulated sample

Infection incidence is externally parameterised. The model does not dynamically
generate infection incidence from the prevalence of infected people. Its
dynamic population interaction concerns resistance profiles and treatment
pressure rather than explicit person-to-person transmission.

The framework deliberately does not model organism-drug MIC values or detailed
drug-specific PK/PD. Resistance severity is represented by bounded model
quantities such as `any_r`, while Figure 2 calibration prevalence is based on
whether `any_r > 0`.

See the [technical model description](model_description/MODEL_DESCRIPTION.md)
for the complete scientific and implementation specification.

## Inventory Authority

The executable inventories are defined in:

- `BACTERIA_LIST` in `src/simulation/population.rs`
- `DRUG_SHORT_NAMES` in `src/simulation/population.rs`
- `ResistanceMechanism::all()` in `src/simulation/population.rs`

Current parameter values are defined in `src/config.rs`. Files under `archive/`
are historical evidence only and are not model or analysis inputs. Repository
tests protect the dimensions and relationships between the executable
inventories, parameters, targets, and output schema.

## Requirements

- Rust stable with Cargo; the crate uses Rust edition 2021
- Python 3.10 or later for analysis
- Substantial RAM and runtime for research-scale simulations

The Python environment is currently specified by lower bounds in
`requirements.txt`; a frozen release environment is still to be added.

## Build and Test

```bash
cargo build --release
cargo test --all-targets
```

The default binary is `executable_amr`.

## Run the Simulation

```bash
cargo run --release
```

Run settings are currently selected near the top of `main()` in
`src/main.rs`. The checked-in configuration uses a population of 10,000,000,
`CalibrationMode::Full`, random seeding, and no individual or infection-journey
logging. This is a long research run, not a quick installation test.

For a smoke test only, temporarily use a much smaller `population_size`.
Outputs from a small population must not be interpreted as calibrated model
results.

### Run Modes

| Mode | Years and output |
|---|---|
| `FullMinimal` | Sparse 2022-2025 output containing drug share and bacteria-drug resistance fields |
| `Full` | Sparse 2022-2025 output containing all fields required by `calibration_summary.py` |
| `Partial` | Daily 1930-2025 output for historical time-series analysis |
| `Partial25Counterfactual` | Daily policy-0 rows for 1930-2025 plus no-resistance policy-2 rows for 2022-2025 |
| `Full25Counterfactual` | Sparse 2022-2025 output for policy 0 and the no-resistance policy 2 |
| `None` | Full 1930-2035 run with selected policy branches from 2027 |

`time_steps` is selected from the mode: 35,040 days for calibration and 2025
counterfactual modes, and 38,325 days for the full policy horizon. The two 2025
counterfactual modes checkpoint immediately before the first 2022 timestep,
retain the completed policy-0 trajectory, then restore that checkpoint and run
the resistance-suppressed policy 2 through the end of 2025.

Branch-enabled runs use disk-backed checkpoints by default. Checkpoint capture
streams borrowed population and mechanism-cache state directly to a temporary
file, then flushes, checksums, and atomically publishes it without cloning the
population in memory. After the baseline summaries have been retained, the
completed active population is released before a checkpoint is deserialized;
the restored state is moved into the policy branch. Multi-policy runs reread
the checkpoint for each branch to keep memory bounded. These modes therefore
still require space for the serialized checkpoint and retained summary rows,
but disk checkpointing does not require another simultaneous full population.
The explicitly selectable in-memory checkpoint path remains suitable only when
the additional population copy is known to fit. Before rerunning at 10 million,
use fixed-seed 300k, 1M, and 3M trials to review the emitted checkpoint size,
write/read duration, phase RSS/PSS, and runner cgroup peak. Provision checkpoint
disk capacity with explicit headroom; an interrupted process can leave a uniquely
named stale file in the configured checkpoint directory for deliberate cleanup.

To compare model-scope infection mortality between those two branches, set
`SIMULATION_CSV` near the top of
`amr_simulation_output_analysis/counterfactual_2025_death_rates.py`, then run:

```powershell
python amr_simulation_output_analysis/counterfactual_2025_death_rates.py
```

The script reports mean annual model-scope infection deaths, population-scaled
to the same world-population target used by `calibration_summary.py`, over
2022-2025 for policy 0 and policy 2. It also reports deaths per 100,000
person-years. It accepts output from either `Full25Counterfactual` or
`Partial25Counterfactual`. A CSV path supplied on the command line overrides
`SIMULATION_CSV`.

### Reproducible Seeds

The launcher generates and records a random `u64` seed unless fixed seeding is
enabled. `AMR_RNG_SEED` overrides the source setting and is the preferred way to
replay a run.

PowerShell:

```powershell
$env:AMR_RNG_SEED = "1234567890"
cargo run --release
```

Bash:

```bash
AMR_RNG_SEED=1234567890 cargo run --release
```

Fixed-seed runs use named ChaCha RNG streams and deterministic population
chunks. For the same source, configuration, and seed, summary output is
expected to be reproducible across repeated runs and different
`RAYON_NUM_THREADS` settings.

### CPU Threads

Rayon uses the available logical cores unless `RAYON_NUM_THREADS` is set.

PowerShell:

```powershell
$env:RAYON_NUM_THREADS = "4"
cargo run --release
```

Bash:

```bash
RAYON_NUM_THREADS=4 cargo run --release
```

## Outputs and Run Provenance

Simulation outputs are written under
`amr_simulation_output_analysis_outputs/`. A completed run normally produces:

- `simulation_summary_NNNNNN.csv`
- `run_metadata_<timestamp>_seed_<seed>.txt`
- `config_validation_<timestamp>.txt`

The summary CSV uses output schema version 1. Its fields depend on the selected
run mode and can number in the tens of thousands. The metadata records the
source hash, seed and seed source, run ID, population, time steps, mode,
policies, thread count, duration, output path, CSV SHA-256 hash, validation
status, and completion or failure state.

The source hash can be supplied by `AMR_SOURCE_HASH` or `source_hash.txt`.
Otherwise the launcher uses the current Git commit and marks a dirty worktree.
For formal analyses, retain the metadata file and exact source snapshot with
the CSV.

Parameter validation is strict by default. `AMR_CONFIG_VALIDATION=warn` permits
a diagnostic run to continue despite validation errors, but such a run should
not be used as a calibrated research result.

## Python Analysis

Create an environment and install the current analysis dependencies:

```bash
python -m venv .venv
python -m pip install -r requirements.txt
```

Activate the environment using the command appropriate for the operating
system. Select the input CSV through `DataConfig.simulation_file` in
`amr_simulation_output_analysis/config.py`, then run:

```bash
python -m amr_simulation_output_analysis.amr_analysis
```

The analysis writes calibration summaries and configured plots under
`output_graphs/`. Plot selection, policies, output format, caching, and memory
settings are controlled by `PlotConfig` in
`amr_simulation_output_analysis/config.py`.

Run the Python regression tests from the repository root with:

```bash
python -m unittest discover -s tests -p "test_*.py"
```

## Calibration Evidence

The current calibration combines sourced estimates, transformed comparisons,
evidence-informed benchmarks, and transparent expert-informed placeholders.
These categories must not be treated as interchangeable observations.

Key provenance documents are:

- [Resistance target provenance](data/RESISTANCE_TARGETS.md)
- [Comparison overlay provenance](data/empirical/OVERLAY_PROVENANCE.md)
- `data/resistance_targets_v1.manifest.json`

Best-guess placeholder overlays are disabled by default and are not calibration
score inputs.

## Policy Branches

`CalibrationMode::None` can run five independent branches from 2027:

| ID | Branch |
|---:|---|
| 0 | Baseline continuation |
| 1 | Antimicrobial stewardship example |
| 2 | Resistance-suppressed AMR counterfactual |
| 3 | Near-complete diagnostics bound |
| 4 | Equal global access example |

`CalibrationMode::Partial25Counterfactual` and
`CalibrationMode::Full25Counterfactual` instead compare policy 0 with policy 2
from the start of 2022 through the end of 2025. The former retains the complete
policy-0 history from 1930; the latter retains only the 2022-2025 calibration
window. In both modes, calibration summaries select policy 0 and do not mix the
counterfactual rows into baseline fit statistics.

These branches are research scenarios and should not be interpreted as
validated policy forecasts without scenario-specific calibration, uncertainty
analysis, and suitable comparison across well-fitting parameter sets.

## Repository Layout

```text
src/
  main.rs                         Run launcher and provenance
  config.rs                       Model parameters
  config_validation.rs            Parameter validation
  observability.rs                Run/source observability
  rules/mod.rs                    Daily individual-level rules
  simulation/
    population.rs                 State, inventories, and enums
    simulation.rs                 Simulation loop, caches, branches, CSV export
    journey_logger.rs             Optional sampled infection journeys

amr_simulation_output_analysis/   Python analysis package
data/                             Calibration targets and comparison data
model_description/                Technical model description
tests/                            Rust integration and Python regression tests
paper_tables/                     Generated manuscript tables and figures
archive/                          Historical, non-executable material
```

## Infection Journeys

Sampled infection-journey logging can be enabled in `src/main.rs` for
illustrative or diagnostic work. It records much denser individual traces and
can materially slow a run, so it is disabled for routine calibration.

## Contribution and Release Status

Contribution guidance, citation metadata, a stable run configuration
interface, a locked Python environment, and a formal release archive are still
being prepared. Until those are available, please treat the repository as an
active research workspace rather than a stable public API.

## Licence

No software licence has yet been selected. A standard open-source licence and
the appropriate UCL/contributor copyright notice must be added before formal
public release and reuse.
