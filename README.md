# AMR Simulation

Agent-based simulation of antimicrobial resistance dynamics across 42 bacterial species, 58 antibiotics, 35 resistance mechanisms, and 6 world regions over the period 1930–2035.

## Quick Start

### Prerequisites

- **Rust** (stable, edition 2021)
- **Cargo** (comes with Rust)
- Optional: Python 3.10+ for output analysis

### Build and Run

```bash
cargo build --release
cargo run --release
```

The executable is named `executable_amr`. All configuration is hardcoded in `src/main.rs` — there are no command-line arguments. To adjust the simulation, edit the constants at the top of `main()`:

| Variable | Default | Purpose |
|----------|---------|---------|
| `population_size` | 100,000 | Number of simulated individuals |
| `time_steps` | 38,325 | Days from 1930 to ~2035 |
| `log_individuals` | `false` | Per-individual logging (very verbose) |
| `log_infection_journeys` | `false` | Infection lifecycle logging |
| `infection_journey_sample_rate` | 1.0 | Fraction of journeys to log |
| `use_fixed_seed` | `false` | Deterministic RNG for reproducibility |
| `fixed_seed_value` | 1,234,567,890 | Seed when fixed seeding is on |
| `infection_journey_bacteria_filter` | `None` | Filter journeys to one species |
| `use_disk_branch_checkpointing` | `false` | Serialise branch checkpoints to disk |

### Limiting CPU Cores

```powershell
$env:RAYON_NUM_THREADS = "4"
cargo run --release
```

### Remote Submission

```powershell
.\submit.ps1
```

Copies the project to a remote Linux build runner via SSH, compiles, and runs there.

## Project Structure

```
src/
├── main.rs                      Entry point, run configuration, validation
├── config.rs                    All parameters (~11,700 lines), HashMap-based
├── rules/
│   └── mod.rs                   Daily update rules (~5,800 lines)
└── simulation/
    ├── mod.rs                   Module declarations
    ├── simulation.rs            Main simulation loop, CSV export (~5,200 lines)
    ├── population.rs            Individual struct, enums, constants (~1,700 lines)
    └── journey_logger.rs        Infection journey CSV logger

amr_simulation_output_analysis/  Python analysis package
├── amr_analysis.py              Main analysis entry point
├── calibration_summary.py       Calibration metric computation
├── column_selector.py           Dynamic column selection from wide CSVs
├── config.py                    Analysis configuration
├── data_loader.py               Pandas-based CSV loader
├── polars_loader.py             Polars-based fast loader
├── multi_run_activity_r_plot.py Multi-run resistance comparison
├── utils.py                     Shared utilities
├── empirical/                   Empirical data loading & normalisation
│   ├── acquire_empirical_data.py
│   ├── data_loader.py
│   ├── enhanced_empirical_loader.py
│   └── normalizers.py
└── plotting/                    Visualisation modules
    ├── base_plots.py
    ├── detail_plots.py
    └── grouped_plots.py

amr_simulation_output_analysis_outputs/
    simulation_summary_NNNNNN.csv   Output CSVs (one per run)

archive/                         Historical code and process documentation
data/                            Empirical reference data
infection_journeys/              Infection journey log output
output_graphs/                   Generated plots
```

## Architecture

### Simulation Loop

`Simulation::run()` iterates over `time_steps` days. Each day:

1. The population is partitioned across threads with **Rayon** (`par_chunks_mut`).
2. Each thread applies the daily rule sequence to its individuals and accumulates per-thread `LocalTotals`.
3. Thread-local totals are merged into a global `TimestepSummary`.
4. Population-level caches (`MajorityRCache`, `MechanismProfileCache`) are rebuilt from the merged totals for the next day.

### Daily Rule Sequence

Applied in order for each individual (see `rules/mod.rs`):

| # | Rule | Description |
|---|------|-------------|
| 1 | Age update | Increment age, assign age category |
| 2 | Background mortality | Non-infection death (age/region/sex-dependent) |
| 3 | Immunodeficiency transitions | Onset/recovery of transient or chronic immunosuppression |
| 4 | Hospitalisation | Admission/discharge based on clinical state |
| 5 | Travel | Region reassignment |
| 6 | Infection acquisition | Community/hospital/carrier-derived new infections |
| 7 | Symptom onset | Infection → symptomatic transition |
| 8 | Sepsis onset | Symptomatic → sepsis progression |
| 9 | Sepsis resolution | Recovery or death from sepsis |
| 10 | Bacterial testing | Culture ordering and result |
| 11 | Resistance testing | AST ordering and result |
| 12 | Drug initiation | Empiric or targeted therapy selection |
| 13 | Drug level update | PK decay, toxicity reservoir accumulation |
| 14 | Drug toxicity check | Sub-lethal discontinuation, lethal toxicity death |
| 15 | Drug efficacy | Active drug reduces bacterial load |
| 16 | Drug treatment stop | Course completion, clinical response |
| 17 | Natural clearance | Immune-mediated infection resolution |
| 18 | Resistance emergence | De novo mutation under drug pressure |
| 19 | Resistance reversion | Fitness-cost-driven loss of resistance |
| 20 | Microbiome dynamics | Resistance promotion/decay in carriage flora |
| 21 | HGT | Horizontal gene transfer (mechanism-driven, compartment-aware) |

### Key Data Structures

#### `Individual` (population.rs)

~60 fields organised into:

- **Demographics**: `age_days`, `sex`, `region`, `age_category`
- **Location**: `hospital_status`, `hospital_days`
- **Infection state** (per-bacteria `Vec`): `cur_bacteria_level`, `cur_infection_duration`, `cur_syndrome`, `symptomatic`, `sepsis`
- **Microbiome** (per-bacteria `Vec`): `is_carrier`, `carriage_days`
- **Treatment** (per-drug `Vec`): `cur_drug_level`, `drug_treatment_days`, `toxicity_reservoir`, `toxicity_stopped_drug_day`
- **Resistance** (bacteria × mechanism `Vec<Vec<bool>>`): `cur_resistance`
- **Computed resistance** (per-bacteria): `Resistance` struct with `microbiome_r`, `test_r`, `activity_r`, `any_r`, `majority_r`
- **Testing**: `test_bacteria_result`, `test_r_result`, `test_day`
- **Mortality/risk**: `immunodeficiency_type`, `not_under_care`
- **Tracking**: `drugs_taken`, `lifetime_infection_count`, `death_cause`

#### `Resistance` Struct

| Field | Type | Meaning |
|-------|------|---------|
| `microbiome_r` | `f64` | Fraction of resistance in gut/carriage flora |
| `test_r` | `f64` | Last AST result (0 or 1) |
| `activity_r` | `f64` | Computed resistance affecting drug efficacy |
| `any_r` | `f64` | Max of microbiome/test/activity |
| `majority_r` | `f64` | Population-level prevalence (from cache) |

#### Parameter System (config.rs)

All ~5,000+ parameters are stored in a global `lazy_static` `HashMap<String, f64>`. Parameters are looked up by string key at runtime:

```rust
let val = get_param("sepsis_death_base_log_odds");
let val = get_bacteria_param("escherichia_coli", "acquisition_log_odds");
let val = get_region_param("africa", "hospitalization_log_odds");
```

Key naming conventions:
- `bacteria_{name}_{param}` — per-bacteria parameters
- `drug_{name}_{param}` — per-drug parameters
- `{region}_{param}` — per-region parameters
- `{bacteria}_{drug}_potency_when_no_r` — potency matrix entries
- `hgt_prob_{source}_to_{target}` — HGT probability matrix

See [MODEL_DESCRIPTION.md](MODEL_DESCRIPTION.md) for the complete parameter reference.

#### Caches

| Cache | Rebuilt | Purpose |
|-------|---------|---------|
| `MajorityRCache` | Every timestep | Rolling window of population-level resistance prevalence per bacteria×drug |
| `MechanismProfileCache` | Every timestep | Reservoir sample (≤200) of resistance mechanism profiles from infected individuals per bacteria |

### Policy System

`PolicyAdjustments` defines three built-in policy branches that diverge at `POLICY_BRANCH_YEAR` (2027):

| Branch | Description |
|--------|-------------|
| `baseline` | No intervention changes |
| `stewardship` | Antibiotic stewardship adjustments |
| `counterfactual` | Alternative scenario |

At the branch point, the simulation serialises (or clones) the full population state and runs each branch independently to the end.

### Key Constants

| Constant | Value | Location |
|----------|-------|----------|
| `INFECTION_EPS` | 0.001 | Minimum meaningful infection level |
| `MICROBIOME_MAJORITY_THRESHOLD` | 0.5 | Threshold for minority→majority resistance promotion |
| `SIMULATION_START_YEAR` | 1930.0 | Calendar year at day 0 |
| `POLICY_BRANCH_YEAR` | 2027.0 | Year policies diverge |
| `MAX_MECHANISM_PROFILES` | 200 | Reservoir sample size per bacteria |
| `REGION_COUNT` | 6 | Number of world regions |
| `BACTERIA_COUNT` | 42 | Number of bacterial species |
| `DrugClass::NUM_CLASSES` | 18 | Number of drug classes |

## Output Format

Each run produces a single CSV file: `amr_simulation_output_analysis_outputs/simulation_summary_NNNNNN.csv`

Each row is one timestep (day). Columns include:

- **Scalar**: `day`, `year`, `total_alive`, `total_infected`, `total_on_treatment`, `total_in_hospital`, `total_sepsis`, `total_died_infection`, `total_died_sepsis`, `total_died_background`, `total_died_toxicity`, `total_new_infections`, `drug_stops_due_to_toxicity`, `policy_name`, ...
- **Per-bacteria** (~42 each): `{bacteria}_infected`, `{bacteria}_carriers`, `{bacteria}_deaths`, `{bacteria}_new_infections`, ...
- **Per-drug** (~58 each): `{drug}_prescribed`, `{drug}_active_treatments`, ...
- **Per-bacteria×drug** (~2,436 each): `{bacteria}_{drug}_activity_r`, `{bacteria}_{drug}_majority_r`, ...
- **Per-region**: `{region}_infected`, `{region}_hospitalized`, ...

The CSV is wide-format with hundreds of columns per row.

## Analysis Tools

### Python Setup

```bash
pip install -r requirements.txt
```

Dependencies: `pandas`, `matplotlib`, `seaborn`, `numpy`, `scipy`, `polars`, `pyarrow`

### Running Analysis

```python
from amr_simulation_output_analysis.amr_analysis import run_analysis
run_analysis("amr_simulation_output_analysis_outputs/simulation_summary_NNNNNN.csv")
```

Key analysis modules:

| Module | Purpose |
|--------|---------|
| `amr_analysis.py` | Main orchestrator |
| `calibration_summary.py` | Compare outputs to calibration targets |
| `multi_run_activity_r_plot.py` | Overlay resistance trends across runs |
| `plotting/base_plots.py` | Time series, prevalence plots |
| `plotting/detail_plots.py` | Per-bacteria/drug detailed views |
| `plotting/grouped_plots.py` | Drug-class and region aggregations |
| `empirical/` | Load and normalise empirical AMR data for comparison |

## Dependencies

| Crate | Version | Purpose |
|-------|---------|---------|
| `rand` | 0.8 | Random number generation (SmallRng) |
| `rayon` | 1.7 | Data-parallel iteration |
| `lazy_static` | 1.4 | Global parameter initialisation |
| `chrono` | 0.4 | Timestamp formatting |
| `csv` | 1.3 | CSV export |
| `serde` | 1.0 | Serialisation for branch checkpointing |
| `bincode` | 1.3 | Binary serialisation |
| `log` | 0.4 | Logging framework |
| `env_logger` | 0.11 | Log output to stderr |
| `plotters` | 0.3 | Optional plotting (feature-gated behind `plots`) |

## License

Not yet specified.
