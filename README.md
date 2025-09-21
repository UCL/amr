
We are developing a stochastic individual-based model for anti-bacterial resistance .  This aims to model populations of individuals being tracked over time for their infection and drug resistance status.  

Note that the variable any_r takes the value 0 for the bacteria/drug combination either if the infection is not present or, if the infection is present and within the bacteria causing the infection there are none with any resistance to the drug.  If there are any bacteria within the population of bacteria causing the infection which have some level of resistance to the drug then the variable takes a non-zero value between 0 and 1.  This value indicates the level of resistance that the bacteria has, with 1 being total resistance meaning the drug has no effect at all.  The variable majority_r is the same as any_r except that it indicates whether the resistant bacteria form the majority, so it only takes non-zero values when this is the case.  When majority_r and any_r have non-zero values these values are always equal.

For any drug that is being taken or has recently been take and some non-zero drug level remains there is a variable indicating the current activity level of the drug against the bacteria.  This is called activity_r.

activity_r depends on the underlying potency of the drug against the bacteria (were no resistance present), any_r and the current level of the drug.  When a person has presence of more than one drug in their system then the activity_r of the two drugs is summed to get the total activity, although there is capacity for introducing interactions between effects of drugs.

Another variable relating to resistance is microbiome_r which indicates whether the person is carrying bacteria with resistance to the given drug.  

The time step is daily and initially we aim to simulate from the date of the introduction of penicillin in 1942 to the present time.  However, the model can be used with a later time zero than this with pre-existing levels of resistance included at time 0.

Age is conveyed in days and to account for children born after time 0 we give some people a negative age.  Nothing happens to them in the simulation until they are born, reaching age 0.

Infection risk is determined by directly specifying these risks by age, region and calendar time.  We specify whether a person is community acquired from another person, acquired from the environment, or acquired in hospital.  The infectious syndrome (site of infection) is assigned at the time of infection with a given bacteria.

There are a number of variables relating to exposure level: sexual_contact_level, airborne_contact_level_with_adults, airborne_contact_level_with_children, oral_exposure_level, mosquito_exposure_level.  Alongisde age and region, these are used as multipliers of the risk of acquisition of a given bacteria from another person in the population. 

If a person is infected directly or indirectly from another person then we randomly sample from people in the same region who have the bacterial infection to assign the value of any_r for the newly infected person.

The level of any antibiotic given is conveyed on a standardized scale of 0-10, with 10 being the daily level of drug on days in which it is taken/administered using the standard age-specific dose.  After stopping the drug the drug level in the persons system decays.  If a double dose is given then the maximum level is 20, etc.  

We account for the fact that testing to identify the bacteria and to identify levels of drug resistance of drugs to the bacteria may take place.  We have variables that indicate whether the bacteria has been identified in a test, narrowing the range of likely antibiotics used, and, separately, whether a test has been done to assess levels of resistance, making choice of a drug to which the bacteria is sensitive more likely. 

We have a variable for the current level of immunity the person carries to each specific bacteria.  If infected with a bacteria the level of immunity grows dependent on bacteria level and days since infection.  We do not have a variable indicating whether the person is currently severely immunosuppressed but we can add this.  If added, this will determine both risk of death if infected and the level of immunity to each bacteria. 

People live in a certain region, but they may visit other regions.  Currently the regions are broadly aligned with continents but there could be flexibility over this.

Currently the model considers 30 bacteria (the ones used in the Global Burden of Disease work on AMR) and 42 antibiotics (so 30 x 42 bacteria-drug combinations) but this can be expanded.

Other variables include whether the person is hospitalized, with consequences for the range of bacteria exposed to.

Mortality risk is separated by (i) background mortality risk, which is age and region-specific (noting that region can be re-coded such that home is a given single country); (ii) mortality risk given sepsis (and possibly according to infection in a person with severe immunosuppression) and (iii) mortality risk specifically due to adverse antibiotic drug effects.


# AMR Simulation and Analysis Project

This project simulates the spread and treatment of antimicrobial resistance (AMR) in a population using an agent-based model written in Rust, with comprehensive Python-based analysis and visualization tools.

## Project Structure

### Rust Simulation Model (`src/`)
The core simulation engine written in Rust for high-performance modeling of large populations.

**Key Components:**
- `simulation/`: Orchestrates the simulation, manages populations, and aggregates statistics
- `population.rs`: Defines Individual and Population structs with comprehensive attributes
- `rules/`: Core logic for updating individuals each time step (infection, resistance, drug initiation)
- `config/`: Model parameters and lookup functions

### Python Analysis System (`amr_simulation_output_analysis/`)
Modular analysis and visualization system for processing simulation outputs.

**Architecture:**
```
amr_simulation_output_analysis/
├── __init__.py
├── config.py              # Configuration management with granular plot controls
├── data_loader.py         # Centralized data loading and caching
├── utils.py              # Common utilities and helper functions
├── plotting/
│   ├── __init__.py
│   ├── grouped_plots.py   # Main grouped figures (Figures 1-9)
│   ├── detail_plots.py    # Individual detailed visualizations
│   └── utils.py          # Plotting utilities and styling
└── empirical/
    ├── __init__.py
    ├── data_config.py     # Empirical data source configuration
    └── parsers.py         # Data parsing and integration functions
```

## Quick Start

### Running the Simulation
```bash
cargo run --release
```

Configure population size and time steps in `main.rs`:
```rust
let population_size = 100_000;
let time_steps = 20;
```

### Analyzing Results
```python
from amr_simulation_output_analysis import create_all_plots, PlotConfig

# Generate all plots with default configuration
create_all_plots()

# Generate specific plot categories
config = PlotConfig(
    grouped_plots=True,
    detail_plots=False,
    age_specific_plots=True,
    drug_failure_plots=True
)
create_all_plots(config)
```

For complete usage examples, see `amr_analysis_examples.py`.

## Analysis Capabilities

### Grouped Visualizations (Figures 1-9)
- **Figure 1**: Infection resolution patterns by bacteria
- **Figure 2**: Mean activity R analysis by bacteria
- **Figure 3**: Day-7 drug initiation patterns
- **Figure 4**: Syndrome distribution patterns  
- **Figure 5**: Population dynamics tracking
- **Figure 6**: Drug initiation timing analysis
- **Figure 7**: Regional syndrome patterns
- **Figure 8**: Comprehensive resistance tracking
- **Figure 9**: Treatment outcome analysis

### Detailed Analysis Categories
- **Age-specific death rates** (7 visualizations)
- **Regional analysis** (age distribution, death rates, drug failure)
- **Bacteria-specific analysis** (death rates, drug failure, infection resolution)
- **Drug analysis** (usage patterns, resistance emergence, MIC analysis)
- **Resistance mechanisms** (source tracking, emergence patterns)
- **Population health metrics** (mortality, syndrome distribution)

### Configuration Options
The system provides granular control over plot generation:

```python
config = PlotConfig(
    # Main categories
    grouped_plots=True,
    detail_plots=True,
    
    # Specific plot types
    drug_failure_rate_by_bacteria_region=True,
    mean_mic_by_drug_for_each_bacteria=True,
    incidence_of_infection_hospital=True,
    proportion_of_people_taking_each_drug=True,
    # ... 20+ additional granular controls
)
```

## Model Features

### Biological Processes
- **Infection Acquisition**: Contact-based transmission with regional, age, and exposure factors
- **Resistance Emergence**: Dynamic resistance development in infections and microbiomes
- **Drug Pharmacodynamics**: Multi-drug interactions with realistic decay and efficacy
- **Testing and Diagnosis**: Laboratory delays, identification accuracy, resistance testing
- **Immunity Development**: Bacteria-specific immunity acquisition and maintenance

### Population Dynamics
- **Demographics**: Age, sex, region with realistic population structures
- **Healthcare**: Hospitalization, healthcare access, treatment protocols
- **Mobility**: Inter-regional travel affecting exposure and transmission
- **Mortality**: Background, sepsis-related, and drug toxicity mortality

### Drug and Resistance Modeling
- **Multi-drug Support**: 357 drugs across major antibiotic classes
- **Resistance Mechanisms**: Genetic and phenotypic resistance evolution
- **Treatment Protocols**: Empirical and targeted therapy based on testing
- **Microbiome Effects**: Commensal bacteria resistance affecting treatment

## Data Integration

The analysis system integrates multiple empirical data sources:
- **CDC AR Threats 2019**: US antimicrobial resistance surveillance
- **ECDC Data 2023**: European resistance and consumption data
- **IQVIA Sales 2023**: Global antibiotic usage patterns
- **WHO Global Data**: International resistance surveillance
- **GBD Study Data**: Global burden of disease metrics

## Output and Visualization

The system generates comprehensive visualizations:
- **2,000+ individual plots** across all analysis categories
- **Multi-format support**: PNG, PDF, SVG output options
- **Publication-ready figures** with consistent styling
- **Interactive dashboards** (planned enhancement)
- **Statistical summaries** with empirical data overlays

## Development and Customization

### Adding New Analysis
1. Implement new functions in `detail_plots.py`
2. Add configuration flags in `config.py`
3. Update plot dispatcher in `detail_plots.py`
4. Add empirical data sources as needed

### Extending the Simulation
- **Parameters**: Edit `src/config.rs` for model parameters
- **Bacteria/Drugs**: Modify lists in `population.rs` and update `config.rs`
- **Rules**: Update `src/rules/mod.rs` for biological or treatment logic
- **Output**: Extend CSV output in simulation for new metrics

## Archive and Legacy
Previous monolithic analysis scripts have been archived in `archive/legacy_scripts/`:
- `analyze_simulation.py` (6,623-line original script)
- `empirical_enhancement.py` (legacy enhancement tools)
- Development and migration scripts in `archive/development_scripts/`

## Requirements

**Rust**: Latest stable version
**Python**: 3.8+ with packages listed in `requirements.txt`

This project is for research and educational use in antimicrobial resistance modeling. 

