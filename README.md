
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

Currently the model considers 21 bacteria (the ones used in the Global Burden of Disease work on AMR) and 41 antibiotics (so 861 bacteria-drug combinations) but this can be expanded.

Other variables include whether the person is hospitalized, with consequences for the range of bacteria exposed to.

Mortality risk is separated by (i) background mortality risk, which is age and region-specific (noting that region can be re-coded such that home is a given single country); (ii) mortality risk given sepsis (and possibly according to infection in a person with severe immunosuppression) and (iii) mortality risk specifically due to adverse antibiotic drug effects.


# Rust AMR Simulation Model

This project simulates the spread and treatment of antimicrobial resistance (AMR) in a population using an agent-based model written in Rust.

## Model Overview

- **Population:** Each individual is modeled with attributes such as age, sex, region, infection status, resistance levels, drug usage, and more.
- **Bacteria & Drugs:** The model supports multiple bacteria and drug types, with customizable parameters for acquisition, resistance, and treatment.
- **Time Steps:** The simulation proceeds in discrete time steps, updating each individual's state according to biological and treatment rules.

## Main Components

- `simulation/`  
  - `simulation.rs`: Orchestrates the simulation, manages the population, and aggregates statistics.
  - `population.rs`: Defines the `Individual` and `Population` structs and their attributes.
- `rules/`  
  - `mod.rs`: Contains the core logic for updating individuals each time step (infection, resistance, drug initiation, etc.).
- `config/`  
  - `config.rs`: Stores all model parameters and provides lookup functions.

## Key Processes

- **Infection Acquisition:** Probability-based, influenced by contact/exposure, vaccination, microbiome, and region.
- **Resistance Emergence:** Modeled for both infection site and microbiome, with parameters for baseline and drug/bacteria-specific rates.
- **Drug Initiation/Stopping:** Drugs are started/stopped based on infection status, resistance test results, and other rules.
- **Testing:** Simulates lab identification of bacteria and resistance, with delays and error probabilities.
- **Hospitalization & Travel:** Individuals may be hospitalized or travel between regions, affecting exposure and risk.
- **Mortality:** Death risk is calculated from background, sepsis, and drug toxicity.

## Running the Simulation

```sh
cargo run --release
```

population size and number of time steps is detrmined in `main.rs`:

```rust
let population_size = 100_000;
let time_steps = 20;
```

## Example Output

The simulation prints detailed information for individual 0 at each time step, as well as summary statistics.

## Customization

- **Parameters:**  
  Edit `src/config.rs` to change model parameters (infection rates, resistance emergence, drug properties, etc.).
- **Bacteria/Drugs:**  
  Add or remove bacteria/drugs by editing the lists in `population.rs` and updating parameters in `config.rs`.
- **Rules:**  
  Modify `src/rules/mod.rs` to change the biological or treatment logic.


This project is for research and educational use. 

