# AMR Simulation Technical Documentation

This directory contains detailed technical documentation for developers new to the AMR (Antimicrobial Resistance) simulation codebase.

## Documentation Structure

The documentation is organized into logical groups corresponding to major subsystems:

| Document | Description |
|----------|-------------|
| [01_individual_state.md](01_individual_state.md) | Per-person state variables in the `Individual` struct |
| [02_resistance_system.md](02_resistance_system.md) | Resistance levels, mechanisms, and how they evolve |
| [03_drug_treatment.md](03_drug_treatment.md) | Drug selection, pharmacokinetics, and treatment flow |
| [04_infection_dynamics.md](04_infection_dynamics.md) | Infection acquisition, progression, and resolution |
| [05_microbiome_system.md](05_microbiome_system.md) | Microbiome carriage and its role in resistance |
| [06_config_parameters.md](06_config_parameters.md) | Configuration system and parameter categories |
| [07_simulation_flow.md](07_simulation_flow.md) | Main simulation loop and timestep processing |
| [08_enums_constants.md](08_enums_constants.md) | Enumerations, constants, and static data |

## Quick Navigation by Task

### "I want to understand how resistance works"
→ Start with [02_resistance_system.md](02_resistance_system.md)

### "I want to add a new drug"
→ See [03_drug_treatment.md](03_drug_treatment.md) and [06_config_parameters.md](06_config_parameters.md)

### "I want to add a new bacteria"
→ See [08_enums_constants.md](08_enums_constants.md) (BACTERIA_LIST) and [06_config_parameters.md](06_config_parameters.md)

### "I want to understand the main loop"
→ Start with [07_simulation_flow.md](07_simulation_flow.md)

### "I want to modify infection acquisition"
→ See [04_infection_dynamics.md](04_infection_dynamics.md)

## Key Source Files

| File | Purpose |
|------|---------|
| `src/simulation/population.rs` | Data structures: `Individual`, enums, constants |
| `src/simulation/simulation.rs` | Main simulation loop, caching, aggregation |
| `src/rules/mod.rs` | Core update logic (5000+ lines of state transitions) |
| `src/config.rs` | All configurable parameters (13000+ lines) |
| `src/main.rs` | Entry point, CLI, output handling |

## Architecture Overview

```
┌─────────────────────────────────────────────────────────────────┐
│                         main.rs                                  │
│  - Parse CLI arguments                                          │
│  - Initialize population                                        │
│  - Run simulation loop                                          │
│  - Write outputs                                                │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                    simulation/simulation.rs                      │
│  - run_simulation(): main time-stepping loop                    │
│  - MajorityRCache: resistance prevalence sampling               │
│  - Aggregation and reporting                                    │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                       rules/mod.rs                               │
│  - apply_rules(): per-individual state updates                  │
│  - Infection acquisition, progression, clearance                │
│  - Drug selection, effects, toxicity                            │
│  - Resistance emergence, HGT, reversion                         │
└─────────────────────────────────────────────────────────────────┘
                              │
                              ▼
┌─────────────────────────────────────────────────────────────────┐
│                         config.rs                                │
│  - PARAMETERS: HashMap of all simulation parameters             │
│  - Reader structs: type-safe parameter access                   │
│  - Defaults organized by category (A-K)                         │
└─────────────────────────────────────────────────────────────────┘
```

## Simulation Time

- **Start Year**: 1930 (day 0)
- **Time Step**: 1 day
- **Typical Run**: 105 years (1930-2035) = 38,325 days
- **Age Convention**: Negative age = not yet born (birth at age 0)

## Key Concepts

### Resistance Levels
- `any_r`: Resistance present in any fraction of bacteria (0.0-1.0)
- `majority_r`: Resistance in >50% of bacteria (always ≤ any_r)
- `microbiome_r`: Resistance level in gut/respiratory flora

### Drug Levels
- Standard dose = level 10.0 per day
- Decays with drug-specific half-life
- Efficacy = potency × level × (1 - resistance)

### Infection States
- `level[b] > INFECTION_EPS`: Active infection with bacteria b
- `presence_microbiome[b]`: Colonization without active infection
- `sepsis[b]`: Life-threatening complication of infection

