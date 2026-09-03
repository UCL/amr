// Simulation-layer modules.
//
// `population` holds the persistent state for individuals and enums/constants shared
// across the model. `simulation` orchestrates timestep execution and summary export.
// `journey_logger` is an opt-in diagnostic path for detailed infection trajectories.
// `rng` provides seed-derived domain-separated streams and an entropy-backed fallback.
pub mod journey_logger;
pub mod population;
pub mod rng;
pub mod simulation;
