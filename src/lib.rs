//! Library crate for the AMR simulation.
//!
//! Most executable behaviour lives in three top-level modules:
//! - `config`: parameter defaults, typed readers, and lookup helpers
//! - `rules`: per-timestep state transitions for individuals and pathogens
//! - `simulation`: population data structures, orchestration, summaries, and logging

pub mod config;
pub mod rules;
pub mod simulation;