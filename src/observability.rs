use std::process::Command;
use std::sync::atomic::{AtomicU32, AtomicUsize, Ordering};

const UNKNOWN_RUN_ID: u32 = 0;
const UNKNOWN_TIMESTEP: usize = usize::MAX;

static CURRENT_RUN_ID: AtomicU32 = AtomicU32::new(UNKNOWN_RUN_ID);
static CURRENT_TIMESTEP: AtomicUsize = AtomicUsize::new(UNKNOWN_TIMESTEP);

pub fn clear_run_context() {
    CURRENT_RUN_ID.store(UNKNOWN_RUN_ID, Ordering::Relaxed);
    CURRENT_TIMESTEP.store(UNKNOWN_TIMESTEP, Ordering::Relaxed);
}

pub fn set_current_run_id(run_id: u32) {
    CURRENT_RUN_ID.store(run_id, Ordering::Relaxed);
}

pub fn current_run_id() -> Option<u32> {
    let run_id = CURRENT_RUN_ID.load(Ordering::Relaxed);
    if run_id == UNKNOWN_RUN_ID {
        None
    } else {
        Some(run_id)
    }
}

pub fn set_current_timestep(timestep: usize) {
    CURRENT_TIMESTEP.store(timestep, Ordering::Relaxed);
}

pub fn current_timestep() -> Option<usize> {
    let timestep = CURRENT_TIMESTEP.load(Ordering::Relaxed);
    if timestep == UNKNOWN_TIMESTEP {
        None
    } else {
        Some(timestep)
    }
}

pub fn resolve_source_hash() -> String {
    if let Ok(value) = std::env::var("AMR_SOURCE_HASH") {
        let trimmed = value.trim();
        if !trimmed.is_empty() {
            return trimmed.to_string();
        }
    }

    if let Ok(value) = std::fs::read_to_string("source_hash.txt") {
        let trimmed = value.trim();
        if !trimmed.is_empty() {
            return trimmed.to_string();
        }
    }

    if let Ok(output) = Command::new("git").args(["rev-parse", "HEAD"]).output() {
        if output.status.success() {
            if let Ok(text) = String::from_utf8(output.stdout) {
                let trimmed = text.trim();
                if !trimmed.is_empty() {
                    return trimmed.to_string();
                }
            }
        }
    }

    "unknown".to_string()
}
