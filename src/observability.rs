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

/// Emit Linux process-memory counters at infrequent lifecycle boundaries.
///
/// Values come from procfs and are intentionally unavailable on other platforms.
/// Runner-level cgroup metrics remain the authority for job memory limits.
pub fn log_process_memory(label: &str) {
    #[cfg(target_os = "linux")]
    {
        let rss_kib = proc_memory_value("/proc/self/status", "VmRSS:");
        let peak_rss_kib = proc_memory_value("/proc/self/status", "VmHWM:");
        let pss_kib = proc_memory_value("/proc/self/smaps_rollup", "Pss:");
        println!(
            "Process memory [{}]: rss_kib={}, pss_kib={}, peak_rss_kib={}",
            label,
            format_optional_counter(rss_kib),
            format_optional_counter(pss_kib),
            format_optional_counter(peak_rss_kib)
        );
    }

    #[cfg(not(target_os = "linux"))]
    let _ = label;
}

#[cfg(target_os = "linux")]
fn proc_memory_value(path: &str, key: &str) -> Option<u64> {
    let contents = std::fs::read_to_string(path).ok()?;
    let line = contents.lines().find(|line| line.starts_with(key))?;
    line.split_whitespace().nth(1)?.parse().ok()
}

#[cfg(target_os = "linux")]
fn format_optional_counter(value: Option<u64>) -> String {
    value
        .map(|counter| counter.to_string())
        .unwrap_or_else(|| "unavailable".to_string())
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

    if let Some(head) = git_output(["rev-parse", "HEAD"]) {
        let dirty = git_output(["status", "--porcelain"]).is_some_and(|status| !status.is_empty());
        if dirty {
            return format!("{head}-dirty");
        }

        return head;
    }

    "unknown".to_string()
}

fn git_output<const N: usize>(args: [&str; N]) -> Option<String> {
    let output = Command::new("git").args(args).output().ok()?;
    if !output.status.success() {
        return None;
    }

    let text = String::from_utf8(output.stdout).ok()?;
    let trimmed = text.trim();
    if trimmed.is_empty() {
        None
    } else {
        Some(trimmed.to_string())
    }
}
