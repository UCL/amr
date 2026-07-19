use chrono::{SecondsFormat, Utc};
use serde::Serialize;
use std::fs::{self, OpenOptions};
use std::io::{self, Write};
use std::path::{Path, PathBuf};
use std::process::Command;
use std::sync::atomic::{AtomicBool, AtomicU32, AtomicUsize, Ordering};
use std::sync::OnceLock;

const UNKNOWN_RUN_ID: u32 = 0;
const UNKNOWN_TIMESTEP: usize = usize::MAX;

static CURRENT_RUN_ID: AtomicU32 = AtomicU32::new(UNKNOWN_RUN_ID);
static CURRENT_TIMESTEP: AtomicUsize = AtomicUsize::new(UNKNOWN_TIMESTEP);
static LAST_PROGRESS_TIMESTEP: AtomicUsize = AtomicUsize::new(UNKNOWN_TIMESTEP);
static PROGRESS_WRITE_ERROR_REPORTED: AtomicBool = AtomicBool::new(false);
static PROGRESS_SINK: OnceLock<ProgressConfig> = OnceLock::new();

const PROGRESS_STEP_INTERVAL: usize = 100;

#[derive(Debug, Clone)]
pub struct ProgressConfig {
    pub snapshot_path: PathBuf,
    pub events_path: PathBuf,
    pub population_size: usize,
    pub total_steps: usize,
    pub calibration_mode: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProgressStatus {
    Starting,
    Running,
    Completed,
    Failed,
}

impl ProgressStatus {
    const fn as_str(self) -> &'static str {
        match self {
            Self::Starting => "starting",
            Self::Running => "running",
            Self::Completed => "completed",
            Self::Failed => "failed",
        }
    }
}

#[derive(Serialize)]
struct ProgressSnapshot<'a> {
    schema_version: &'static str,
    status: &'static str,
    phase: &'static str,
    run_id: Option<u32>,
    current_step: Option<usize>,
    total_steps: usize,
    step_unit: &'static str,
    progress_label: Option<String>,
    population_size: usize,
    calibration_mode: &'a str,
    updated_at_utc: String,
    source: &'static str,
}

pub fn configure_progress(config: ProgressConfig) -> io::Result<()> {
    ensure_parent(&config.snapshot_path)?;
    ensure_parent(&config.events_path)?;
    PROGRESS_SINK.set(config).map_err(|_| {
        io::Error::new(
            io::ErrorKind::AlreadyExists,
            "AMR progress output is already configured",
        )
    })?;
    publish_progress(ProgressStatus::Starting, None)
}

pub fn publish_progress(status: ProgressStatus, current_step: Option<usize>) -> io::Result<()> {
    let Some(config) = PROGRESS_SINK.get() else {
        return Ok(());
    };
    let snapshot = ProgressSnapshot {
        schema_version: "amr-progress/v1",
        status: status.as_str(),
        phase: status.as_str(),
        run_id: current_run_id(),
        current_step,
        total_steps: config.total_steps,
        step_unit: "time_step",
        progress_label: current_step
            .map(|step| format!("time step {step} / {}", config.total_steps)),
        population_size: config.population_size,
        calibration_mode: &config.calibration_mode,
        updated_at_utc: Utc::now().to_rfc3339_opts(SecondsFormat::Millis, true),
        source: "amr_model",
    };
    write_json_atomically(&config.snapshot_path, &snapshot)?;
    append_json_line(&config.events_path, &snapshot)
}

pub fn clear_run_context() {
    CURRENT_RUN_ID.store(UNKNOWN_RUN_ID, Ordering::Relaxed);
    CURRENT_TIMESTEP.store(UNKNOWN_TIMESTEP, Ordering::Relaxed);
    LAST_PROGRESS_TIMESTEP.store(UNKNOWN_TIMESTEP, Ordering::Relaxed);
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
    if timestep.checked_rem(PROGRESS_STEP_INTERVAL) == Some(0)
        && LAST_PROGRESS_TIMESTEP.load(Ordering::Relaxed) != timestep
    {
        match publish_progress(ProgressStatus::Running, Some(timestep)) {
            Ok(()) => {
                LAST_PROGRESS_TIMESTEP.store(timestep, Ordering::Relaxed);
            }
            Err(error) => report_progress_write_error(&error),
        }
    }
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

    if let Some(head) = git_output(["rev-parse", "HEAD"]) {
        let dirty = git_output(["status", "--porcelain"]).is_some_and(|status| !status.is_empty());
        if dirty {
            return format!("{head}-dirty");
        }

        return head;
    }

    "unknown".to_string()
}

pub fn write_json_atomically<T: Serialize>(path: &Path, value: &T) -> io::Result<()> {
    ensure_parent(path)?;
    let filename = path
        .file_name()
        .and_then(|name| name.to_str())
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "JSON output path must name a file",
            )
        })?;
    let temporary_path = path.with_file_name(format!(".{filename}.tmp-{}", std::process::id()));
    let mut encoded = serde_json::to_vec(value)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
    encoded.push(b'\n');
    fs::write(&temporary_path, encoded)?;
    fs::rename(&temporary_path, path).inspect_err(|_| {
        let _ = fs::remove_file(&temporary_path);
    })
}

fn append_json_line<T: Serialize>(path: &Path, value: &T) -> io::Result<()> {
    let mut file = OpenOptions::new().create(true).append(true).open(path)?;
    serde_json::to_writer(&mut file, value)
        .map_err(|error| io::Error::new(io::ErrorKind::InvalidData, error))?;
    file.write_all(b"\n")?;
    file.flush()
}

fn ensure_parent(path: &Path) -> io::Result<()> {
    let parent = path.parent().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "output path has no parent directory",
        )
    })?;
    if !parent.as_os_str().is_empty() {
        fs::create_dir_all(parent)?;
    }
    Ok(())
}

fn report_progress_write_error(error: &io::Error) {
    if !PROGRESS_WRITE_ERROR_REPORTED.swap(true, Ordering::Relaxed) {
        eprintln!("[progress] unable to update AMR progress evidence: {error}");
    }
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
