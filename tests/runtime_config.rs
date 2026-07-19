use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::fs;
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::time::{SystemTime, UNIX_EPOCH};

fn test_root(label: &str) -> PathBuf {
    let unique = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock must be after Unix epoch")
        .as_nanos();
    std::env::temp_dir().join(format!(
        "amr_runtime_config_{label}_{}_{}",
        std::process::id(),
        unique
    ))
}

fn write_config(root: &Path, seed: Option<u64>) -> (PathBuf, PathBuf) {
    let output_dir = root.join("model-outputs");
    fs::create_dir_all(root).expect("test run directory must be creatable");

    let mut overrides = json!({
        "population_size": 64,
        "time_steps": 4,
        "calibration_mode": "partial",
        "config_validation_mode": "strict"
    });
    if let Some(value) = seed {
        overrides["rng_seed"] = json!(value);
    }
    let config = json!({
        "schema_version": "amr-run-config/v1",
        "contract": "typed_adapter_runtime_config",
        "override_semantics": "explicit_keys_only_model_defaults_when_omitted",
        "overrides": overrides,
        "outputs": {
            "run_root": root,
            "report_json": root.join("amr_report.json"),
            "progress_json": root.join("amr_progress.json"),
            "progress_jsonl": root.join("amr_progress.jsonl"),
            "model_output_dir": output_dir
        }
    });
    let config_path = root.join("amr_run_config.json");
    fs::write(
        &config_path,
        serde_json::to_vec_pretty(&config).expect("test config must serialize"),
    )
    .expect("test config must be writable");
    (config_path, output_dir)
}

fn run_binary(root: &Path, config_path: &Path, environment_seed: Option<u64>) -> Output {
    let mut command = Command::new(env!("CARGO_BIN_EXE_executable_amr"));
    command
        .current_dir(root)
        .env("AMR_RUN_CONFIG", config_path)
        .env_remove("AMR_RNG_SEED")
        .env("AMR_CONFIG_VALIDATION", "strict")
        .env("AMR_SOURCE_HASH", "runtime-config-integration-test")
        .env("RAYON_NUM_THREADS", "2");
    if let Some(seed) = environment_seed {
        command.env("AMR_RNG_SEED", seed.to_string());
    }
    command.output().expect("AMR binary must start")
}

fn assert_success(output: &Output) -> (String, String) {
    let stdout = String::from_utf8_lossy(&output.stdout).into_owned();
    let stderr = String::from_utf8_lossy(&output.stderr).into_owned();
    assert!(
        output.status.success(),
        "AMR binary failed\nstdout:\n{stdout}\nstderr:\n{stderr}"
    );
    (stdout, stderr)
}

fn output_file_with_prefix(output_dir: &Path, prefix: &str) -> PathBuf {
    fs::read_dir(output_dir)
        .expect("configured model output directory must exist")
        .filter_map(Result::ok)
        .find(|entry| entry.file_name().to_string_lossy().starts_with(prefix))
        .unwrap_or_else(|| panic!("output with prefix {prefix} must exist"))
        .path()
}

fn read_json(path: &Path) -> Value {
    serde_json::from_slice(&fs::read(path).expect("JSON evidence must be readable"))
        .expect("JSON evidence must parse")
}

fn sha256_file(path: &Path) -> String {
    format!(
        "{:x}",
        Sha256::digest(fs::read(path).expect("hashed test file must be readable"))
    )
}

#[test]
fn binary_consumes_typed_runtime_config_without_source_rewrite() {
    let root = test_root("explicit_seed");
    let (config_path, output_dir) = write_config(&root, Some(1729));
    let output = run_binary(&root, &config_path, Some(1729));
    let (stdout, stderr) = assert_success(&output);

    assert!(stdout.contains("last_timestep=3"));
    assert!(stderr.contains("source=runtime_config"));
    assert!(stderr.contains("AMR_RNG_SEED=1729 source=runtime_config"));
    assert!(output_file_with_prefix(&output_dir, "simulation_summary_").is_file());

    let metadata = fs::read_to_string(output_file_with_prefix(&output_dir, "run_metadata_"))
        .expect("run metadata must be readable");
    for expected in [
        "status=completed",
        "rng_seed=1729",
        "rng_seed_source=runtime_config",
        "population_size=64",
        "time_steps=4",
        "calibration_mode=Partial",
    ] {
        assert!(metadata.contains(expected), "missing metadata: {expected}");
    }

    let progress = read_json(&root.join("amr_progress.json"));
    assert_eq!(progress["schema_version"], "amr-progress/v1");
    assert_eq!(progress["status"], "completed");
    assert_eq!(progress["current_step"], 4);
    assert_eq!(progress["total_steps"], 4);
    assert_eq!(progress["source"], "amr_model");
    let progress_events = fs::read_to_string(root.join("amr_progress.jsonl"))
        .expect("progress event history must be readable");
    assert!(progress_events.lines().count() >= 3);

    let report = read_json(&root.join("amr_report.json"));
    assert_eq!(report["schema_version"], "amr-report/v1");
    assert_eq!(report["status"], "completed");
    assert_eq!(report["exit_code"], 0);
    assert_eq!(report["runtime_config_consumed"], true);
    assert_eq!(report["runtime_config_sha256"], sha256_file(&config_path));
    assert_eq!(report["effective_config"]["rng_seed"]["value"], 1729);
    assert_eq!(
        report["effective_config"]["rng_seed"]["source"],
        "runtime_config"
    );
    assert!(report["summary_sha256"].as_str().is_some());

    fs::remove_dir_all(&root).expect("test run directory must be removable");
}

#[test]
fn binary_generates_its_seed_when_no_override_is_supplied() {
    let root = test_root("generated_seed");
    let (config_path, output_dir) = write_config(&root, None);
    let output = run_binary(&root, &config_path, None);
    let (_, stderr) = assert_success(&output);

    assert!(stderr.contains("source=generated"));
    let metadata = fs::read_to_string(output_file_with_prefix(&output_dir, "run_metadata_"))
        .expect("run metadata must be readable");
    assert!(metadata.contains("rng_seed_source=generated"));

    let report = read_json(&root.join("amr_report.json"));
    assert_eq!(report["runtime_config_consumed"], true);
    assert_eq!(report["runtime_config_sha256"], sha256_file(&config_path));
    assert_eq!(
        report["effective_config"]["rng_seed"]["source"],
        "generated"
    );
    assert!(report["effective_config"]["rng_seed"]["value"]
        .as_u64()
        .is_some());

    fs::remove_dir_all(&root).expect("test run directory must be removable");
}
