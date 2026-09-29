//! Publish completed output files without replacing artifacts from an earlier run.

use std::ffi::OsString;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

/// Publish a complete staged summary under the requested name or a free repeat name.
///
/// Each hard-link creation atomically reserves its destination without replacing an
/// existing file. The staged file must be on the same filesystem as the destination;
/// unsupported hard links and other publication failures leave it available for recovery.
/// Existing directories and symlinks are errors rather than filename reservations.
pub fn publish_summary_no_clobber(staged: &Path, requested: &Path) -> io::Result<PathBuf> {
    let mut destination = requested.to_path_buf();
    let mut repeat = 0_u64;

    loop {
        match fs::hard_link(staged, &destination) {
            Ok(()) => {
                // Publication has succeeded. Failure to remove this second link must not
                // misreport or undo the successfully published, complete destination.
                if let Err(error) = fs::remove_file(staged) {
                    eprintln!(
                        "Warning: summary published to {} but unable to remove staged file {}: {}",
                        destination.display(),
                        staged.display(),
                        error
                    );
                }
                return Ok(destination);
            }
            Err(error) => {
                if error.kind() != io::ErrorKind::AlreadyExists {
                    return Err(error);
                }

                // Inspect only after the atomic create reports a collision. Do not follow
                // symlinks, and retain the original link error if the entry cannot be read.
                match fs::symlink_metadata(&destination) {
                    Ok(metadata) if metadata.file_type().is_file() => {}
                    _ => return Err(error),
                }
                repeat = match repeat.checked_add(1) {
                    Some(next) => next,
                    None => return Err(error),
                };
                destination = repeat_path(requested, repeat)?;
            }
        }
    }
}

fn repeat_path(requested: &Path, repeat: u64) -> io::Result<PathBuf> {
    let stem = requested.file_stem().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "summary output path must have a filename",
        )
    })?;
    let mut filename = OsString::from(stem);
    filename.push(format!("_repeat_{repeat}"));
    if let Some(extension) = requested.extension() {
        filename.push(".");
        filename.push(extension);
    }
    Ok(requested.with_file_name(filename))
}

#[cfg(test)]
mod tests {
    use super::publish_summary_no_clobber;
    use std::collections::HashSet;
    use std::fs;
    use std::io;
    use std::path::{Path, PathBuf};
    use std::sync::atomic::{AtomicU64, Ordering};
    use std::sync::{Arc, Barrier};
    use std::time::{SystemTime, UNIX_EPOCH};

    static NEXT_DIRECTORY: AtomicU64 = AtomicU64::new(0);

    struct TestDirectory(PathBuf);

    impl TestDirectory {
        fn new() -> Self {
            let timestamp = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .unwrap()
                .as_nanos();
            let counter = NEXT_DIRECTORY.fetch_add(1, Ordering::Relaxed);
            let path = std::env::temp_dir().join(format!(
                "amr_output_files_{}_{}_{}",
                std::process::id(),
                timestamp,
                counter
            ));
            fs::create_dir(&path).expect("test directory should be created");
            Self(path)
        }

        fn path(&self) -> &Path {
            &self.0
        }
    }

    impl Drop for TestDirectory {
        fn drop(&mut self) {
            if self.0.parent() == Some(std::env::temp_dir().as_path())
                && self
                    .0
                    .file_name()
                    .is_some_and(|name| name.to_string_lossy().starts_with("amr_output_files_"))
            {
                let _ = fs::remove_dir_all(&self.0);
            }
        }
    }

    #[test]
    fn first_publication_uses_requested_name_and_preserves_bytes() {
        let directory = TestDirectory::new();
        let staged = directory.path().join("summary.incomplete");
        let requested = directory.path().join("simulation_summary_123456.csv");
        let payload = b"time_step,run_id\n0,123456\n1,123456\n";
        fs::write(&staged, payload).unwrap();

        let published = publish_summary_no_clobber(&staged, &requested).unwrap();

        assert_eq!(published, requested);
        assert_eq!(fs::read(&published).unwrap(), payload);
        assert!(!staged.exists());
        assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 1);
    }

    #[test]
    fn collisions_preserve_every_existing_payload_and_choose_next_repeat() {
        let directory = TestDirectory::new();
        let requested = directory.path().join("simulation_summary_123456.csv");
        let repeat_one = directory
            .path()
            .join("simulation_summary_123456_repeat_1.csv");
        let staged = directory.path().join("new_summary.incomplete");
        fs::write(&requested, b"original run").unwrap();
        fs::write(&repeat_one, b"earlier replay").unwrap();
        fs::write(&staged, b"new replay").unwrap();

        let published = publish_summary_no_clobber(&staged, &requested).unwrap();

        assert_eq!(
            published,
            directory
                .path()
                .join("simulation_summary_123456_repeat_2.csv")
        );
        assert_eq!(fs::read(&requested).unwrap(), b"original run");
        assert_eq!(fs::read(&repeat_one).unwrap(), b"earlier replay");
        assert_eq!(fs::read(&published).unwrap(), b"new replay");
        assert!(!staged.exists());
        assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 3);
    }

    #[test]
    fn identical_replays_still_get_distinct_artifacts_without_changing_content() {
        let directory = TestDirectory::new();
        let requested = directory.path().join("simulation_summary_654321.csv");
        let payload = b"time_step,run_id\n0,654321\n";
        let mut published = Vec::new();
        for repeat in 0..3 {
            let staged = directory.path().join(format!("replay_{repeat}.incomplete"));
            fs::write(&staged, payload).unwrap();
            published.push(publish_summary_no_clobber(&staged, &requested).unwrap());
            assert!(!staged.exists());
        }

        assert_eq!(published.iter().collect::<HashSet<_>>().len(), 3);
        assert_eq!(published[0], requested);
        for path in published {
            assert_eq!(fs::read(path).unwrap(), payload);
        }
        assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 3);
    }

    #[test]
    fn concurrent_collisions_publish_each_payload_once_without_replacing_files() {
        const WRITERS: usize = 8;
        let directory = TestDirectory::new();
        let requested = directory.path().join("simulation_summary_123456.csv");
        fs::write(&requested, b"original run").unwrap();
        let barrier = Arc::new(Barrier::new(WRITERS));
        let mut workers = Vec::new();
        for writer in 0..WRITERS {
            let staged = directory.path().join(format!("writer_{writer}.incomplete"));
            let payload = format!("complete payload from writer {writer}\n");
            fs::write(&staged, payload.as_bytes()).unwrap();
            let requested = requested.clone();
            let barrier = Arc::clone(&barrier);
            workers.push(std::thread::spawn(move || {
                barrier.wait();
                let published = publish_summary_no_clobber(&staged, &requested).unwrap();
                assert!(!staged.exists());
                (published, payload)
            }));
        }

        let mut published_paths = HashSet::new();
        for worker in workers {
            let (published, payload) = worker.join().expect("publisher should not panic");
            assert_eq!(fs::read(&published).unwrap(), payload.as_bytes());
            assert!(published_paths.insert(published));
        }
        let expected_paths = (1..=WRITERS)
            .map(|repeat| {
                directory
                    .path()
                    .join(format!("simulation_summary_123456_repeat_{repeat}.csv"))
            })
            .collect::<HashSet<_>>();
        assert_eq!(published_paths, expected_paths);
        assert_eq!(fs::read(&requested).unwrap(), b"original run");
        assert_eq!(fs::read_dir(directory.path()).unwrap().count(), WRITERS + 1);
    }

    #[test]
    fn existing_directory_is_an_error_and_retains_staged_data() {
        let directory = TestDirectory::new();
        let staged = directory.path().join("summary.incomplete");
        let requested = directory.path().join("simulation_summary_123456.csv");
        fs::write(&staged, b"complete staged data").unwrap();
        fs::create_dir(&requested).unwrap();
        let existing = requested.join("keep.txt");
        fs::write(&existing, b"existing directory contents").unwrap();

        assert!(publish_summary_no_clobber(&staged, &requested).is_err());

        assert_eq!(fs::read(&staged).unwrap(), b"complete staged data");
        assert_eq!(fs::read(&existing).unwrap(), b"existing directory contents");
        assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 2);
    }

    #[test]
    fn directory_at_repeat_name_does_not_get_skipped_or_modified() {
        let directory = TestDirectory::new();
        let staged = directory.path().join("summary.incomplete");
        let requested = directory.path().join("summary.csv");
        let collision = directory.path().join("summary_repeat_1.csv");
        fs::write(&staged, b"complete staged data").unwrap();
        fs::write(&requested, b"original run").unwrap();
        fs::create_dir(&collision).unwrap();

        assert!(publish_summary_no_clobber(&staged, &requested).is_err());

        assert_eq!(fs::read(&staged).unwrap(), b"complete staged data");
        assert_eq!(fs::read(&requested).unwrap(), b"original run");
        assert!(collision.is_dir());
        assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 3);
    }

    #[test]
    fn missing_staged_file_returns_original_io_kind_without_creating_output() {
        let directory = TestDirectory::new();
        let staged = directory.path().join("missing.incomplete");
        let requested = directory.path().join("summary.csv");

        let error = publish_summary_no_clobber(&staged, &requested).unwrap_err();

        assert_eq!(error.kind(), io::ErrorKind::NotFound);
        assert_eq!(fs::read_dir(directory.path()).unwrap().count(), 0);
    }

    #[cfg(unix)]
    #[test]
    fn existing_symlinks_are_errors_and_remain_unchanged() {
        use std::os::unix::fs::symlink;

        let directory = TestDirectory::new();
        let staged = directory.path().join("summary.incomplete");
        fs::write(&staged, b"complete staged data").unwrap();
        let target = directory.path().join("original.csv");
        fs::write(&target, b"original run").unwrap();

        for (name, link_target) in [
            ("linked.csv", target.clone()),
            ("dangling.csv", directory.path().join("absent.csv")),
        ] {
            let requested = directory.path().join(name);
            symlink(&link_target, &requested).unwrap();
            assert!(publish_summary_no_clobber(&staged, &requested).is_err());
            assert_eq!(fs::read_link(&requested).unwrap(), link_target);
            assert_eq!(fs::read(&staged).unwrap(), b"complete staged data");
        }
        assert_eq!(fs::read(&target).unwrap(), b"original run");
    }
}
