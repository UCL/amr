"""Regenerate and publish a consistent resistance-target set from Rust potencies.

Run from the repository root with::

    python -m amr_simulation_output_analysis.refresh_resistance_targets

Compilation and generation happen in a staging directory. Cooperative readers
wait only during the short publication step, when the update marker is present.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import shutil
import subprocess
import tempfile
import uuid
from pathlib import Path
from typing import Dict, Optional, Tuple

from .build_resistance_targets_v1 import (
    MANIFEST_FILENAME,
    TARGET_SET_VERSION,
    UPDATE_LOCK_FILENAME,
    build_resistance_targets_v1,
)


RESISTANCE_TARGET_UPDATE_LOCK_FILENAME = UPDATE_LOCK_FILENAME
INPUT_FILENAMES = (
    "resistance_prevalence_values.csv",
    "resistance_average_resistant_values.csv",
    "resistance_targets_v1.schema.json",
)
# The manifest is the commit record and must be published last.
GENERATED_FILENAMES = (
    "model_potency_matrix.csv",
    "model_resistance_reachability_matrix.csv",
    "resistance_targets_v1.csv",
    "resistance_target_sources_v1.csv",
    MANIFEST_FILENAME,
)


def _rust_source_hashes(root: Path) -> Dict[str, str]:
    paths = list((root / "src").rglob("*.rs"))
    paths.extend(root / name for name in ("Cargo.toml", "Cargo.lock") if (root / name).is_file())
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(paths)
    }


def _verify_rust_sources(root: Path, expected: Dict[str, str]) -> None:
    if _rust_source_hashes(root) != expected:
        raise RuntimeError(
            "Rust source inputs changed during resistance-target refresh; rerun the refresh"
        )


def _export_rust_projections(root: Path, staged_data: Path) -> None:
    for binary, filename in (
        ("export_potency_matrix", "model_potency_matrix.csv"),
        (
            "export_resistance_reachability_matrix",
            "model_resistance_reachability_matrix.csv",
        ),
    ):
        subprocess.run(
            ["cargo", "run", "--quiet", "--bin", binary, "--", str(staged_data / filename)],
            cwd=root,
            check=True,
        )


def _verify_staged_manifest(data_dir: Path) -> None:
    manifest = json.loads((data_dir / MANIFEST_FILENAME).read_text(encoding="utf-8"))
    expected = set(INPUT_FILENAMES + GENERATED_FILENAMES[:-1])
    if (
        manifest.get("target_set_version") != TARGET_SET_VERSION
        or manifest.get("hash_algorithm") != "sha256"
        or set(manifest.get("artifacts", {})) != expected
    ):
        raise ValueError("Staged resistance-target manifest has an invalid artifact contract")
    for filename, entry in manifest["artifacts"].items():
        content = (data_dir / filename).read_bytes()
        if (
            len(content) != entry.get("bytes")
            or hashlib.sha256(content).hexdigest() != entry.get("sha256")
        ):
            raise ValueError(f"Staged resistance-target manifest mismatch: {filename}")


def _publish_staged(
    staged_data: Path, data_dir: Path, root: Path, rust_sources: Dict[str, str]
) -> None:
    marker = data_dir / RESISTANCE_TARGET_UPDATE_LOCK_FILENAME
    token = f"pid={os.getpid()} token={uuid.uuid4().hex}\n".encode("ascii")
    try:
        handle = marker.open("xb")
    except FileExistsError as error:
        raise RuntimeError(
            f"Resistance-target refresh already in progress: {marker}. "
            "If no refresh is running, investigate the retained marker before retrying."
        ) from error

    release_marker = True
    replaced = []
    backup_dir = staged_data.parent / "backup"
    originals = {}
    try:
        with handle:
            handle.write(token)
        _verify_rust_sources(root, rust_sources)
        # A researcher may edit these while Rust is compiling. Do not publish a
        # manifest for old copies of the editable sources.
        for filename in INPUT_FILENAMES:
            if (data_dir / filename).read_bytes() != (staged_data / filename).read_bytes():
                raise RuntimeError(
                    f"Resistance-target input changed during refresh: {filename}; rerun the refresh"
                )

        backup_dir.mkdir()
        for filename in GENERATED_FILENAMES:
            original = data_dir / filename
            originals[filename] = original.exists()
            if originals[filename]:
                shutil.copyfile(original, backup_dir / filename)

        try:
            for filename in GENERATED_FILENAMES:
                os.replace(staged_data / filename, data_dir / filename)
                replaced.append(filename)
        except BaseException as publication_error:
            rollback_errors = []
            for filename in reversed(replaced):
                try:
                    if originals[filename]:
                        os.replace(backup_dir / filename, data_dir / filename)
                    else:
                        (data_dir / filename).unlink()
                except BaseException as error:
                    rollback_errors.append(f"{filename}: {error}")
            if rollback_errors:
                release_marker = False
                raise RuntimeError(
                    "Resistance-target publication and rollback failed. "
                    f"The update marker {marker} and recovery files in {backup_dir} "
                    "were retained; investigate and restore the target set before "
                    "removing the marker. Rollback errors: " + "; ".join(rollback_errors)
                ) from publication_error
            raise
    finally:
        # Exclusive creation gives this writer ownership. The token also stops
        # cleanup from deleting a marker that somebody has replaced meanwhile.
        if release_marker and marker.exists() and marker.read_bytes() == token:
            marker.unlink()


def refresh_resistance_targets(root: Optional[Path] = None) -> Tuple[Path, ...]:
    """Export Rust matrices, rebuild targets, and publish the five generated files.

    Existing live artifacts are untouched if staging fails. Publication failures
    restore the previous artifacts before allowing readers to resume. A failed
    restoration retains the marker and staging directory for manual recovery.
    """
    project_root = Path(root).resolve() if root is not None else Path(__file__).resolve().parents[1]
    data_dir = (project_root / "data").resolve()
    # Same filesystem as the live files, so os.replace remains atomic per file.
    staged_root = Path(tempfile.mkdtemp(prefix=".resistance-targets-", dir=data_dir))
    staged_data = staged_root / "data"
    staged_data.mkdir()
    try:
        for filename in INPUT_FILENAMES:
            shutil.copyfile(data_dir / filename, staged_data / filename)
        rust_sources = _rust_source_hashes(project_root)
        _export_rust_projections(project_root, staged_data)
        _verify_rust_sources(project_root, rust_sources)
        build_resistance_targets_v1(staged_root)
        _verify_staged_manifest(staged_data)
        _publish_staged(staged_data, data_dir, project_root, rust_sources)
    finally:
        # A backup directory left behind with the marker may be needed to repair
        # an unsuccessful rollback. Other failures have no published changes.
        marker = data_dir / RESISTANCE_TARGET_UPDATE_LOCK_FILENAME
        if not ((staged_root / "backup").exists() and marker.exists()):
            resolved_stage = staged_root.resolve()
            if resolved_stage.parent != data_dir or not resolved_stage.name.startswith(".resistance-targets-"):
                raise RuntimeError(f"Refusing to remove unexpected staging directory: {resolved_stage}")
            shutil.rmtree(resolved_stage)
    return tuple(data_dir / filename for filename in GENERATED_FILENAMES)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, help="Repository root (defaults to this module's repository)")
    args = parser.parse_args()
    for path in refresh_resistance_targets(args.root):
        print(f"Wrote {path}")


if __name__ == "__main__":
    main()
