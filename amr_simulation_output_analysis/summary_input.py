"""Resolve AMR summary inputs without treating model-local IDs as run identity."""

from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Mapping, Sequence


_MODEL_RUN_ID_RE = re.compile(r"^simulation_summary_(\d{6})(?:--[^.]*)?\.csv$", re.IGNORECASE)
_EXPORT_ID_RE = re.compile(r"--([A-Za-z0-9][A-Za-z0-9_-]*)\.csv$", re.IGNORECASE)


class SummaryInputError(ValueError):
    """Raised when a summary input is missing, unsafe, or ambiguous."""


def model_run_id_from_filename(path: str | Path) -> str | None:
    """Return the model-local six-digit ID as evidence, never as global identity."""

    match = _MODEL_RUN_ID_RE.match(Path(path).name)
    return match.group(1) if match else None


def canonical_summary_identity(path: str | Path) -> str:
    """Return a collision-safe local identity for one exported summary."""

    resolved = Path(path).resolve()
    match = _EXPORT_ID_RE.search(resolved.name)
    if match:
        return f"submission:{match.group(1)}"
    return f"path:{resolved.as_posix()}"


def _is_summary_filename(path: Path) -> bool:
    name = path.name.lower()
    return (
        name == "summary.csv"
        or name.startswith("summary--")
        or name.startswith("simulation_summary_")
    ) and name.endswith(".csv")


def discover_summary_csvs(directory: str | Path) -> tuple[Path, ...]:
    """Return all recognizable summary CSVs below a directory in stable order."""

    root = Path(directory)
    if not root.is_dir():
        raise SummaryInputError(f"summary input directory does not exist: {root}")
    return tuple(
        sorted(
            {path.resolve() for path in root.rglob("*.csv") if path.is_file() and _is_summary_filename(path)},
            key=lambda path: path.as_posix(),
        )
    )


def _contained_existing_file(candidate: Path, root: Path) -> Path | None:
    try:
        resolved = candidate.resolve()
        resolved.relative_to(root.resolve())
    except (OSError, ValueError):
        return None
    return resolved if resolved.is_file() else None


def _manifest_summary_candidates(manifest_path: Path, manifest: Mapping[str, object]) -> tuple[Path, ...]:
    root = manifest_path.parent.resolve()
    for field in ("summary_csv", "summary_upload_csv"):
        value = str(manifest.get(field) or "").strip()
        if not value:
            continue
        raw = Path(value)
        if raw.is_absolute():
            accepted = _contained_existing_file(raw, root)
            if accepted is not None:
                return (accepted,)
            continue
        accepted = _contained_existing_file(root / raw, root)
        if accepted is None and ".." in raw.parts:
            raise SummaryInputError(f"run manifest {field} escapes its export directory: {value}")
        if accepted is not None:
            return (accepted,)

    original_name = str(manifest.get("summary_original_filename") or "").strip()
    if original_name:
        original = Path(original_name)
        if original.name != original_name or original_name in {".", ".."}:
            raise SummaryInputError("run manifest summary_original_filename must be a basename")
        candidates = tuple(
            path
            for path in discover_summary_csvs(root)
            if path.name == original_name
            or (
                path.suffix.lower() == original.suffix.lower()
                and path.stem.startswith(f"{original.stem}--")
            )
        )
        return candidates
    return ()


def _read_manifest(path: Path) -> Mapping[str, object]:
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise SummaryInputError(f"run manifest is not readable JSON: {path}") from exc
    if not isinstance(payload, Mapping):
        raise SummaryInputError(f"run manifest must contain a JSON object: {path}")
    return payload


def _require_one(candidates: Sequence[Path], *, source: Path) -> Path:
    unique = tuple(dict.fromkeys(path.resolve() for path in candidates))
    if len(unique) == 1:
        return unique[0]
    if not unique:
        raise SummaryInputError(f"no AMR summary CSV was found for {source}")
    names = ", ".join(path.name for path in unique[:5])
    raise SummaryInputError(
        f"multiple AMR summary CSVs were found for {source}: {names}; "
        "select one explicitly or use a multi-run tool"
    )


def resolve_summary_csv(source: str | Path) -> Path:
    """Resolve one explicit CSV, run manifest, or unambiguous output directory."""

    input_path = Path(source)
    if input_path.is_file() and input_path.suffix.lower() == ".csv":
        if not _is_summary_filename(input_path):
            raise SummaryInputError(f"file is not a recognized AMR summary CSV: {input_path}")
        return input_path.resolve()
    if input_path.is_file() and input_path.suffix.lower() == ".json":
        manifest = _read_manifest(input_path)
        candidates = _manifest_summary_candidates(input_path, manifest)
        if not candidates:
            candidates = discover_summary_csvs(input_path.parent)
        return _require_one(candidates, source=input_path)
    if input_path.is_dir():
        return _require_one(discover_summary_csvs(input_path), source=input_path)
    raise SummaryInputError(f"summary input does not exist: {input_path}")
