"""Identify completed output artifacts separately from deterministic numeric run IDs.

The first output keeps its historical numeric name; later publications retain a
positive ``_repeat_N`` suffix. That suffix belongs in derived output names even
though CSV rows continue to contain the original numeric ``run_id``.
"""

from __future__ import annotations

import os
from pathlib import PurePath
import re


_ARTIFACT_TOKEN = r"[0-9]{6,7}(?:_repeat_[1-9][0-9]*)?"
_BARE_ARTIFACT = re.compile(rf"^{_ARTIFACT_TOKEN}$")
_SIMULATION_NAME = re.compile(rf"^simulation_summary_({_ARTIFACT_TOKEN})$")
_CALIBRATION_NAME = re.compile(
    rf"^calibration_summary_[A-Za-z_\-]*({_ARTIFACT_TOKEN})$"
)
_LEGACY_IDENTIFIER = re.compile(r"(?<![0-9])([0-9]{6,7})(?![0-9])")


def extract_run_artifact_id(path: str | os.PathLike[str] | None) -> str | None:
    """Return the complete numeric/repeat token from a summary name or bare token.

    Canonical simulation CSVs and calibration TXT reports are supported, along
    with historical names such as ``calibration_summary_fff123456.txt``. Numeric
    IDs have six or seven digits, including the model's maximum ID ``1000000``.
    Nonrepeat legacy filenames retain the historical embedded-ID fallback.

    The repeat marker is reserved: an invalid repeat name never falls back to a
    base ID. Cache, temporary, and incomplete-file extensions are not artifacts.
    """
    if path is None:
        return None
    # Handle recorded Windows source paths when reports are analyzed on Unix too.
    name = re.split(r"[/\\]", os.fspath(path))[-1]
    if not name:
        return None
    suffix = PurePath(name).suffix.lower()
    if suffix not in {"", ".csv", ".txt"}:
        return None
    stem = PurePath(name).stem if suffix else name

    if not suffix and _BARE_ARTIFACT.fullmatch(stem):
        return stem
    if suffix in {"", ".csv"}:
        match = _SIMULATION_NAME.fullmatch(stem)
        if match:
            return match.group(1)
    if suffix in {"", ".txt"}:
        match = _CALIBRATION_NAME.fullmatch(stem)
        if match:
            return match.group(1)

    if "repeat" in stem.casefold():
        return None
    match = _LEGACY_IDENTIFIER.search(stem)
    return match.group(1) if match else None
