"""Provenance contract for optional model-comparison overlays.

The historical overlay files contain useful source-patterned values, but their
row-level observational provenance is not recoverable.  This module keeps
those values available as best-guess placeholders while ensuring that an
authoritative-looking source label can never make a row observational evidence.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Optional

import pandas as pd


OBSERVED_COMPARISON = "observed_comparison"
GENERATED_BEST_GUESS_PLACEHOLDER = "generated_best_guess_placeholder"
SOURCE_INFORMED_BEST_GUESS_PLACEHOLDER = (
    "source_informed_best_guess_placeholder_provenance_unverified"
)

ALLOWED_PROVENANCE_CLASSES = frozenset(
    {
        OBSERVED_COMPARISON,
        GENERATED_BEST_GUESS_PLACEHOLDER,
        SOURCE_INFORMED_BEST_GUESS_PLACEHOLDER,
    }
)

PROVENANCE_COLUMNS = (
    "overlay_provenance_class",
    "generated",
    "generation_method",
    "source_id",
    "source_url_or_doi",
    "reference_year",
    "uncertainty",
    "rationale",
    "last_reviewed",
)

REGISTRY_COLUMNS = (
    "relative_path",
    "source_quality",
    *PROVENANCE_COLUMNS,
)

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DEFAULT_REGISTRY_PATH = PROJECT_ROOT / "data" / "empirical" / "overlay_provenance_v1.csv"


class OverlayProvenanceError(ValueError):
    """Raised when overlay provenance metadata violates the contract."""


def _text(series: pd.Series) -> pd.Series:
    return series.astype("object").fillna("").astype(str).str.strip()


def _nullable_bool(series: pd.Series, context: str) -> pd.Series:
    true_values = {"true", "1", "yes"}
    false_values = {"false", "0", "no"}
    values = _text(series).str.lower()
    invalid = ~values.isin(true_values | false_values | {""})
    if invalid.any():
        examples = sorted(values[invalid].unique())[:3]
        raise OverlayProvenanceError(
            f"{context}: generated must be true, false, or blank; found {examples}"
        )

    converted = pd.Series(pd.NA, index=series.index, dtype="boolean")
    converted.loc[values.isin(true_values)] = True
    converted.loc[values.isin(false_values)] = False
    return converted


def _validate_metadata(df: pd.DataFrame, context: str) -> pd.DataFrame:
    missing = [column for column in PROVENANCE_COLUMNS if column not in df.columns]
    if missing:
        raise OverlayProvenanceError(
            f"{context}: incomplete overlay provenance metadata; missing {missing}"
        )

    result = df.copy()
    result["overlay_provenance_class"] = _text(result["overlay_provenance_class"])
    unknown = sorted(
        set(result["overlay_provenance_class"]) - ALLOWED_PROVENANCE_CLASSES
    )
    if unknown:
        raise OverlayProvenanceError(
            f"{context}: unknown overlay provenance classes {unknown}"
        )

    result["generated"] = _nullable_bool(result["generated"], context)
    for column in PROVENANCE_COLUMNS:
        if column not in {"overlay_provenance_class", "generated"}:
            result[column] = _text(result[column])

    observed = result["overlay_provenance_class"].eq(OBSERVED_COMPARISON)
    observed_generated_ok = result["generated"].eq(False).fillna(False)
    if (observed & ~observed_generated_ok).any():
        raise OverlayProvenanceError(
            f"{context}: observed comparison rows must explicitly set generated=false"
        )
    for column in (
        "generation_method",
        "source_id",
        "source_url_or_doi",
        "reference_year",
        "uncertainty",
        "rationale",
        "last_reviewed",
    ):
        if (observed & result[column].eq("")).any():
            raise OverlayProvenanceError(
                f"{context}: observed comparison rows require non-blank {column}"
            )

    generated = result["overlay_provenance_class"].eq(
        GENERATED_BEST_GUESS_PLACEHOLDER
    )
    generated_flag_ok = result["generated"].eq(True).fillna(False)
    if (generated & ~generated_flag_ok).any():
        raise OverlayProvenanceError(
            f"{context}: generated placeholders must explicitly set generated=true"
        )
    for column in ("generation_method", "rationale", "last_reviewed"):
        if (generated & result[column].eq("")).any():
            raise OverlayProvenanceError(
                f"{context}: generated placeholders require non-blank {column}"
            )

    unverified = result["overlay_provenance_class"].eq(
        SOURCE_INFORMED_BEST_GUESS_PLACEHOLDER
    )
    if (unverified & result["rationale"].eq("")).any():
        raise OverlayProvenanceError(
            f"{context}: source-informed placeholders require a rationale"
        )
    if (unverified & result["last_reviewed"].eq("")).any():
        raise OverlayProvenanceError(
            f"{context}: source-informed placeholders require last_reviewed"
        )

    result["eligible_as_observed_comparison"] = observed
    return result


@lru_cache(maxsize=4)
def load_overlay_provenance_registry(
    registry_path: Path = DEFAULT_REGISTRY_PATH,
) -> pd.DataFrame:
    """Load and validate the file/source-quality provenance registry."""
    path = Path(registry_path)
    if not path.exists():
        raise OverlayProvenanceError(f"Overlay provenance registry not found: {path}")

    registry = pd.read_csv(path, dtype=str, keep_default_na=False)
    missing = [column for column in REGISTRY_COLUMNS if column not in registry.columns]
    if missing:
        raise OverlayProvenanceError(
            f"{path}: overlay provenance registry is missing columns {missing}"
        )

    registry = registry.loc[:, REGISTRY_COLUMNS].copy()
    registry["relative_path"] = _text(registry["relative_path"]).str.replace(
        "\\", "/", regex=False
    )
    registry["source_quality"] = _text(registry["source_quality"])
    duplicate = registry.duplicated(["relative_path", "source_quality"], keep=False)
    if duplicate.any():
        keys = registry.loc[duplicate, ["relative_path", "source_quality"]]
        raise OverlayProvenanceError(
            f"{path}: duplicate registry selectors {keys.to_dict('records')[:3]}"
        )

    _validate_metadata(registry, str(path))
    return registry


def _relative_source_path(source_path: Optional[Path]) -> Optional[str]:
    if source_path is None:
        return None
    path = Path(source_path).resolve()
    try:
        return path.relative_to(PROJECT_ROOT).as_posix()
    except ValueError:
        return path.as_posix()


def _fallback_metadata(df: pd.DataFrame, source_path: Optional[Path]) -> pd.DataFrame:
    source_quality = (
        _text(df["source_quality"])
        if "source_quality" in df.columns
        else pd.Series("", index=df.index, dtype="object")
    )
    path_label = _relative_source_path(source_path) or "in_memory_dataframe"
    result = df.copy()
    result["overlay_provenance_class"] = SOURCE_INFORMED_BEST_GUESS_PLACEHOLDER
    result["generated"] = pd.Series("", index=result.index, dtype="object")
    result["generation_method"] = "legacy_overlay_without_row_provenance"
    result["source_id"] = source_quality.where(
        source_quality.ne(""), f"unverified:{path_label}"
    )
    result["source_url_or_doi"] = ""
    result["reference_year"] = ""
    result["uncertainty"] = "not_quantified"
    result["rationale"] = (
        "Retained as a best-guess comparison placeholder; row-level observational "
        "provenance has not been recovered."
    )
    result["last_reviewed"] = "2026-07-22"
    return result


def annotate_overlay_provenance(
    df: pd.DataFrame,
    source_path: Optional[Path] = None,
    registry_path: Path = DEFAULT_REGISTRY_PATH,
) -> pd.DataFrame:
    """Attach and validate provenance without changing any numerical values."""
    present = [column for column in PROVENANCE_COLUMNS if column in df.columns]
    if present:
        if len(present) != len(PROVENANCE_COLUMNS):
            missing = [column for column in PROVENANCE_COLUMNS if column not in df.columns]
            raise OverlayProvenanceError(
                f"{source_path or 'dataframe'}: partial provenance metadata; missing {missing}"
            )
        return _validate_metadata(df, str(source_path or "dataframe"))

    result = _fallback_metadata(df, source_path)
    relative_path = _relative_source_path(source_path)
    if relative_path is None:
        return _validate_metadata(result, "in-memory overlay dataframe")

    registry = load_overlay_provenance_registry(Path(registry_path))
    matches = registry.loc[registry["relative_path"].eq(relative_path)]
    if matches.empty:
        return _validate_metadata(result, relative_path)

    source_quality = (
        _text(result["source_quality"])
        if "source_quality" in result.columns
        else pd.Series("", index=result.index, dtype="object")
    )

    # Apply the file-wide rule first, then source-quality-specific overrides.
    ordered = pd.concat(
        [
            matches.loc[matches["source_quality"].eq("*")],
            matches.loc[~matches["source_quality"].eq("*")],
        ],
        ignore_index=True,
    )
    for _, metadata in ordered.iterrows():
        selector = metadata["source_quality"]
        mask = pd.Series(True, index=result.index) if selector == "*" else source_quality.eq(selector)
        for column in PROVENANCE_COLUMNS:
            value = metadata[column]
            if column == "source_id" and value == "{source_quality}":
                result.loc[mask, column] = source_quality.loc[mask]
            else:
                result.loc[mask, column] = value

    return _validate_metadata(result, relative_path)


def filter_overlay_rows(
    df: Optional[pd.DataFrame],
    *,
    include_best_guess_placeholders: bool = False,
    source_path: Optional[Path] = None,
) -> Optional[pd.DataFrame]:
    """Return observed rows by default, or all rows after explicit placeholder opt-in."""
    if df is None:
        return None
    annotated = annotate_overlay_provenance(df, source_path=source_path)
    if include_best_guess_placeholders:
        return annotated
    return annotated.loc[annotated["eligible_as_observed_comparison"]].copy()


def overlay_display_label(df: Optional[pd.DataFrame]) -> str:
    """Return an honest plot label for the provenance classes present."""
    if df is None or df.empty or "overlay_provenance_class" not in df.columns:
        return "Comparison overlay"
    classes = set(_text(df["overlay_provenance_class"]))
    if classes == {OBSERVED_COMPARISON}:
        return "Observed comparison"
    if classes == {GENERATED_BEST_GUESS_PLACEHOLDER}:
        return "Generated best-guess placeholder"
    if classes == {SOURCE_INFORMED_BEST_GUESS_PLACEHOLDER}:
        return "Source-informed best-guess placeholder"
    return "Comparison overlay (includes best-guess placeholders)"
