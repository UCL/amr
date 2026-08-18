"""Build the versioned display-only plausible-range registry for calibration targets.

The central values remain authoritative in their existing JSON/CSV files. This
builder copies those values into the registry so tests can detect drift, then
adds either an explicit range or a transparent tier-based expert range.

These ranges are not calibration-score tolerances and do not enter the score.
Only rows explicitly labelled ``published_uncertainty_range`` reproduce a
published interval. All other bounds are plausible ranges for communication
and sensitivity review.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = REPO_ROOT / "data"
OUTPUT_PATH = DATA_ROOT / "calibration_target_ranges_v1.csv"
TARGET_YEAR = 2025
TARGET_SET_VERSION = "calibration-target-ranges-v1"
LAST_REVIEWED = "2026-07-28"

FIELDNAMES = [
    "target_set_version",
    "target_family",
    "target_key",
    "target_label",
    "target_year",
    "central_value",
    "plausible_lower",
    "plausible_upper",
    "unit",
    "interval_kind",
    "provenance_class",
    "source_id",
    "source_url_or_doi",
    "range_method",
    "rationale",
    "last_reviewed",
]

GBD_BACTERIAL_MORTALITY_URL = "https://pubmed.ncbi.nlm.nih.gov/36423648/"
GBD_2021_SEPSIS_URL = "https://pubmed.ncbi.nlm.nih.gov/41135560/"
WHO_ANTIBIOTIC_USE_URL = (
    "https://www.who.int/data/gho/data/indicators/indicator-details/GHO/"
    "antibacterial-consumption--total-consumption-of-antibacterials-expressed-as-"
    "ddd-per-1000-inhabitants-per-day"
)
GLOBAL_ANTIBIOTIC_USE_URL = "https://pmc.ncbi.nlm.nih.gov/articles/PMC8654683/"
WHO_GLASS_2022_URL = "https://www.who.int/publications/i/item/9789240108127"
WHO_CHOLERA_URL = "https://www.who.int/news-room/fact-sheets/detail/cholera"
WHO_FERG_2021_DIARRHOEAL_URL = "https://pubmed.ncbi.nlm.nih.gov/42296981/"


def _round_significant(value: float, digits: int = 2) -> float:
    if value == 0.0:
        return 0.0
    places = digits - int(math.floor(math.log10(abs(value)))) - 1
    return round(value, places)


def _make_row(
    *,
    family: str,
    key: str,
    label: str,
    central: float,
    lower: float,
    upper: float,
    unit: str,
    interval_kind: str,
    provenance_class: str,
    source_id: str,
    source_url: str,
    range_method: str,
    rationale: str,
    last_reviewed: str = LAST_REVIEWED,
) -> dict[str, Any]:
    if not math.isfinite(central) or not math.isfinite(lower) or not math.isfinite(upper):
        raise ValueError(f"Non-finite target range for {family}/{key}")
    if lower > central or central > upper:
        raise ValueError(
            f"Range does not contain central value for {family}/{key}: "
            f"{lower} <= {central} <= {upper}"
        )
    if lower < 0:
        raise ValueError(f"Negative lower bound for {family}/{key}")
    if unit == "proportion" and upper > 1:
        raise ValueError(f"Proportion upper bound exceeds one for {family}/{key}")

    return {
        "target_set_version": TARGET_SET_VERSION,
        "target_family": family,
        "target_key": key,
        "target_label": label,
        "target_year": TARGET_YEAR,
        "central_value": central,
        "plausible_lower": lower,
        "plausible_upper": upper,
        "unit": unit,
        "interval_kind": interval_kind,
        "provenance_class": provenance_class,
        "source_id": source_id,
        "source_url_or_doi": source_url,
        "range_method": range_method,
        "rationale": rationale,
        "last_reviewed": last_reviewed,
    }


HEADLINE_RANGES = {
    "infection_deaths_millions": {
        "lower": 5.7,
        "upper": 10.2,
        "interval_kind": "published_uncertainty_range",
        "provenance_class": "published_source_range",
        "source_id": "gbd_33_bacterial_pathogens_2019",
        "source_url": GBD_BACTERIAL_MORTALITY_URL,
        "range_method": "Published GBD 2019 95% uncertainty interval.",
        "rationale": (
            "The published 7.7 million estimate is used unchanged as a pragmatic "
            "calibration benchmark. The model's organism set and target year differ "
            "slightly; those differences are documented rather than encoded as a "
            "precision-implying numerical adjustment."
        ),
    },
    "people_on_antibiotics_millions": {
        "lower": 70.0,
        "upper": 150.0,
        "interval_kind": "expert_plausible_range",
        "provenance_class": "source_informed_transformed",
        "source_id": "who_glass_amu_and_global_spatial_model",
        "source_url": WHO_ANTIBIOTIC_USE_URL,
        "range_method": "Expert conversion range from DDD/day to unique current users.",
        "rationale": (
            "WHO describes DDD/1000/day as only a rough proxy for people treated. The "
            "range spans uncertainty in dose intensity, combination treatment, data "
            "coverage, stock-versus-use measurement, and conversion to unique users."
        ),
    },
    "annual_infection_incidence_percent": {
        "lower": 15.0,
        "upper": 30.0,
        "interval_kind": "derived_plausible_range",
        "provenance_class": "source_informed_transformed",
        "source_id": "model_scope_pathogen_incidence_synthesis_who_ferg_2021",
        "source_url": WHO_FERG_2021_DIARRHOEAL_URL,
        "last_reviewed": "2026-08-18",
        "range_method": (
            "Derived model-scope range around selectively updated organism targets, "
            "allowing for same-day polymicrobial acquisition and cross-source uncertainty."
        ),
        "rationale": (
            "The 42 organism targets sum to 22.34% after updating E. coli, Shigella, and "
            "Campylobacter from WHO FERG 2021 estimates. The 20% headline counts each "
            "person at most once per day and therefore allows for same-day polymicrobial "
            "acquisition. It is a derived benchmark, not a published all-bacteria estimate."
        ),
    },
    "sepsis_incident_cases_millions": {
        "lower": 50.0,
        "upper": 100.0,
        "interval_kind": "derived_plausible_range",
        "provenance_class": "source_informed_transformed",
        "source_id": "gbd_2021_global_sepsis_2025_bacterial_subset",
        "source_url": GBD_2021_SEPSIS_URL,
        "last_reviewed": "2026-08-18",
        "range_method": (
            "Use a broad 50-100 million bacterial-subset range within the published "
            "135-201 million all-cause sepsis UI for 2021."
        ),
        "rationale": (
            "GBD estimates 166 million all-cause sepsis cases in 2021. The model is "
            "bacteria-only and excludes viral, fungal, and parasitic sepsis while not "
            "representing the full burden of sepsis complicating non-infectious causes. "
            "The 70 million central target and 50-100 million range are explicit "
            "model-scope assumptions, not a published GBD point estimate or interval."
        ),
    },
}


DRUG_CLASS_RANGES = {
    "Penicillins (J01C)": (12.0, 24.0),
    "Beta-lactamase combinations (J01CR)": (14.0, 25.0),
    "Cephalosporins 1-2G": (6.0, 14.0),
    "Cephalosporins 3G": (3.0, 8.0),
    "Cephalosporins 3G/BLI": (0.1, 1.2),
    "Cephalosporins 4G": (1.0, 4.0),
    "Anti-MRSA Cephalosporins (5G)": (0.03, 0.6),
    "Siderophore Cephalosporins": (0.01, 0.4),
    "Novel BL/BLI": (0.03, 0.8),
    "Monobactams": (0.1, 1.5),
    "Macrolides (J01F)": (7.0, 15.0),
    "Lincosamides (J01FF)": (1.0, 4.0),
    "Quinolones/Fluoroquinolones (J01M)": (6.0, 15.0),
    "Tetracyclines (J01A)": (3.0, 10.0),
    "Carbapenems (J01DH)": (1.0, 4.0),
    "Aminoglycosides (J01G)": (1.0, 4.0),
    "Sulfonamides (J01E)": (2.0, 8.0),
    "Nitrofurans (J01XE)": (1.5, 7.0),
    "Fosfomycin (J01XX01)": (0.3, 2.5),
    "Glycopeptides (J01XA)": (0.8, 4.0),
    "Lipoglycopeptides & Oxazolidinones": (0.3, 2.5),
    "Lipopeptides (J01XX09)": (0.1, 1.5),
    "Fidaxomicin": (0.01, 0.4),
    "Nitroimidazoles": (0.3, 3.0),
    "Polymyxins (J01XB)": (0.005, 0.2),
    "Rifamycins (J04AB)": (0.03, 0.6),
    "Chloramphenicol (J01BA)": (0.01, 0.5),
    "Other Antibiotics": (0.01, 0.6),
}


CARRIAGE_RANGE_OVERRIDES = {
    "acinetobacter baumannii": (0.003, 0.03),
    "enterococcus faecalis": (0.6, 0.95),
    "enterococcus faecium": (0.1, 0.5),
    "escherichia coli": (0.85, 0.99),
    "klebsiella pneumoniae": (0.05, 0.4),
    "pseudomonas aeruginosa": (0.01, 0.1),
    "staphylococcus aureus": (0.2, 0.4),
    "staphylococcus epidermidis": (0.85, 0.99),
    "streptococcus pneumoniae": (0.1, 0.6),
    "salmonella enterica serovar typhi": (0.00005, 0.0006),
    "neisseria gonorrhoeae": (0.002, 0.01),
    "streptococcus pyogenes": (0.05, 0.2),
    "streptococcus agalactiae": (0.15, 0.35),
    "haemophilus influenzae": (0.25, 0.7),
    "chlamydia trachomatis": (0.008, 0.03),
    "neisseria meningitidis": (0.05, 0.2),
    "listeria monocytogenes": (0.01, 0.1),
    "clostridioides difficile": (0.02, 0.2),
    "moraxella catarrhalis": (0.2, 0.6),
    "mdr mycobacterium tuberculosis": (0.001, 0.005),
    "bacteroides fragilis": (0.7, 0.95),
    "mycoplasma genitalium": (0.005, 0.02),
    "mycoplasma pneumoniae": (0.01, 0.05),
}


INCIDENCE_RANGE_OVERRIDES = {
    "escherichia coli": (0.0185, 0.0516),
    "salmonella enterica serovar typhi": (0.0008, 0.0015),
    "shigella spp.": (0.0285, 0.0868),
    "neisseria gonorrhoeae": (0.007, 0.013),
    "chlamydia trachomatis": (0.012, 0.02),
    "vibrio cholerae": (0.00015, 0.0006),
    "treponema pallidum": (0.0007, 0.0013),
    "bordetella pertussis": (0.001, 0.003),
    "mdr mycobacterium tuberculosis": (0.00004, 0.00007),
    "campylobacter jejuni": (0.0206, 0.0560),
}


WHO_FERG_2021_INCIDENCE_KEYS = {
    "escherichia coli",
    "shigella spp.",
    "campylobacter jejuni",
}


DEATH_RANGE_OVERRIDES = {
    "salmonella enterica serovar typhi": (0.08, 0.15),
    "vibrio cholerae": (0.021, 0.143),
    "bordetella pertussis": (0.08, 0.25),
    "mdr mycobacterium tuberculosis": (0.15, 0.25),
}


PLACEHOLDER_MARKERS = (
    "placeholder",
    "uncertain",
    "extrapolat",
    "rare",
    "adjusted",
)
SOURCE_MARKERS = (
    "who",
    "gbd",
    "lancet",
    "cdc",
    "microbiome project",
    "efsa",
    "iarc",
)


def _tier_range(
    central: float,
    notes: str,
    *,
    unit: str,
) -> tuple[float, float, str, str, str]:
    if central == 0.0:
        return (
            0.0,
            0.0,
            "design_constraint",
            "structural_design",
            "Fixed zero representing an explicit model-scope or ecological design constraint.",
        )

    notes_lower = notes.lower()
    if any(marker in notes_lower for marker in PLACEHOLDER_MARKERS):
        lower_factor, upper_factor = 0.25, 4.0
        interval_kind = "expert_plausible_range"
        provenance_class = "expert_best_guess_placeholder"
        method = "Broad 0.25x-4x multiplicative range for an explicitly uncertain placeholder."
    elif any(marker in notes_lower for marker in SOURCE_MARKERS):
        lower_factor, upper_factor = 0.5, 2.0
        interval_kind = "expert_plausible_range"
        provenance_class = "source_informed_transformed"
        method = "0.5x-2x multiplicative range around a source-informed central estimate."
    else:
        lower_factor, upper_factor = 1.0 / 3.0, 3.0
        interval_kind = "expert_plausible_range"
        provenance_class = "expert_best_guess_placeholder"
        method = "Broad one-third-to-threefold range where source precision is unresolved."

    lower = _round_significant(central * lower_factor)
    upper = _round_significant(central * upper_factor)
    if unit == "proportion":
        upper = min(1.0, upper)
    return lower, upper, interval_kind, provenance_class, method


def _headline_rows() -> list[dict[str, Any]]:
    payload = json.loads((DATA_ROOT / "calibration_targets.json").read_text(encoding="utf-8"))
    rows: list[dict[str, Any]] = []
    for metric in payload["headline_metrics"]:
        key = str(metric["key"])
        central = float(metric["target"])
        spec = HEADLINE_RANGES[key]
        rows.append(
            _make_row(
                family="headline",
                key=key,
                label=str(metric["label"]),
                central=central,
                lower=float(spec["lower"]),
                upper=float(spec["upper"]),
                unit=str(metric["unit"]),
                interval_kind=str(spec["interval_kind"]),
                provenance_class=str(spec["provenance_class"]),
                source_id=str(spec["source_id"]),
                source_url=str(spec["source_url"]),
                range_method=str(spec["range_method"]),
                rationale=str(spec["rationale"]),
                last_reviewed=str(spec.get("last_reviewed", LAST_REVIEWED)),
            )
        )
    return rows


def _drug_class_rows() -> list[dict[str, Any]]:
    target_path = DATA_ROOT / "drug_class_share_history_targets.csv"
    targets = pd.read_csv(target_path)
    rows: list[dict[str, Any]] = []
    for _, target in targets.iterrows():
        key = str(target["Class"])
        central = float(target["Share_2025 (%)"])
        if key not in DRUG_CLASS_RANGES:
            raise KeyError(f"Missing drug-class plausible range: {key}")
        lower, upper = DRUG_CLASS_RANGES[key]
        rows.append(
            _make_row(
                family="drug_class_share",
                key=key,
                label=key,
                central=central,
                lower=lower,
                upper=upper,
                unit="percent",
                interval_kind="expert_plausible_range",
                provenance_class="source_informed_transformed",
                source_id="who_glass_amu_2022_and_existing_coarse_class_ranges",
                source_url=WHO_GLASS_2022_URL,
                range_method=(
                    "Marginal expert range informed by heterogeneous national class-use "
                    "data and data/number_on_drug_by_class.csv."
                ),
                rationale=(
                    "The 28 class shares form one composition, but these marginal ranges "
                    "are not a joint confidence region and are not required to sum to 100%."
                ),
            )
        )
    return rows


def _burden_rows(
    *,
    family: str,
    filename: str,
    value_column: str,
    unit: str,
    overrides: dict[str, tuple[float, float]],
) -> list[dict[str, Any]]:
    targets = pd.read_csv(DATA_ROOT / filename)
    rows: list[dict[str, Any]] = []
    for _, target in targets.iterrows():
        key = str(target["Bacteria"]).strip()
        central = float(target[value_column])
        notes = str(target.get("notes", "") or "")
        lower, upper, interval_kind, provenance_class, method = _tier_range(
            central,
            notes,
            unit=unit,
        )
        source_id = "canonical_target_note_source_unresolved"
        source_url = ""
        rationale = (
            "Range tier follows the provenance wording in the canonical target note; "
            "it is a display-only plausible range, not a statistical confidence interval."
        )

        explicit_lower = pd.to_numeric(target.get("plausible_lower"), errors="coerce")
        explicit_upper = pd.to_numeric(target.get("plausible_upper"), errors="coerce")
        has_explicit_range = bool(
            pd.notna(explicit_lower) and pd.notna(explicit_upper)
        )
        if has_explicit_range:
            lower = float(explicit_lower)
            upper = float(explicit_upper)
            interval_kind = str(target.get("interval_kind") or interval_kind)
            provenance_class = str(
                target.get("provenance_class") or provenance_class
            )
            source_id = str(target.get("source_id") or source_id)
            source_url_value = target.get("source_url_or_doi")
            source_url = "" if pd.isna(source_url_value) else str(source_url_value).strip()
            mapping_method = str(target.get("mapping_method") or "").strip()
            if interval_kind == "published_uncertainty_range":
                method = "Published GBD 2019 95% uncertainty interval."
            elif source_id == "gbd_33_bacterial_pathogens_2019":
                method = (
                    "GBD 2019 estimate and interval transferred using the documented "
                    f"model mapping ({mapping_method})."
                )
            else:
                method = (
                    "Explicit review-informed plausible range from the canonical "
                    "target file."
                )
            rationale = notes

        if not has_explicit_range and key in overrides:
            lower, upper = overrides[key]
            interval_kind = "derived_plausible_range"
            provenance_class = "source_informed_transformed"
            method = "Explicit source-informed range override for this target."
            source_id = "canonical_target_note_explicit_range"
            rationale = (
                "The canonical target note contains a quantitative range or a sufficiently "
                "specific source estimate to use an explicit, rounded plausible interval."
            )

        if family == "infection_incidence" and key in WHO_FERG_2021_INCIDENCE_KEYS:
            source_id = "who_ferg_2021_diarrhoeal_hazards"
            source_url = WHO_FERG_2021_DIARRHOEAL_URL
            method = (
                "Published 2021 illness uncertainty bounds divided by 8.2 billion and "
                "treated as a source-informed model-category range."
            )
            rationale = notes

        if (
            not has_explicit_range
            and family == "infection_deaths"
            and key == "vibrio cholerae"
        ):
            interval_kind = "published_uncertainty_range"
            provenance_class = "published_source_range"
            source_id = "who_cholera_global_deaths"
            source_url = WHO_CHOLERA_URL
            method = "WHO-published global annual death range."
            rationale = "WHO reports 21,000-143,000 global cholera deaths per year."
        elif family == "infection_incidence" and key == "vibrio cholerae":
            source_id = "who_cholera_global_cases"
            source_url = WHO_CHOLERA_URL
            method = (
                "WHO 1.3-4.0 million annual case range divided by 8.2 billion, "
                "with a rounded upper bound that contains the existing 4.1 million target."
            )

        rows.append(
            _make_row(
                family=family,
                key=key,
                label=key,
                central=central,
                lower=float(lower),
                upper=float(upper),
                unit=unit,
                interval_kind=interval_kind,
                provenance_class=provenance_class,
                source_id=source_id,
                source_url=source_url,
                range_method=method,
                rationale=rationale,
                last_reviewed=(
                    "2026-08-18"
                    if family == "infection_incidence"
                    and key in WHO_FERG_2021_INCIDENCE_KEYS
                    else LAST_REVIEWED
                ),
            )
        )
    return rows


def build_rows() -> list[dict[str, Any]]:
    rows = _headline_rows()
    rows.extend(_drug_class_rows())
    rows.extend(
        _burden_rows(
            family="infection_incidence",
            filename="infection_incidence_by_bacteria.csv",
            value_column="annual_infection_proportion",
            unit="proportion",
            overrides=INCIDENCE_RANGE_OVERRIDES,
        )
    )
    rows.extend(
        _burden_rows(
            family="carriage",
            filename="microbiome_carriage_by_bacteria.csv",
            value_column="carriage_proportion",
            unit="proportion",
            overrides=CARRIAGE_RANGE_OVERRIDES,
        )
    )
    rows.extend(
        _burden_rows(
            family="infection_deaths",
            filename="deaths_by_bacteria.csv",
            value_column="annual_deaths_millions",
            unit="millions",
            overrides=DEATH_RANGE_OVERRIDES,
        )
    )
    return rows


def main() -> None:
    rows = build_rows()
    expected = 4 + 28 + 42 + 42 + 42
    if len(rows) != expected:
        raise ValueError(f"Expected {expected} target-range rows, found {len(rows)}")

    keys = [(row["target_family"], row["target_key"], row["target_year"]) for row in rows]
    if len(keys) != len(set(keys)):
        raise ValueError("Duplicate target-family/key/year rows in plausible-range registry")

    with OUTPUT_PATH.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=FIELDNAMES)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} target ranges to {OUTPUT_PATH}")
    print(f"Additional antibiotic-use context: {GLOBAL_ANTIBIOTIC_USE_URL}")


if __name__ == "__main__":
    main()
