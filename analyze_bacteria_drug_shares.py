from __future__ import annotations

import argparse
from collections import defaultdict
from pathlib import Path
from typing import Iterable

import pandas as pd


DRUG_CLASS_OVERRIDES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("Penicillins (J01C)", ("penicillin_g", "ampicillin", "amoxicillin", "piperacillin", "ticarcillin", "flucloxacillin")),
    (
        "Beta-lactamase combinations (J01CR)",
        (
            "amoxicillin_clavulanate",
            "piperacillin_tazobactam",
            "ampicillin_sulbactam",
            "ticarcillin_clavulanate",
        ),
    ),
    ("Cephalosporins 1-2G", ("cephalexin", "cefazolin", "cefuroxime")),
    ("Cephalosporins 3G", ("ceftriaxone", "ceftazidime", "cefixime")),
    ("Cephalosporins 3G/BLI", ("ceftolozane_tazobactam",)),
    ("Cephalosporins 4G", ("cefepime",)),
    ("Anti-MRSA Cephalosporins (5G)", ("ceftaroline",)),
    ("Siderophore Cephalosporins", ("cefiderocol",)),
    ("Novel BL/BLI", ("ceftazidime_avibactam", "meropenem_vaborbactam", "aztreonam_avibactam")),
    ("Monobactams", ("aztreonam",)),
    ("Carbapenems (J01DH)", ("meropenem", "imipenem_c", "ertapenem")),
    ("Macrolides (J01F)", ("erythromycin", "azithromycin", "clarithromycin")),
    ("Lincosamides (J01FF)", ("clindamycin",)),
    ("Fluoroquinolones (J01M)", ("ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin")),
    ("Aminoglycosides (J01G)", ("gentamicin", "tobramycin", "amikacin")),
    ("Tetracyclines (J01A)", ("tetracycline", "doxycycline", "minocycline", "tigecycline")),
    ("Sulfonamides (J01E)", ("trim_sulf", "sulfanilamide")),
    ("Glycopeptides (J01XA)", ("vancomycin", "teicoplanin")),
    ("Lipoglycopeptides", ("dalbavancin",)),
    ("Oxazolidinones (J01XX)", ("linezolid", "tedizolid")),
    ("Lipopeptides (J01XX09)", ("daptomycin",)),
    ("Polymyxins (J01XB)", ("colistin",)),
    ("Rifamycins (J04AB)", ("rifampicin",)),
    ("Chloramphenicol (J01BA)", ("chloramphenicol",)),
    ("Nitrofurans (J01XE)", ("nitrofurantoin", "furazolidone")),
    ("Fosfomycin (J01XX01)", ("fosfomycin",)),
    ("Fusidic acid (J01XC)", ("fusidic_a",)),
    ("Pleuromutilins", ("retapamulin",)),
    ("Streptogramins (J01FG)", ("quinu_dalfo",)),
    ("Nitroimidazoles", ("metronidazole",)),
    ("Fidaxomicin", ("fidaxomicin",)),
)

DEFAULT_BACTERIA = (
    "escherichia_coli",
    "klebsiella_pneumoniae",
    "enterobacter_spp.",
    "enterobacter_cloacae",
    "citrobacter_spp.",
    "morganella_spp.",
)

SULFONAMIDE_DRUGS = ("trim_sulf", "sulfanilamide")
MONOBACTAM_DRUGS = ("aztreonam",)
CONTEXT_GROUPS: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("Sulfonamides", SULFONAMIDE_DRUGS),
    ("Monobactams", MONOBACTAM_DRUGS),
)


def build_drug_class_lookup() -> dict[str, str]:
    lookup: dict[str, str] = {}
    for label, slugs in DRUG_CLASS_OVERRIDES:
        for slug in slugs:
            lookup[slug] = label
    return lookup


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Summarize per-bacteria drug and drug-class shares from a simulation_summary CSV. "
            "Uses the per-bacteria currently_on_drug columns already exported by the Rust model."
        )
    )
    parser.add_argument("csv_path", type=Path, help="Path to simulation_summary_*.csv")
    parser.add_argument(
        "--start-year",
        type=float,
        default=2022.0,
        help="Inclusive start year for the averaging window (default: 2022)",
    )
    parser.add_argument(
        "--end-year",
        type=float,
        default=2026.0,
        help="Exclusive end year for the averaging window (default: 2026, i.e. 2022-2025 inclusive)",
    )
    parser.add_argument(
        "--top-drugs",
        type=int,
        default=10,
        help="How many drugs to show per organism (default: 10)",
    )
    parser.add_argument(
        "--bacteria",
        nargs="*",
        default=list(DEFAULT_BACTERIA),
        help="Canonical bacteria slugs to summarize (default: main Enterobacterales culprits)",
    )
    return parser.parse_args()


def select_year_window(df: pd.DataFrame, start_year: float, end_year: float) -> pd.DataFrame:
    if "time_in_years" in df.columns:
        year_values = pd.to_numeric(df["time_in_years"], errors="coerce")
        # Rust export stores elapsed years since SIMULATION_START_YEAR (1930.0), not absolute year.
        if float(year_values.max()) < 500.0:
            df = df.copy()
            df["calendar_year"] = 1930.0 + year_values
            year_col = "calendar_year"
        else:
            year_col = "time_in_years"
    elif "year" in df.columns:
        year_col = "year"
    else:
        raise ValueError("CSV is missing both 'time_in_years' and 'year' columns")

    window_df = df.loc[(df[year_col] >= start_year) & (df[year_col] < end_year)].copy()
    if window_df.empty:
        raise ValueError(
            f"No rows found in requested window {start_year:.0f}-{end_year:.0f} using column '{year_col}'"
        )
    return window_df


def get_global_drug_columns(df: pd.DataFrame) -> list[str]:
    return [col for col in df.columns if col.endswith("_currently_on_drug") and "_currently_on_drug_" not in col]


def get_bacteria_from_columns(df: pd.DataFrame) -> list[str]:
    suffix = "_currently_infected"
    excluded_prefixes = {"total"}
    bacteria = []
    for col in df.columns:
        if not col.endswith(suffix):
            continue
        slug = col[: -len(suffix)]
        if slug in excluded_prefixes:
            continue
        bacteria.append(slug)
    return sorted(bacteria)


def normalize_display(name: str) -> str:
    return name.replace("_", " ").title()


def compute_drug_share_table(
    df: pd.DataFrame,
    bacteria: str,
    drugs: Iterable[str],
) -> pd.DataFrame:
    records: list[dict[str, object]] = []
    for drug in drugs:
        column = f"{bacteria}_currently_on_drug_{drug}"
        if column not in df.columns:
            continue
        mean_count = float(df[column].mean())
        records.append({"drug": drug, "mean_count": mean_count})

    result = pd.DataFrame(records)
    if result.empty:
        return result

    total = float(result["mean_count"].sum())
    if total <= 0.0:
        result["share_pct"] = 0.0
    else:
        result["share_pct"] = result["mean_count"] / total * 100.0
    return result.sort_values(["share_pct", "drug"], ascending=[False, True]).reset_index(drop=True)


def compute_class_share_table(drug_table: pd.DataFrame, class_lookup: dict[str, str]) -> pd.DataFrame:
    if drug_table.empty:
        return pd.DataFrame(columns=["class_label", "share_pct"])

    class_totals: dict[str, float] = defaultdict(float)
    for row in drug_table.itertuples(index=False):
        class_label = class_lookup.get(row.drug, "Other / unspecified")
        class_totals[class_label] += float(row.share_pct)

    class_df = pd.DataFrame(
        [{"class_label": label, "share_pct": share_pct} for label, share_pct in class_totals.items()]
    )
    return class_df.sort_values(["share_pct", "class_label"], ascending=[False, True]).reset_index(drop=True)


def print_table(title: str, rows: pd.DataFrame, name_col: str, top_n: int | None = None) -> None:
    print(title)
    if rows.empty:
        print("  no matching data")
        return

    subset = rows if top_n is None else rows.head(top_n)
    for row in subset.itertuples(index=False):
        row_map = row._asdict()
        label = row_map[name_col]
        print(f"  {label:<35} {row_map['share_pct']:6.2f}%")


def build_sulfonamide_by_bacteria_table(df: pd.DataFrame, drugs: Iterable[str]) -> pd.DataFrame:
    available_drugs = [drug for drug in SULFONAMIDE_DRUGS if drug in set(drugs)]
    if not available_drugs:
        return pd.DataFrame(
            columns=[
                "bacteria",
                "mean_infected_count",
                "mean_infected_and_on_drug_count",
                "sulfonamide_mean_count",
                "share_pct",
            ]
        )

    records: list[dict[str, object]] = []
    for bacteria in get_bacteria_from_columns(df):
        infected_col = f"{bacteria}_currently_infected"
        if infected_col not in df.columns:
            continue

        drug_table = compute_drug_share_table(df, bacteria, drugs)
        total_mean_on_drug = float(drug_table["mean_count"].sum()) if not drug_table.empty else 0.0
        sulfonamide_mean_count = 0.0
        if not drug_table.empty:
            sulfonamide_mean_count = float(
                drug_table.loc[drug_table["drug"].isin(available_drugs), "mean_count"].sum()
            )

        share_pct = 0.0 if total_mean_on_drug <= 0.0 else sulfonamide_mean_count / total_mean_on_drug * 100.0
        records.append(
            {
                "bacteria": bacteria,
                "mean_infected_count": float(df[infected_col].mean()),
                "mean_infected_and_on_drug_count": total_mean_on_drug,
                "sulfonamide_mean_count": sulfonamide_mean_count,
                "share_pct": share_pct,
            }
        )

    result = pd.DataFrame(records)
    if result.empty:
        return result
    return result.sort_values(
        ["share_pct", "sulfonamide_mean_count", "mean_infected_and_on_drug_count", "bacteria"],
        ascending=[False, False, False, True],
    ).reset_index(drop=True)


def print_sulfonamide_by_bacteria_table(rows: pd.DataFrame) -> None:
    print("Sulfonamide share by bacteria (all bacteria)\n")
    if rows.empty:
        print("no sulfonamide-capable drugs found in CSV")
        return

    print(
        f"{'bacteria':<32} {'sulfonamide_share':>18} {'mean_sulfa_on_drug':>20} {'mean_on_any_drug':>18} {'mean_infected':>16}"
    )
    for row in rows.itertuples(index=False):
        row_map = row._asdict()
        print(
            f"{row_map['bacteria']:<32} "
            f"{row_map['share_pct']:17.2f}% "
            f"{row_map['sulfonamide_mean_count']:20.2f} "
            f"{row_map['mean_infected_and_on_drug_count']:18.2f} "
            f"{row_map['mean_infected_count']:16.2f}"
        )


def build_context_split_table(df: pd.DataFrame, drugs: Iterable[str]) -> pd.DataFrame:
    available_drugs = set(drugs)
    bacteria = get_bacteria_from_columns(df)
    records: list[dict[str, object]] = []

    mean_currently_taking_any_drug = float(df["currently_taking_drug_count"].mean())
    mean_currently_infected_and_on_drug = float(df["currently_infected_and_on_drug_count"].mean())

    for label, group_drugs in CONTEXT_GROUPS:
        present_group_drugs = [drug for drug in group_drugs if drug in available_drugs]
        global_cols = [f"{drug}_currently_on_drug" for drug in present_group_drugs if f"{drug}_currently_on_drug" in df.columns]
        infected_cols = [
            f"{bacteria_slug}_currently_on_drug_{drug}"
            for bacteria_slug in bacteria
            for drug in present_group_drugs
            if f"{bacteria_slug}_currently_on_drug_{drug}" in df.columns
        ]

        global_mean_count = float(df[global_cols].mean().sum()) if global_cols else 0.0
        infected_mean_count = float(df[infected_cols].mean().sum()) if infected_cols else 0.0
        residual_outside_infected_context = global_mean_count - infected_mean_count

        records.append(
            {
                "group_label": label,
                "drugs_present": ", ".join(present_group_drugs) if present_group_drugs else "none",
                "global_mean_count": global_mean_count,
                "infected_context_mean_count": infected_mean_count,
                "residual_outside_infected_context": residual_outside_infected_context,
                "global_share_of_all_drug_use_pct": 0.0 if mean_currently_taking_any_drug <= 0.0 else global_mean_count / mean_currently_taking_any_drug * 100.0,
                "infected_share_of_infected_drug_use_pct": 0.0 if mean_currently_infected_and_on_drug <= 0.0 else infected_mean_count / mean_currently_infected_and_on_drug * 100.0,
            }
        )

    return pd.DataFrame(records)


def build_noninfected_proxy_table(df: pd.DataFrame, drugs: Iterable[str]) -> tuple[pd.DataFrame, dict[str, float]]:
    available_drugs = set(drugs)
    bacteria = get_bacteria_from_columns(df)

    mean_currently_taking_any_drug = float(df["currently_taking_drug_count"].mean())
    mean_currently_infected_and_on_drug = float(df["currently_infected_and_on_drug_count"].mean())
    mean_noninfected_on_drug_proxy = mean_currently_taking_any_drug - mean_currently_infected_and_on_drug
    mean_new_drug_initiations = float(df["new_drug_initiations_count"].mean())
    mean_new_drug_initiations_infected = float(df["new_drug_initiations_count_infected"].mean())
    mean_new_drug_initiations_noninfected_proxy = (
        mean_new_drug_initiations - mean_new_drug_initiations_infected
    )
    mean_immunosuppressed = float(df["number_severely_immunosuppressed"].mean())

    records: list[dict[str, object]] = []
    for label, group_drugs in CONTEXT_GROUPS:
        present_group_drugs = [drug for drug in group_drugs if drug in available_drugs]
        global_cols = [
            f"{drug}_currently_on_drug" for drug in present_group_drugs if f"{drug}_currently_on_drug" in df.columns
        ]
        infected_cols = [
            f"{bacteria_slug}_currently_on_drug_{drug}"
            for bacteria_slug in bacteria
            for drug in present_group_drugs
            if f"{bacteria_slug}_currently_on_drug_{drug}" in df.columns
        ]

        global_mean_count = float(df[global_cols].mean().sum()) if global_cols else 0.0
        infected_mean_count = float(df[infected_cols].mean().sum()) if infected_cols else 0.0
        noninfected_proxy_count = global_mean_count - infected_mean_count

        records.append(
            {
                "group_label": label,
                "drugs_present": ", ".join(present_group_drugs) if present_group_drugs else "none",
                "noninfected_proxy_mean_count": noninfected_proxy_count,
                "share_of_noninfected_drug_stock_pct": 0.0
                if mean_noninfected_on_drug_proxy <= 0.0
                else noninfected_proxy_count / mean_noninfected_on_drug_proxy * 100.0,
                "share_of_immunosuppressed_pool_pct": 0.0
                if mean_immunosuppressed <= 0.0
                else noninfected_proxy_count / mean_immunosuppressed * 100.0,
            }
        )

    summary_metrics = {
        "mean_immunosuppressed": mean_immunosuppressed,
        "mean_currently_taking_any_drug": mean_currently_taking_any_drug,
        "mean_currently_infected_and_on_drug": mean_currently_infected_and_on_drug,
        "mean_noninfected_on_drug_proxy": mean_noninfected_on_drug_proxy,
        "mean_new_drug_initiations": mean_new_drug_initiations,
        "mean_new_drug_initiations_infected": mean_new_drug_initiations_infected,
        "mean_new_drug_initiations_noninfected_proxy": mean_new_drug_initiations_noninfected_proxy,
    }

    return pd.DataFrame(records), summary_metrics


def print_context_split_table(rows: pd.DataFrame) -> None:
    print("Global vs infected-context use\n")
    if rows.empty:
        print("no context split rows available")
        return

    for row in rows.itertuples(index=False):
        row_map = row._asdict()
        print(f"{row_map['group_label']} ({row_map['drugs_present']})")
        print(f"  global mean currently_on_drug: {row_map['global_mean_count']:.2f}")
        print(f"  infected-context mean summed over bacteria: {row_map['infected_context_mean_count']:.2f}")
        print(f"  residual outside infected-context accounting: {row_map['residual_outside_infected_context']:.2f}")
        print(f"  global share of all drug use: {row_map['global_share_of_all_drug_use_pct']:.2f}%")
        print(f"  infected-context share of infected-and-on-drug use: {row_map['infected_share_of_infected_drug_use_pct']:.2f}%")
        print()


def print_noninfected_proxy_table(rows: pd.DataFrame, metrics: dict[str, float]) -> None:
    print("Non-infected / prophylaxis proxy\n")
    print(f"  mean severely immunosuppressed: {metrics['mean_immunosuppressed']:.2f}")
    print(f"  mean currently taking any drug: {metrics['mean_currently_taking_any_drug']:.2f}")
    print(f"  mean currently infected and on drug: {metrics['mean_currently_infected_and_on_drug']:.2f}")
    print(f"  mean non-infected on-drug proxy: {metrics['mean_noninfected_on_drug_proxy']:.2f}")
    print(f"  mean new drug initiations: {metrics['mean_new_drug_initiations']:.2f}")
    print(f"  mean new drug initiations in infected people: {metrics['mean_new_drug_initiations_infected']:.2f}")
    print(
        "  mean new drug initiations outside infected people (proxy): "
        f"{metrics['mean_new_drug_initiations_noninfected_proxy']:.2f}"
    )
    print()

    if rows.empty:
        print("no non-infected proxy rows available")
        return

    for row in rows.itertuples(index=False):
        row_map = row._asdict()
        print(f"{row_map['group_label']} ({row_map['drugs_present']})")
        print(f"  non-infected stock proxy: {row_map['noninfected_proxy_mean_count']:.2f}")
        print(
            "  share of all non-infected on-drug stock: "
            f"{row_map['share_of_noninfected_drug_stock_pct']:.2f}%"
        )
        print(
            "  share of severely immunosuppressed pool: "
            f"{row_map['share_of_immunosuppressed_pool_pct']:.2f}%"
        )
        print()


def main() -> None:
    args = parse_args()
    if not args.csv_path.exists():
        raise FileNotFoundError(f"CSV not found: {args.csv_path}")

    df = pd.read_csv(args.csv_path, low_memory=False)
    df = select_year_window(df, args.start_year, args.end_year)

    global_drug_columns = get_global_drug_columns(df)
    drugs = [col[: -len("_currently_on_drug")] for col in global_drug_columns]
    if not drugs:
        raise ValueError("No '*_currently_on_drug' columns found in the CSV")

    class_lookup = build_drug_class_lookup()
    sulfonamide_table = build_sulfonamide_by_bacteria_table(df, drugs)
    context_split_table = build_context_split_table(df, drugs)
    noninfected_proxy_table, noninfected_proxy_metrics = build_noninfected_proxy_table(df, drugs)

    print(
        f"Per-bacteria prescribing shares from {args.csv_path.name} "
        f"for {args.start_year:.0f}-{args.end_year - 1:.0f}\n"
    )

    print_context_split_table(context_split_table)
    print_noninfected_proxy_table(noninfected_proxy_table, noninfected_proxy_metrics)

    print_sulfonamide_by_bacteria_table(sulfonamide_table)
    print()

    for bacteria in args.bacteria:
        infected_col = f"{bacteria}_currently_infected"
        if infected_col not in df.columns:
            print(f"{bacteria}: missing infected column, skipped\n")
            continue

        drug_table = compute_drug_share_table(df, bacteria, drugs)
        class_table = compute_class_share_table(drug_table, class_lookup)

        mean_infected = float(df[infected_col].mean())
        total_mean_on_drug = float(drug_table["mean_count"].sum()) if not drug_table.empty else 0.0

        print(f"{bacteria} ({normalize_display(bacteria)})")
        print(f"  mean infected count in window: {mean_infected:.2f}")
        print(f"  mean infected-and-on-drug count: {total_mean_on_drug:.2f}")
        print_table("  drug shares:", drug_table, "drug", top_n=args.top_drugs)
        print_table("  class shares:", class_table, "class_label")
        print()


if __name__ == "__main__":
    main()