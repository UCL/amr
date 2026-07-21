"""
make_paper_tables.py

Generate the paper-facing HTML outputs from one or more calibration_summary_*.txt
files: Table 1, Supplementary Table S2, main Figures 1-13,
Supplementary Figures S1-S3 and S5-S8, and diagnostic Supplementary Figure SX.

Usage
-----
Single run:
    python amr_simulation_output_analysis/make_paper_tables.py output_graphs/calibration_summary_958282.txt

Multiple runs (pass explicit paths or a glob):
    python amr_simulation_output_analysis/make_paper_tables.py output_graphs/calibration_summary_*.txt

Package/module form:
    python -m amr_simulation_output_analysis.make_paper_tables output_graphs/calibration_summary_*.txt

Figure 2 summary mode:
    Edit FIGURE2_SUMMARY_MODE below: "median_range" or "mean_ci".

Output
------
paper_tables/
    index.html
    Tables/
        T1__model_summary.html
        Supplementary_Table_S2__detailed_bacterium_drug_resistance_benchmarks.html
    Figures/
        Figure_1__calibration_headline_metrics.html/.png/.svg
        Figure_2__calibration_resistance_fit_by_bacteria_drug_class.html/.png/.svg
        Figure_3__calibration_drug_class_share.html/.png/.svg
        Figure_4__calibration_infection_deaths_by_bacteria.html/.png/.svg
        Figure_5__calibration_carriage_prevalence_by_bacteria.html/.png/.svg
        Figure_6A__resistance_trends.html/.png/.svg
        Figure_6B__resistance_trends_by_bacterium.html/.png/.svg
        Figure_6C__serious_r_trends_by_bacterium.html/.png/.svg
        Figure_7__serious_r_by_hospital_community.html/.png/.svg
        Figure_8__infection_death_rate_by_region.html/.png/.svg
        Figure_9__antibiotic_use_by_treatment_context.html/.png/.svg
        Figure_10__sepsis_context_effective_therapy.html/.png/.svg
        Figure_11__activity_retained_by_bacterium.html/.png/.svg
        Figure_12__distribution_drug_use_by_bacteria.html/.png/.svg
        Figure_13__resistance_pathway_counterfactuals.html/.png/.svg
        Supplementary_Figure_S1__potential_activity_retained.html/.png/.svg
        Supplementary_Figure_S2__microbiome_resistance_reservoir.html/.png/.svg
        Supplementary_Figure_S3__carrier_vs_non_carrier_infection_incidence.html/.png/.svg
        Supplementary_Figure_S5__diagnostic_testing_targeted_treatment_cascade.html/.png/.svg
        Supplementary_Figure_S6__new_active_infection_denominators_by_bacterium.html/.png/.svg
        Supplementary_Figure_S7__active_infection_incidence_by_bacterium.html/.png/.svg
        Supplementary_Figure_S8__infection_outcome_pathway_by_bacterium.html/.png/.svg
        Supplementary_Figure_SX__modelled_resistance_mechanisms_by_bacterium.html/.png/.svg
"""

from __future__ import annotations

import glob
import io
import json
import math
import re
import shutil
import sys
from pathlib import Path
from typing import Callable, Union

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.ticker as mticker
import numpy as np
import pandas as pd

try:
    from .parse_calibration import aggregate, parse_files
except ImportError:  # Allows direct script execution from this folder.
    from parse_calibration import aggregate, parse_files

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parent
OUT_DIR = REPO_ROOT / "paper_tables"
CALIBRATION_TARGETS_PATH = REPO_ROOT / "data" / "calibration_targets.json"
SIMULATION_OUTPUTS_DIR = REPO_ROOT / "amr_simulation_output_analysis_outputs"

# Figure 2 toggle. Options:
#   "median_range" - simulation median with 5th-95th percentile range
#   "mean_ci"      - simulation mean with 95% confidence interval
FIGURE2_SUMMARY_MODE = "mean_ci"

TABLES_DIRNAME = "Tables"
FIGURES_DIRNAME = "Figures"
GENERATED_SUBDIRS = (
    "main",
    "supplementary",
    "Supplementary",
    "figures",
    TABLES_DIRNAME,
    FIGURES_DIRNAME,
)

# ---------------------------------------------------------------------------
# Shared HTML / CSS
# ---------------------------------------------------------------------------

_CSS = """
body {
  font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
  line-height: 1.5; max-width: 1500px; margin: 0 auto;
  padding: 24px; color: #222; font-size: 13px;
}
h1  { font-size: 1.25em; border-bottom: 2px solid #2c3e50;
      padding-bottom: 8px; color: #2c3e50; margin-bottom: 4px; }
h2  { font-size: 1.05em; margin-top: 26px; color: #2c3e50;
      border-bottom: 1px solid #ddd; padding-bottom: 4px; }
.subtitle { color: #555; font-size: 0.9em; margin-bottom: 16px; }
table { border-collapse: collapse; width: 100%; margin: 10px 0 18px 0;
        font-size: 12px; box-shadow: 0 1px 3px rgba(0,0,0,0.08); }
th, td { border: 1px solid #ccc; padding: 5px 9px; text-align: right; white-space: nowrap; }
th:first-child, td:first-child { text-align: left; white-space: normal; min-width: 160px; }
th:nth-child(2), td:nth-child(2) { text-align: left; white-space: normal; }
th { background: #eef1f5; font-weight: 600; color: #2c3e50; }
tr:nth-child(even) td { background: #fafbfc; }
tr:hover td { background: #eef4fb; }
.total-row td { font-weight: 600; background: #e2e8f0 !important; }
.note  { font-style: italic; color: #666; font-size: 0.87em; margin: 4px 0; }
.meta-box {
  background: #f4f7fb; border-left: 4px solid #3498db;
  padding: 9px 14px; margin-bottom: 14px;
  border-radius: 0 4px 4px 0; font-size: 0.88em; color: #333;
}
ol.footnotes { font-size: 0.87em; color: #444; padding-left: 22px; margin-top: 20px; }
ol.footnotes li { margin-bottom: 7px; line-height: 1.5; }
.back-link { font-size: 0.83em; margin-bottom: 12px; }
a { color: #2980b9; }
"""


def _html_head(title: str) -> str:
    return (
        f'<!DOCTYPE html>\n<html lang="en">\n<head>'
        f'<meta charset="UTF-8"><title>{title}</title>'
        f"<style>{_CSS}</style></head>\n<body>\n"
    )


def _html_table(df: pd.DataFrame, total_marker: str | None = None) -> str:
    """Render a DataFrame as an HTML <table>."""
    if df is None or df.empty:
        return "<p class='note'>No data available.</p>"
    lines = ["<table>", "<thead><tr>"]
    for col in df.columns:
        lines.append(f"<th>{col}</th>")
    lines.append("</tr></thead><tbody>")
    for _, row in df.iterrows():
        vals = list(row)
        first_str = str(vals[0]).lower() if vals else ""
        is_total = total_marker and total_marker.lower() in first_str
        cls = " class='total-row'" if is_total else ""
        lines.append(f"<tr{cls}>")
        for v in vals:
            cell = str(v) if v not in (None, np.nan) and str(v) not in ("nan", "None") else "—"
            lines.append(f"<td>{cell}</td>")
        lines.append("</tr>")
    lines.append("</tbody></table>")
    return "\n".join(lines)


def _html_footnotes(notes: list[str]) -> str:
    if not notes:
        return ""
    items = "".join(f"<li>{n}</li>" for n in notes)
    return f"<ol class='footnotes'>\n{items}\n</ol>"


def _meta_box(agg: dict) -> str:
    m = agg.get("meta", {})
    n = agg.get("n_runs", 1)
    run_label = f"{n} accepted calibration run{'s' if n > 1 else ''}"
    parts = [
        f"<strong>Target year:</strong> {m.get('target_year', '—')}",
        f"<strong>Window:</strong> {m.get('window_duration', '—')}",
        f"<strong>Simulated population (mean):</strong> {m.get('mean_pop', '—')}",
        f"<strong>Population scale factor:</strong> {m.get('scale_factor', '—')}",
        f"<strong>Runs:</strong> {run_label}",
    ]
    return "<div class='meta-box'>" + " &nbsp;|&nbsp; ".join(parts) + "</div>\n"


def _meta_footnote(agg: dict) -> str:
    m = agg.get("meta", {})
    n = agg.get("n_runs", 1)
    run_label = f"{n} accepted calibration run{'s' if n > 1 else ''}"
    return (
        f"<strong>Run/window:</strong> Target year: {m.get('target_year', '—')}; "
        f"window: {m.get('window_duration', '—')}; "
        f"simulated population (mean): {m.get('mean_pop', '—')}; "
        f"population scale factor: {m.get('scale_factor', '—')}; "
        f"runs: {run_label}."
    )


def _back_link() -> str:
    return "<p class='back-link'><a href='../index.html'>← Back to index</a></p>\n"


def _save(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    print(f"  Saved: {path}")


def _resolve_project_path(path: Union[str, Path]) -> Path:
    """Resolve relative paths from cwd first, then from the repository root."""
    candidate = Path(path)
    if candidate.is_absolute():
        return candidate
    if candidate.exists():
        return candidate
    root_candidate = REPO_ROOT / candidate
    if root_candidate.exists():
        return root_candidate
    return candidate


def _prepare_output_dirs(out_dir: Path) -> None:
    """Remove only known generated paper output folders, then recreate the new layout."""
    out_dir.mkdir(parents=True, exist_ok=True)
    root = out_dir.resolve()
    for name in GENERATED_SUBDIRS:
        target = (out_dir / name).resolve()
        if target.parent == root and target.exists() and target.is_dir():
            shutil.rmtree(target)
            print(f"  Removed old generated folder: {out_dir / name}")
    (out_dir / TABLES_DIRNAME).mkdir(parents=True, exist_ok=True)
    (out_dir / FIGURES_DIRNAME).mkdir(parents=True, exist_ok=True)


def _save_figure(
    fig: "plt.Figure",
    out_dir: Path,
    stem: str,
    title: str,
    note: str,
    footnotes: list[str],
    subfolder: str = "main",
    agg: dict | None = None,
    extra_html: str = "",
) -> None:
    """Save a matplotlib figure as PNG/SVG and write an HTML wrapper page."""
    fig_dir = out_dir / FIGURES_DIRNAME
    fig_dir.mkdir(parents=True, exist_ok=True)
    png_path = fig_dir / f"{stem}.png"
    svg_path = fig_dir / f"{stem}.svg"
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {png_path}")
    print(f"  Saved: {svg_path}")
    html_rel = f"{stem}.png"
    svg_rel = f"{stem}.svg"
    body  = _html_head(title)
    body += _back_link()
    body += (
        f"<img src='{html_rel}' alt='{stem}' "
        f"style='max-width:100%; border:1px solid #ddd; border-radius:4px;'>\n"
    )
    body += f"<p class='note'>Download: <a href='{html_rel}'>PNG</a> | <a href='{svg_rel}'>SVG</a></p>\n"
    if extra_html:
        body += extra_html
    notes = []
    if agg is not None:
        notes.append(_meta_footnote(agg))
    if note:
        notes.append(note)
    notes.extend(footnotes)
    body += _html_footnotes(notes)
    body += "</body></html>"
    html_path = fig_dir / f"{stem}.html"
    _save(html_path, body)


# Best-estimate targets for % of new infections that are hospital-acquired.
# Sources: ECDC EARS-Net, WHO HAI reports, INICC, published HAI epidemiology ~2015–2024.
# Ranges are broad; central estimates used as calibration targets.
_HA_PCT_TARGETS: dict[str, float] = {
    "acinetobacter baumannii":                   65.0,  # ESKAPE; ICU/VAP ~60–80%
    "bacteroides fragilis":                       30.0,  # post-surgical intra-abdominal ~20–40%
    "bordetella pertussis":                        5.0,  # occasionally nosocomial in neonates
    "burkholderia cepacia complex":               65.0,  # CF centres / CGD ~50–80%
    "campylobacter jejuni":                        2.0,  # foodborne; very rare HA
    "chlamydia trachomatis":                       1.0,  # STI; negligible HA
    "citrobacter spp.":                           45.0,  # opportunistic; device-associated ~40–55%
    "clostridioides difficile":                   50.0,  # classic HAI ~40–60%
    "enterobacter cloacae":                       45.0,  # nosocomial ~40–55%
    "enterobacter spp.":                          45.0,
    "enterococcus faecalis":                      30.0,  # UTI/wound HA ~20–40%
    "enterococcus faecium":                       50.0,  # VRE; BSI/UTI ~40–60%
    "escherichia coli":                           15.0,  # CAUTI/BSI ~10–20%
    "haemophilus influenzae":                     10.0,  # mainly community; HA neonates/elderly
    "helicobacter pylori":                        10.0,  # endoscopy-related seeding possible
    "invasive non-typhoidal salmonella spp.":     10.0,
    "klebsiella pneumoniae":                      40.0,  # HAI ~30–50%
    "legionella pneumophila":                     20.0,  # hospital water systems ~15–30%
    "listeria monocytogenes":                     10.0,  # foodborne; HA in immunocompromised
    "mdr mycobacterium tuberculosis":              8.0,
    "moraxella catarrhalis":                      10.0,
    "morganella spp.":                            40.0,  # UTI/wound HA ~30–50%
    "mycoplasma genitalium":                       1.0,  # STI
    "mycoplasma pneumoniae":                       8.0,  # community; HA in elderly/outbreaks
    "neisseria gonorrhoeae":                       1.0,  # STI
    "neisseria meningitidis":                     20.0,  # HA in infants/elderly ~15–25%
    "proteus spp.":                               30.0,  # catheter-associated ~25–35%
    "providencia stuartii":                       65.0,  # long-term-care catheter ~60–75%
    "pseudomonas aeruginosa":                     45.0,  # VAP/wound HA ~35–55%
    "salmonella enterica serovar paratyphi a":     3.0,
    "salmonella enterica serovar typhi":           3.0,
    "serratia spp.":                              50.0,  # ICU/NICU ~40–60%
    "shigella spp.":                               2.0,  # foodborne/waterborne
    "staphylococcus aureus":                      25.0,  # MRSA/SSTI HA ~20–30%
    "staphylococcus epidermidis":                 75.0,  # device/implant ~70–85%
    "stenotrophomonas maltophilia":               70.0,  # ventilated/immunocompromised ~60–80%
    "streptococcus agalactiae":                   30.0,  # neonatal/obstetric HA ~20–35%
    "streptococcus pneumoniae":                   10.0,  # mostly community; HA in elderly
    "streptococcus pyogenes":                     10.0,
    "treponema pallidum":                          1.0,  # STI
    "vibrio cholerae":                             2.0,  # waterborne
    "yersinia enterocolitica":                     3.0,
}


def _clean_df(
    df: pd.DataFrame,
    *,
    target_label: str = "Observed estimate",
) -> pd.DataFrame:
    """Drop delta columns and replace 'Target' with the requested display label."""
    if df is None or df.empty:
        return df
    df = df.copy()
    drop = [c for c in df.columns if re.search(r'\bDelta\b|\bΔ\b', c, re.IGNORECASE)]
    df = df.drop(columns=drop, errors='ignore')
    def _rename_col(c: str) -> str:
        def repl(m):
            return target_label if m.group(0)[0].isupper() else target_label.lower()
        return re.sub(r'(?i)\btarget\b', repl, c)
    df.columns = [_rename_col(c) for c in df.columns]
    return df


def _window_note(n: int) -> str:
    """Standard simulation window / run-count note."""
    interval = (
        f"Values are medians (5th\u201395th percentile) across {n} accepted calibration runs."
        if n > 1 else "Values are from a single accepted calibration run."
    )
    return (
        "Simulation outputs represent 2025, averaged over a 4-year calibration window "
        f"(2022\u20132025). {interval}"
    )


_HEADLINE_TARGET_SOURCE_NOTES = [
    "Infection-death target: 6.4 million model-scope bacterial infection deaths per year. "
    "This is the rounded sum of per-organism mortality targets after excluding H. pylori "
    "and MDR-TB, matching the Figure 1 simulation numerator. It sits below the 7.7 million "
    "GBD/33-pathogen central estimate and is not the AMR-specific mortality estimate.",
    "Antibiotic-use target: Klein EY et al. (2018). Global increase and geographic "
    "convergence in antibiotic consumption between 2000 and 2015. <em>PNAS</em> "
    "115:E3463-E3470. The model uses a person-prevalence proxy for daily users, rather "
    "than DDD-equivalent consumption.",
    "Bacterial-infection incidence target: Vos T et al. (2020). Global burden of 369 "
    "diseases and injuries in 204 countries and territories, 1990-2019. "
    "<em>Lancet</em> 396:1204-1222.",
    "Sepsis target: Rudd KE et al. (2020). Global, regional, and national sepsis incidence "
    "and mortality, 1990-2017. <em>Lancet</em> 395:200-211. The model target is a "
    "bacterial-subset anchor rather than the full all-cause sepsis estimate.",
]

_RESISTANCE_TARGET_SOURCE_NOTES = [
    "Resistance-prevalence values are evidence-informed calibration benchmarks. ECDC "
    "EARS-Net, WHO GLASS reports, and organism-drug literature informed the legacy matrix, "
    "but cell-level citations and harmonised denominator definitions were not retained.",
    "Conditional mean any_r values are expert-assigned model benchmarks on the model's "
    "unitless resistance scale; they are not MIC values or direct surveillance estimates.",
    "Both benchmark families compare with simulated active-infection person-days and should "
    "not be interpreted as a harmonised global clinical-isolate surveillance dataset.",
]

_DRUG_CLASS_TARGET_SOURCE_NOTES = [
    "Global drug-class share estimates are derived from the WHO AWaRe classification "
    "database, IQVIA MIDAS market data, and ECDC ESAC-Net surveillance, adjusted to "
    "represent the global population mix.",
    "Drug class codes follow the WHO ATC classification system (J01 antibacterials for "
    "systemic use). Classes without a standard J01 code are assigned their closest "
    "available ATC grouping.",
]

_MORTALITY_TARGET_SOURCE_NOTES = [
    "Observed infection-death estimates are based on GBD 2019/2020 cause-of-death "
    "attributions, WHO mortality data, and organism-specific published literature.",
    "Bacterium-level death estimates can sum to more than the headline all-cause bacterial "
    "death estimate because polymicrobial deaths may be attributed to all contributing "
    "pathogens.",
]

_CARRIAGE_TARGET_SOURCE_NOTES = [
    "Carriage target/observed-estimate values are drawn from published cross-sectional "
    "carriage surveys; individual source details are given in the model description.",
    "Carriage values are percentages of the world population carrying the organism "
    "asymptomatically in the modelled microbiome/carriage compartment.",
]

_INFECTION_DEATH_EXCLUDED_BACTERIA_SLUGS = {
    "helicobacter_pylori",
    "mdr_mycobacterium_tuberculosis",
}


def _bacteria_slug_for_filter(value: object) -> str:
    clean_value = re.sub(r"\s+\*$", "", str(value or "").strip().lower())
    return clean_value.replace(" ", "_")


# ---------------------------------------------------------------------------
# Table T1 — Model Summary (hand-written; no agg data required)
# ---------------------------------------------------------------------------

_T1_ROWS: list[tuple[str, str, str]] = [
    # (Section heading, Feature, Detail)
    # ── Model framework ──────────────────────────────────────────────────────
    ("Model framework", "Model type",
     "Individual-based model (IBM / agent-based model). Each simulated person is "
     "an independent agent with their own state (age, sex, region, immune status, "
     "active infections, microbiome, current antibiotics)."),

    ("Model framework", "Simulation time-span",
     "1930–2035 (105 years; 38,325 daily time steps). Starting before the antibiotic "
     "era allows the model to reproduce the full historical arc of drug introduction, "
     "rising consumption, and accumulating resistance."),

    ("Model framework", "Time step",
     "One calendar day. Every living individual is processed through 21 ordered "
     "mechanistic rules each day."),

    ("Model framework", "Population size",
     "100,000 synthetic individuals per run. Results are rescaled to the global population "
     "using a run-specific scale factor derived from calibration targets."),

    ("Model framework", "Geographic scope",
     "Six world regions: North America, Europe, Asia, Oceania, South America, Africa. "
     "Regions differ in antibiotic availability, hospital capacity, diagnostic testing "
     "rates, and pathogen epidemiology."),

    ("Model framework", "Stochasticity",
     "All events (infection, testing, treatment initiation, resistance mutation, death) "
     "are sampled from daily Bernoulli probabilities. Multiple independent runs characterise "
     "the distribution of outcomes; accepted calibration runs form the uncertainty ensemble."),

    # ── Biological scope ─────────────────────────────────────────────────────
    ("Biological scope", "Bacteria modelled",
     "42 species, including all ESKAPE pathogens and IHME 2019 Global Burden of "
     "Antimicrobial Resistance priority organisms."),

    ("Biological scope", "Antibiotics modelled",
     "61 drugs across 31 antibiotic classes (ATC J01 hierarchy plus key non-J01 agents "
     "such as metronidazole, fidaxomicin, and polymyxins)."),

    ("Biological scope", "Resistance mechanisms",
     "40 distinct biochemical mechanisms including β-lactamases (TEM, SHV, CTX-M, NDM, "
     "OXA-type), efflux pumps (MexAB-OprM, AcrAB-TolC, norA), target-site modifications "
     "(PBP2a, GyrA/ParC, rpsL), and porin loss."),

    ("Biological scope", "Drug–bacteria potency matrix",
     "Full 61 × 42 matrix of minimum inhibitory concentration (MIC) shifts, intrinsic "
     "susceptibility flags, and mechanism-specific potency overrides. Each cell is "
     "individually parameterised."),

    # ── Disease processes ────────────────────────────────────────────────────
    ("Disease processes", "Infection acquisition",
     "Three routes: (i) community acquisition (region-, age-, and sex-specific daily "
     "incidence); (ii) hospital acquisition (admission probability × nosocomial hazard "
     "× length-of-stay); (iii) endogenous infection seeded from the microbiome carriage "
     "compartment."),

    ("Disease processes", "Clinical presentation",
     "Nine syndrome categories (respiratory, urinary tract, bloodstream, "
     "skin and soft tissue, gastrointestinal, sexually transmitted, bone and joint, "
     "CNS, other). Syndrome drives empiric drug choice and sepsis risk."),

    ("Disease processes", "Sepsis",
     "Log-odds model combining bacteraemia probability, syndrome severity, age, "
     "immune status, and treatment adequacy. Sepsis substantially increases "
     "mortality risk and triggers escalation of antibiotic therapy."),

    ("Disease processes", "Diagnostic testing",
     "Microbiology culture followed by antimicrobial susceptibility testing (AST). "
     "Testing probability depends on region, syndrome, hospitalisation status, and "
     "clinician ordering behaviour. Region-specific multipliers range from ×0.3 "
     "(Africa) to ×1.2 (Europe)."),

    ("Disease processes", "Antibiotic treatment",
     "Two-phase prescribing: (i) empiric initiation based on syndrome and local "
     "guideline preferences; (ii) culture-guided de-escalation or switching once "
     "AST results arrive. Drug selection factors in availability tier, formulary "
     "restrictions, and organism-specific access rules. Course duration, toxicity, "
     "and treatment failure are explicitly modelled."),

    ("Disease processes", "Microbiome carriage",
     "Individuals carry a microbiome compartment for each of the 42 bacteria. "
     "Colonisation is acquired from the environment and lost at organism-specific "
     "clearance rates. Carriage drives endogenous infection and acts as a reservoir "
     "for resistance gene amplification."),

    ("Disease processes", "Horizontal gene transfer (HGT)",
     "Inter-species plasmid transfer of resistance genes, modelled as a per-contact "
     "probability matrix across organism pairs. Transfer can occur during co-colonisation "
     "in the microbiome compartment. A 42 × 42 × mechanism HGT matrix governs "
     "species-pair permissibility."),

    ("Disease processes", "Mortality",
     "Three causes: (i) background age- and region-specific all-cause mortality; "
     "(ii) infection-attributable death (case fatality rates by syndrome, severity, "
     "organism, and treatment outcome); (iii) drug toxicity mortality for selected agents."),

    # ── Resistance biology ────────────────────────────────────────────────────
    ("Resistance biology", "De novo emergence",
     "Each infected individual has a daily probability of generating a resistant mutant, "
     "parameterised as organism- and mechanism-specific mutation rate × antibiotic "
     "selection pressure (measured by active treatment). Rates are calibrated effective "
     "hazards rather than literal point-mutation rates."),

    ("Resistance biology", "Fitness cost and reversion",
     "Resistance mechanisms carry a fitness cost. In the absence of antibiotic pressure, "
     "resistance reverts at mechanism-specific daily rates, scaled by a global reversion "
     "multiplier. This allows resistance to decline when drug use is reduced."),

    ("Resistance biology", "Resistance profile inheritance",
     "On infection, a resistance profile is sampled from a running cache of observed "
     "profiles for that organism (weighted by recency and current prevalence), ensuring "
     "that transmitted infections reflect the current circulating resistance landscape "
     "rather than drawing each mechanism independently."),

    # ── Calibration ─────────────────────────────────────────────────────────
    ("Calibration", "Calibration window",
     "2022–2025 (4-year window; ~1,461 simulated days). Summary statistics are averaged "
     "over this window to smooth stochastic noise before comparison with surveillance "
     "point estimates."),

    ("Calibration", "Calibration targets",
     "Antibiotic consumption (ECDC ESAC-Net, WHO AWaRe consumption data); infection "
     "resistance benchmarks informed by ECDC EARS-Net, WHO GLASS reports, and "
     "organism-specific literature; infection incidence and "
     "deaths (IHME Global Burden of Disease 2019); sepsis incidence (Rudd et al. 2020); "
     "carriage and hospital-acquired infection rates (published HAI epidemiology)."),

    ("Calibration", "Acceptance criteria",
     "There is an established set of criteria for a calibration being considered adequate — "
     "in initial work we generate 5 diverse sets of parameter values each leading to adequate "
     "calibration and evaluate policy comparisons in the context of each — in this way we begin "
     "to take account of parameter uncertainty in the answer to the policy question."),

    ("Calibration", "Uncertainty quantification",
     "Accepted runs differ in their random-number seeds, capturing stochastic variability. "
     "Results are reported as median (5th–95th percentile) across accepted runs."),

    # ── Counterfactual and burden estimation ────────────────────────────────
    ("Counterfactual", "Counterfactual design",
     "A resistance-free counterfactual is constructed by replaying each accepted run with "
     "resistance emergence disabled. Individuals still acquire infections and may die of "
     "them, but all bacteria remain fully susceptible. AMR-attributable deaths = observed "
     "deaths (with resistance) − counterfactual deaths (without resistance)."),

    ("Counterfactual", "AMR-attributable mortality",
     "Separable into directly attributable (deaths where resistance caused treatment "
     "failure) and associated (deaths in resistant-infection patients who may have died "
     "even with full susceptibility). Both components are reported."),

    # ── Software ──────────────────────────────────────────────────────────────
    ("Software and reproducibility", "Implementation language",
     "Rust (edition 2021, stable toolchain). Parallelised across simulation runs using "
     "Rayon. Analysis and paper tables produced in Python (pandas, matplotlib)."),

    ("Software and reproducibility", "Reproducibility",
     "Fixed-seed mode available for deterministic replication. All parameters are stored "
     "in a single configuration module (<code>src/config.rs</code>, ~11,700 lines). "
     "Code and configuration will be archived at acceptance."),
]


def make_t1(out_dir: Path) -> None:
    """
    Generate Table 1 — Model summary. This is a static 'hand-written' table;
    it draws no data from the calibration output files.
    """
    # Group rows by section
    sections: dict[str, list[tuple[str, str]]] = {}
    for section, feature, detail in _T1_ROWS:
        sections.setdefault(section, []).append((feature, detail))

    body  = _html_head("Table 1 — Model Summary")
    body += _back_link()
    body += "<h1>Table 1. Summary of simulation model design</h1>\n"
    body += (
        "<p class='note'>"
        "This table summarises the structure and scope of the individual-based "
        "antimicrobial resistance simulation. Full methodological detail is provided "
        "in the accompanying Technical Model Description."
        "</p>\n"
    )

    # Build a single merged table with section-header rows
    lines = ["<table>"]
    lines.append(
        "<thead><tr>"
        "<th style='width:14%'>Section</th>"
        "<th style='width:22%'>Feature</th>"
        "<th>Detail</th>"
        "</tr></thead><tbody>"
    )

    for section, rows in sections.items():
        n = len(rows)
        for i, (feature, detail) in enumerate(rows):
            tr_cls = " class='total-row'" if i == 0 else ""
            lines.append("<tr>")
            if i == 0:
                lines.append(
                    f"<td rowspan='{n}' style='font-weight:600; vertical-align:top; "
                    f"background:#eef1f5; white-space:normal;'>{section}</td>"
                )
            lines.append(
                f"<td style='font-weight:500; vertical-align:top; white-space:normal; "
                f"background:#f8f9fb;'>{feature}</td>"
            )
            lines.append(
                f"<td style='vertical-align:top; white-space:normal; text-align:left;'>{detail}</td>"
            )
            lines.append("</tr>")

    lines.append("</tbody></table>")
    body += "\n".join(lines) + "\n"

    body += _html_footnotes([
        "ESKAPE pathogens: <em>Enterococcus faecium</em>, <em>Staphylococcus aureus</em>, "
        "<em>Klebsiella pneumoniae</em>, <em>Acinetobacter baumannii</em>, "
        "<em>Pseudomonas aeruginosa</em>, and <em>Enterobacter</em> spp.",
        "IHME priority organisms: those listed in Murray CJ et al. (2022). Global burden "
        "of bacterial antimicrobial resistance in 2019: a systematic analysis. "
        "<em>Lancet</em> 399:629–655.",
        "ATC = Anatomical Therapeutic Chemical classification. AWaRe = Access, Watch, "
        "Reserve (WHO antibiotic categorisation). EARS-Net = European Antimicrobial "
        "Resistance Surveillance Network. GLASS = Global Antimicrobial Resistance and "
        "Use Surveillance System.",
    ])
    body += "</body></html>"
    _save(out_dir / TABLES_DIRNAME / "T1__model_summary.html", body)


# ---------------------------------------------------------------------------
# Table T2 — Headline Calibration Metrics + Block Scores
# ---------------------------------------------------------------------------

def _load_current_headline_targets() -> dict[str, float]:
    path = CALIBRATION_TARGETS_PATH
    if not path.exists():
        return {}

    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return {}

    key_to_metric = {
        "infection_deaths_millions": "Annual infection deaths (millions per year)",
        "people_on_antibiotics_millions": "People on antibiotics on an average day (millions)",
        "annual_infection_incidence_percent": "Incidence of bacterial infection per year (%)",
        "sepsis_incident_cases_millions": "Incident cases of sepsis per year (millions)",
    }

    targets: dict[str, float] = {}
    for metric in payload.get("headline_metrics", []):
        if not isinstance(metric, dict):
            continue
        key = metric.get("key")
        target = metric.get("target")
        if key in key_to_metric and isinstance(target, (int, float)):
            targets[key_to_metric[key]] = float(target)
    return targets


# Legacy table generators below are retained for reference only. They are not
# called by main(); the current paper output includes Table 1 and main Figures 1-13.
def make_t2(agg: dict, out_dir: Path) -> None:
    hm = agg.get("headline_metrics", pd.DataFrame()).copy()
    n  = agg.get("n_runs", 1)
    configured_targets = _load_current_headline_targets()

    if not hm.empty:
        import re

        def _rename_metric(raw: str) -> str:
            name = re.sub(r'\s*\(\d+\)\s*$', '', str(raw).strip())
            lo = name.lower()
            if 'infection deaths' in lo:
                return 'Annual infection deaths (millions per year)'
            if 'antibiotics' in lo:
                return 'People on antibiotics on an average day (millions)'
            if 'incidence' in lo and 'infection' in lo:
                return 'Incidence of bacterial infection per year (%)'
            if 'sepsis' in lo:
                return 'Incident cases of sepsis per year (millions)'
            return name

        def _ref_for_metric(renamed: str) -> str:
            lo = renamed.lower()
            if 'infection deaths' in lo:
                return 'Murray et al. 2022 (ref 1)'
            if 'antibiotics' in lo:
                return 'Klein et al. 2018 (ref 2)'
            if 'incidence' in lo and 'infection' in lo:
                return 'Vos et al. 2020 (ref 3)'
            if 'sepsis' in lo:
                return 'Rudd et al. 2020 (ref 4)'
            return '—'

        hm["Metric"] = hm["Metric"].apply(_rename_metric)
        if configured_targets:
            hm["Target"] = hm.apply(
                lambda row: configured_targets.get(str(row.get("Metric", "")), row.get("Target")),
                axis=1,
            )

        # Drop Delta and Unit columns; rename Target; add References column
        hm = hm.drop(columns=[c for c in hm.columns if c.lower().startswith("delta")], errors="ignore")
        hm = hm.drop(columns=["Unit"], errors="ignore")
        hm = hm.rename(columns={"Target": "Estimate from observed data"})

        hm["References"] = hm["Metric"].apply(_ref_for_metric)

    interval_note = (
        f"Median (5th–95th percentile) across {n} accepted calibration runs."
        if n > 1 else "Single accepted calibration run."
    )

    methods_note = (
        f"Simulation values are means over the 2022–2025 calibration window "
        f"(1,461 simulated days), scaled to the global population using the "
        f"run-specific population scale factor. {interval_note}"
    )

    footnotes = [
        "Infection deaths: the headline target is set to 6.4 million model-scope bacterial "
        "infection deaths per year, matching the simulation numerator by excluding H. pylori "
        "and MDR-TB. The broader GBD/33-pathogen central estimate is approximately 7.7 "
        "million bacterial-pathogen-associated deaths; AMR-specific deaths are a separate "
        "subset.",
        "Klein EY et al. (2018). Global increase and geographic convergence in antibiotic "
        "consumption between 2000 and 2015. <em>PNAS</em> 115:E3463–E3470. The antibiotic "
        "headline target is set to 100 million daily users, not 130 million, because Klein "
        "reports DDD-based consumption rather than unique people on treatment. The revised "
        "target is a person-prevalence proxy that sits below the DDD-equivalent total after "
        "allowing for dose intensity, stock-sales mismatch, and wastage.",
        "Vos T et al. (2020). Global burden of 369 diseases and injuries in 204 countries "
        "and territories, 1990–2019. <em>Lancet</em> 396:1204–1222.",
        "Rudd KE et al. (2020). Global, regional, and national sepsis incidence and mortality, "
        "1990–2017. <em>Lancet</em> 395:200–211. The sepsis headline target is set to 30 "
        "million, not 35 million, because Rudd reports all-cause sepsis whereas the model is "
        "bacteria-only; the revised target preserves a large bacterial burden without forcing "
        "the simulation up to an all-cause benchmark.",
    ]

    body  = _html_head("Table 2 — Comparison of simulation outputs with observed data")
    body += _back_link()
    body += "<h1>Table 2. Comparison of simulation outputs with estimates based on observed data</h1>\n"
    body += f"<p class='note'>{methods_note}</p>\n"

    if not hm.empty:
        body += _html_table(hm)

    body += _html_footnotes(footnotes)
    body += "</body></html>"
    _save(out_dir / "main" / "T2_headline_metrics.html", body)


# ---------------------------------------------------------------------------
# Table T3 — Drug Class Share
# ---------------------------------------------------------------------------

def make_t3(agg: dict, out_dir: Path) -> None:
    dc = agg.get("drug_class_share", pd.DataFrame()).copy()
    n  = agg.get("n_runs", 1)
    if dc is None or dc.empty:
        return
    dc = _clean_df(dc)

    footnotes = [
        _window_note(n),
        "Global drug-class share estimates are derived from the WHO AWaRe classification "
        "database, IQVIA MIDAS market data, and ECDC ESAC-Net surveillance, adjusted to "
        "represent the global population mix. Class-level global shares carry substantial "
        "uncertainty; estimates should be interpreted as approximate calibration anchors.",
        "Drug class codes follow the WHO ATC classification system (J01 antibacterials for "
        "systemic use). Classes without a standard J01 code (e.g. fidaxomicin, rifamycins) "
        "are given their closest available ATC code.",
    ]

    body  = _html_head("Table 3 — Drug Class Share")
    body += _back_link()
    body += "<h1>Table 3. Antibiotic use by drug class: simulation vs. global estimates, 2025</h1>\n"
    body += _html_table(dc)
    body += _html_footnotes(footnotes)
    body += "</body></html>"
    _save(out_dir / "main" / "T3_drug_class_share.html", body)


# ---------------------------------------------------------------------------
# Table T4 — Bacterial Infection Prevalence and Carriage
# ---------------------------------------------------------------------------

def make_t4(agg: dict, out_dir: Path) -> None:
    bi = agg.get("bacteria_infections", pd.DataFrame())
    n  = agg.get("n_runs", 1)
    if bi is None or bi.empty:
        return

    # Select the columns we want in the paper table
    display_cols = [
        "Bacteria",
        "Infection target (%)", "Infection simulation (%)",
        "Hospital Acquired (%)",
        "Carriage target (%)", "Carriage simulation (%)",
    ]
    present = [c for c in display_cols if c in bi.columns]
    table_df = bi[present].copy()

    interval_note = (
        f"Median (5th–95th percentile) across {n} accepted calibration runs."
        if n > 1 else "Single accepted calibration run."
    )

    footnotes = [
        f"Infection prevalence is the percentage of the global population with an active "
        f"bacterial infection from the specified organism on an average day during the "
        f"calibration window. {interval_note}",
        "Hospital-acquired (%) is the proportion of simulated new infections where acquisition "
        "occurred during hospitalisation. Patients who were admitted with a community-acquired "
        "infection also appear in hospital but are not counted here as hospital-acquired.",
        "Carriage target and simulation values represent the percentage of the global population "
        "colonised asymptomatically with the organism in the gut or upper respiratory microbiome. "
        "Targets are drawn from published cross-sectional carriage surveys; individual sources are "
        "given in the model description (Section 8).",
        "Organisms where the simulated infection rate exceeds twice or falls below half of the "
        "target value in the reference calibration run are flagged in the full calibration output "
        "(Supplementary Table S4).",
    ]

def make_t4(agg: dict, out_dir: Path) -> None:
    bi = agg.get("bacteria_infections", pd.DataFrame()).copy()
    n  = agg.get("n_runs", 1)
    if bi is None or bi.empty:
        return
    bi = _clean_df(bi)

    display_cols = [
        "Bacteria",
        "Infection observed estimate (%)", "Infection simulation (%)",
        "Carriage observed estimate (%)", "Carriage simulation (%)",
    ]
    present = [c for c in display_cols if c in bi.columns]
    table_df = bi[present] if present else bi

    # Summary row: sum each numeric column across all bacteria.
    # Because each prevalence figure is an independent % of the world population
    # (a person can carry multiple organisms simultaneously), the sum gives the
    # pooled cross-organism burden, not a probability — hence the label "pooled sum".
    numeric_cols = [c for c in table_df.columns if c != "Bacteria"]
    summary: dict = {"Bacteria": "All organisms (pooled sum)"}
    for col in numeric_cols:
        vals = pd.to_numeric(table_df[col], errors="coerce")
        summary[col] = round(float(vals.sum()), 2)
    summary_row = pd.DataFrame([summary], columns=table_df.columns)
    table_df = pd.concat([table_df, summary_row], ignore_index=True)

    footnotes = [
        _window_note(n),
        "Infection prevalence is the percentage of the global population with an active "
        "bacterial infection from the specified organism on an average day.",
        "Carriage estimates represent the percentage of the global population colonised "
        "asymptomatically in the gut or upper respiratory microbiome. Sources are given in "
        "the model description (Section 8).",
        "The bottom row pools across all 42 organisms by summing each column. Because a "
        "single person may carry or be infected by multiple organisms simultaneously, the "
        "pooled sum can exceed 100% and should be read as a total burden figure rather "
        "than a population prevalence.",
    ]

    body  = _html_head("Table 4 — Bacterial Infection and Carriage Prevalence")
    body += _back_link()
    body += (
        "<h1>Table 4. Bacterial infection prevalence and microbiome carriage "
        "(% world population), 2025</h1>\n"
    )
    body += _html_table(table_df, total_marker="all organisms")
    body += _html_footnotes(footnotes)
    body += "</body></html>"
    _save(out_dir / "main" / "T4_bacteria_burden.html", body)


# ---------------------------------------------------------------------------
# Table T5 — Hospital vs. Community Resistance by Organism
# ---------------------------------------------------------------------------

def make_t5(agg: dict, out_dir: Path) -> None:
    srl = agg.get("serious_resistance_locus", pd.DataFrame())
    ril = agg.get("resistance_incidence_locus", pd.DataFrame())
    bi  = agg.get("bacteria_infections", pd.DataFrame())
    n   = agg.get("n_runs", 1)

    if srl is None or srl.empty:
        return

    # Filter on the any-R H:C benchmark BEFORE _clean_df, because _clean_df renames every
    # column containing "target" → "observed estimate", which would break the lookup.
    # Select on the any-R structural benchmark; it is not a serious-R target.
    target_col = "Target H:C ratio"
    first_col = srl.columns[0]
    summary_mask = srl[first_col].astype(str).str.match(
        r"^\s*(-|Resistance Locus|Serious Resistance|Mean |H:C)", na=False
    )
    srl = srl[~summary_mask].copy()

    ril_raw = ril.copy() if ril is not None and not ril.empty else pd.DataFrame()
    if not ril_raw.empty:
        ril_first_col = ril_raw.columns[0]
        ril_summary_mask = ril_raw[ril_first_col].astype(str).str.match(
            r"^\s*(-|Resistance Locus|Serious Resistance|Mean |H:C)", na=False
        )
        ril_raw = ril_raw[~ril_summary_mask].copy()

    hc_col_present = target_col in ril_raw.columns and "Bacteria" in ril_raw.columns
    if hc_col_present:
        ril_raw[target_col] = pd.to_numeric(ril_raw[target_col], errors="coerce")
        included_names = set(ril_raw.loc[ril_raw[target_col] > 1.0, "Bacteria"].astype(str))
        included_raw = srl[srl["Bacteria"].astype(str).isin(included_names)].copy()
        excluded_raw = ril_raw[
            ril_raw[target_col].notna() & (ril_raw[target_col] <= 1.0)
        ].copy()
    else:
        included_raw = srl.copy()
        excluded_raw = pd.DataFrame()

    srl = _clean_df(included_raw)
    ril = _clean_df(ril_raw) if not ril_raw.empty else pd.DataFrame()
    bi  = _clean_df(bi.copy())  if bi  is not None and not bi.empty  else pd.DataFrame()

    # `excluded` only needed for the footnote — keep the Bacteria column from the raw slice.
    excluded = _clean_df(excluded_raw) if not excluded_raw.empty else pd.DataFrame()
    included = srl  # already filtered and cleaned

    # Build output table — start from serious-R columns
    keep_srl = ["Bacteria", "Hospital Serious-R (%)", "Community Serious-R (%)"]
    out = included[[c for c in keep_srl if c in included.columns]].copy()

    # Merge any-R columns from resistance_incidence_locus
    if not ril.empty and "Bacteria" in ril.columns:
        any_r_cols = ["Bacteria", "Hospital any-R (%)", "Community any-R (%)"]
        ril_sub = ril[[c for c in any_r_cols if c in ril.columns]].copy()
        out = out.merge(ril_sub, on="Bacteria", how="left")

    # Merge hospital-acquired % from bacteria_infections
    if not bi.empty and "Bacteria" in bi.columns:
        ha_col = "Hospital Acquired (%)"
        if ha_col in bi.columns:
            out = out.merge(bi[["Bacteria", ha_col]], on="Bacteria", how="left")

    # Reorder columns: Bacteria | HA% | H any-R | C any-R | H serious-R | C serious-R
    desired_order = [
        "Bacteria",
        "Hospital Acquired (%)",
        "Hospital any-R (%)",
        "Community any-R (%)",
        "Hospital Serious-R (%)",
        "Community Serious-R (%)",
    ]
    out = out[[c for c in desired_order if c in out.columns]]

    # Summary row: unweighted mean across all included organisms for each % column.
    # A plain sum would be meaningless here because columns are percentages on different
    # denominators (hospital-acquired vs community-acquired sub-populations).
    # An unweighted mean gives a single-number cross-organism comparison consistent with
    # how multi-organism resistance summaries are reported in surveillance literature.
    numeric_cols = [c for c in out.columns if c != "Bacteria"]
    summary: dict = {"Bacteria": "All included organisms (mean)"}
    for col in numeric_cols:
        vals = pd.to_numeric(out[col], errors="coerce")
        summary[col] = round(float(vals.mean()), 1)
    summary_row = pd.DataFrame([summary], columns=out.columns)
    out = pd.concat([out, summary_row], ignore_index=True)

    # Footnote: list organisms whose structural any-R benchmark is 1.0.
    if not excluded.empty and "Bacteria" in excluded.columns:
        excl_names = sorted(str(v) for v in excluded["Bacteria"].dropna().unique())
        excl_list = "; ".join(excl_names)
    else:
        excl_list = "none"

    footnotes = [
        _window_note(n),
        "This table includes organisms whose expert-assigned structural any-R benchmark indicates "
        "higher resistance among hospital-acquired than community-acquired cases. These benchmarks "
        "encode a qualitative expected setting gradient and are not direct harmonised empirical "
        "estimates or marker-drug serious-R targets. Operationally, inclusion requires an any-R "
        "hospital:community (H:C) benchmark greater than 1.0. "
        f"Organisms with an any-R H:C benchmark of 1.0 are excluded ({excl_list}).",
        "Hospital-acquired (%) is the simulated proportion of new infections acquired during "
        "hospitalisation.",
        "Hospital any-R (%) and Community any-R (%) are the percentages of new hospital- and "
        "community-acquired infections, respectively, carrying any resistance mechanism, "
        "averaged across all drugs with non-negligible potency for that organism.",
        "Hospital serious-R (%) and Community serious-R (%) use a single clinically important "
        "marker drug per organism (e.g. meropenem for Gram-negatives, flucloxacillin for "
        "<em>S. aureus</em>, vancomycin for enterococci) to give a focused hospital vs. "
        "community resistance comparison.",
        "The bottom row gives the unweighted mean across all included organisms. Because each "
        "column is a percentage on a different denominator (hospital-acquired vs. community-acquired "
        "infections for that specific organism), summation is not meaningful; the unweighted mean "
        "is consistent with multi-organism summaries in AMR surveillance literature.",
    ]

    body  = _html_head("Table 5 — Hospital vs. Community Resistance by Organism")
    body += _back_link()
    body += (
        "<h1>Table 5. Hospital- vs. community-acquired infection resistance "
        "by organism, 2025</h1>\n"
    )
    body += _html_table(out, total_marker="all included organisms")
    body += _html_footnotes(footnotes)
    body += "</body></html>"
    _save(out_dir / "main" / "T5_resistance_fit.html", body)


# ---------------------------------------------------------------------------
# Table T6 — AMR-Attributable Deaths (placeholder)
# ---------------------------------------------------------------------------

def make_t6_placeholder(out_dir: Path) -> None:
    body  = _html_head("Table 6 — AMR-Attributable Deaths (placeholder)")
    body += _back_link()
    body += "<h1>Table 6. AMR-attributable deaths, 2022–2025 (placeholder)</h1>\n"
    body += (
        "<p>This table will present AMR-attributable deaths derived from the counterfactual "
        "branch (resistance-free scenario from 2022). It requires a completed set of paired "
        "baseline and counterfactual simulation runs for each accepted parameter set.</p>"
        "<p>Planned columns:</p>"
        "<ul>"
        "<li>Organism</li>"
        "<li>Baseline deaths (millions/year): median (5th–95th percentile)</li>"
        "<li>Counterfactual deaths (millions/year): median (5th–95th percentile)</li>"
        "<li>AMR-attributable deaths (millions/year): median (5th–95th percentile)</li>"
        "<li>Attributable fraction (%)</li>"
        "</ul>"
        "<p>See Section&nbsp;11.1 of the model description for the counterfactual design. "
        "The branch point will be set to 2022 once calibration is finalised.</p>"
    )
    body += "</body></html>"
    _save(out_dir / "main" / "T6_amr_attributable_deaths_PLACEHOLDER.html", body)


# ---------------------------------------------------------------------------
# Supplementary S1 — Infection Deaths per Organism
# ---------------------------------------------------------------------------

def make_s1(agg: dict, out_dir: Path) -> None:
    bm = agg.get("bacteria_mortality", pd.DataFrame()).copy()
    n  = agg.get("n_runs", 1)
    if bm is None or bm.empty:
        return
    bm = _clean_df(bm)
    bm = bm.drop(columns=["Mortality Hospital Acquired (%)"], errors="ignore")

    footnotes = [
        _window_note(n),
        "Deaths are scaled to the global population using the run-specific population scale "
        "factor and annualised to yearly equivalents.",
        "Death estimates are based on GBD 2019/2020 cause-of-death attributions, WHO mortality "
        "data, and organism-specific published literature. Organism-level totals sum to more "
        "than the headline all-cause bacterial death estimate because polymicrobial deaths are "
        "attributed to all contributing pathogens.",
    ]

    body  = _html_head("S1 — Infection Deaths per Organism")
    body += _back_link()
    body += (
        "<h1>Supplementary Table S1. Infection deaths per organism: "
        "simulation vs. observed estimate (2025)</h1>\n"
    )
    body += _html_table(bm)
    body += _html_footnotes(footnotes)
    body += "</body></html>"
    _save(out_dir / "supplementary" / "S1_infection_deaths.html", body)


# ---------------------------------------------------------------------------
# Supplementary S2 — Syndrome Incidence
# ---------------------------------------------------------------------------

def make_s2(agg: dict, out_dir: Path) -> None:
    si = agg.get("syndrome_incidence", pd.DataFrame())
    n  = agg.get("n_runs", 1)

    footnotes = [
        _window_note(n),
        "Syndrome incidence is the simulated annual rate per 100,000 population. "
        "Separate syndrome-specific calibration targets have not been defined; "
        "the syndrome distribution is an emergent product of organism-specific "
        "infection rates and syndrome-assignment probabilities.",
    ]

    body  = _html_head("S2 — Syndrome Incidence")
    body += _back_link()
    body += "<h1>Supplementary Table S2. Syndrome incidence, 2025</h1>\n"

    if si is not None and not si.empty:
        body += _html_table(si, total_marker="TOTAL")

    body += _html_footnotes(footnotes)
    body += "</body></html>"
    _save(out_dir / "supplementary" / "S2_syndrome_incidence.html", body)


# ---------------------------------------------------------------------------
# Supplementary S3 — % of New Infections with Resistance by Organism and Hospital Status
# ---------------------------------------------------------------------------

def make_s3(agg: dict, out_dir: Path) -> None:
    ril = agg.get("resistance_incidence_locus", pd.DataFrame())
    n   = agg.get("n_runs", 1)

    if ril is not None and not ril.empty:
        ril = ril.drop(columns=["Total New Infections"], errors="ignore").copy()
        ril = ril.rename(columns={
            "Hospital Infections with Any Resistance (%)":  "Hospital-acquired new infections with any resistance (%)",
            "Community Infections with Any Resistance (%)": "Community-acquired new infections with any resistance (%)",
        })
        # Drop spurious summary-stat rows that the parser picks up after the
        # data table (e.g. "Resistance Locus Fit Summary", "- Bacteria with...",
        # "Serious Resistance Locus Summary"). These have non-bacteria text in
        # the first column and NaN in all numeric columns.
        first_col = ril.columns[0]
        summary_mask = ril[first_col].astype(str).str.match(
            r"^\s*(-|Resistance Locus|Serious Resistance|Mean |H:C)", na=False
        )
        ril = ril[~summary_mask].copy()

    footnotes = [
        _window_note(n),
        "Shows the percentage of new infections carrying any resistance mechanism, "
        "split by hospital-acquired and community-acquired origin. "
        "This is a descriptive output, not a calibration target.",
    ]

    body  = _html_head("S3 — Percentage of New Infections with Resistance by Organism and Hospital Status")
    body += _back_link()
    body += (
        "<h1>Supplementary Table S3. Percentage of new infections with resistance, "
        "by organism and hospital status, 2025</h1>\n"
    )

    if ril is not None and not ril.empty:
        body += _html_table(ril)

    body += _html_footnotes(footnotes)
    body += "</body></html>"
    _save(out_dir / "supplementary" / "S3_resistance_by_acquisition_route.html", body)


# ---------------------------------------------------------------------------
# Supplementary S4 — Full Resistance Benchmarks per Organism
# ---------------------------------------------------------------------------

def make_s4(agg: dict, out_dir: Path) -> None:
    rb_all = agg.get("resistance_benchmarks", pd.DataFrame())
    n      = agg.get("n_runs", 1)
    if rb_all is None or rb_all.empty:
        return

    # Exclude negligible-potency rows
    if "Flags" in rb_all.columns:
        rb_nonneg = rb_all[
            ~rb_all["Flags"].str.contains("negligible", case=False, na=False)
        ].copy()
    else:
        rb_nonneg = rb_all.copy()

    rb_nonneg = _clean_df(rb_nonneg, target_label="Calibration benchmark")
    rb_nonneg = rb_nonneg.rename(columns={
        "Inf sim (%)":               "Percent of infections with resistance — simulation (%)",
        "Inf calibration benchmark (%)": "Percent of infections with resistance — evidence-informed calibration benchmark (%)",
        "Avg sim (%)":               "Average resistance level among resistant infection-days — simulation (%)",
        "Avg calibration benchmark (%)": "Average resistance level among resistant infection-days — expert-assigned model benchmark (%)",
        "Micro sim (%)":             "Percent of people carrying the bacterium in whom a resistant strain is present (%)",
    })

    display_cols = [c for c in [
        "Drug", "Class",
        "Percent of infections with resistance — simulation (%)",
        "Percent of infections with resistance — evidence-informed calibration benchmark (%)",
        "Average resistance level among resistant infection-days — simulation (%)",
        "Average resistance level among resistant infection-days — expert-assigned model benchmark (%)",
        "Percent of people carrying the bacterium in whom a resistant strain is present (%)",
    ] if c in rb_nonneg.columns]

    footnotes = [
        _window_note(n),
        "Only organism–drug combinations where the drug has non-negligible potency "
        "(baseline potency >= 0.15) are shown.",
        "<em>Percent of infections with resistance</em>: percentage of active infections "
        "carrying any resistance to this drug at a point in time "
        "(simulation vs. evidence-informed calibration benchmark).",
        "<em>Average resistance level among resistant infection-days</em>: among infection-days "
        "where any resistance is present, the mean resistance level expressed as a percentage (0–100%). "
        "A value near 100% indicates that resistance, when present, is essentially complete; "
        "lower values indicate partial resistance. This is distinct from the prevalence column above, "
        "which measures the proportion of infection-days with any resistance. Its comparison value is "
        "an expert-assigned model benchmark, not a direct surveillance estimate.",
        "<em>Percent of people carrying the bacterium in whom a resistant strain is present</em>: "
        "percentage of the global population carrying a resistant strain of this organism "
        "in the gut or upper respiratory microbiome.",
    ]

    body  = _html_head("S4 — Full Resistance Benchmarks per Organism")
    body += _back_link()
    body += (
        "<h1>Supplementary Table S4. Full resistance benchmarks per organism and drug, "
        "2025</h1>\n"
    )
    body += "<p class='note'>Only drug–organism combinations with non-negligible potency are shown. "
    body += "Rows are grouped by organism.</p>\n"

    if "Bacteria" in rb_nonneg.columns:
        organisms = rb_nonneg["Bacteria"].unique()
        for org in organisms:
            subset = rb_nonneg[rb_nonneg["Bacteria"] == org][display_cols].copy()
            body += f"<h2>{org}</h2>\n"
            body += _html_table(subset)
    else:
        body += _html_table(rb_nonneg[display_cols])

    body += _html_footnotes(footnotes)
    body += "</body></html>"
    _save(out_dir / "supplementary" / "S4_resistance_benchmarks.html", body)


# ---------------------------------------------------------------------------
# Figure 6A/B/C — Resistance trends (1930–2025)
# ---------------------------------------------------------------------------

#: Calendar year corresponding to simulation time_in_years == 0.
_F1_SIM_EPOCH_YEAR: int = 1930

#: Colours for the trend figure.
_F1_TREND_COLOUR_MEAN   = "#1565C0"   # dark blue — mean line
_F1_TREND_COLOUR_CLOUD  = "#90CAF9"   # light blue — 90% CI band

def _discover_f1_simulation_csvs(input_paths: list[Union[str, Path]]) -> list[Path]:
    """
    Return simulation_summary CSVs matching the supplied calibration files.

    Calibration summary names are not always exactly calibration_summary_{seed}.txt;
    some accepted-run files carry prefixes such as calibration_summary_abc574337.txt.
    For F1 we need the numeric run id, so extract the trailing six digits and look
    for simulation_summary_{run_id}.csv in the standard output directory.
    """
    csv_dir = SIMULATION_OUTPUTS_DIR
    csv_paths: list[Path] = []
    seen: set[Path] = set()

    for input_path in input_paths:
        path = _resolve_project_path(input_path)

        candidates: list[Path] = []
        if path.name.startswith("simulation_summary_") and path.suffix.lower() == ".csv":
            candidates.append(path)

        seed_token = path.stem.split("_")[-1]
        candidates.append(csv_dir / f"simulation_summary_{seed_token}.csv")

        run_id_match = re.search(r"(\d{6})$", path.stem)
        if run_id_match:
            run_id = run_id_match.group(1)
            candidates.append(csv_dir / f"simulation_summary_{run_id}.csv")

        for candidate in candidates:
            if candidate.exists():
                resolved = candidate.resolve()
                if resolved not in seen:
                    seen.add(resolved)
                    csv_paths.append(candidate)
                break

    return csv_paths


def _population_scale_factor_from_calibration(path: Union[str, Path]) -> float | None:
    try:
        text = _resolve_project_path(path).read_text(encoding="utf-8", errors="replace")
    except OSError:
        return None
    match = re.search(
        r"Population scale factor relative to calibration targets:\s*([0-9,.\-+eE]+)",
        text,
    )
    if not match:
        return None
    try:
        value = float(match.group(1).replace(",", ""))
    except ValueError:
        return None
    return value if np.isfinite(value) and value > 0.0 else None


def _discover_simulation_csvs_with_scale(
    input_paths: list[Union[str, Path]],
) -> list[tuple[Path, float | None]]:
    """
    Return matching simulation_summary CSVs paired with their calibration scale factor.
    """
    csv_dir = SIMULATION_OUTPUTS_DIR
    rows: list[tuple[Path, float | None]] = []
    seen: set[Path] = set()

    for input_path in input_paths:
        path = _resolve_project_path(input_path)
        scale_factor = (
            None
            if path.name.startswith("simulation_summary_") and path.suffix.lower() == ".csv"
            else _population_scale_factor_from_calibration(path)
        )

        candidates: list[Path] = []
        if path.name.startswith("simulation_summary_") and path.suffix.lower() == ".csv":
            candidates.append(path)

        seed_token = path.stem.split("_")[-1]
        candidates.append(csv_dir / f"simulation_summary_{seed_token}.csv")

        run_id_match = re.search(r"(\d{6})$", path.stem)
        if run_id_match:
            run_id = run_id_match.group(1)
            candidates.append(csv_dir / f"simulation_summary_{run_id}.csv")

        for candidate in candidates:
            if candidate.exists():
                resolved = candidate.resolve()
                if resolved not in seen:
                    seen.add(resolved)
                    rows.append((candidate, scale_factor))
                break

    return rows


_CSV_HEADER_CACHE: dict[Path, tuple[str, ...] | None] = {}


def _simulation_csv_column_names(csv_path: Path) -> tuple[str, ...] | None:
    resolved = csv_path.resolve()
    if resolved not in _CSV_HEADER_CACHE:
        try:
            _CSV_HEADER_CACHE[resolved] = tuple(pd.read_csv(csv_path, nrows=0).columns)
        except (FileNotFoundError, pd.errors.EmptyDataError, OSError, ValueError):
            _CSV_HEADER_CACHE[resolved] = None
    return _CSV_HEADER_CACHE[resolved]


def _simulation_csv_columns(csv_path: Path) -> set[str] | None:
    names = _simulation_csv_column_names(csv_path)
    return set(names) if names is not None else None


def _filter_simulation_csvs_with_columns(
    csv_paths: list[Path],
    required_columns: list[str],
    label: str,
) -> list[Path]:
    filtered: list[Path] = []
    for csv_path in csv_paths:
        columns = _simulation_csv_columns(csv_path)
        if columns is not None and all(column in columns for column in required_columns):
            filtered.append(csv_path)
    if csv_paths:
        print(
            f"  {label}: {len(filtered)} of {len(csv_paths)} matching simulation CSV(s) "
            "contain the required aggregate columns."
        )
    return filtered


def _read_csv_selected(csv_path: Path, usecols: list[str] | set[str]) -> pd.DataFrame:
    selected = list(dict.fromkeys(usecols))
    available = _simulation_csv_columns(csv_path)
    if available is not None:
        selected = [column for column in selected if column in available]
    try:
        return pd.read_csv(csv_path, usecols=selected, engine="pyarrow")
    except Exception:
        return pd.read_csv(csv_path, usecols=selected)


_F6_SERIOUS_RESISTANCE_COLUMN = "newly_infected_with_serious_resistance_count"
_F6_MARKER_ELIGIBLE_COLUMN = "newly_infected_serious_resistance_marker_eligible_count"
_F6_SERIOUS_R_MISSING_NOTE = (
    "Serious-R trend requires simulation_summary column "
    "newly_infected_with_serious_resistance_count. Re-run the Rust simulation after "
    "adding this field to show the serious-R line."
)
_F6B_DENOMINATOR_COLUMN = "new_active_infections_by_bacteria"
_F6_TOP_SELECTION_YEAR = 2025
_F6B_TOP_N_BACTERIA = 15
_F6_HOSPITAL_REGIONS = [
    "north_america",
    "south_america",
    "africa",
    "asia",
    "europe",
    "oceania",
]
_F6C_SERIOUS_R_VECTOR_COLUMNS = (
    "newly_infected_serious_r_by_bacteria",
    "newly_infected_with_serious_resistance_by_bacteria",
    "newly_infected_with_serious_resistance_count_by_bacteria",
)
_F6C_SERIOUS_R_PER_BACTERIUM_SUFFIXES = (
    "_newly_infected_serious_r",
    "_newly_infected_with_serious_resistance",
    "_newly_infected_serious_resistance",
    "_newly_infected_serious_resistance_count",
)


def _f6_hospital_columns_for_slug(slug: str, available: set[str]) -> list[str]:
    return [
        f"{slug}_newly_infected_hospital_{region}"
        for region in _F6_HOSPITAL_REGIONS
        if f"{slug}_newly_infected_hospital_{region}" in available
    ]


def _load_resistance_series(csv_path: Path) -> tuple[pd.DataFrame | None, bool, bool]:
    """
    Load a simulation_summary CSV and return annual-resistance inputs.

    The serious-R and marker-eligible columns are optional so older CSVs still
    generate the original any-resistance trend.
    """
    needed = ["time_in_years", "newly_infected_count", "newly_infected_with_resistance_count"]
    available_columns = _simulation_csv_columns(csv_path)
    if available_columns is None:
        return None, False, False

    has_serious_column = _F6_SERIOUS_RESISTANCE_COLUMN in available_columns
    has_marker_eligible_column = _F6_MARKER_ELIGIBLE_COLUMN in available_columns
    if any(column not in available_columns for column in needed):
        return None, has_serious_column, has_marker_eligible_column

    hospital_denominator_columns: list[str] = []
    hospital_any_r_columns: list[str] = []
    community_any_r_columns: list[str] = []
    for slug in _F15_KNOWN_BACTERIA_SLUGS:
        hospital_denominator_columns.extend(_f6_hospital_columns_for_slug(slug, available_columns))
        hospital_any_r_col = f"{slug}_newly_infected_any_r_hospital"
        community_any_r_col = f"{slug}_newly_infected_any_r_community"
        if hospital_any_r_col in available_columns:
            hospital_any_r_columns.append(hospital_any_r_col)
        if community_any_r_col in available_columns:
            community_any_r_columns.append(community_any_r_col)
    has_setting_columns = (
        _F6B_DENOMINATOR_COLUMN in available_columns
        and bool(hospital_denominator_columns)
        and bool(hospital_any_r_columns)
        and bool(community_any_r_columns)
    )

    usecols = list(needed)
    if has_serious_column:
        usecols.append(_F6_SERIOUS_RESISTANCE_COLUMN)
    if has_setting_columns:
        usecols.extend([
            _F6B_DENOMINATOR_COLUMN,
            *hospital_denominator_columns,
            *hospital_any_r_columns,
            *community_any_r_columns,
        ])
    try:
        df = _read_csv_selected(csv_path, usecols)
    except (FileNotFoundError, ValueError):
        return None, has_serious_column, has_marker_eligible_column

    df = df.dropna(subset=needed)
    df = df[df["newly_infected_count"] > 0].copy()
    if df.empty:
        return None, has_serious_column, has_marker_eligible_column
    df["year"] = _F1_SIM_EPOCH_YEAR + df["time_in_years"]
    df["pct_resistant"] = (
        df["newly_infected_with_resistance_count"] / df["newly_infected_count"] * 100.0
    )
    columns = ["year", "pct_resistant"]
    if has_setting_columns:
        numeric_hospital_denom = pd.DataFrame({
            col: pd.to_numeric(df[col], errors="coerce").fillna(0.0)
            for col in hospital_denominator_columns
        })
        numeric_hospital_num = pd.DataFrame({
            col: pd.to_numeric(df[col], errors="coerce").fillna(0.0)
            for col in hospital_any_r_columns
        })
        numeric_community_num = pd.DataFrame({
            col: pd.to_numeric(df[col], errors="coerce").fillna(0.0)
            for col in community_any_r_columns
        })
        target_len = len(_F15_KNOWN_BACTERIA_SLUGS)
        total_denominator = np.array([
            np.nansum(
                _figure_15_extend_array(
                    np.array(_figure_15_parse_vector_cell(value), dtype=float),
                    target_len,
                )
            )
            for value in df[_F6B_DENOMINATOR_COLUMN]
        ], dtype=float)
        hospital_denominator = numeric_hospital_denom.sum(axis=1).to_numpy(dtype=float)
        community_denominator = np.clip(total_denominator - hospital_denominator, 0.0, None)
        hospital_numerator = numeric_hospital_num.sum(axis=1).to_numpy(dtype=float)
        community_numerator = numeric_community_num.sum(axis=1).to_numpy(dtype=float)
        hospital_pct = np.full_like(hospital_denominator, np.nan, dtype=float)
        community_pct = np.full_like(community_denominator, np.nan, dtype=float)
        np.divide(
            hospital_numerator,
            hospital_denominator,
            out=hospital_pct,
            where=hospital_denominator > 0.0,
        )
        np.divide(
            community_numerator,
            community_denominator,
            out=community_pct,
            where=community_denominator > 0.0,
        )
        df["pct_any_resistant_hospital"] = hospital_pct * 100.0
        df["pct_any_resistant_community"] = community_pct * 100.0
        columns.extend(["pct_any_resistant_hospital", "pct_any_resistant_community"])
    if has_serious_column:
        df["pct_serious_resistant"] = (
            df[_F6_SERIOUS_RESISTANCE_COLUMN] / df["newly_infected_count"] * 100.0
        )
        columns.append("pct_serious_resistant")
    df["year_int"] = df["year"].apply(int)
    annual = df.groupby("year_int", as_index=False)[columns[1:]].mean()
    annual["year"] = annual["year_int"].astype(float)
    return annual[columns], has_serious_column, has_marker_eligible_column


def _combine_annual_resistance_series(
    series_list: list[pd.DataFrame],
    value_column: str,
) -> pd.DataFrame:
    annual_frames: list[pd.Series] = []
    for series in series_list:
        if value_column not in series.columns:
            continue
        one_run = series.dropna(subset=[value_column]).copy()
        if one_run.empty:
            continue
        one_run["year_int"] = one_run["year"].apply(int)
        annual_frames.append(one_run.groupby("year_int")[value_column].mean())
    if not annual_frames:
        return pd.DataFrame()
    return pd.concat(annual_frames, axis=1)


def make_figure_6_resistance_trend(csv_paths: list[Path], out_dir: Path) -> None:
    """
    Figure 6A: time trend of resistance among new active bacterial infections.

    Each run in *csv_paths* contributes one time series.  The figure shows:
      - solid mean line across all runs
      - shaded 90% credible-interval cloud (5th–95th percentile across runs)

    If *csv_paths* is empty, or the available data cover fewer than 2 calendar
    years, a clearly labelled placeholder panel is saved instead.
    """
    fig_dir = out_dir / FIGURES_DIRNAME
    fig_dir.mkdir(parents=True, exist_ok=True)
    stem = "Figure_6A__resistance_trends"
    png_path = fig_dir / f"{stem}.png"
    svg_path = fig_dir / f"{stem}.svg"
    html_path = fig_dir / f"{stem}.html"

    # ------------------------------------------------------------------ #
    # Load all run series                                                   #
    # ------------------------------------------------------------------ #
    series_list: list[pd.DataFrame] = []
    serious_series_list: list[pd.DataFrame] = []
    setting_series_list: list[pd.DataFrame] = []
    marker_eligible_column_available = False
    for p in csv_paths:
        s, has_serious_column, has_marker_eligible_column = _load_resistance_series(p)
        marker_eligible_column_available = (
            marker_eligible_column_available or has_marker_eligible_column
        )
        if s is not None and not s.empty:
            series_list.append(s)
            if has_serious_column and "pct_serious_resistant" in s.columns:
                serious_series_list.append(s)
            if {
                "pct_any_resistant_hospital",
                "pct_any_resistant_community",
            }.issubset(s.columns):
                setting_series_list.append(s)

    n_runs = len(series_list)
    n_serious_runs = len(serious_series_list)
    n_setting_runs = len(setting_series_list)
    serious_available = n_serious_runs > 0
    serious_column_missing = n_serious_runs < n_runs
    setting_available = n_setting_runs > 0
    setting_column_missing = n_setting_runs < n_runs
    data_available = n_runs > 0

    # Determine calendar year range covered by the data
    if data_available:
        year_min_data = min(s["year"].min() for s in series_list)
        year_max_data = max(s["year"].max() for s in series_list)
        data_span_years = year_max_data - year_min_data
        if data_span_years < 1.0:
            data_available = False   # too short to be meaningful

    # ------------------------------------------------------------------ #
    # Build the figure                                                      #
    # ------------------------------------------------------------------ #
    fig, ax = plt.subplots(figsize=(10, 5))

    if data_available:
        combined = _combine_annual_resistance_series(series_list, "pct_resistant")
        years = combined.index.values
        mean_vals = combined.mean(axis=1).values
        p5_vals = combined.quantile(0.05, axis=1).values
        p95_vals = combined.quantile(0.95, axis=1).values

        ax.fill_between(
            years,
            p5_vals,
            p95_vals,
            color=_F1_TREND_COLOUR_CLOUD,
            alpha=0.45,
        )
        ax.plot(
            years,
            mean_vals,
            color=_F1_TREND_COLOUR_MEAN,
            linewidth=1.8,
            label="Any resistance, all new infections",
        )

        if setting_available:
            for value_column, label, colour in [
                ("pct_any_resistant_hospital", "Any resistance, hospital-acquired", "#7B3294"),
                ("pct_any_resistant_community", "Any resistance, community-acquired", "#008837"),
            ]:
                setting_combined = _combine_annual_resistance_series(
                    setting_series_list,
                    value_column,
                )
                if setting_combined.empty:
                    continue
                setting_years = setting_combined.index.values
                setting_mean = setting_combined.mean(axis=1).values
                setting_p5 = setting_combined.quantile(0.05, axis=1).values
                setting_p95 = setting_combined.quantile(0.95, axis=1).values
                ax.fill_between(
                    setting_years,
                    setting_p5,
                    setting_p95,
                    color=colour,
                    alpha=0.12,
                )
                ax.plot(
                    setting_years,
                    setting_mean,
                    color=colour,
                    linewidth=1.55,
                    label=label,
                )

        if serious_available:
            serious_combined = _combine_annual_resistance_series(
                serious_series_list,
                "pct_serious_resistant",
            )
            if not serious_combined.empty:
                serious_years = serious_combined.index.values
                serious_mean = serious_combined.mean(axis=1).values
                serious_p5 = serious_combined.quantile(0.05, axis=1).values
                serious_p95 = serious_combined.quantile(0.95, axis=1).values
                ax.fill_between(
                    serious_years,
                    serious_p5,
                    serious_p95,
                    color="#F4A261",
                    alpha=0.28,
                )
                ax.plot(
                    serious_years,
                    serious_mean,
                    color="#C65D00",
                    linewidth=1.8,
                    label="Serious-R marker resistance",
                )

        # Annotate if data are sparse (calibration window only)
        if year_min_data > _F1_SIM_EPOCH_YEAR + 80:
            ax.text(
                0.02, 0.96,
                "Note: data currently cover calibration window only "
                f"({int(year_min_data)}\u2013{int(year_max_data)}).\n"
                "Full 1930\u20132025 trend will be available after a full-period simulation run.",
                transform=ax.transAxes, fontsize=8, va="top",
                color="#c0392b",
                bbox=dict(boxstyle="round,pad=0.3", fc="#fff3cd", ec="#f0c040", alpha=0.9),
            )
    else:
        # Placeholder panel
        ax.text(0.5, 0.5,
                "Figure 6A. Resistance trends\n\n"
                "Data not yet available.\n"
                "Re-run with full-period simulation output\n"
                "(simulation_summary_*.csv covering 1930-2025)\n"
                "to generate this figure.",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=11, color="#555",
                bbox=dict(boxstyle="round,pad=0.6", fc="#f5f5f5", ec="#bbb"))
        ax.set_axis_off()

    if data_available:
        ax.set_xlim(_F1_SIM_EPOCH_YEAR, _F1_SIM_EPOCH_YEAR + 96)
        ax.set_ylim(0, 100)
        ax.set_xlabel("Year", fontsize=11)
        ax.set_ylabel("New infections with resistance (%)", fontsize=11)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", linewidth=0.4, alpha=0.5)
        if serious_available and n_serious_runs != n_runs:
            n_label = f"Any n = {n_runs}; Serious-R n = {n_serious_runs}"
        elif setting_available and n_setting_runs != n_runs:
            n_label = f"Any n = {n_runs}; setting n = {n_setting_runs}"
        else:
            n_label = f"n = {n_runs} run{'s' if n_runs > 1 else ''}"
        ax.legend(fontsize=9, frameon=False, title=n_label, title_fontsize=8)

    fig.suptitle("Figure 6A. Resistance trends", fontsize=11, fontweight="bold")
    fig.tight_layout()

    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {png_path}")
    print(f"  Saved: {svg_path}")

    # HTML wrapper
    html_rel_img = f"{stem}.png"
    html_rel_svg = f"{stem}.svg"
    body  = _html_head("Figure 6A. Resistance Trends")
    body += _back_link()
    if data_available:
        figure_note = (
            f"Resistance among new active bacterial infections, plotted annually "
            f"from {int(year_min_data)} to {int(year_max_data)}. "
            f"Solid lines: means across accepted simulation run"
            f"{'s' if n_runs > 1 else ''}. "
            f"Shaded bands: 5th-95th percentile (90% interval) across runs."
        )
        footnotes = [
            figure_note,
            "Any resistance is the percentage of new active bacterial infections whose "
            "infection-level any-R value is positive for at least one drug. Serious-R marker "
            "resistance is the percentage whose bacterium-specific serious-R marker drug has "
            "infection-level any-R above the infection threshold at infection establishment. "
            "Both metrics use the same denominator: all new active bacterial infections.",
            "Values are averaged within each calendar year across all daily time steps.",
            "Full 1930-2025 data require a simulation run spanning the complete historical period, "
            "not just the 2022-2025 calibration window.",
        ]
        if serious_column_missing:
            footnotes.insert(1, _F6_SERIOUS_R_MISSING_NOTE)
        if setting_available:
            footnotes.append(
                "Hospital-acquired and community-acquired any-R lines use the same event "
                "definition as the overall any-R line. Hospital denominators are summed from "
                "per-bacterium newly_infected_hospital_<region> columns; community denominators "
                f"are {_F6B_DENOMINATOR_COLUMN} minus hospital-acquired new infections."
            )
        if setting_column_missing:
            footnotes.append(
                "Hospital/community any-R trend lines are shown only for simulation summary "
                "files containing the required per-bacterium hospital/community any-R columns."
            )
        if marker_eligible_column_available:
            footnotes.append(
                "The simulation summary also includes "
                "newly_infected_serious_resistance_marker_eligible_count as a QC field for "
                "new infections whose bacterium has a serious-R marker; the plotted denominator "
                "remains all new active bacterial infections."
            )
    else:
        figure_note = (
            "Placeholder - full-period simulation data (1930-2025) not yet available. "
            "Run the simulation in non-calibration mode and supply the resulting "
            "<code>simulation_summary_*.csv</code> files to generate this figure."
        )
        footnotes = [figure_note]
    body += (
        f"<img src='{html_rel_img}' alt='Figure 6A' "
        f"style='max-width:100%; border:1px solid #ddd; border-radius:4px;'>\n"
    )
    body += f"<p class='note'>Download: <a href='{html_rel_img}'>PNG</a> | <a href='{html_rel_svg}'>SVG</a></p>\n"
    body += _html_footnotes(footnotes)
    body += "</body></html>"
    _save(html_path, body)


def _f6b_any_r_columns_for_slug(slug: str, available: set[str]) -> list[str]:
    return [
        column
        for column in (
            f"{slug}_newly_infected_any_r_hospital",
            f"{slug}_newly_infected_any_r_community",
        )
        if column in available
    ]


def _load_bacteria_resistance_trend_rows(
    csv_path: Path,
    numerator_columns_for_slug: Callable[[str, set[str]], list[str]],
    numerator_field: str,
    metric_label: str,
) -> tuple[list[dict[str, object]], str | None]:
    columns = _simulation_csv_columns(csv_path)
    if columns is None:
        return [], f"{csv_path.name}: unable to read simulation summary header."
    if _F6B_DENOMINATOR_COLUMN not in columns:
        return [], f"{csv_path.name}: missing {_F6B_DENOMINATOR_COLUMN}."

    optional = ["policy_option", "run_id", "simulation_year", "year", "time_in_years", "time_step"]
    numerator_columns: list[str] = []
    per_bacterium_columns: dict[str, list[str]] = {}
    for slug in _F15_KNOWN_BACTERIA_SLUGS:
        slug_columns = numerator_columns_for_slug(slug, columns)
        per_bacterium_columns[slug] = slug_columns
        numerator_columns.extend(slug_columns)
    if not numerator_columns:
        return [], f"{csv_path.name}: missing per-bacterium {metric_label} newly infected columns."

    usecols = [_F6B_DENOMINATOR_COLUMN, *optional, *numerator_columns]
    try:
        df = _read_csv_selected(csv_path, usecols)
    except (FileNotFoundError, ValueError, OSError) as exc:
        return [], f"{csv_path.name}: unable to load Figure 6B columns ({exc})."
    if df.empty:
        return [], f"{csv_path.name}: no rows available."

    years = _simulation_year_series(df)
    valid_year = years.notna()
    df = df.loc[valid_year].copy()
    years = years.loc[valid_year]
    if df.empty:
        return [], f"{csv_path.name}: no rows with valid simulation year."

    target_len = len(_F15_KNOWN_BACTERIA_SLUGS)
    denominator_values = [
        _figure_15_extend_array(
            np.array(_figure_15_parse_vector_cell(value), dtype=float),
            target_len,
        )
        for value in df[_F6B_DENOMINATOR_COLUMN]
    ]
    if not denominator_values:
        return [], f"{csv_path.name}: no denominator values found."
    denominator_matrix = np.vstack(denominator_values)

    run_series = (
        df["run_id"].astype(str)
        if "run_id" in df.columns
        else pd.Series(csv_path.stem, index=df.index)
    )
    work = pd.DataFrame({
        "year_int": years.astype(int).to_numpy(),
        "run_id": run_series.to_numpy(),
    })
    records: list[dict[str, object]] = []
    for b_idx, slug in enumerate(_F15_KNOWN_BACTERIA_SLUGS):
        numerator_cols = per_bacterium_columns.get(slug, [])
        if not numerator_cols:
            continue
        numeric_num = pd.DataFrame({
            col: pd.to_numeric(df[col], errors="coerce").fillna(0.0)
            for col in numerator_cols
        })
        bacterium_work = work.copy()
        bacterium_work["denominator"] = denominator_matrix[:, b_idx]
        bacterium_work["numerator"] = numeric_num.sum(axis=1).to_numpy(dtype=float)
        annual = (
            bacterium_work
            .groupby(["run_id", "year_int"], as_index=False)[["denominator", "numerator"]]
            .sum()
        )
        annual = annual[annual["denominator"] > 0.0].copy()
        for _, row in annual.iterrows():
            records.append({
                "source": csv_path.name,
                "run_id": row["run_id"],
                "year": int(row["year_int"]),
                "bacterium_slug": slug,
                "bacterium": _figure_15_bacterium_label(slug),
                "new_active_infections": float(row["denominator"]),
                numerator_field: float(row["numerator"]),
            })

    if not records:
        return [], f"{csv_path.name}: no positive per-bacterium new-infection denominators."
    return records, None


def _f6b_placeholder(
    out_dir: Path,
    message: str,
    problems: list[str] | None = None,
) -> None:
    fig_dir = out_dir / FIGURES_DIRNAME
    fig_dir.mkdir(parents=True, exist_ok=True)
    stem = "Figure_6B__resistance_trends_by_bacterium"
    fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.text(
        0.5,
        0.5,
        f"Figure 6B. Resistance trends by bacterium\n\n{message}",
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=10.5,
        color="#555",
        bbox=dict(boxstyle="round,pad=0.6", fc="#f5f5f5", ec="#bbb"),
        wrap=True,
    )
    ax.set_axis_off()
    fig.subplots_adjust(left=0.03, right=0.97, top=0.90, bottom=0.08)
    png_path = fig_dir / f"{stem}.png"
    svg_path = fig_dir / f"{stem}.svg"
    html_path = fig_dir / f"{stem}.html"
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)

    body = _html_head("Figure 6B. Resistance Trends by Bacterium")
    body += _back_link()
    body += (
        f"<img src='{stem}.png' alt='Figure 6B' "
        f"style='max-width:100%; border:1px solid #ddd; border-radius:4px;'>\n"
    )
    body += f"<p class='note'>Download: <a href='{stem}.png'>PNG</a> | <a href='{stem}.svg'>SVG</a></p>\n"
    footnotes = [message]
    if problems:
        footnotes.append("Parser notes: " + " ".join(problems[:5]))
    body += _html_footnotes(footnotes)
    body += "</body></html>"
    _save(html_path, body)
    print(f"  Saved: {png_path}")
    print(f"  Saved: {svg_path}")


def make_figure_6b_resistance_trend_by_bacterium(
    csv_paths: list[Path],
    out_dir: Path,
) -> None:
    """Figure 6B: any-resistance trend for the highest-resistance bacteria in 2025."""
    rows: list[dict[str, object]] = []
    problems: list[str] = []
    for csv_path in csv_paths:
        csv_rows, problem = _load_bacteria_resistance_trend_rows(
            csv_path,
            _f6b_any_r_columns_for_slug,
            "new_active_infections_any_r",
            "any-R",
        )
        rows.extend(csv_rows)
        if problem:
            problems.append(problem)

    if not rows:
        _f6b_placeholder(
            out_dir,
            "Per-bacterium resistance trend data are not available. Required columns are "
            f"{_F6B_DENOMINATOR_COLUMN} and per-bacterium "
            "newly_infected_any_r_hospital/community columns.",
            problems,
        )
        return

    data = pd.DataFrame(rows)
    data["pct_any_resistant"] = np.where(
        data["new_active_infections"] > 0.0,
        data["new_active_infections_any_r"] / data["new_active_infections"] * 100.0,
        np.nan,
    )
    selection_year_data = data[data["year"] == _F6_TOP_SELECTION_YEAR].copy()
    if selection_year_data.empty:
        _f6b_placeholder(
            out_dir,
            f"No per-bacterium any-R rows were found for {_F6_TOP_SELECTION_YEAR}.",
            problems,
        )
        return
    top = (
        selection_year_data
        .groupby(["bacterium_slug", "bacterium"], as_index=False)[
            ["new_active_infections", "new_active_infections_any_r"]
        ]
        .sum()
    )
    top["selection_any_r_pct"] = np.where(
        top["new_active_infections"] > 0.0,
        top["new_active_infections_any_r"] / top["new_active_infections"] * 100.0,
        np.nan,
    )
    top = (
        top.dropna(subset=["selection_any_r_pct"])
        .sort_values(
            ["selection_any_r_pct", "new_active_infections"],
            ascending=[False, False],
        )
        .head(_F6B_TOP_N_BACTERIA)
    )
    top_slugs = top["bacterium_slug"].tolist()
    if not top_slugs:
        _f6b_placeholder(
            out_dir,
            f"No positive per-bacterium any-R denominators were found for {_F6_TOP_SELECTION_YEAR}.",
            problems,
        )
        return

    plot_data = data[data["bacterium_slug"].isin(top_slugs)].copy()
    annual = (
        plot_data
        .groupby(["bacterium_slug", "bacterium", "year"], as_index=False)[
            ["new_active_infections", "new_active_infections_any_r"]
        ]
        .sum()
    )
    annual["pct_any_resistant"] = np.where(
        annual["new_active_infections"] > 0.0,
        annual["new_active_infections_any_r"] / annual["new_active_infections"] * 100.0,
        np.nan,
    )
    year_min_data = int(annual["year"].min())
    year_max_data = int(annual["year"].max())

    fig_dir = out_dir / FIGURES_DIRNAME
    fig_dir.mkdir(parents=True, exist_ok=True)
    stem = "Figure_6B__resistance_trends_by_bacterium"
    png_path = fig_dir / f"{stem}.png"
    svg_path = fig_dir / f"{stem}.svg"
    html_path = fig_dir / f"{stem}.html"

    fig, ax = plt.subplots(figsize=(11.5, 6.2))
    colour_cycle = plt.cm.tab20(np.linspace(0, 1, len(top_slugs)))
    rank_lookup = {slug: idx + 1 for idx, slug in enumerate(top_slugs)}
    for colour, slug in zip(colour_cycle, top_slugs):
        bacterium_rows = annual[annual["bacterium_slug"] == slug].sort_values("year")
        if bacterium_rows.empty:
            continue
        label = f"{rank_lookup[slug]}. {bacterium_rows['bacterium'].iloc[0]}"
        ax.plot(
            bacterium_rows["year"].to_numpy(dtype=float),
            bacterium_rows["pct_any_resistant"].to_numpy(dtype=float),
            linewidth=1.7,
            marker="o",
            markersize=2.4,
            color=colour,
            label=label,
        )

    if year_min_data > _F1_SIM_EPOCH_YEAR + 80:
        ax.text(
            0.02,
            0.97,
            "Note: data currently cover calibration window only "
            f"({year_min_data}\u2013{year_max_data}).",
            transform=ax.transAxes,
            fontsize=8,
            va="top",
            color="#c0392b",
            bbox=dict(boxstyle="round,pad=0.3", fc="#fff3cd", ec="#f0c040", alpha=0.9),
        )

    ax.set_xlim(_F1_SIM_EPOCH_YEAR, _F1_SIM_EPOCH_YEAR + 96)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Year", fontsize=11)
    ax.set_ylabel("New infections with any resistance (%)", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", linewidth=0.4, alpha=0.5)
    ax.legend(
        fontsize=7.5,
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        title=f"Top {_F6B_TOP_N_BACTERIA} by {_F6_TOP_SELECTION_YEAR} any-R",
        title_fontsize=8,
    )
    fig.suptitle(
        "Figure 6B. Resistance trends by bacterium",
        fontsize=11,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 0.80, 1])
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {png_path}")
    print(f"  Saved: {svg_path}")

    top_table = top.copy()
    top_table["Rank"] = np.arange(1, len(top_table) + 1)
    top_table[f"{_F6_TOP_SELECTION_YEAR} any-R (%)"] = top_table["selection_any_r_pct"].map(
        lambda value: f"{float(value):.1f}"
    )
    top_table[f"{_F6_TOP_SELECTION_YEAR} new active infections"] = top_table["new_active_infections"].map(
        lambda value: f"{float(value):,.0f}"
    )
    top_table = top_table[[
        "Rank",
        "bacterium",
        f"{_F6_TOP_SELECTION_YEAR} any-R (%)",
        f"{_F6_TOP_SELECTION_YEAR} new active infections",
    ]].rename(
        columns={"bacterium": "Bacterium"}
    )

    body = _html_head("Figure 6B. Resistance Trends by Bacterium")
    body += _back_link()
    body += (
        f"<img src='{stem}.png' alt='Figure 6B' "
        f"style='max-width:100%; border:1px solid #ddd; border-radius:4px;'>\n"
    )
    body += f"<p class='note'>Download: <a href='{stem}.png'>PNG</a> | <a href='{stem}.svg'>SVG</a></p>\n"
    body += "<h2>Included bacteria</h2>\n"
    body += _html_table(top_table)
    footnotes = [
        f"Lines show the {_F6B_TOP_N_BACTERIA} bacteria with the highest count-weighted "
        f"any-R percentage among new active infections in {_F6_TOP_SELECTION_YEAR} across "
        "the supplied simulation summary files.",
        "For each bacterium and calendar year, the numerator is new active infections with any-R "
        "recorded in either hospital-acquired or community-acquired infection columns. The denominator "
        f"is {_F6B_DENOMINATOR_COLUMN}.",
        "Values are annual count-weighted percentages pooled across supplied runs to keep the multi-line "
        "comparison readable.",
        "Figure 6A remains the overall population-level trend; Figure 6B decomposes the any-resistance "
        "component by bacterium.",
    ]
    if problems:
        footnotes.append("Parser notes: " + " ".join(problems[:5]))
    body += _html_footnotes(footnotes)
    body += "</body></html>"
    _save(html_path, body)


def _f6c_serious_r_columns_for_slug(slug: str, available: set[str]) -> list[str]:
    return [
        f"{slug}{suffix}"
        for suffix in _F6C_SERIOUS_R_PER_BACTERIUM_SUFFIXES
        if f"{slug}{suffix}" in available
    ]


def _f6c_placeholder(
    out_dir: Path,
    message: str,
    problems: list[str] | None = None,
) -> None:
    fig_dir = out_dir / FIGURES_DIRNAME
    fig_dir.mkdir(parents=True, exist_ok=True)
    stem = "Figure_6C__serious_r_trends_by_bacterium"
    fig, ax = plt.subplots(figsize=(10, 4.2))
    ax.text(
        0.5,
        0.5,
        f"Figure 6C. Serious-R trends by bacterium\n\n{message}",
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=10.5,
        color="#555",
        bbox=dict(boxstyle="round,pad=0.6", fc="#f5f5f5", ec="#bbb"),
        wrap=True,
    )
    ax.set_axis_off()
    fig.subplots_adjust(left=0.03, right=0.97, top=0.90, bottom=0.08)
    png_path = fig_dir / f"{stem}.png"
    svg_path = fig_dir / f"{stem}.svg"
    html_path = fig_dir / f"{stem}.html"
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)

    body = _html_head("Figure 6C. Serious-R Trends by Bacterium")
    body += _back_link()
    body += (
        f"<img src='{stem}.png' alt='Figure 6C' "
        f"style='max-width:100%; border:1px solid #ddd; border-radius:4px;'>\n"
    )
    body += f"<p class='note'>Download: <a href='{stem}.png'>PNG</a> | <a href='{stem}.svg'>SVG</a></p>\n"
    footnotes = [
        message,
        "The current simulation_summary schema contains aggregate serious-R event fields, "
        f"including {_F6_SERIOUS_RESISTANCE_COLUMN}, but not bacterium-specific serious-R "
        "event counts. This placeholder is generated so the paper-output build remains explicit.",
        "Supported future numerator fields are either a vector column named "
        + ", ".join(_F6C_SERIOUS_R_VECTOR_COLUMNS)
        + " or per-bacterium columns ending "
        + ", ".join(_F6C_SERIOUS_R_PER_BACTERIUM_SUFFIXES)
        + ".",
    ]
    if problems:
        footnotes.append("Parser notes: " + " ".join(problems[:5]))
    body += _html_footnotes(footnotes)
    body += "</body></html>"
    _save(html_path, body)
    print(f"  Saved: {png_path}")
    print(f"  Saved: {svg_path}")


def _load_bacteria_serious_r_trend_rows(
    csv_path: Path,
) -> tuple[list[dict[str, object]], str | None]:
    columns = _simulation_csv_columns(csv_path)
    if columns is None:
        return [], f"{csv_path.name}: unable to read simulation summary header."
    vector_col = next((col for col in _F6C_SERIOUS_R_VECTOR_COLUMNS if col in columns), None)
    if vector_col is None:
        return _load_bacteria_resistance_trend_rows(
            csv_path,
            _f6c_serious_r_columns_for_slug,
            "new_active_infections_serious_r",
            "serious-R",
        )
    if _F6B_DENOMINATOR_COLUMN not in columns:
        return [], f"{csv_path.name}: missing {_F6B_DENOMINATOR_COLUMN}."

    optional = ["policy_option", "run_id", "simulation_year", "year", "time_in_years", "time_step"]
    usecols = [_F6B_DENOMINATOR_COLUMN, vector_col, *optional]
    try:
        df = _read_csv_selected(csv_path, usecols)
    except (FileNotFoundError, ValueError, OSError) as exc:
        return [], f"{csv_path.name}: unable to load Figure 6C columns ({exc})."
    if df.empty:
        return [], f"{csv_path.name}: no rows available."

    years = _simulation_year_series(df)
    valid_year = years.notna()
    df = df.loc[valid_year].copy()
    years = years.loc[valid_year]
    if df.empty:
        return [], f"{csv_path.name}: no rows with valid simulation year."

    target_len = len(_F15_KNOWN_BACTERIA_SLUGS)
    denominator_values = [
        _figure_15_extend_array(
            np.array(_figure_15_parse_vector_cell(value), dtype=float),
            target_len,
        )
        for value in df[_F6B_DENOMINATOR_COLUMN]
    ]
    numerator_values = [
        _figure_15_extend_array(
            np.array(_figure_15_parse_vector_cell(value), dtype=float),
            target_len,
        )
        for value in df[vector_col]
    ]
    if not denominator_values or not numerator_values:
        return [], f"{csv_path.name}: no serious-R vector values found."
    denominator_matrix = np.vstack(denominator_values)
    numerator_matrix = np.vstack(numerator_values)

    run_series = (
        df["run_id"].astype(str)
        if "run_id" in df.columns
        else pd.Series(csv_path.stem, index=df.index)
    )
    work = pd.DataFrame({
        "year_int": years.astype(int).to_numpy(),
        "run_id": run_series.to_numpy(),
    })
    records: list[dict[str, object]] = []
    for b_idx, slug in enumerate(_F15_KNOWN_BACTERIA_SLUGS):
        bacterium_work = work.copy()
        bacterium_work["denominator"] = denominator_matrix[:, b_idx]
        bacterium_work["numerator"] = numerator_matrix[:, b_idx]
        annual = (
            bacterium_work
            .groupby(["run_id", "year_int"], as_index=False)[["denominator", "numerator"]]
            .sum()
        )
        annual = annual[annual["denominator"] > 0.0].copy()
        for _, row in annual.iterrows():
            records.append({
                "source": csv_path.name,
                "run_id": row["run_id"],
                "year": int(row["year_int"]),
                "bacterium_slug": slug,
                "bacterium": _figure_15_bacterium_label(slug),
                "new_active_infections": float(row["denominator"]),
                "new_active_infections_serious_r": float(row["numerator"]),
            })

    if not records:
        return [], f"{csv_path.name}: no positive per-bacterium new-infection denominators."
    return records, None


def make_figure_6c_serious_r_trend_by_bacterium(
    csv_paths: list[Path],
    out_dir: Path,
) -> None:
    """Figure 6C: serious-R trend for the highest serious-R bacteria in 2025."""
    rows: list[dict[str, object]] = []
    problems: list[str] = []
    for csv_path in csv_paths:
        csv_rows, problem = _load_bacteria_serious_r_trend_rows(csv_path)
        rows.extend(csv_rows)
        if problem:
            problems.append(problem)

    if not rows:
        _f6c_placeholder(
            out_dir,
            "Per-bacterium serious-R trend data are not available in the supplied "
            "simulation_summary files.",
            problems,
        )
        return

    data = pd.DataFrame(rows)
    data["pct_serious_r"] = np.where(
        data["new_active_infections"] > 0.0,
        data["new_active_infections_serious_r"] / data["new_active_infections"] * 100.0,
        np.nan,
    )
    selection_year_data = data[data["year"] == _F6_TOP_SELECTION_YEAR].copy()
    if selection_year_data.empty:
        _f6c_placeholder(
            out_dir,
            f"No per-bacterium serious-R rows were found for {_F6_TOP_SELECTION_YEAR}.",
            problems,
        )
        return
    top = (
        selection_year_data
        .groupby(["bacterium_slug", "bacterium"], as_index=False)[
            ["new_active_infections", "new_active_infections_serious_r"]
        ]
        .sum()
    )
    top["selection_serious_r_pct"] = np.where(
        top["new_active_infections"] > 0.0,
        top["new_active_infections_serious_r"] / top["new_active_infections"] * 100.0,
        np.nan,
    )
    top = (
        top.dropna(subset=["selection_serious_r_pct"])
        .sort_values(
            ["selection_serious_r_pct", "new_active_infections"],
            ascending=[False, False],
        )
        .head(_F6B_TOP_N_BACTERIA)
    )
    top_slugs = top["bacterium_slug"].tolist()
    if not top_slugs:
        _f6c_placeholder(
            out_dir,
            f"No positive per-bacterium serious-R denominators were found for {_F6_TOP_SELECTION_YEAR}.",
            problems,
        )
        return

    plot_data = data[data["bacterium_slug"].isin(top_slugs)].copy()
    annual = (
        plot_data
        .groupby(["bacterium_slug", "bacterium", "year"], as_index=False)[
            ["new_active_infections", "new_active_infections_serious_r"]
        ]
        .sum()
    )
    annual["pct_serious_r"] = np.where(
        annual["new_active_infections"] > 0.0,
        annual["new_active_infections_serious_r"] / annual["new_active_infections"] * 100.0,
        np.nan,
    )
    year_min_data = int(annual["year"].min())
    year_max_data = int(annual["year"].max())

    fig_dir = out_dir / FIGURES_DIRNAME
    fig_dir.mkdir(parents=True, exist_ok=True)
    stem = "Figure_6C__serious_r_trends_by_bacterium"
    png_path = fig_dir / f"{stem}.png"
    svg_path = fig_dir / f"{stem}.svg"
    html_path = fig_dir / f"{stem}.html"

    fig, ax = plt.subplots(figsize=(11.5, 6.2))
    colour_cycle = plt.cm.tab20(np.linspace(0, 1, len(top_slugs)))
    rank_lookup = {slug: idx + 1 for idx, slug in enumerate(top_slugs)}
    for colour, slug in zip(colour_cycle, top_slugs):
        bacterium_rows = annual[annual["bacterium_slug"] == slug].sort_values("year")
        if bacterium_rows.empty:
            continue
        label = f"{rank_lookup[slug]}. {bacterium_rows['bacterium'].iloc[0]}"
        ax.plot(
            bacterium_rows["year"].to_numpy(dtype=float),
            bacterium_rows["pct_serious_r"].to_numpy(dtype=float),
            linewidth=1.7,
            marker="o",
            markersize=2.4,
            color=colour,
            label=label,
        )

    if year_min_data > _F1_SIM_EPOCH_YEAR + 80:
        ax.text(
            0.02,
            0.97,
            "Note: data currently cover calibration window only "
            f"({year_min_data}\u2013{year_max_data}).",
            transform=ax.transAxes,
            fontsize=8,
            va="top",
            color="#c0392b",
            bbox=dict(boxstyle="round,pad=0.3", fc="#fff3cd", ec="#f0c040", alpha=0.9),
        )

    ax.set_xlim(_F1_SIM_EPOCH_YEAR, _F1_SIM_EPOCH_YEAR + 96)
    ax.set_ylim(0, 100)
    ax.set_xlabel("Year", fontsize=11)
    ax.set_ylabel("New infections with serious-R marker resistance (%)", fontsize=11)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", linewidth=0.4, alpha=0.5)
    ax.legend(
        fontsize=7.5,
        frameon=False,
        loc="center left",
        bbox_to_anchor=(1.01, 0.5),
        title=f"Top {_F6B_TOP_N_BACTERIA} by {_F6_TOP_SELECTION_YEAR} serious-R",
        title_fontsize=8,
    )
    fig.suptitle(
        "Figure 6C. Serious-R trends by bacterium",
        fontsize=11,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 0.80, 1])
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {png_path}")
    print(f"  Saved: {svg_path}")

    top_table = top.copy()
    top_table["Rank"] = np.arange(1, len(top_table) + 1)
    top_table[f"{_F6_TOP_SELECTION_YEAR} serious-R (%)"] = top_table[
        "selection_serious_r_pct"
    ].map(lambda value: f"{float(value):.1f}")
    top_table[f"{_F6_TOP_SELECTION_YEAR} new active infections"] = top_table[
        "new_active_infections"
    ].map(lambda value: f"{float(value):,.0f}")
    top_table = top_table[[
        "Rank",
        "bacterium",
        f"{_F6_TOP_SELECTION_YEAR} serious-R (%)",
        f"{_F6_TOP_SELECTION_YEAR} new active infections",
    ]].rename(columns={"bacterium": "Bacterium"})

    body = _html_head("Figure 6C. Serious-R Trends by Bacterium")
    body += _back_link()
    body += (
        f"<img src='{stem}.png' alt='Figure 6C' "
        f"style='max-width:100%; border:1px solid #ddd; border-radius:4px;'>\n"
    )
    body += f"<p class='note'>Download: <a href='{stem}.png'>PNG</a> | <a href='{stem}.svg'>SVG</a></p>\n"
    body += "<h2>Included bacteria</h2>\n"
    body += _html_table(top_table)
    footnotes = [
        f"Lines show the {_F6B_TOP_N_BACTERIA} bacteria with the highest count-weighted "
        f"serious-R percentage among new active infections in {_F6_TOP_SELECTION_YEAR} across "
        "the supplied simulation summary files.",
        "For each bacterium and calendar year, the numerator is new active infections with "
        "bacterium-specific serious-R marker resistance. The denominator is "
        f"{_F6B_DENOMINATOR_COLUMN}.",
        "Values are annual count-weighted percentages pooled across supplied runs to keep the multi-line "
        "comparison readable.",
        "Figure 6C uses serious-R marker resistance, not any-R. It therefore needs "
        "bacterium-specific serious-R event numerator columns in simulation_summary files.",
    ]
    if problems:
        footnotes.append("Parser notes: " + " ".join(problems[:5]))
    body += _html_footnotes(footnotes)
    body += "</body></html>"
    _save(html_path, body)


# ---------------------------------------------------------------------------
# Figure F2 — Resistance calibration fit by organism (all organisms, dynamic grid)
# ---------------------------------------------------------------------------

# Preferred display order: IHME/WHO-ESKAPE priority organisms first, then the
# remainder alphabetically.  Any organism present in the data but not listed
# here will be appended alphabetically at the end.
_F2_ORGANISM_ORDER: list[str] = [
    "Escherichia coli",
    "Klebsiella pneumoniae",
    "Staphylococcus aureus",
    "Streptococcus pneumoniae",
    "Acinetobacter baumannii",
    "Pseudomonas aeruginosa",
    "Enterococcus faecium",
    "Enterococcus faecalis",
    "Haemophilus influenzae",
    "Enterobacter cloacae",
    "Shigella spp.",
    "Streptococcus pyogenes",
]

# Canonical display order for drug classes on the x-axis.
_F2_CLASS_ORDER: list[str] = [
    "Penicillins (J01C)",
    "Beta-lactamase combinations (J01CR)",
    "Cephalosporins 1-2G",
    "Cephalosporins 3G",
    "Cephalosporins 3G/BLI",
    "Cephalosporins 4G",
    "Anti-MRSA Cephalosporins (5G)",
    "Siderophore Cephalosporins",
    "Monobactams",
    "Carbapenems (J01DH)",
    "Novel BL/BLI",
    "Fluoroquinolones (J01M)",
    "Aminoglycosides (J01G)",
    "Tetracyclines (J01A)",
    "Macrolides (J01F)",
    "Sulfonamides (J01E)",
    "Glycopeptides (J01XA)",
    "Lincosamides (J01FF)",
    "Oxazolidinones (J01XX)",
    "Rifamycins (J04AB)",
    "Chloramphenicol (J01BA)",
    "Nitroimidazoles",
    "Polymyxins (J01XB)",
    "Lipopeptides (J01XX09)",
    "Fosfomycin (J01XX01)",
    "Nitrofurans (J01XE)",
    "Fusidic acid (J01XC)",
    "Streptogramins (J01FG)",
    "Lipoglycopeptides",
    "Pleuromutilins",
    "Fidaxomicin",
]

# Short x-axis labels for each class.
_F2_CLASS_SHORT: dict[str, str] = {
    "Penicillins (J01C)":                "Penicillins",
    "Beta-lactamase combinations (J01CR)": "BLIs",
    "Cephalosporins 1-2G":               "Ceph 1-2G",
    "Cephalosporins 3G":                 "Ceph 3G",
    "Cephalosporins 3G/BLI":             "Ceph 3G/BLI",
    "Cephalosporins 4G":                 "Ceph 4G",
    "Anti-MRSA Cephalosporins (5G)":     "Ceph 5G",
    "Siderophore Cephalosporins":        "Sider. Ceph",
    "Monobactams":                       "Monobactams",
    "Carbapenems (J01DH)":               "Carbapenems",
    "Novel BL/BLI":                      "Novel BL/BLI",
    "Fluoroquinolones (J01M)":           "FQs",
    "Aminoglycosides (J01G)":            "Aminoglycosides",
    "Tetracyclines (J01A)":              "Tetracyclines",
    "Macrolides (J01F)":                 "Macrolides",
    "Sulfonamides (J01E)":               "Sulfonamides",
    "Glycopeptides (J01XA)":             "Glycopeptides",
    "Lincosamides (J01FF)":              "Lincosamides",
    "Oxazolidinones (J01XX)":            "Oxazolidinones",
    "Rifamycins (J04AB)":                "Rifamycins",
    "Chloramphenicol (J01BA)":           "Chloramphenicol",
    "Nitroimidazoles":                   "Nitroimidazoles",
    "Polymyxins (J01XB)":                "Polymyxins",
    "Lipopeptides (J01XX09)":            "Lipopeptides",
    "Fosfomycin (J01XX01)":              "Fosfomycin",
    "Nitrofurans (J01XE)":               "Nitrofurans",
    "Fusidic acid (J01XC)":              "Fusidic acid",
    "Streptogramins (J01FG)":            "Streptogramins",
    "Lipoglycopeptides":                 "Lipoglycopeptides",
    "Pleuromutilins":                    "Pleuromutilins",
    "Fidaxomicin":                       "Fidaxomicin",
}

_F2_DISPLAY_POTENCY_THRESHOLD = 0.25
_F2_EXCLUDED_CLASSES = {"Rifamycins (J04AB)"}
_F2_EXCLUDED_BACTERIA_SLUGS = {"mdr_mycobacterium_tuberculosis"}
_F2_SUMMARY_MEDIAN_RANGE = "median_range"
_F2_SUMMARY_MEAN_CI = "mean_ci"
_F2_BACTERIA_SLUG_NORMALIZATION_OVERRIDES = {
    "p_stuartii": "providencia_stuartii",
}
_F2_DRUG_SLUG_NORMALIZATION_OVERRIDES = {
    "doxyclycline": "doxycycline",
}
_F2_BASELINE_POTENCY_LOOKUP_CACHE: dict[tuple[str, str], float] | None = None

# Colour scheme
_F2_COLOUR_SIM    = "#2196F3"   # blue — simulation
_F2_COLOUR_TARGET = "#FF7043"   # deep orange - calibration benchmark


def _parse_interval_val(v: object) -> tuple[float, float, float] | None:
    """Parse an aggregated calibration value into (median, lo, hi)."""
    if v is None:
        return None
    if isinstance(v, float) and np.isnan(v):
        return None
    if isinstance(v, (int, float)):
        f = float(v)
        return (f, f, f)
    s = str(v).strip()
    if s in ("-", "", "N/A", "nan", "\u2014", "\u00e2\u20ac\u201d"):
        return None
    for dash in ("\u2013", "\u2014", "\u2212", "\u00e2\u20ac\u201c", "\u00e2\u20ac\u201d"):
        s = s.replace(dash, "-")
    # In strings such as "19.1 (16.2-20.5)", the hyphen between interval
    # bounds is a range separator, not a negative sign for the upper bound.
    nums = re.findall(r"(?<![\d.,])[-+]?(?:\d[\d,]*\.?\d*|\.\d+)", s)
    if len(nums) >= 3:
        parsed = [float(x.replace(",", "")) for x in nums[:3]]
        return (parsed[0], parsed[1], parsed[2])
    if len(nums) == 1:
        f = float(nums[0].replace(",", ""))
        return (f, f, f)
    try:
        f = float(s.replace(",", ""))
        return (f, f, f)
    except ValueError:
        return None


def _first_numeric_value(v: object) -> float | None:
    """Return the first numeric value from a possibly mixed text cell."""
    if v is None:
        return None
    if isinstance(v, float) and np.isnan(v):
        return None
    if isinstance(v, (int, float)):
        return float(v)
    s = str(v).strip()
    if s in ("-", "", "N/A", "nan", "\u2014", "\u00e2\u20ac\u201d"):
        return None
    nums = re.findall(r"(?<![\d.,])[-+]?(?:\d[\d,]*\.?\d*|\.\d+)", s)
    if not nums:
        return None
    return float(nums[0].replace(",", ""))


def _add_interval_columns(df: pd.DataFrame, source_col: str, prefix: str) -> pd.DataFrame:
    """Add median/lo/hi numeric columns parsed from an aggregated value column."""
    df = df.copy()
    parsed = df[source_col].apply(_parse_interval_val)
    df[f"{prefix}_med"] = parsed.apply(lambda x: x[0] if x else np.nan)
    df[f"{prefix}_lo"] = parsed.apply(lambda x: x[1] if x else np.nan)
    df[f"{prefix}_hi"] = parsed.apply(lambda x: x[2] if x else np.nan)
    return df


def _asymmetric_errors(df: pd.DataFrame, prefix: str) -> tuple[np.ndarray, np.ndarray]:
    """Return non-negative asymmetric lower/upper errors for parsed interval columns."""
    med = df[f"{prefix}_med"].to_numpy(dtype=float)
    lo = df[f"{prefix}_lo"].to_numpy(dtype=float)
    hi = df[f"{prefix}_hi"].to_numpy(dtype=float)
    return np.clip(med - lo, 0, None), np.clip(hi - med, 0, None)


def _parse_resistance_val(v: object) -> tuple[float, float, float] | None:
    """
    Parse a resistance value from resistance_benchmarks into (median, lo, hi).

    Handles:
      - float/int (single run, no CI)
      - "12.3 (10.1–14.5)" or "12.3 (10.1-14.5)"  (aggregated multi-run)
      - "—", "", None, NaN → return None
    """
    return _parse_interval_val(v)


def _f2_normalize_summary_mode(value: object) -> str:
    raw = str(value or "").strip().lower().replace("-", "_").replace(" ", "_")
    aliases = {
        "": _F2_SUMMARY_MEDIAN_RANGE,
        "median": _F2_SUMMARY_MEDIAN_RANGE,
        "median_range": _F2_SUMMARY_MEDIAN_RANGE,
        "median_5_95": _F2_SUMMARY_MEDIAN_RANGE,
        "median_p5_p95": _F2_SUMMARY_MEDIAN_RANGE,
        "mean": _F2_SUMMARY_MEAN_CI,
        "mean_ci": _F2_SUMMARY_MEAN_CI,
        "mean_95ci": _F2_SUMMARY_MEAN_CI,
        "mean_95_ci": _F2_SUMMARY_MEAN_CI,
    }
    if raw not in aliases:
        valid = ", ".join(sorted({_F2_SUMMARY_MEDIAN_RANGE, _F2_SUMMARY_MEAN_CI}))
        raise ValueError(f"Unknown Figure 2 summary mode '{value}'. Use one of: {valid}.")
    return aliases[raw]


def _f2_default_summary_mode() -> str:
    return _f2_normalize_summary_mode(FIGURE2_SUMMARY_MODE)


def _f2_t_critical_975(n_values: int) -> float:
    """Two-sided 95% t critical value for n run-level values."""
    df = max(1, int(n_values) - 1)
    table = {
        1: 12.706, 2: 4.303, 3: 3.182, 4: 2.776, 5: 2.571,
        6: 2.447, 7: 2.365, 8: 2.306, 9: 2.262, 10: 2.228,
        11: 2.201, 12: 2.179, 13: 2.160, 14: 2.145, 15: 2.131,
        16: 2.120, 17: 2.110, 18: 2.101, 19: 2.093, 20: 2.086,
        21: 2.080, 22: 2.074, 23: 2.069, 24: 2.064, 25: 2.060,
        26: 2.056, 27: 2.052, 28: 2.048, 29: 2.045, 30: 2.042,
    }
    return table.get(df, 1.960 if df > 120 else 2.000)


def _f2_ci95(values: list[float]) -> tuple[float | None, float | None, float | None]:
    arr = np.array([v for v in values if np.isfinite(v)], dtype=float)
    if len(arr) == 0:
        return None, None, None
    center = float(np.mean(arr))
    if len(arr) == 1:
        return center, center, center
    se = float(np.std(arr, ddof=1) / math.sqrt(len(arr)))
    half_width = _f2_t_critical_975(len(arr)) * se
    return center, max(0.0, center - half_width), min(100.0, center + half_width)


def _f2_class_summary_from_aggregated_rows(
    rows: pd.DataFrame,
    sim_col: str,
    tgt_col: str,
) -> tuple[float | None, float | None, float | None, float | None]:
    """
    Current Figure 2 behaviour: mean of per-drug run medians, with the
    displayed interval spanning the drug-level 5th-95th percentile limits.
    """
    sims: list[tuple[float, float, float]] = []
    tgts: list[float] = []
    for _, row in rows.iterrows():
        parsed = _parse_resistance_val(row.get(sim_col))
        if parsed is not None:
            sims.append(parsed)
        tv = _parse_resistance_val(row.get(tgt_col))
        if tv is not None:
            tgts.append(tv[0])
    sim_med = float(np.mean([s[0] for s in sims])) if sims else None
    sim_lo  = float(np.min([s[1] for s in sims]))  if sims else None
    sim_hi  = float(np.max([s[2] for s in sims]))  if sims else None
    tgt     = float(np.mean(tgts)) if tgts else None
    return sim_med, sim_lo, sim_hi, tgt


def _f2_run_label(run: dict, index: int) -> str:
    meta = run.get("meta", {}) if isinstance(run, dict) else {}
    return str(meta.get("run_id") or meta.get("source_file") or f"run_{index + 1}")


def _f2_raw_resistance_rows(runs: list[dict] | None) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for idx, run in enumerate(runs or []):
        rb = run.get("resistance_benchmarks", pd.DataFrame()) if isinstance(run, dict) else pd.DataFrame()
        if rb is None or rb.empty:
            continue
        frame = rb.copy()
        frame["_f2_run"] = _f2_run_label(run, idx)
        frames.append(frame)
    if not frames:
        return pd.DataFrame()
    return pd.concat(frames, ignore_index=True, sort=False)


def _f2_build_median_range_class_table(rb: pd.DataFrame, sim_col: str, tgt_col: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for (bacterium, cls), group in rb.groupby(["Bacteria", "Class"], dropna=False, sort=False):
        sim, lo, hi, target = _f2_class_summary_from_aggregated_rows(group, sim_col, tgt_col)
        if sim is None and target is None:
            continue
        rows.append({
            "Bacteria": bacterium,
            "Class": cls,
            "sim": sim,
            "lo": lo,
            "hi": hi,
            "target": target,
        })
    return pd.DataFrame(rows)


def _f2_build_mean_ci_class_table(runs: list[dict] | None, sim_col: str, tgt_col: str) -> pd.DataFrame:
    rb = _f2_raw_resistance_rows(runs)
    if rb.empty:
        return pd.DataFrame()
    if "Flags" in rb.columns:
        rb = rb[~rb["Flags"].astype(str).str.contains("negligible", case=False, na=False)].copy()
    if rb.empty:
        return pd.DataFrame()
    rb = rb.loc[rb["Bacteria"].apply(_f2_is_valid_organism_label)].copy()
    rb = _f2_apply_display_filters(rb)
    if rb.empty:
        return pd.DataFrame()

    rb["_f2_sim"] = rb[sim_col].apply(_first_numeric_value) if sim_col in rb.columns else np.nan
    rb["_f2_target"] = rb[tgt_col].apply(_first_numeric_value) if tgt_col in rb.columns else np.nan

    run_class = (
        rb.dropna(subset=["_f2_sim"])
        .groupby(["Bacteria", "Class", "_f2_run"], as_index=False)["_f2_sim"]
        .mean()
    )
    target_by_class = (
        rb.dropna(subset=["_f2_target"])
        .groupby(["Bacteria", "Class"], as_index=False)["_f2_target"]
        .mean()
        .rename(columns={"_f2_target": "target"})
    )

    keys = set(zip(run_class["Bacteria"], run_class["Class"])) | set(
        zip(target_by_class["Bacteria"], target_by_class["Class"])
    )
    rows: list[dict[str, object]] = []
    target_lookup = {
        (row["Bacteria"], row["Class"]): float(row["target"])
        for _, row in target_by_class.iterrows()
        if pd.notna(row["target"])
    }
    for bacterium, cls in sorted(keys, key=lambda item: (str(item[0]).lower(), str(item[1]).lower())):
        values = run_class.loc[
            (run_class["Bacteria"] == bacterium) & (run_class["Class"] == cls),
            "_f2_sim",
        ].astype(float).tolist()
        center, lo, hi = _f2_ci95(values)
        target = target_lookup.get((bacterium, cls))
        if center is None and target is None:
            continue
        rows.append({
            "Bacteria": bacterium,
            "Class": cls,
            "sim": center,
            "lo": lo,
            "hi": hi,
            "target": target,
            "n_runs": len([v for v in values if np.isfinite(v)]),
        })
    return pd.DataFrame(rows)


def _f2_slugify_value(name: object) -> str:
    return str(name or "").strip().lower().replace(" ", "_")


def _f2_slugify_bacteria_value(name: object) -> str:
    slug = _f2_slugify_value(name)
    return _F2_BACTERIA_SLUG_NORMALIZATION_OVERRIDES.get(slug, slug)


def _f2_normalize_drug_slug(name: object) -> str:
    slug = _f2_slugify_value(name)
    return _F2_DRUG_SLUG_NORMALIZATION_OVERRIDES.get(slug, slug)


def _f2_load_baseline_potency_lookup() -> dict[tuple[str, str], float]:
    global _F2_BASELINE_POTENCY_LOOKUP_CACHE

    if _F2_BASELINE_POTENCY_LOOKUP_CACHE is not None:
        return _F2_BASELINE_POTENCY_LOOKUP_CACHE

    path = REPO_ROOT / "data" / "model_potency_matrix.csv"
    if not path.exists():
        print(
            f"  Figure 2: {path} not found; potency display threshold "
            "will not be applied."
        )
        _F2_BASELINE_POTENCY_LOOKUP_CACHE = {}
        return _F2_BASELINE_POTENCY_LOOKUP_CACHE

    try:
        potency_df = pd.read_csv(path)
    except Exception as exc:
        print(
            f"  Figure 2: could not read {path.name} ({exc}); potency display "
            "threshold will not be applied."
        )
        _F2_BASELINE_POTENCY_LOOKUP_CACHE = {}
        return _F2_BASELINE_POTENCY_LOOKUP_CACHE

    required = {"bacteria", "drug", "potency_when_no_r"}
    if potency_df.empty or not required.issubset(potency_df.columns):
        print(
            f"  Figure 2: {path.name} is missing required columns; potency "
            "display threshold will not be applied."
        )
        _F2_BASELINE_POTENCY_LOOKUP_CACHE = {}
        return _F2_BASELINE_POTENCY_LOOKUP_CACHE

    lookup: dict[tuple[str, str], float] = {}
    for _, row in potency_df.iterrows():
        bacteria = row.get("bacteria")
        drug = row.get("drug")
        potency = row.get("potency_when_no_r")
        if pd.isna(bacteria) or pd.isna(drug) or pd.isna(potency):
            continue
        lookup[(_f2_slugify_bacteria_value(bacteria), _f2_normalize_drug_slug(drug))] = float(potency)

    _F2_BASELINE_POTENCY_LOOKUP_CACHE = lookup
    return _F2_BASELINE_POTENCY_LOOKUP_CACHE


def _f2_class_has_display_potency(rows: pd.DataFrame) -> bool:
    lookup = _f2_load_baseline_potency_lookup()
    if not lookup:
        return True

    for _, row in rows.iterrows():
        potency = lookup.get(
            (
                _f2_slugify_bacteria_value(row.get("Bacteria")),
                _f2_normalize_drug_slug(row.get("Drug")),
            )
        )
        if potency is not None and potency >= _F2_DISPLAY_POTENCY_THRESHOLD:
            return True
    return False


def _f2_apply_display_filters(rb: pd.DataFrame) -> pd.DataFrame:
    """Apply Figure 2 presentation exclusions without changing upstream metrics."""
    if rb.empty:
        return rb

    kept_indices: list[object] = []
    dropped_excluded = 0
    dropped_potency = 0
    for (_, cls), rows in rb.groupby(["Bacteria", "Class"], dropna=False, sort=False):
        bacterium = rows["Bacteria"].iloc[0]
        bacterium_slug = _f2_slugify_bacteria_value(bacterium)
        cls_name = str(cls)

        if bacterium_slug in _F2_EXCLUDED_BACTERIA_SLUGS or cls_name in _F2_EXCLUDED_CLASSES:
            dropped_excluded += 1
            continue
        if not _f2_class_has_display_potency(rows):
            dropped_potency += 1
            continue
        kept_indices.extend(rows.index.tolist())

    if dropped_excluded:
        print(f"  Figure 2: excluded {dropped_excluded} special-case organism/class cell(s).")
    if dropped_potency:
        print(
            f"  Figure 2: excluded {dropped_potency} organism/class cell(s) with no drug "
            f"at baseline potency >= {_F2_DISPLAY_POTENCY_THRESHOLD:.2f}."
        )

    if not kept_indices:
        return rb.iloc[0:0].copy()
    return rb.loc[kept_indices].copy()


def _f2_is_valid_organism_label(v: object) -> bool:
    """Return False for separator/header artifacts that are not organism names."""
    if v is None:
        return False
    if isinstance(v, float) and np.isnan(v):
        return False
    s = str(v).strip()
    if not s or s.lower() == "nan":
        return False
    dash_normalised = (
        s.replace("\u2013", "-")
        .replace("\u2014", "-")
        .replace("\u2212", "-")
        .replace("\u00e2\u20ac\u201c", "-")
        .replace("\u00e2\u20ac\u201d", "-")
    )
    compact = "".join(ch for ch in dash_normalised if not ch.isspace())
    return set(compact) != {"-"}


def _make_figure_2_calibration_resistance_fit_legacy(agg: dict, out_dir: Path) -> None:
    """
    Create Figure 2: dynamic grid showing infection resistance calibration fit
    for all bacteria present in the resistance_benchmarks data.

    Each panel:
      x-axis — drug classes present for that organism
      y-axis — % resistant infections
      Blue bar + error bars — simulation (median ± 5th–95th percentile range)
      Orange bar           — calibration benchmark
    """
    rb = agg.get("resistance_benchmarks", pd.DataFrame())
    if rb is None or rb.empty:
        print("  F2: no resistance_benchmarks data — skipping figure.")
        return

    if "Flags" in rb.columns:
        rb = rb[~rb["Flags"].astype(str).str.contains("negligible", case=False, na=False)].copy()
    if rb.empty:
        print("  Figure 2: no non-negligible resistance_benchmarks rows — skipping figure.")
        return

    valid_organism_mask = rb["Bacteria"].apply(_f2_is_valid_organism_label)
    if not bool(valid_organism_mask.all()):
        skipped = sorted({str(v).strip() for v in rb.loc[~valid_organism_mask, "Bacteria"].dropna()})
        print(
            "  Figure 2: omitted non-organism resistance benchmark label(s): "
            f"{', '.join(repr(v) for v in skipped)}."
        )
        rb = rb.loc[valid_organism_mask].copy()
    if rb.empty:
        print("  Figure 2: no plottable organism resistance_benchmarks rows — skipping figure.")
        return

    rb = _f2_apply_display_filters(rb)
    if rb.empty:
        print("  Figure 2: no rows remain after display potency/special-case filters — skipping figure.")
        return

    n_runs = agg.get("n_runs", 1)
    sim_col = "Inf sim (%)"
    tgt_col = "Inf target (%)"

    # Build full organism list: priority order first, then remaining alphabetically.
    all_organisms_in_data = sorted(rb["Bacteria"].dropna().unique().tolist())
    ordered = [o for o in _F2_ORGANISM_ORDER if o in all_organisms_in_data]
    ordered += sorted(o for o in all_organisms_in_data if o not in ordered)

    ncols = 4
    nrows = math.ceil(len(ordered) / ncols)
    fig_height = max(14, 3.8 * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=(22, fig_height))
    axes_flat = axes.flatten()

    for panel_idx, organism in enumerate(ordered):
        ax = axes_flat[panel_idx]
        org_rows = rb[rb["Bacteria"] == organism].copy()

        if org_rows.empty:
            ax.text(0.5, 0.5, "No data", ha="center", va="center",
                    transform=ax.transAxes, fontsize=9, color="#888")
            org_escaped = organism.replace(" ", r"\ ")
            ax.set_title(f"$\\it{{{org_escaped}}}$",
                         fontsize=9, pad=4)
            ax.axis("off")
            continue

        # Collect data per class in canonical order
        classes_in_data = set(org_rows["Class"].dropna().unique())
        ordered_classes = [c for c in _F2_CLASS_ORDER if c in classes_in_data]
        # Append any unexpected classes alphabetically
        extra = sorted(c for c in classes_in_data if c not in _F2_CLASS_ORDER)
        ordered_classes += extra

        bar_classes: list[str] = []
        sim_meds:    list[float] = []
        sim_los:     list[float] = []
        sim_his:     list[float] = []
        tgt_means:   list[float | None] = []

        for cls in ordered_classes:
            cls_rows = org_rows[org_rows["Class"] == cls]
            sim_med, sim_lo, sim_hi, tgt = _class_summary(cls_rows, sim_col, tgt_col)
            # Only include if we have at least sim OR target
            if sim_med is None and tgt is None:
                continue
            bar_classes.append(cls)
            sim_meds.append(sim_med if sim_med is not None else 0.0)
            sim_los.append(sim_lo  if sim_lo  is not None else sim_meds[-1])
            sim_his.append(sim_hi  if sim_hi  is not None else sim_meds[-1])
            tgt_means.append(tgt)

        if not bar_classes:
            ax.text(0.5, 0.5, "No data", ha="center", va="center",
                    transform=ax.transAxes, fontsize=9, color="#888")
            org_escaped = organism.replace(" ", r"\ ")
            ax.set_title(f"$\\it{{{org_escaped}}}$",
                         fontsize=9, pad=4)
            continue

        x      = np.arange(len(bar_classes))
        width  = 0.38
        labels = [_F2_CLASS_SHORT.get(c, c) for c in bar_classes]

        # Sim bars with error bars (asymmetric)
        err_lo = np.clip(np.array(sim_meds) - np.array(sim_los), 0, None)
        err_hi = np.clip(np.array(sim_his) - np.array(sim_meds), 0, None)
        ax.bar(x - width / 2, sim_meds, width,
               color=_F2_COLOUR_SIM, alpha=0.85, label="Simulation",
               yerr=[err_lo, err_hi], capsize=2.5,
               error_kw={"elinewidth": 0.8, "ecolor": "#0D47A1", "capthick": 0.8})

        # Target bars (no error bar)
        tgt_vals = [t if t is not None else 0.0 for t in tgt_means]
        tgt_alpha = [1.0 if t is not None else 0.0 for t in tgt_means]
        # Draw target bars, making bars with no target invisible
        for i, (tv, ta) in enumerate(zip(tgt_vals, tgt_alpha)):
            ax.bar(x[i] + width / 2, tv, width,
                   color=_F2_COLOUR_TARGET, alpha=0.85 * ta, label=None)

        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=5.5)
        ax.set_ylim(0, 105)
        ax.set_ylabel("Resistance (%)", fontsize=7, labelpad=2)
        ax.yaxis.set_tick_params(labelsize=6)
        ax.axhline(y=100, color="#ccc", linewidth=0.4, linestyle="--")
        ax.grid(axis="y", linewidth=0.35, alpha=0.5)
        ax.spines[["top", "right"]].set_visible(False)

        # Italicised title with genus abbreviation for long names
        title_str = organism
        ax.set_title(title_str, fontsize=8.5, fontstyle="italic", pad=3)

    # Hide any unused panels
    for idx in range(len(ordered), nrows * ncols):
        axes_flat[idx].axis("off")

    # Figure-level legend and title
    sim_patch = mpatches.Patch(color=_F2_COLOUR_SIM,   alpha=0.85, label="Simulation")
    tgt_patch = mpatches.Patch(color=_F2_COLOUR_TARGET, alpha=0.85, label="Calibration benchmark")
    ci_note   = (
        f"Error bars: 5th–95th percentile across {n_runs} run{'s' if n_runs > 1 else ''}."
        if n_runs > 1 else "Single run; no uncertainty interval shown."
    )
    fig.legend(
        handles=[sim_patch, tgt_patch], loc="lower center",
        ncol=2, fontsize=9, frameon=False, bbox_to_anchor=(0.5, 0.0),
    )
    fig.suptitle(
        "Figure 2. Calibration: resistance fit by bacterium and drug class",
        fontsize=11, fontweight="bold", y=1.01,
    )
    fig.text(0.5, -0.01, ci_note, ha="center", fontsize=7.5, color="#555")

    fig.tight_layout(rect=[0, 0.04, 1, 1])

    # Save PNG
    fig_dir = out_dir / FIGURES_DIRNAME
    fig_dir.mkdir(parents=True, exist_ok=True)
    stem = "Figure_2__calibration_resistance_fit_by_bacteria_drug_class"
    png_path = fig_dir / f"{stem}.png"
    svg_path = fig_dir / f"{stem}.svg"
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {png_path}")
    print(f"  Saved: {svg_path}")

    # HTML wrapper
    html_rel_img = f"{stem}.png"
    html_rel_svg = f"{stem}.svg"
    body  = _html_head("Figure 2. Calibration: Resistance Fit")
    body += _back_link()
    figure_note = (
        "Each panel shows the simulated (blue) and calibration-benchmark (orange) "
        "infection resistance percentage by drug class for one bacterium. "
        "Simulation bars show the mean across all drugs in the class; "
        f"error bars span the 5th–95th percentile range across {n_runs} accepted run"
        f"{'s' if n_runs > 1 else ''}. "
        "Classes without data for a given bacterium are omitted."
    )
    body += (
        f"<img src='{html_rel_img}' alt='Figure 2' "
        f"style='max-width:100%; border:1px solid #ddd; border-radius:4px;'>\n"
    )
    body += f"<p class='note'>Download: <a href='{html_rel_img}'>PNG</a> | <a href='{html_rel_svg}'>SVG</a></p>\n"
    body += _html_footnotes([
        _meta_footnote(agg),
        figure_note,
        "Drug class resistance within a panel is averaged across all specific drugs in that class.",
        "Drugs flagged as negligible potency (baseline potency < 0.15) are excluded from class averages.",
        f"Drug classes are shown only where at least one drug in the class has baseline potency >= "
        f"{_F2_DISPLAY_POTENCY_THRESHOLD:.2f} for that bacterium.",
        "Rifamycins and MDR Mycobacterium tuberculosis are omitted from this broad calibration figure "
        "because they are special-case TB/rifampicin-resistance outputs rather than comparable "
        "general drug-class calibration cells.",
        "Eligible bacteria with resistance benchmark data in the simulation output are included. "
        "IHME/WHO-ESKAPE priority bacteria are shown first, remainder alphabetically.",
    ] + _RESISTANCE_TARGET_SOURCE_NOTES)
    body += "</body></html>"
    html_path = fig_dir / f"{stem}.html"
    _save(html_path, body)


# ---------------------------------------------------------------------------
# Active Figure 2 implementation with selectable uncertainty summary
# ---------------------------------------------------------------------------

def make_figure_2_calibration_resistance_fit(
    agg: dict,
    out_dir: Path,
    *,
    runs: list[dict] | None = None,
    summary_mode: str | None = None,
) -> None:
    """
    Create Figure 2 with a selectable uncertainty summary.

    Default mode preserves the existing median/range display from the
    aggregated calibration summary. Mean-CI mode uses the per-run parsed
    resistance benchmark tables to compute one class mean per run, then
    plots the mean and a two-sided 95% t confidence interval.
    """
    mode = _f2_normalize_summary_mode(summary_mode or _f2_default_summary_mode())
    n_runs = int(agg.get("n_runs", 1) or 1)
    sim_col = "Inf sim (%)"
    tgt_col = "Inf target (%)"

    class_summary = pd.DataFrame()
    if mode == _F2_SUMMARY_MEAN_CI:
        class_summary = _f2_build_mean_ci_class_table(runs, sim_col, tgt_col)
        if class_summary.empty:
            print(
                "  Figure 2: mean-CI mode needs per-run resistance_benchmarks data; "
                "falling back to median/range mode."
            )
            mode = _F2_SUMMARY_MEDIAN_RANGE

    if mode == _F2_SUMMARY_MEDIAN_RANGE:
        rb = agg.get("resistance_benchmarks", pd.DataFrame())
        if rb is None or rb.empty:
            print("  F2: no resistance_benchmarks data - skipping figure.")
            return

        if "Flags" in rb.columns:
            rb = rb[~rb["Flags"].astype(str).str.contains("negligible", case=False, na=False)].copy()
        if rb.empty:
            print("  Figure 2: no non-negligible resistance_benchmarks rows - skipping figure.")
            return

        valid_organism_mask = rb["Bacteria"].apply(_f2_is_valid_organism_label)
        if not bool(valid_organism_mask.all()):
            skipped = sorted({str(v).strip() for v in rb.loc[~valid_organism_mask, "Bacteria"].dropna()})
            print(
                "  Figure 2: omitted non-organism resistance benchmark label(s): "
                f"{', '.join(repr(v) for v in skipped)}."
            )
            rb = rb.loc[valid_organism_mask].copy()
        if rb.empty:
            print("  Figure 2: no plottable organism resistance_benchmarks rows - skipping figure.")
            return

        rb = _f2_apply_display_filters(rb)
        if rb.empty:
            print("  Figure 2: no rows remain after display potency/special-case filters - skipping figure.")
            return

        class_summary = _f2_build_median_range_class_table(rb, sim_col, tgt_col)

    if class_summary.empty:
        print("  Figure 2: no plottable resistance_benchmarks class summaries - skipping figure.")
        return

    print(f"  Figure 2: summary mode = {mode}.")

    all_organisms_in_data = sorted(class_summary["Bacteria"].dropna().unique().tolist())
    ordered = [o for o in _F2_ORGANISM_ORDER if o in all_organisms_in_data]
    ordered += sorted(o for o in all_organisms_in_data if o not in ordered)

    ncols = 4
    nrows = math.ceil(len(ordered) / ncols)
    fig_height = max(14, 3.8 * nrows)
    fig, axes = plt.subplots(nrows, ncols, figsize=(22, fig_height))
    axes_flat = np.array(axes).flatten()

    for panel_idx, organism in enumerate(ordered):
        ax = axes_flat[panel_idx]
        org_rows = class_summary[class_summary["Bacteria"] == organism].copy()

        if org_rows.empty:
            ax.text(
                0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes, fontsize=9, color="#888",
            )
            ax.set_title(organism, fontsize=8.5, fontstyle="italic", pad=3)
            ax.axis("off")
            continue

        classes_in_data = set(org_rows["Class"].dropna().unique())
        ordered_classes = [c for c in _F2_CLASS_ORDER if c in classes_in_data]
        ordered_classes += sorted(c for c in classes_in_data if c not in _F2_CLASS_ORDER)

        bar_classes: list[str] = []
        sim_centers: list[float] = []
        sim_los: list[float] = []
        sim_his: list[float] = []
        tgt_means: list[float | None] = []

        for cls in ordered_classes:
            cls_rows = org_rows[org_rows["Class"] == cls]
            if cls_rows.empty:
                continue
            row = cls_rows.iloc[0]
            sim = row.get("sim")
            lo = row.get("lo")
            hi = row.get("hi")
            target = row.get("target")
            sim_val = float(sim) if pd.notna(sim) else None
            lo_val = float(lo) if pd.notna(lo) else sim_val
            hi_val = float(hi) if pd.notna(hi) else sim_val
            target_val = float(target) if pd.notna(target) else None

            if sim_val is None and target_val is None:
                continue
            bar_classes.append(cls)
            sim_centers.append(sim_val if sim_val is not None else 0.0)
            sim_los.append(lo_val if lo_val is not None else sim_centers[-1])
            sim_his.append(hi_val if hi_val is not None else sim_centers[-1])
            tgt_means.append(target_val)

        if not bar_classes:
            ax.text(
                0.5, 0.5, "No data", ha="center", va="center",
                transform=ax.transAxes, fontsize=9, color="#888",
            )
            ax.set_title(organism, fontsize=8.5, fontstyle="italic", pad=3)
            continue

        x = np.arange(len(bar_classes))
        width = 0.38
        labels = [_F2_CLASS_SHORT.get(c, c) for c in bar_classes]

        err_lo = np.clip(np.array(sim_centers) - np.array(sim_los), 0, None)
        err_hi = np.clip(np.array(sim_his) - np.array(sim_centers), 0, None)
        ax.bar(
            x - width / 2,
            sim_centers,
            width,
            color=_F2_COLOUR_SIM,
            alpha=0.85,
            label=None,
            yerr=[err_lo, err_hi],
            capsize=2.5,
            error_kw={"elinewidth": 0.8, "ecolor": "#0D47A1", "capthick": 0.8},
        )

        tgt_vals = [t if t is not None else 0.0 for t in tgt_means]
        tgt_alpha = [1.0 if t is not None else 0.0 for t in tgt_means]
        for i, (tv, ta) in enumerate(zip(tgt_vals, tgt_alpha)):
            ax.bar(
                x[i] + width / 2,
                tv,
                width,
                color=_F2_COLOUR_TARGET,
                alpha=0.85 * ta,
                label=None,
            )

        ax.set_xticks(x)
        ax.set_xticklabels(labels, rotation=45, ha="right", fontsize=5.5)
        ax.set_ylim(0, 105)
        ax.set_ylabel("Resistance (%)", fontsize=7, labelpad=2)
        ax.yaxis.set_tick_params(labelsize=6)
        ax.axhline(y=100, color="#ccc", linewidth=0.4, linestyle="--")
        ax.grid(axis="y", linewidth=0.35, alpha=0.5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(organism, fontsize=8.5, fontstyle="italic", pad=3)

    for idx in range(len(ordered), nrows * ncols):
        axes_flat[idx].axis("off")

    sim_label = "Simulation mean" if mode == _F2_SUMMARY_MEAN_CI else "Simulation median"
    sim_patch = mpatches.Patch(color=_F2_COLOUR_SIM, alpha=0.85, label=sim_label)
    tgt_patch = mpatches.Patch(color=_F2_COLOUR_TARGET, alpha=0.85, label="Calibration benchmark")
    if mode == _F2_SUMMARY_MEAN_CI:
        ci_note = (
            f"Error bars: 95% confidence interval for the mean across {n_runs} stochastic run"
            f"{'s' if n_runs != 1 else ''}."
            if n_runs > 1 else "Single run; mean equals the single run value and no confidence interval is shown."
        )
    else:
        ci_note = (
            f"Error bars: aggregated 5th-95th percentile range across {n_runs} accepted run"
            f"{'s' if n_runs != 1 else ''}."
            if n_runs > 1 else "Single run; no uncertainty interval shown."
        )

    fig.legend(
        handles=[sim_patch, tgt_patch],
        loc="lower center",
        ncol=2,
        fontsize=9,
        frameon=False,
        bbox_to_anchor=(0.5, 0.0),
    )
    fig.suptitle(
        "Figure 2. Calibration: resistance fit by bacterium and drug class",
        fontsize=11,
        fontweight="bold",
        y=1.01,
    )
    fig.text(0.5, -0.01, ci_note, ha="center", fontsize=7.5, color="#555")
    fig.tight_layout(rect=[0, 0.04, 1, 1])

    fig_dir = out_dir / FIGURES_DIRNAME
    fig_dir.mkdir(parents=True, exist_ok=True)
    stem = "Figure_2__calibration_resistance_fit_by_bacteria_drug_class"
    png_path = fig_dir / f"{stem}.png"
    svg_path = fig_dir / f"{stem}.svg"
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {png_path}")
    print(f"  Saved: {svg_path}")

    html_rel_img = f"{stem}.png"
    html_rel_svg = f"{stem}.svg"
    body = _html_head("Figure 2. Calibration: Resistance Fit")
    body += _back_link()
    if mode == _F2_SUMMARY_MEAN_CI:
        figure_note = (
            "Each panel shows the simulated (blue) and calibration-benchmark (orange) "
            "infection resistance percentage by drug class for one bacterium. "
            "Simulation bars show the mean of run-level class means; error bars show "
            "a two-sided 95% t confidence interval across stochastic runs. "
            "Classes without data for a given bacterium are omitted."
        )
    else:
        figure_note = (
            "Each panel shows the simulated (blue) and calibration-benchmark (orange) "
            "infection resistance percentage by drug class for one bacterium. "
            "Simulation bars show the mean across all drugs in the class using the "
            "aggregated calibration-summary median; error bars retain the aggregated "
            "5th-95th percentile range. Classes without data for a given bacterium "
            "are omitted."
        )
    body += (
        f"<img src='{html_rel_img}' alt='Figure 2' "
        f"style='max-width:100%; border:1px solid #ddd; border-radius:4px;'>\n"
    )
    body += f"<p class='note'>Download: <a href='{html_rel_img}'>PNG</a> | <a href='{html_rel_svg}'>SVG</a></p>\n"
    body += _html_footnotes([
        _meta_footnote(agg),
        figure_note,
        (
            f"Figure 2 summary mode: {mode}. To switch modes, edit "
            "FIGURE2_SUMMARY_MODE in make_paper_tables.py."
        ),
        "Drug class resistance within a panel is averaged across all specific drugs in that class.",
        "Drugs flagged as negligible potency (baseline potency < 0.15) are excluded from class averages.",
        f"Drug classes are shown only where at least one drug in the class has baseline potency >= "
        f"{_F2_DISPLAY_POTENCY_THRESHOLD:.2f} for that bacterium.",
        "Rifamycins and MDR Mycobacterium tuberculosis are omitted from this broad calibration figure "
        "because they are special-case TB/rifampicin-resistance outputs rather than comparable "
        "general drug-class calibration cells.",
        "Eligible bacteria with resistance benchmark data in the simulation output are included. "
        "IHME/WHO-ESKAPE priority bacteria are shown first, remainder alphabetically.",
    ] + _RESISTANCE_TARGET_SOURCE_NOTES)
    body += "</body></html>"
    html_path = fig_dir / f"{stem}.html"
    _save(html_path, body)


# ---------------------------------------------------------------------------
# Figure F3 - Calibration block scores
# ---------------------------------------------------------------------------

# Legacy calibration-score figure retained for reference only; not called by main().
def make_f3_calibration_scores(agg: dict, out_dir: Path) -> None:
    """
    Figure F3: horizontal bar chart of calibration block scores.
    Score ≤ 1.0 = accepted; > 1.0 = failed.
    """
    bs = agg.get("block_scores", pd.DataFrame())
    n  = agg.get("n_runs", 1)
    if bs is None or bs.empty:
        print("  F3: no block_scores data — skipping.")
        return
    block_col = bs.columns[0]
    score_col = next((c for c in bs.columns if c.lower() == "score"), None)
    if score_col is None:
        print("  F3: no 'Score' column — skipping.")
        return
    bs = bs[[block_col, score_col]].copy()
    bs = _add_interval_columns(bs, score_col, "_score")
    bs = bs.dropna(subset=["_score_med"]).sort_values("_score_med", ascending=True)
    colors = ["#EF5350" if v > 1.0 else "#42A5F5" for v in bs["_score_med"]]
    fig, ax = plt.subplots(figsize=(8, max(2.5, 0.7 * len(bs))))
    err_lo, err_hi = _asymmetric_errors(bs, "_score")
    bars = ax.barh(
        range(len(bs)),
        bs["_score_med"].values,
        color=colors,
        edgecolor="none",
        height=0.6,
        xerr=[err_lo, err_hi] if n > 1 else None,
        error_kw={"elinewidth": 1.0, "ecolor": "#333", "capthick": 1.0, "capsize": 3},
    )
    ax.set_yticks(range(len(bs)))
    ax.set_yticklabels(bs[block_col].values, fontsize=10)
    ax.axvline(1.0, color="#333", linewidth=1.2, linestyle="--", alpha=0.8)
    ax.set_xlabel("Block score (lower = better fit)", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", linewidth=0.4, alpha=0.5)
    for bar, val in zip(bars, bs["_score_med"]):
        ax.text(val + 0.02, bar.get_y() + bar.get_height() / 2,
                f"{val:.3f}", va="center", ha="left", fontsize=9)
    accepted = mpatches.Patch(color="#42A5F5", label="Score \u2264 1.0 (accepted)")
    failed   = mpatches.Patch(color="#EF5350", label="Score > 1.0 (failed)")
    thresh   = plt.Line2D([0], [0], color="#333", linewidth=1.2,
                          linestyle="--", label="Acceptance threshold (1.0)")
    ax.legend(handles=[accepted, failed, thresh], fontsize=8, frameon=False)
    fig.suptitle("Figure F3 \u2014 Calibration block scores",
                 fontsize=11, fontweight="bold")
    fig.tight_layout()
    _save_figure(
        fig, out_dir, "F3_calibration_scores",
        "Figure F3 \u2014 Calibration Block Scores",
        f"Normalised block scores for the accepted calibration run{'s' if n > 1 else ''} "
        f"(n\u2009=\u2009{n}). Scores \u2264 1.0 (dashed line) indicate all targets within "
        f"the block are within tolerance. Blue = accepted; red = failed. "
        f"{'Error bars show 5th-95th percentile ranges across accepted runs.' if n > 1 else ''}",
        [
            "Block scores are normalised so that 1.0 represents the acceptance boundary. "
            "A score below 1.0 means all targets in that calibration block are within the "
            "specified absolute or relative tolerance.",
            "The overall calibration score is the weighted sum across blocks. Calibration "
            "is considered accepted when the overall weighted score \u2264 1.0.",
        ],
    )


# ---------------------------------------------------------------------------
# Figure F4 — Headline metrics (figure version of T2)
# ---------------------------------------------------------------------------

def make_figure_1_calibration_headline_metrics(agg: dict, out_dir: Path) -> None:
    """
    Figure 1: 4-panel grouped bar chart for headline calibration metrics
    (figure version of Table T2).
    """
    hm = agg.get("headline_metrics", pd.DataFrame()).copy()
    n  = agg.get("n_runs", 1)
    if hm is None or hm.empty:
        print("  F4: no headline_metrics data — skipping.")
        return
    metric_col = hm.columns[0]
    sim_col    = next((c for c in hm.columns if "simulation" in c.lower()), None)
    tgt_col    = next((c for c in hm.columns
                       if "target" in c.lower() or "observed estimate" in c.lower()
                       or "observed data" in c.lower()), None)
    if sim_col is None:
        return

    wanted = [
        "infection deaths",
        "people on antibiotics",
        "incidence of bacterial infection",
        "incident cases of sepsis",
    ]
    hm["_metric_order"] = hm[metric_col].astype(str).str.lower().apply(
        lambda s: next((i for i, needle in enumerate(wanted) if needle in s), np.nan)
    )
    hm = hm.dropna(subset=["_metric_order"]).sort_values("_metric_order").drop(columns=["_metric_order"])
    if hm.empty:
        print("  Figure 2: requested headline metrics not found — skipping.")
        return

    def _short(raw: str) -> str:
        lo = re.sub(r'\s*\(\d+\)\s*$', '', str(raw).strip()).lower()
        if 'infection deaths'    in lo: return 'Infection\ndeaths (M/yr)'
        if 'antibiotics'         in lo: return 'People on\nantibiotics (M)'
        if 'incidence'           in lo: return 'Infection\nincidence (%/yr)'
        if 'sepsis'              in lo: return 'Sepsis\ncases (M/yr)'
        return re.sub(r'\s*\(\d+\)\s*$', '', str(raw).strip())

    hm["_raw_metric"] = hm[metric_col].astype(str)
    hm[metric_col] = hm[metric_col].apply(_short)
    ncols = len(hm)
    fig, axes = plt.subplots(1, ncols, figsize=(3.8 * ncols, 4.5))
    if ncols == 1:
        axes = [axes]

    for idx, (_, row) in enumerate(hm.iterrows()):
        ax   = axes[idx]
        name = row[metric_col]
        sim_p = _parse_resistance_val(row.get(sim_col))
        tgt_p = _parse_resistance_val(row.get(tgt_col)) if tgt_col else None
        if tgt_p is None and tgt_col:
            tgt_first = _first_numeric_value(row.get(tgt_col))
            if tgt_first is not None:
                tgt_p = (tgt_first, tgt_first, tgt_first)

        vals, labels, colors, err_lo, err_hi = [], [], [], [], []
        if sim_p:
            vals.append(sim_p[0]); labels.append("Simulation"); colors.append("#2196F3")
            err_lo.append(max(0, sim_p[0] - sim_p[1]))
            err_hi.append(max(0, sim_p[2] - sim_p[0]))
        if tgt_p:
            vals.append(tgt_p[0]); labels.append("Target/\nobserved"); colors.append("#FF7043")
            err_lo.append(0.0); err_hi.append(0.0)
        if not vals:
            ax.axis("off"); continue

        x = np.arange(len(vals))
        ax.bar(x, vals, color=colors, width=0.55, edgecolor="none", alpha=0.88)
        if sim_p and err_lo[0] + err_hi[0] > 0:
            ax.errorbar([0], [vals[0]], yerr=[[err_lo[0]], [err_hi[0]]],
                        fmt="none", color="#0D47A1", capsize=5, linewidth=1.5)
        ax.set_xticks(x)
        ax.set_xticklabels(labels, fontsize=9)
        for i, v in enumerate(vals):
            ax.text(i, v * 1.04, f"{v:.1f}", ha="center", va="bottom", fontsize=9)
        ax.set_title(name, fontsize=9.5, pad=5)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", linewidth=0.4, alpha=0.5)
        ax.set_ylabel("Value", fontsize=8)

    fig.suptitle(
        "Figure 1. Calibration: 2025 headline health and antibiotic-use metrics",
        fontsize=11, fontweight="bold",
    )
    fig.tight_layout()
    _save_figure(
        fig, out_dir, "Figure_1__calibration_headline_metrics",
        "Figure 1. Calibration: 2025 headline health and antibiotic-use metrics",
        f"Simulation compared with target/observed estimates. Error bars show 5th\u201395th "
        f"percentile range across {n} accepted run{'s' if n > 1 else ''}.",
        _HEADLINE_TARGET_SOURCE_NOTES,
        agg=agg,
    )


# ---------------------------------------------------------------------------
# Figure F5 — Drug class share (figure version of T3)
# ---------------------------------------------------------------------------

_F5_AWARE: dict[str, str] = {
    "Penicillins (J01C)":                      "#4CAF50",   # Access
    "Beta-lactamase combinations (J01CR)":      "#4CAF50",
    "Cephalosporins 1-2G":                      "#4CAF50",
    "Sulfonamides (J01E)":                      "#4CAF50",
    "Nitrofurans (J01XE)":                      "#4CAF50",
    "Fosfomycin (J01XX01)":                     "#4CAF50",
    "Chloramphenicol (J01BA)":                  "#4CAF50",
    "Cephalosporins 3G":                        "#FF9800",   # Watch
    "Cephalosporins 3G/BLI":                    "#FF9800",
    "Cephalosporins 4G":                        "#FF9800",
    "Fluoroquinolones (J01M)":                  "#FF9800",
    "Macrolides (J01F)":                        "#FF9800",
    "Glycopeptides (J01XA)":                    "#FF9800",
    "Aminoglycosides (J01G)":                   "#FF9800",
    "Tetracyclines (J01A)":                     "#FF9800",
    "Carbapenems (J01DH)":                      "#FF9800",
    "Rifamycins (J04AB)":                       "#FF9800",
    "Nitroimidazoles":                          "#FF9800",
    "Lincosamides (J01FF)":                     "#FF9800",
    "Anti-MRSA Cephalosporins (5G)":            "#F44336",   # Reserve
    "Siderophore Cephalosporins":               "#F44336",
    "Novel BL/BLI":                             "#F44336",
    "Monobactams":                              "#F44336",
    "Polymyxins (J01XB)":                       "#F44336",
    "Lipopeptides (J01XX09)":                   "#F44336",
    "Oxazolidinones (J01XX)":                   "#F44336",
    "Streptogramins (J01FG)":                   "#F44336",
    "Lipoglycopeptides":                        "#F44336",
    "Pleuromutilins":                           "#F44336",
    "Fidaxomicin":                              "#F44336",
    "Fusidic acid (J01XC)":                     "#9E9E9E",   # Not classified
}


def make_figure_3_calibration_drug_class_share(agg: dict, out_dir: Path) -> None:
    """
    Figure 3: paired horizontal bars of drug-class share.
    Bars coloured by WHO AWaRe category.
    """
    dc_raw = agg.get("drug_class_share", pd.DataFrame()).copy()
    n  = agg.get("n_runs", 1)
    if dc_raw is None or dc_raw.empty:
        print("  F5: no drug_class_share data — skipping.")
        return
    dc = _clean_df(dc_raw)
    class_col = dc.columns[0]
    sim_col   = next((c for c in dc.columns
                      if "share" in c.lower() and "%" in c
                      and "observed estimate" not in c.lower()
                      and "target" not in c.lower()), None)
    tgt_col   = next((c for c in dc.columns
                      if ("target" in c.lower() or "observed estimate" in c.lower()) and "%" in c), None)
    if sim_col is None:
        return

    plot_df = dc.copy()

    plot_df = _add_interval_columns(plot_df, sim_col, "_sim")
    if tgt_col:
        plot_df = _add_interval_columns(plot_df, tgt_col, "_tgt")
    sort_col = "_tgt_med" if tgt_col else "_sim_med"
    order = plot_df.sort_values(sort_col, ascending=True, na_position="first").index
    plot_df = plot_df.loc[order].reset_index(drop=True)
    dc = dc.loc[order].reset_index(drop=True)
    n_cls = len(plot_df)
    fig, ax = plt.subplots(figsize=(9, max(4, 0.45 * n_cls)))
    y     = np.arange(n_cls)
    bar_h = 0.36
    sim_colors = [_F5_AWARE.get(str(c), "#78909C") for c in plot_df[class_col]]
    sim_err_lo, sim_err_hi = _asymmetric_errors(plot_df, "_sim")
    ax.barh(
        y + bar_h / 2,
        plot_df["_sim_med"].fillna(0),
        bar_h,
        color=sim_colors,
        alpha=0.85,
        label="Simulation",
        xerr=[sim_err_lo, sim_err_hi] if n > 1 else None,
        error_kw={"elinewidth": 0.9, "ecolor": "#263238", "capthick": 0.9, "capsize": 2.5},
    )
    if tgt_col:
        ax.barh(y - bar_h / 2, plot_df["_tgt_med"].fillna(0), bar_h,
                color="none", edgecolor="#444", linewidth=0.9, hatch="///",
                label="Target (estimate)")
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df[class_col].values, fontsize=7.5)
    ax.set_xlabel("Share of active antibiotic drug exposure (%)", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", linewidth=0.4, alpha=0.5)
    access_p  = mpatches.Patch(color="#4CAF50", alpha=0.85, label="Access")
    watch_p   = mpatches.Patch(color="#FF9800", alpha=0.85, label="Watch")
    reserve_p = mpatches.Patch(color="#F44336", alpha=0.85, label="Reserve")
    other_p   = mpatches.Patch(color="#9E9E9E", alpha=0.85, label="Not classified")
    tgt_p     = mpatches.Patch(facecolor="none", edgecolor="#444",
                                hatch="///", label="Target (estimate)")
    handles = [access_p, watch_p, reserve_p, other_p, tgt_p]
    if n > 1:
        handles.append(plt.Line2D([0], [0], color="#263238", linewidth=1.0,
                                  label="Simulation 5th-95th percentile"))
    ax.legend(handles=handles,
              fontsize=7.5, frameon=False, loc="lower right")
    fig.suptitle("Figure 3. Calibration: 2025 antibiotic use by drug class", fontsize=10, fontweight="bold")
    fig.tight_layout()
    _save_figure(
        fig, out_dir, "Figure_3__calibration_drug_class_share",
        "Figure 3. Calibration: 2025 antibiotic use by drug class",
        f"Bars coloured by WHO AWaRe category "
        f"(green = Access, orange = Watch, red = Reserve). "
        f"Simulation shares use total active antibiotic drug-days across configured classes as the denominator. "
        f"Hatched outlines = global target/observed estimates. "
        f"{'Error bars show simulation 5th-95th percentile ranges. ' if n > 1 else ''}"
        f"n\u2009=\u2009{n} run{'s' if n > 1 else ''}.",
        [
            "WHO AWaRe classification: Access, Watch, Reserve (WHO 2023 AWaRe antibiotic book).",
            "Class-level global shares are approximate calibration anchors and carry substantial uncertainty.",
        ] + _DRUG_CLASS_TARGET_SOURCE_NOTES,
        agg=agg,
    )


# ---------------------------------------------------------------------------
# Figure F6 — Bacteria infection calibration scatter (figure version of T4)
# ---------------------------------------------------------------------------

# Legacy infection-prevalence scatter; retained for reference, not called by main().
def make_f6_bacteria_scatter(agg: dict, out_dir: Path) -> None:
    """
    Figure F6: scatter of simulated vs. target infection prevalence for all 42
    organisms (calibration fit overview — figure version of Table T4).
    """
    bi = agg.get("bacteria_infections", pd.DataFrame()).copy()
    n  = agg.get("n_runs", 1)
    if bi is None or bi.empty:
        print("  F6: no bacteria_infections data — skipping.")
        return
    bact_col = bi.columns[0]
    tgt_col  = next((c for c in bi.columns
                     if "infection" in c.lower()
                     and ("target" in c.lower() or "observed" in c.lower())), None)
    sim_col  = next((c for c in bi.columns
                     if "infection simulation" in c.lower()
                     or "infection sim" in c.lower()), None)
    if tgt_col is None or sim_col is None:
        print("  F6: cannot identify infection target/simulation columns — skipping.")
        return
    bi = _add_interval_columns(bi, tgt_col, "_tgt")
    bi = _add_interval_columns(bi, sim_col, "_sim")
    bi["_tgt"] = bi["_tgt_med"]
    bi["_sim"] = bi["_sim_med"]
    bi = bi.dropna(subset=["_tgt", "_sim"])
    bi = bi[bi["_tgt"] > 0].copy()
    if bi.empty:
        print("  F6: no valid bacteria_infections rows after filtering — skipping.")
        return
    bi["_ratio"] = bi["_sim"] / bi["_tgt"]
    raw_max = max(bi["_tgt"].max(), bi["_sim"].max())
    if not np.isfinite(raw_max) or raw_max <= 0:
        print("  F6: bacteria_infections values are NaN/Inf/zero — skipping.")
        return
    max_val = raw_max * 1.15
    colors  = ["#EF5350" if r > 2.0
               else "#FFA726" if r > 1.25
               else "#42A5F5" if r < 0.5
               else "#66BB6A"
               for r in bi["_ratio"]]
    fig, ax = plt.subplots(figsize=(8, 7))
    tgt_err_lo, tgt_err_hi = _asymmetric_errors(bi, "_tgt")
    sim_err_lo, sim_err_hi = _asymmetric_errors(bi, "_sim")
    if n > 1:
        ax.errorbar(
            bi["_tgt"],
            bi["_sim"],
            xerr=[tgt_err_lo, tgt_err_hi],
            yerr=[sim_err_lo, sim_err_hi],
            fmt="none",
            ecolor="#455A64",
            elinewidth=0.7,
            capsize=2.0,
            alpha=0.55,
            zorder=2,
        )
    ax.scatter(bi["_tgt"], bi["_sim"], c=colors, s=50,
               edgecolors="white", linewidths=0.4, zorder=3)
    ax.plot([0, max_val], [0, max_val], color="#555", linewidth=1.0,
            linestyle="--", label="Perfect fit (y = x)", zorder=2)
    ax.plot([0, max_val], [0, 2 * max_val], color="#FF9800", linewidth=0.7,
            linestyle=":", alpha=0.7, label="\xd72 / \xf70.5 tolerance")
    ax.plot([0, max_val / 2], [0, max_val], color="#FF9800", linewidth=0.7,
            linestyle=":", alpha=0.7)
    # Label outliers
    for _, row in bi[(bi["_ratio"] > 2.0) | (bi["_ratio"] < 0.5)].iterrows():
        parts = str(row[bact_col]).split()
        abbr  = (parts[0][0] + ". " + " ".join(parts[1:])) if len(parts) > 1 else str(row[bact_col])
        ax.annotate(abbr, (row["_tgt"], row["_sim"]),
                    fontsize=6, xytext=(4, 2), textcoords="offset points",
                    color="#333", clip_on=True)
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ax.set_xlabel("Infection prevalence — target (% world population)", fontsize=10)
    ax.set_ylabel("Infection prevalence — simulation (% world population)", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(linewidth=0.35, alpha=0.5)
    good_p   = mpatches.Patch(color="#66BB6A", label="0.5\xd7 \u2013 1.25\xd7 target")
    over_p   = mpatches.Patch(color="#FFA726", label="1.25\xd7 \u2013 2\xd7 target")
    way_over = mpatches.Patch(color="#EF5350", label=">2\xd7 target")
    under_p  = mpatches.Patch(color="#42A5F5", label="<0.5\xd7 target")
    ax.legend(handles=[good_p, over_p, way_over, under_p], fontsize=8, frameon=False)
    fig.suptitle(
        "Figure F6 \u2014 Bacterial infection prevalence: simulation vs. calibration target",
        fontsize=10, fontweight="bold")
    fig.tight_layout()
    _save_figure(
        fig, out_dir, "F6_bacteria_scatter",
        "Figure F6 \u2014 Bacterial Infection Prevalence: Simulation vs. Calibration Target",
        f"Each point represents one of 42 organisms. Dashed diagonal = perfect fit. "
        f"Orange dotted lines = \xd72/\xf70.5 tolerance. "
        f"Outliers (>2\xd7 or <0.5\xd7 target) labelled with abbreviated names. "
        f"{'Error bars show 5th-95th percentile ranges across accepted runs. ' if n > 1 else ''}"
        f"n\u2009=\u2009{n} run{'s' if n > 1 else ''}.",
        ["Infection prevalence is the percentage of the world population with an active "
         "infection from the specified organism on an average day during the "
         "2022\u20132025 calibration window.",
         "Figure version of Table T4."],
        agg=agg,
    )


# ---------------------------------------------------------------------------
# Figure F7 — Hospital vs. community resistance heatmap (figure version of T5)
# ---------------------------------------------------------------------------

# Legacy hospital/community heatmap; retained for reference, not called by main().
def make_f7_hc_resistance_heatmap(agg: dict, out_dir: Path) -> None:
    """
    Figure F7: heatmap of hospital vs. community acquisition and resistance rates
    (figure version of Table T5).
    """
    srl = agg.get("serious_resistance_locus", pd.DataFrame())
    ril = agg.get("resistance_incidence_locus", pd.DataFrame())
    bi  = agg.get("bacteria_infections",       pd.DataFrame())
    n   = agg.get("n_runs", 1)
    if srl is None or srl.empty:
        print("  F7: no serious_resistance_locus data — skipping.")
        return
    # Filter on the any-R structural benchmark (same as T5).
    first_col = srl.columns[0]
    summary_mask = srl[first_col].astype(str).str.match(
        r"^\s*(-|Resistance Locus|Serious Resistance|Mean |H:C)", na=False)
    srl = srl[~summary_mask].copy()
    ril_c = pd.DataFrame()
    if ril is not None and not ril.empty:
        ril_mask = ril[ril.columns[0]].astype(str).str.match(
            r"^\s*(-|Resistance Locus|Serious Resistance|Mean |H:C)", na=False
        )
        ril_c = ril[~ril_mask].copy()
    target_col = "Target H:C ratio"
    if target_col in ril_c.columns and "Bacteria" in ril_c.columns:
        ril_c = _add_interval_columns(ril_c, target_col, "_target_hc")
        ril_c[target_col] = ril_c["_target_hc_med"]
        included_names = set(ril_c.loc[ril_c[target_col] > 1.0, "Bacteria"].astype(str))
        srl = srl[srl["Bacteria"].astype(str).isin(included_names)].copy()
    # Build merged matrix
    keep_srl = ["Bacteria"] + [c for c in ["Hospital Serious-R (%)", "Community Serious-R (%)"]
                                if c in srl.columns]
    out = srl[keep_srl].copy()
    if not ril_c.empty:
        any_cols = ["Bacteria"] + [c for c in [
            "Hospital Infections with Any Resistance (%)",
            "Community Infections with Any Resistance (%)"] if c in ril_c.columns]
        out = out.merge(ril_c[any_cols], on="Bacteria", how="left")
    if bi is not None and not bi.empty and "Hospital Acquired (%)" in bi.columns:
        out = out.merge(bi[["Bacteria", "Hospital Acquired (%)"]], on="Bacteria", how="left")
    col_renames = {
        "Hospital Acquired (%)":                        "Hosp\nacquired",
        "Hospital Infections with Any Resistance (%)":  "Hosp\nany-R",
        "Community Infections with Any Resistance (%)": "Comm\nany-R",
        "Hospital Serious-R (%)":                       "Hosp\nserious-R",
        "Community Serious-R (%)":                      "Comm\nserious-R",
    }
    out = out.rename(columns=col_renames)
    value_cols = [v for v in col_renames.values() if v in out.columns]
    for c in value_cols:
        out = _add_interval_columns(out, c, f"_{c.replace(chr(10), '_').lower()}")
        out[c] = out[f"_{c.replace(chr(10), '_').lower()}_med"]
    matrix  = out[value_cols].values.astype(float)
    yticks  = out["Bacteria"].values
    fig_h   = max(5, 0.40 * len(yticks))
    fig_w   = max(6, 1.85 * len(value_cols))
    fig, ax = plt.subplots(figsize=(fig_w, fig_h))
    norm    = mcolors.Normalize(vmin=0, vmax=100)
    im      = ax.imshow(matrix, cmap=plt.cm.YlOrRd, norm=norm, aspect="auto")
    ax.set_xticks(range(len(value_cols)))
    ax.set_xticklabels(value_cols, fontsize=9)
    ax.set_yticks(range(len(yticks)))
    ax.set_yticklabels(yticks, fontsize=7.5, fontstyle="italic")
    ax.tick_params(top=True, bottom=False, labeltop=True, labelbottom=False)
    for i in range(len(yticks)):
        for j in range(len(value_cols)):
            v = matrix[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.0f}", ha="center", va="center",
                        fontsize=6.5, color="white" if v > 55 else "black")
    cbar = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02)
    cbar.set_label("Percentage (%)", fontsize=8)
    fig.suptitle(
        "Figure F7 \u2014 Hospital vs. community resistance and acquisition rates",
        fontsize=10, fontweight="bold", y=1.02)
    fig.tight_layout()
    _save_figure(
        fig, out_dir, "F7_hc_resistance_heatmap",
        "Figure F7 \u2014 Hospital vs. Community Resistance and Acquisition Rates",
        "Figure version of Table T5. Colour scale: 0\u2013100%. "
        "Only organisms with an expert-assigned any-R structural H:C benchmark > 1.0 are shown; "
        "the benchmark is not a marker-drug serious-R target.",
        ["Hosp/Comm any-R: percentage of new hospital/community-acquired infections "
         "carrying any resistance mechanism.",
         "Hosp/Comm serious-R: percentage with resistance to the marker drug for that organism "
         "(e.g. meropenem for Gram-negatives, flucloxacillin for S. aureus).",
         "Hosp acquired: percentage of new infections acquired during hospitalisation.",
         "Figure version of Table T5."],
        agg=agg,
    )


# ---------------------------------------------------------------------------
# Supplementary Figure FS1 — Infection deaths per organism (figure version of S1)
# ---------------------------------------------------------------------------

def make_figure_4_calibration_infection_deaths(agg: dict, out_dir: Path) -> None:
    """
    Figure 4: horizontal bar chart of infection deaths per bacterium,
    simulation vs. target/observed estimate.
    """
    bm = agg.get("bacteria_mortality", pd.DataFrame()).copy()
    n  = agg.get("n_runs", 1)
    if bm is None or bm.empty:
        print("  FS1: no bacteria_mortality data — skipping.")
        return
    bact_col = bm.columns[0]
    tgt_col  = next((c for c in bm.columns
                     if "target" in c.lower() and "death" in c.lower()), None)
    sim_col  = next((c for c in bm.columns
                     if "simulation" in c.lower() and "death" in c.lower()), None)
    if sim_col is None:
        return
    bm = bm[
        ~bm[bact_col].apply(
            lambda value: _bacteria_slug_for_filter(value) in _INFECTION_DEATH_EXCLUDED_BACTERIA_SLUGS
        )
    ].copy()
    if bm.empty:
        print("  Figure 4: no non-excluded bacteria mortality data — skipping.")
        return
    if tgt_col:
        bm = _add_interval_columns(bm, tgt_col, "_tgt")
        bm["_tgt"] = bm["_tgt_med"]
    else:
        bm["_tgt"] = np.nan
    bm = _add_interval_columns(bm, sim_col, "_sim")
    bm["_sim"] = bm["_sim_med"]
    bm = bm.dropna(subset=["_sim"])
    bm = bm[(bm["_tgt"].fillna(0) > 0) | (bm["_sim"].fillna(0) > 0)].copy()
    sort_key = "_tgt" if bm["_tgt"].notna().any() else "_sim"
    bm = bm.sort_values(sort_key, ascending=True, na_position="first")
    n_orgs = len(bm)
    fig, ax = plt.subplots(figsize=(9, max(4, 0.42 * n_orgs)))
    y = np.arange(n_orgs)
    bar_h = 0.38
    sim_err_lo, sim_err_hi = _asymmetric_errors(bm, "_sim")
    ax.barh(
        y + bar_h / 2,
        bm["_sim"].fillna(0),
        bar_h,
        color="#2196F3",
        alpha=0.85,
        label="Simulation",
        xerr=[sim_err_lo, sim_err_hi] if n > 1 else None,
        error_kw={"elinewidth": 0.9, "ecolor": "#0D47A1", "capthick": 0.9, "capsize": 2.5},
    )
    if bm["_tgt"].notna().any():
        ax.barh(y - bar_h / 2, bm["_tgt"].fillna(0), bar_h,
                color="#FF7043", alpha=0.85, label="Observed estimate")
    ax.set_yticks(y)
    ax.set_yticklabels(bm[bact_col].values, fontsize=7.5, fontstyle="italic")
    ax.set_xlabel("Deaths (millions per year)", fontsize=10)
    ax.legend(fontsize=9, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", linewidth=0.4, alpha=0.5)
    fig.suptitle("Figure 4. Calibration: 2025 infection deaths by bacterium", fontsize=10, fontweight="bold")
    fig.tight_layout()
    _save_figure(
        fig, out_dir, "Figure_4__calibration_infection_deaths_by_bacteria",
        "Figure 4. Calibration: 2025 infection deaths by bacterium",
        f"Bacteria with zero modelled deaths are omitted. Sorted by target/observed estimate "
        f"where available, otherwise by simulation. "
        f"{'Error bars show simulation 5th-95th percentile ranges. ' if n > 1 else ''}"
        f"n\u2009=\u2009{n} run{'s' if n > 1 else ''}.",
        ["Deaths are scaled to the global population using the run-specific population scale factor "
         "and annualised to yearly equivalents."] + _MORTALITY_TARGET_SOURCE_NOTES,
        agg=agg,
    )


# ---------------------------------------------------------------------------
# Supplementary Figure FS2 — Syndrome incidence (figure version of S2)
# ---------------------------------------------------------------------------

def make_figure_5_calibration_carriage_prevalence(agg: dict, out_dir: Path) -> None:
    """
    Figure 5: horizontal bar chart of asymptomatic carriage prevalence by bacterium,
    simulation vs. target/observed estimate where available.
    """
    bi_raw = agg.get("bacteria_infections", pd.DataFrame()).copy()
    n = agg.get("n_runs", 1)
    if bi_raw is None or bi_raw.empty:
        print("  Figure 5: no bacteria_infections data — skipping.")
        return

    bi = _clean_df(bi_raw)
    bact_col = bi.columns[0]
    sim_col = next(
        (c for c in bi.columns if "carriage" in c.lower() and "simulation" in c.lower() and "%" in c),
        None,
    )
    tgt_col = next(
        (
            c for c in bi.columns
            if "carriage" in c.lower()
            and ("observed estimate" in c.lower() or "target" in c.lower())
            and "%" in c
        ),
        None,
    )
    if sim_col is None:
        print("  Figure 5: no carriage simulation column found — skipping.")
        return

    bi = _add_interval_columns(bi, sim_col, "_sim")
    if tgt_col:
        bi = _add_interval_columns(bi, tgt_col, "_tgt")
    else:
        bi["_tgt_med"] = np.nan

    bi = bi.dropna(subset=["_sim_med"]).copy()
    if bi.empty:
        print("  Figure 5: no valid carriage simulation values — skipping.")
        return

    sort_col = "_tgt_med" if bi["_tgt_med"].notna().any() else "_sim_med"
    bi = bi.sort_values(sort_col, ascending=True, na_position="first").reset_index(drop=True)

    n_bacteria = len(bi)
    fig, ax = plt.subplots(figsize=(9, max(4, 0.42 * n_bacteria)))
    y = np.arange(n_bacteria)
    bar_h = 0.38
    sim_err_lo, sim_err_hi = _asymmetric_errors(bi, "_sim")

    ax.barh(
        y + bar_h / 2,
        bi["_sim_med"].fillna(0),
        bar_h,
        color="#2196F3",
        alpha=0.85,
        label="Simulation",
        xerr=[sim_err_lo, sim_err_hi] if n > 1 else None,
        error_kw={"elinewidth": 0.9, "ecolor": "#0D47A1", "capthick": 0.9, "capsize": 2.5},
    )
    if bi["_tgt_med"].notna().any():
        ax.barh(
            y - bar_h / 2,
            bi["_tgt_med"].fillna(0),
            bar_h,
            color="#FF7043",
            alpha=0.85,
            label="Target/observed estimate",
        )

    ax.set_yticks(y)
    ax.set_yticklabels(bi[bact_col].values, fontsize=7.5, fontstyle="italic")
    ax.set_xlabel("Carriage prevalence (% of world population)", fontsize=10)
    ax.legend(fontsize=9, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", linewidth=0.4, alpha=0.5)
    fig.suptitle("Figure 5. Calibration: 2025 prevalence of carriage by bacterium", fontsize=10, fontweight="bold")
    fig.tight_layout()

    _save_figure(
        fig,
        out_dir,
        "Figure_5__calibration_carriage_prevalence_by_bacteria",
        "Figure 5. Calibration: 2025 prevalence of carriage by bacterium",
        f"Horizontal paired bars show simulated asymptomatic carriage prevalence and "
        f"target/observed estimates where available. "
        f"{'Error bars show simulation 5th-95th percentile ranges. ' if n > 1 else ''}"
        f"n\u2009=\u2009{n} run{'s' if n > 1 else ''}.",
        ["Zero target/observed-estimate values are retained where they represent design zeros."]
        + _CARRIAGE_TARGET_SOURCE_NOTES,
        agg=agg,
    )


# Legacy simulation-only carriage builder retained for reference only; not called by main().
def make_legacy_carriage_prevalence_simulation_only(agg: dict, out_dir: Path) -> None:
    """
    Figure 6: simulation-only asymptomatic carriage prevalence by bacterium.
    """
    bi_raw = agg.get("bacteria_infections", pd.DataFrame()).copy()
    n = agg.get("n_runs", 1)
    if bi_raw is None or bi_raw.empty:
        print("  Figure 6: no bacteria_infections data - skipping.")
        return

    bi = _clean_df(bi_raw)
    bact_col = bi.columns[0]
    sim_col = next(
        (c for c in bi.columns if "carriage" in c.lower() and "simulation" in c.lower() and "%" in c),
        None,
    )
    if sim_col is None:
        print("  Figure 6: no carriage simulation column found - skipping.")
        return

    bi = _add_interval_columns(bi, sim_col, "_sim")
    bi = bi.dropna(subset=["_sim_med"]).copy()
    if bi.empty:
        print("  Figure 6: no valid carriage simulation values - skipping.")
        return

    bi = bi.sort_values("_sim_med", ascending=True, na_position="first").reset_index(drop=True)

    n_bacteria = len(bi)
    fig, ax = plt.subplots(figsize=(9, max(4, 0.42 * n_bacteria)))
    y = np.arange(n_bacteria)
    sim_err_lo, sim_err_hi = _asymmetric_errors(bi, "_sim")
    ax.barh(
        y,
        bi["_sim_med"].fillna(0),
        0.55,
        color="#2196F3",
        alpha=0.85,
        label="Simulation",
        xerr=[sim_err_lo, sim_err_hi] if n > 1 else None,
        error_kw={"elinewidth": 0.9, "ecolor": "#0D47A1", "capthick": 0.9, "capsize": 2.5},
    )
    ax.set_yticks(y)
    ax.set_yticklabels(bi[bact_col].values, fontsize=7.5, fontstyle="italic")
    ax.set_xlabel("Carriage prevalence (% of world population)", fontsize=10)
    ax.legend(fontsize=9, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", linewidth=0.4, alpha=0.5)
    fig.suptitle("Figure 6. Prevalence of carriage by bacterium, 2022-2025", fontsize=10, fontweight="bold")
    fig.tight_layout()

    _save_figure(
        fig,
        out_dir,
        "Figure_6__carriage_prevalence_by_bacteria",
        "Figure 6. Prevalence of carriage by bacterium, 2022-2025",
        f"Simulation output for the 2022-2025 calibration window; target carriage values "
        f"are not shown in this figure. "
        f"{'Error bars show simulation 5th-95th percentile ranges. ' if n > 1 else ''}"
        f"n = {n} run{'s' if n > 1 else ''}.",
        [
            "Carriage prevalence is the percentage of the world population carrying the bacterium "
            "asymptomatically in the modelled microbiome/carriage compartment."
        ],
        agg=agg,
    )


_REGION_LABELS = {
    "africa": "Africa",
    "asia": "Asia",
    "europe": "Europe",
    "north_america": "North America",
    "oceania": "Oceania",
    "south_america": "South America",
    "home": "Home",
}

_REGION_ORDER = (
    "africa",
    "asia",
    "europe",
    "north_america",
    "oceania",
    "south_america",
    "home",
)


def _figure_7_rows_from_simulation_csv(csv_path: Path) -> list[dict[str, object]]:
    """
    Compute all-age regional infection death rates from one simulation_summary CSV.

    The rate is deaths from sepsis plus infection_non_sepsis per 100,000 alive
    per year during calendar years 2022-2025. The denominator is regional
    person-years alive, not an unweighted mean of age-specific rates.
    """
    columns = _simulation_csv_columns(csv_path)
    if columns is None:
        return []

    if "time_in_years" not in columns:
        return []

    available_regions = [
        region for region in _REGION_ORDER
        if (
            f"{region}_population" in columns
            and f"{region}_deaths_sepsis" in columns
            and f"{region}_deaths_infection_non_sepsis" in columns
        )
    ]
    if not available_regions:
        return []

    usecols = ["time_in_years"]
    if "run_id" in columns:
        usecols.append("run_id")
    for region in available_regions:
        usecols.extend([
            f"{region}_population",
            f"{region}_deaths_sepsis",
            f"{region}_deaths_infection_non_sepsis",
        ])

    try:
        df = _read_csv_selected(csv_path, usecols)
    except (FileNotFoundError, ValueError, OSError):
        return []

    df["calendar_year"] = _F1_SIM_EPOCH_YEAR + pd.to_numeric(df["time_in_years"], errors="coerce")
    df = df[(df["calendar_year"] >= 2022.0) & (df["calendar_year"] < 2026.0)].copy()
    if df.empty:
        return []

    time_values = pd.to_numeric(df["time_in_years"], errors="coerce").dropna().sort_values()
    diffs = time_values.diff().dropna()
    diffs = diffs[diffs > 0]
    step_years = float(diffs.median()) if not diffs.empty else 1.0 / 365.0
    if not np.isfinite(step_years) or step_years <= 0:
        step_years = 1.0 / 365.0

    group_cols = ["run_id"] if "run_id" in df.columns else []
    grouped = df.groupby(group_cols, dropna=False) if group_cols else [(csv_path.stem, df)]

    rows: list[dict[str, object]] = []
    for run_key, run_df in grouped:
        for region in available_regions:
            pop_col = f"{region}_population"
            sepsis_col = f"{region}_deaths_sepsis"
            non_sepsis_col = f"{region}_deaths_infection_non_sepsis"
            population = pd.to_numeric(run_df[pop_col], errors="coerce")
            deaths = (
                pd.to_numeric(run_df[sepsis_col], errors="coerce").fillna(0)
                + pd.to_numeric(run_df[non_sepsis_col], errors="coerce").fillna(0)
            )
            person_years = float(population.fillna(0).sum() * step_years)
            infection_deaths = float(deaths.sum())
            if person_years <= 0 or not np.isfinite(person_years):
                continue
            rows.append({
                "source": csv_path.name,
                "run": str(run_key),
                "region": region,
                "region_label": _REGION_LABELS.get(region, region.replace("_", " ").title()),
                "rate": 100000.0 * infection_deaths / person_years,
                "infection_deaths": infection_deaths,
                "person_years": person_years,
            })
    return rows


def make_figure_7_infection_death_rate_by_region(csv_paths: list[Path], out_dir: Path, agg: dict | None = None) -> None:
    """
    Figure 8: all-age regional infection death rates from simulation_summary CSVs.
    """
    rows: list[dict[str, object]] = []
    for csv_path in csv_paths:
        rows.extend(_figure_7_rows_from_simulation_csv(csv_path))

    title = "Figure 8. Infection death rate by region, 2022-2025"
    stem = "Figure_8__infection_death_rate_by_region"

    if not rows:
        fig, ax = plt.subplots(figsize=(9, 4.8))
        ax.text(
            0.5,
            0.5,
            "Figure 8. Infection death rate by region\n\n"
            "All-age regional infection death rates require regional infection-death counts\n"
            "and regional population denominators for the calibration window.\n\n"
            "The available inputs did not contain the required numerator and denominator\n"
            "columns in matching simulation_summary_*.csv files.",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=10.5,
            color="#555",
            bbox=dict(boxstyle="round,pad=0.6", fc="#f5f5f5", ec="#bbb"),
        )
        ax.set_axis_off()
        fig.tight_layout()
        _save_figure(
            fig,
            out_dir,
            stem,
            title,
            "All-age regional infection death rates require regional infection-death counts and "
            "regional population denominators for the calibration window. The available calibration "
            "summary contains age-specific rates but not the denominators needed to combine them "
            "into all-age regional rates.",
            [],
            agg=agg,
        )
        return

    rate_df = pd.DataFrame(rows)
    summary = (
        rate_df.groupby(["region", "region_label"], as_index=False)["rate"]
        .agg(
            median="median",
            p5=lambda s: float(np.percentile(s, 5)),
            p95=lambda s: float(np.percentile(s, 95)),
            n="count",
        )
    )
    summary = summary.sort_values("median", ascending=True).reset_index(drop=True)
    n_runs = int(rate_df[["source", "run"]].drop_duplicates().shape[0])

    fig, ax = plt.subplots(figsize=(8.5, max(3.5, 0.55 * len(summary))))
    y = np.arange(len(summary))
    err_lo = np.clip(summary["median"].to_numpy(float) - summary["p5"].to_numpy(float), 0, None)
    err_hi = np.clip(summary["p95"].to_numpy(float) - summary["median"].to_numpy(float), 0, None)
    ax.barh(
        y,
        summary["median"].to_numpy(float),
        0.55,
        color="#546E7A",
        alpha=0.88,
        xerr=[err_lo, err_hi] if n_runs > 1 else None,
        error_kw={"elinewidth": 0.9, "ecolor": "#263238", "capthick": 0.9, "capsize": 3},
    )
    ax.set_yticks(y)
    ax.set_yticklabels(summary["region_label"].values, fontsize=9)
    ax.set_xlabel("Infection deaths per 100,000 alive per year", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", linewidth=0.4, alpha=0.5)
    fig.suptitle(title, fontsize=10.5, fontweight="bold")
    fig.tight_layout()

    _save_figure(
        fig,
        out_dir,
        stem,
        title,
        "Deaths include sepsis and infection_non_sepsis deaths. Rates are all-age rates "
        "for the 2022-2025 calibration window. "
        f"{'Error bars show 5th-95th percentile ranges across runs. ' if n_runs > 1 else ''}"
        f"n = {n_runs} simulation run{'s' if n_runs > 1 else ''}.",
        [
            "Source: matched simulation_summary_*.csv files discovered from the supplied "
            "calibration_summary_*.txt paths.",
            "Calculation: 100,000 x regional infection deaths divided by regional person-years "
            "alive. Person-years are summed from regional alive population counts across "
            "daily/timestep simulation rows.",
        ],
        agg=agg,
    )


_F8_OLD_CONTEXT_COLUMNS: list[tuple[str, str]] = [
    ("Empiric", "currently_taking_drug_count_empiric"),
    ("Targeted", "currently_taking_drug_count_targeted"),
    ("Prophylaxis", "currently_taking_drug_count_prophylaxis"),
    ("Other", "currently_taking_drug_count_other"),
]

_F8_DETAILED_CONTEXT_COLUMNS: list[tuple[str, str]] = [
    ("Empiric", "currently_taking_drug_count_empiric"),
    ("Targeted", "currently_taking_drug_count_targeted"),
    ("Prophylaxis", "currently_taking_drug_count_prophylaxis"),
    (
        "Other: active asymptomatic infection",
        "currently_taking_drug_count_other_active_asymptomatic_modelled_bacterial_infection",
    ),
    (
        "Other: no active infection",
        "currently_taking_drug_count_other_no_active_modelled_infection",
    ),
    ("Unknown / legacy", "currently_taking_drug_count_other_unknown_or_legacy"),
]

_F8_CONTEXT_COLOURS: dict[str, str] = {
    "Empiric": "#4E79A7",
    "Targeted": "#59A14F",
    "Prophylaxis": "#F28E2B",
    "Other": "#9C755F",
    "Other: active asymptomatic infection": "#B07AA1",
    "Other: no active infection": "#9C755F",
    "Unknown / legacy": "#BAB0AC",
}


def _figure_8_placeholder(out_dir: Path, agg: dict | None, message: str) -> None:
    title = "Figure 9. Antibiotic use by treatment context, 2025"
    stem = "Figure_9__antibiotic_use_by_treatment_context"
    fig, ax = plt.subplots(figsize=(9, 3.2))
    ax.text(
        0.5,
        0.5,
        f"{title}\n\n{message}",
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=10.5,
        color="#555",
        bbox=dict(boxstyle="round,pad=0.6", fc="#f5f5f5", ec="#bbb"),
    )
    ax.set_axis_off()
    fig.subplots_adjust(left=0.03, right=0.97, top=0.92, bottom=0.08)
    _save_figure(fig, out_dir, stem, title, message, [], agg=agg)


def _simulation_year_series(df: pd.DataFrame) -> pd.Series:
    for col in ("simulation_year", "year"):
        if col in df.columns:
            return pd.to_numeric(df[col], errors="coerce")
    if "time_in_years" in df.columns:
        return _F1_SIM_EPOCH_YEAR + pd.to_numeric(df["time_in_years"], errors="coerce")
    if "time_step" in df.columns:
        return _F1_SIM_EPOCH_YEAR + pd.to_numeric(df["time_step"], errors="coerce") / 365.0
    return pd.Series(np.nan, index=df.index)


def _figure_8_rows_from_simulation_csv(
    csv_path: Path,
    scale_factor: float | None,
) -> tuple[list[dict[str, object]], str | None, str]:
    old_required = [column for _, column in _F8_OLD_CONTEXT_COLUMNS]
    detailed_required = [column for _, column in _F8_DETAILED_CONTEXT_COLUMNS]
    optional = ["policy_option", "run_id", "simulation_year", "year", "time_in_years", "time_step"]
    wanted = set(old_required + detailed_required + optional)
    try:
        df = _read_csv_selected(csv_path, wanted)
    except (FileNotFoundError, ValueError, OSError) as exc:
        return [], f"{csv_path.name}: could not read simulation CSV ({exc}).", "missing"

    has_detailed = all(column in df.columns for column in detailed_required)
    has_old = all(column in df.columns for column in old_required)
    if has_detailed:
        context_columns = _F8_DETAILED_CONTEXT_COLUMNS
        mode = "detailed"
    elif has_old:
        context_columns = _F8_OLD_CONTEXT_COLUMNS
        mode = "old"
    else:
        missing = [column for column in old_required if column not in df.columns]
        return [], f"{csv_path.name}: missing columns {', '.join(missing)}.", "missing"

    if scale_factor is None or not np.isfinite(scale_factor) or scale_factor <= 0.0:
        return [], f"{csv_path.name}: missing population scale factor from calibration summary.", mode

    if "policy_option" in df.columns:
        policy = pd.to_numeric(df["policy_option"], errors="coerce")
        df = df[policy == 0].copy()

    df["simulation_year_for_f8"] = _simulation_year_series(df)
    df = df[(df["simulation_year_for_f8"] >= 2025.0) & (df["simulation_year_for_f8"] < 2026.0)].copy()
    if df.empty:
        return [], f"{csv_path.name}: no baseline 2025 rows available.", mode

    for _, column in context_columns:
        df[column] = pd.to_numeric(df[column], errors="coerce")

    grouped = df.groupby("run_id", dropna=False) if "run_id" in df.columns else [(csv_path.stem, df)]
    rows: list[dict[str, object]] = []
    for run_key, run_df in grouped:
        row: dict[str, object] = {"source": csv_path.name, "run": str(run_key)}
        has_value = False
        for label, column in context_columns:
            mean_count = float(run_df[column].mean(skipna=True))
            if np.isfinite(mean_count):
                row[label] = mean_count * scale_factor / 1_000_000.0
                has_value = True
            else:
                row[label] = np.nan
        if has_value:
            rows.append(row)
    if not rows:
        return [], f"{csv_path.name}: context columns were present but contained no usable 2025 values.", mode
    return rows, None, mode


def make_figure_8_antibiotic_use_by_context(
    csv_runs: list[tuple[Path, float | None]],
    out_dir: Path,
    agg: dict | None = None,
) -> None:
    title = "Figure 9. Antibiotic use by treatment context, 2025"
    stem = "Figure_9__antibiotic_use_by_treatment_context"
    missing_columns_message = (
        "Figure 9 requires simulation_summary CSV columns for antibiotic-use context. "
        "Re-run the Rust simulation after adding currently_taking_drug_count_empiric, "
        "currently_taking_drug_count_targeted, currently_taking_drug_count_prophylaxis, "
        "and currently_taking_drug_count_other."
    )

    if not csv_runs:
        _figure_8_placeholder(
            out_dir,
            agg,
            "Figure 9 requires matching simulation_summary_*.csv files with antibiotic-use context columns.",
        )
        return

    rows: list[dict[str, object]] = []
    detailed_rows: list[dict[str, object]] = []
    old_rows: list[dict[str, object]] = []
    problems: list[str] = []
    saw_missing_required_columns = False
    for csv_path, scale_factor in csv_runs:
        run_rows, problem, mode = _figure_8_rows_from_simulation_csv(csv_path, scale_factor)
        if mode == "detailed":
            detailed_rows.extend(run_rows)
        elif mode == "old":
            old_rows.extend(run_rows)
        if problem:
            problems.append(problem)
            if "missing columns" in problem:
                saw_missing_required_columns = True

    mode_note = ""
    context_columns = _F8_DETAILED_CONTEXT_COLUMNS
    if detailed_rows:
        rows = detailed_rows
        context_columns = _F8_DETAILED_CONTEXT_COLUMNS
        if old_rows:
            mode_note = (
                "Some supplied simulation CSVs predate the detailed Other split and were not "
                "included in the detailed Figure 9 summary."
            )
    elif old_rows:
        rows = old_rows
        context_columns = _F8_OLD_CONTEXT_COLUMNS
        mode_note = (
            "This run predates the detailed Other split; rerun the Rust simulation to split "
            "Other into no-active-infection, active-asymptomatic-infection, and unknown/legacy "
            "components."
        )

    if not rows:
        if saw_missing_required_columns:
            _figure_8_placeholder(out_dir, agg, missing_columns_message)
        else:
            detail = " ".join(problems) if problems else "No usable 2025 baseline context rows were found."
            _figure_8_placeholder(out_dir, agg, detail)
        return

    df = pd.DataFrame(rows)
    labels = [label for label, _ in context_columns]
    medians = [float(df[label].median(skipna=True)) for label in labels]
    p5s = [float(np.nanpercentile(df[label], 5)) for label in labels]
    p95s = [float(np.nanpercentile(df[label], 95)) for label in labels]
    total_median = float(np.nansum(medians))
    n_runs = int(df[["source", "run"]].drop_duplicates().shape[0])

    fig, ax = plt.subplots(figsize=(9, 3.2))
    left = 0.0
    for label, value in zip(labels, medians):
        safe_value = 0.0 if not np.isfinite(value) else value
        if label == "Unknown / legacy" and safe_value <= 0.0:
            continue
        ax.barh(
            [0],
            [safe_value],
            left=left,
            height=0.46,
            color=_F8_CONTEXT_COLOURS[label],
            label=label,
        )
        if safe_value > 0.0 and total_median > 0.0 and safe_value / total_median >= 0.08:
            pct = 100.0 * safe_value / total_median
            ax.text(
                left + safe_value / 2.0,
                0,
                f"{label}\n{safe_value:.2f}M\n{pct:.0f}%",
                ha="center",
                va="center",
                fontsize=8,
                color="white",
                fontweight="bold",
            )
        left += safe_value

    ax.set_yticks([])
    ax.set_xlabel("People on antibiotics on an average day in 2025 (millions)", fontsize=10)
    ax.set_xlim(0, max(total_median * 1.08, 0.1))
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", linewidth=0.4, alpha=0.5)
    ax.legend(loc="lower center", bbox_to_anchor=(0.5, -0.36), ncol=4, frameon=False, fontsize=9)
    fig.suptitle(title, fontsize=10.5, fontweight="bold")
    fig.tight_layout(rect=[0, 0.05, 1, 1])

    table_rows = []
    for label, value in zip(labels, medians):
        safe_value = 0.0 if not np.isfinite(value) else value
        table_rows.append({
            "Context": label,
            "People on antibiotics (millions)": f"{safe_value:.3f}",
            "% of total": f"{(100.0 * safe_value / total_median):.1f}%" if total_median > 0.0 else "0.0%",
        })
    summary_table_html = "<h2>Figure 9 Summary</h2>\n" + _html_table(pd.DataFrame(table_rows))

    interval_note = ""
    if n_runs > 1:
        intervals = [
            f"{label}: {lo:.2f}-{hi:.2f}M"
            for label, lo, hi in zip(labels, p5s, p95s)
            if np.isfinite(lo) and np.isfinite(hi)
        ]
        interval_note = " Per-context 5th-95th percentile ranges: " + "; ".join(intervals) + "."

    footnotes = [
        "Values are mean daily people on antibiotics during baseline-policy rows in calendar "
        "year 2025, scaled to millions using the population scale factor from the matching "
        "calibration_summary file.",
        f"Stacked segments show medians across {n_runs} simulation run"
        f"{'s' if n_runs > 1 else ''}.{interval_note}",
    ]
    if mode_note:
        footnotes.append(mode_note)

    _save_figure(
        fig,
        out_dir,
        stem,
        title,
        "Context is assigned when each antibiotic course starts and is retained until that "
        "course stops. People taking multiple antibiotics are assigned to one category using "
        "priority: targeted, empiric, prophylaxis, other active asymptomatic modelled bacterial "
        "infection, other no active modelled infection, then unknown/legacy. 'Other: no active "
        "modelled infection' is a proxy for non-bacterial, non-modelled, or background prescribing, "
        "not a direct viral diagnosis. 'Unknown/legacy' indicates missing or legacy context labels "
        "and should be near zero in new runs.",
        footnotes,
        agg=agg,
        extra_html=summary_table_html,
    )


_F10_TITLE = "Figure 13. Counterfactual resistance-acquisition pathway comparisons"
_F10_STEM = "Figure_13__resistance_pathway_counterfactuals"
_F10_NOTE = (
    "Figure 13 is intentionally shown as a placeholder. The planned analysis requires a set "
    "of counterfactual model runs in which resistance-acquisition or resistance-spread pathways "
    "are individually disabled or modified. The current baseline simulation_summary and "
    "calibration_summary files do not contain the required scenario comparisons."
)


def make_figure_10_resistance_pathway_counterfactuals(
    out_dir: Path,
    agg: dict | None = None,
) -> None:
    placeholder_text = (
        "Figure 13 placeholder\n\n"
        "This figure will compare resistance and health outcomes across counterfactual\n"
        "model runs in which individual resistance-acquisition or resistance-spread\n"
        "pathways are disabled or modified.\n\n"
        "Required future inputs:\n"
        "- baseline run with all pathways active\n"
        "- pathway-ablation runs\n"
        "- consistent run settings and calibration window\n"
        "- agreed outcome metrics, such as any-R, serious-R,\n"
        "  resistance-adjusted activity, antibiotic use, sepsis, and infection deaths\n\n"
        "This figure will be completed after the counterfactual run set is defined\n"
        "and generated."
    )

    fig, ax = plt.subplots(figsize=(10.5, 5.2))
    ax.text(
        0.5,
        0.5,
        placeholder_text,
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=11,
        color="#333",
        linespacing=1.35,
        bbox=dict(boxstyle="round,pad=0.7", fc="#f7f8fa", ec="#b7c0cc"),
    )
    ax.set_axis_off()
    fig.suptitle(_F10_TITLE, fontsize=11.5, fontweight="bold")
    fig.subplots_adjust(left=0.03, right=0.97, top=0.88, bottom=0.08)

    planned_scenarios = pd.DataFrame([
        {
            "Planned scenario": "All pathways active",
            "Resistance pathway affected": "Reference scenario",
            "Status": "Not yet generated for this figure",
        },
        {
            "Planned scenario": "No infection de novo resistance",
            "Resistance pathway affected": "Within-infection emergence",
            "Status": "Pending",
        },
        {
            "Planned scenario": "No microbiome de novo resistance",
            "Resistance pathway affected": "Carriage/bystander selection",
            "Status": "Pending",
        },
        {
            "Planned scenario": "No HGT",
            "Resistance pathway affected": "Horizontal gene transfer",
            "Status": "Pending",
        },
        {
            "Planned scenario": "No carrier-to-infection inheritance",
            "Resistance pathway affected": "Carriage-to-infection transfer",
            "Status": "Pending",
        },
        {
            "Planned scenario": "No hospital enrichment",
            "Resistance pathway affected": "Hospital resistance-profile sampling",
            "Status": "Pending",
        },
        {
            "Planned scenario": "No environmental/static floors",
            "Resistance pathway affected": "Background reseeding",
            "Status": "Pending",
        },
        {
            "Planned scenario": "No persistence/ratchet floors",
            "Resistance pathway affected": "Long-term persistence",
            "Status": "Pending",
        },
        {
            "Planned scenario": "No late de novo resistance after drug maturity window",
            "Resistance pathway affected": "Post-introduction emergence",
            "Status": "Pending",
        },
    ])
    extra_html = "<h2>Planned Scenarios</h2>\n" + _html_table(planned_scenarios)
    footnotes = [
        "The final Figure 13 is intended to compare resistance and health outcomes across "
        "dedicated counterfactual model runs with matched random seeds and run settings.",
        "No counterfactual results are displayed here, and no pathway-ablation settings are "
        "inferred from current calibration summaries or baseline simulation summaries.",
    ]
    _save_figure(
        fig,
        out_dir,
        _F10_STEM,
        _F10_TITLE,
        _F10_NOTE,
        footnotes,
        agg=agg,
        extra_html=extra_html,
    )
    print("  Figure 13: placeholder for future resistance-pathway counterfactual runs.")


_F11_TITLE = "Figure 10. Underlying sepsis onset context and time to effective therapy, 2022\u20132025"
_F11_STEM = "Figure_10__sepsis_context_effective_therapy"
_F11_REQUIRED_MESSAGE = (
    "Figure 10 requires simulation_summary columns for sepsis-onset context counts and "
    "effective-therapy delay buckets. Re-run the Rust simulation after adding these aggregate "
    "columns to generate this figure."
)
_F11_CONTEXT_COLUMNS = [
    ("Targeted, effective", "sepsis_onset_targeted_effective_count"),
    ("Targeted, not effective", "sepsis_onset_targeted_not_effective_count"),
    ("Empiric, effective", "sepsis_onset_empiric_effective_count"),
    ("Empiric, not effective", "sepsis_onset_empiric_not_effective_count"),
    ("Other / prophylaxis only", "sepsis_onset_other_or_prophylaxis_only_count"),
    ("No antibiotic active", "sepsis_onset_no_antibiotic_count"),
    ("Unknown / legacy", "sepsis_onset_unknown_legacy_count"),
]
_F11_DELAY_EFFECTIVE_COLUMNS = [
    (
        "Effective on or before sepsis onset",
        "sepsis_effective_therapy_on_or_before_onset_count",
    ),
    (
        "Effective later in the same model day",
        "sepsis_effective_therapy_later_same_day_count",
    ),
    ("1 day after onset", "sepsis_effective_therapy_1_day_count"),
    ("2-3 days after onset", "sepsis_effective_therapy_2_3_days_count"),
    ("4+ days after onset", "sepsis_effective_therapy_4plus_days_count"),
]
_F11_NO_EFFECTIVE_COMBINED_COLUMN = (
    "No effective therapy before resolution, death, or censoring",
    "sepsis_no_effective_therapy_before_resolution_death_or_censoring_count",
)
_F11_NO_EFFECTIVE_SPLIT_COLUMNS = [
    (
        "No effective therapy before recovery",
        "sepsis_no_effective_therapy_before_recovery_count",
    ),
    (
        "No effective therapy before death",
        "sepsis_no_effective_therapy_before_death_count",
    ),
    (
        "No effective therapy before censoring/end of episode",
        "sepsis_no_effective_therapy_before_censoring_count",
    ),
]
_F11_NO_EFFECTIVE_UNKNOWN_COLUMN = (
    "No effective therapy outcome unknown",
    "sepsis_no_effective_therapy_unknown_count",
)
_F11_DELAY_UNKNOWN_COLUMN = (
    "Unknown / censored",
    "sepsis_effective_therapy_unknown_or_censored_count",
)
_F11_REQUIRED_DELAY_COLUMNS = (
    _F11_DELAY_EFFECTIVE_COLUMNS
    + [_F11_NO_EFFECTIVE_COMBINED_COLUMN, _F11_DELAY_UNKNOWN_COLUMN]
)
_F11_REQUIRED_COLUMNS = (
    ["time_in_years", "policy_option"]
    + [column for _, column in _F11_CONTEXT_COLUMNS]
    + [column for _, column in _F11_REQUIRED_DELAY_COLUMNS]
)
_F11_OPTIONAL_COLUMNS = (
    ["run_id", "time_step"]
    + [column for _, column in _F11_NO_EFFECTIVE_SPLIT_COLUMNS]
    + [_F11_NO_EFFECTIVE_UNKNOWN_COLUMN[1]]
)
_F11_CONTEXT_COLOURS = {
    "Targeted, effective": "#2A9D8F",
    "Targeted, not effective": "#D1495B",
    "Empiric, effective": "#4C78A8",
    "Empiric, not effective": "#F28E2B",
    "Other / prophylaxis only": "#8D99AE",
    "No antibiotic active": "#5C677D",
    "Unknown / legacy": "#6D597A",
}
_F11_DELAY_COLOURS = {
    "Effective on or before sepsis onset": "#2A9D8F",
    "Effective later in the same model day": "#3A86FF",
    "1 day after onset": "#59A14F",
    "2-3 days after onset": "#EDC948",
    "4+ days after onset": "#F28E2B",
    "No effective therapy before resolution, death, or censoring": "#D1495B",
    "No effective therapy before recovery": "#B56576",
    "No effective therapy before death": "#D1495B",
    "No effective therapy before censoring/end of episode": "#8D99AE",
    "No effective therapy outcome unknown": "#7A7A7A",
    "Unknown / censored": "#8D99AE",
}


def _figure_11_placeholder(out_dir: Path, agg: dict | None, message: str) -> None:
    fig, ax = plt.subplots(figsize=(9.5, 3.4))
    ax.text(
        0.5,
        0.5,
        f"{_F11_TITLE}\n\n{message}",
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=10.5,
        color="#555",
        bbox=dict(boxstyle="round,pad=0.6", fc="#f5f5f5", ec="#bbb"),
    )
    ax.set_axis_off()
    fig.subplots_adjust(left=0.03, right=0.97, top=0.92, bottom=0.08)
    _save_figure(fig, out_dir, _F11_STEM, _F11_TITLE, message, [], agg=agg)


def _figure_11_load_summary(csv_path: Path) -> tuple[pd.DataFrame | None, str | None]:
    available = _simulation_csv_columns(csv_path)
    if available is None:
        return None, f"{csv_path.name}: unable to read simulation summary CSV header."

    missing = [column for column in _F11_REQUIRED_COLUMNS if column not in available]
    if missing:
        return None, f"{csv_path.name}: missing columns {', '.join(missing)}."
    split_available = all(
        column in available
        for _, column in _F11_NO_EFFECTIVE_SPLIT_COLUMNS + [_F11_NO_EFFECTIVE_UNKNOWN_COLUMN]
    )

    usecols = [
        column
        for column in _F11_REQUIRED_COLUMNS + _F11_OPTIONAL_COLUMNS
        if column in available
    ]
    try:
        df = _read_csv_selected(csv_path, usecols)
    except (ValueError, OSError) as exc:
        return None, f"{csv_path.name}: unable to load summary rows ({exc})."

    df["source_file"] = csv_path.name
    df["figure11_no_effective_split_available"] = split_available
    return df, None


def _figure_11_counts_from_columns(
    df: pd.DataFrame,
    column_specs: list[tuple[str, str]],
) -> list[tuple[str, int]]:
    counts: list[tuple[str, int]] = []
    for label, column in column_specs:
        values = pd.to_numeric(df[column], errors="coerce").fillna(0)
        counts.append((label, int(values.sum())))
    return counts


def _figure_11_delay_counts(df: pd.DataFrame) -> tuple[list[tuple[str, int]], bool, bool]:
    counts = _figure_11_counts_from_columns(df, _F11_DELAY_EFFECTIVE_COLUMNS)
    split_flag = df.get("figure11_no_effective_split_available")
    if split_flag is None:
        split_mask = pd.Series(False, index=df.index)
    else:
        split_mask = split_flag.fillna(False).astype(bool)

    split_rows = df[split_mask]
    legacy_rows = df[~split_mask]
    split_available = not split_rows.empty
    legacy_combined_present = not legacy_rows.empty

    if split_available:
        counts.extend(_figure_11_counts_from_columns(split_rows, _F11_NO_EFFECTIVE_SPLIT_COLUMNS))
        unknown_count = _figure_11_counts_from_columns(
            split_rows,
            [_F11_NO_EFFECTIVE_UNKNOWN_COLUMN],
        )[0][1]
        if unknown_count:
            counts.append((_F11_NO_EFFECTIVE_UNKNOWN_COLUMN[0], unknown_count))
        if legacy_combined_present:
            counts.extend(
                _figure_11_counts_from_columns(legacy_rows, [_F11_NO_EFFECTIVE_COMBINED_COLUMN])
            )
    else:
        counts.extend(_figure_11_counts_from_columns(df, [_F11_NO_EFFECTIVE_COMBINED_COLUMN]))

    unknown_count = _figure_11_counts_from_columns(df, [_F11_DELAY_UNKNOWN_COLUMN])[0][1]
    if unknown_count:
        counts.append((_F11_DELAY_UNKNOWN_COLUMN[0], unknown_count))

    return counts, split_available, legacy_combined_present


def _figure_11_stacked_bar(
    ax,
    counts: list[tuple[str, int]],
    total: int,
    colours: dict[str, str],
    title: str,
) -> None:
    left = 0.0
    handles = []
    for label, count in counts:
        pct = 100.0 * count / total if total else 0.0
        if pct <= 0.0:
            continue
        color = colours.get(label, "#777777")
        ax.barh([0], [pct], left=left, height=0.48, color=color, edgecolor="white", linewidth=0.8)
        if pct >= 7.0:
            ax.text(
                left + pct / 2.0,
                0,
                f"{pct:.0f}%",
                ha="center",
                va="center",
                fontsize=8,
                color="white",
                fontweight="bold",
            )
        handles.append(mpatches.Patch(color=color, label=f"{label} ({pct:.1f}%)"))
        left += pct

    ax.set_xlim(0, 100)
    ax.set_yticks([])
    ax.set_xlabel("Incident sepsis episodes (%)", fontsize=9.5)
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.grid(axis="x", linewidth=0.35, alpha=0.45)
    if handles:
        ax.legend(
            handles=handles,
            loc="center left",
            bbox_to_anchor=(1.01, 0.5),
            frameon=False,
            fontsize=8,
        )


def make_figure_11_sepsis_context_effective_therapy(
    csv_paths: list[Path],
    out_dir: Path,
    agg: dict | None = None,
) -> None:
    if not csv_paths:
        _figure_11_placeholder(out_dir, agg, _F11_REQUIRED_MESSAGE)
        print("  Figure 10: placeholder (no simulation summary CSVs found).")
        return

    frames: list[pd.DataFrame] = []
    problems: list[str] = []
    for csv_path in csv_paths:
        frame, problem = _figure_11_load_summary(csv_path)
        if problem:
            problems.append(problem)
        if frame is not None and not frame.empty:
            frames.append(frame)

    if not frames:
        _figure_11_placeholder(out_dir, agg, _F11_REQUIRED_MESSAGE)
        if problems:
            print("  Figure 10: placeholder; " + " ".join(problems[:3]))
        else:
            print("  Figure 10: placeholder (simulation summary CSVs contained no rows).")
        return

    df = pd.concat(frames, ignore_index=True)
    df["policy_option"] = pd.to_numeric(df["policy_option"], errors="coerce")
    df["sepsis_onset_year"] = _F1_SIM_EPOCH_YEAR + pd.to_numeric(
        df["time_in_years"],
        errors="coerce",
    )
    df = df[
        (df["policy_option"] == 0)
        & (df["sepsis_onset_year"] >= 2022.0)
        & (df["sepsis_onset_year"] < 2026.0)
    ].copy()
    if df.empty:
        _figure_11_placeholder(
            out_dir,
            agg,
            "No baseline-policy incident sepsis episodes with onset years 2022-2025 were found.",
        )
        print("  Figure 10: placeholder (no baseline 2022-2025 sepsis episode rows).")
        return

    context_counts_all = _figure_11_counts_from_columns(df, _F11_CONTEXT_COLUMNS)
    delay_counts_all, no_effective_split_available, legacy_combined_present = (
        _figure_11_delay_counts(df)
    )
    context_total = sum(count for _, count in context_counts_all)
    delay_total = sum(count for _, count in delay_counts_all)

    if context_total == 0 or delay_total == 0:
        _figure_11_placeholder(
            out_dir,
            agg,
            "No baseline-policy sepsis-onset aggregate counts were found for onset years 2022-2025.",
        )
        print("  Figure 10: placeholder (no baseline 2022-2025 aggregate sepsis counts).")
        return

    context_counts = [(label, count) for label, count in context_counts_all if count]
    delay_counts = [(label, count) for label, count in delay_counts_all if count]
    threshold_note = "Effective therapy threshold: activity_r >= 0.500."

    fig, axes = plt.subplots(2, 1, figsize=(11.0, 4.9), sharex=True)
    _figure_11_stacked_bar(
        axes[0],
        context_counts,
        context_total,
        _F11_CONTEXT_COLOURS,
        "A. Antibiotic context at underlying sepsis onset",
    )
    _figure_11_stacked_bar(
        axes[1],
        delay_counts,
        delay_total,
        _F11_DELAY_COLOURS,
        "B. Time from underlying sepsis onset to first effective therapy",
    )
    fig.suptitle(_F11_TITLE, fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0.02, 0.78, 0.96])

    detail_rows: list[dict[str, object]] = []
    for panel, counts, denominator in [
        ("A. Antibiotic context at underlying sepsis onset", context_counts_all, context_total),
        ("B. Time from underlying sepsis onset to first effective therapy", delay_counts_all, delay_total),
    ]:
        for label, count in counts:
            pct = 100.0 * count / denominator if denominator else 0.0
            note = threshold_note if panel.startswith("B.") or "effective" in label.lower() else ""
            detail_rows.append({
                "Panel": panel,
                "Category": label,
                "Count": f"{count:,}",
                "Percent of panel denominator": f"{pct:.1f}%",
                "Notes": note,
            })

    details = pd.DataFrame(detail_rows)
    n_runs = int(df["run_id"].nunique()) if "run_id" in df.columns else len(csv_paths)
    summary_html = (
        "<div class='meta-box'>"
        f"Panel A denominator: <strong>{context_total:,}</strong> incident sepsis onsets. "
        f"Panel B denominator: <strong>{delay_total:,}</strong> sepsis episodes with an assigned delay bucket. "
        f"Runs: <strong>{n_runs}</strong>. "
        f"{threshold_note}"
        " Counts are raw simulated sepsis-episode counts from the accepted run(s); "
        "percentages are the primary quantities for interpretation. Counts are not population-scaled."
        "</div>\n"
        "<h2>Figure 10 Details</h2>\n"
        + _html_table(details)
        + "<p class='note'>Recovery/resolution, death, and censoring categories are assigned only "
        "for episodes that did not receive effective therapy before episode closure or end of observation.</p>\n"
    )

    footnotes = [
        "Figure 10 reads only aggregate columns in simulation_summary_run#.csv; no event-level "
        "sepsis CSV is required or produced.",
        "Sepsis onset is the modelled onset of underlying sepsis physiology, not the time of "
        "clinical recognition, healthcare presentation, or diagnosis. Some modelled sepsis "
        "episodes occur in people who never reach care or whose sepsis is not clinically "
        "recognised. \u2018No antibiotic active\u2019 means no active antibiotic in the model at "
        "underlying sepsis onset; it should not be interpreted as failure to administer "
        "antibiotics after recognised sepsis. Clinical recognition and attendance are modelled "
        "separately from underlying sepsis onset. Timing is measured in model days, so "
        "\u2018later in the same model day\u2019 is not a clinical hour-level delay.",
        "Panel A uses course-start context labels in priority order: targeted, empiric, "
        "other/prophylaxis, no active antibiotic, then unknown/legacy. Effective categories "
        "mean at least one currently active antibiotic met the activity threshold for that bacterium.",
        "Panel A classifies sepsis episodes by course-start treatment context at onset using priority "
        "order. Panel B classifies whether any active antibiotic was effective by the specified activity "
        "threshold, regardless of whether that antibiotic was empiric, targeted, prophylaxis, or "
        "other/background.",
        "Panel A counts are incremented on the timestep when a new sepsis episode begins. "
        "The denominator is incident sepsis onsets in baseline-policy rows with onset years 2022-2025.",
        "Panel B internally tracks each open sepsis episode until first effective therapy, resolution, "
        "death, or censoring, then assigns the aggregate delay bucket back to that episode's sepsis-onset "
        "timestep. The 2022-2025 filter is therefore also based on onset year.",
        "Activity uses the model's existing resistance-adjusted activity_r stored for each "
        "bacterium-drug pair. The Figure 10 threshold is output-only and does not change "
        "treatment selection or infection dynamics.",
        "Panel B separates antibiotics already effective at the pre-onset snapshot from therapy "
        "that first becomes effective later in the same model day.",
        "Because the model uses daily timesteps, 'later in the same model day' reflects event order "
        "within the daily rule update and should not be interpreted as a precise sub-day clinical timestamp.",
    ]
    if legacy_combined_present:
        footnotes.append(
            "This run predates the no-effective-therapy outcome split; no-effective episodes are "
            "shown as one combined recovery/death/censoring category for legacy rows."
        )
    elif no_effective_split_available:
        footnotes.append(
            "No-effective-therapy episodes are split by whether the episode closed by recovery, death, "
            "or censoring/end of observation."
        )
    if problems:
        footnotes.append("Some simulation summary CSVs were skipped: " + " ".join(problems[:3]))

    _save_figure(
        fig,
        out_dir,
        _F11_STEM,
        _F11_TITLE,
        "Baseline-policy sepsis-onset aggregate counts with onset years 2022-2025.",
        footnotes,
        agg=agg,
        extra_html=summary_html,
    )
    print(
        "  Figure 10: "
        f"{context_total:,} context count{'s' if context_total != 1 else ''} and "
        f"{delay_total:,} delay count{'s' if delay_total != 1 else ''} "
        f"from {len(csv_paths)} simulation CSV{'s' if len(csv_paths) != 1 else ''}."
    )


_F15_TITLE = "Figure 11. Resistance-adjusted antibiotic activity retained by bacterium, 2022–2025"
_F15_STEM = "Figure_11__activity_retained_by_bacterium"
_F15_REQUIRED_MESSAGE = (
    "Figure 11 requires simulation_summary CSV columns for activity_r_sum_by_bacteria "
    "and max_possible_activity_r_sum_by_bacteria."
)

_F15_KNOWN_BACTERIA_SLUGS = [
    "acinetobacter_baumannii",
    "citrobacter_spp.",
    "enterobacter_spp.",
    "enterococcus_faecalis",
    "enterococcus_faecium",
    "escherichia_coli",
    "klebsiella_pneumoniae",
    "morganella_spp.",
    "proteus_spp.",
    "serratia_spp.",
    "p_stuartii",
    "pseudomonas_aeruginosa",
    "stenotrophomonas_maltophilia",
    "staphylococcus_aureus",
    "staphylococcus_epidermidis",
    "streptococcus_pneumoniae",
    "salmonella_enterica_serovar_typhi",
    "salmonella_enterica_serovar_paratyphi_a",
    "invasive_non-typhoidal_salmonella_spp.",
    "shigella_spp.",
    "neisseria_gonorrhoeae",
    "streptococcus_pyogenes",
    "streptococcus_agalactiae",
    "haemophilus_influenzae",
    "chlamydia_trachomatis",
    "mycoplasma_genitalium",
    "vibrio_cholerae",
    "neisseria_meningitidis",
    "listeria_monocytogenes",
    "clostridioides_difficile",
    "bacteroides_fragilis",
    "campylobacter_jejuni",
    "enterobacter_cloacae",
    "yersinia_enterocolitica",
    "moraxella_catarrhalis",
    "treponema_pallidum",
    "bordetella_pertussis",
    "helicobacter_pylori",
    "mdr_mycobacterium_tuberculosis",
    "mycoplasma_pneumoniae",
    "legionella_pneumophila",
    "burkholderia_cepacia_complex",
]


def _figure_15_placeholder(out_dir: Path, agg: dict | None, message: str) -> None:
    fig, ax = plt.subplots(figsize=(9, 3.6))
    ax.text(
        0.5,
        0.5,
        f"{_F15_TITLE}\n\n{message}",
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=10.5,
        color="#555",
        bbox=dict(boxstyle="round,pad=0.6", fc="#f5f5f5", ec="#bbb"),
    )
    ax.set_axis_off()
    fig.subplots_adjust(left=0.03, right=0.97, top=0.92, bottom=0.08)
    _save_figure(fig, out_dir, _F15_STEM, _F15_TITLE, message, [], agg=agg)


def _figure_15_bacterium_label(slug: str) -> str:
    label = str(slug or "").strip().replace("_", " ")
    special = {
        "p stuartii": "P. stuartii",
        "mdr mycobacterium tuberculosis": "MDR Mycobacterium tuberculosis",
        "invasive non-typhoidal salmonella spp.": "Invasive non-typhoidal Salmonella spp.",
        "salmonella enterica serovar typhi": "Salmonella enterica serovar Typhi",
        "salmonella enterica serovar paratyphi a": "Salmonella enterica serovar Paratyphi A",
    }
    lowered = label.lower()
    if lowered in special:
        return special[lowered]
    words = lowered.split()
    if not words:
        return str(slug)
    formatted = [words[0].capitalize()]
    for word in words[1:]:
        formatted.append("spp." if word == "spp." else word)
    return " ".join(formatted)


def _figure_15_parse_vector_cell(value: object) -> list[float]:
    if value is None or (isinstance(value, float) and np.isnan(value)):
        return []
    text = str(value).strip()
    if not text or text.lower() in {"nan", "none", "null"}:
        return []
    text = text.strip("[]()")
    if not text:
        return []
    parts = [part.strip() for part in re.split(r"[;,]", text) if part.strip()]
    if len(parts) == 1 and " " in parts[0]:
        parts = [part.strip() for part in re.split(r"\s+", parts[0]) if part.strip()]
    values: list[float] = []
    for part in parts:
        try:
            values.append(float(part))
        except ValueError:
            values.append(float("nan"))
    return values


def _figure_15_extend_array(values: np.ndarray, target_len: int) -> np.ndarray:
    if len(values) >= target_len:
        return values
    extended = np.zeros(target_len, dtype=float)
    if len(values):
        extended[: len(values)] = values
    return extended


def _figure_15_detect_activity_columns(
    columns: list[str],
) -> tuple[dict[str, tuple[str, str]], tuple[str, str] | None, list[str], bool, str]:
    column_set = set(columns)
    wide_pairs: dict[str, tuple[str, str]] = {}
    missing_denominator: list[str] = []

    for col in columns:
        if (
            col.endswith("_activity_r_sum")
            and not col.endswith("_activity_r_pure_sum")
            and not col.endswith("_max_possible_activity_r_sum")
            and not col.startswith("max_possible_")
        ):
            slug = col[: -len("_activity_r_sum")]
            denom_col = f"{slug}_max_possible_activity_r_sum"
            if denom_col in column_set:
                wide_pairs[slug] = (col, denom_col)
            else:
                missing_denominator.append(slug)

    prefix_patterns = [
        ("activity_r_sum_by_bacteria__", "max_possible_activity_r_sum_by_bacteria__"),
        ("activity_r_sum_by_bacteria_", "max_possible_activity_r_sum_by_bacteria_"),
    ]
    for num_prefix, denom_prefix in prefix_patterns:
        for col in columns:
            if num_prefix.endswith("_") and not num_prefix.endswith("__") and col.startswith(
                "activity_r_sum_by_bacteria__"
            ):
                continue
            if not col.startswith(num_prefix):
                continue
            slug = col[len(num_prefix):]
            denom_col = f"{denom_prefix}{slug}"
            if denom_col in column_set:
                wide_pairs[slug] = (col, denom_col)
            else:
                missing_denominator.append(slug)

    vector_pair = None
    if (
        "activity_r_sum_by_bacteria" in column_set
        and "max_possible_activity_r_sum_by_bacteria" in column_set
    ):
        vector_pair = ("activity_r_sum_by_bacteria", "max_possible_activity_r_sum_by_bacteria")

    pure_present = any(
        col.endswith("_activity_r_pure_sum")
        or col.endswith("_max_possible_activity_r_pure_sum")
        or col in {
            "activity_r_pure_sum_by_bacteria",
            "max_possible_activity_r_pure_sum_by_bacteria",
        }
        for col in columns
    )
    detection = "wide per-bacterium suffix columns" if wide_pairs else "vector by-bacteria columns"
    return wide_pairs, vector_pair, missing_denominator, pure_present, detection


_FIGURE_15_ROWS_CACHE: dict[
    Path,
    tuple[list[dict[str, object]], str | None, set[str], bool, str | None],
] = {}


def _copy_figure_15_rows_result(
    result: tuple[list[dict[str, object]], str | None, set[str], bool, str | None],
) -> tuple[list[dict[str, object]], str | None, set[str], bool, str | None]:
    rows, problem, no_denominator, saw_pure, detection = result
    return [dict(row) for row in rows], problem, set(no_denominator), saw_pure, detection


def _figure_15_rows_from_simulation_csv_uncached(
    csv_path: Path,
) -> tuple[list[dict[str, object]], str | None, set[str], bool, str | None]:
    columns = _simulation_csv_column_names(csv_path)
    if columns is None:
        return [], f"{csv_path.name}: could not read simulation CSV header.", set(), False, None

    wide_pairs, vector_pair, missing_denominator, pure_present, detection = (
        _figure_15_detect_activity_columns(columns)
    )
    if not wide_pairs and vector_pair is None:
        return [], f"{csv_path.name}: missing activity_r_sum/max_possible activity columns.", set(), pure_present, None

    optional = ["policy_option", "run_id", "simulation_year", "year", "time_in_years", "time_step"]
    wanted = set(optional)
    if wide_pairs:
        for num_col, denom_col in wide_pairs.values():
            wanted.add(num_col)
            wanted.add(denom_col)
    elif vector_pair is not None:
        wanted.update(vector_pair)

    try:
        df = _read_csv_selected(csv_path, wanted)
    except (FileNotFoundError, ValueError, OSError) as exc:
        return [], f"{csv_path.name}: could not load Figure 11 columns ({exc}).", set(), pure_present, detection

    if "policy_option" in df.columns:
        policy = pd.to_numeric(df["policy_option"], errors="coerce")
        df = df[policy == 0].copy()

    df["simulation_year_for_f15"] = _simulation_year_series(df)
    df = df[(df["simulation_year_for_f15"] >= 2022.0) & (df["simulation_year_for_f15"] < 2026.0)].copy()
    if df.empty:
        return [], f"{csv_path.name}: no baseline-policy rows in 2022-2025.", set(missing_denominator), pure_present, detection

    grouped = df.groupby("run_id", dropna=False) if "run_id" in df.columns else [(csv_path.stem, df)]
    rows: list[dict[str, object]] = []
    no_denominator = set(missing_denominator)

    if wide_pairs:
        for num_col, denom_col in wide_pairs.values():
            df[num_col] = pd.to_numeric(df[num_col], errors="coerce")
            df[denom_col] = pd.to_numeric(df[denom_col], errors="coerce")
        for run_key, run_df in grouped:
            for slug, (num_col, denom_col) in wide_pairs.items():
                numerator = float(run_df[num_col].sum(skipna=True))
                denominator = float(run_df[denom_col].sum(skipna=True))
                if not np.isfinite(denominator) or denominator <= 0.0:
                    no_denominator.add(slug)
                    continue
                rows.append({
                    "source": csv_path.name,
                    "run": str(run_key),
                    "bacterium_slug": slug,
                    "bacterium": _figure_15_bacterium_label(slug),
                    "mean_activity_r_percent": 100.0 * numerator / denominator,
                    "numerator_sum": numerator,
                    "denominator_sum": denominator,
                    "detection": detection,
                })
    elif vector_pair is not None:
        num_col, denom_col = vector_pair
        for run_key, run_df in grouped:
            numerator_sums = np.zeros(0, dtype=float)
            denominator_sums = np.zeros(0, dtype=float)
            for num_value, denom_value in zip(run_df[num_col], run_df[denom_col]):
                numerator_values = np.array(_figure_15_parse_vector_cell(num_value), dtype=float)
                denominator_values = np.array(_figure_15_parse_vector_cell(denom_value), dtype=float)
                target_len = max(len(numerator_sums), len(numerator_values), len(denominator_values))
                numerator_sums = _figure_15_extend_array(numerator_sums, target_len)
                denominator_sums = _figure_15_extend_array(denominator_sums, target_len)
                numerator_values = _figure_15_extend_array(numerator_values, target_len)
                denominator_values = _figure_15_extend_array(denominator_values, target_len)
                numerator_sums += np.nan_to_num(numerator_values, nan=0.0)
                denominator_sums += np.nan_to_num(denominator_values, nan=0.0)

            for idx, denominator in enumerate(denominator_sums):
                slug = (
                    _F15_KNOWN_BACTERIA_SLUGS[idx]
                    if idx < len(_F15_KNOWN_BACTERIA_SLUGS)
                    else f"bacterium_{idx + 1}"
                )
                if not np.isfinite(denominator) or denominator <= 0.0:
                    no_denominator.add(slug)
                    continue
                numerator = float(numerator_sums[idx])
                rows.append({
                    "source": csv_path.name,
                    "run": str(run_key),
                    "bacterium_slug": slug,
                    "bacterium": _figure_15_bacterium_label(slug),
                    "mean_activity_r_percent": 100.0 * numerator / float(denominator),
                    "numerator_sum": numerator,
                    "denominator_sum": float(denominator),
                    "detection": detection,
                })

    return rows, None, no_denominator, pure_present, detection


def _figure_15_rows_from_simulation_csv(
    csv_path: Path,
) -> tuple[list[dict[str, object]], str | None, set[str], bool, str | None]:
    resolved = csv_path.resolve()
    if resolved not in _FIGURE_15_ROWS_CACHE:
        _FIGURE_15_ROWS_CACHE[resolved] = _figure_15_rows_from_simulation_csv_uncached(csv_path)
    return _copy_figure_15_rows_result(_FIGURE_15_ROWS_CACHE[resolved])


def make_figure_15_mean_activity_by_bacteria(
    csv_paths: list[Path],
    out_dir: Path,
    agg: dict | None = None,
) -> None:
    if not csv_paths:
        _figure_15_placeholder(out_dir, agg, _F15_REQUIRED_MESSAGE)
        return

    rows: list[dict[str, object]] = []
    problems: list[str] = []
    denominatorless: set[str] = set()
    pure_present = False
    detections: set[str] = set()

    for csv_path in csv_paths:
        run_rows, problem, no_denominator, saw_pure, detection = _figure_15_rows_from_simulation_csv(csv_path)
        rows.extend(run_rows)
        denominatorless.update(no_denominator)
        pure_present = pure_present or saw_pure
        if detection:
            detections.add(detection)
        if problem:
            problems.append(problem)

    if not rows:
        detail = _F15_REQUIRED_MESSAGE
        if problems and not any("missing activity" in p for p in problems):
            detail = " ".join(problems)
        _figure_15_placeholder(out_dir, agg, detail)
        return

    df = pd.DataFrame(rows)
    included_slugs = set(df["bacterium_slug"].dropna().astype(str))
    excluded_denominatorless = sorted(denominatorless - included_slugs)
    summary = (
        df.groupby(["bacterium_slug", "bacterium"], as_index=False)
        .agg(
            activity_median=("mean_activity_r_percent", "median"),
            activity_p5=("mean_activity_r_percent", lambda s: float(np.nanpercentile(s, 5))),
            activity_p95=("mean_activity_r_percent", lambda s: float(np.nanpercentile(s, 95))),
            numerator_median=("numerator_sum", "median"),
            denominator_median=("denominator_sum", "median"),
            n=("mean_activity_r_percent", "count"),
        )
        .sort_values("activity_median", ascending=True)
        .reset_index(drop=True)
    )
    n_runs = int(df[["source", "run"]].drop_duplicates().shape[0])

    fig, ax = plt.subplots(figsize=(9, max(5.0, 0.34 * len(summary))))
    y = np.arange(len(summary))
    medians = summary["activity_median"].to_numpy(float)
    err_lo = np.clip(medians - summary["activity_p5"].to_numpy(float), 0, None)
    err_hi = np.clip(summary["activity_p95"].to_numpy(float) - medians, 0, None)
    ax.barh(
        y,
        medians,
        0.58,
        color="#3F7F93",
        alpha=0.9,
        xerr=[err_lo, err_hi] if n_runs > 1 else None,
        error_kw={"elinewidth": 0.9, "ecolor": "#263238", "capthick": 0.9, "capsize": 2.5},
    )
    ax.set_yticks(y)
    ax.set_yticklabels(summary["bacterium"].values, fontsize=7.5, fontstyle="italic")
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xlabel("Resistance-adjusted activity retained (% of no-resistance activity)", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", linewidth=0.4, alpha=0.5)
    fig.suptitle(_F15_TITLE, fontsize=10.5, fontweight="bold")
    fig.tight_layout()

    denominator_total = float(summary["denominator_median"].sum())

    def _fmt_activity(value: object) -> str:
        return f"{float(value):.1f}" if pd.notna(value) and np.isfinite(float(value)) else "—"

    def _fmt_share(value: object) -> str:
        return f"{float(value):.2f}" if pd.notna(value) and np.isfinite(float(value)) else "—"

    def _fmt_compact(value: object) -> str:
        if value is None or pd.isna(value):
            return "—"
        value_f = float(value)
        if not np.isfinite(value_f):
            return "—"
        return f"{value_f:.3g}"

    denominator_table = summary.copy()
    denominator_table["denominator_share_percent"] = (
        100.0 * denominator_table["denominator_median"] / denominator_total
        if denominator_total > 0.0
        else np.nan
    )
    denominator_table["denominator_flag"] = np.where(
        denominator_table["denominator_share_percent"] < 0.1,
        "low denominator",
        "",
    )
    denominator_table = pd.DataFrame({
        "Bacterium": denominator_table["bacterium"],
        "Activity retained (%)": denominator_table["activity_median"].map(_fmt_activity),
        "Activity numerator": denominator_table["numerator_median"].map(_fmt_compact),
        "No-resistance activity denominator": denominator_table["denominator_median"].map(_fmt_compact),
        "Denominator share (%)": denominator_table["denominator_share_percent"].map(_fmt_share),
        "Denominator flag": denominator_table["denominator_flag"],
    })
    lowest_names = summary.head(5)["bacterium"].tolist()
    lowest_context = "this run" if n_runs == 1 else "these runs"
    lowest_sentence = (
        f"<p class='note'>The lowest retained-activity bacteria in {lowest_context} were: "
        + ", ".join(lowest_names)
        + ".</p>\n"
        if lowest_names
        else ""
    )
    extra_html = (
        "<h2>Figure 11 Denominator Table</h2>\n"
        + _html_table(denominator_table)
        + lowest_sentence
    )

    footnotes = [
        "Activity retained is calculated as sum(activity_r) divided by sum(max_possible_activity_r) "
        "across baseline-policy rows in the 2022-2025 calibration window. The denominator "
        "represents expected activity under the same drug exposure if resistance were absent. "
        "The metric is weighted by actual antibiotic exposure and therefore can differ from "
        "simple resistance prevalence. Low values indicate that the drugs being used for that "
        "bacterium have little retained activity after resistance.",
        "This is not an absolute measure of treatment adequacy, and it is not the percentage "
        "of infections resistant. It measures resistance-related loss of activity conditional "
        "on the antibiotics used.",
        "Because this metric is weighted by drug exposure, organisms can have lower retained "
        "activity than their simple resistance prevalence would suggest if treated person-days "
        "are concentrated among resistant or poorly responding infections.",
        "Bacteria with no treated-infection activity denominator are excluded from the plot "
        "and denominator table.",
        f"Source: matched simulation_summary_*.csv files. Values are medians across {n_runs} "
        f"simulation run{'s' if n_runs > 1 else ''}; "
        f"{'error bars show 5th-95th percentile ranges. ' if n_runs > 1 else 'no error bars are shown for a single run. '}",
    ]
    if detections:
        footnotes.append("Activity columns detected as: " + "; ".join(sorted(detections)) + ".")
    if pure_present:
        footnotes.append(
            "Pure activity columns were present, but this figure uses the drug-exposure-weighted "
            "activity metric rather than the pure potency-only metric."
        )
    if excluded_denominatorless:
        labels = [_figure_15_bacterium_label(slug) for slug in excluded_denominatorless]
        footnotes.append(
            "Excluded due to no treated-infection activity denominator: "
            + ", ".join(labels)
            + "."
        )

    _save_figure(
        fig,
        out_dir,
        _F15_STEM,
        _F15_TITLE,
        "Higher values mean better retained antibiotic activity; lower values mean resistance "
        "is reducing the activity of the drugs being used for that bacterium.",
        footnotes,
        agg=agg,
        extra_html=extra_html,
    )


_F19_TITLE = "Figure 12. Antibiotic exposure distribution by bacterium, 2022\u20132025"
_F19_STEM = "Figure_12__distribution_drug_use_by_bacteria"
_F19_REQUIRED_MESSAGE = (
    "Figure 12 requires simulation_summary CSV output for "
    "currently_on_drug_by_bacteria_drug. Re-run the Rust simulation with full "
    "per-bacterium/per-drug summary content enabled, or provide simulation_summary "
    "files containing this field."
)

_F19_OTHER_UNMAPPED_CLASS = "Other / unmapped"
_F19_OTHER_CLASSES = "Other classes"
_F19_MAX_DISPLAY_CLASSES = 16

_F19_KNOWN_DRUG_SLUGS: list[str] = [
    "sulfanilamide",
    "penicillin_g",
    "ampicillin",
    "amoxicillin",
    "piperacillin",
    "ticarcillin",
    "cephalexin",
    "cefazolin",
    "cefuroxime",
    "ceftriaxone",
    "ceftazidime",
    "cefepime",
    "ceftaroline",
    "ceftolozane_tazobactam",
    "cefiderocol",
    "meropenem",
    "imipenem_c",
    "ertapenem",
    "aztreonam",
    "erythromycin",
    "azithromycin",
    "clarithromycin",
    "clindamycin",
    "gentamicin",
    "tobramycin",
    "amikacin",
    "ciprofloxacin",
    "levofloxacin",
    "moxifloxacin",
    "ofloxacin",
    "tetracycline",
    "doxycycline",
    "minocycline",
    "tigecycline",
    "vancomycin",
    "teicoplanin",
    "dalbavancin",
    "linezolid",
    "tedizolid",
    "daptomycin",
    "quinu_dalfo",
    "trim_sulf",
    "chloramphenicol",
    "nitrofurantoin",
    "fosfomycin",
    "retapamulin",
    "fusidic_a",
    "metronidazole",
    "fidaxomicin",
    "furazolidone",
    "rifampicin",
    "amoxicillin_clavulanate",
    "piperacillin_tazobactam",
    "ampicillin_sulbactam",
    "ticarcillin_clavulanate",
    "ceftazidime_avibactam",
    "meropenem_vaborbactam",
    "colistin",
    "flucloxacillin",
    "aztreonam_avibactam",
    "cefixime",
    "nalidixic_acid",
]


def _figure_19_fallback_drug_class_map() -> dict[str, str]:
    groups = {
        "Sulfonamides (J01E)": ["sulfanilamide", "trim_sulf"],
        "Penicillins (J01C)": [
            "penicillin_g", "ampicillin", "amoxicillin", "piperacillin",
            "ticarcillin", "flucloxacillin",
        ],
        "Beta-lactamase combinations (J01CR)": [
            "amoxicillin_clavulanate", "piperacillin_tazobactam",
            "ampicillin_sulbactam", "ticarcillin_clavulanate",
        ],
        "Cephalosporins 1-2G": ["cephalexin", "cefazolin", "cefuroxime"],
        "Cephalosporins 3G": ["ceftriaxone", "ceftazidime", "cefixime"],
        "Cephalosporins 3G/BLI": ["ceftolozane_tazobactam"],
        "Cephalosporins 4G": ["cefepime"],
        "Anti-MRSA Cephalosporins (5G)": ["ceftaroline"],
        "Siderophore Cephalosporins": ["cefiderocol"],
        "Carbapenems (J01DH)": ["meropenem", "imipenem_c", "ertapenem"],
        "Novel BL/BLI": [
            "ceftazidime_avibactam", "meropenem_vaborbactam",
            "aztreonam_avibactam",
        ],
        "Monobactams": ["aztreonam"],
        "Macrolides (J01F)": ["erythromycin", "azithromycin", "clarithromycin"],
        "Lincosamides (J01FF)": ["clindamycin"],
        "Aminoglycosides (J01G)": ["gentamicin", "tobramycin", "amikacin"],
        "Fluoroquinolones (J01M)": [
            "ciprofloxacin", "levofloxacin", "moxifloxacin", "ofloxacin",
            "nalidixic_acid",
        ],
        "Tetracyclines (J01A)": ["tetracycline", "doxycycline", "minocycline", "tigecycline"],
        "Glycopeptides (J01XA)": ["vancomycin", "teicoplanin"],
        "Lipoglycopeptides": ["dalbavancin"],
        "Oxazolidinones (J01XX)": ["linezolid", "tedizolid"],
        "Lipopeptides (J01XX09)": ["daptomycin"],
        "Streptogramins (J01FG)": ["quinu_dalfo"],
        "Chloramphenicol (J01BA)": ["chloramphenicol"],
        "Nitrofurans (J01XE)": ["nitrofurantoin", "furazolidone"],
        "Fosfomycin (J01XX01)": ["fosfomycin"],
        "Pleuromutilins": ["retapamulin"],
        "Fusidic acid (J01XC)": ["fusidic_a"],
        "Nitroimidazoles": ["metronidazole"],
        "Fidaxomicin": ["fidaxomicin"],
        "Rifamycins (J04AB)": ["rifampicin"],
        "Polymyxins (J01XB)": ["colistin"],
    }
    mapping: dict[str, str] = {}
    for drug_class, drugs in groups.items():
        for drug in drugs:
            mapping[drug] = drug_class
    return mapping


def _figure_19_normalize_drug_name(value: object) -> str:
    text = str(value or "").strip().lower()
    text = re.sub(r"\s*/\s*", "_", text)
    text = re.sub(r"[^a-z0-9]+", "_", text)
    text = re.sub(r"_+", "_", text).strip("_")
    aliases = {
        "co_trimoxazole": "trim_sulf",
        "trimethoprim_sulfamethoxazole": "trim_sulf",
        "tmp_smx": "trim_sulf",
        "quinupristin_dalfopristin": "quinu_dalfo",
        "imipenem_cilastatin": "imipenem_c",
        "fusidic_acid": "fusidic_a",
    }
    return aliases.get(text, text)


def _figure_19_drug_class_map(agg: dict | None) -> tuple[dict[str, str], str]:
    mapping = _figure_19_fallback_drug_class_map()
    source = "fallback drug-class map"
    rb = agg.get("resistance_benchmarks", pd.DataFrame()) if agg is not None else pd.DataFrame()
    if rb is not None and not rb.empty and {"Drug", "Class"}.issubset(rb.columns):
        rb_mapping: dict[str, str] = {}
        for _, row in rb[["Drug", "Class"]].dropna().iterrows():
            drug = _figure_19_normalize_drug_name(row["Drug"])
            drug_class = str(row["Class"]).strip()
            if drug and drug_class and drug_class.lower() not in {"nan", "none", "---"}:
                rb_mapping[drug] = drug_class
        if rb_mapping:
            mapping.update(rb_mapping)
            source = "resistance_benchmarks Drug/Class columns plus fallback drug-class map"
    return mapping, source


def _figure_19_placeholder(out_dir: Path, agg: dict | None, message: str) -> None:
    fig, ax = plt.subplots(figsize=(10, 3.8))
    ax.text(
        0.5,
        0.5,
        f"{_F19_TITLE}\n\n{message}",
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=10.5,
        color="#555",
        bbox=dict(boxstyle="round,pad=0.6", fc="#f5f5f5", ec="#bbb"),
    )
    ax.set_axis_off()
    fig.subplots_adjust(left=0.03, right=0.97, top=0.92, bottom=0.08)
    _save_figure(
        fig,
        out_dir,
        _F19_STEM,
        _F19_TITLE,
        message,
        [
            "Cells show the percentage distribution of active antibiotic exposure among people "
            "infected with each bacterium during the 2022-2025 calibration window. The numerator "
            "is summed currently_on_drug_by_bacteria_drug for each bacterium-drug pair, aggregated "
            "to drug class. This is not necessarily the drug prescribed specifically for that "
            "bacterium; it includes empiric therapy, targeted therapy, combination therapy, and "
            "bystander exposure from antibiotics used for other concurrent infections.",
            "Rows are normalised within bacterium and sum to 100% across displayed classes plus "
            "'Other classes'.",
        ],
        agg=agg,
    )


def _figure_19_detect_exposure_columns(
    columns: list[str],
    known_drugs: list[str],
) -> tuple[dict[tuple[str, str], str], str | None, str | None]:
    column_set = set(columns)
    wide: dict[tuple[str, str], str] = {}
    double_prefix = "currently_on_drug_by_bacteria_drug__"
    single_prefix = "currently_on_drug_by_bacteria_drug_"

    for col in columns:
        if not col.startswith(double_prefix):
            continue
        remainder = col[len(double_prefix):]
        if "__" not in remainder:
            continue
        bacterium, drug = remainder.split("__", 1)
        if bacterium and drug:
            wide[(bacterium, _figure_19_normalize_drug_name(drug))] = col

    sorted_drugs = sorted(set(known_drugs), key=len, reverse=True)
    for col in columns:
        if not col.startswith(single_prefix) or col.startswith(double_prefix):
            continue
        remainder = col[len(single_prefix):]
        for drug in sorted_drugs:
            suffix = f"_{drug}"
            if remainder.endswith(suffix) and len(remainder) > len(suffix):
                bacterium = remainder[: -len(suffix)]
                wide[(bacterium, drug)] = col
                break

    legacy_matches = 0
    for bacterium in _F15_KNOWN_BACTERIA_SLUGS:
        for drug in sorted_drugs:
            col = f"{bacterium}_currently_on_drug_{drug}"
            if col in column_set:
                wide[(bacterium, drug)] = col
                legacy_matches += 1

    vector_col = (
        "currently_on_drug_by_bacteria_drug"
        if "currently_on_drug_by_bacteria_drug" in column_set
        else None
    )
    if wide:
        if legacy_matches:
            return wide, None, "legacy Rust <bacterium>_currently_on_drug_<drug> columns"
        return wide, None, "wide per-bacterium/per-drug columns"
    if vector_col:
        return {}, vector_col, "flat vector column"
    return {}, None, None


def _figure_19_rows_from_simulation_csv(
    csv_path: Path,
    drug_class_map: dict[str, str],
) -> tuple[list[dict[str, object]], str | None, set[str], str | None]:
    columns = _simulation_csv_column_names(csv_path)
    if columns is None:
        return [], f"{csv_path.name}: could not read simulation CSV header.", set(), None

    known_drugs = list(dict.fromkeys(_F19_KNOWN_DRUG_SLUGS + sorted(drug_class_map)))
    wide_cols, vector_col, detection = _figure_19_detect_exposure_columns(columns, known_drugs)
    if not wide_cols and vector_col is None:
        return [], f"{csv_path.name}: missing currently_on_drug_by_bacteria_drug columns.", set(), None

    optional = ["policy_option", "run_id", "simulation_year", "year", "time_in_years", "time_step"]
    wanted = set(optional)
    if wide_cols:
        wanted.update(wide_cols.values())
    elif vector_col is not None:
        wanted.add(vector_col)

    try:
        df = _read_csv_selected(csv_path, wanted)
    except (FileNotFoundError, ValueError, OSError) as exc:
        return [], f"{csv_path.name}: could not load Figure 12 columns ({exc}).", set(), detection

    if "policy_option" in df.columns:
        policy = pd.to_numeric(df["policy_option"], errors="coerce")
        df = df[policy == 0].copy()

    df["simulation_year_for_f19"] = _simulation_year_series(df)
    df = df[(df["simulation_year_for_f19"] >= 2022.0) & (df["simulation_year_for_f19"] < 2026.0)].copy()
    if df.empty:
        return [], f"{csv_path.name}: no baseline-policy rows in 2022-2025.", set(), detection

    grouped = df.groupby("run_id", dropna=False) if "run_id" in df.columns else [(csv_path.stem, df)]
    rows: list[dict[str, object]] = []
    unmapped: set[str] = set()

    def _record_count(run_key: object, bacterium_slug: str, drug_slug: str, count: float) -> None:
        if not np.isfinite(count) or count <= 0.0:
            return
        drug_class = drug_class_map.get(drug_slug)
        if drug_class is None:
            drug_class = _F19_OTHER_UNMAPPED_CLASS
            unmapped.add(drug_slug)
        rows.append({
            "source": csv_path.name,
            "run": str(run_key),
            "bacterium_slug": bacterium_slug,
            "bacterium": _figure_15_bacterium_label(bacterium_slug),
            "drug": drug_slug,
            "drug_class": drug_class,
            "exposure_count": float(count),
            "detection": detection,
        })

    if wide_cols:
        for col in wide_cols.values():
            df[col] = pd.to_numeric(df[col], errors="coerce")
        for run_key, run_df in grouped:
            for (bacterium_slug, drug_slug), col in wide_cols.items():
                _record_count(run_key, bacterium_slug, drug_slug, float(run_df[col].sum(skipna=True)))
    elif vector_col is not None:
        parsed_lengths = [
            len(values)
            for values in (_figure_15_parse_vector_cell(value) for value in df[vector_col])
            if values
        ]
        if not parsed_lengths:
            return [], f"{csv_path.name}: Figure 12 vector column contained no numeric values.", unmapped, detection
        vector_len = max(parsed_lengths)
        n_bacteria = len(_F15_KNOWN_BACTERIA_SLUGS)
        if vector_len % n_bacteria != 0:
            return (
                [],
                f"{csv_path.name}: Figure 12 vector length {vector_len} is not divisible "
                f"by known bacteria count {n_bacteria}.",
                unmapped,
                detection,
            )
        n_drugs = vector_len // n_bacteria
        drug_slugs = list(_F19_KNOWN_DRUG_SLUGS[:n_drugs])
        if n_drugs > len(drug_slugs):
            drug_slugs.extend(f"drug_{idx + 1}" for idx in range(len(drug_slugs), n_drugs))

        for run_key, run_df in grouped:
            sums = np.zeros(vector_len, dtype=float)
            for value in run_df[vector_col]:
                values = np.array(_figure_15_parse_vector_cell(value), dtype=float)
                values = _figure_15_extend_array(values, vector_len)
                sums += np.nan_to_num(values, nan=0.0)
            for bacterium_idx, bacterium_slug in enumerate(_F15_KNOWN_BACTERIA_SLUGS):
                offset = bacterium_idx * n_drugs
                for drug_idx, drug_slug in enumerate(drug_slugs):
                    _record_count(run_key, bacterium_slug, drug_slug, float(sums[offset + drug_idx]))

    return rows, None, unmapped, detection


def _figure_19_activity_sort_values(csv_paths: list[Path]) -> dict[str, float]:
    rows: list[dict[str, object]] = []
    for csv_path in csv_paths:
        run_rows, _, _, _, _ = _figure_15_rows_from_simulation_csv(csv_path)
        rows.extend(run_rows)
    if not rows:
        return {}
    df = pd.DataFrame(rows)
    return (
        df.groupby("bacterium_slug")["mean_activity_r_percent"]
        .median()
        .dropna()
        .to_dict()
    )


def _figure_19_compact_count(value: object) -> str:
    if value is None or pd.isna(value):
        return "-"
    value_f = float(value)
    if not np.isfinite(value_f):
        return "-"
    abs_value = abs(value_f)
    if abs_value >= 1_000_000_000:
        return f"{value_f / 1_000_000_000:.2f}B"
    if abs_value >= 1_000_000:
        return f"{value_f / 1_000_000:.2f}M"
    if abs_value >= 1_000:
        return f"{value_f / 1_000:.1f}K"
    return f"{value_f:.0f}"


def make_figure_19_antibiotic_exposure_distribution(
    csv_paths: list[Path],
    out_dir: Path,
    agg: dict | None = None,
) -> None:
    if not csv_paths:
        _figure_19_placeholder(out_dir, agg, _F19_REQUIRED_MESSAGE)
        return

    drug_class_map, class_mapping_source = _figure_19_drug_class_map(agg)
    rows: list[dict[str, object]] = []
    problems: list[str] = []
    unmapped: set[str] = set()
    detections: set[str] = set()

    for csv_path in csv_paths:
        run_rows, problem, run_unmapped, detection = _figure_19_rows_from_simulation_csv(
            csv_path,
            drug_class_map,
        )
        rows.extend(run_rows)
        unmapped.update(run_unmapped)
        if detection:
            detections.add(detection)
        if problem:
            problems.append(problem)

    if not rows:
        detail = _F19_REQUIRED_MESSAGE
        if problems and not any("missing currently_on_drug_by_bacteria_drug" in p for p in problems):
            detail = " ".join(problems)
        _figure_19_placeholder(out_dir, agg, detail)
        return

    count_df = pd.DataFrame(rows)
    class_counts = (
        count_df.groupby(["source", "run", "bacterium_slug", "bacterium", "drug_class"], as_index=False)
        ["exposure_count"]
        .sum()
    )
    class_totals = class_counts.groupby("drug_class")["exposure_count"].sum().sort_values(ascending=False)
    top_classes = [
        cls
        for cls in class_totals.index.tolist()
        if cls != _F19_OTHER_CLASSES
    ][: _F19_MAX_DISPLAY_CLASSES]
    if len(class_totals) > len(top_classes):
        display_classes = top_classes + [_F19_OTHER_CLASSES]
    else:
        display_classes = top_classes

    class_counts["display_class"] = np.where(
        class_counts["drug_class"].isin(top_classes),
        class_counts["drug_class"],
        _F19_OTHER_CLASSES,
    )
    display_counts = (
        class_counts.groupby(
            ["source", "run", "bacterium_slug", "bacterium", "display_class"],
            as_index=False,
        )["exposure_count"]
        .sum()
    )
    display_counts["denominator"] = display_counts.groupby(
        ["source", "run", "bacterium_slug"]
    )["exposure_count"].transform("sum")
    display_counts = display_counts[display_counts["denominator"] > 0.0].copy()
    display_counts["share_percent"] = (
        100.0 * display_counts["exposure_count"] / display_counts["denominator"]
    )

    summary = (
        display_counts.groupby(["bacterium_slug", "bacterium", "display_class"], as_index=False)
        .agg(share_percent=("share_percent", "median"))
    )
    row_totals = summary.groupby("bacterium_slug")["share_percent"].transform("sum")
    summary["share_percent"] = np.where(
        row_totals > 0.0,
        100.0 * summary["share_percent"] / row_totals,
        summary["share_percent"],
    )

    denominators = (
        display_counts[["source", "run", "bacterium_slug", "bacterium", "denominator"]]
        .drop_duplicates()
        .groupby(["bacterium_slug", "bacterium"], as_index=False)
        .agg(denominator_median=("denominator", "median"))
    )
    activity_sort = _figure_19_activity_sort_values(csv_paths)
    if activity_sort:
        denominators["_sort_activity"] = denominators["bacterium_slug"].map(activity_sort)
        denominators = denominators.sort_values(
            ["_sort_activity", "denominator_median"],
            ascending=[True, False],
            na_position="last",
        )
    else:
        denominators = denominators.sort_values("denominator_median", ascending=False)
    bacteria_order = denominators["bacterium_slug"].tolist()

    pivot = (
        summary.pivot(index="bacterium_slug", columns="display_class", values="share_percent")
        .reindex(index=bacteria_order, columns=display_classes)
        .fillna(0.0)
    )
    labels = (
        denominators.set_index("bacterium_slug")
        .reindex(bacteria_order)["bacterium"]
        .fillna(pd.Series(bacteria_order, index=bacteria_order))
        .tolist()
    )
    class_labels = [_F2_CLASS_SHORT.get(cls, cls) for cls in display_classes]

    fig_width = max(10.5, 3.5 + 0.62 * len(display_classes))
    fig_height = max(6.0, 2.2 + 0.34 * len(bacteria_order))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height))
    matrix = pivot.to_numpy(dtype=float)
    observed_max = float(np.nanmax(matrix)) if matrix.size else 0.0
    vmax = min(100.0, max(10.0, math.ceil(observed_max / 5.0) * 5.0)) if observed_max > 0 else 100.0
    im = ax.imshow(matrix, aspect="auto", cmap="YlGnBu", vmin=0.0, vmax=vmax)
    cbar = fig.colorbar(im, ax=ax, fraction=0.028, pad=0.015)
    cbar.set_label("Share of antibiotic exposure while infected (%)", fontsize=9)
    ax.set_xticks(np.arange(len(display_classes)))
    ax.set_xticklabels(class_labels, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(np.arange(len(bacteria_order)))
    ax.set_yticklabels(labels, fontsize=7.5, fontstyle="italic")
    ax.set_xlabel("Drug class", fontsize=10)
    ax.set_ylabel("Bacterium", fontsize=10)
    ax.set_title(_F19_TITLE, fontsize=10.5, fontweight="bold", pad=10)
    ax.set_xticks(np.arange(-0.5, len(display_classes), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(bacteria_order), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=0.45)
    ax.tick_params(which="minor", bottom=False, left=False)
    fig.tight_layout()

    top_by_bacterium = (
        summary.sort_values("share_percent", ascending=False)
        .groupby(["bacterium_slug", "bacterium"], as_index=False)
        .first()[["bacterium_slug", "display_class", "share_percent"]]
        .rename(columns={"display_class": "top_class", "share_percent": "top_share"})
    )
    nonzero_counts = (
        summary[summary["share_percent"] > 0.0]
        .groupby("bacterium_slug")["display_class"]
        .nunique()
        .rename("nonzero_classes")
        .reset_index()
    )
    table_df = (
        denominators.merge(top_by_bacterium, on="bacterium_slug", how="left")
        .merge(nonzero_counts, on="bacterium_slug", how="left")
    )
    denom_table = pd.DataFrame({
        "Bacterium": table_df["bacterium"],
        "Total antibiotic-exposure denominator": table_df["denominator_median"].map(_figure_19_compact_count),
        "Top drug class": table_df["top_class"].fillna("-"),
        "Top drug-class share (%)": table_df["top_share"].map(
            lambda v: f"{float(v):.1f}" if pd.notna(v) and np.isfinite(float(v)) else "-"
        ),
        "Number of displayed classes with nonzero exposure": table_df["nonzero_classes"].fillna(0).astype(int),
    })
    extra_html = "<h2>Figure 12 Denominator Table</h2>\n" + _html_table(denom_table)

    n_runs = int(display_counts[["source", "run"]].drop_duplicates().shape[0])
    footnotes = [
        "Cells show the percentage distribution of active antibiotic exposure among people "
        "infected with each bacterium during the 2022-2025 calibration window. The numerator "
        "is summed currently_on_drug_by_bacteria_drug for each bacterium-drug pair, aggregated "
        "to drug class. This is not necessarily the drug prescribed specifically for that "
        "bacterium; it includes empiric therapy, targeted therapy, combination therapy, and "
        "bystander exposure from antibiotics used for other concurrent infections.",
        "Rows are normalised within bacterium and sum to 100% across displayed classes plus "
        "'Other classes'.",
        f"Values are median within-bacterium class shares across {n_runs} simulation run"
        f"{'s' if n_runs > 1 else ''}; exposure counts are summed first within each run, "
        "then converted to percentages.",
        "Drug-to-class mapping source: " + class_mapping_source + ".",
    ]
    if detections:
        footnotes.append("Exposure columns detected as: " + "; ".join(sorted(detections)) + ".")
    if activity_sort:
        footnotes.append("Rows are sorted by Figure 11 retained activity, lowest first.")
    else:
        footnotes.append("Rows are sorted by median antibiotic-exposure denominator, highest first.")
    if unmapped:
        footnotes.append("Unmapped drugs assigned to 'Other / unmapped': " + ", ".join(sorted(unmapped)) + ".")
    if len(class_totals) > len(top_classes):
        collapsed = [cls for cls in class_totals.index if cls not in top_classes]
        footnotes.append("Collapsed into 'Other classes': " + ", ".join(collapsed) + ".")

    _save_figure(
        fig,
        out_dir,
        _F19_STEM,
        _F19_TITLE,
        "Heatmap of antibiotic exposure distribution among infected people, normalised "
        "within each bacterium.",
        footnotes,
        agg=agg,
        extra_html=extra_html,
    )


# ---------------------------------------------------------------------------
# Supplementary Table S1. Infection outcomes by bacterium
# ---------------------------------------------------------------------------

_ST1_TITLE = "Supplementary Table S1. Infection outcomes by bacterium, 2022\u20132025"
_ST1_STEM = "Supplementary_Table_S1__infection_outcomes_by_bacterium"
_S8_TITLE = "Supplementary Figure S8. Infection outcome pathway by bacterium, 2022\u20132025"
_S8_STEM = "Supplementary_Figure_S8__infection_outcome_pathway_by_bacterium"
_ST1_REQUIRED_VECTOR_COLUMNS = [
    "new_active_infections_by_bacteria",
    "active_infection_days_by_bacteria",
    "treated_infection_days_by_bacteria",
    "effective_treated_infection_days_by_bacteria",
    "infection_resolution_count_by_bacteria",
    "sepsis_onset_count_by_bacteria",
    "infection_death_count_by_bacteria",
]
_ST1_OPTIONAL_VECTOR_COLUMNS = [
    "drug_failure_count_by_bacteria",
]
_ST1_VECTOR_COLUMNS = _ST1_REQUIRED_VECTOR_COLUMNS + _ST1_OPTIONAL_VECTOR_COLUMNS
_ST1_REQUIRED_MESSAGE = (
    "Supplementary Figure S8 requires aggregate per-bacterium vector columns in "
    "simulation_summary_*.csv: " + ", ".join(_ST1_REQUIRED_VECTOR_COLUMNS) + ". "
    "The optional drug_failure_count_by_bacteria column is used for the treatment-failure panel "
    "when available. Re-run the Rust simulation with the required aggregate observability "
    "fields to generate the real pathway figure."
)


def _supplementary_table_s1_placeholder(
    out_dir: Path,
    agg: dict | None,
    message: str,
    problems: list[str] | None = None,
) -> None:
    path = out_dir / TABLES_DIRNAME / f"{_ST1_STEM}.html"
    body = _html_head(_ST1_TITLE)
    body += _back_link()
    body += f"<h1>{_ST1_TITLE}</h1>\n"
    if agg is not None:
        body += _meta_box(agg)
    body += f"<p class='note'>{message}</p>\n"
    body += "<h2>Required simulation_summary columns</h2>\n<ul>\n"
    for column in _ST1_VECTOR_COLUMNS:
        body += f"<li><code>{column}</code></li>\n"
    body += "</ul>\n"
    if problems:
        body += "<h2>Parser notes</h2>\n<ul>\n"
        for problem in problems[:8]:
            body += f"<li>{problem}</li>\n"
        body += "</ul>\n"
    body += _html_footnotes([
        "Old simulation_summary files remain readable: when the required aggregate columns are "
        "absent, this placeholder is generated rather than failing the paper-output build.",
        "No supplementary-table-specific CSV or event-level output is required.",
    ])
    body += "</body></html>"
    _save(path, body)


def _st1_sum_vector_column(run_df: pd.DataFrame, column: str, target_len: int) -> np.ndarray:
    totals = np.zeros(target_len, dtype=float)
    for value in run_df[column]:
        values = np.array(_figure_15_parse_vector_cell(value), dtype=float)
        if len(values) > len(totals):
            totals = _figure_15_extend_array(totals, len(values))
        values = _figure_15_extend_array(values, len(totals))
        totals += np.nan_to_num(values, nan=0.0)
    return totals


def _st1_rate_per_1000(numerator: float, denominator: float) -> float:
    if not np.isfinite(denominator) or denominator <= 0.0:
        return np.nan
    return 1000.0 * numerator / denominator


def _st1_percent(numerator: float, denominator: float) -> float:
    if not np.isfinite(denominator) or denominator <= 0.0:
        return np.nan
    return 100.0 * numerator / denominator


def _st1_reliability_flags(
    new_infections: float,
    treated_days: float,
    sepsis_onsets: float,
    deaths: float,
) -> str:
    flags: list[str] = []
    if not np.isfinite(new_infections) or new_infections < 100:
        flags.append("low infection denominator")
    if not np.isfinite(treated_days) or treated_days < 100:
        flags.append("low treated denominator")
    if not np.isfinite(sepsis_onsets) or sepsis_onsets < 20:
        flags.append("low sepsis count")
    if not np.isfinite(deaths) or deaths < 20:
        flags.append("low death count")
    if flags:
        flags.append("unstable rate")
    return "; ".join(flags)


def _st1_rows_from_simulation_csv(csv_path: Path) -> tuple[list[dict[str, object]], str | None]:
    header = _simulation_csv_column_names(csv_path)
    if header is None:
        return [], f"{csv_path.name}: unable to read simulation summary CSV header."

    available = set(header)
    missing = [column for column in _ST1_REQUIRED_VECTOR_COLUMNS if column not in available]
    if missing:
        return [], f"{csv_path.name}: missing columns {', '.join(missing)}."
    drug_failure_available = "drug_failure_count_by_bacteria" in available

    optional = ["policy_option", "run_id", "simulation_year", "year", "time_in_years", "time_step"]
    usecols = [
        column
        for column in _ST1_REQUIRED_VECTOR_COLUMNS + _ST1_OPTIONAL_VECTOR_COLUMNS + optional
        if column in available
    ]
    try:
        df = _read_csv_selected(csv_path, usecols)
    except (ValueError, OSError) as exc:
        return [], f"{csv_path.name}: unable to load Supplementary Table S1 columns ({exc})."

    if "policy_option" in df.columns:
        df = df[pd.to_numeric(df["policy_option"], errors="coerce") == 0].copy()
    df["st1_year"] = _simulation_year_series(df)
    df = df[(df["st1_year"] >= 2022.0) & (df["st1_year"] < 2026.0)].copy()
    if df.empty:
        return [], f"{csv_path.name}: no baseline-policy rows in 2022-2025."

    grouped = df.groupby("run_id", dropna=False) if "run_id" in df.columns else [(csv_path.stem, df)]
    target_len = len(_F15_KNOWN_BACTERIA_SLUGS)
    rows: list[dict[str, object]] = []

    for run_key, run_df in grouped:
        totals = {
            column: _st1_sum_vector_column(run_df, column, target_len)
            for column in _ST1_REQUIRED_VECTOR_COLUMNS
        }
        if drug_failure_available:
            totals["drug_failure_count_by_bacteria"] = _st1_sum_vector_column(
                run_df,
                "drug_failure_count_by_bacteria",
                target_len,
            )
        else:
            totals["drug_failure_count_by_bacteria"] = np.full(target_len, np.nan)
        run_len = max(len(values) for values in totals.values())
        for column in totals:
            totals[column] = _figure_15_extend_array(totals[column], run_len)

        for b_idx in range(run_len):
            slug = (
                _F15_KNOWN_BACTERIA_SLUGS[b_idx]
                if b_idx < len(_F15_KNOWN_BACTERIA_SLUGS)
                else f"bacterium_{b_idx + 1}"
            )
            new_infections = float(totals["new_active_infections_by_bacteria"][b_idx])
            active_days = float(totals["active_infection_days_by_bacteria"][b_idx])
            treated_days = float(totals["treated_infection_days_by_bacteria"][b_idx])
            effective_days = float(totals["effective_treated_infection_days_by_bacteria"][b_idx])
            resolutions = float(totals["infection_resolution_count_by_bacteria"][b_idx])
            sepsis_onsets = float(totals["sepsis_onset_count_by_bacteria"][b_idx])
            deaths = float(totals["infection_death_count_by_bacteria"][b_idx])
            drug_failures = float(totals["drug_failure_count_by_bacteria"][b_idx])

            rows.append({
                "source": csv_path.name,
                "run": str(run_key),
                "bacterium_slug": slug,
                "bacterium": _figure_15_bacterium_label(slug),
                "new_active_infections": new_infections,
                "active_infection_days": active_days,
                "treated_infection_days": treated_days,
                "treated_infection_percent": _st1_percent(treated_days, active_days),
                "effective_treated_infection_days": effective_days,
                "effective_treated_percent": _st1_percent(effective_days, treated_days),
                "infection_resolutions": resolutions,
                "resolution_rate_per_1000": _st1_rate_per_1000(resolutions, new_infections),
                "drug_failure_events": drug_failures,
                "drug_failure_available": drug_failure_available,
                "drug_failure_rate_per_1000_treated_days": _st1_rate_per_1000(
                    drug_failures,
                    treated_days,
                ),
                "sepsis_onsets": sepsis_onsets,
                "sepsis_rate_per_1000": _st1_rate_per_1000(sepsis_onsets, new_infections),
                "infection_deaths": deaths,
                "fatality_per_1000": _st1_rate_per_1000(deaths, new_infections),
            })

    if not rows:
        return [], f"{csv_path.name}: Supplementary Table S1 columns contained no usable values."
    return rows, None


def _st1_format_count(value: object) -> str:
    if value is None or pd.isna(value):
        return "\u2014"
    value_f = float(value)
    if not np.isfinite(value_f):
        return "\u2014"
    return f"{int(np.rint(value_f)):,}"


def _st1_format_percent(value: object) -> str:
    if value is None or pd.isna(value):
        return "\u2014"
    value_f = float(value)
    if not np.isfinite(value_f):
        return "\u2014"
    return f"{value_f:.1f}"


def _st1_format_rate(value: object) -> str:
    if value is None or pd.isna(value):
        return "\u2014"
    value_f = float(value)
    if not np.isfinite(value_f):
        return "\u2014"
    return f"{value_f:.2f}" if abs(value_f) < 10 else f"{value_f:.1f}"


def _st1_build_summary(rows: list[dict[str, object]]) -> tuple[pd.DataFrame, int]:
    if not rows:
        return pd.DataFrame(), 0
    df = pd.DataFrame(rows)
    n_runs = int(df[["source", "run"]].drop_duplicates().shape[0])
    summary = (
        df.groupby(["bacterium_slug", "bacterium"], as_index=False)
        .agg(
            new_active_infections=("new_active_infections", "median"),
            active_infection_days=("active_infection_days", "median"),
            treated_infection_days=("treated_infection_days", "median"),
            treated_infection_percent=("treated_infection_percent", "median"),
            effective_treated_infection_days=("effective_treated_infection_days", "median"),
            effective_treated_percent=("effective_treated_percent", "median"),
            infection_resolutions=("infection_resolutions", "median"),
            resolution_rate_per_1000=("resolution_rate_per_1000", "median"),
            drug_failure_events=("drug_failure_events", "median"),
            drug_failure_available=("drug_failure_available", "any"),
            drug_failure_rate_per_1000_treated_days=(
                "drug_failure_rate_per_1000_treated_days",
                "median",
            ),
            sepsis_onsets=("sepsis_onsets", "median"),
            sepsis_rate_per_1000=("sepsis_rate_per_1000", "median"),
            infection_deaths=("infection_deaths", "median"),
            fatality_per_1000=("fatality_per_1000", "median"),
        )
    )
    summary["reliability_flag"] = summary.apply(
        lambda row: _st1_reliability_flags(
            float(row["new_active_infections"]),
            float(row["treated_infection_days"]),
            float(row["sepsis_onsets"]),
            float(row["infection_deaths"]),
        ),
        axis=1,
    )
    summary = summary.sort_values(
        ["infection_deaths", "new_active_infections"],
        ascending=[False, False],
    ).reset_index(drop=True)
    return summary, n_runs


def _st1_detail_table_from_summary(summary: pd.DataFrame) -> pd.DataFrame:
    return pd.DataFrame({
        "Bacterium": summary["bacterium"],
        "New active infections": summary["new_active_infections"].map(_st1_format_count),
        "Active infection-days": summary["active_infection_days"].map(_st1_format_count),
        "Treated infection-days": summary["treated_infection_days"].map(_st1_format_count),
        "Treated infection-days (% of active)": summary["treated_infection_percent"].map(_st1_format_percent),
        "Effective treated infection-days": summary["effective_treated_infection_days"].map(_st1_format_count),
        "Effective treated infection-days (% of treated)": summary["effective_treated_percent"].map(_st1_format_percent),
        "Infection resolutions": summary["infection_resolutions"].map(_st1_format_count),
        "Resolution rate per 1,000 new active infections": summary["resolution_rate_per_1000"].map(_st1_format_rate),
        "Drug / treatment failure events": summary["drug_failure_events"].map(_st1_format_count),
        "Failure rate per 1,000 treated infection-days": summary[
            "drug_failure_rate_per_1000_treated_days"
        ].map(_st1_format_rate),
        "Sepsis onsets": summary["sepsis_onsets"].map(_st1_format_count),
        "Sepsis onset rate per 1,000 new active infections": summary["sepsis_rate_per_1000"].map(_st1_format_rate),
        "Infection deaths": summary["infection_deaths"].map(_st1_format_count),
        "Infection fatality per 1,000 new active infections": summary["fatality_per_1000"].map(_st1_format_rate),
        "Reliability flag": summary["reliability_flag"],
    })


def _s8_reliability_legend_html() -> str:
    legend = pd.DataFrame({
        "Reliability flag": [
            "low infection denominator",
            "low treated denominator",
            "low sepsis count",
            "low death count",
            "unstable rate",
        ],
        "Meaning": [
            "Fewer than 100 new active infections in the run-level median denominator.",
            "Fewer than 100 treated infection-days in the run-level median denominator.",
            "Fewer than 20 sepsis-onset events in the run-level median numerator.",
            "Fewer than 20 infection-death events in the run-level median numerator.",
            "At least one small-count flag applies; pathway rates should be interpreted cautiously.",
        ],
    })
    return "<h2>Reliability flag legend</h2>\n" + _html_table(legend)


def _s8_placeholder(
    out_dir: Path,
    agg: dict | None,
    message: str,
    problems: list[str] | None = None,
) -> dict[str, object]:
    fig, ax = plt.subplots(figsize=(10.5, 3.8))
    ax.text(
        0.5,
        0.5,
        f"{_S8_TITLE}\n\n{message}",
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=10.2,
        color="#555",
        bbox=dict(boxstyle="round,pad=0.6", fc="#f5f5f5", ec="#bbb"),
        wrap=True,
    )
    ax.set_axis_off()
    fig.subplots_adjust(left=0.03, right=0.97, top=0.88, bottom=0.08)

    required = pd.DataFrame({
        "Required simulation_summary field": _ST1_REQUIRED_VECTOR_COLUMNS,
    })
    extra_html = (
        "<p class='note'>This placeholder replaces the prior standalone Supplementary Table S1 "
        "layout. A model run with the required aggregate simulation_summary fields is needed "
        "to render the pathway figure and detailed table.</p>\n"
        "<h2>Required fields</h2>\n"
        + _html_table(required)
    )
    if problems:
        extra_html += "<h2>Parser notes</h2>\n<ul>\n"
        for problem in problems[:10]:
            extra_html += f"<li>{problem}</li>\n"
        extra_html += "</ul>\n"

    _save_figure(
        fig,
        out_dir,
        _S8_STEM,
        _S8_TITLE,
        message,
        [
            "Old simulation_summary files remain readable: when required aggregate columns are absent, "
            "this placeholder is generated rather than failing the paper-output build.",
            "No new model-output CSV file is required or produced.",
        ],
        agg=agg,
        extra_html=extra_html,
    )
    return {"generated": "placeholder", "bacteria_included": 0, "n_runs": 0}


def _s8_apply_count_axis(ax: "plt.Axes", values: pd.Series | np.ndarray) -> None:
    finite = [
        float(value)
        for value in pd.Series(values).dropna().to_numpy(dtype=float)
        if np.isfinite(float(value)) and float(value) > 0.0
    ]
    if not finite:
        ax.set_xlim(0, 1)
        return
    max_value = max(finite)
    min_value = min(finite)
    if max_value / max(min_value, 1.0) >= 100.0:
        ax.set_xscale("symlog", linthresh=max(1.0, min_value))
    ax.set_xlim(0, max_value * 1.18)


def _s8_plot_rate_panel(
    ax: "plt.Axes",
    y: np.ndarray,
    values: pd.Series,
    title: str,
    color: str,
    xlabel: str,
) -> None:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    mask = np.isfinite(numeric)
    if mask.any():
        ax.scatter(numeric[mask], y[mask], s=20, color=color, edgecolors="white", linewidths=0.35, zorder=3)
        ax.set_xlim(0, max(1.0, float(np.nanmax(numeric[mask])) * 1.15))
    else:
        ax.text(
            0.5,
            0.5,
            "No usable values",
            ha="center",
            va="center",
            transform=ax.transAxes,
            fontsize=8.2,
            color="#555",
        )
        ax.set_xlim(0, 1)
    ax.set_title(title, loc="left", fontsize=9.2, fontweight="bold")
    ax.set_xlabel(xlabel, fontsize=8.1)
    ax.grid(axis="x", linewidth=0.35, alpha=0.45)
    ax.spines[["top", "right"]].set_visible(False)


def make_supplementary_figure_s8_infection_outcome_pathway(
    csv_paths: list[Path],
    out_dir: Path,
    agg: dict | None = None,
) -> dict[str, object]:
    if not csv_paths:
        print("  Supplementary Figure S8: placeholder (no simulation summary CSVs found).")
        return _s8_placeholder(
            out_dir,
            agg,
            _ST1_REQUIRED_MESSAGE,
            ["No matching simulation_summary_*.csv files were found."],
        )

    rows: list[dict[str, object]] = []
    problems: list[str] = []
    for csv_path in csv_paths:
        run_rows, problem = _st1_rows_from_simulation_csv(csv_path)
        rows.extend(run_rows)
        if problem:
            problems.append(problem)

    if not rows:
        print("  Supplementary Figure S8: placeholder (no usable ST1/S8 rows).")
        return _s8_placeholder(out_dir, agg, _ST1_REQUIRED_MESSAGE, problems)

    summary, n_runs = _st1_build_summary(rows)
    if summary.empty:
        print("  Supplementary Figure S8: placeholder (empty pathway summary).")
        return _s8_placeholder(out_dir, agg, _ST1_REQUIRED_MESSAGE, problems)

    table = _st1_detail_table_from_summary(summary)
    has_drug_failure = bool(summary["drug_failure_available"].any())
    plot = summary.reset_index(drop=True)
    y = np.arange(len(plot))
    fig_height = max(8.5, 2.2 + 0.25 * len(plot))
    fig, axes = plt.subplots(
        1,
        5,
        figsize=(18.2, fig_height),
        sharey=True,
        gridspec_kw={"width_ratios": [1.1, 1.35, 1.05, 1.0, 1.0], "wspace": 0.1},
    )

    axes[0].barh(y, plot["new_active_infections"].to_numpy(dtype=float), height=0.58, color="#4E79A7", alpha=0.85)
    axes[0].set_title("A. New active infections", loc="left", fontsize=9.2, fontweight="bold")
    axes[0].set_xlabel("Raw simulated n", fontsize=8.1)
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(plot["bacterium"].values, fontsize=6.5, fontstyle="italic")
    axes[0].invert_yaxis()
    _s8_apply_count_axis(axes[0], plot["new_active_infections"])
    axes[0].grid(axis="x", linewidth=0.35, alpha=0.45)
    axes[0].spines[["top", "right"]].set_visible(False)

    treated_pct = plot["treated_infection_percent"].to_numpy(dtype=float)
    effective_pct = plot["effective_treated_percent"].to_numpy(dtype=float)
    for idx, (treated, effective) in enumerate(zip(treated_pct, effective_pct)):
        if np.isfinite(treated) and np.isfinite(effective):
            axes[1].plot([treated, effective], [idx, idx], color="#b0bec5", linewidth=0.75, zorder=1)
    treated_mask = np.isfinite(treated_pct)
    effective_mask = np.isfinite(effective_pct)
    axes[1].scatter(
        treated_pct[treated_mask],
        y[treated_mask] - 0.09,
        s=18,
        color="#59A14F",
        edgecolors="white",
        linewidths=0.3,
        label="Treated / active",
        zorder=3,
    )
    axes[1].scatter(
        effective_pct[effective_mask],
        y[effective_mask] + 0.09,
        s=18,
        color="#F28E2B",
        edgecolors="white",
        linewidths=0.3,
        label="Effective / treated",
        zorder=3,
    )
    max_pct = np.nanmax([np.nanmax(treated_pct) if treated_mask.any() else 0, np.nanmax(effective_pct) if effective_mask.any() else 0])
    axes[1].set_xlim(0, max(100.0, float(max_pct) * 1.08))
    axes[1].set_title("B. Treatment exposure and effectiveness", loc="left", fontsize=9.2, fontweight="bold")
    axes[1].set_xlabel("% of denominator days", fontsize=8.1)
    axes[1].grid(axis="x", linewidth=0.35, alpha=0.45)
    axes[1].legend(fontsize=6.8, frameon=False, loc="lower right")
    axes[1].spines[["top", "right"]].set_visible(False)

    if has_drug_failure:
        _s8_plot_rate_panel(
            axes[2],
            y,
            plot["drug_failure_rate_per_1000_treated_days"],
            "C. Treatment failure",
            "#E15759",
            "Events per 1,000 treated days",
        )
    else:
        axes[2].text(
            0.5,
            0.5,
            "drug_failure_count_by_bacteria\nwas not available",
            ha="center",
            va="center",
            transform=axes[2].transAxes,
            fontsize=8.2,
            color="#555",
            bbox=dict(boxstyle="round,pad=0.4", fc="#f5f5f5", ec="#bbb"),
        )
        axes[2].set_title("C. Treatment failure", loc="left", fontsize=9.2, fontweight="bold")
        axes[2].set_xlabel("Unavailable", fontsize=8.1)
        axes[2].set_xlim(0, 1)
        axes[2].spines[["top", "right"]].set_visible(False)

    _s8_plot_rate_panel(
        axes[3],
        y,
        plot["sepsis_rate_per_1000"],
        "D. Sepsis onset",
        "#B07AA1",
        "Onsets per 1,000 infections",
    )
    _s8_plot_rate_panel(
        axes[4],
        y,
        plot["fatality_per_1000"],
        "E. Infection fatality",
        "#7F3C8D",
        "Deaths per 1,000 infections",
    )
    for ax in axes[1:]:
        ax.tick_params(axis="y", left=False, labelleft=False)
    fig.suptitle(_S8_TITLE, fontsize=11, fontweight="bold", y=0.995)
    fig.subplots_adjust(left=0.16, right=0.995, top=0.955, bottom=0.065, wspace=0.12)

    extra_html = (
        "<p class='note'><strong>Layout note:</strong> This figure replaces the prior standalone "
        "Supplementary Table S1 layout. The detailed numeric Supplementary Table S1 content is "
        "preserved below the multipanel figure.</p>\n"
        "<p class='note'>Panel A provides denominator context for the pathway rates in this figure. "
        "Supplementary Figure S6 focuses specifically on serious-R denominators, and Supplementary "
        "Figure S7 compares active infection incidence with calibration targets.</p>\n"
        "<h2>Detailed infection outcome table</h2>\n"
        + _html_table(table)
        + _s8_reliability_legend_html()
    )
    note = (
        "Integrated pathway view of baseline-policy infection outcomes by bacterium; rows are sorted "
        "by median infection deaths, then median new active infections."
    )
    footnotes = [
        "Values summarise baseline-policy infection outcomes by bacterium in the 2022-2025 window.",
        "Counts are raw simulated counts unless otherwise stated.",
        "Rates are calculated from summed run-level numerators and denominators over the window.",
        "Treatment exposure is counted while a person is infected with the bacterium and may include empiric therapy, targeted therapy, combination therapy, and bystander antibiotic exposure.",
        "Effective treated infection-days use the same resistance-adjusted activity threshold as the sepsis effective-therapy summary, expected activity_r >= 0.500.",
        "Infection resolutions count active infection episodes ending by clearance/resolution rather than death.",
        "Sepsis onsets use the existing simulation_summary aggregate semantics.",
        "Deaths are counted per pathogen involved; in polymicrobial cases, the same death can appear under more than one bacterium.",
        "Drug / treatment failure events use the existing model-defined failure-event semantics.",
        "Reliability flags mark small denominators and unstable rates.",
        f"Metrics are calculated within each run first and displayed as medians across {n_runs} accepted simulation run{'s' if n_runs != 1 else ''}.",
        "Source: existing simulation_summary_*.csv aggregate columns only. No event-level or supplementary-output CSV file is required or produced.",
    ]
    if not has_drug_failure:
        footnotes.append("drug_failure_count_by_bacteria was unavailable, so Panel C and failure-rate table cells are placeholders.")
    if problems:
        footnotes.append("Some simulation summary CSVs were skipped: " + " ".join(problems[:3]))

    _save_figure(
        fig,
        out_dir,
        _S8_STEM,
        _S8_TITLE,
        note,
        footnotes,
        agg=agg,
        extra_html=extra_html,
    )
    print(
        "  Supplementary Figure S8: real data; "
        f"{len(summary)} bacterium row{'s' if len(summary) != 1 else ''} from "
        f"{n_runs} run{'s' if n_runs != 1 else ''}; "
        f"drug-failure panel {'available' if has_drug_failure else 'placeholder'}."
    )
    return {
        "generated": "real data",
        "bacteria_included": int(len(summary)),
        "n_runs": n_runs,
        "drug_failure_available": has_drug_failure,
    }


def make_supplementary_table_s1_infection_outcomes(
    csv_paths: list[Path],
    out_dir: Path,
    agg: dict | None = None,
) -> None:
    if not csv_paths:
        _supplementary_table_s1_placeholder(
            out_dir,
            agg,
            _ST1_REQUIRED_MESSAGE,
            ["No matching simulation_summary_*.csv files were found."],
        )
        print("  Supplementary Table S1: placeholder (no simulation summary CSVs found).")
        return

    rows: list[dict[str, object]] = []
    problems: list[str] = []
    for csv_path in csv_paths:
        run_rows, problem = _st1_rows_from_simulation_csv(csv_path)
        rows.extend(run_rows)
        if problem:
            problems.append(problem)

    if not rows:
        _supplementary_table_s1_placeholder(out_dir, agg, _ST1_REQUIRED_MESSAGE, problems)
        if problems:
            print("  Supplementary Table S1: placeholder; " + " ".join(problems[:3]))
        else:
            print("  Supplementary Table S1: placeholder (no usable rows).")
        return

    df = pd.DataFrame(rows)
    n_runs = int(df[["source", "run"]].drop_duplicates().shape[0])
    summary = (
        df.groupby(["bacterium_slug", "bacterium"], as_index=False)
        .agg(
            new_active_infections=("new_active_infections", "median"),
            active_infection_days=("active_infection_days", "median"),
            treated_infection_days=("treated_infection_days", "median"),
            effective_treated_infection_days=("effective_treated_infection_days", "median"),
            effective_treated_percent=("effective_treated_percent", "median"),
            infection_resolutions=("infection_resolutions", "median"),
            resolution_rate_per_1000=("resolution_rate_per_1000", "median"),
            drug_failure_events=("drug_failure_events", "median"),
            drug_failure_rate_per_1000_treated_days=(
                "drug_failure_rate_per_1000_treated_days",
                "median",
            ),
            sepsis_onsets=("sepsis_onsets", "median"),
            sepsis_rate_per_1000=("sepsis_rate_per_1000", "median"),
            infection_deaths=("infection_deaths", "median"),
            fatality_per_1000=("fatality_per_1000", "median"),
        )
    )
    summary["reliability_flag"] = summary.apply(
        lambda row: _st1_reliability_flags(
            float(row["new_active_infections"]),
            float(row["treated_infection_days"]),
            float(row["sepsis_onsets"]),
            float(row["infection_deaths"]),
        ),
        axis=1,
    )
    summary = summary.sort_values(
        ["infection_deaths", "new_active_infections"],
        ascending=[False, False],
    ).reset_index(drop=True)

    table = pd.DataFrame({
        "Bacterium": summary["bacterium"],
        "New active infections": summary["new_active_infections"].map(_st1_format_count),
        "Active infection-days": summary["active_infection_days"].map(_st1_format_count),
        "Treated infection-days": summary["treated_infection_days"].map(_st1_format_count),
        "Effective treated infection-days": summary["effective_treated_infection_days"].map(_st1_format_count),
        "Effective treated infection-days (% of treated)": summary["effective_treated_percent"].map(_st1_format_percent),
        "Infection resolutions": summary["infection_resolutions"].map(_st1_format_count),
        "Resolution rate per 1,000 new active infections": summary["resolution_rate_per_1000"].map(_st1_format_rate),
        "Drug / treatment failure events": summary["drug_failure_events"].map(_st1_format_count),
        "Failure rate per 1,000 treated infection-days": summary[
            "drug_failure_rate_per_1000_treated_days"
        ].map(_st1_format_rate),
        "Sepsis onsets": summary["sepsis_onsets"].map(_st1_format_count),
        "Sepsis onset rate per 1,000 new active infections": summary["sepsis_rate_per_1000"].map(_st1_format_rate),
        "Infection deaths": summary["infection_deaths"].map(_st1_format_count),
        "Infection fatality per 1,000 new active infections": summary["fatality_per_1000"].map(_st1_format_rate),
        "Reliability flag": summary["reliability_flag"],
    })

    body = _html_head(_ST1_TITLE)
    body += _back_link()
    body += f"<h1>{_ST1_TITLE}</h1>\n"
    if agg is not None:
        body += _meta_box(agg)
    body += (
        "<p class='note'>Baseline-policy rows with 2022 <= simulation year < 2026. "
        f"Values are median run-level counts and rates across {n_runs} accepted simulation run"
        f"{'s' if n_runs != 1 else ''}; each run is summed across the window before rates are calculated.</p>\n"
    )
    body += _html_table(table)
    footnotes = [
        "Supplementary Table S1 summarises baseline-policy infection outcomes by bacterium in the 2022-2025 window. "
        "Counts are raw simulated counts unless otherwise stated. Rates are calculated from summed numerators and "
        "denominators over the window within each run before cross-run summarisation.",
        "Treatment exposure is counted while a person is infected with the bacterium and may include empiric therapy, "
        "targeted therapy, combination therapy, and bystander antibiotic exposure from other concurrent infections.",
        "Effective treated infection-days use the same resistance-adjusted activity threshold as Figure 10 / the sepsis "
        "effective-therapy summary: activity_r >= 0.500.",
        "Infection resolutions count immune-clearance plus drug-assisted-clearance resolution records, excluding death "
        "as the episode-ending event.",
        "Sepsis onsets use the same per-bacterium incident sepsis semantics as the sepsis effective-therapy summary: "
        "a bacterium contributes on the timestep its sepsis episode begins while the infection is active.",
        "Deaths are counted per pathogen involved. In polymicrobial cases, the same death can appear under more than "
        "one bacterium, so the sum across bacteria can exceed the headline infection-death total.",
        "Drug / treatment failure events are the model's existing day-5 post-drug-initiation failure events, aggregated "
        "to bacterium. They are not inferred from ineffective treated infection-days.",
        "Reliability flags mark rows with small denominators; rates in those rows should be interpreted cautiously.",
        "Source: existing simulation_summary_run#.csv aggregate columns only. No event-level CSV is required or produced.",
        _meta_footnote(agg) if agg is not None else "",
    ]
    if n_runs == 1:
        footnotes.append("Single run: no cross-run interval is shown.")
    else:
        footnotes.append("Intervals are omitted to keep the supplementary table readable; displayed values are medians.")
    if problems:
        footnotes.append("Some simulation summary CSVs were skipped: " + " ".join(problems[:3]))
    body += _html_footnotes([note for note in footnotes if note])
    body += "</body></html>"
    _save(out_dir / TABLES_DIRNAME / f"{_ST1_STEM}.html", body)
    print(
        "  Supplementary Table S1: "
        f"{len(table)} bacterium row{'s' if len(table) != 1 else ''} from "
        f"{n_runs} run{'s' if n_runs != 1 else ''}."
    )


# ---------------------------------------------------------------------------
# Supplementary Table S2. Detailed bacterium-drug resistance benchmarks
# ---------------------------------------------------------------------------

_ST2_TITLE = "Supplementary Table S2. Detailed bacterium-drug resistance benchmarks, 2022\u20132025"
_ST2_STEM = "Supplementary_Table_S2__detailed_bacterium_drug_resistance_benchmarks"
_ST2_REQUIRED_COLUMNS = [
    "Bacteria",
    "Drug",
    "Class",
    "Inf sim (%)",
    "Inf target (%)",
    "Avg sim (%)",
    "Avg target (%)",
    "Micro sim (%)",
    "Inf days",
    "Res days",
    "Carrier days",
    "Flags",
]
_ST2_PLACEHOLDER_MESSAGE = (
    "Supplementary Table S2 requires the 'Resistance Benchmarks (percent resistant)' table "
    "in calibration_summary_*.txt, including Bacteria, Drug, Class, Inf sim (%), Inf target (%), "
    "Avg sim (%), Avg target (%), Micro sim (%), Inf days, Res days, Carrier days, and Flags."
)
_ST2_MISSING_SENTINELS = {"", "-", "---", "\u2014", "NaN", "nan", "N/A", "None", "none"}


def _st2_numeric(value: object) -> float:
    if value is None:
        return np.nan
    if isinstance(value, (int, float, np.integer, np.floating)):
        value_f = float(value)
        return value_f if np.isfinite(value_f) else np.nan
    text = str(value).strip()
    if text in _ST2_MISSING_SENTINELS:
        return np.nan
    text = text.replace(",", "").replace("%", "")
    try:
        value_f = float(text)
    except ValueError:
        return np.nan
    return value_f if np.isfinite(value_f) else np.nan


def _st2_text(value: object) -> str:
    if value is None:
        return ""
    if isinstance(value, float) and np.isnan(value):
        return ""
    text = str(value).strip()
    if text in _ST2_MISSING_SENTINELS:
        return ""
    return text


def _st2_numeric_values(values: pd.Series) -> np.ndarray:
    arr = np.array([_st2_numeric(value) for value in values], dtype=float)
    return arr[np.isfinite(arr)]


def _st2_first_nonmissing(values: pd.Series) -> float:
    arr = _st2_numeric_values(values)
    if arr.size == 0:
        return np.nan
    return float(arr[0])


def _st2_target_inconsistent(values: pd.Series) -> bool:
    arr = _st2_numeric_values(values)
    if arr.size <= 1:
        return False
    return bool(np.nanmax(arr) - np.nanmin(arr) > 1e-9)


def _st2_median(values: pd.Series) -> float:
    arr = _st2_numeric_values(values)
    if arr.size == 0:
        return np.nan
    return float(np.nanmedian(arr))


def _st2_percentile(values: pd.Series, q: float) -> float:
    arr = _st2_numeric_values(values)
    if arr.size == 0:
        return np.nan
    return float(np.nanpercentile(arr, q))


def _st2_delta(sim_value: float, target_value: float) -> float:
    if not np.isfinite(sim_value) or not np.isfinite(target_value):
        return np.nan
    return float(sim_value - target_value)


def _st2_format_number(value: object, decimals: int = 2) -> str:
    value_f = _st2_numeric(value)
    if not np.isfinite(value_f):
        return "\u2014"
    if abs(value_f) >= 10_000:
        return f"{value_f:,.0f}"
    return f"{value_f:.{decimals}f}"


def _st2_format_count(value: object) -> str:
    value_f = _st2_numeric(value)
    if not np.isfinite(value_f):
        return "\u2014"
    return f"{int(np.rint(value_f)):,}"


def _st2_format_interval(lo: object, hi: object) -> str:
    lo_f = _st2_numeric(lo)
    hi_f = _st2_numeric(hi)
    if not np.isfinite(lo_f) or not np.isfinite(hi_f):
        return "\u2014"
    return f"{lo_f:.2f}\u2013{hi_f:.2f}"


def _st2_split_flags(flags: object) -> list[str]:
    text = _st2_text(flags)
    if not text:
        return []
    parts = [part.strip() for part in re.split(r";|\|", text) if part.strip()]
    return parts if parts else [text]


def _st2_combine_flags(values: pd.Series) -> str:
    seen: set[str] = set()
    ordered: list[str] = []
    for value in values:
        for flag in _st2_split_flags(value):
            if flag not in seen:
                seen.add(flag)
                ordered.append(flag)
    return "; ".join(ordered)


def _st2_derived_flags(row: pd.Series) -> str:
    flags: list[str] = []
    original = str(row.get("original_flags", "")).lower()
    inf_days = _st2_numeric(row.get("inf_days_median"))
    res_days = _st2_numeric(row.get("res_days_median"))
    carrier_days = _st2_numeric(row.get("carrier_days_median"))

    if "negligible potency" in original:
        flags.append("negligible potency")
    if not np.isfinite(_st2_numeric(row.get("inf_target"))):
        flags.append("infection-resistance benchmark not assigned")
    if not np.isfinite(_st2_numeric(row.get("inf_sim_median"))):
        flags.append("missing infection-resistance simulation")
    if not np.isfinite(_st2_numeric(row.get("avg_target"))):
        flags.append("average resistant-level benchmark not assigned")
    if "no resistant infections" in original or (np.isfinite(res_days) and res_days <= 0):
        flags.append("no resistant infections for average metric")
    if "low resistant sample" in original:
        flags.append("low resistant sample")
    elif np.isfinite(res_days) and 0 < res_days < 50:
        flags.append("low resistant sample")
    if np.isfinite(inf_days) and inf_days < 100:
        flags.append("low infection-days")
    if not np.isfinite(carrier_days):
        flags.append("missing carrier-days")
        flags.append("missing microbiome denominator")
    elif carrier_days < 1000:
        flags.append("low carrier-days")
    if "expanded window" in original:
        flags.append("expanded window")
    if bool(row.get("target_inconsistent", False)):
        flags.append("benchmark value varied across runs")

    seen: set[str] = set()
    unique = [flag for flag in flags if not (flag in seen or seen.add(flag))]
    return "; ".join(unique) if unique else "no caveat"


def _st2_rows_from_runs(runs: list[dict]) -> tuple[pd.DataFrame, list[str], list[str]]:
    rows: list[dict[str, object]] = []
    problems: list[str] = []
    missing_columns: list[str] = []

    for run_idx, run in enumerate(runs, start=1):
        rb = run.get("resistance_benchmarks", pd.DataFrame())
        if rb is None or rb.empty:
            problems.append(f"run {run_idx}: missing Resistance Benchmarks table.")
            continue
        rb = rb.copy()
        missing = [column for column in _ST2_REQUIRED_COLUMNS if column not in rb.columns]
        if missing:
            problems.append(f"run {run_idx}: missing columns {', '.join(missing)}.")
            for column in missing:
                if column not in missing_columns:
                    missing_columns.append(column)
                rb[column] = np.nan if column not in {"Bacteria", "Drug", "Class", "Flags"} else ""

        for _, source_row in rb.iterrows():
            rows.append({
                "run": f"run {run_idx}",
                "Bacteria": _st2_text(source_row.get("Bacteria")),
                "Drug": _st2_text(source_row.get("Drug")),
                "Class": _st2_text(source_row.get("Class")),
                "Inf sim (%)": _st2_numeric(source_row.get("Inf sim (%)")),
                "Inf target (%)": _st2_numeric(source_row.get("Inf target (%)")),
                "Avg sim (%)": _st2_numeric(source_row.get("Avg sim (%)")),
                "Avg target (%)": _st2_numeric(source_row.get("Avg target (%)")),
                "Micro sim (%)": _st2_numeric(source_row.get("Micro sim (%)")),
                "Inf days": _st2_numeric(source_row.get("Inf days")),
                "Res days": _st2_numeric(source_row.get("Res days")),
                "Carrier days": _st2_numeric(source_row.get("Carrier days")),
                "Flags": _st2_text(source_row.get("Flags")),
            })

    return pd.DataFrame(rows), problems, missing_columns


def _st2_summarise_rows(raw: pd.DataFrame) -> pd.DataFrame:
    if raw.empty:
        return pd.DataFrame()

    records: list[dict[str, object]] = []
    for (bacterium, drug, drug_class), group in raw.groupby(
        ["Bacteria", "Drug", "Class"],
        dropna=False,
        sort=False,
    ):
        inf_target = _st2_first_nonmissing(group["Inf target (%)"])
        avg_target = _st2_first_nonmissing(group["Avg target (%)"])
        record = {
            "bacterium": str(bacterium),
            "drug": str(drug),
            "drug_class": str(drug_class),
            "inf_sim_median": _st2_median(group["Inf sim (%)"]),
            "inf_sim_p5": _st2_percentile(group["Inf sim (%)"], 5),
            "inf_sim_p95": _st2_percentile(group["Inf sim (%)"], 95),
            "inf_target": inf_target,
            "avg_sim_median": _st2_median(group["Avg sim (%)"]),
            "avg_sim_p5": _st2_percentile(group["Avg sim (%)"], 5),
            "avg_sim_p95": _st2_percentile(group["Avg sim (%)"], 95),
            "avg_target": avg_target,
            "micro_sim_median": _st2_median(group["Micro sim (%)"]),
            "micro_sim_p5": _st2_percentile(group["Micro sim (%)"], 5),
            "micro_sim_p95": _st2_percentile(group["Micro sim (%)"], 95),
            "inf_days_median": _st2_median(group["Inf days"]),
            "res_days_median": _st2_median(group["Res days"]),
            "carrier_days_median": _st2_median(group["Carrier days"]),
            "original_flags": _st2_combine_flags(group["Flags"]),
            "runs_contributing": int(group["run"].nunique()),
            "target_inconsistent": (
                _st2_target_inconsistent(group["Inf target (%)"])
                or _st2_target_inconsistent(group["Avg target (%)"])
            ),
        }
        record["infection_resistance_delta_pp"] = _st2_delta(
            float(record["inf_sim_median"]),
            inf_target,
        )
        record["abs_infection_resistance_delta_pp"] = (
            abs(float(record["infection_resistance_delta_pp"]))
            if np.isfinite(float(record["infection_resistance_delta_pp"]))
            else np.nan
        )
        record["avg_resistant_level_delta_pp"] = _st2_delta(
            float(record["avg_sim_median"]),
            avg_target,
        )
        records.append(record)

    summary = pd.DataFrame(records)
    if summary.empty:
        return summary
    summary["derived_flags"] = summary.apply(_st2_derived_flags, axis=1)
    return summary.sort_values(["bacterium", "drug_class", "drug"]).reset_index(drop=True)


def _st2_display_table(summary: pd.DataFrame, multiple_runs: bool) -> pd.DataFrame:
    if multiple_runs:
        return pd.DataFrame({
            "Bacterium": summary["bacterium"],
            "Drug": summary["drug"],
            "Drug class": summary["drug_class"],
            "Infection resistance simulation median (%)": summary["inf_sim_median"].map(_st2_format_number),
            "Infection resistance simulation 5th-95th percentile (%)": summary.apply(
                lambda row: _st2_format_interval(row["inf_sim_p5"], row["inf_sim_p95"]),
                axis=1,
            ),
            "Infection resistance calibration benchmark (%)": summary["inf_target"].map(_st2_format_number),
            "Infection resistance difference, median simulation minus benchmark (pp)": summary[
                "infection_resistance_delta_pp"
            ].map(_st2_format_number),
            "Absolute infection-resistance difference (pp)": summary[
                "abs_infection_resistance_delta_pp"
            ].map(_st2_format_number),
            "Average resistant level simulation median (%)": summary["avg_sim_median"].map(_st2_format_number),
            "Average resistant level simulation 5th-95th percentile (%)": summary.apply(
                lambda row: _st2_format_interval(row["avg_sim_p5"], row["avg_sim_p95"]),
                axis=1,
            ),
            "Average resistant level expert-assigned model benchmark (%)": summary["avg_target"].map(_st2_format_number),
            "Average resistant-level difference, median simulation minus benchmark (pp)": summary[
                "avg_resistant_level_delta_pp"
            ].map(_st2_format_number),
            "Microbiome/carriage resistance simulation median (%)": summary["micro_sim_median"].map(_st2_format_number),
            "Microbiome/carriage resistance simulation 5th-95th percentile (%)": summary.apply(
                lambda row: _st2_format_interval(row["micro_sim_p5"], row["micro_sim_p95"]),
                axis=1,
            ),
            "Median infection-days": summary["inf_days_median"].map(_st2_format_count),
            "Median resistant infection-days": summary["res_days_median"].map(_st2_format_count),
            "Median carrier-days": summary["carrier_days_median"].map(_st2_format_count),
            "Number of runs contributing": summary["runs_contributing"],
            "Original calibration flags": summary["original_flags"],
            "Derived reliability / interpretation flag": summary["derived_flags"],
        })

    return pd.DataFrame({
        "Bacterium": summary["bacterium"],
        "Drug": summary["drug"],
        "Drug class": summary["drug_class"],
        "Infection resistance, simulation (%)": summary["inf_sim_median"].map(_st2_format_number),
        "Infection resistance, calibration benchmark (%)": summary["inf_target"].map(_st2_format_number),
        "Infection resistance difference, simulation minus benchmark (pp)": summary[
            "infection_resistance_delta_pp"
        ].map(_st2_format_number),
        "Absolute infection-resistance difference (pp)": summary[
            "abs_infection_resistance_delta_pp"
        ].map(_st2_format_number),
        "Average resistant level, simulation (%)": summary["avg_sim_median"].map(_st2_format_number),
        "Average resistant level, expert-assigned model benchmark (%)": summary["avg_target"].map(_st2_format_number),
        "Average resistant-level difference, simulation minus benchmark (pp)": summary[
            "avg_resistant_level_delta_pp"
        ].map(_st2_format_number),
        "Microbiome/carriage resistance, simulation (%)": summary["micro_sim_median"].map(_st2_format_number),
        "Infection-days": summary["inf_days_median"].map(_st2_format_count),
        "Resistant infection-days": summary["res_days_median"].map(_st2_format_count),
        "Carrier-days": summary["carrier_days_median"].map(_st2_format_count),
        "Original calibration flags": summary["original_flags"],
        "Derived reliability / interpretation flag": summary["derived_flags"],
    })


def _st2_largest_differences_table(summary: pd.DataFrame) -> pd.DataFrame:
    eligible = summary.dropna(
        subset=[
            "inf_sim_median",
            "inf_target",
            "infection_resistance_delta_pp",
            "abs_infection_resistance_delta_pp",
        ]
    ).copy()
    if eligible.empty:
        return pd.DataFrame()
    eligible = eligible.sort_values(
        "abs_infection_resistance_delta_pp",
        ascending=False,
    ).head(20)
    return pd.DataFrame({
        "Rank": np.arange(1, len(eligible) + 1),
        "Bacterium": eligible["bacterium"],
        "Drug": eligible["drug"],
        "Drug class": eligible["drug_class"],
        "Infection resistance simulation (%)": eligible["inf_sim_median"].map(_st2_format_number),
        "Infection resistance calibration benchmark (%)": eligible["inf_target"].map(_st2_format_number),
        "Difference (pp)": eligible["infection_resistance_delta_pp"].map(_st2_format_number),
        "Absolute difference (pp)": eligible["abs_infection_resistance_delta_pp"].map(_st2_format_number),
        "Flags": eligible["original_flags"],
    })


def _st2_flag_legend_html() -> str:
    legend = pd.DataFrame([
        ("negligible potency", "Baseline potency was too low for the benchmark row to be meaningful."),
        ("infection-resistance benchmark not assigned", "No infection-resistance calibration benchmark was assigned to the row."),
        ("missing infection-resistance simulation", "The simulation infection-resistance metric was not defined or calculable."),
        ("average resistant-level benchmark not assigned", "No expert-assigned model benchmark was assigned to the average resistant-level metric."),
        ("no resistant infections for average metric", "Average resistant level is not meaningful because no resistant infection-days were observed."),
        ("low resistant sample", "Resistant infection-days were below 50 or the original summary flagged a low sample."),
        ("low infection-days", "Infection-days were below 100."),
        ("low carrier-days", "Carrier-days were below 1,000."),
        ("missing microbiome denominator", "Carrier-days were missing, so microbiome/carriage resistance is not interpretable."),
        ("missing carrier-days", "Carrier-day denominator was missing."),
        ("expanded window", "The original calibration summary flagged an expanded observation window."),
        ("benchmark value varied across runs", "Benchmarks differed across supplied calibration summaries; the first non-missing benchmark is displayed."),
        ("no caveat", "No derived caveat was added."),
    ], columns=["Derived flag", "Meaning"])
    return "<h2>Flag Legend</h2>\n" + _html_table(legend)


def _st2_notes_html(multiple_runs: bool, missing_columns: list[str], target_inconsistency_count: int) -> str:
    notes = [
        "Data source: Resistance Benchmarks table parsed from calibration_summary_*.txt. "
        "calibration_summary_*.txt is an internal diagnostic file and is not required to interpret this page.",
        "Values refer to the shared calibration window reported in the source calibration summaries.",
        "Infection resistance simulation (%) is the simulated percentage of infection-days for the bacterium-drug combination classified as resistant.",
        "Infection resistance calibration benchmark (%) is an evidence-informed comparison value, not a direct harmonised surveillance estimate.",
        "Average resistant-level comparison values are expert-assigned model benchmarks for mean any_r conditional on any_r > 0; they are not MIC values or direct surveillance estimates.",
        "Average resistant level is summarised among resistant positives where defined; rows with no resistant infections do not have a meaningful average resistant level.",
        "Microbiome resistance simulation (%) describes simulated resistance in the microbiome/carriage reservoir and is not clinical isolate resistance.",
        "Infection-days, resistant infection-days, and carrier-days are the denominators used for infection resistance, average resistant-level, and microbiome-resistance summaries respectively.",
        "Flags are copied from the calibration summary and supplemented with simple display reliability flags.",
        "Rows flagged as negligible potency correspond to bacterium-drug combinations where baseline potency was too low for the benchmark to be meaningful.",
        "Missing values indicate that a benchmark was not assigned or that the simulated metric was not defined or calculable for that bacterium-drug combination.",
        "This table is detailed calibration support. It should not be interpreted as an independent surveillance dataset.",
        "Detailed bibliographic source attribution for each benchmark value is not encoded in the generated calibration outputs and is therefore not shown here.",
        "Rows excluded from the largest-differences summary because simulation or benchmark values are missing are still retained in the full table.",
    ]
    if multiple_runs:
        notes.append(
            "For multiple calibration summaries, simulation values are calculated within run first and displayed as medians; "
            "5th-95th percentile intervals are shown for simulated percentage columns. Denominators are reported as medians."
        )
    else:
        notes.append("Single calibration summary supplied: no cross-run interval is shown.")
    if missing_columns:
        notes.append("Missing expected source columns in at least one parsed table: " + ", ".join(missing_columns) + ".")
    if target_inconsistency_count:
        notes.append(
            f"{target_inconsistency_count} row{'s' if target_inconsistency_count != 1 else ''} had benchmark values that varied across runs; "
            "the first non-missing benchmark is displayed and the row is flagged."
        )
    return _html_footnotes(notes)


def _st2_meta_box(agg: dict | None, stats: dict[str, object]) -> str:
    meta = agg.get("meta", {}) if agg is not None else {}
    target_year = meta.get("target_year", "\u2014")
    window = meta.get("window_duration", "2022\u20132025 calibration window")
    parts = [
        f"<strong>Target year:</strong> {target_year}",
        f"<strong>Window:</strong> {window}",
        f"<strong>Calibration summaries parsed:</strong> {stats['n_runs']}",
        f"<strong>Bacterium-drug rows:</strong> {stats['n_rows']:,}",
        f"<strong>Bacteria:</strong> {stats['n_bacteria']:,}",
        f"<strong>Drugs:</strong> {stats['n_drugs']:,}",
        f"<strong>Rows with infection simulation and benchmark:</strong> {stats['n_complete_inf']:,}",
        f"<strong>Negligible potency rows:</strong> {stats['n_negligible']:,}",
        f"<strong>Low resistant-sample rows:</strong> {stats['n_low_resistant_sample']:,}",
    ]
    return "<div class='meta-box'>" + " &nbsp;|&nbsp; ".join(parts) + "</div>\n"


def _st2_placeholder(
    out_dir: Path,
    agg: dict | None,
    message: str,
    problems: list[str] | None = None,
) -> None:
    body = _html_head(_ST2_TITLE)
    body += _back_link()
    body += f"<h1>{_ST2_TITLE}</h1>\n"
    if agg is not None:
        body += _meta_box(agg)
    body += f"<p class='note'>{message}</p>\n"
    body += "<h2>Required source columns</h2>\n<ul>\n"
    for column in _ST2_REQUIRED_COLUMNS:
        body += f"<li><code>{column}</code></li>\n"
    body += "</ul>\n"
    if problems:
        body += "<h2>Parser Notes</h2>\n<ul>\n"
        for problem in problems[:8]:
            body += f"<li>{problem}</li>\n"
        body += "</ul>\n"
    body += _st2_notes_html(False, [], 0)
    body += "</body></html>"
    _save(out_dir / TABLES_DIRNAME / f"{_ST2_STEM}.html", body)


def make_supplementary_table_s2_resistance_benchmarks(
    runs: list[dict],
    out_dir: Path,
    agg: dict | None = None,
) -> None:
    raw, problems, missing_columns = _st2_rows_from_runs(runs)
    summary = _st2_summarise_rows(raw)
    if summary.empty:
        _st2_placeholder(out_dir, agg, _ST2_PLACEHOLDER_MESSAGE, problems)
        if problems:
            print("  Supplementary Table S2: placeholder; " + " ".join(problems[:3]))
        else:
            print("  Supplementary Table S2: placeholder (Resistance Benchmarks table not found).")
        return

    n_runs = len(runs)
    multiple_runs = n_runs > 1
    full_table = _st2_display_table(summary, multiple_runs)
    largest_table = _st2_largest_differences_table(summary)
    target_inconsistency_count = int(summary["target_inconsistent"].sum())
    stats = {
        "n_runs": n_runs,
        "n_rows": len(summary),
        "n_bacteria": int(summary["bacterium"].nunique()),
        "n_drugs": int(summary["drug"].nunique()),
        "n_complete_inf": int(
            (summary["inf_sim_median"].notna() & summary["inf_target"].notna()).sum()
        ),
        "n_negligible": int(summary["derived_flags"].str.contains("negligible potency", na=False).sum()),
        "n_low_resistant_sample": int(summary["derived_flags"].str.contains("low resistant sample", na=False).sum()),
    }

    body = _html_head(_ST2_TITLE)
    body += _back_link()
    body += f"<h1>{_ST2_TITLE}</h1>\n"
    body += (
        "<p class='subtitle'>Detailed bacterium-drug resistance benchmark rows parsed "
        "from accepted calibration summaries.</p>\n"
    )
    body += _st2_meta_box(agg, stats)
    if missing_columns:
        body += (
            "<p class='note'>Some expected source columns were missing in at least one parsed "
            "Resistance Benchmarks table: "
            + ", ".join(missing_columns)
            + ". Available columns are still shown where possible.</p>\n"
        )

    body += "<h2>Largest Absolute Infection-Resistance Differences</h2>\n"
    body += (
        "<p class='note'>Rows in this compact summary require both non-missing simulation and "
        "target infection-resistance values. All rows remain included in the full table below.</p>\n"
    )
    body += _html_table(largest_table)

    body += "<h2>Full Detailed Resistance Benchmark Table</h2>\n"
    body += (
        "<div style='overflow-x:auto; max-height:78vh; border:1px solid #ddd; "
        "padding:0 4px; margin-bottom:14px;'>\n"
    )
    body += _html_table(full_table)
    body += "\n</div>\n"
    body += _st2_flag_legend_html()
    body += _st2_notes_html(multiple_runs, missing_columns, target_inconsistency_count)
    if problems:
        body += "<h2>Parser Notes</h2>\n<ul>\n"
        for problem in problems[:8]:
            body += f"<li>{problem}</li>\n"
        body += "</ul>\n"
    body += "</body></html>"
    _save(out_dir / TABLES_DIRNAME / f"{_ST2_STEM}.html", body)
    print(
        "  Supplementary Table S2: "
        f"{len(summary)} bacterium-drug row{'s' if len(summary) != 1 else ''}, "
        f"{stats['n_bacteria']} bacteria, {stats['n_drugs']} drugs, "
        f"{stats['n_complete_inf']} row{'s' if stats['n_complete_inf'] != 1 else ''} with infection simulation and target, "
        f"from {n_runs} calibration summar{'ies' if n_runs != 1 else 'y'}."
    )


_SF1_TITLE = "Supplementary Figure S1. Antibiotic activity retained after resistance among treated active infections"
_SF1_STEM = "Supplementary_Figure_S1__potential_activity_retained"
_SF1_NUMERATOR_COLUMN = "potential_activity_existing_drugs_sum_by_bacteria"
_SF1_DENOMINATOR_COLUMN = "max_possible_potential_activity_existing_drugs_sum_by_bacteria"
_SF1_REQUIRED_MESSAGE = (
    "Supplementary Figure S1 requires simulation_summary columns "
    "`potential_activity_existing_drugs_sum_by_bacteria` and "
    "`max_possible_potential_activity_existing_drugs_sum_by_bacteria`. Re-run the Rust "
    "simulation after adding these aggregate outputs."
)
_SF1_NEGLIGIBLE_POTENCY_THRESHOLD = 0.15
_SF1_SENTINEL_SLUGS = [
    "escherichia_coli",
    "klebsiella_pneumoniae",
    "staphylococcus_aureus",
    "pseudomonas_aeruginosa",
    "acinetobacter_baumannii",
    "neisseria_gonorrhoeae",
    "mycoplasma_genitalium",
    "streptococcus_pneumoniae",
]


def _sf1_placeholder(out_dir: Path, agg: dict | None, message: str) -> None:
    fig, ax = plt.subplots(figsize=(10, 3.8))
    ax.text(
        0.5,
        0.5,
        f"{_SF1_TITLE}\n\n{message}",
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=10.5,
        color="#555",
        bbox=dict(boxstyle="round,pad=0.6", fc="#f5f5f5", ec="#bbb"),
    )
    ax.set_axis_off()
    fig.subplots_adjust(left=0.03, right=0.97, top=0.92, bottom=0.08)
    _save_figure(fig, out_dir, _SF1_STEM, _SF1_TITLE, message, [], agg=agg)


def _sf1_rows_from_simulation_csv(csv_path: Path) -> tuple[list[dict[str, object]], str | None]:
    header = _simulation_csv_column_names(csv_path)
    if header is None:
        return [], f"{csv_path.name}: unable to read simulation summary CSV header."

    required = {_SF1_NUMERATOR_COLUMN, _SF1_DENOMINATOR_COLUMN}
    available = set(header)
    missing = sorted(required - available)
    if missing:
        return [], f"{csv_path.name}: missing columns {', '.join(missing)}."

    optional = ["policy_option", "run_id", "simulation_year", "year", "time_in_years", "time_step"]
    usecols = [column for column in list(required) + optional if column in available]
    try:
        df = _read_csv_selected(csv_path, usecols)
    except (ValueError, OSError) as exc:
        return [], f"{csv_path.name}: unable to load Supplementary Figure S1 columns ({exc})."

    if "policy_option" in df.columns:
        policy = pd.to_numeric(df["policy_option"], errors="coerce")
        df = df[policy == 0].copy()
    df["sf1_year"] = _simulation_year_series(df)
    df = df.dropna(subset=["sf1_year"]).copy()
    if df.empty:
        return [], f"{csv_path.name}: no baseline-policy rows with usable years."

    grouped = df.groupby("run_id", dropna=False) if "run_id" in df.columns else [(csv_path.stem, df)]
    n_known_bacteria = len(_F15_KNOWN_BACTERIA_SLUGS)
    rows: list[dict[str, object]] = []

    for run_key, run_df in grouped:
        totals: dict[tuple[int, int], list[float]] = {}
        for _, row in run_df.iterrows():
            year = int(math.floor(float(row["sf1_year"])))
            numerator_values = np.array(
                _figure_15_parse_vector_cell(row[_SF1_NUMERATOR_COLUMN]),
                dtype=float,
            )
            denominator_values = np.array(
                _figure_15_parse_vector_cell(row[_SF1_DENOMINATOR_COLUMN]),
                dtype=float,
            )
            target_len = max(n_known_bacteria, len(numerator_values), len(denominator_values))
            numerator_values = _figure_15_extend_array(numerator_values, target_len)
            denominator_values = _figure_15_extend_array(denominator_values, target_len)
            for b_idx in range(target_len):
                numerator = float(numerator_values[b_idx])
                denominator = float(denominator_values[b_idx])
                if not np.isfinite(numerator):
                    numerator = 0.0
                if not np.isfinite(denominator):
                    denominator = 0.0
                if denominator <= 0.0 and numerator <= 0.0:
                    continue
                bucket = totals.setdefault((year, b_idx), [0.0, 0.0])
                bucket[0] += numerator
                bucket[1] += denominator

        for (year, b_idx), (numerator, denominator) in totals.items():
            slug = (
                _F15_KNOWN_BACTERIA_SLUGS[b_idx]
                if b_idx < n_known_bacteria
                else f"bacterium_{b_idx + 1}"
            )
            rows.append({
                "source": csv_path.name,
                "run": str(run_key),
                "year": year,
                "bacterium_slug": slug,
                "bacterium": _figure_15_bacterium_label(slug),
                "numerator": numerator,
                "denominator": denominator,
            })

    if not rows:
        return [], f"{csv_path.name}: Supplementary Figure S1 columns contained no usable values."
    return rows, None


def _sf1_percent(numerator: object, denominator: object) -> float:
    denominator_f = float(denominator)
    if not np.isfinite(denominator_f) or denominator_f <= 0.0:
        return np.nan
    return 100.0 * float(numerator) / denominator_f


def make_supplementary_figure_s1_potential_activity_retained(
    csv_paths: list[Path],
    out_dir: Path,
    agg: dict | None = None,
) -> None:
    if not csv_paths:
        _sf1_placeholder(out_dir, agg, _SF1_REQUIRED_MESSAGE)
        print("  Supplementary Figure S1: placeholder (no simulation summary CSVs found).")
        return

    rows: list[dict[str, object]] = []
    problems: list[str] = []
    for csv_path in csv_paths:
        run_rows, problem = _sf1_rows_from_simulation_csv(csv_path)
        rows.extend(run_rows)
        if problem:
            problems.append(problem)

    if not rows:
        _sf1_placeholder(out_dir, agg, _SF1_REQUIRED_MESSAGE)
        if problems:
            print("  Supplementary Figure S1: placeholder; " + " ".join(problems[:3]))
        else:
            print("  Supplementary Figure S1: placeholder (no usable rows).")
        return

    df = pd.DataFrame(rows)
    per_run_overall = (
        df.groupby(["source", "run", "year"], as_index=False)[["numerator", "denominator"]]
        .sum()
    )
    per_run_overall["retained_percent"] = per_run_overall.apply(
        lambda row: _sf1_percent(row["numerator"], row["denominator"]),
        axis=1,
    )
    per_run_overall = per_run_overall.dropna(subset=["retained_percent"])
    if per_run_overall.empty:
        _sf1_placeholder(out_dir, agg, _SF1_REQUIRED_MESSAGE)
        print("  Supplementary Figure S1: placeholder (all denominator sums were zero).")
        return

    overall_summary = (
        per_run_overall.groupby("year", as_index=False)["retained_percent"]
        .agg(
            median="median",
            p5=lambda s: float(np.nanpercentile(s, 5)),
            p95=lambda s: float(np.nanpercentile(s, 95)),
            n="count",
        )
        .sort_values("year")
    )

    per_run_bacteria = df[df["denominator"] > 0.0].copy()
    per_run_bacteria["retained_percent"] = 100.0 * per_run_bacteria["numerator"] / per_run_bacteria["denominator"]
    sentinel_summary = (
        per_run_bacteria[per_run_bacteria["bacterium_slug"].isin(_SF1_SENTINEL_SLUGS)]
        .groupby(["bacterium_slug", "bacterium", "year"], as_index=False)["retained_percent"]
        .median()
    )

    n_runs = int(per_run_overall[["source", "run"]].drop_duplicates().shape[0])
    fig, axes = plt.subplots(1, 2, figsize=(13.0, 4.8))

    years = overall_summary["year"].to_numpy(dtype=float)
    axes[0].plot(years, overall_summary["median"].to_numpy(float), color="#2A9D8F", linewidth=2.0)
    if n_runs > 1:
        axes[0].fill_between(
            years,
            overall_summary["p5"].to_numpy(float),
            overall_summary["p95"].to_numpy(float),
            color="#2A9D8F",
            alpha=0.18,
        )
    axes[0].set_title("A. Overall active-infection weighted", loc="left", fontsize=10, fontweight="bold")
    axes[0].set_xlabel("Year", fontsize=9.5)
    axes[0].set_ylabel("Potential activity retained (%)", fontsize=9.5)
    axes[0].set_ylim(0, 100)
    axes[0].spines[["top", "right"]].set_visible(False)
    axes[0].grid(axis="y", linewidth=0.35, alpha=0.45)

    if sentinel_summary.empty:
        axes[1].text(
            0.5,
            0.5,
            "No sentinel bacterium rows had nonzero denominators.",
            ha="center",
            va="center",
            transform=axes[1].transAxes,
            color="#555",
        )
    else:
        for slug in _SF1_SENTINEL_SLUGS:
            one = sentinel_summary[sentinel_summary["bacterium_slug"] == slug].sort_values("year")
            if one.empty:
                continue
            axes[1].plot(
                one["year"].to_numpy(dtype=float),
                one["retained_percent"].to_numpy(dtype=float),
                linewidth=1.5,
                label=_figure_15_bacterium_label(slug),
            )
    axes[1].set_title("B. Sentinel bacteria", loc="left", fontsize=10, fontweight="bold")
    axes[1].set_xlabel("Year", fontsize=9.5)
    axes[1].set_ylabel("Potential activity retained (%)", fontsize=9.5)
    axes[1].set_ylim(0, 100)
    axes[1].spines[["top", "right"]].set_visible(False)
    axes[1].grid(axis="y", linewidth=0.35, alpha=0.45)
    axes[1].legend(fontsize=7.2, frameon=False, loc="center left", bbox_to_anchor=(1.01, 0.5))
    fig.suptitle(_SF1_TITLE, fontsize=11, fontweight="bold")
    fig.tight_layout(rect=[0, 0.02, 0.86, 0.94])

    available_years = sorted(int(year) for year in per_run_bacteria["year"].dropna().unique())
    table_year = 2025 if 2025 in available_years else max(available_years)
    table_base = per_run_bacteria[per_run_bacteria["year"] == table_year].copy()
    run_totals = (
        table_base.groupby(["source", "run"], as_index=False)["denominator"]
        .sum()
        .rename(columns={"denominator": "run_total_denominator"})
    )
    table_base = table_base.merge(run_totals, on=["source", "run"], how="left")
    table_base["denominator_share_percent"] = np.where(
        table_base["run_total_denominator"] > 0.0,
        100.0 * table_base["denominator"] / table_base["run_total_denominator"],
        np.nan,
    )
    table_summary = (
        table_base.groupby(["bacterium_slug", "bacterium"], as_index=False)
        .agg(
            retained_percent=("retained_percent", "median"),
            denominator_share_percent=("denominator_share_percent", "median"),
            denominator=("denominator", "median"),
        )
        .sort_values("retained_percent", ascending=True)
    )

    table = pd.DataFrame({
        "Bacterium": table_summary["bacterium"],
        f"Potential activity retained in {table_year} (%)": table_summary["retained_percent"].map(
            lambda v: f"{float(v):.1f}" if pd.notna(v) and np.isfinite(float(v)) else "-",
        ),
        f"Denominator share in {table_year} (%)": table_summary["denominator_share_percent"].map(
            lambda v: f"{float(v):.2f}" if pd.notna(v) and np.isfinite(float(v)) else "-",
        ),
        "No-resistance activity denominator": table_summary["denominator"].map(_figure_19_compact_count),
        "Denominator flag": np.where(
            table_summary["denominator_share_percent"] < 0.1,
            "low denominator",
            "",
        ),
    })
    extra_html = (
        f"<h2>Supplementary Figure S1 {table_year} Bacterium Table</h2>\n"
        + _html_table(table)
    )
    if table_year != 2025:
        extra_html += (
            f"<p class='note'>Calendar year 2025 was not available; table uses latest available "
            f"year {table_year}.</p>\n"
        )

    interval_note = (
        "Panel A shaded band shows 5th-95th percentile across runs. "
        if n_runs > 1
        else "Single run; no uncertainty interval is shown. "
    )
    footnotes = [
        "Supplementary Figure S1 estimates the potential treatment-option landscape, not realised "
        "prescribing. For each active infection, the numerator sums baseline potency x "
        "resistance-adjusted retained activity across drugs that existed at that time and had "
        "non-negligible baseline potency against the bacterium. The denominator sums the same "
        "baseline potencies without resistance. The metric therefore estimates the fraction of "
        "potential activity retained after resistance across available modelled drugs.",
        "This differs from the main activity-retained figure, which is weighted by antibiotics "
        "actually being used. Supplementary Figure S1 does not measure treatment adequacy and "
        "does not imply that all included drugs were clinically available, locally accessible, "
        "or prescribed.",
        "Drugs with baseline potency below the negligible-potency threshold are excluded from "
        f"the denominator. Negligible potency is defined as baseline potency < "
        f"{_SF1_NEGLIGIBLE_POTENCY_THRESHOLD:.2f}.",
        "Drug existence is based on the model's configured drug introduction day; regional access "
        "or person-level prescribing constraints are not applied to this first treatment-option "
        "landscape metric.",
        "Annual percentages are calculated by summing numerator and denominator first within each "
        "run-year, then dividing. " + interval_note + f"Runs: {n_runs}.",
    ]
    if problems:
        footnotes.append("Some simulation summary CSVs were skipped: " + " ".join(problems[:3]))

    _save_figure(
        fig,
        out_dir,
        _SF1_STEM,
        _SF1_TITLE,
        "Baseline-policy active-infection weighted potential activity retained across modelled drugs "
        "that had been introduced and had non-negligible baseline potency.",
        footnotes,
        agg=agg,
        extra_html=extra_html,
    )
    print(
        "  Supplementary Figure S1: "
        f"{len(overall_summary)} annual point{'s' if len(overall_summary) != 1 else ''} "
        f"from {n_runs} run{'s' if n_runs != 1 else ''}."
    )


_SF2_TITLE = "Supplementary Figure S2. Microbiome resistance reservoir, 2022\u20132025"
_SF2_STEM = "Supplementary_Figure_S2__microbiome_resistance_reservoir"
_SF2_REQUIRED_MESSAGE = (
    "Supplementary Figure S2 requires microbiome-resistance fields from "
    "calibration_summary_*.txt, including Microbiome Resistance Prevalence (%) "
    "and Resistance Benchmarks columns such as Micro sim (%), Inf sim (%), "
    "Carrier days, Inf days, and Class."
)


def _sf2_numeric(value: object) -> float:
    parsed = _first_numeric_value(value)
    if parsed is None or not np.isfinite(parsed):
        return np.nan
    return float(parsed)


def _sf2_clean_bacterium_name(value: object) -> str:
    text = str(value or "").strip()
    text = re.sub(r"\s*\*$", "", text)
    slug = text.lower().replace(" ", "_")
    return _figure_15_bacterium_label(slug)


def _sf2_axis_placeholder(ax, title: str, message: str) -> None:
    ax.set_title(title, loc="left", fontsize=10, fontweight="bold")
    ax.text(
        0.5,
        0.5,
        message,
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=9.2,
        color="#555",
        wrap=True,
        bbox=dict(boxstyle="round,pad=0.45", fc="#f5f5f5", ec="#bbb"),
    )
    ax.set_axis_off()


def _sf2_placeholder(out_dir: Path, agg: dict | None, message: str) -> None:
    fig, ax = plt.subplots(figsize=(10, 3.8))
    _sf2_axis_placeholder(ax, _SF2_TITLE, message)
    fig.subplots_adjust(left=0.03, right=0.97, top=0.88, bottom=0.08)
    _save_figure(fig, out_dir, _SF2_STEM, _SF2_TITLE, message, [], agg=agg)


def _sf2_panel_a_rows_from_runs(
    runs: list[dict],
) -> tuple[list[dict[str, object]], list[str]]:
    rows: list[dict[str, object]] = []
    problems: list[str] = []
    required = ["Bacteria", "Microbiome Resistance Prevalence (%)"]
    for run in runs:
        run_id = run.get("meta", {}).get("run_id", "run")
        bi = run.get("bacteria_infections", pd.DataFrame())
        if bi is None or bi.empty:
            problems.append(f"{run_id}: missing Bacteria Burden Benchmarks table.")
            continue
        missing = [column for column in required if column not in bi.columns]
        if missing:
            problems.append(f"{run_id}: bacteria burden table missing {', '.join(missing)}.")
            continue
        for _, row in bi.iterrows():
            bacterium = _sf2_clean_bacterium_name(row.get("Bacteria"))
            micro_prev = _sf2_numeric(row.get("Microbiome Resistance Prevalence (%)"))
            carriage = _sf2_numeric(row.get("Carriage simulation (%)"))
            if not np.isfinite(micro_prev):
                continue
            rows.append({
                "run": str(run_id),
                "bacterium": bacterium,
                "microbiome_resistance_percent": micro_prev,
                "carriage_simulation_percent": carriage,
            })
    return rows, problems


def _sf2_class_rows_from_runs(
    runs: list[dict],
) -> tuple[list[dict[str, object]], list[str]]:
    rows: list[dict[str, object]] = []
    problems: list[str] = []
    required = ["Class", "Micro sim (%)", "Inf sim (%)", "Carrier days", "Inf days"]
    for run in runs:
        run_id = run.get("meta", {}).get("run_id", "run")
        rb = run.get("resistance_benchmarks", pd.DataFrame())
        if rb is None or rb.empty:
            problems.append(f"{run_id}: missing Resistance Benchmarks table.")
            continue
        missing = [column for column in required if column not in rb.columns]
        if missing:
            problems.append(f"{run_id}: resistance benchmarks missing {', '.join(missing)}.")
            continue
        if "Flags" not in rb.columns:
            problems.append(f"{run_id}: resistance benchmarks missing Flags; negligible-potency exclusion unavailable.")

        rb = rb.copy()
        rb["sf2_class"] = rb["Class"].astype(str).str.strip()
        rb = rb[rb["sf2_class"].ne("") & rb["sf2_class"].ne("nan")].copy()
        for drug_class, class_rows in rb.groupby("sf2_class", dropna=False):
            total_rows = int(len(class_rows))
            carrier_num = 0.0
            carrier_den = 0.0
            infection_num = 0.0
            infection_den = 0.0
            included_rows = 0
            excluded_rows = 0
            for _, row in class_rows.iterrows():
                flags = str(row.get("Flags", "")).lower()
                negligible = "negligible" in flags
                micro = _sf2_numeric(row.get("Micro sim (%)"))
                inf = _sf2_numeric(row.get("Inf sim (%)"))
                carrier_days = _sf2_numeric(row.get("Carrier days"))
                inf_days = _sf2_numeric(row.get("Inf days"))
                include = (
                    not negligible
                    and np.isfinite(micro)
                    and np.isfinite(inf)
                    and np.isfinite(carrier_days)
                    and np.isfinite(inf_days)
                    and carrier_days > 0.0
                    and inf_days > 0.0
                )
                if include:
                    carrier_num += micro * carrier_days
                    carrier_den += carrier_days
                    infection_num += inf * inf_days
                    infection_den += inf_days
                    included_rows += 1
                else:
                    excluded_rows += 1

            if included_rows == 0 or carrier_den <= 0.0 or infection_den <= 0.0:
                continue
            micro_percent = carrier_num / carrier_den
            infection_percent = infection_num / infection_den
            rows.append({
                "run": str(run_id),
                "drug_class": str(drug_class),
                "microbiome_resistance_percent": micro_percent,
                "infection_resistance_percent": infection_percent,
                "difference_pp": micro_percent - infection_percent,
                "carrier_days": carrier_den,
                "infection_days": infection_den,
                "included_rows": included_rows,
                "excluded_rows": excluded_rows,
                "total_rows": total_rows,
            })
    return rows, problems


def _sf2_summarise_panel_a(rows: list[dict[str, object]]) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    return (
        df.groupby("bacterium", as_index=False)
        .agg(
            microbiome_resistance_percent=("microbiome_resistance_percent", "median"),
            microbiome_resistance_p5=("microbiome_resistance_percent", lambda s: float(np.nanpercentile(s, 5))),
            microbiome_resistance_p95=("microbiome_resistance_percent", lambda s: float(np.nanpercentile(s, 95))),
            carriage_simulation_percent=("carriage_simulation_percent", "median"),
            n=("microbiome_resistance_percent", "count"),
        )
        .sort_values("microbiome_resistance_percent", ascending=False)
        .reset_index(drop=True)
    )


def _sf2_summarise_classes(rows: list[dict[str, object]]) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    summary = (
        df.groupby("drug_class", as_index=False)
        .agg(
            microbiome_resistance_percent=("microbiome_resistance_percent", "median"),
            infection_resistance_percent=("infection_resistance_percent", "median"),
            difference_pp=("difference_pp", "median"),
            carrier_days=("carrier_days", "median"),
            infection_days=("infection_days", "median"),
            included_rows=("included_rows", "median"),
            excluded_rows=("excluded_rows", "median"),
            total_rows=("total_rows", "median"),
            n=("microbiome_resistance_percent", "count"),
        )
    )
    total_carrier = float(summary["carrier_days"].sum())
    total_infection = float(summary["infection_days"].sum())

    def reliability(row: pd.Series) -> str:
        flags: list[str] = []
        carrier_share = (
            float(row["carrier_days"]) / total_carrier
            if total_carrier > 0.0 and np.isfinite(float(row["carrier_days"]))
            else np.nan
        )
        infection_share = (
            float(row["infection_days"]) / total_infection
            if total_infection > 0.0 and np.isfinite(float(row["infection_days"]))
            else np.nan
        )
        total_rows_f = float(row["total_rows"])
        excluded_rows_f = float(row["excluded_rows"])
        if not np.isfinite(carrier_share) or carrier_share < 0.001:
            flags.append("low carrier denominator")
        if not np.isfinite(infection_share) or infection_share < 0.001:
            flags.append("low infection denominator")
        if total_rows_f > 0 and excluded_rows_f / total_rows_f >= 0.25 and excluded_rows_f >= 2:
            flags.append("many excluded rows")
        return "; ".join(flags)

    summary["reliability_flag"] = summary.apply(reliability, axis=1)
    return summary.sort_values("microbiome_resistance_percent", ascending=False).reset_index(drop=True)


def _sf2_format_percent(value: object) -> str:
    if value is None or pd.isna(value):
        return "\u2014"
    value_f = float(value)
    if not np.isfinite(value_f):
        return "\u2014"
    return f"{value_f:.1f}"


def _sf2_format_days(value: object) -> str:
    if value is None or pd.isna(value):
        return "\u2014"
    value_f = float(value)
    if not np.isfinite(value_f):
        return "\u2014"
    return f"{int(np.rint(value_f)):,}"


def _sf2_format_int(value: object) -> str:
    if value is None or pd.isna(value):
        return "\u2014"
    value_f = float(value)
    if not np.isfinite(value_f):
        return "\u2014"
    return f"{int(np.rint(value_f)):,}"


def make_supplementary_figure_s2_microbiome_resistance_reservoir(
    runs: list[dict],
    out_dir: Path,
    agg: dict | None = None,
) -> None:
    panel_a_rows, panel_a_problems = _sf2_panel_a_rows_from_runs(runs)
    class_rows, class_problems = _sf2_class_rows_from_runs(runs)
    panel_a = _sf2_summarise_panel_a(panel_a_rows)
    class_summary = _sf2_summarise_classes(class_rows)
    n_runs = len(runs)

    if panel_a.empty and class_summary.empty:
        _sf2_placeholder(out_dir, agg, _SF2_REQUIRED_MESSAGE)
        print("  Supplementary Figure S2: placeholder (required microbiome resistance tables missing).")
        return

    panel_a_count = int(len(panel_a))
    class_count = int(len(class_summary))
    fig_height = max(8.2, 2.2 + 0.25 * max(panel_a_count, class_count, 16))
    fig, axes = plt.subplots(
        1,
        3,
        figsize=(16.5, fig_height),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.35, 1.05, 1.05], "wspace": 0.55},
    )

    if panel_a.empty:
        _sf2_axis_placeholder(
            axes[0],
            "A. Microbiome resistance prevalence by bacterium",
            "Bacteria burden table or Microbiome Resistance Prevalence (%) column missing.",
        )
    else:
        plot_a = panel_a.sort_values("microbiome_resistance_percent", ascending=True)
        y = np.arange(len(plot_a))
        axes[0].barh(y, plot_a["microbiome_resistance_percent"].to_numpy(float), color="#3F7F93")
        axes[0].set_yticks(y)
        axes[0].set_yticklabels(plot_a["bacterium"].values, fontsize=7.0, fontstyle="italic")
        axes[0].set_xlabel("Microbiome resistance prevalence (%)", fontsize=9)
        axes[0].set_title(
            "A. Microbiome resistance prevalence by bacterium",
            loc="left",
            fontsize=10,
            fontweight="bold",
        )
        axes[0].set_xlim(0, max(100.0, float(plot_a["microbiome_resistance_percent"].max()) * 1.05))
        axes[0].spines[["top", "right"]].set_visible(False)
        axes[0].grid(axis="x", linewidth=0.35, alpha=0.45)

    if class_summary.empty:
        _sf2_axis_placeholder(
            axes[1],
            "B. Microbiome resistance by drug class",
            "Resistance Benchmarks columns for Micro sim (%), Carrier days, Inf sim (%), Inf days, and Class were not all available.",
        )
        _sf2_axis_placeholder(
            axes[2],
            "C. Microbiome versus active-infection resistance",
            "Drug-class weighted resistance summaries could not be calculated.",
        )
    else:
        plot_b = class_summary.sort_values("microbiome_resistance_percent", ascending=True)
        y = np.arange(len(plot_b))
        axes[1].barh(y, plot_b["microbiome_resistance_percent"].to_numpy(float), color="#2A9D8F")
        axes[1].set_yticks(y)
        axes[1].set_yticklabels(plot_b["drug_class"].values, fontsize=7.0)
        axes[1].set_xlabel("Carrier-day-weighted microbiome resistance (%)", fontsize=9)
        axes[1].set_title(
            "B. Microbiome resistance by drug class",
            loc="left",
            fontsize=10,
            fontweight="bold",
        )
        axes[1].set_xlim(0, max(100.0, float(plot_b["microbiome_resistance_percent"].max()) * 1.05))
        axes[1].spines[["top", "right"]].set_visible(False)
        axes[1].grid(axis="x", linewidth=0.35, alpha=0.45)

        x = class_summary["microbiome_resistance_percent"].to_numpy(float)
        y_vals = class_summary["infection_resistance_percent"].to_numpy(float)
        sizes = np.clip(
            24.0 + 12.0 * np.log10(class_summary["carrier_days"].to_numpy(float) + 1.0),
            28.0,
            120.0,
        )
        axes[2].scatter(x, y_vals, s=sizes, color="#D1495B", alpha=0.78, edgecolor="white", linewidth=0.6)
        lim = max(100.0, float(np.nanmax([np.nanmax(x), np.nanmax(y_vals)])) * 1.05)
        axes[2].plot([0, lim], [0, lim], color="#555", linestyle="--", linewidth=0.9)
        axes[2].set_xlim(0, lim)
        axes[2].set_ylim(0, lim)
        axes[2].set_xlabel("Microbiome resistance, carrier-day weighted (%)", fontsize=9)
        axes[2].set_ylabel("Active-infection resistance, infection-day weighted (%)", fontsize=9)
        axes[2].set_title(
            "C. Microbiome versus active-infection resistance",
            loc="left",
            fontsize=10,
            fontweight="bold",
        )
        axes[2].spines[["top", "right"]].set_visible(False)
        axes[2].grid(linewidth=0.35, alpha=0.45)
        outliers = class_summary.reindex(
            class_summary["difference_pp"].abs().sort_values(ascending=False).index
        ).head(5)
        for _, row in outliers.iterrows():
            axes[2].annotate(
                _F2_CLASS_SHORT.get(str(row["drug_class"]), str(row["drug_class"])),
                (float(row["microbiome_resistance_percent"]), float(row["infection_resistance_percent"])),
                xytext=(4, 4),
                textcoords="offset points",
                fontsize=7.0,
            )

    fig.suptitle(_SF2_TITLE, fontsize=11, fontweight="bold")

    extra_html = ""
    if not class_summary.empty:
        class_table_base = class_summary.copy()
        class_table_base["abs_difference"] = class_table_base["difference_pp"].abs()
        class_table_base = class_table_base.sort_values("abs_difference", ascending=False)
        class_table = pd.DataFrame({
            "Drug class": class_table_base["drug_class"],
            "Microbiome resistance (%)": class_table_base["microbiome_resistance_percent"].map(_sf2_format_percent),
            "Infection resistance (%)": class_table_base["infection_resistance_percent"].map(_sf2_format_percent),
            "Difference: microbiome minus infection (pp)": class_table_base["difference_pp"].map(_sf2_format_percent),
            "Carrier days": class_table_base["carrier_days"].map(_sf2_format_days),
            "Infection days": class_table_base["infection_days"].map(_sf2_format_days),
            "Number of bacterium-drug rows included": class_table_base["included_rows"].map(_sf2_format_int),
            "Reliability flag": class_table_base["reliability_flag"],
        })
        extra_html += "<h2>Drug-Class Microbiome vs Infection Resistance</h2>\n"
        extra_html += _html_table(class_table)
    if not panel_a.empty:
        bacterium_table = pd.DataFrame({
            "Bacterium": panel_a["bacterium"],
            "Microbiome resistance prevalence (%)": panel_a["microbiome_resistance_percent"].map(_sf2_format_percent),
            "Carriage simulation (%)": panel_a["carriage_simulation_percent"].map(_sf2_format_percent),
            "Reliability flag": np.where(
                pd.to_numeric(panel_a["carriage_simulation_percent"], errors="coerce") < 0.1,
                "low carriage prevalence",
                "",
            ),
        })
        extra_html += "<h2>Bacterium-Level Microbiome Reservoir Table</h2>\n"
        extra_html += _html_table(bacterium_table)

    footnotes = [
        "Microbiome resistance describes resistance in the simulated microbiome / carriage reservoir. "
        "It is not restricted to active infections and should not be interpreted as clinical isolate resistance.",
        "Panel B is carrier-day weighted. Panel C compares carrier-day-weighted microbiome resistance with "
        "infection-day-weighted active-infection resistance by drug class.",
        "Rows flagged as negligible potency or with missing denominators are excluded from the weighted drug-class summaries.",
        "Because this figure uses resistance benchmarks from the calibration summary, the values reflect the 2022-2025 "
        "calibration window unless otherwise stated.",
        "This figure is descriptive simulation output and does not introduce new calibration targets.",
        f"Values are medians across {n_runs} calibration summary run{'s' if n_runs != 1 else ''}; "
        "run-level weighted numerators and denominators are calculated before cross-run summarisation.",
    ]
    problems = panel_a_problems + class_problems
    if problems:
        footnotes.append("Parser notes: " + " ".join(problems[:4]))
    if n_runs > 1:
        footnotes.append("5th-95th intervals are omitted to keep the three-panel figure readable.")

    _save_figure(
        fig,
        out_dir,
        _SF2_STEM,
        _SF2_TITLE,
        "Microbiome resistance reservoir summaries from calibration_summary_*.txt.",
        footnotes,
        agg=agg,
        extra_html=extra_html,
    )
    status = "real data"
    if panel_a.empty or class_summary.empty:
        status = "partial data"
    print(
        "  Supplementary Figure S2: "
        f"{status}; {panel_a_count} bacterium row{'s' if panel_a_count != 1 else ''} "
        f"and {class_count} drug class{'es' if class_count != 1 else ''}."
    )


_SF3_TITLE = "Supplementary Figure S3. Carrier versus non-carrier infection incidence, 2022\u20132025"
_SF3_STEM = "Supplementary_Figure_S3__carrier_vs_non_carrier_infection_incidence"
_SF3_REQUIRED_MESSAGE = (
    "Supplementary Figure S3 requires simulation_summary columns for carrier/non-carrier "
    "at-risk person-days and new infections by bacterium."
)
_SF3_REQUIRED_COLUMNS = [
    "carrier_at_risk_person_days_by_bacteria",
    "non_carrier_at_risk_person_days_by_bacteria",
    "new_infections_in_carriers_by_bacteria",
    "new_infections_in_non_carriers_by_bacteria",
]
_SF3_ANY_R_COLUMNS = [
    "new_any_r_infections_in_carriers_by_bacteria",
    "new_any_r_infections_in_non_carriers_by_bacteria",
]


def _sf3_placeholder(out_dir: Path, agg: dict | None, message: str) -> None:
    fig, ax = plt.subplots(figsize=(10, 3.8))
    _sf2_axis_placeholder(ax, _SF3_TITLE, message)
    fig.subplots_adjust(left=0.03, right=0.97, top=0.88, bottom=0.08)
    _save_figure(fig, out_dir, _SF3_STEM, _SF3_TITLE, message, [], agg=agg)


def _sf3_vector_cell_to_array(value: object, target_len: int) -> np.ndarray:
    arr = np.zeros(target_len, dtype=float)
    values = _figure_15_parse_vector_cell(value)
    limit = min(len(values), target_len)
    if limit:
        arr[:limit] = np.asarray(values[:limit], dtype=float)
    return arr


def _sf3_sum_vector_column(df: pd.DataFrame, column: str, target_len: int) -> np.ndarray:
    total = np.zeros(target_len, dtype=float)
    if column not in df.columns:
        return total
    for value in df[column].values:
        total += _sf3_vector_cell_to_array(value, target_len)
    return total


def _sf3_rate_per_100k(events: float, person_days: float) -> float:
    if not np.isfinite(events) or not np.isfinite(person_days) or person_days <= 0.0:
        return np.nan
    return 100000.0 * float(events) / (float(person_days) / 365.0)


def _sf3_ratio(numerator_rate: float, denominator_rate: float) -> float:
    if not np.isfinite(numerator_rate) or not np.isfinite(denominator_rate):
        return np.nan
    if denominator_rate == 0.0:
        if numerator_rate > 0.0:
            return np.inf
        return np.nan
    return numerator_rate / denominator_rate


def _sf3_median_with_inf(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce").to_numpy(dtype=float)
    finite = numeric[np.isfinite(numeric)]
    if finite.size:
        return float(np.nanmedian(finite))
    if np.isposinf(numeric).any():
        return np.inf
    return np.nan


def _sf3_rows_from_csvs(csv_paths: list[Path]) -> tuple[list[dict[str, object]], list[str], bool]:
    rows: list[dict[str, object]] = []
    problems: list[str] = []
    any_r_seen = False
    target_len = len(_F15_KNOWN_BACTERIA_SLUGS)
    optional = ["policy_option", "run_id", "simulation_year", "year", "time_in_years", "time_step"]

    for csv_path in csv_paths:
        header = _simulation_csv_column_names(csv_path)
        if header is None:
            problems.append(f"{csv_path.name}: could not read simulation CSV header.")
            continue

        missing = [column for column in _SF3_REQUIRED_COLUMNS if column not in header]
        if missing:
            problems.append(f"{csv_path.name}: missing {', '.join(missing)}.")
            continue

        has_any_r = all(column in header for column in _SF3_ANY_R_COLUMNS)
        any_r_seen = any_r_seen or has_any_r
        wanted = set(_SF3_REQUIRED_COLUMNS + optional)
        if has_any_r:
            wanted.update(_SF3_ANY_R_COLUMNS)

        try:
            df = _read_csv_selected(csv_path, wanted)
        except (FileNotFoundError, ValueError, OSError) as exc:
            problems.append(f"{csv_path.name}: could not read required columns ({exc}).")
            continue
        if df.empty:
            problems.append(f"{csv_path.name}: simulation CSV has no rows.")
            continue

        if "policy_option" in df.columns:
            policy = pd.to_numeric(df["policy_option"], errors="coerce")
            df = df[policy.eq(0)].copy()
        years = _simulation_year_series(df)
        df = df[(years >= 2022.0) & (years < 2026.0)].copy()
        if df.empty:
            problems.append(f"{csv_path.name}: no baseline-policy rows in 2022-2025.")
            continue

        if "run_id" not in df.columns:
            df["run_id"] = csv_path.stem
        for run_id, group in df.groupby("run_id", dropna=False):
            carrier_days = _sf3_sum_vector_column(
                group, "carrier_at_risk_person_days_by_bacteria", target_len
            )
            non_carrier_days = _sf3_sum_vector_column(
                group, "non_carrier_at_risk_person_days_by_bacteria", target_len
            )
            carrier_events = _sf3_sum_vector_column(
                group, "new_infections_in_carriers_by_bacteria", target_len
            )
            non_carrier_events = _sf3_sum_vector_column(
                group, "new_infections_in_non_carriers_by_bacteria", target_len
            )
            if has_any_r:
                carrier_any_r_events = _sf3_sum_vector_column(
                    group, "new_any_r_infections_in_carriers_by_bacteria", target_len
                )
                non_carrier_any_r_events = _sf3_sum_vector_column(
                    group, "new_any_r_infections_in_non_carriers_by_bacteria", target_len
                )
            else:
                carrier_any_r_events = np.full(target_len, np.nan)
                non_carrier_any_r_events = np.full(target_len, np.nan)

            for b_idx, slug in enumerate(_F15_KNOWN_BACTERIA_SLUGS):
                carrier_rate = _sf3_rate_per_100k(carrier_events[b_idx], carrier_days[b_idx])
                non_carrier_rate = _sf3_rate_per_100k(
                    non_carrier_events[b_idx], non_carrier_days[b_idx]
                )
                any_r_carrier_rate = _sf3_rate_per_100k(
                    carrier_any_r_events[b_idx], carrier_days[b_idx]
                )
                any_r_non_carrier_rate = _sf3_rate_per_100k(
                    non_carrier_any_r_events[b_idx], non_carrier_days[b_idx]
                )
                rows.append({
                    "source": csv_path.name,
                    "run": str(run_id),
                    "bacterium_idx": b_idx,
                    "bacterium": _figure_15_bacterium_label(slug),
                    "carrier_py": carrier_days[b_idx] / 365.0,
                    "non_carrier_py": non_carrier_days[b_idx] / 365.0,
                    "carrier_events": carrier_events[b_idx],
                    "non_carrier_events": non_carrier_events[b_idx],
                    "carrier_rate": carrier_rate,
                    "non_carrier_rate": non_carrier_rate,
                    "rate_ratio": _sf3_ratio(carrier_rate, non_carrier_rate),
                    "carrier_any_r_events": carrier_any_r_events[b_idx],
                    "non_carrier_any_r_events": non_carrier_any_r_events[b_idx],
                    "carrier_any_r_rate": any_r_carrier_rate,
                    "non_carrier_any_r_rate": any_r_non_carrier_rate,
                    "any_r_rate_ratio": _sf3_ratio(any_r_carrier_rate, any_r_non_carrier_rate),
                    "any_r_available": has_any_r,
                })
    return rows, problems, any_r_seen


def _sf3_reliability_flag(row: pd.Series) -> str:
    flags: list[str] = []
    carrier_py = float(row.get("carrier_py", np.nan))
    non_carrier_py = float(row.get("non_carrier_py", np.nan))
    carrier_events = float(row.get("carrier_events", np.nan))
    non_carrier_events = float(row.get("non_carrier_events", np.nan))
    carrier_rate = float(row.get("carrier_rate", np.nan))
    non_carrier_rate = float(row.get("non_carrier_rate", np.nan))

    if not np.isfinite(carrier_py) or carrier_py <= 0.0 or not np.isfinite(non_carrier_py) or non_carrier_py <= 0.0:
        flags.append("zero denominator")
    if np.isfinite(carrier_py) and 0.0 < carrier_py < 100.0:
        flags.append("low carrier denominator")
    if np.isfinite(non_carrier_py) and 0.0 < non_carrier_py < 100.0:
        flags.append("low non-carrier denominator")
    if np.isfinite(carrier_events) and carrier_events < 20.0:
        flags.append("low carrier events")
    if np.isfinite(non_carrier_events) and non_carrier_events < 20.0:
        flags.append("low non-carrier events")
    if np.isfinite(non_carrier_rate) and non_carrier_rate == 0.0:
        if np.isfinite(carrier_rate) and carrier_rate > 0.0:
            flags.append("zero comparison rate")
        elif np.isfinite(carrier_rate) and carrier_rate == 0.0:
            flags.append("zero comparison rate")
    return "; ".join(dict.fromkeys(flags))


def _sf3_summarise(rows: list[dict[str, object]]) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    summary = (
        df.groupby(["bacterium_idx", "bacterium"], as_index=False)
        .agg(
            carrier_py=("carrier_py", "median"),
            non_carrier_py=("non_carrier_py", "median"),
            carrier_events=("carrier_events", "median"),
            non_carrier_events=("non_carrier_events", "median"),
            carrier_rate=("carrier_rate", "median"),
            non_carrier_rate=("non_carrier_rate", "median"),
            rate_ratio=("rate_ratio", _sf3_median_with_inf),
            carrier_any_r_events=("carrier_any_r_events", "median"),
            non_carrier_any_r_events=("non_carrier_any_r_events", "median"),
            carrier_any_r_rate=("carrier_any_r_rate", "median"),
            non_carrier_any_r_rate=("non_carrier_any_r_rate", "median"),
            any_r_rate_ratio=("any_r_rate_ratio", _sf3_median_with_inf),
            n_runs=("run", "nunique"),
        )
        .reset_index(drop=True)
    )
    summary["reliability_flag"] = summary.apply(_sf3_reliability_flag, axis=1)
    finite_ratio = pd.to_numeric(summary["rate_ratio"], errors="coerce").replace([np.inf, -np.inf], np.nan)
    max_finite = float(finite_ratio.max()) if finite_ratio.notna().any() else 1.0
    summary["ratio_sort"] = np.where(
        np.isposinf(summary["rate_ratio"]),
        max(max_finite * 1.5, 10.0),
        pd.to_numeric(summary["rate_ratio"], errors="coerce").fillna(-1.0),
    )
    summary = summary.sort_values(
        ["ratio_sort", "carrier_rate", "non_carrier_rate"],
        ascending=[False, False, False],
    ).reset_index(drop=True)
    return summary


def _sf3_format_number(value: object, digits: int = 1) -> str:
    if value is None or pd.isna(value):
        return "\u2014"
    value_f = float(value)
    if np.isposinf(value_f):
        return "\u221e"
    if not np.isfinite(value_f):
        return "\u2014"
    if abs(value_f) >= 1000.0:
        return f"{value_f:,.0f}"
    return f"{value_f:,.{digits}f}"


def _sf3_format_count(value: object) -> str:
    if value is None or pd.isna(value):
        return "\u2014"
    value_f = float(value)
    if not np.isfinite(value_f):
        return "\u2014"
    return f"{int(np.rint(value_f)):,}"


def _sf3_plot_rate_panel(ax, plot_df: pd.DataFrame) -> None:
    y = np.arange(len(plot_df))
    carrier = pd.to_numeric(plot_df["carrier_rate"], errors="coerce").to_numpy(dtype=float)
    non_carrier = pd.to_numeric(plot_df["non_carrier_rate"], errors="coerce").to_numpy(dtype=float)
    for idx, (x_carrier, x_non_carrier) in enumerate(zip(carrier, non_carrier)):
        if np.isfinite(x_carrier) and np.isfinite(x_non_carrier):
            ax.plot([x_non_carrier, x_carrier], [idx, idx], color="#B8B8B8", linewidth=0.8, zorder=1)
    ax.scatter(non_carrier, y, color="#4E79A7", s=26, label="Non-carrier", zorder=3)
    ax.scatter(carrier, y, color="#D1495B", s=26, label="Carrier", zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels(plot_df["bacterium"].values, fontsize=6.7, fontstyle="italic")
    ax.invert_yaxis()
    ax.set_xlabel("New active infection incidence per 100,000 at-risk person-years", fontsize=8.7)
    ax.set_title("A. Incidence among carriers and non-carriers", loc="left", fontsize=10, fontweight="bold")
    all_rates = np.concatenate([carrier[np.isfinite(carrier)], non_carrier[np.isfinite(non_carrier)]])
    positive = all_rates[all_rates > 0.0]
    if positive.size and float(np.nanmax(positive) / max(np.nanmin(positive), 1e-12)) > 100.0:
        ax.set_xscale("symlog", linthresh=1.0)
    ax.set_xlim(left=0.0)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", linewidth=0.35, alpha=0.45)
    ax.legend(loc="lower right", fontsize=8, frameon=False)


def _sf3_plot_ratio_panel(ax, plot_df: pd.DataFrame) -> None:
    y = np.arange(len(plot_df))
    ratios = pd.to_numeric(plot_df["rate_ratio"], errors="coerce").to_numpy(dtype=float)
    finite_positive = ratios[np.isfinite(ratios) & (ratios > 0.0)]
    cap = max(float(np.nanmax(finite_positive)) * 1.5, 10.0) if finite_positive.size else 10.0
    plot_x = np.full(len(plot_df), np.nan)
    marker = np.full(len(plot_df), "o", dtype=object)
    for idx, ratio in enumerate(ratios):
        if np.isposinf(ratio):
            plot_x[idx] = cap
            marker[idx] = ">"
        elif np.isfinite(ratio) and ratio > 0.0:
            plot_x[idx] = ratio
    ax.axvline(1.0, color="#555", linestyle="--", linewidth=0.9)
    for symbol, colour, label in [("o", "#6A4C93", "Finite ratio"), (">", "#D1495B", "Infinite, capped for display")]:
        mask = (marker == symbol) & np.isfinite(plot_x)
        if mask.any():
            ax.scatter(plot_x[mask], y[mask], marker=symbol, s=28, color=colour, label=label, zorder=3)
    ax.set_yticks(y)
    ax.set_yticklabels([])
    ax.invert_yaxis()
    ax.set_xscale("log")
    lower = min(float(np.nanmin(finite_positive)) / 1.5, 1.0 / 3.0) if finite_positive.size else 1.0 / 3.0
    upper = max(float(np.nanmax(plot_x[np.isfinite(plot_x)])) * 1.5, 3.0) if np.isfinite(plot_x).any() else 3.0
    ax.set_xlim(max(lower, 1e-3), upper)
    ax.set_xlabel("Carrier / non-carrier incidence rate ratio", fontsize=8.7)
    ax.set_title("B. Incidence rate ratio", loc="left", fontsize=10, fontweight="bold")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", linewidth=0.35, alpha=0.45)
    if np.isposinf(ratios).any():
        ax.legend(loc="lower right", fontsize=8, frameon=False)


def make_supplementary_figure_s3_carrier_vs_non_carrier_incidence(
    csv_paths: list[Path],
    out_dir: Path,
    agg: dict | None = None,
) -> None:
    rows, problems, any_r_available = _sf3_rows_from_csvs(csv_paths)
    summary = _sf3_summarise(rows)
    if summary.empty:
        message = _SF3_REQUIRED_MESSAGE
        if problems:
            message += " Parser notes: " + " ".join(problems[:3])
        _sf3_placeholder(out_dir, agg, message)
        print("  Supplementary Figure S3: placeholder (required simulation summary fields missing).")
        return

    fig_height = max(8.0, 2.0 + 0.24 * len(summary))
    fig, axes = plt.subplots(
        1,
        2,
        figsize=(15.5, fig_height),
        constrained_layout=True,
        gridspec_kw={"width_ratios": [1.25, 0.9], "wspace": 0.12},
    )
    plot_df = summary.copy()
    _sf3_plot_rate_panel(axes[0], plot_df)
    _sf3_plot_ratio_panel(axes[1], plot_df)
    fig.suptitle(_SF3_TITLE, fontsize=11, fontweight="bold")

    table_data = {
        "Bacterium": summary["bacterium"],
        "Carrier person-years at risk": summary["carrier_py"].map(lambda v: _sf3_format_number(v, 1)),
        "Non-carrier person-years at risk": summary["non_carrier_py"].map(lambda v: _sf3_format_number(v, 1)),
        "New infections among carriers": summary["carrier_events"].map(_sf3_format_count),
        "New infections among non-carriers": summary["non_carrier_events"].map(_sf3_format_count),
        "Carrier incidence per 100,000 person-years": summary["carrier_rate"].map(lambda v: _sf3_format_number(v, 1)),
        "Non-carrier incidence per 100,000 person-years": summary["non_carrier_rate"].map(lambda v: _sf3_format_number(v, 1)),
        "Carrier/non-carrier rate ratio": summary["rate_ratio"].map(lambda v: _sf3_format_number(v, 2)),
        "Reliability flag": summary["reliability_flag"],
    }
    if any_r_available:
        table_data.update({
            "New any-R infections among carriers": summary["carrier_any_r_events"].map(_sf3_format_count),
            "New any-R infections among non-carriers": summary["non_carrier_any_r_events"].map(_sf3_format_count),
            "Carrier any-R incidence per 100,000 person-years": summary["carrier_any_r_rate"].map(lambda v: _sf3_format_number(v, 1)),
            "Non-carrier any-R incidence per 100,000 person-years": summary["non_carrier_any_r_rate"].map(lambda v: _sf3_format_number(v, 1)),
            "Any-R carrier/non-carrier rate ratio": summary["any_r_rate_ratio"].map(lambda v: _sf3_format_number(v, 2)),
        })
    extra_html = "<h2>Bacterium-Level Carrier vs Non-Carrier Incidence</h2>\n"
    extra_html += _html_table(pd.DataFrame(table_data))

    n_runs = int(max(summary["n_runs"].max(), 1))
    footnotes = [
        "Carrier status is assessed before infection acquisition in each daily timestep. "
        "Carriers are people with microbiome/carriage presence for the bacterium and no active "
        "infection with that bacterium at the start of the timestep. Non-carriers have no "
        "microbiome/carriage presence and no active infection with that bacterium. Incident "
        "infections are new active infections occurring later in that timestep.",
        "Rates use at-risk person-years as the denominator. People already actively infected "
        "with the bacterium are excluded from both carrier and non-carrier denominators for that bacterium.",
        "This figure is a structural diagnostic of the carriage-to-infection pathway. It should "
        "not be interpreted as an external epidemiological estimate unless separately calibrated.",
        "Rates and rate ratios are calculated within each run from summed 2022-2025 baseline-policy "
        "numerators and denominators, then summarised as medians across runs.",
        f"Values are raw simulated counts and rates from {n_runs} simulation run{'s' if n_runs != 1 else ''}; "
        "5th-95th intervals are omitted to keep the paired-bacterium figure readable.",
    ]
    if any_r_available:
        footnotes.append(
            "Any-R incident infections are classified using the same infection-level any-R threshold "
            "used elsewhere in the simulation summaries."
        )
    else:
        footnotes.append("Any-R carrier/non-carrier rates were not available for this run.")
    if problems:
        footnotes.append("Parser notes: " + " ".join(problems[:4]))

    _save_figure(
        fig,
        out_dir,
        _SF3_STEM,
        _SF3_TITLE,
        "Baseline-policy carrier versus non-carrier infection incidence using simulation_summary_*.csv aggregate fields.",
        footnotes,
        agg=agg,
        extra_html=extra_html,
    )
    print(
        "  Supplementary Figure S3: real data; "
        f"{len(summary)} bacterium row{'s' if len(summary) != 1 else ''}; "
        f"any-R {'included' if any_r_available else 'not available'}."
    )


_SF4_TITLE = "Supplementary Figure SX. Modelled resistance mechanisms by bacterium, 2022\u20132025"
_SF4_STEM = "Supplementary_Figure_SX__modelled_resistance_mechanisms_by_bacterium"
_SF4_REQUIRED_MESSAGE = (
    "Supplementary Figure SX requires simulation_summary aggregate columns for active "
    "infection-days by bacterium and exact per-bacterium/per-ResistanceMechanism active "
    "infection-day counts named <bacterium>_infected_with_<mechanism>."
)

# Keep this mapping aligned with resistance_mechanism_family_idx() in src/simulation/simulation.rs.
_SF4_MECHANISM_FAMILIES: list[dict[str, object]] = [
    {
        "slug": "beta_lactamase_esbl_or_broad",
        "label": "ESBL / broad beta-lactam",
        "short_label": "ESBL /\nbroad BL",
        "variants": [
            "ResistanceMechanism::EnzymeEsblCtxM",
            "ResistanceMechanism::EnzymeEsblTem",
            "ResistanceMechanism::EnzymeEsblShv",
            "ResistanceMechanism::TargetSitePbp2aMecA",
            "ResistanceMechanism::EnzymeBlaZ",
            "ResistanceMechanism::EnzymeTem1",
            "ResistanceMechanism::MutationPbpMosaic",
        ],
        "notes": "Includes ESBL enzymes and broad beta-lactam target/penicillinase mechanisms.",
    },
    {
        "slug": "ampc",
        "label": "AmpC",
        "short_label": "AmpC",
        "variants": [
            "ResistanceMechanism::EnzymeAmpcCmy",
            "ResistanceMechanism::EnzymeAmpcDha",
            "ResistanceMechanism::MutationAmpCDerepression",
        ],
        "notes": "Plasmid AmpC enzymes and chromosomal AmpC derepression.",
    },
    {
        "slug": "carbapenemase",
        "label": "Carbapenemase",
        "short_label": "Carbapen-\nemase",
        "variants": [
            "ResistanceMechanism::EnzymeKpc",
            "ResistanceMechanism::EnzymeNdmVim",
            "ResistanceMechanism::EnzymeOxa48",
            "ResistanceMechanism::EnzymeOxaAcinetobacter",
        ],
        "notes": "KPC, metallo-beta-lactamase, OXA-48, and Acinetobacter OXA carbapenemases.",
    },
    {
        "slug": "porin_loss",
        "label": "Porin loss",
        "short_label": "Porin\nloss",
        "variants": [
            "ResistanceMechanism::PorinLossOmpk35_36",
            "ResistanceMechanism::PorinLossOprd",
            "ResistanceMechanism::GlobalPorinLoss",
        ],
        "notes": "Specific and global porin-loss mechanisms.",
    },
    {
        "slug": "efflux",
        "label": "Efflux",
        "short_label": "Efflux",
        "variants": [
            "ResistanceMechanism::EffluxAcrabTolc",
            "ResistanceMechanism::EffluxMexxyOprm",
            "ResistanceMechanism::GlobalEffluxPump",
            "ResistanceMechanism::EffluxMtrCde",
        ],
        "notes": "Specific and global efflux mechanisms.",
    },
    {
        "slug": "fluoroquinolone_target_or_qnr",
        "label": "Fluoroquinolone target / qnr",
        "short_label": "FQ target\n/ qnr",
        "variants": [
            "ResistanceMechanism::MutationGyrAPrimary",
            "ResistanceMechanism::MutationGyrAParCSecondary",
            "ResistanceMechanism::ProtectionQnr",
        ],
        "notes": "Chromosomal quinolone target mutations and qnr protection.",
    },
    {
        "slug": "macrolide_lincosamide_ribosomal",
        "label": "Macrolide / lincosamide / ribosomal",
        "short_label": "MLS /\nribosomal",
        "variants": [
            "ResistanceMechanism::TargetSiteErmB",
            "ResistanceMechanism::EnzymeMphA",
            "ResistanceMechanism::Mutation23sRrna",
        ],
        "notes": "Ribosomal methylation, macrolide phosphotransferase, and 23S rRNA macrolide mechanisms.",
    },
    {
        "slug": "aminoglycoside_ribosomal_or_enzyme",
        "label": "Aminoglycoside",
        "short_label": "Amino-\nglycoside",
        "variants": [
            "ResistanceMechanism::Enzyme16sRrmt",
            "ResistanceMechanism::EnzymeAacAph",
        ],
        "notes": "16S rRNA methyltransferase and aminoglycoside-modifying enzymes.",
    },
    {
        "slug": "phenicol_oxazolidinone",
        "label": "Phenicol / oxazolidinone",
        "short_label": "Phenicol /\noxazolid.",
        "variants": [
            "ResistanceMechanism::TargetSiteCfr",
            "ResistanceMechanism::EnzymeCat",
            "ResistanceMechanism::Mutation23sRrnaOxazolidinone",
        ],
        "notes": "Phenicol acetyltransferase, cfr, and oxazolidinone 23S rRNA mechanisms.",
    },
    {
        "slug": "tetracycline",
        "label": "Tetracycline",
        "short_label": "Tetra-\ncycline",
        "variants": [
            "ResistanceMechanism::ProtectionTetM",
            "ResistanceMechanism::EffluxTetAbc",
        ],
        "notes": "Ribosomal protection and tetracycline efflux.",
    },
    {
        "slug": "folate_pathway",
        "label": "Folate pathway",
        "short_label": "Folate\npathway",
        "variants": ["ResistanceMechanism::MutationFolatePathway"],
        "notes": "Sulfonamide/trimethoprim pathway marker.",
    },
    {
        "slug": "colistin",
        "label": "Colistin",
        "short_label": "Colistin",
        "variants": [
            "ResistanceMechanism::ModificationMcr1",
            "ResistanceMechanism::MutationPolymyxinRegulatory",
        ],
        "notes": "mcr-1 and polymyxin regulatory/lipid-A pathways.",
    },
    {
        "slug": "rifampicin",
        "label": "Rifampicin",
        "short_label": "Rifampicin",
        "variants": ["ResistanceMechanism::MutationRpoB"],
        "notes": "rpoB target mutation family.",
    },
    {
        "slug": "fosfomycin_nitrofuran",
        "label": "Fosfomycin / nitrofuran",
        "short_label": "Fosfo. /\nnitrofuran",
        "variants": [
            "ResistanceMechanism::MutationNitroreductase",
            "ResistanceMechanism::EnzymeFos",
        ],
        "notes": "Nitroreductase and fosfomycin enzyme mechanisms.",
    },
    {
        "slug": "daptomycin_fusidic",
        "label": "Daptomycin / fusidic",
        "short_label": "Dapto. /\nfusidic",
        "variants": [
            "ResistanceMechanism::MutationMprF",
            "ResistanceMechanism::MutationLiafsrCls",
            "ResistanceMechanism::ProtectionFusB",
        ],
        "notes": "Membrane/cell-envelope daptomycin mechanisms and fusidic-acid protection.",
    },
    {
        "slug": "other_unknown",
        "label": "Other / unknown",
        "short_label": "Other /\nunknown",
        "variants": [
            "ResistanceMechanism::TargetSiteVanA",
            "ResistanceMechanism::TargetSiteVanB",
            "ResistanceMechanism::AsYetUnknown",
        ],
        "notes": "Compact bucket for glycopeptide targets and calibration placeholder mechanisms.",
    },
]


_SF4_EXACT_MECHANISMS: list[dict[str, str]] = [
    {"variant": "ResistanceMechanism::EnzymeEsblCtxM", "slug": "enzyme_esbl_ctx_m", "label": "ESBL CTX-M"},
    {"variant": "ResistanceMechanism::EnzymeEsblTem", "slug": "enzyme_esbl_tem", "label": "ESBL TEM"},
    {"variant": "ResistanceMechanism::EnzymeEsblShv", "slug": "enzyme_esbl_shv", "label": "ESBL SHV"},
    {"variant": "ResistanceMechanism::EnzymeKpc", "slug": "enzyme_kpc", "label": "KPC"},
    {"variant": "ResistanceMechanism::EnzymeNdmVim", "slug": "enzyme_ndm_vim", "label": "NDM/VIM"},
    {"variant": "ResistanceMechanism::EnzymeOxa48", "slug": "enzyme_oxa_48", "label": "OXA-48"},
    {"variant": "ResistanceMechanism::EnzymeAmpcCmy", "slug": "enzyme_ampc_cmy", "label": "AmpC CMY"},
    {"variant": "ResistanceMechanism::EnzymeAmpcDha", "slug": "enzyme_ampc_dha", "label": "AmpC DHA"},
    {"variant": "ResistanceMechanism::MutationAmpCDerepression", "slug": "mutation_ampc_derepression", "label": "AmpC derepression"},
    {"variant": "ResistanceMechanism::TargetSitePbp2aMecA", "slug": "target_site_pbp2a_meca", "label": "PBP2a mecA"},
    {"variant": "ResistanceMechanism::TargetSiteVanA", "slug": "target_site_van_a", "label": "vanA"},
    {"variant": "ResistanceMechanism::TargetSiteVanB", "slug": "target_site_van_b", "label": "vanB"},
    {"variant": "ResistanceMechanism::MutationGyrAPrimary", "slug": "mutation_gyra_primary", "label": "gyrA primary"},
    {"variant": "ResistanceMechanism::MutationGyrAParCSecondary", "slug": "mutation_gyra_parc_secondary", "label": "gyrA/parC secondary"},
    {"variant": "ResistanceMechanism::ProtectionQnr", "slug": "protection_qnr", "label": "qnr protection"},
    {"variant": "ResistanceMechanism::Enzyme16sRrmt", "slug": "enzyme_16s_rrmt", "label": "16S rRNA methyltransferase"},
    {"variant": "ResistanceMechanism::TargetSiteErmB", "slug": "target_site_erm_b", "label": "ermB"},
    {"variant": "ResistanceMechanism::TargetSiteCfr", "slug": "target_site_cfr", "label": "cfr"},
    {"variant": "ResistanceMechanism::EnzymeCat", "slug": "enzyme_cat", "label": "CAT"},
    {"variant": "ResistanceMechanism::EffluxAcrabTolc", "slug": "efflux_acrab_tolc", "label": "AcrAB-TolC"},
    {"variant": "ResistanceMechanism::EffluxMexxyOprm", "slug": "efflux_mexxy_oprm", "label": "MexXY-OprM"},
    {"variant": "ResistanceMechanism::PorinLossOmpk35_36", "slug": "porin_loss_ompk35_36", "label": "OmpK35/36 loss"},
    {"variant": "ResistanceMechanism::PorinLossOprd", "slug": "porin_loss_oprd", "label": "OprD loss"},
    {"variant": "ResistanceMechanism::ModificationMcr1", "slug": "modification_mcr_1", "label": "mcr-1"},
    {"variant": "ResistanceMechanism::MutationPolymyxinRegulatory", "slug": "mutation_polymyxin_regulatory", "label": "polymyxin regulatory"},
    {"variant": "ResistanceMechanism::GlobalEffluxPump", "slug": "global_efflux_pump", "label": "global efflux pump"},
    {"variant": "ResistanceMechanism::GlobalPorinLoss", "slug": "global_porin_loss", "label": "global porin loss"},
    {"variant": "ResistanceMechanism::MutationFolatePathway", "slug": "mutation_folate_pathway", "label": "folate pathway"},
    {"variant": "ResistanceMechanism::MutationNitroreductase", "slug": "mutation_nitroreductase", "label": "nitroreductase"},
    {"variant": "ResistanceMechanism::EnzymeFos", "slug": "enzyme_fos", "label": "Fos enzyme"},
    {"variant": "ResistanceMechanism::MutationMprF", "slug": "mutation_mpr_f", "label": "mprF"},
    {"variant": "ResistanceMechanism::MutationLiafsrCls", "slug": "mutation_liafsr_cls", "label": "liaFSR/cls"},
    {"variant": "ResistanceMechanism::MutationRpoB", "slug": "mutation_rpo_b", "label": "rpoB"},
    {"variant": "ResistanceMechanism::ProtectionFusB", "slug": "protection_fus_b", "label": "fusB protection"},
    {"variant": "ResistanceMechanism::ProtectionTetM", "slug": "protection_tet_m", "label": "tetM protection"},
    {"variant": "ResistanceMechanism::EnzymeAacAph", "slug": "enzyme_aac_aph", "label": "AAC/APH"},
    {"variant": "ResistanceMechanism::EnzymeBlaZ", "slug": "enzyme_bla_z", "label": "blaZ"},
    {"variant": "ResistanceMechanism::EnzymeTem1", "slug": "enzyme_tem_1", "label": "TEM-1"},
    {"variant": "ResistanceMechanism::EnzymeMphA", "slug": "enzyme_mph_a", "label": "mphA"},
    {"variant": "ResistanceMechanism::EnzymeOxaAcinetobacter", "slug": "enzyme_oxa_acinetobacter", "label": "Acinetobacter OXA"},
    {"variant": "ResistanceMechanism::Mutation23sRrna", "slug": "mutation_23s_rrna", "label": "23S rRNA macrolide"},
    {"variant": "ResistanceMechanism::Mutation23sRrnaOxazolidinone", "slug": "mutation_23s_rrna_oxazolidinone", "label": "23S rRNA oxazolidinone"},
    {"variant": "ResistanceMechanism::EffluxTetAbc", "slug": "efflux_tet_abc", "label": "TetABC efflux"},
    {"variant": "ResistanceMechanism::MutationPbpMosaic", "slug": "mutation_pbp_mosaic", "label": "PBP mosaic"},
    {"variant": "ResistanceMechanism::EffluxMtrCde", "slug": "efflux_mtr_cde", "label": "MtrCDE efflux"},
    {"variant": "ResistanceMechanism::AsYetUnknown", "slug": "as_yet_unknown", "label": "as-yet-unknown"},
]


def _sf4_family_column(slug: str) -> str:
    return f"infection_days_with_mechanism_family_{slug}_by_bacteria"


def _sf4_family_lookup() -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    for family in _SF4_MECHANISM_FAMILIES:
        for variant in family["variants"]:
            lookup[str(variant)] = {
                "family": str(family["label"]),
                "notes": str(family["notes"]),
            }
    return lookup


def _sf4_mechanism_definitions_table_html() -> str:
    family_lookup = _sf4_family_lookup()
    rows: list[dict[str, str]] = []
    for mechanism in _SF4_EXACT_MECHANISMS:
        info = family_lookup.get(str(mechanism["variant"]), {})
        rows.append({
            "Mechanism enum variant": str(mechanism["variant"]),
            "Short display label": str(mechanism["label"]),
            "Broad interpretive family": info.get("family", "\u2014"),
            "Notes": info.get("notes", "\u2014"),
        })
    return "<h2>ResistanceMechanism variant definitions</h2>\n" + _html_table(pd.DataFrame(rows))


def _sf4_placeholder(out_dir: Path, agg: dict | None, message: str) -> dict[str, object]:
    fig, ax = plt.subplots(figsize=(10.5, 3.8))
    _sf2_axis_placeholder(ax, _SF4_TITLE, message)
    fig.subplots_adjust(left=0.03, right=0.97, top=0.88, bottom=0.08)
    footnotes = [
        "A real Supplementary Figure SX requires aggregate simulation_summary_*.csv fields with exact "
        "per-bacterium/per-ResistanceMechanism active infection-day counts.",
        "Grouped mechanism-family columns are not used as a fallback for the SX heatmap.",
    ]
    _save_figure(
        fig,
        out_dir,
        _SF4_STEM,
        _SF4_TITLE,
        message,
        footnotes,
        agg=agg,
        extra_html=_sf4_mechanism_definitions_table_html(),
    )
    return {"generated": "placeholder", "bacteria_included": 0, "mechanisms_included": 0}


def _sf4_sum_vector_column(df: pd.DataFrame, column: str, target_len: int) -> np.ndarray:
    total = np.zeros(target_len, dtype=float)
    if column not in df.columns:
        return total
    for value in df[column].values:
        total += _sf3_vector_cell_to_array(value, target_len)
    return total


def _sf4_exact_column(bacterium_slug: str, mechanism_slug: str) -> str:
    return f"{bacterium_slug}_infected_with_{mechanism_slug}"


def _sf4_exact_required_columns() -> list[str]:
    return [
        _sf4_exact_column(bacterium_slug, str(mechanism["slug"]))
        for bacterium_slug in _F15_KNOWN_BACTERIA_SLUGS
        for mechanism in _SF4_EXACT_MECHANISMS
    ]


def _sf4_percent(numerator: float, denominator: float) -> float:
    if not np.isfinite(numerator) or not np.isfinite(denominator) or denominator <= 0.0:
        return np.nan
    return 100.0 * numerator / denominator


def _sf4_rows_from_csvs(csv_paths: list[Path]) -> tuple[list[dict[str, object]], list[str], bool]:
    rows: list[dict[str, object]] = []
    problems: list[str] = []
    target_len = len(_F15_KNOWN_BACTERIA_SLUGS)
    exact_columns = _sf4_exact_required_columns()
    required = [
        "active_infection_days_by_bacteria",
        "infection_days_with_any_resistance_mechanism_by_bacteria",
        *exact_columns,
    ]
    optional = ["policy_option", "run_id", "simulation_year", "year", "time_in_years", "time_step"]
    any_new_active_available = False

    if not csv_paths:
        problems.append("No matching simulation_summary_*.csv files were found.")
        return rows, problems, any_new_active_available

    for csv_path in csv_paths:
        header = _simulation_csv_column_names(csv_path)
        if header is None:
            problems.append(f"{csv_path.name}: could not read simulation CSV header.")
            continue

        header_set = set(header)
        missing = [column for column in required if column not in header_set]
        if missing:
            problems.append(
                f"{csv_path.name}: missing exact SX fields; first missing columns: "
                + ", ".join(missing[:6])
            )
            continue

        new_active_available = "new_active_infections_by_bacteria" in header_set
        any_new_active_available = any_new_active_available or new_active_available
        wanted = [
            column
            for column in dict.fromkeys(
                required + optional + (["new_active_infections_by_bacteria"] if new_active_available else [])
            )
            if column in header_set
        ]
        n_mechanisms = len(_SF4_EXACT_MECHANISMS)
        run_states: dict[str, dict[str, object]] = {}
        matched_rows = 0
        try:
            chunk_iter = pd.read_csv(csv_path, usecols=wanted, chunksize=4000)
            for chunk in chunk_iter:
                if chunk.empty:
                    continue
                if "policy_option" in chunk.columns:
                    policy = pd.to_numeric(chunk["policy_option"], errors="coerce")
                    chunk = chunk[policy.eq(0)].copy()
                years = _simulation_year_series(chunk)
                chunk = chunk[(years >= 2022.0) & (years < 2026.0)].copy()
                if chunk.empty:
                    continue
                matched_rows += int(len(chunk))
                if "run_id" not in chunk.columns:
                    chunk["run_id"] = csv_path.stem

                for run_id, group in chunk.groupby("run_id", dropna=False):
                    run_key = str(run_id)
                    state = run_states.setdefault(
                        run_key,
                        {
                            "active_days": np.zeros(target_len, dtype=float),
                            "any_days": np.zeros(target_len, dtype=float),
                            "new_active": (
                                np.zeros(target_len, dtype=float)
                                if new_active_available
                                else np.full(target_len, np.nan)
                            ),
                            "mechanism_days": np.zeros((target_len, n_mechanisms), dtype=float),
                        },
                    )
                    state["active_days"] = state["active_days"] + _sf4_sum_vector_column(
                        group,
                        "active_infection_days_by_bacteria",
                        target_len,
                    )
                    state["any_days"] = state["any_days"] + _sf4_sum_vector_column(
                        group,
                        "infection_days_with_any_resistance_mechanism_by_bacteria",
                        target_len,
                    )
                    if new_active_available:
                        state["new_active"] = state["new_active"] + _sf4_sum_vector_column(
                            group,
                            "new_active_infections_by_bacteria",
                            target_len,
                        )
                    try:
                        exact_sums = group[exact_columns].sum(axis=0, skipna=True).to_numpy(dtype=float)
                    except (TypeError, ValueError):
                        exact_sums = (
                            group[exact_columns]
                            .apply(pd.to_numeric, errors="coerce")
                            .sum(axis=0, skipna=True)
                            .to_numpy(dtype=float)
                        )
                    state["mechanism_days"] = state["mechanism_days"] + exact_sums.reshape(
                        (target_len, n_mechanisms)
                    )
        except (FileNotFoundError, ValueError, OSError, pd.errors.EmptyDataError) as exc:
            problems.append(f"{csv_path.name}: could not read required SX columns ({exc}).")
            continue

        if matched_rows == 0 or not run_states:
            problems.append(f"{csv_path.name}: no baseline-policy rows in 2022-2025.")
            continue

        for run_id, state in run_states.items():
            active_days = np.asarray(state["active_days"], dtype=float)
            any_days = np.asarray(state["any_days"], dtype=float)
            new_active = np.asarray(state["new_active"], dtype=float)
            mechanism_days = np.asarray(state["mechanism_days"], dtype=float)
            for b_idx, bacterium_slug in enumerate(_F15_KNOWN_BACTERIA_SLUGS):
                row: dict[str, object] = {
                    "source": csv_path.name,
                    "run": str(run_id),
                    "bacterium_idx": b_idx,
                    "bacterium_slug": bacterium_slug,
                    "bacterium": _figure_15_bacterium_label(bacterium_slug),
                    "active_infection_days": active_days[b_idx],
                    "new_active_infections": new_active[b_idx],
                    "new_active_infections_available": new_active_available,
                    "any_mechanism_days": any_days[b_idx],
                    "any_mechanism_percent": _sf4_percent(any_days[b_idx], active_days[b_idx]),
                }
                for mech_idx, mechanism in enumerate(_SF4_EXACT_MECHANISMS):
                    slug = str(mechanism["slug"])
                    days = float(mechanism_days[b_idx, mech_idx])
                    row[f"{slug}_days"] = days
                    row[f"{slug}_percent"] = _sf4_percent(days, active_days[b_idx])
                rows.append(row)
    return rows, problems, any_new_active_available


def _sf4_median_numeric(values: pd.Series) -> float:
    numeric = pd.to_numeric(values, errors="coerce").dropna().to_numpy(dtype=float)
    if len(numeric) == 0:
        return np.nan
    return float(np.nanmedian(numeric))


def _sf4_reliability_flag(row: dict[str, object]) -> str:
    flags: list[str] = []
    active_days = float(row.get("active_infection_days", np.nan))
    new_active = float(row.get("new_active_infections", np.nan))
    any_days = float(row.get("any_mechanism_days", np.nan))
    exact_counts = [
        float(row.get(f"{mechanism['slug']}_days", np.nan))
        for mechanism in _SF4_EXACT_MECHANISMS
    ]
    finite_exact_counts = [value for value in exact_counts if np.isfinite(value)]
    max_exact_count = max(finite_exact_counts) if finite_exact_counts else np.nan

    if not np.isfinite(active_days) or active_days <= 0.0:
        flags.append("no active infection-days")
    elif active_days < 1000.0:
        flags.append("low active infection-days")

    if not bool(row.get("new_active_infections_available", False)) or not np.isfinite(new_active):
        flags.append("incident denominator unavailable")
    elif new_active < 100.0:
        flags.append("low incident infection denominator")
    elif new_active < 500.0:
        flags.append("moderate incident infection denominator")

    if np.isfinite(any_days) and any_days > 0.0 and np.isfinite(max_exact_count) and max_exact_count < 20.0:
        flags.append("sparse mechanism counts")

    return "; ".join(flags) if flags else "\u2014"


def _sf4_summarise(rows: list[dict[str, object]]) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    records: list[dict[str, object]] = []
    for (bacterium_idx, bacterium), group in df.groupby(["bacterium_idx", "bacterium"], dropna=False):
        run_pairs = group[["source", "run"]].drop_duplicates()
        record: dict[str, object] = {
            "bacterium_idx": int(bacterium_idx),
            "bacterium": str(bacterium),
            "active_infection_days": _sf4_median_numeric(group["active_infection_days"]),
            "new_active_infections": _sf4_median_numeric(group["new_active_infections"]),
            "new_active_infections_available": bool(group["new_active_infections_available"].any()),
            "any_mechanism_days": _sf4_median_numeric(group["any_mechanism_days"]),
            "any_mechanism_percent": _sf4_median_numeric(group["any_mechanism_percent"]),
            "n_runs": int(run_pairs.shape[0]),
        }
        for mechanism in _SF4_EXACT_MECHANISMS:
            slug = str(mechanism["slug"])
            record[f"{slug}_days"] = _sf4_median_numeric(group[f"{slug}_days"])
            record[f"{slug}_percent"] = _sf4_median_numeric(group[f"{slug}_percent"])
        record["reliability_flag"] = _sf4_reliability_flag(record)
        records.append(record)
    summary = pd.DataFrame(records)
    return summary.sort_values(
        ["any_mechanism_percent", "active_infection_days"],
        ascending=[False, False],
        na_position="last",
    ).reset_index(drop=True)


def _sf4_format_percent(value: object) -> str:
    if value is None or pd.isna(value):
        return "\u2014"
    value_f = float(value)
    if not np.isfinite(value_f):
        return "\u2014"
    if value_f > 0.0 and value_f < 0.1:
        return f"{value_f:.3f}"
    return f"{value_f:.1f}"


def _sf4_format_count(value: object) -> str:
    if value is None or pd.isna(value):
        return "\u2014"
    value_f = float(value)
    if not np.isfinite(value_f):
        return "\u2014"
    return f"{int(np.rint(value_f)):,}"


def _sf4_axis_short_label(label: str) -> str:
    replacements = {
        "16S rRNA methyltransferase": "16S rRNA\nmethyl.",
        "23S rRNA oxazolidinone": "23S rRNA\noxazol.",
        "23S rRNA macrolide": "23S rRNA\nmacrol.",
        "polymyxin regulatory": "polymyxin\nreg.",
        "global efflux pump": "global\nefflux",
        "global porin loss": "global\nporin",
        "gyrA/parC secondary": "gyrA/parC\nsecondary",
        "AmpC derepression": "AmpC\nderepr.",
        "Acinetobacter OXA": "Acinet.\nOXA",
        "as-yet-unknown": "as-yet\nunknown",
    }
    return replacements.get(label, label.replace(" ", "\n", 1))


def make_supplementary_figure_s4_resistance_mechanisms_by_bacterium(
    csv_paths: list[Path],
    out_dir: Path,
    agg: dict | None = None,
) -> dict[str, object]:
    rows, problems, any_new_active_available = _sf4_rows_from_csvs(csv_paths)
    summary = _sf4_summarise(rows)
    if summary.empty:
        message = _SF4_REQUIRED_MESSAGE
        if problems:
            message += " Parser notes: " + " ".join(problems[:3])
        print("  Supplementary Figure SX: placeholder (required exact mechanism fields missing).")
        return _sf4_placeholder(out_dir, agg, message)

    mechanism_slugs = [str(mechanism["slug"]) for mechanism in _SF4_EXACT_MECHANISMS]
    mechanism_labels = [str(mechanism["label"]) for mechanism in _SF4_EXACT_MECHANISMS]
    heatmap = summary[[f"{slug}_percent" for slug in mechanism_slugs]].to_numpy(dtype=float)
    finite_values = heatmap[np.isfinite(heatmap)]
    positive_values = finite_values[finite_values > 0.0]
    if positive_values.size:
        vmax = max(1.0, min(100.0, float(np.nanpercentile(positive_values, 95))))
    else:
        vmax = 1.0

    fig_height = max(9.0, 2.0 + 0.24 * len(summary))
    fig_width = max(18.0, 0.42 * len(mechanism_slugs))
    fig, ax = plt.subplots(figsize=(fig_width, fig_height), constrained_layout=True)
    cmap = plt.cm.YlOrRd.copy()
    cmap.set_bad("#F2F2F2")
    image = ax.imshow(heatmap, aspect="auto", cmap=cmap, vmin=0.0, vmax=vmax)
    ax.set_xticks(np.arange(len(mechanism_labels)))
    ax.set_xticklabels([_sf4_axis_short_label(label) for label in mechanism_labels], rotation=70, ha="right", fontsize=6.2)
    ax.set_yticks(np.arange(len(summary)))
    ax.set_yticklabels(summary["bacterium"].values, fontsize=6.8, fontstyle="italic")
    ax.set_title(
        "Exact ResistanceMechanism prevalence among active infection-days",
        loc="left",
        fontsize=10.5,
        fontweight="bold",
    )
    ax.set_xlabel("ResistanceMechanism enum variant", fontsize=9)
    ax.set_ylabel("Bacterium", fontsize=9)
    cbar = fig.colorbar(image, ax=ax, shrink=0.72)
    cbar.set_label("% of active infection-days", fontsize=8.5)
    fig.suptitle(_SF4_TITLE, fontsize=11, fontweight="bold")

    table_data: dict[str, object] = {
        "Bacterium": summary["bacterium"],
    }
    if bool(summary["new_active_infections_available"].any()):
        table_data["New active infections"] = summary["new_active_infections"].map(_sf4_format_count)
    table_data["Active infection-days"] = summary["active_infection_days"].map(_sf4_format_count)
    table_data["Infection-days with any recorded resistance mechanism"] = summary[
        "any_mechanism_days"
    ].map(_sf4_format_count)
    table_data["Any recorded mechanism (%)"] = summary["any_mechanism_percent"].map(_sf4_format_percent)
    for mechanism in _SF4_EXACT_MECHANISMS:
        slug = str(mechanism["slug"])
        label = str(mechanism["label"])
        table_data[f"{label} (%)"] = summary[f"{slug}_percent"].map(_sf4_format_percent)
    table_data["Reliability flag"] = summary["reliability_flag"]
    table_data["Runs contributing"] = summary["n_runs"].astype(int).astype(str)

    extra_html = "<h2>Exact mechanism prevalence by bacterium</h2>\n"
    extra_html += _html_table(pd.DataFrame(table_data))
    extra_html += _sf4_mechanism_definitions_table_html()

    n_runs = int(max(summary["n_runs"].max(), 1))
    footnotes = [
        "This figure shows modelled, biologically motivated resistance-mechanism state in the simulation. It is not an externally observed genomic surveillance estimate.",
        "Mechanisms are the exact ResistanceMechanism enum variants used by the model, not grouped mechanism families.",
        "Percentages are calculated among active infection-days for each bacterium in the 2022-2025 baseline-policy window.",
        "Active infection-days are repeated model observations within infection episodes, not independent clinical isolates.",
        "Mechanism columns are not mutually exclusive; one infection-day can carry more than one mechanism.",
        "Mechanism state can arise through the model's resistance pathways, including acquisition from circulating profiles, microbiome/carriage inheritance, HGT where allowed, environmental or ratchet floors, and de novo emergence under drug pressure.",
        "The figure is intended as model-transparency output. The primary resistance calibration result remains the phenotypic bacterium-drug benchmark shown in Figure 2 and Supplementary Table S2.",
        "Rare bacteria may have small incident infection denominators even when active infection-days are nonzero; those rows should be interpreted cautiously.",
        f"Percentages are calculated within each run first and shown as medians across {n_runs} simulation "
        f"run{'s' if n_runs != 1 else ''}; intervals are omitted to keep the heatmap readable.",
        "Reliability flags identify rows with low active infection-days, low or moderate incident infection denominators, sparse exact mechanism counts, or unavailable incident denominators.",
    ]
    if not any_new_active_available:
        footnotes.append("new_active_infections_by_bacteria was unavailable; incident infection denominator flags were therefore unavailable.")
    if problems:
        footnotes.append("Parser notes: " + " ".join(problems[:4]))

    _save_figure(
        fig,
        out_dir,
        _SF4_STEM,
        _SF4_TITLE,
        "Baseline-policy exact resistance-mechanism prevalence among active infection-days.",
        footnotes,
        agg=agg,
        extra_html=extra_html,
    )
    print(
        "  Supplementary Figure SX: real data; "
        f"{len(summary)} bacterium row{'s' if len(summary) != 1 else ''}, "
        f"{len(_SF4_EXACT_MECHANISMS)} exact mechanism columns, "
        f"new active infection denominators {'available' if any_new_active_available else 'unavailable'}."
    )
    return {
        "generated": "real data",
        "bacteria_included": int(len(summary)),
        "mechanisms_included": int(len(_SF4_EXACT_MECHANISMS)),
        "new_active_infections_available": bool(any_new_active_available),
        "n_runs": n_runs,
    }


_SF5_TITLE = "Supplementary Figure S5. Diagnostic testing and targeted-treatment cascade, 2022\u20132025"
_SF5_STEM = "Supplementary_Figure_S5__diagnostic_testing_targeted_treatment_cascade"
_SF5_REQUIRED_MESSAGE = (
    "Supplementary Figure S5 requires simulation_summary aggregate columns for eligible "
    "symptomatic infections, bacterial identification, resistance testing, targeted treatment, "
    "and effective targeted treatment."
)
_SF5_STAGES: list[dict[str, str]] = [
    {
        "key": "eligible",
        "label": "Eligible symptomatic infection",
        "column": "diagnostic_cascade_eligible_symptomatic_infections",
        "definition": "First model day when an active bacterium-specific infection is symptomatic and eligible for diagnostic testing logic.",
    },
    {
        "key": "id",
        "label": "Bacterial identification done",
        "column": "diagnostic_cascade_bacterial_identification_done",
        "definition": "The model's bacterium-specific identification flag has become true for the infection episode.",
    },
    {
        "key": "ast",
        "label": "Resistance testing initiated",
        "column": "diagnostic_cascade_resistance_testing_done",
        "definition": "The Rust flag records resistance/susceptibility testing initiation; this is labelled as initiation even though the compatibility column name says done.",
    },
    {
        "key": "targeted",
        "label": "Targeted antibiotic treatment started",
        "column": "diagnostic_cascade_targeted_treatment_started",
        "definition": "A course-start targeted-context antibiotic is observed for the identified active infection.",
    },
    {
        "key": "effective",
        "label": "Effective targeted antibiotic treatment started",
        "column": "diagnostic_cascade_effective_targeted_treatment_started",
        "definition": "A targeted-context active antibiotic reaches activity_r >= 0.500 for the bacterium.",
    },
]
_SF5_SETTINGS = [
    ("Overall", ""),
    ("Community-onset / cascade-entry", "_community"),
    ("Hospital-onset / cascade-entry", "_hospital"),
]


def _sf5_stage_columns_for_suffix(suffix: str) -> list[str]:
    return [f"{stage['column']}{suffix}" for stage in _SF5_STAGES]


def _sf5_definitions_table_html() -> str:
    definitions = pd.DataFrame({
        "Item": [
            "Denominator",
            "Time window",
            "Data source",
            "Community/hospital split",
            "Resistance testing stage",
            "Effective targeted therapy",
            "Timing semantics",
            "Organism exclusions",
            "Exclusions and caveats",
        ],
        "Definition": [
            "Eligible symptomatic infection episodes entering the diagnostic cascade.",
            "Baseline-policy cascade-entry years 2022-2025, using rows from simulation_summary_run#.csv.",
            "Aggregate simulation_summary_run#.csv columns only; calibration_summary_*.txt is not reader-facing for this figure.",
            "Hospital status is captured at cascade entry. Community means not hospitalized; hospital means hospitalized.",
            "The Rust flag records resistance/susceptibility testing initiation, not confirmed result availability.",
            "A targeted-context active antibiotic with resistance-adjusted activity_r >= 0.500 for the bacterium.",
            "All downstream stages are assigned back to the episode's cascade-entry timestep/year. The model uses daily timesteps, so same-day order is not a precise sub-day clinical timestamp.",
            "The Rust aggregate excludes organisms returned by is_microbiome_excluded (currently treponema_pallidum). H. pylori is not separately excluded by that helper in the current code.",
            "Counts are raw simulated counts. Episodes that resolve, die, or are censored before a later stage remain in the earlier-stage denominator and are not imputed into later stages.",
        ],
    })
    return "<h2>Definitions, Denominators, and Exclusions</h2>\n" + _html_table(definitions)


def _sf5_stage_definition_table_html() -> str:
    stages = pd.DataFrame({
        "Cascade stage": [str(stage["label"]) for stage in _SF5_STAGES],
        "simulation_summary column": [str(stage["column"]) for stage in _SF5_STAGES],
        "Definition": [str(stage["definition"]) for stage in _SF5_STAGES],
    })
    return "<h2>Stage Definitions</h2>\n" + _html_table(stages)


def _sf5_placeholder(out_dir: Path, agg: dict | None, message: str) -> None:
    fig, ax = plt.subplots(figsize=(10.5, 3.8))
    _sf2_axis_placeholder(ax, _SF5_TITLE, message)
    fig.subplots_adjust(left=0.03, right=0.97, top=0.88, bottom=0.08)
    footnotes = [
        "Supplementary Figure S5 is generated from aggregate simulation_summary_run#.csv columns only.",
        "Old simulation_summary files without the diagnostic-cascade aggregate columns intentionally render this placeholder.",
    ]
    _save_figure(
        fig,
        out_dir,
        _SF5_STEM,
        _SF5_TITLE,
        message,
        footnotes,
        agg=agg,
        extra_html=_sf5_stage_definition_table_html() + _sf5_definitions_table_html(),
    )


def _sf5_percent(numerator: float, denominator: float) -> float:
    if not np.isfinite(numerator) or not np.isfinite(denominator) or denominator <= 0.0:
        return np.nan
    return 100.0 * float(numerator) / float(denominator)


def _sf5_format_count(value: object) -> str:
    if value is None or pd.isna(value):
        return "\u2014"
    value_f = float(value)
    if not np.isfinite(value_f):
        return "\u2014"
    return f"{int(np.rint(value_f)):,}"


def _sf5_format_percent(value: object) -> str:
    if value is None or pd.isna(value):
        return "\u2014"
    value_f = float(value)
    if not np.isfinite(value_f):
        return "\u2014"
    return f"{value_f:.1f}"


def _sf5_reliability_flag(eligible: float, count: float, previous_count: float | None) -> str:
    flags: list[str] = []
    if not np.isfinite(eligible) or eligible <= 0.0:
        flags.append("zero eligible denominator")
    elif eligible < 100.0:
        flags.append("low eligible denominator")
    if previous_count is not None and np.isfinite(previous_count) and np.isfinite(count):
        if count > previous_count:
            flags.append("non-monotonic cascade count")
    return "; ".join(flags)


def _sf5_rows_from_csvs(csv_paths: list[Path]) -> tuple[list[dict[str, object]], list[str]]:
    rows: list[dict[str, object]] = []
    problems: list[str] = []
    required = ["time_in_years", "policy_option"] + [
        column
        for _, suffix in _SF5_SETTINGS
        for column in _sf5_stage_columns_for_suffix(suffix)
    ]
    optional = ["run_id", "simulation_year", "year", "time_step"]

    for csv_path in csv_paths:
        header = _simulation_csv_column_names(csv_path)
        if header is None:
            problems.append(f"{csv_path.name}: could not read simulation CSV header.")
            continue
        missing = [column for column in required if column not in header]
        if missing:
            problems.append(f"{csv_path.name}: missing {', '.join(missing[:4])}.")
            continue

        wanted = set(required + optional)
        try:
            df = _read_csv_selected(csv_path, wanted)
        except (FileNotFoundError, ValueError, OSError) as exc:
            problems.append(f"{csv_path.name}: could not read required S5 columns ({exc}).")
            continue
        if df.empty:
            problems.append(f"{csv_path.name}: simulation CSV has no rows.")
            continue

        if "policy_option" in df.columns:
            policy = pd.to_numeric(df["policy_option"], errors="coerce")
            df = df[policy.eq(0)].copy()
        years = _simulation_year_series(df)
        df = df[(years >= 2022.0) & (years < 2026.0)].copy()
        if df.empty:
            problems.append(f"{csv_path.name}: no baseline-policy cascade-entry rows in 2022-2025.")
            continue
        if "run_id" not in df.columns:
            df["run_id"] = csv_path.stem

        for run_id, group in df.groupby("run_id", dropna=False):
            for setting_label, suffix in _SF5_SETTINGS:
                stage_columns = _sf5_stage_columns_for_suffix(suffix)
                counts = [
                    float(pd.to_numeric(group[column], errors="coerce").fillna(0.0).sum())
                    for column in stage_columns
                ]
                eligible = counts[0] if counts else 0.0
                for stage_idx, stage in enumerate(_SF5_STAGES):
                    previous_count = counts[stage_idx - 1] if stage_idx > 0 else None
                    count = counts[stage_idx]
                    rows.append({
                        "source": csv_path.name,
                        "run": str(run_id),
                        "setting": setting_label,
                        "stage_idx": stage_idx,
                        "stage": str(stage["label"]),
                        "definition": str(stage["definition"]),
                        "count": count,
                        "eligible_denominator": eligible,
                        "previous_denominator": np.nan if previous_count is None else previous_count,
                        "pct_of_eligible": _sf5_percent(count, eligible),
                        "pct_of_previous_stage": (
                            np.nan if previous_count is None else _sf5_percent(count, previous_count)
                        ),
                        "reliability_flag": _sf5_reliability_flag(eligible, count, previous_count),
                    })
    return rows, problems


def _sf5_summarise(rows: list[dict[str, object]]) -> pd.DataFrame:
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    records: list[dict[str, object]] = []
    for (setting, stage_idx, stage), group in df.groupby(
        ["setting", "stage_idx", "stage"], dropna=False
    ):
        record = {
            "setting": str(setting),
            "stage_idx": int(stage_idx),
            "stage": str(stage),
            "definition": str(group["definition"].iloc[0]),
            "count": float(np.nanmedian(pd.to_numeric(group["count"], errors="coerce"))),
            "eligible_denominator": float(np.nanmedian(pd.to_numeric(group["eligible_denominator"], errors="coerce"))),
            "previous_denominator": float(np.nanmedian(pd.to_numeric(group["previous_denominator"], errors="coerce"))),
            "pct_of_eligible": float(np.nanmedian(pd.to_numeric(group["pct_of_eligible"], errors="coerce"))),
            "pct_of_previous_stage": float(np.nanmedian(pd.to_numeric(group["pct_of_previous_stage"], errors="coerce"))),
            "n_runs": int(group["run"].nunique()),
        }
        flags = sorted(
            {
                str(flag)
                for flag in group["reliability_flag"].dropna()
                if str(flag).strip()
            }
        )
        record["reliability_flag"] = "; ".join(flags)
        records.append(record)
    return pd.DataFrame(records).sort_values(["setting", "stage_idx"]).reset_index(drop=True)


def make_supplementary_figure_s5_diagnostic_testing_targeted_treatment_cascade(
    csv_paths: list[Path],
    out_dir: Path,
    agg: dict | None = None,
) -> None:
    rows, problems = _sf5_rows_from_csvs(csv_paths)
    summary = _sf5_summarise(rows)
    if summary.empty:
        message = _SF5_REQUIRED_MESSAGE
        if problems:
            message += " Parser notes: " + " ".join(problems[:3])
        _sf5_placeholder(out_dir, agg, message)
        print("  Supplementary Figure S5: placeholder (required diagnostic-cascade fields missing).")
        return

    overall = summary[summary["setting"] == "Overall"].sort_values("stage_idx")
    by_setting = summary[summary["setting"] != "Overall"].sort_values(["setting", "stage_idx"])
    fig, axes = plt.subplots(1, 2, figsize=(15.5, 6.4), constrained_layout=True)

    y = np.arange(len(overall))
    axes[0].barh(y, overall["pct_of_eligible"].to_numpy(float), color="#2A9D8F")
    axes[0].set_yticks(y)
    axes[0].set_yticklabels(overall["stage"].values, fontsize=8.5)
    axes[0].invert_yaxis()
    axes[0].set_xlim(0, 100)
    axes[0].set_xlabel("% of eligible symptomatic infection episodes", fontsize=9)
    axes[0].set_title("A. Overall cascade", loc="left", fontsize=10, fontweight="bold")
    axes[0].spines[["top", "right"]].set_visible(False)
    axes[0].grid(axis="x", linewidth=0.35, alpha=0.45)
    for idx, value in enumerate(overall["pct_of_eligible"].to_numpy(float)):
        if np.isfinite(value):
            axes[0].text(min(value + 1.0, 99.0), idx, f"{value:.1f}%", va="center", fontsize=7.5)

    stage_labels = [str(stage["label"]) for stage in _SF5_STAGES]
    x = np.arange(len(stage_labels))
    width = 0.36
    colours = {
        "Community-onset / cascade-entry": "#4C78A8",
        "Hospital-onset / cascade-entry": "#F28E2B",
    }
    for offset_idx, setting_label in enumerate(
        ["Community-onset / cascade-entry", "Hospital-onset / cascade-entry"]
    ):
        setting_rows = by_setting[by_setting["setting"] == setting_label].sort_values("stage_idx")
        values = setting_rows["pct_of_eligible"].to_numpy(float)
        axes[1].bar(
            x + (offset_idx - 0.5) * width,
            values,
            width=width,
            label=setting_label.split(" / ")[0],
            color=colours[setting_label],
        )
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(stage_labels, rotation=35, ha="right", fontsize=8)
    axes[1].set_ylim(0, 100)
    axes[1].set_ylabel("% of setting-specific eligible episodes", fontsize=9)
    axes[1].set_title("B. Cascade by entry setting", loc="left", fontsize=10, fontweight="bold")
    axes[1].spines[["top", "right"]].set_visible(False)
    axes[1].grid(axis="y", linewidth=0.35, alpha=0.45)
    axes[1].legend(frameon=False, fontsize=8)
    fig.suptitle(_SF5_TITLE, fontsize=11, fontweight="bold")

    overall_table = pd.DataFrame({
        "Stage": overall["stage"],
        "Definition": overall["definition"],
        "Raw simulated count": overall["count"].map(_sf5_format_count),
        "% of eligible episodes": overall["pct_of_eligible"].map(_sf5_format_percent),
        "% of previous stage": overall["pct_of_previous_stage"].map(_sf5_format_percent),
        "Denominator": np.where(
            overall["stage_idx"].eq(0),
            "Eligible symptomatic infection episodes",
            "Previous cascade stage",
        ),
        "Reliability flag": overall["reliability_flag"],
    })
    setting_table = pd.DataFrame({
        "Setting at cascade entry": by_setting["setting"],
        "Stage": by_setting["stage"],
        "Raw simulated count": by_setting["count"].map(_sf5_format_count),
        "% of setting-specific eligible episodes": by_setting["pct_of_eligible"].map(_sf5_format_percent),
        "% of previous stage": by_setting["pct_of_previous_stage"].map(_sf5_format_percent),
        "Eligible denominator": by_setting["eligible_denominator"].map(_sf5_format_count),
        "Reliability flag": by_setting["reliability_flag"],
    })
    extra_html = "<h2>Overall Cascade</h2>\n" + _html_table(overall_table)
    extra_html += "<h2>Cascade by Hospital/Community Status</h2>\n" + _html_table(setting_table)
    extra_html += _sf5_stage_definition_table_html()
    extra_html += _sf5_definitions_table_html()

    n_runs = int(max(summary["n_runs"].max(), 1))
    footnotes = [
        "Denominator: eligible symptomatic bacterium-specific infection episodes entering the diagnostic cascade under baseline policy with cascade-entry years 2022-2025.",
        "Counts are raw simulated counts. When multiple runs are supplied, counts and percentages are calculated within each run first and shown as medians across runs.",
        "All stages after eligibility are assigned back to the episode's cascade-entry timestep/year, so the 2022-2025 filter is based on cascade-entry year.",
        "Community and hospital groups use hospital status at cascade entry; community means not hospitalized and hospital means hospitalized.",
        "The resistance-testing stage is labelled as initiation because the Rust flag records AST/resistance testing initiation, not final result availability.",
        "Effective targeted therapy means a targeted-context active antibiotic with activity_r >= 0.500 for the bacterium, matching the sepsis effective-therapy threshold.",
        "The model uses daily timesteps, so same-model-day ordering is not a precise sub-day clinical timestamp.",
        "Data source: aggregate simulation_summary_run#.csv fields only. calibration_summary_*.txt is not used as a reader-facing substitute for this figure.",
        "The Rust aggregate excludes organisms returned by is_microbiome_excluded (currently treponema_pallidum). H. pylori is not separately excluded by that helper in the current code.",
        f"Values shown are medians across {n_runs} simulation run{'s' if n_runs != 1 else ''}.",
    ]
    if problems:
        footnotes.append("Parser notes: " + " ".join(problems[:4]))

    _save_figure(
        fig,
        out_dir,
        _SF5_STEM,
        _SF5_TITLE,
        "Diagnostic testing and targeted-treatment cascade from aggregate simulation summary columns.",
        footnotes,
        agg=agg,
        extra_html=extra_html,
    )
    print(
        "  Supplementary Figure S5: real data; "
        f"{len(overall)} overall stage row{'s' if len(overall) != 1 else ''} "
        f"from {n_runs} run{'s' if n_runs != 1 else ''}."
    )


# ---------------------------------------------------------------------------
# Supplementary Figure S6. New active infection denominators by bacterium
# ---------------------------------------------------------------------------

_SF6_TITLE = "Supplementary Figure S6. New active infection denominators by bacterium, 2022\u20132025"
_SF6_STEM = "Supplementary_Figure_S6__new_active_infection_denominators_by_bacterium"
_SF6_TOTAL_COLUMN = "new_active_infections_by_bacteria"
_SF6_HOSPITAL_REGIONS = _F6_HOSPITAL_REGIONS
_SF6_REQUIRED_MESSAGE = (
    "Supplementary Figure S6 requires total, hospital, and community new active "
    "infection denominators by bacterium for the 2022-2025 serious-R window. "
    "The total denominator column is new_active_infections_by_bacteria; "
    "hospital/community denominators are derived from existing per-region "
    "newly_infected_hospital columns when available."
)
_SF6_SPECIAL_ORGANISM_NOTES = {
    "mdr_mycobacterium_tuberculosis": (
        "special organism; MDR-TB is excluded from serious-R summaries because "
        "rifampicin resistance is definitional in the MDR-TB model"
    ),
}


def _sf6_placeholder(
    out_dir: Path,
    agg: dict | None,
    message: str,
    problems: list[str] | None = None,
) -> None:
    fig, ax = plt.subplots(figsize=(10, 3.8))
    ax.text(
        0.5,
        0.5,
        f"{_SF6_TITLE}\n\n{message}",
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=10.2,
        color="#555",
        bbox=dict(boxstyle="round,pad=0.6", fc="#f5f5f5", ec="#bbb"),
    )
    ax.set_axis_off()
    fig.subplots_adjust(left=0.03, right=0.97, top=0.92, bottom=0.08)
    extra_html = ""
    if problems:
        extra_html += "<h2>Parser notes</h2>\n<ul>\n"
        for problem in problems[:10]:
            extra_html += f"<li>{problem}</li>\n"
        extra_html += "</ul>\n"
    _save_figure(
        fig,
        out_dir,
        _SF6_STEM,
        _SF6_TITLE,
        message,
        [
            "Supplementary Figure S6 requires total, hospital, and community new active "
            "infection denominators by bacterium for the 2022-2025 serious-R window.",
            "calibration_summary_*.txt is an internal diagnostic file and is not required "
            "to interpret this reader-facing page.",
            "Old simulation_summary files remain readable: when required aggregate columns "
            "are absent, this placeholder is generated rather than failing the paper-output build.",
        ],
        agg=agg,
        extra_html=extra_html,
    )


def _sf6_hospital_columns_for_slug(slug: str, available: set[str]) -> list[str]:
    return _f6_hospital_columns_for_slug(slug, available)


def _sf6_sum_series(df: pd.DataFrame, columns: list[str]) -> float:
    if not columns:
        return np.nan
    total = 0.0
    seen = False
    for column in columns:
        if column not in df.columns:
            continue
        total += float(pd.to_numeric(df[column], errors="coerce").sum(skipna=True))
        seen = True
    return total if seen else np.nan


def _sf6_rows_from_simulation_csv(csv_path: Path) -> tuple[list[dict[str, object]], str | None]:
    columns = _simulation_csv_columns(csv_path)
    if columns is None:
        return [], f"{csv_path.name}: unable to read simulation summary header."
    if _SF6_TOTAL_COLUMN not in columns:
        return [], f"{csv_path.name}: missing {_SF6_TOTAL_COLUMN}."

    optional = ["policy_option", "run_id", "simulation_year", "year", "time_in_years", "time_step"]
    hospital_columns: list[str] = []
    for slug in _F15_KNOWN_BACTERIA_SLUGS:
        hospital_columns.extend(_sf6_hospital_columns_for_slug(slug, columns))

    usecols = [_SF6_TOTAL_COLUMN, *optional, *hospital_columns]
    try:
        df = _read_csv_selected(csv_path, usecols)
    except (FileNotFoundError, ValueError, OSError) as exc:
        return [], f"{csv_path.name}: unable to load S6 denominator columns ({exc})."

    if "policy_option" in df.columns:
        df = df[pd.to_numeric(df["policy_option"], errors="coerce") == 0].copy()
    df["sf6_year"] = _simulation_year_series(df)
    df = df[(df["sf6_year"] >= 2022.0) & (df["sf6_year"] < 2026.0)].copy()
    if df.empty:
        return [], f"{csv_path.name}: no baseline-policy rows in 2022-2025."

    grouped = df.groupby("run_id", dropna=False) if "run_id" in df.columns else [(csv_path.stem, df)]
    target_len = len(_F15_KNOWN_BACTERIA_SLUGS)
    rows: list[dict[str, object]] = []
    hospital_source_available = bool(hospital_columns)

    for run_key, run_df in grouped:
        total_values = _st1_sum_vector_column(run_df, _SF6_TOTAL_COLUMN, target_len)
        run_len = max(target_len, len(total_values))
        total_values = _figure_15_extend_array(total_values, run_len)

        for b_idx in range(run_len):
            slug = (
                _F15_KNOWN_BACTERIA_SLUGS[b_idx]
                if b_idx < len(_F15_KNOWN_BACTERIA_SLUGS)
                else f"bacterium_{b_idx + 1}"
            )
            total = float(total_values[b_idx])
            hospital_columns_for_bacterium = _sf6_hospital_columns_for_slug(slug, set(run_df.columns))
            hospital = _sf6_sum_series(run_df, hospital_columns_for_bacterium)
            if np.isfinite(hospital):
                community = total - hospital
                if community < 0 and abs(community) < 1e-9:
                    community = 0.0
                elif community < 0:
                    community = np.nan
            else:
                community = np.nan

            rows.append({
                "source": csv_path.name,
                "run": str(run_key),
                "bacterium_slug": slug,
                "bacterium": _figure_15_bacterium_label(slug),
                "total_new_active_infections": total,
                "hospital_new_active_infections": hospital,
                "community_new_active_infections": community,
                "hospital_source_available": bool(hospital_columns_for_bacterium),
                "all_hospital_source_available": hospital_source_available,
                "sum_mismatch": (
                    bool(np.isfinite(hospital) and np.isfinite(community) and abs((hospital + community) - total) > 0.5)
                ),
            })

    if not rows:
        return [], f"{csv_path.name}: S6 denominator columns contained no usable values."
    return rows, None


def _sf6_bacterium_key(label: object) -> str:
    return str(label or "").strip().lower()


def _sf6_serious_r_context(calibration_paths: list[Union[str, Path]]) -> pd.DataFrame:
    frames: list[pd.DataFrame] = []
    for path in calibration_paths:
        frame, _summary_values = _figure_20_parse_calibration_summary(path)
        if not frame.empty:
            frames.append(frame)
    if not frames:
        return pd.DataFrame()
    all_rows = pd.concat(frames, ignore_index=True)
    summary, _notes = _figure_20_summarise_rows(all_rows)
    if summary.empty:
        return pd.DataFrame()
    return pd.DataFrame({
        "bacterium_key": summary["bacterium"].map(_sf6_bacterium_key),
        "Marker drug(s)": summary["marker_drugs"],
        "Overall serious-R (%)": summary["overall_median"],
        "Hospital serious-R (%)": summary["hospital_median"],
        "Community serious-R (%)": summary["community_median"],
    })


def _sf6_adequacy_label(value: object) -> str:
    if value is None or pd.isna(value):
        return "not available"
    n = float(value)
    if not np.isfinite(n):
        return "not available"
    if n <= 0:
        return "no denominator"
    if n < 20:
        return "very sparse"
    if n < 100:
        return "sparse"
    if n < 500:
        return "moderate"
    return "larger denominator"


def _sf6_adequacy_code(value: object) -> int:
    label = _sf6_adequacy_label(value)
    return {
        "not available": 0,
        "no denominator": 1,
        "very sparse": 2,
        "sparse": 3,
        "moderate": 4,
        "larger denominator": 5,
    }[label]


def _sf6_format_count(value: object) -> str:
    if value is None or pd.isna(value):
        return "\u2014"
    value_f = float(value)
    if not np.isfinite(value_f):
        return "\u2014"
    return f"{value_f:,.0f}"


def _sf6_format_interval(p5: object, p95: object) -> str:
    if p5 is None or p95 is None or pd.isna(p5) or pd.isna(p95):
        return "\u2014"
    p5_f = float(p5)
    p95_f = float(p95)
    if not np.isfinite(p5_f) or not np.isfinite(p95_f):
        return "\u2014"
    return f"{p5_f:,.0f}-{p95_f:,.0f}"


def _sf6_format_percent(value: object) -> str:
    if value is None or pd.isna(value):
        return "\u2014"
    value_f = float(value)
    if not np.isfinite(value_f):
        return "\u2014"
    return f"{value_f:.1f}"


def _sf6_run_notes(row: pd.Series) -> str:
    notes: list[str] = []
    slug = str(row.get("bacterium_slug", ""))
    if slug in _SF6_SPECIAL_ORGANISM_NOTES:
        notes.append(_SF6_SPECIAL_ORGANISM_NOTES[slug])
    if bool(row.get("sum_mismatch", False)):
        notes.append("hospital/community sum mismatch")
    if not bool(row.get("hospital_source_available_any", False)):
        notes.append("hospital/community new-infection denominators not available")
    for stratum, pct_col, n_col in [
        ("hospital", "Hospital serious-R (%)", "hospital_new_active_infections_median"),
        ("community", "Community serious-R (%)", "community_new_active_infections_median"),
    ]:
        pct = row.get(pct_col)
        n = row.get(n_col)
        if (
            pct is not None
            and pd.notna(pct)
            and np.isfinite(float(pct))
            and (float(pct) <= 0.0 or float(pct) >= 100.0)
            and _sf6_adequacy_label(n) in {"no denominator", "very sparse", "sparse"}
        ):
            notes.append(f"apparent {stratum} {float(pct):.0f}% serious-R with sparse denominator")
    return "; ".join(dict.fromkeys(notes)) or "\u2014"


def _sf6_summarise_rows(rows: list[dict[str, object]], calibration_paths: list[Union[str, Path]]) -> pd.DataFrame:
    raw = pd.DataFrame(rows)
    if raw.empty:
        return pd.DataFrame()

    summary_rows: list[dict[str, object]] = []
    for (slug, bacterium), group in raw.groupby(["bacterium_slug", "bacterium"], sort=False):
        row: dict[str, object] = {
            "bacterium_slug": slug,
            "bacterium": bacterium,
            "bacterium_key": _sf6_bacterium_key(bacterium),
            "runs_contributing": int(group[["source", "run"]].drop_duplicates().shape[0]),
            "hospital_source_available_any": bool(group["hospital_source_available"].any()),
            "sum_mismatch": bool(group["sum_mismatch"].any()),
        }
        for col in [
            "total_new_active_infections",
            "hospital_new_active_infections",
            "community_new_active_infections",
        ]:
            values = pd.to_numeric(group[col], errors="coerce").dropna().astype(float).to_numpy()
            if len(values) == 0:
                row[f"{col}_median"] = np.nan
                row[f"{col}_p5"] = np.nan
                row[f"{col}_p95"] = np.nan
            else:
                row[f"{col}_median"] = float(np.nanmedian(values))
                row[f"{col}_p5"] = float(np.nanpercentile(values, 5))
                row[f"{col}_p95"] = float(np.nanpercentile(values, 95))
        total = row["total_new_active_infections_median"]
        hospital = row["hospital_new_active_infections_median"]
        if pd.notna(total) and np.isfinite(float(total)) and float(total) > 0 and pd.notna(hospital):
            row["percent_hospital_median"] = 100.0 * float(hospital) / float(total)
        else:
            row["percent_hospital_median"] = np.nan
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    context = _sf6_serious_r_context(calibration_paths)
    if not context.empty:
        summary = summary.merge(context, on="bacterium_key", how="left")
    else:
        summary["Marker drug(s)"] = "\u2014"
        summary["Overall serious-R (%)"] = np.nan
        summary["Hospital serious-R (%)"] = np.nan
        summary["Community serious-R (%)"] = np.nan

    summary["Marker drug(s)"] = summary["Marker drug(s)"].fillna("\u2014")
    summary["Notes"] = summary.apply(_sf6_run_notes, axis=1)
    summary = summary.sort_values(
        ["total_new_active_infections_median", "bacterium"],
        ascending=[False, True],
        na_position="last",
    ).reset_index(drop=True)
    return summary


def _sf6_add_errorbar(
    ax: "plt.Axes",
    summary: pd.DataFrame,
    y: np.ndarray,
    prefix: str,
    label: str,
    color: str,
    marker: str,
    n_runs: int,
) -> None:
    values = summary[f"{prefix}_median"].astype(float).to_numpy()
    mask = np.isfinite(values)
    if not mask.any():
        return
    xerr = None
    if n_runs > 1:
        p5 = summary[f"{prefix}_p5"].astype(float).to_numpy()
        p95 = summary[f"{prefix}_p95"].astype(float).to_numpy()
        err_low = np.where(np.isfinite(p5), np.clip(values - p5, 0.0, None), 0.0)
        err_high = np.where(np.isfinite(p95), np.clip(p95 - values, 0.0, None), 0.0)
        xerr = np.vstack([err_low[mask], err_high[mask]])
    ax.errorbar(
        values[mask],
        y[mask],
        xerr=xerr,
        fmt=marker,
        markersize=4.6,
        color=color,
        ecolor=color,
        elinewidth=0.75,
        capsize=2.2 if n_runs > 1 else 0,
        linestyle="none",
        label=label,
        zorder=3,
    )


def _sf6_apply_count_axis(ax: "plt.Axes", values: list[float]) -> None:
    finite = [float(value) for value in values if np.isfinite(float(value)) and float(value) > 0]
    if not finite:
        ax.set_xlim(0, 1)
        return
    max_value = max(finite)
    min_value = min(finite)
    if max_value / max(min_value, 1.0) >= 100.0:
        ax.set_xscale("symlog", linthresh=1.0)
    ax.set_xlim(0, max_value * 1.18)
    ax.grid(axis="x", linewidth=0.4, alpha=0.45)


def make_supplementary_figure_s6_new_active_infection_denominators(
    csv_paths: list[Path],
    calibration_paths: list[Union[str, Path]],
    out_dir: Path,
    agg: dict | None = None,
) -> dict[str, object]:
    if not csv_paths:
        _sf6_placeholder(
            out_dir,
            agg,
            _SF6_REQUIRED_MESSAGE,
            ["No matching simulation_summary_*.csv files with new_active_infections_by_bacteria were found."],
        )
        print("  Supplementary Figure S6: placeholder (no matching denominator CSVs).")
        return {"generated": "placeholder", "bacteria_included": 0, "n_runs": 0}

    rows: list[dict[str, object]] = []
    problems: list[str] = []
    for csv_path in csv_paths:
        run_rows, problem = _sf6_rows_from_simulation_csv(csv_path)
        rows.extend(run_rows)
        if problem:
            problems.append(problem)

    if not rows:
        _sf6_placeholder(out_dir, agg, _SF6_REQUIRED_MESSAGE, problems)
        print("  Supplementary Figure S6: placeholder (no usable denominator rows).")
        return {"generated": "placeholder", "bacteria_included": 0, "n_runs": 0}

    summary = _sf6_summarise_rows(rows, calibration_paths)
    if summary.empty:
        _sf6_placeholder(out_dir, agg, _SF6_REQUIRED_MESSAGE, problems)
        print("  Supplementary Figure S6: placeholder (empty denominator summary).")
        return {"generated": "placeholder", "bacteria_included": 0, "n_runs": 0}

    n_runs = int(pd.DataFrame(rows)[["source", "run"]].drop_duplicates().shape[0])
    has_hospital = bool(summary["hospital_source_available_any"].any())
    plot_summary = summary.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(plot_summary))

    fig_height = max(8.0, 2.2 + 0.28 * len(plot_summary))
    fig_width = 15.5 if has_hospital else 11.5
    ncols = 3 if has_hospital else 2
    width_ratios = [1.15, 1.35, 0.62] if has_hospital else [1.25, 0.9]
    fig, axes = plt.subplots(
        1,
        ncols,
        figsize=(fig_width, fig_height),
        sharey=True,
        gridspec_kw={"width_ratios": width_ratios, "wspace": 0.08},
    )
    if ncols == 2:
        ax_total, ax_placeholder = axes
        ax_split = None
        ax_heat = ax_placeholder
    else:
        ax_total, ax_split, ax_heat = axes

    _sf6_add_errorbar(
        ax_total,
        plot_summary,
        y,
        "total_new_active_infections",
        "Total",
        "#1565C0",
        "o",
        n_runs,
    )
    _sf6_apply_count_axis(
        ax_total,
        plot_summary["total_new_active_infections_median"].fillna(0).astype(float).tolist(),
    )
    ax_total.set_title("A. Total new active infections", fontsize=9.5, fontweight="bold")
    ax_total.set_xlabel("Raw simulated n, 2022-2025", fontsize=8.8)
    ax_total.set_yticks(y)
    ax_total.set_yticklabels(plot_summary["bacterium"].values, fontsize=7.1, fontstyle="italic")
    ax_total.spines[["top", "right"]].set_visible(False)

    if has_hospital and ax_split is not None:
        for idx, row in plot_summary.iterrows():
            h = row["hospital_new_active_infections_median"]
            c = row["community_new_active_infections_median"]
            if pd.notna(h) and pd.notna(c) and np.isfinite(float(h)) and np.isfinite(float(c)):
                ax_split.plot([float(h), float(c)], [idx, idx], color="#b0bec5", linewidth=0.75, alpha=0.6)
        _sf6_add_errorbar(
            ax_split,
            plot_summary,
            y,
            "hospital_new_active_infections",
            "Hospital-acquired",
            "#C65D00",
            "D",
            n_runs,
        )
        _sf6_add_errorbar(
            ax_split,
            plot_summary,
            y,
            "community_new_active_infections",
            "Community-acquired",
            "#2A9D8F",
            "o",
            n_runs,
        )
        _sf6_apply_count_axis(
            ax_split,
            pd.concat(
                [
                    plot_summary["hospital_new_active_infections_median"],
                    plot_summary["community_new_active_infections_median"],
                ],
                ignore_index=True,
            )
            .fillna(0)
            .astype(float)
            .tolist(),
        )
        ax_split.set_title("B. Hospital and community denominators", fontsize=9.5, fontweight="bold")
        ax_split.set_xlabel("Raw simulated n, 2022-2025", fontsize=8.8)
        ax_split.spines[["top", "right"]].set_visible(False)
        ax_split.tick_params(axis="y", left=False, labelleft=False)
        ax_split.legend(fontsize=7.6, frameon=False, loc="lower right")
    else:
        ax_placeholder.text(
            0.5,
            0.5,
            "Hospital/community denominator\ncolumns were not available.\nPanel A uses total n only.",
            ha="center",
            va="center",
            transform=ax_placeholder.transAxes,
            fontsize=9.2,
            color="#555",
            bbox=dict(boxstyle="round,pad=0.5", fc="#f5f5f5", ec="#bbb"),
        )
        ax_placeholder.set_axis_off()

    if has_hospital:
        heat_values = np.column_stack([
            plot_summary["hospital_new_active_infections_median"].map(_sf6_adequacy_code).to_numpy(),
            plot_summary["community_new_active_infections_median"].map(_sf6_adequacy_code).to_numpy(),
        ])
        cmap = mcolors.ListedColormap([
            "#e0e0e0",
            "#f7f7f7",
            "#d73027",
            "#fc8d59",
            "#fee08b",
            "#1a9850",
        ])
        ax_heat.imshow(heat_values, aspect="auto", interpolation="nearest", cmap=cmap, vmin=0, vmax=5)
        ax_heat.set_title("C. Denominator adequacy", fontsize=9.5, fontweight="bold")
        ax_heat.set_xticks([0, 1])
        ax_heat.set_xticklabels(["Hospital", "Community"], rotation=45, ha="right", fontsize=7.5)
        ax_heat.tick_params(axis="y", left=False, labelleft=False)
        ax_heat.set_yticks(y)
        ax_heat.set_yticklabels([])
        for spine in ax_heat.spines.values():
            spine.set_visible(False)
        legend_items = [
            mpatches.Patch(color="#e0e0e0", label="not available"),
            mpatches.Patch(color="#f7f7f7", label="n = 0"),
            mpatches.Patch(color="#d73027", label="1-19"),
            mpatches.Patch(color="#fc8d59", label="20-99"),
            mpatches.Patch(color="#fee08b", label="100-499"),
            mpatches.Patch(color="#1a9850", label=">=500"),
        ]
        ax_heat.legend(
            handles=legend_items,
            title="Raw n",
            fontsize=6.6,
            title_fontsize=7.2,
            frameon=False,
            loc="lower center",
            bbox_to_anchor=(0.5, -0.01),
        )

    fig.suptitle(_SF6_TITLE, fontsize=10.8, fontweight="bold", y=0.995)
    fig.tight_layout(rect=[0, 0.01, 1, 0.985])

    if n_runs > 1:
        table_df = pd.DataFrame({
            "Bacterium": summary["bacterium"],
            "Serious-R marker drug(s)": summary["Marker drug(s)"],
            "Total n, median": summary["total_new_active_infections_median"].map(_sf6_format_count),
            "Total n, 5th-95th": [
                _sf6_format_interval(p5, p95)
                for p5, p95 in zip(
                    summary["total_new_active_infections_p5"],
                    summary["total_new_active_infections_p95"],
                )
            ],
            "Hospital n, median": summary["hospital_new_active_infections_median"].map(_sf6_format_count),
            "Hospital n, 5th-95th": [
                _sf6_format_interval(p5, p95)
                for p5, p95 in zip(
                    summary["hospital_new_active_infections_p5"],
                    summary["hospital_new_active_infections_p95"],
                )
            ],
            "Community n, median": summary["community_new_active_infections_median"].map(_sf6_format_count),
            "Community n, 5th-95th": [
                _sf6_format_interval(p5, p95)
                for p5, p95 in zip(
                    summary["community_new_active_infections_p5"],
                    summary["community_new_active_infections_p95"],
                )
            ],
        })
    else:
        table_df = pd.DataFrame({
            "Bacterium": summary["bacterium"],
            "Serious-R marker drug(s)": summary["Marker drug(s)"],
            "Total new active infections, raw simulated n": summary[
                "total_new_active_infections_median"
            ].map(_sf6_format_count),
            "Hospital new active infections, raw simulated n": summary[
                "hospital_new_active_infections_median"
            ].map(_sf6_format_count),
            "Community new active infections, raw simulated n": summary[
                "community_new_active_infections_median"
            ].map(_sf6_format_count),
        })

    table_df["Percent hospital"] = summary["percent_hospital_median"].map(_sf6_format_percent)
    table_df["Denominator adequacy flag, hospital"] = summary[
        "hospital_new_active_infections_median"
    ].map(_sf6_adequacy_label)
    table_df["Denominator adequacy flag, community"] = summary[
        "community_new_active_infections_median"
    ].map(_sf6_adequacy_label)
    table_df["Overall serious-R (%)"] = summary["Overall serious-R (%)"].map(_sf6_format_percent)
    table_df["Hospital serious-R (%)"] = summary["Hospital serious-R (%)"].map(_sf6_format_percent)
    table_df["Community serious-R (%)"] = summary["Community serious-R (%)"].map(_sf6_format_percent)
    table_df["Runs contributing"] = summary["runs_contributing"].astype(int).astype(str)
    table_df["Notes"] = summary["Notes"]

    extra_html = (
        "<h2>Detailed new active infection denominators for serious-R interpretation</h2>\n"
        + _html_table(table_df)
    )
    if problems:
        extra_html += "<h2>Parser notes</h2>\n<ul>\n"
        for problem in problems[:10]:
            extra_html += f"<li>{problem}</li>\n"
        extra_html += "</ul>\n"

    run_note = (
        f"Counts are summarised across {n_runs} accepted runs; displayed values are "
        "run-level medians and intervals are 5th-95th percentiles."
        if n_runs > 1
        else "Counts are from the supplied run."
    )
    footnotes = [
        "Counts are summarised over baseline-policy new active infections in the 2022-2025 calibration window.",
        "Data source: aggregate model outputs parsed into paper_tables. calibration_summary_*.txt is not required to interpret this page.",
        "Counts are raw simulated bacterium-specific new active infection events and are provided to interpret bacterium-specific serious-R percentages.",
        "Hospital and community new-infection strata use the infection_hospital_acquired event flag exported by the model: hospital-acquired events are hospital, and all other new active infections are community. The existing Figure 7 serious-R percentages preferentially use current hospital/community infected-stock denominators when those stock-split columns are available; this page shows new active infection event denominators for interpretability and labels that distinction directly.",
        "Counts are raw simulated n values, not population-scaled incidence estimates.",
        "Small denominators can make serious-R percentages unstable; apparent 0% or 100% serious-R values should be interpreted cautiously when denominators are sparse.",
        "Serious-R is a model-defined marker-drug resistance endpoint and is not a directly observed surveillance category.",
        "If multiple bacteria can be involved in one infection episode, bacterium-specific infection events may count separately by bacterium.",
        run_note,
        "Denominator adequacy categories are interpretation aids only: no denominator = n = 0; very sparse = 1-19; sparse = 20-99; moderate = 100-499; larger denominator = >=500.",
    ]

    _save_figure(
        fig,
        out_dir,
        _SF6_STEM,
        _SF6_TITLE,
        "Total and hospital/community raw simulated new active infection denominators by bacterium.",
        footnotes,
        agg=agg,
        extra_html=extra_html,
    )
    generated_status = "real data" if has_hospital else "partial data"
    print(
        "  Supplementary Figure S6: "
        f"{generated_status}; {len(summary)} bacteria included from {n_runs} run"
        f"{'s' if n_runs != 1 else ''}."
    )
    return {
        "generated": generated_status,
        "bacteria_included": len(summary),
        "n_runs": n_runs,
        "hospital_community_available": has_hospital,
    }


# ---------------------------------------------------------------------------
# Supplementary Figure S7. Active infection incidence by bacterium
# ---------------------------------------------------------------------------

_SF7_TITLE = (
    "Supplementary Figure S7. Active infection incidence by bacterium: "
    "simulation versus target, 2022\u20132025"
)
_SF7_STEM = "Supplementary_Figure_S7__active_infection_incidence_by_bacterium"
_SF7_REQUIRED_MESSAGE = (
    "Supplementary Figure S7 requires the Bacteria Burden Benchmarks infection "
    "table in calibration_summary_*.txt, including Bacteria, Infection target "
    "(%), and Infection simulation (%)."
)


def _sf7_placeholder(out_dir: Path, agg: dict | None, message: str) -> dict[str, object]:
    fig, ax = plt.subplots(figsize=(10, 3.8))
    ax.text(
        0.5,
        0.5,
        f"{_SF7_TITLE}\n\n{message}",
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=10.2,
        color="#555",
        bbox=dict(boxstyle="round,pad=0.6", fc="#f5f5f5", ec="#bbb"),
        wrap=True,
    )
    ax.set_axis_off()
    fig.subplots_adjust(left=0.03, right=0.97, top=0.92, bottom=0.08)
    _save_figure(
        fig,
        out_dir,
        _SF7_STEM,
        _SF7_TITLE,
        message,
        [
            "The placeholder is generated so the paper-output build remains complete "
            "when older calibration summaries lack the required burden-benchmark columns.",
        ],
        agg=agg,
    )
    return {"generated": "placeholder", "bacteria_included": 0}


def _sf7_find_column(columns: pd.Index, exact: str, include: list[str], exclude: list[str] | None = None) -> str | None:
    lower_exact = exact.lower()
    for column in columns:
        if str(column).strip().lower() == lower_exact:
            return str(column)
    exclude_terms = [term.lower() for term in (exclude or [])]
    for column in columns:
        text = str(column).strip().lower()
        if all(term.lower() in text for term in include) and not any(term in text for term in exclude_terms):
            return str(column)
    return None


def _sf7_format_percent(value: object) -> str:
    parsed = _first_numeric_value(value)
    if parsed is None or not np.isfinite(parsed):
        return "\u2014"
    value_f = float(parsed)
    if value_f == 0:
        return "0"
    if abs(value_f) < 0.01:
        return f"{value_f:.4f}"
    if abs(value_f) < 1:
        return f"{value_f:.3f}"
    return f"{value_f:.2f}"


def _sf7_format_interval(lo: object, hi: object) -> str:
    lo_text = _sf7_format_percent(lo)
    hi_text = _sf7_format_percent(hi)
    if lo_text == "\u2014" or hi_text == "\u2014":
        return "\u2014"
    return f"{lo_text}-{hi_text}"


def _sf7_format_ratio(value: object) -> str:
    if value is None or pd.isna(value):
        return "\u2014"
    value_f = float(value)
    if not np.isfinite(value_f):
        return "\u2014"
    return f"{value_f:.2f}"


def _sf7_axis_label(value: float, _pos: int) -> str:
    if value == 0:
        return "0"
    if abs(value) < 0.01:
        return f"{value:.4f}".rstrip("0").rstrip(".")
    if abs(value) < 1:
        return f"{value:.2f}".rstrip("0").rstrip(".")
    return f"{value:.1f}".rstrip("0").rstrip(".")


def make_supplementary_figure_s7_active_infection_incidence(
    agg: dict,
    out_dir: Path,
) -> dict[str, object]:
    bi = agg.get("bacteria_infections", pd.DataFrame()).copy()
    n_runs = int(agg.get("n_runs", 1) or 1)
    if bi is None or bi.empty:
        print("  Supplementary Figure S7: placeholder (no bacteria burden table).")
        return _sf7_placeholder(out_dir, agg, _SF7_REQUIRED_MESSAGE)

    bacterium_col = "Bacteria" if "Bacteria" in bi.columns else str(bi.columns[0])
    target_col = _sf7_find_column(
        bi.columns,
        "Infection target (%)",
        ["infection", "target", "%"],
        ["hospital", "carriage", "<5", "65+"],
    )
    simulation_col = _sf7_find_column(
        bi.columns,
        "Infection simulation (%)",
        ["infection", "simulation", "%"],
        ["hospital", "carriage", "<5", "65+"],
    )
    if target_col is None or simulation_col is None:
        missing = []
        if target_col is None:
            missing.append("Infection target (%)")
        if simulation_col is None:
            missing.append("Infection simulation (%)")
        message = _SF7_REQUIRED_MESSAGE + " Missing: " + ", ".join(missing) + "."
        print("  Supplementary Figure S7: placeholder (missing infection columns).")
        return _sf7_placeholder(out_dir, agg, message)

    work = bi[[bacterium_col, target_col, simulation_col]].copy()
    work = _add_interval_columns(work, target_col, "_target")
    work = _add_interval_columns(work, simulation_col, "_simulation")
    work = work.dropna(subset=["_target_med", "_simulation_med"], how="any").copy()
    if work.empty:
        print("  Supplementary Figure S7: placeholder (no numeric infection rows).")
        return _sf7_placeholder(out_dir, agg, "No numeric target/simulation infection rows were found.")

    work["_ratio"] = np.where(
        work["_target_med"].astype(float) > 0.0,
        work["_simulation_med"].astype(float) / work["_target_med"].astype(float),
        np.nan,
    )
    work["_difference"] = work["_simulation_med"].astype(float) - work["_target_med"].astype(float)
    work["_label"] = work[bacterium_col].astype(str).str.strip()
    work["_summary_flag"] = np.where(
        work["_label"].str.endswith("*"),
        "simulation >2x or <0.5x target in calibration summary",
        "\u2014",
    )
    work = work.sort_values(
        ["_target_med", "_simulation_med", "_label"],
        ascending=[False, False, True],
        na_position="last",
    ).reset_index(drop=True)

    plot = work.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(plot))
    fig_height = max(6.5, 2.2 + 0.28 * len(plot))
    fig, ax = plt.subplots(figsize=(10.4, fig_height))

    for idx, row in plot.iterrows():
        target = float(row["_target_med"])
        simulation = float(row["_simulation_med"])
        if np.isfinite(target) and np.isfinite(simulation):
            ax.plot([target, simulation], [idx, idx], color="#b0bec5", linewidth=0.85, zorder=1)

    target_values = plot["_target_med"].to_numpy(dtype=float)
    simulation_values = plot["_simulation_med"].to_numpy(dtype=float)
    target_mask = np.isfinite(target_values)
    simulation_mask = np.isfinite(simulation_values)

    ax.scatter(
        target_values[target_mask],
        y[target_mask] - 0.11,
        marker="D",
        s=25,
        color="#FF7043",
        edgecolors="white",
        linewidths=0.35,
        label="Target",
        zorder=3,
    )

    sim_err_low, sim_err_high = _asymmetric_errors(plot, "_simulation")
    xerr = None
    if n_runs > 1:
        xerr = np.vstack([sim_err_low[simulation_mask], sim_err_high[simulation_mask]])
    ax.errorbar(
        simulation_values[simulation_mask],
        y[simulation_mask] + 0.11,
        xerr=xerr,
        fmt="o",
        markersize=4.5,
        color="#1565C0",
        ecolor="#1565C0",
        elinewidth=0.75,
        capsize=2.2 if n_runs > 1 else 0,
        linestyle="none",
        label="Simulation median" if n_runs > 1 else "Simulation",
        zorder=4,
    )

    finite_values = [
        float(value)
        for value in [*target_values.tolist(), *simulation_values.tolist()]
        if np.isfinite(float(value)) and float(value) >= 0.0
    ]
    positive_values = [value for value in finite_values if value > 0.0]
    axis_note = ""
    if positive_values:
        max_value = max(positive_values)
        min_value = min(positive_values)
        if max_value / max(min_value, 1e-12) >= 100.0:
            linthresh = max(min_value, max_value / 10_000.0)
            ax.set_xscale("symlog", linthresh=linthresh, linscale=0.8)
            axis_note = (
                "The x-axis uses a symmetric log scale near zero because target "
                "infection percentages vary by several orders of magnitude across bacteria."
            )
        ax.set_xlim(left=0.0, right=max_value * 1.25)
    else:
        ax.set_xlim(0.0, 1.0)

    ax.xaxis.set_major_formatter(mticker.FuncFormatter(_sf7_axis_label))
    ax.set_yticks(y)
    ax.set_yticklabels(plot["_label"].values, fontsize=7.1, fontstyle="italic")
    ax.set_xlabel("Active infection incidence (% of world population)", fontsize=9.5)
    ax.set_title(_SF7_TITLE, fontsize=10.8, fontweight="bold", pad=10)
    ax.grid(axis="x", linewidth=0.4, alpha=0.45)
    ax.legend(fontsize=8.2, frameon=False, loc="lower right")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()

    table_cols: dict[str, object] = {
        "Bacterium": work["_label"],
        "Infection target (%)": work["_target_med"].map(_sf7_format_percent),
        "Infection simulation, median (%)": work["_simulation_med"].map(_sf7_format_percent),
    }
    if n_runs > 1:
        table_cols["Infection simulation, 5th-95th (%)"] = [
            _sf7_format_interval(lo, hi)
            for lo, hi in zip(work["_simulation_lo"], work["_simulation_hi"])
        ]
    table_cols["Simulation/target ratio"] = work["_ratio"].map(_sf7_format_ratio)
    table_cols["Simulation minus target (percentage points)"] = work["_difference"].map(_sf7_format_percent)
    table_cols["Calibration summary flag"] = work["_summary_flag"]
    table_df = pd.DataFrame(table_cols)

    extra_html = (
        "<h2>Active infection incidence by bacterium</h2>\n"
        + _html_table(table_df)
    )
    interval_note = (
        "Simulation points are medians across accepted calibration summaries; horizontal "
        "error bars show 5th-95th percentiles."
        if n_runs > 1
        else "Simulation points are from the supplied calibration summary."
    )
    footnotes = [
        "Data source: Bacteria Burden Benchmarks \u2014 Infections & Carriage in "
        "calibration_summary_*.txt. This figure displays Infection target (%) and "
        "Infection simulation (%) only.",
        "The plotted values are calibration-summary active-infection percentages by "
        "bacterium. No new model-output CSV files are used.",
        interval_note,
        "Asterisks in bacterium labels are carried over from the calibration summary and "
        "indicate simulated infection rates greater than 2x or less than 0.5x the target.",
    ]
    if axis_note:
        footnotes.append(axis_note)

    _save_figure(
        fig,
        out_dir,
        _SF7_STEM,
        _SF7_TITLE,
        "Target and simulated active-infection incidence by bacterium from the "
        "Bacteria Burden Benchmarks infection table.",
        footnotes,
        agg=agg,
        extra_html=extra_html,
    )
    print(
        "  Supplementary Figure S7: real data; "
        f"{len(work)} bacteria included from {n_runs} calibration summary "
        f"{'files' if n_runs != 1 else 'file'}."
    )
    return {"generated": "real data", "bacteria_included": len(work), "n_runs": n_runs}


_F20_TITLE = "Figure 7. Serious-R by hospital and community, 2022\u20132025"
_F20_STEM = "Figure_7__serious_r_by_hospital_community"
_F20_REQUIRED_MESSAGE = (
    "Figure 7 requires the Serious Resistance Locus table in calibration_summary_*.txt, "
    "including hospital and community serious-R percentages."
)
_F20_VALUE_COLUMNS = [
    "Overall Serious-R (%)",
    "Hospital Serious-R (%)",
    "Community Serious-R (%)",
    "Sim H:C ratio",
    "Total New Infections",
]


def _figure_20_split_row(line: str) -> list[str]:
    return re.split(r"\s{2,}", line.strip())


def _figure_20_parse_number(value: object) -> float:
    if value is None:
        return np.nan
    text = str(value).strip()
    if text.lower() in {"", "nan", "none", "null", "---", "-", "\u2014"}:
        return np.nan
    try:
        return float(text.replace(",", ""))
    except ValueError:
        return np.nan


def _figure_20_bacterium_label(value: object) -> str:
    text = str(value or "").strip()
    if not text:
        return ""
    slug = re.sub(r"\s+", "_", text.lower().replace("_", " "))
    return _figure_15_bacterium_label(slug)


def _figure_20_parse_calibration_summary(path: Union[str, Path]) -> tuple[pd.DataFrame, dict[str, float]]:
    try:
        lines = _resolve_project_path(path).read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return pd.DataFrame(), {}

    summary: dict[str, float] = {}
    in_summary = False
    table_title_idx: int | None = None

    for idx, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("Serious Resistance Locus Summary"):
            in_summary = True
            continue
        if stripped.startswith("Serious Resistance Locus (marker-drug hospital vs community resistance gap)"):
            table_title_idx = idx
            in_summary = False
            break
        if not in_summary:
            continue
        patterns = [
            ("mean_overall_serious_r", r"Mean overall serious-R:\s*([0-9.+\-eE]+)%"),
            ("mean_hospital_serious_r", r"Mean hospital serious-R:\s*([0-9.+\-eE]+)%"),
            ("mean_community_serious_r", r"Mean community serious-R:\s*([0-9.+\-eE]+)%"),
        ]
        for key, pattern in patterns:
            match = re.search(pattern, stripped)
            if match:
                summary[key] = _figure_20_parse_number(match.group(1))

    if table_title_idx is None:
        return pd.DataFrame(), summary

    header_idx: int | None = None
    for idx in range(table_title_idx + 1, len(lines)):
        parts = _figure_20_split_row(lines[idx])
        if "Bacteria" in parts and "Marker drug(s)" in parts:
            header_idx = idx
            break
    if header_idx is None:
        return pd.DataFrame(), summary

    headers = _figure_20_split_row(lines[header_idx])
    rows: list[list[str]] = []
    for line in lines[header_idx + 1:]:
        stripped = line.strip()
        if not stripped:
            break
        parts = _figure_20_split_row(line)
        if len(parts) < 3:
            break
        while len(parts) < len(headers):
            parts.append("")
        rows.append(parts[:len(headers)])

    if not rows:
        return pd.DataFrame(columns=headers), summary

    df = pd.DataFrame(rows, columns=headers)
    required_columns = {
        "Bacteria",
        "Marker drug(s)",
        "Total New Infections",
        "Hospital Serious-R (%)",
        "Community Serious-R (%)",
        "Sim H:C ratio",
    }
    if not required_columns.issubset(df.columns):
        return pd.DataFrame(), summary
    for col in _F20_VALUE_COLUMNS:
        if col in df.columns:
            df[col] = df[col].map(_figure_20_parse_number)
    if "Overall Serious-R (%)" not in df.columns:
        df["Overall Serious-R (%)"] = np.nan
    if "Marker drug(s)" not in df.columns:
        df["Marker drug(s)"] = ""
    if "Total New Infections" in df.columns:
        df = df[df["Total New Infections"] > 0].copy()
    else:
        df = pd.DataFrame()
    if df.empty:
        return df, summary

    df["source_file"] = _resolve_project_path(path).name
    df["bacterium_key"] = df["Bacteria"].astype(str).str.strip().str.lower()
    df["Bacterium display"] = df["Bacteria"].map(_figure_20_bacterium_label)
    return df, summary


def _figure_20_summarise_rows(rows: pd.DataFrame) -> tuple[pd.DataFrame, list[str]]:
    notes: list[str] = []
    if rows.empty:
        return pd.DataFrame(), notes

    summary_rows: list[dict[str, object]] = []
    for key, group in rows.groupby("bacterium_key", sort=False):
        marker_values = [
            str(value).strip()
            for value in group["Marker drug(s)"].dropna().tolist()
            if str(value).strip()
        ]
        row: dict[str, object] = {
            "bacterium_key": key,
            "bacterium": group["Bacterium display"].iloc[0],
            "marker_drugs": marker_values[0] if marker_values else "\u2014",
            "n_runs": int(group["source_file"].nunique()),
            "missing_hospital_any_run": bool(group["Hospital Serious-R (%)"].isna().any()),
            "missing_community_any_run": bool(group["Community Serious-R (%)"].isna().any()),
        }

        for source_col, prefix in [
            ("Overall Serious-R (%)", "overall"),
            ("Hospital Serious-R (%)", "hospital"),
            ("Community Serious-R (%)", "community"),
            ("Sim H:C ratio", "sim_hc_ratio"),
            ("Total New Infections", "total_new_infections"),
        ]:
            values = group[source_col].dropna().astype(float).to_numpy()
            if len(values) == 0:
                row[f"{prefix}_median"] = np.nan
                row[f"{prefix}_p5"] = np.nan
                row[f"{prefix}_p95"] = np.nan
            else:
                row[f"{prefix}_median"] = float(np.nanmedian(values))
                row[f"{prefix}_p5"] = float(np.nanpercentile(values, 5))
                row[f"{prefix}_p95"] = float(np.nanpercentile(values, 95))
        summary_rows.append(row)

    summary = pd.DataFrame(summary_rows)
    summary["_overall_sort"] = summary["overall_median"].fillna(-1.0)
    summary = (
        summary.sort_values(
            ["_overall_sort", "total_new_infections_median"],
            ascending=[False, False],
        )
        .drop(columns=["_overall_sort"])
        .reset_index(drop=True)
    )
    return summary, notes


def _figure_20_format_percent(value: object) -> str:
    return f"{float(value):.1f}" if pd.notna(value) and np.isfinite(float(value)) else "\u2014"


def _figure_20_format_ratio(value: object) -> str:
    return f"{float(value):.2f}" if pd.notna(value) and np.isfinite(float(value)) else "\u2014"


def _figure_20_format_count(value: object) -> str:
    if value is None or pd.isna(value):
        return "\u2014"
    value_f = float(value)
    if not np.isfinite(value_f):
        return "\u2014"
    return f"{value_f:,.0f}"


def _figure_20_summary_box(summary_values: list[dict[str, float]], n_runs: int) -> str:
    if not summary_values:
        return ""
    parts: list[str] = []
    labels = [
        ("mean_overall_serious_r", "overall", ".1f", "%"),
        ("mean_hospital_serious_r", "hospital", ".1f", "%"),
        ("mean_community_serious_r", "community", ".1f", "%"),
    ]
    for key, label, fmt, suffix in labels:
        values = [
            float(item[key])
            for item in summary_values
            if key in item and np.isfinite(float(item[key]))
        ]
        if not values:
            continue
        value = float(np.nanmedian(values))
        parts.append(f"{label} {value:{fmt}}{suffix}")
    if not parts:
        return ""
    scope = "this run" if n_runs == 1 else f"{n_runs} runs (median summary values)"
    return (
        "<div class='meta-box'><strong>Mean serious-R in "
        + scope
        + ":</strong> "
        + ", ".join(parts)
        + ".</div>\n"
    )


def _figure_20_placeholder(out_dir: Path, agg: dict | None, message: str) -> None:
    fig, ax = plt.subplots(figsize=(10, 3.8))
    ax.text(
        0.5,
        0.5,
        f"{_F20_TITLE}\n\n{message}",
        ha="center",
        va="center",
        transform=ax.transAxes,
        fontsize=10.5,
        color="#555",
        bbox=dict(boxstyle="round,pad=0.6", fc="#f5f5f5", ec="#bbb"),
    )
    ax.set_axis_off()
    fig.subplots_adjust(left=0.03, right=0.97, top=0.92, bottom=0.08)
    _save_figure(
        fig,
        out_dir,
        _F20_STEM,
        _F20_TITLE,
        message,
        [
            "Serious-R is defined using the bacterium-specific marker drug(s) listed in "
            "the calibration summary. Hospital and community percentages are calculated "
            "among new active infections in the 2022-2025 calibration window. Missing "
            "hospital or community points indicate that the corresponding denominator was "
            "unavailable or zero in the summary.",
            "This figure uses serious-R marker resistance, not any-R. It should not be read "
            "as resistance to all drugs.",
        ],
        agg=agg,
    )


def make_figure_20_serious_r_by_hospital_community(
    calibration_paths: list[Union[str, Path]],
    out_dir: Path,
    agg: dict | None = None,
) -> dict[str, object]:
    parsed_frames: list[pd.DataFrame] = []
    parsed_summaries: list[dict[str, float]] = []
    for path in calibration_paths:
        frame, summary_values = _figure_20_parse_calibration_summary(path)
        if not frame.empty:
            parsed_frames.append(frame)
        if summary_values:
            parsed_summaries.append(summary_values)

    if not parsed_frames:
        _figure_20_placeholder(out_dir, agg, _F20_REQUIRED_MESSAGE)
        return {
            "generated": "placeholder",
            "bacteria_included": 0,
            "n_runs": 0,
            "missing_hospital_or_community_rows": 0,
            "summary_block_parsed": bool(parsed_summaries),
        }

    all_rows = pd.concat(parsed_frames, ignore_index=True)
    summary, _notes = _figure_20_summarise_rows(all_rows)
    if summary.empty:
        _figure_20_placeholder(out_dir, agg, _F20_REQUIRED_MESSAGE)
        return {
            "generated": "placeholder",
            "bacteria_included": 0,
            "n_runs": int(all_rows["source_file"].nunique()),
            "missing_hospital_or_community_rows": 0,
            "summary_block_parsed": bool(parsed_summaries),
        }

    n_runs = int(all_rows["source_file"].nunique())
    missing_rows = int(
        (
            summary["missing_hospital_any_run"]
            | summary["missing_community_any_run"]
            | summary["hospital_median"].isna()
            | summary["community_median"].isna()
        ).sum()
    )

    plot_summary = summary.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(plot_summary))
    fig_height = max(6.0, 2.2 + 0.34 * len(plot_summary))
    fig, ax = plt.subplots(figsize=(10.5, fig_height))

    for idx, row in plot_summary.iterrows():
        hospital = row["hospital_median"]
        community = row["community_median"]
        if pd.notna(hospital) and pd.notna(community):
            ax.plot(
                [float(community), float(hospital)],
                [idx, idx],
                color="#b0bec5",
                linewidth=0.9,
                alpha=0.65,
                zorder=1,
            )

    def _plot_metric(prefix: str, label: str, color: str, marker: str) -> None:
        values = plot_summary[f"{prefix}_median"].astype(float).to_numpy()
        mask = np.isfinite(values)
        if not mask.any():
            return
        xerr = None
        if n_runs > 1:
            p5 = plot_summary[f"{prefix}_p5"].astype(float).to_numpy()
            p95 = plot_summary[f"{prefix}_p95"].astype(float).to_numpy()
            err_low = np.where(np.isfinite(p5), np.clip(values - p5, 0.0, None), 0.0)
            err_high = np.where(np.isfinite(p95), np.clip(p95 - values, 0.0, None), 0.0)
            xerr = np.vstack([err_low[mask], err_high[mask]])
        ax.errorbar(
            values[mask],
            y[mask],
            xerr=xerr,
            fmt=marker,
            markersize=5.5,
            color=color,
            ecolor=color,
            elinewidth=0.8,
            capsize=2.5 if n_runs > 1 else 0,
            linestyle="none",
            label=label,
            zorder=3,
        )

    _plot_metric("community", "Community Serious-R (%)", "#2A9D8F", "o")
    _plot_metric("hospital", "Hospital Serious-R (%)", "#C65D00", "D")

    ax.set_yticks(y)
    ax.set_yticklabels(plot_summary["bacterium"].values, fontsize=7.5, fontstyle="italic")
    ax.set_xlim(0, 100)
    ax.set_xlabel("Serious-R among new active infections (%)", fontsize=10)
    ax.set_ylabel("Bacteria", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", linewidth=0.4, alpha=0.5)
    ax.legend(fontsize=8.5, frameon=False, loc="lower right")
    fig.suptitle(_F20_TITLE, fontsize=10.5, fontweight="bold")
    fig.tight_layout()

    details = pd.DataFrame({
        "Bacterium": summary["bacterium"],
        "Marker drug(s)": summary["marker_drugs"],
        "Total new infections": summary["total_new_infections_median"].map(_figure_20_format_count),
        "Overall serious-R (%)": summary["overall_median"].map(_figure_20_format_percent),
        "Hospital serious-R (%)": summary["hospital_median"].map(_figure_20_format_percent),
        "Community serious-R (%)": summary["community_median"].map(_figure_20_format_percent),
        "Sim H:C ratio": summary["sim_hc_ratio_median"].map(_figure_20_format_ratio),
    })

    summary_box = _figure_20_summary_box(parsed_summaries, n_runs)
    s6_link = (
        "<p class='note'>For the raw new active infection denominator counts used to "
        "interpret bacterium-specific serious-R estimates, see "
        "<a href='Supplementary_Figure_S6__new_active_infection_denominators_by_bacterium.html'>"
        "Supplementary Figure S6</a>.</p>\n"
    )
    extra_html = summary_box + s6_link + "<h2>Figure 7 Details</h2>\n" + _html_table(details)

    run_note = (
        f"Values are medians across {n_runs} calibration summary file"
        f"{'s' if n_runs > 1 else ''}; horizontal intervals show 5th-95th percentile ranges. "
        "Total new infections in the table is the median across runs."
        if n_runs > 1
        else "Values are taken directly from the supplied calibration summary file."
    )
    footnotes = [
        "Serious-R is defined using the bacterium-specific marker drug(s) listed in the "
        "calibration summary. Hospital and community percentages are calculated among new "
        "active infections in the 2022-2025 calibration window. Missing hospital or community "
        "points indicate that the corresponding denominator was unavailable or zero in the summary.",
        "This figure uses serious-R marker resistance, not any-R. It should not be read as "
        "resistance to all drugs.",
        run_note,
    ]
    _save_figure(
        fig,
        out_dir,
        _F20_STEM,
        _F20_TITLE,
        "Horizontal paired-dot plot of hospital-associated and community-associated serious-R "
        "percentages from the Serious Resistance Locus calibration-summary table.",
        footnotes,
        agg=agg,
        extra_html=extra_html,
    )

    print(
        "  Figure 7: "
        f"{len(summary)} bacteria included from {n_runs} calibration summary file"
        f"{'s' if n_runs > 1 else ''}; "
        f"{missing_rows} row{'s' if missing_rows != 1 else ''} had missing hospital or community values; "
        f"summary block {'parsed' if parsed_summaries else 'not found'}."
    )
    return {
        "generated": "real data",
        "bacteria_included": len(summary),
        "n_runs": n_runs,
        "missing_hospital_or_community_rows": missing_rows,
        "summary_block_parsed": bool(parsed_summaries),
    }


# Legacy supplementary/additional figures below are retained for reference only;
# they are not called by main().
def make_fs2_syndrome_bars(agg: dict, out_dir: Path) -> None:
    """
    Figure FS2: horizontal bar chart of syndrome incidence per 100,000
    (figure version of Supplementary Table S2).
    """
    si = agg.get("syndrome_incidence", pd.DataFrame()).copy()
    n  = agg.get("n_runs", 1)
    if si is None or si.empty:
        print("  FS2: no syndrome_incidence data — skipping.")
        return
    first_col = si.columns[0]
    si = si[si[first_col].astype(str).str.upper() != "TOTAL"].copy()
    inc_col = next((c for c in si.columns if "incidence" in c.lower()), None)
    if inc_col is None:
        return
    si = _add_interval_columns(si, inc_col, "_inc")
    si[inc_col] = si["_inc_med"]
    si = si.dropna(subset=[inc_col]).sort_values(inc_col, ascending=True)
    share_col = next((c for c in si.columns if "share" in c.lower()), None)
    if share_col:
        si = _add_interval_columns(si, share_col, "_share")
        si[share_col] = si["_share_med"]
    tab_colors = plt.cm.tab10(np.linspace(0, 0.9, len(si)))
    fig, ax = plt.subplots(figsize=(8, max(3, 0.55 * len(si))))
    inc_err_lo, inc_err_hi = _asymmetric_errors(si, "_inc")
    ax.barh(
        range(len(si)),
        si[inc_col].values,
        color=tab_colors,
        edgecolor="none",
        height=0.6,
        xerr=[inc_err_lo, inc_err_hi] if n > 1 else None,
        error_kw={"elinewidth": 0.9, "ecolor": "#333", "capthick": 0.9, "capsize": 2.5},
    )
    ax.set_yticks(range(len(si)))
    ax.set_yticklabels(si[first_col].values, fontsize=9)
    ax.set_xlabel("Incidence per 100,000 population per year", fontsize=10)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="x", linewidth=0.4, alpha=0.5)
    for i, v in enumerate(si[inc_col].values):
        ax.text(v + max(si[inc_col].max() * 0.01, 5),
                i, f"{v:.0f}", va="center", ha="left", fontsize=8)
    fig.suptitle("Figure FS2 \u2014 Syndrome incidence: simulated annual rates",
                 fontsize=10, fontweight="bold")
    fig.tight_layout()
    _save_figure(
        fig, out_dir, "FS2_syndrome_bars",
        "Figure FS2 \u2014 Syndrome Incidence: Simulated Annual Rates",
        f"Figure version of Supplementary Table S2. "
        f"{'Error bars show 5th-95th percentile ranges across accepted runs. ' if n > 1 else ''}"
        f"n\u2009=\u2009{n} run{'s' if n > 1 else ''}.",
        ["Syndrome incidence is the simulated annual rate per 100,000 population during "
         "the 2022\u20132025 calibration window.",
         "No external calibration targets are defined for individual syndromes; the "
         "distribution is an emergent product of organism-specific infection rates and "
         "syndrome-assignment probabilities.",
         "Figure version of Supplementary Table S2."],
        subfolder="supplementary",
        agg=agg,
    )


# ---------------------------------------------------------------------------
# Supplementary Figure FS3 — Resistance by acquisition route (figure version of S3)
# ---------------------------------------------------------------------------

def make_fs3_resistance_scatter(agg: dict, out_dir: Path) -> None:
    """
    Figure FS3: scatter of hospital any-R% vs. community any-R% per organism
    (figure version of Supplementary Table S3).
    """
    ril = agg.get("resistance_incidence_locus", pd.DataFrame())
    n   = agg.get("n_runs", 1)
    if ril is None or ril.empty:
        print("  FS3: no resistance_incidence_locus data — skipping.")
        return
    first_col = ril.columns[0]
    summary_mask = ril[first_col].astype(str).str.match(
        r"^\s*(-|Resistance Locus|Serious Resistance|Mean |H:C)", na=False)
    ril = ril[~summary_mask].copy()
    hosp_col = next((c for c in ril.columns
                     if "hospital" in c.lower() and "%" in c
                     and "total" not in c.lower()), None)
    comm_col = next((c for c in ril.columns
                     if "community" in c.lower() and "%" in c
                     and "total" not in c.lower()), None)
    if hosp_col is None or comm_col is None:
        print("  FS3: cannot find hospital/community resistance columns — skipping.")
        return
    ril = _add_interval_columns(ril, hosp_col, "_hosp")
    ril = _add_interval_columns(ril, comm_col, "_comm")
    ril[hosp_col] = ril["_hosp_med"]
    ril[comm_col] = ril["_comm_med"]
    ril = ril.dropna(subset=[hosp_col, comm_col])
    if ril.empty:
        print("  FS3: no valid resistance rows after filtering — skipping.")
        return
    raw_max = max(ril[hosp_col].max(), ril[comm_col].max())
    if not np.isfinite(raw_max) or raw_max <= 0:
        print("  FS3: resistance values are NaN/Inf/zero — skipping.")
        return
    max_val = raw_max * 1.12
    fig, ax = plt.subplots(figsize=(8, 7))
    ax.fill_between([0, max_val], [0, max_val], [max_val, max_val],
                    alpha=0.05, color="#5C6BC0")
    comm_err_lo, comm_err_hi = _asymmetric_errors(ril, "_comm")
    hosp_err_lo, hosp_err_hi = _asymmetric_errors(ril, "_hosp")
    if n > 1:
        ax.errorbar(
            ril[comm_col],
            ril[hosp_col],
            xerr=[comm_err_lo, comm_err_hi],
            yerr=[hosp_err_lo, hosp_err_hi],
            fmt="none",
            ecolor="#455A64",
            elinewidth=0.7,
            capsize=2.0,
            alpha=0.55,
            zorder=2,
        )
    ax.scatter(ril[comm_col], ril[hosp_col], s=52,
               color="#5C6BC0", edgecolors="white", linewidths=0.4, zorder=3, alpha=0.85)
    ax.plot([0, max_val], [0, max_val], color="#555", linewidth=1.0, linestyle="--",
            label="Hospital = community", zorder=2)
    for _, row in ril.iterrows():
        parts = str(row[first_col]).split()
        abbr  = (parts[0][0] + ". " + " ".join(parts[1:])) if len(parts) > 1 else str(row[first_col])
        ax.annotate(abbr, (row[comm_col], row[hosp_col]),
                    fontsize=6, xytext=(3, 2), textcoords="offset points",
                    color="#333", clip_on=True)
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ax.set_xlabel("Community-acquired new infections with any resistance (%)", fontsize=10)
    ax.set_ylabel("Hospital-acquired new infections with any resistance (%)", fontsize=10)
    ax.legend(fontsize=9, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(linewidth=0.35, alpha=0.5)
    ax.text(max_val * 0.55, max_val * 0.88,
            "Higher hospital\nresistance \u2191", fontsize=8.5,
            color="#5C6BC0", alpha=0.7, ha="center")
    fig.suptitle(
        "Figure FS3 \u2014 Resistance by acquisition route: hospital vs. community",
        fontsize=10, fontweight="bold")
    fig.tight_layout()
    _save_figure(
        fig, out_dir, "FS3_resistance_scatter",
        "Figure FS3 \u2014 Resistance by Acquisition Route: Hospital vs. Community",
        f"Each point is one organism. Points above the dashed diagonal carry higher resistance "
        f"in hospital-acquired than community-acquired infections (blue shaded region). "
        f"{'Error bars show 5th-95th percentile ranges across accepted runs. ' if n > 1 else ''}"
        f"n\u2009=\u2009{n} run{'s' if n > 1 else ''}.",
        ["Each axis shows the percentage of new infections of that organism carrying any "
         "resistance mechanism, averaged across all drugs with non-negligible potency for "
         "that organism.",
         "Figure version of Supplementary Table S3."],
        subfolder="supplementary",
        agg=agg,
    )


# ---------------------------------------------------------------------------
# Additional Figure FA1 — Infection vs. carriage ecology (42 organisms)
# ---------------------------------------------------------------------------

def make_fa1_infection_carriage(agg: dict, out_dir: Path) -> None:
    """
    Figure FA1 (additional): scatter/bubble of simulated infection prevalence vs.
    carriage prevalence for all 42 organisms, bubble-sized by infection deaths.
    Shows the pathogen ecology spectrum from obligate invaders to pure colonisers.
    """
    bi  = agg.get("bacteria_infections", pd.DataFrame()).copy()
    bm  = agg.get("bacteria_mortality",  pd.DataFrame()).copy()
    n   = agg.get("n_runs", 1)
    if bi is None or bi.empty:
        return
    bact_col = bi.columns[0]
    inf_col  = next((c for c in bi.columns if "infection simulation" in c.lower()), None)
    carr_col = next((c for c in bi.columns if "carriage simulation" in c.lower()), None)
    if inf_col is None or carr_col is None:
        print("  FA1: infection/carriage simulation columns not found — skipping.")
        return
    bi = _add_interval_columns(bi, inf_col, "_inf")
    bi = _add_interval_columns(bi, carr_col, "_carr")
    bi[inf_col] = bi["_inf_med"]
    bi[carr_col] = bi["_carr_med"]
    df = bi[[bact_col, inf_col, carr_col]].dropna(subset=[inf_col, carr_col]).copy()
    if df.empty:
        print("  FA1: no valid infection/carriage rows after filtering — skipping.")
        return
    if not bm.empty:
        death_col = next((c for c in bm.columns
                          if "simulation" in c.lower() and "death" in c.lower()), None)
        if death_col:
            bm = _add_interval_columns(bm, death_col, "_death")
            bm[death_col] = bm["_death_med"]
            df = df.merge(bm[[bm.columns[0], death_col]],
                          left_on=bact_col, right_on=bm.columns[0], how="left")
            df["_deaths"] = df[death_col].fillna(0)
        else:
            df["_deaths"] = 0.0
    else:
        df["_deaths"] = 0.0
    max_d = df["_deaths"].max()
    df["_size"] = 20 + (df["_deaths"] / max_d * 600 if max_d > 0 else 0)
    # Colour by infection:carriage ratio (capped to finite range)
    df["_ratio"] = np.where(df[carr_col] > 0, df[inf_col] / df[carr_col], 100.0)
    log_ratio = np.clip(np.log10(df["_ratio"].clip(lower=1e-4)), -2, 2)
    max_v = max(df[inf_col].max(), df[carr_col].max()) * 1.2
    if not np.isfinite(max_v) or max_v <= 0:
        print("  FA1: infection/carriage values are NaN/Inf/zero — skipping.")
        return
    fig, ax = plt.subplots(figsize=(9.5, 8))
    sc = ax.scatter(
        df[carr_col], df[inf_col],
        s=df["_size"], c=log_ratio,
        cmap="RdYlGn_r", vmin=-2, vmax=2,
        alpha=0.82, edgecolors="white", linewidths=0.5, zorder=3,
    )
    cbar = fig.colorbar(sc, ax=ax, fraction=0.030, pad=0.02)
    cbar.set_label("log\u2081\u2080(infection / carriage ratio)", fontsize=8)
    cbar.set_ticks([-2, -1, 0, 1, 2])
    cbar.set_ticklabels(["0.01\xd7", "0.1\xd7", "1\xd7", "10\xd7", "100\xd7"])
    ax.plot([0, max_v], [0, max_v], color="#888", linewidth=0.8, linestyle="--",
            label="Infection = carriage prevalence", zorder=2)
    for _, row in df.iterrows():
        parts = str(row[bact_col]).split()
        abbr  = (parts[0][0] + ". " + " ".join(parts[1:])) if len(parts) > 1 else str(row[bact_col])
        ax.annotate(abbr, (row[carr_col], row[inf_col]),
                    fontsize=5.5, xytext=(3, 2), textcoords="offset points",
                    color="#333", clip_on=True)
    ax.set_xlim(0, max_v)
    ax.set_ylim(0, max_v)
    ax.set_xscale("symlog", linthresh=0.005)
    ax.set_yscale("symlog", linthresh=0.005)
    ax.set_xlabel("Carriage prevalence — simulation (% world population)", fontsize=10)
    ax.set_ylabel("Infection prevalence — simulation (% world population)", fontsize=10)
    ax.legend(fontsize=9, frameon=False)
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(linewidth=0.35, alpha=0.4)
    # Size legend
    if max_d > 0:
        for deaths_val, label in [(0.1, "0.1 M"), (0.5, "0.5 M"), (2.0, "2 M")]:
            size = 20 + deaths_val / max_d * 600
            ax.scatter([], [], s=size, color="#888", alpha=0.5, label=f"Deaths: {label}/yr")
        ax.legend(fontsize=7.5, frameon=False, loc="upper left")
    fig.suptitle(
        "Figure FA1 \u2014 Infection vs. carriage ecology across 42 organisms",
        fontsize=10, fontweight="bold")
    fig.tight_layout()
    _save_figure(
        fig, out_dir, "FA1_infection_carriage",
        "Figure FA1 \u2014 Infection vs. Carriage Ecology Across 42 Organisms",
        "Each point is one organism; bubble area is proportional to simulated annual "
        "infection deaths. Colour: red = infection >> carriage (obligate/invasive pathogens); "
        "green = carriage >> infection (mainly commensals). Both axes are symmetric log scale.",
        ["Points above the dashed diagonal have higher infection than carriage prevalence, "
         "characteristic of invasive obligate pathogens (e.g. Neisseria gonorrhoeae, "
         "Neisseria meningitidis, Mycobacterium tuberculosis).",
         "Points below the diagonal have higher carriage than infection prevalence, "
         "characteristic of opportunistic commensals (e.g. E. coli, S. epidermidis, "
         "Bacteroides fragilis).",
         "Bubble area is proportional to simulated infection deaths. Organisms with zero "
         "deaths are shown at a fixed minimum size.",
         "This figure is not derived from an existing table; it provides a cross-organism "
         "ecological overview of the simulated pathogen landscape."],
        agg=agg,
    )


# ---------------------------------------------------------------------------
# Additional Figure FA2 — Age distribution of infections and deaths
# ---------------------------------------------------------------------------

def make_fa2_age_distribution(agg: dict, out_dir: Path) -> None:
    """
    Figure FA2 (additional): stacked horizontal bars showing the simulated
    age distribution (<5 yr, 5–64 yr, ≥65 yr) of infections and deaths
    for all 42 organisms.
    """
    bi  = agg.get("bacteria_infections", pd.DataFrame()).copy()
    bm  = agg.get("bacteria_mortality",  pd.DataFrame()).copy()
    n   = agg.get("n_runs", 1)

    bact_inf  = bi.columns[0]  if not bi.empty else None
    u5_inf    = next((c for c in bi.columns if "<5"   in c and "infection" in c.lower()), None) \
                if not bi.empty else None
    o65_inf   = next((c for c in bi.columns if "65+"  in c and "infection" in c.lower()
                      or "65" in c and "infection" in c.lower()), None) \
                if not bi.empty else None
    bact_mort = bm.columns[0]  if not bm.empty else None
    u5_mort   = next((c for c in bm.columns if "<5"   in c and "mortalit" in c.lower()), None) \
                if not bm.empty else None
    o65_mort  = next((c for c in bm.columns if "65+"  in c and "mortalit" in c.lower()
                      or "65" in c and "mortalit" in c.lower()), None) \
                if not bm.empty else None

    if u5_inf is None and u5_mort is None:
        print("  FA2: no age-distribution columns found — skipping.")
        return

    def _plot_age(ax, df, bc, u5c, o65c, title, xlabel):
        if df is None or df.empty or u5c is None or o65c is None:
            ax.text(0.5, 0.5, "No data", ha="center", va="center",
                    transform=ax.transAxes, fontsize=11, color="#888")
            ax.set_title(title, fontsize=10); return
        df = df.copy()
        df[u5c]  = pd.to_numeric(df[u5c],  errors="coerce").fillna(0)
        df[o65c] = pd.to_numeric(df[o65c], errors="coerce").fillna(0)
        df["_mid"] = (100 - df[u5c] - df[o65c]).clip(lower=0)
        df = df.sort_values(u5c, ascending=True)
        y = np.arange(len(df))
        ax.barh(y, df[u5c].values,   0.6, color="#42A5F5", label="<5 years")
        ax.barh(y, df["_mid"].values, 0.6, left=df[u5c].values,
                color="#66BB6A", label="5\u201364 years")
        ax.barh(y, df[o65c].values,  0.6,
                left=(df[u5c] + df["_mid"]).values,
                color="#FF7043", label="\u226565 years")
        ax.set_yticks(y)
        ax.set_yticklabels(df[bc].values, fontsize=6.5, fontstyle="italic")
        ax.set_xlabel(xlabel, fontsize=9)
        ax.axvline(100, color="#aaa", linewidth=0.5, linestyle="--")
        ax.set_xlim(0, 105)
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_title(title, fontsize=10, fontweight="bold")
        ax.legend(fontsize=8, frameon=False, loc="lower right")

    n_orgs = max(len(bi) if not bi.empty else 0,
                 len(bm) if not bm.empty else 0)
    fig, axes = plt.subplots(1, 2, figsize=(14, max(5, 0.33 * n_orgs)))
    _plot_age(axes[0], bi, bact_inf,  u5_inf,  o65_inf,
              "(A) Age distribution of infections",
              "Share of infections in age group (%)")
    _plot_age(axes[1], bm, bact_mort, u5_mort, o65_mort,
              "(B) Age distribution of deaths",
              "Share of deaths in age group (%)")
    fig.suptitle(
        "Figure FA2 \u2014 Age distribution of bacterial infections and deaths by organism",
        fontsize=10, fontweight="bold")
    fig.tight_layout()
    _save_figure(
        fig, out_dir, "FA2_age_distribution",
        "Figure FA2 \u2014 Age Distribution of Bacterial Infections and Deaths",
        f"Stacked bars showing the simulated proportion of infections (left panel) and deaths "
        f"(right panel) in the <5, 5\u201364, and \u226565 age groups, for each of 42 organisms. "
        f"Sorted by percentage in the under-5 age group. n\u2009=\u2009{n} run{'s' if n > 1 else ''}.",
        ["Values are percentages of that organism\u2019s total infections or deaths "
         "attributable to each age group, during the 2022\u20132025 calibration window.",
         "The \u20185\u201364 years\u2019 bar is derived as 100% minus the <5 and \u226565 shares.",
         "This figure is not derived from an existing table; it provides an age-structured "
         "overview of infection and mortality burden across all modelled organisms."],
        agg=agg,
    )


# ---------------------------------------------------------------------------
# Index page
# ---------------------------------------------------------------------------

def make_legacy_index(agg: dict, out_dir: Path) -> None:
    m  = agg.get("meta", {})
    n  = agg.get("n_runs", 1)
    cs = agg.get("calibration_score", {})

    body  = _html_head("Paper Tables — AMR Simulation")
    body += "<h1>Paper Tables — AMR Simulation Calibration</h1>\n"
    body += f"<p class='note'>Accepted runs: {n}</p>\n"

    body += "<h2>Main manuscript tables</h2>\n<ul>\n"
    for fname, label, desc in [
        ("main/T1_model_summary.html",
         "Table 1 \u2014 Model summary",
         "Individual-based model design: scope, biological detail, calibration approach, and counterfactual"),
        ("main/T2_headline_metrics.html",
         "Table 2 — Headline calibration metrics",
         "4 headline metrics + calibration block scores"),
        ("main/T3_drug_class_share.html",
         "Table 3 — Drug class share",
         "28 antibiotic classes: simulation vs. global use targets"),
        ("main/T4_bacteria_burden.html",
         "Table 4 — Bacterial infection &amp; carriage prevalence",
         "42 organisms: infection%, carriage%"),
        ("main/T5_resistance_fit.html",
         "Table 5 — Hospital vs. community resistance by organism",
         "Organisms with any-R structural H:C benchmark &gt; 1: hospital-acquired%, any-R%, serious-R% by locus"),
        ("main/T6_amr_attributable_deaths_PLACEHOLDER.html",
         "Table 6 — AMR-attributable deaths",
         "Placeholder — requires completed counterfactual runs"),
    ]:
        body += f"<li><a href='{fname}'><strong>{label}</strong></a> — {desc}</li>\n"
    body += "</ul>\n"

    body += "<h2>Main manuscript figures</h2>\n<ul>\n"
    for fname, label, desc in [
        ("main/F1_resistance_trend.html",
         "Figure F1 \u2014 Proportion of new infections with any resistance, 1930\u20132025",
         "Historical resistance trend: mean + 90% interval across accepted runs"),
        ("main/F2_resistance_fit.html",
         "Figure F2 \u2014 Resistance calibration fit by organism",
         "Multi-panel bar chart: simulated vs. calibration-benchmark infection resistance "
         "by drug class for IHME-priority pathogens"),
        ("main/F3_calibration_scores.html",
         "Figure F3 \u2014 Calibration block scores",
         "Horizontal bar chart: normalised block scores vs. acceptance threshold (1.0)"),
        ("main/F4_headline_metrics.html",
         "Figure F4 \u2014 Headline calibration metrics: simulation vs. target",
         "4-panel grouped bars (simulation vs. target) for headline calibration metrics; figure version of T2"),
        ("main/F5_drug_class_share.html",
         "Figure F5 \u2014 Antibiotic use by drug class: simulation vs. global estimates",
         "Horizontal paired bars coloured by WHO AWaRe category; figure version of T3"),
        ("main/F6_bacteria_scatter.html",
         "Figure F6 \u2014 Bacterial infection prevalence: simulation vs. calibration target",
         "Scatter plot of simulated vs. target infection % for all 42 organisms; figure version of T4"),
        ("main/F7_hc_resistance_heatmap.html",
         "Figure F7 \u2014 Hospital vs. community resistance and acquisition rates",
         "Heatmap of hospital/community resistance and acquisition rates; figure version of T5"),
    ]:
        body += f"<li><a href='{fname}'><strong>{label}</strong></a> \u2014 {desc}</li>\n"
    body += "</ul>\n"

    body += "<h2>Supplementary tables</h2>\n<ul>\n"
    for fname, label, desc in [
        ("supplementary/S1_infection_deaths.html",
         "S1 — Infection deaths per organism",
         "42 organisms: deaths target and simulation (millions/year)"),
        ("supplementary/S2_syndrome_incidence.html",
         "S2 — Syndrome incidence",
         "Simulated annual incidence per 100,000 population by syndrome"),
        ("supplementary/S3_resistance_by_acquisition_route.html",
         "S3 — Percentage of new infections with resistance by organism and hospital status",
         "Hospital-acquired vs. community-acquired new infections carrying any resistance"),
        ("supplementary/S4_resistance_benchmarks.html",
         "S4 — Full resistance benchmarks per organism",
         "All non-negligible organism × drug combinations: infection and carriage resistance"),
    ]:
        body += f"<li><a href='{fname}'><strong>{label}</strong></a> \u2014 {desc}</li>\n"
    body += "</ul>\n"

    body += "<h2>Supplementary figures</h2>\n<ul>\n"
    for fname, label, desc in [
        ("supplementary/FS1_mortality_bars.html",
         "Figure FS1 \u2014 Infection deaths per organism: simulation vs. observed estimate",
         "Horizontal paired bars sorted by burden; figure version of S1"),
        ("supplementary/FS2_syndrome_bars.html",
         "Figure FS2 \u2014 Syndrome incidence: simulated annual rates",
         "Horizontal bar chart of syndrome incidence per 100,000; figure version of S2"),
        ("supplementary/FS3_resistance_scatter.html",
         "Figure FS3 \u2014 Resistance by acquisition route: hospital vs. community",
         "Scatter of hospital any-R% vs. community any-R% per organism; figure version of S3"),
    ]:
        body += f"<li><a href='{fname}'><strong>{label}</strong></a> \u2014 {desc}</li>\n"
    body += "</ul>\n"

    body += "<h2>Additional figures</h2>\n<ul>\n"
    for fname, label, desc in [
        ("main/FA1_infection_carriage.html",
         "Figure FA1 \u2014 Infection vs. carriage ecology across 42 organisms",
         "Bubble chart: infection% vs. carriage%, sized by deaths, coloured by infection:carriage ratio"),
        ("main/FA2_age_distribution.html",
         "Figure FA2 \u2014 Age distribution of bacterial infections and deaths",
         "Stacked bars: <5, 5\u201364, \u226565 age-group shares of infections and deaths for all 42 organisms"),
    ]:
        body += f"<li><a href='{fname}'><strong>{label}</strong></a> \u2014 {desc}</li>\n"
    body += "</ul>\n"
    body += "</body></html>"
    _save(out_dir / "index.html", body)


# ---------------------------------------------------------------------------
# Focused paper-output index
# ---------------------------------------------------------------------------

def _index_existing_link_items(
    out_dir: Path,
    items: list[tuple[str, str]],
    empty_message: str,
) -> str:
    body = ""
    for fname, label in items:
        if (out_dir / fname).exists():
            body += f"<li><a href='{fname}'><strong>{label}</strong></a></li>\n"
        else:
            print(f"  Index: skipped missing output {fname}")
    if not body:
        body = f"<li><em>{empty_message}</em></li>\n"
    return body


def make_index(agg: dict, out_dir: Path) -> None:
    n = agg.get("n_runs", 1)

    body  = _html_head("Paper Outputs — AMR Simulation")
    body += "<h1>Paper outputs — AMR Simulation Calibration</h1>\n"
    body += f"<p class='note'>Accepted runs: {n}</p>\n"

    body += "<h2>Tables</h2>\n<ul>\n"
    body += (
        "<li><a href='Tables/T1__model_summary.html'>"
        "<strong>Table 1. Model summary</strong></a></li>\n"
    )
    body += "</ul>\n"

    body += "<h2>Supplementary Tables</h2>\n<ul>\n"
    body += _index_existing_link_items(
        out_dir,
        [
            (
                "Tables/Supplementary_Table_S2__detailed_bacterium_drug_resistance_benchmarks.html",
                "Supplementary Table S2. Detailed bacterium-drug resistance benchmarks, "
                "2022\u20132025",
            ),
        ],
        "No supplementary tables were generated.",
    )
    body += "</ul>\n"

    body += "<h2>Figures</h2>\n<ul>\n"
    body += _index_existing_link_items(
        out_dir,
        [
            (
                "Figures/Figure_1__calibration_headline_metrics.html",
                "Figure 1. Calibration: 2025 headline health and antibiotic-use metrics",
            ),
            (
                "Figures/Figure_2__calibration_resistance_fit_by_bacteria_drug_class.html",
                "Figure 2. Calibration: resistance fit by bacterium and drug class",
            ),
            (
                "Figures/Figure_3__calibration_drug_class_share.html",
                "Figure 3. Calibration: 2025 antibiotic use by drug class",
            ),
            (
                "Figures/Figure_4__calibration_infection_deaths_by_bacteria.html",
                "Figure 4. Calibration: 2025 infection deaths by bacterium",
            ),
            (
                "Figures/Figure_5__calibration_carriage_prevalence_by_bacteria.html",
                "Figure 5. Calibration: 2025 prevalence of carriage by bacterium",
            ),
            ("Figures/Figure_6A__resistance_trends.html", "Figure 6A. Resistance trends"),
            (
                "Figures/Figure_6B__resistance_trends_by_bacterium.html",
                "Figure 6B. Resistance trends by bacterium",
            ),
            (
                "Figures/Figure_6C__serious_r_trends_by_bacterium.html",
                "Figure 6C. Serious-R trends by bacterium",
            ),
            (
                "Figures/Figure_7__serious_r_by_hospital_community.html",
                "Figure 7. Serious-R by hospital and community, 2022\u20132025",
            ),
            (
                "Figures/Figure_8__infection_death_rate_by_region.html",
                "Figure 8. Infection death rate by region, 2022\u20132025",
            ),
            (
                "Figures/Figure_9__antibiotic_use_by_treatment_context.html",
                "Figure 9. Antibiotic use by treatment context, 2025",
            ),
            (
                "Figures/Figure_10__sepsis_context_effective_therapy.html",
                "Figure 10. Underlying sepsis onset context and time to effective therapy, "
                "2022\u20132025",
            ),
            (
                "Figures/Figure_11__activity_retained_by_bacterium.html",
                "Figure 11. Resistance-adjusted antibiotic activity retained by bacterium, 2022\u20132025",
            ),
            (
                "Figures/Figure_12__distribution_drug_use_by_bacteria.html",
                "Figure 12. Antibiotic exposure distribution by bacterium, 2022\u20132025",
            ),
            (
                "Figures/Figure_13__resistance_pathway_counterfactuals.html",
                "Figure 13. Counterfactual resistance-acquisition pathway comparisons",
            ),
        ],
        "No main figures were generated.",
    )
    body += "</ul>\n"

    body += "<h2>Supplementary Figures</h2>\n<ul>\n"
    body += _index_existing_link_items(
        out_dir,
        [
            (
                "Figures/Supplementary_Figure_S1__potential_activity_retained.html",
                "Supplementary Figure S1. Potential antibiotic activity retained across "
                "available drugs",
            ),
            (
                "Figures/Supplementary_Figure_S2__microbiome_resistance_reservoir.html",
                "Supplementary Figure S2. Microbiome resistance reservoir, 2022\u20132025",
            ),
            (
                "Figures/Supplementary_Figure_S3__carrier_vs_non_carrier_infection_incidence.html",
                "Supplementary Figure S3. Carrier versus non-carrier infection incidence, "
                "2022\u20132025",
            ),
            (
                "Figures/Supplementary_Figure_S5__diagnostic_testing_targeted_treatment_cascade.html",
                "Supplementary Figure S5. Diagnostic testing and targeted-treatment cascade, "
                "2022\u20132025",
            ),
            (
                "Figures/Supplementary_Figure_S6__new_active_infection_denominators_by_bacterium.html",
                "Supplementary Figure S6. New active infection denominators by bacterium, "
                "2022\u20132025",
            ),
            (
                "Figures/Supplementary_Figure_S7__active_infection_incidence_by_bacterium.html",
                "Supplementary Figure S7. Active infection incidence by bacterium: "
                "simulation versus target, 2022\u20132025",
            ),
            (
                "Figures/Supplementary_Figure_S8__infection_outcome_pathway_by_bacterium.html",
                "Supplementary Figure S8. Infection outcome pathway by bacterium, "
                "2022\u20132025",
            ),
        ],
        "No supplementary figures were generated.",
    )
    body += "</ul>\n"
    body += "<h2>Supplementary / Diagnostic Figures</h2>\n<ul>\n"
    body += _index_existing_link_items(
        out_dir,
        [
            (
                "Figures/Supplementary_Figure_SX__modelled_resistance_mechanisms_by_bacterium.html",
                "Supplementary Figure SX. Modelled resistance mechanisms by bacterium, "
                "2022\u20132025",
            ),
        ],
        "No supplementary diagnostic figures were generated.",
    )
    body += "</ul>\n"
    body += "</body></html>"
    _save(out_dir / "index.html", body)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(input_args: list[str]) -> None:
    if not input_args:
        print(
            "Usage: python -m amr_simulation_output_analysis.make_paper_tables "
            "<calibration_summary_*.txt> [...]"
        )
        sys.exit(1)

    # Expand any glob patterns
    paths: list[Path] = []
    seen_paths: set[Path] = set()
    for arg in input_args:
        candidates = [_resolve_project_path(item) for item in glob.glob(arg)]
        if not candidates and not Path(arg).is_absolute():
            candidates = [Path(item) for item in glob.glob(str(REPO_ROOT / arg))]
        if not candidates:
            candidate = _resolve_project_path(arg)
            if candidate.exists():
                candidates.append(candidate)

        if not candidates:
            print(f"  Warning: no file found for '{arg}'")
            continue

        for candidate in candidates:
            resolved = candidate.resolve()
            if resolved in seen_paths:
                continue
            seen_paths.add(resolved)
            paths.append(candidate)

    if not paths:
        print("No input files found.")
        sys.exit(1)

    print(f"Parsing {len(paths)} calibration file(s)...")
    runs = parse_files(paths)
    agg  = aggregate(runs)
    n    = len(runs)
    print(f"  -> {n} run(s) parsed and aggregated.")

    # Auto-discover matching simulation_summary CSVs for simulation-output figures.
    csv_paths = _discover_f1_simulation_csvs(paths)
    csv_runs_with_scale = _discover_simulation_csvs_with_scale(paths)
    if csv_paths:
        print(
            f"  Found {len(csv_paths)} simulation CSV(s) for Figures 6A, 6B, 8, 9, 10, 11, 12, "
            "Supplementary Figures S1, S3, S5, S6, and S8, "
            "and diagnostic Supplementary Figure SX. "
            "Supplementary Table S2 and Supplementary Figures S2 and S7 use calibration summary tables."
        )
    else:
        print(
            "  No matching simulation CSVs found; Figures 6A, 6B, 8, 9, 10, 11, 12, and "
            "Supplementary Figures S1, S3, S5, S6, S8, and diagnostic SX may render as placeholders. "
            "Supplementary Table S2 and Supplementary Figures S2 and S7 use calibration summary tables."
        )

    st1_csv_paths = _filter_simulation_csvs_with_columns(
        csv_paths,
        _ST1_REQUIRED_VECTOR_COLUMNS,
        "Supplementary Figure S8",
    )
    sf1_csv_paths = _filter_simulation_csvs_with_columns(
        csv_paths,
        [_SF1_NUMERATOR_COLUMN, _SF1_DENOMINATOR_COLUMN],
        "Supplementary Figure S1",
    )
    sf3_csv_paths = _filter_simulation_csvs_with_columns(
        csv_paths,
        _SF3_REQUIRED_COLUMNS,
        "Supplementary Figure S3",
    )
    sf4_csv_paths = csv_paths
    sf5_required_columns = [
        column
        for _, suffix in _SF5_SETTINGS
        for column in _sf5_stage_columns_for_suffix(suffix)
    ]
    sf5_csv_paths = _filter_simulation_csvs_with_columns(
        csv_paths,
        ["time_in_years", "policy_option", *sf5_required_columns],
        "Supplementary Figure S5",
    )
    sf6_csv_paths = _filter_simulation_csvs_with_columns(
        csv_paths,
        [_SF6_TOTAL_COLUMN],
        "Supplementary Figure S6",
    )

    out = OUT_DIR
    print(f"\nGenerating paper outputs in {out.absolute()} ...")
    _prepare_output_dirs(out)
    make_t1(out)
    make_supplementary_table_s2_resistance_benchmarks(runs, out, agg=agg)
    make_figure_1_calibration_headline_metrics(agg, out)
    make_figure_2_calibration_resistance_fit(
        agg,
        out,
        runs=runs,
        summary_mode=FIGURE2_SUMMARY_MODE,
    )
    make_figure_3_calibration_drug_class_share(agg, out)
    make_figure_4_calibration_infection_deaths(agg, out)
    make_figure_5_calibration_carriage_prevalence(agg, out)
    make_supplementary_figure_s1_potential_activity_retained(sf1_csv_paths, out, agg=agg)
    make_supplementary_figure_s2_microbiome_resistance_reservoir(runs, out, agg=agg)
    make_supplementary_figure_s3_carrier_vs_non_carrier_incidence(sf3_csv_paths, out, agg=agg)
    make_supplementary_figure_s4_resistance_mechanisms_by_bacterium(sf4_csv_paths, out, agg=agg)
    make_supplementary_figure_s5_diagnostic_testing_targeted_treatment_cascade(sf5_csv_paths, out, agg=agg)
    make_supplementary_figure_s6_new_active_infection_denominators(sf6_csv_paths, paths, out, agg=agg)
    make_supplementary_figure_s7_active_infection_incidence(agg, out)
    make_supplementary_figure_s8_infection_outcome_pathway(st1_csv_paths, out, agg=agg)
    make_figure_6_resistance_trend(csv_paths, out)
    make_figure_6b_resistance_trend_by_bacterium(csv_paths, out)
    make_figure_6c_serious_r_trend_by_bacterium(csv_paths, out)
    make_figure_20_serious_r_by_hospital_community(paths, out, agg=agg)
    make_figure_7_infection_death_rate_by_region(csv_paths, out, agg=agg)
    make_figure_8_antibiotic_use_by_context(csv_runs_with_scale, out, agg=agg)
    make_figure_11_sepsis_context_effective_therapy(csv_paths, out, agg=agg)
    make_figure_15_mean_activity_by_bacteria(csv_paths, out, agg=agg)
    make_figure_19_antibiotic_exposure_distribution(csv_paths, out, agg=agg)
    make_figure_10_resistance_pathway_counterfactuals(out, agg=agg)
    make_index(agg, out)

    print(f"\nDone. Open {out / 'index.html'} to browse paper outputs.")


if __name__ == "__main__":
    main(sys.argv[1:])
