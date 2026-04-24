"""
make_paper_tables.py

Generate HTML paper tables (T2–T6, S1–S4) from one or more
calibration_summary_*.txt files.

Usage
-----
Single run:
    python make_paper_tables.py output_graphs/calibration_summary_958282.txt

Multiple runs (pass explicit paths or a glob):
    python make_paper_tables.py output_graphs/calibration_summary_*.txt

Output
------
paper_tables/
    index.html                              — navigation page
    main/
        T2_headline_metrics.html
        T3_drug_class_share.html
        T4_bacteria_burden.html
        T5_resistance_fit.html
        T6_amr_attributable_deaths.html     — placeholder until counterfactual runs ready
    supplementary/
        S1_infection_deaths.html
        S2_syndrome_and_resistance_locus.html
        S3_resistance_benchmarks.html
        S4_calibration_score_details.html
"""

from __future__ import annotations

import glob
import io
import json
import math
import re
import sys
from pathlib import Path
from typing import Union

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
import pandas as pd

from parse_calibration import aggregate, parse_files

OUT_DIR = Path("paper_tables")

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


def _back_link() -> str:
    return "<p class='back-link'><a href='../index.html'>← Back to index</a></p>\n"


def _save(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content, encoding="utf-8")
    print(f"  Saved: {path}")


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


def _clean_df(df: pd.DataFrame) -> pd.DataFrame:
    """Drop delta columns; rename 'Target' → 'Observed estimate' in column names."""
    if df is None or df.empty:
        return df
    df = df.copy()
    drop = [c for c in df.columns if re.search(r'\bDelta\b|\bΔ\b', c, re.IGNORECASE)]
    df = df.drop(columns=drop, errors='ignore')
    def _rename_col(c: str) -> str:
        def repl(m):
            return 'Observed estimate' if m.group(0)[0].isupper() else 'observed estimate'
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
     "resistance prevalence (ECDC EARS-Net, WHO GLASS 2026); infection incidence and "
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
    _save(out_dir / "main" / "T1_model_summary.html", body)


# ---------------------------------------------------------------------------
# Table T2 — Headline Calibration Metrics + Block Scores
# ---------------------------------------------------------------------------

def _load_current_headline_targets() -> dict[str, float]:
    path = Path("data") / "calibration_targets.json"
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
        "Murray CJ et al. (2022). Global burden of bacterial antimicrobial resistance "
        "in 2019: a systematic analysis. <em>Lancet</em> 399:629–655.",
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

    # Filter on the H:C target BEFORE _clean_df, because _clean_df renames every
    # column containing "target" → "observed estimate", which would break the lookup.
    target_col = "Target H:C ratio"
    hc_col_present = target_col in srl.columns
    if hc_col_present:
        # Drop summary-stat rows (non-bacteria text in first column)
        first_col = srl.columns[0]
        summary_mask = srl[first_col].astype(str).str.match(
            r"^\s*(-|Resistance Locus|Serious Resistance|Mean |H:C)", na=False
        )
        srl = srl[~summary_mask].copy()

        srl[target_col] = pd.to_numeric(srl[target_col], errors="coerce")
        included_raw = srl[srl[target_col] > 1.0].copy()
        excluded_raw = srl[srl[target_col].notna() & (srl[target_col] <= 1.0)].copy()
    else:
        included_raw = srl.copy()
        excluded_raw = pd.DataFrame()

    srl = _clean_df(included_raw)
    ril = _clean_df(ril.copy()) if ril is not None and not ril.empty else pd.DataFrame()
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

    # Footnote: list the excluded organisms (target H:C = 1.0)
    if not excluded.empty and "Bacteria" in excluded.columns:
        excl_names = sorted(str(v) for v in excluded["Bacteria"].dropna().unique())
        excl_list = "; ".join(excl_names)
    else:
        excl_list = "none"

    footnotes = [
        _window_note(n),
        "This table includes only organisms for which there is at least some empirical suggestion "
        "that hospital-acquired cases carry higher resistance levels than community-acquired cases. "
        "Operationally, inclusion requires a literature-based target hospital:community (H:C) "
        "resistance ratio greater than 1.0 for the organism's marker drug. "
        f"Organisms with target H:C&nbsp;=&nbsp;1.0 are excluded ({excl_list}).",
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

    rb_nonneg = _clean_df(rb_nonneg)
    rb_nonneg = rb_nonneg.rename(columns={
        "Inf sim (%)":               "Percent of infections with resistance — simulation (%)",
        "Inf observed estimate (%)": "Percent of infections with resistance — observed estimate (%)",
        "Avg sim (%)":               "Average resistance level among resistant infection-days — simulation (%)",
        "Avg observed estimate (%)": "Average resistance level among resistant infection-days — observed estimate (%)",
        "Micro sim (%)":             "Percent of people carrying the bacterium in whom a resistant strain is present (%)",
    })

    display_cols = [c for c in [
        "Drug", "Class",
        "Percent of infections with resistance — simulation (%)",
        "Percent of infections with resistance — observed estimate (%)",
        "Average resistance level among resistant infection-days — simulation (%)",
        "Average resistance level among resistant infection-days — observed estimate (%)",
        "Percent of people carrying the bacterium in whom a resistant strain is present (%)",
    ] if c in rb_nonneg.columns]

    footnotes = [
        _window_note(n),
        "Only organism–drug combinations where the drug has non-negligible potency "
        "(baseline potency > 0.1) are shown.",
        "<em>Percent of infections with resistance</em>: percentage of active infections "
        "carrying any resistance to this drug at a point in time "
        "(simulated vs. surveillance estimate).",
        "<em>Average resistance level among resistant infection-days</em>: among infection-days "
        "where any resistance is present, the mean resistance level expressed as a percentage (0–100%). "
        "A value near 100% indicates that resistance, when present, is essentially complete; "
        "lower values indicate partial resistance. This is distinct from the prevalence column above, "
        "which measures the proportion of infection-days with any resistance.",
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
# Figure F1 — Historical resistance trend (1930–2025)
# ---------------------------------------------------------------------------

#: Calendar year corresponding to simulation time_in_years == 0.
_F1_SIM_EPOCH_YEAR: int = 1930

#: Colours for the trend figure.
_F1_TREND_COLOUR_MEAN   = "#1565C0"   # dark blue — mean line
_F1_TREND_COLOUR_CLOUD  = "#90CAF9"   # light blue — 90% CI band


def _load_resistance_series(csv_path: Path) -> pd.DataFrame | None:
    """
    Load a simulation_summary CSV and return a DataFrame with columns:
      year          — calendar year (float, 1930-based)
      pct_resistant — % of newly infected individuals carrying any resistance

    Returns None if the required columns are absent.
    """
    needed = ["time_in_years", "newly_infected_count", "newly_infected_with_resistance_count"]
    try:
        df = pd.read_csv(csv_path, usecols=needed)
    except (FileNotFoundError, ValueError):
        return None
    df = df.dropna(subset=needed)
    df = df[df["newly_infected_count"] > 0].copy()
    df["year"] = _F1_SIM_EPOCH_YEAR + df["time_in_years"]
    df["pct_resistant"] = (
        df["newly_infected_with_resistance_count"] / df["newly_infected_count"] * 100.0
    )
    return df[["year", "pct_resistant"]]


def make_f1_resistance_trend(csv_paths: list[Path], out_dir: Path) -> None:
    """
    Figure F1: time trend of the proportion of new bacterial infections in which
    resistance to any antibiotic was present at infection, 1930–2025.

    Each run in *csv_paths* contributes one time series.  The figure shows:
      - solid mean line across all runs
      - shaded 90% credible-interval cloud (5th–95th percentile across runs)

    If *csv_paths* is empty, or the available data cover fewer than 2 calendar
    years, a clearly labelled placeholder panel is saved instead.
    """
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    png_path = fig_dir / "F1_resistance_trend.png"
    svg_path = fig_dir / "F1_resistance_trend.svg"
    html_path = out_dir / "main" / "F1_resistance_trend.html"

    # ------------------------------------------------------------------ #
    # Load all run series                                                   #
    # ------------------------------------------------------------------ #
    series_list: list[pd.DataFrame] = []
    for p in csv_paths:
        s = _load_resistance_series(p)
        if s is not None and not s.empty:
            series_list.append(s)

    n_runs = len(series_list)
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
        # Resample every run to annual means, then stack
        annual_frames: list[pd.Series] = []
        for s in series_list:
            s = s.copy()
            s["year_int"] = s["year"].apply(int)
            annual = s.groupby("year_int")["pct_resistant"].mean()
            annual_frames.append(annual)

        combined = pd.concat(annual_frames, axis=1)   # rows=years, cols=runs
        years     = combined.index.values
        mean_vals = combined.mean(axis=1).values
        p5_vals   = combined.quantile(0.05, axis=1).values
        p95_vals  = combined.quantile(0.95, axis=1).values

        ax.fill_between(years, p5_vals, p95_vals,
                        color=_F1_TREND_COLOUR_CLOUD, alpha=0.55,
                        label="90% interval across runs")
        ax.plot(years, mean_vals,
                color=_F1_TREND_COLOUR_MEAN, linewidth=1.8, label="Mean")

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
                "Figure F1 — Historical resistance trend\n\n"
                "Data not yet available.\n"
                "Re-run with full-period simulation output\n"
                "(simulation_summary_*.csv covering 1930–2025)\n"
                "to generate this figure.",
                ha="center", va="center", transform=ax.transAxes,
                fontsize=11, color="#555",
                bbox=dict(boxstyle="round,pad=0.6", fc="#f5f5f5", ec="#bbb"))
        ax.set_axis_off()

    if data_available:
        ax.set_xlim(_F1_SIM_EPOCH_YEAR, _F1_SIM_EPOCH_YEAR + 96)
        ax.set_ylim(0, 100)
        ax.set_xlabel("Year", fontsize=11)
        ax.set_ylabel("New infections with any resistance (%)", fontsize=11)
        ax.spines[["top", "right"]].set_visible(False)
        ax.grid(axis="y", linewidth=0.4, alpha=0.5)
        n_label = f"n\u2009=\u2009{n_runs} run{'s' if n_runs > 1 else ''}"
        ax.legend(fontsize=9, frameon=False, title=n_label, title_fontsize=8)

    fig.suptitle(
        "Figure F1 \u2014 Proportion of new bacterial infections with any resistance, 1930\u20132025",
        fontsize=11, fontweight="bold",
    )
    fig.tight_layout()

    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {png_path}")
    print(f"  Saved: {svg_path}")

    # HTML wrapper
    html_rel_img     = "../figures/F1_resistance_trend.png"
    html_rel_img_svg = "../figures/F1_resistance_trend.svg"
    body  = _html_head("Figure F1 \u2014 Historical Resistance Trend")
    body += _back_link()
    body += "<h1>Figure F1 \u2014 Proportion of New Infections with Any Resistance, 1930\u20132025</h1>\n"
    if data_available:
        body += (
            f"<p class='note'>Proportion of newly infected individuals carrying any antibiotic-"
            f"resistant strain, plotted annually from {int(year_min_data)} to "
            f"{int(year_max_data)}. "
            f"Solid line: mean across {n_runs} accepted simulation run"
            f"{'s' if n_runs > 1 else ''}. "
            f"Shaded band: 5th\u201395th percentile (90% interval) across runs.</p>\n"
        )
    else:
        body += (
            "<p class='note' style='color:#c0392b;'>"
            "Placeholder \u2014 full-period simulation data (1930\u20132025) not yet available. "
            "Run the simulation in non-calibration mode and supply the resulting "
            "<code>simulation_summary_*.csv</code> files to generate this figure.</p>\n"
        )
    body += (
        f"<p><a href='{html_rel_img_svg}' target='_blank'>[Download SVG]</a></p>\n"
        f"<img src='{html_rel_img}' alt='Figure F1' "
        f"style='max-width:100%; border:1px solid #ddd; border-radius:4px;'>\n"
    )
    body += _html_footnotes([
        "Resistance presence is defined as the individual carrying at least one bacterial strain "
        "with any resistance at the time of infection acquisition.",
        "Values are averaged within each calendar year across all daily time steps.",
        "Full 1930\u20132025 data require a simulation run spanning the complete historical period, "
        "not just the 2022\u20132025 calibration window.",
    ])
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

# Colour scheme
_F2_COLOUR_SIM    = "#2196F3"   # blue — simulation
_F2_COLOUR_TARGET = "#FF7043"   # deep orange — surveillance target


def _parse_resistance_val(v: object) -> tuple[float, float, float] | None:
    """
    Parse a resistance value from resistance_benchmarks into (median, lo, hi).

    Handles:
      - float/int (single run, no CI)
      - "12.3 (10.1–14.5)" or "12.3 (10.1-14.5)"  (aggregated multi-run)
      - "—", "", None, NaN → return None
    """
    if v is None:
        return None
    if isinstance(v, float) and np.isnan(v):
        return None
    if isinstance(v, (int, float)):
        f = float(v)
        return (f, f, f)
    s = str(v).strip()
    if s in ("—", "-", "", "N/A", "nan"):
        return None
    # "12.3 (10.1–14.5)" with em-dash or hyphen
    m = re.match(r"([\d.]+)\s*\(\s*([\d.]+)\s*[–\-]\s*([\d.]+)\s*\)", s)
    if m:
        return (float(m.group(1)), float(m.group(2)), float(m.group(3)))
    # Plain numeric string
    try:
        f = float(s.replace(",", ""))
        return (f, f, f)
    except ValueError:
        return None


def _class_summary(
    rows: pd.DataFrame,
    sim_col: str,
    tgt_col: str,
) -> tuple[float | None, float | None, float | None, float | None]:
    """
    Summarise all drug rows for one class:
    Returns (sim_median, sim_lo, sim_hi, tgt_mean).
    Uses mean of per-drug medians for bar height; [min-lo, max-hi] for error span.
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


def make_f2_resistance_barplot(agg: dict, out_dir: Path) -> None:
    """
    Create Figure F2: dynamic grid showing infection resistance calibration fit
    for all organisms present in the resistance_benchmarks data.

    Each panel:
      x-axis — drug classes present for that organism
      y-axis — % resistant infections
      Blue bar + error bars — simulation (median ± 5th–95th percentile range)
      Orange bar           — surveillance target
    """
    rb = agg.get("resistance_benchmarks", pd.DataFrame())
    if rb is None or rb.empty:
        print("  F2: no resistance_benchmarks data — skipping figure.")
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
    tgt_patch = mpatches.Patch(color=_F2_COLOUR_TARGET, alpha=0.85, label="Surveillance target")
    ci_note   = (
        f"Error bars: 5th–95th percentile across {n_runs} run{'s' if n_runs > 1 else ''}."
        if n_runs > 1 else "Single run; no uncertainty interval shown."
    )
    fig.legend(
        handles=[sim_patch, tgt_patch], loc="lower center",
        ncol=2, fontsize=9, frameon=False, bbox_to_anchor=(0.5, 0.0),
    )
    fig.suptitle(
        "Figure F2 — Infection resistance calibration fit by organism",
        fontsize=11, fontweight="bold", y=1.01,
    )
    fig.text(0.5, -0.01, ci_note, ha="center", fontsize=7.5, color="#555")

    fig.tight_layout(rect=[0, 0.04, 1, 1])

    # Save PNG and SVG
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    png_path = fig_dir / "F2_resistance_fit.png"
    svg_path = fig_dir / "F2_resistance_fit.svg"
    fig.savefig(png_path, dpi=150, bbox_inches="tight")
    fig.savefig(svg_path, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {png_path}")
    print(f"  Saved: {svg_path}")

    # HTML wrapper
    html_rel_img = "../figures/F2_resistance_fit.png"
    html_rel_img_svg = "../figures/F2_resistance_fit.svg"
    body  = _html_head("Figure F2 — Resistance Calibration Fit")
    body += _back_link()
    body += "<h1>Figure F2 — Infection Resistance Calibration Fit by Organism</h1>\n"
    body += _meta_box(agg)
    body += (
        "<p class='note'>Each panel shows the simulated (blue) and surveillance-target (orange) "
        "infection resistance percentage by drug class for one priority organism. "
        "Simulation bars show the mean across all drugs in the class; "
        f"error bars span the 5th–95th percentile range across {n_runs} accepted run"
        f"{'s' if n_runs > 1 else ''}. "
        "Classes without data for a given organism are omitted.</p>\n"
    )
    body += (
        f"<p><a href='{html_rel_img_svg}' target='_blank'>[Download SVG]</a></p>\n"
        f"<img src='{html_rel_img}' alt='Figure F2' "
        f"style='max-width:100%; border:1px solid #ddd; border-radius:4px;'>\n"
    )
    body += _html_footnotes([
        "Drug class resistance within a panel is averaged across all specific drugs in that class.",
        "Drugs marked 'negligible potency' in the simulation are excluded from class averages.",
        "All organisms with resistance benchmark data in the simulation output are included. "
        "IHME/WHO-ESKAPE priority organisms are shown first, remainder alphabetically.",
    ])
    body += "</body></html>"
    html_path = out_dir / "main" / "F2_resistance_fit.html"
    _save(html_path, body)


# ---------------------------------------------------------------------------
# Index page
# ---------------------------------------------------------------------------

def make_index(agg: dict, out_dir: Path) -> None:
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
         "Organisms with H:C target &gt; 1: hospital-acquired%, any-R%, serious-R% by locus"),
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
         "12-panel bar chart: simulated vs. surveillance-target infection resistance "
         "by drug class for IHME-priority pathogens"),
    ]:
        body += f"<li><a href='{fname}'><strong>{label}</strong></a> — {desc}</li>\n"
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
        body += f"<li><a href='{fname}'><strong>{label}</strong></a> — {desc}</li>\n"
    body += "</ul>\n"
    body += "</body></html>"
    _save(out_dir / "index.html", body)


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------

def main(input_args: list[str]) -> None:
    if not input_args:
        print("Usage: python make_paper_tables.py <calibration_summary_*.txt> [...]")
        sys.exit(1)

    # Expand any glob patterns
    paths: list[str] = []
    for arg in input_args:
        expanded = glob.glob(arg)
        if expanded:
            paths.extend(expanded)
        elif Path(arg).exists():
            paths.append(arg)
        else:
            print(f"  Warning: no file found for '{arg}'")

    if not paths:
        print("No input files found.")
        sys.exit(1)

    print(f"Parsing {len(paths)} calibration file(s)...")
    runs = parse_files(paths)
    agg  = aggregate(runs)
    n    = len(runs)
    print(f"  -> {n} run(s) parsed and aggregated.")

    # Auto-discover matching simulation_summary CSVs for F1 trend figure.
    # Calibration files are named calibration_summary_{seed}.txt;
    # matching CSVs are simulation_summary_{seed}.csv in the output analysis folder.
    _csv_dir = Path("amr_simulation_output_analysis_outputs")
    csv_paths: list[Path] = []
    for cal_path in paths:
        stem = Path(cal_path).stem                 # e.g. "calibration_summary_958282"
        seed = stem.split("_")[-1]                 # e.g. "958282"
        candidate = _csv_dir / f"simulation_summary_{seed}.csv"
        if candidate.exists():
            csv_paths.append(candidate)
    if csv_paths:
        print(f"  Found {len(csv_paths)} simulation CSV(s) for F1 trend figure.")
    else:
        print("  No matching simulation CSVs found; F1 will render as placeholder.")

    out = OUT_DIR
    print(f"\nGenerating tables in {out.absolute()} ...")
    make_t1(out)
    make_t2(agg, out)
    make_t3(agg, out)
    make_t4(agg, out)
    make_t5(agg, out)
    make_t6_placeholder(out)
    make_s1(agg, out)
    make_s2(agg, out)
    make_s3(agg, out)
    make_s4(agg, out)
    make_f1_resistance_trend(csv_paths, out)
    make_f2_resistance_barplot(agg, out)
    make_index(agg, out)

    print(f"\nDone. Open {out / 'index.html'} to browse all tables.")


if __name__ == "__main__":
    main(sys.argv[1:])
