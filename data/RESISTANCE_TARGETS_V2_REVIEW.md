# Resistance reference review: resistance_targets_v2

Review date: 28 September 2026. The active reference version is `resistance_targets_v2`.
This is a reference correction, with no change to biological model parameters or existing
simulation measurements. A score change caused by revising or withdrawing a reference
does not demonstrate improved simulated biology.

## Decisions and units

Prevalence is stored as a proportion: `0.40` means 40%, and `0.0` is an assigned numeric
zero. Blank long-form values are unassigned; historical wide matrices use `.` for missing.
Conditional severity is a separate, unitless model-scale quantity, conditional on
`any_r > 0`. Prevalence values must not be copied into severity.

The three T. pallidum tetracycline-family references now have a pragmatic near-zero
structural prior centred at `0.0`. This is an explicit modelling convention for the
contemporary reference window, not a measured global prevalence of exactly zero.
Supporting evidence primarily concerns doxycycline/tetracycline; minocycline uses a
class-based prior without independent surveillance. Uncertainty limits and evidence
weights remain unassigned. The preserved conditional-severity numbers are labelled
rare-positive structural priors: they describe severity if simulated positive infections
occur. Existing potency and component eligibility restrictions still apply, and zero
uses ordinary finite scoring rather than a hard constraint.

E. faecium tigecycline and linezolid are unassigned pending recovery of defensible
all-region resistant proportions and susceptibility criteria. No new positive numerical
reference is assigned. The inspected ATLAS abstract does not establish exact percentages;
neither susceptibility above 90%, its complement, nor a selected regional extreme is
used as an exact resistance estimate. Country-specific evidence is not transformed into
a global centre. The existing tedizolid prevalence `0.10` remains unchanged and unresolved
for separate review. Conventional tetracycline `0.50` is unchanged. Ampicillin and
amoxicillin remain missing, with no transfer between aminopenicillins.

All 13 positive C. trachomatis entries had only the generic legacy WHO GLASS/GRAM note,
without recovered compatible species-drug provenance. Each was withdrawn individually;
the inspected local files contained no newer cell-level source to preserve for these
entries. They are missing, not zero. Treatment failure, persistence, and experimental
inhibitory concentrations do not establish globally representative acquired-resistance
prevalence. The Korean study provides context, not a basis for extrapolation across drugs.

The source observation is a clinical isolate or genomic sample. The model observation
is an active-infection person-day with `any_r > 0`. These denominators and definitions
differ; a small positive model effect is not equivalent to clinical breakpoint resistance.
No observation-model or resistance-mechanism redesign is made by this patch.

## Row-level changes

The following prevalence values are proportions. “Unassigned” means blank and excluded
from scoring and paired reference averages. All listed prevalence rows were previously
included. “Prior” means structural prior, not direct empirical estimate.

| Organism | Drug | Previous prevalence | V2 prevalence | Revision | Preserved severity | V2 severity eligibility |
|---|---|---:|---:|---|---:|---|
| T. pallidum | tetracycline | 0.08 | 0.0 | Near-zero prior | 0.62 | Included; rare-positive prior |
| T. pallidum | doxycycline | 0.06 | 0.0 | Near-zero prior | 0.60 | Included; rare-positive prior |
| T. pallidum | minocycline | 0.05 | 0.0 | Class-based near-zero prior | 0.58 | Included; rare-positive prior |
| E. faecium | tigecycline | 0.40 | Unassigned | Withdrawn pending review | 0.78 | Unpaired; already excluded before review |
| E. faecium | linezolid | 0.10 | Unassigned | Withdrawn pending review | 0.62 | Unpaired; newly excluded |
| C. trachomatis | erythromycin | 0.10 | Unassigned | Unsupported legacy reference withdrawn | 0.62 | Unpaired; newly excluded |
| C. trachomatis | azithromycin | 0.05 | Unassigned | Unsupported legacy reference withdrawn | 0.60 | Unpaired; newly excluded |
| C. trachomatis | clarithromycin | 0.10 | Unassigned | Unsupported legacy reference withdrawn | 0.62 | Unpaired; newly excluded |
| C. trachomatis | clindamycin | 0.20 | Unassigned | Unsupported legacy reference withdrawn | 0.68 | Unpaired; newly excluded |
| C. trachomatis | ciprofloxacin | 0.10 | Unassigned | Unsupported legacy reference withdrawn | 0.62 | Unpaired; newly excluded |
| C. trachomatis | levofloxacin | 0.05 | Unassigned | Unsupported legacy reference withdrawn | 0.60 | Unpaired; newly excluded |
| C. trachomatis | moxifloxacin | 0.05 | Unassigned | Unsupported legacy reference withdrawn | 0.60 | Unpaired; newly excluded |
| C. trachomatis | ofloxacin | 0.10 | Unassigned | Unsupported legacy reference withdrawn | 0.62 | Unpaired; newly excluded |
| C. trachomatis | tetracycline | 0.05 | Unassigned | Unsupported legacy reference withdrawn | 0.60 | Unpaired; newly excluded |
| C. trachomatis | doxycycline | 0.05 | Unassigned | Unsupported legacy reference withdrawn | 0.60 | Unpaired; newly excluded |
| C. trachomatis | minocycline | 0.05 | Unassigned | Unsupported legacy reference withdrawn | 0.60 | Unpaired; newly excluded |
| C. trachomatis | tigecycline | 0.05 | Unassigned | Unsupported legacy reference withdrawn | 0.60 | Unpaired; already excluded before review |
| C. trachomatis | chloramphenicol | 0.10 | Unassigned | Unsupported legacy reference withdrawn | 0.62 | Unpaired; newly excluded |

The machine-readable `resistance_targets_v2.revisions.csv` contains all 36 changed rows:
18 prevalence revisions and 18 severity status/provenance revisions. It records previous
and new values, statuses, inclusion flags, exclusion reasons, provenance classes, source
identifiers, rationale, and checked-source identifiers. Every severity number is preserved.
No other organism-drug-component metadata or values changed beyond the common version tag.
This includes the existing S. aureus glycopeptide references.

Static scoring inclusion changes from 1,176 to 1,161 prevalence rows, and from 1,056 to
1,043 severity rows. Two of the 15 newly unpaired severity rows were already excluded
because their values exceeded the model-representable maximum. They retain those
restrictions and now also use `inactive_unpaired_legacy_benchmark`. Missing prevalence
continues to use the existing `legacy_unclassified_missing` status; withdrawal is made
explicit through `not_assigned`, the review source, rationale, and audit rather than an
unsupported new status.

Paired class comparisons use the same eligible drugs for simulation and reference means.
Thus E. faecium tetracycline-family comparisons lose tigecycline while retaining
eligible conventional members; oxazolidinone comparisons lose linezolid but do not alter
tedizolid. All C. trachomatis paired class comparisons become unavailable because no
assigned prevalence references remain. T. pallidum tetracycline-family references retain
their eligible membership and show explicit zero markers. Per-run membership and aggregate
changes are recorded in the version-labelled comparison output; simulation-only summaries
remain separately labelled where available.

## Checked sources and limitations

Machine-readable source observations are in `resistance_targets_v2.review.json`, with
publication year separated from observation years, geography, setting, denominator,
breakpoint criteria, access level, and transformation decision. Missing information is
recorded as unknown, not inferred. No confidence interval or evidence-quality weight is
invented.

The BASHH full guideline HTML was checked. It reports no observed doxycycline resistance
in T. pallidum at publication. The CDC full article and Table 1 caption were checked; its
genomic observations support surveillance context, not an exact phenotypic global zero.
The prior decision remains a modelling judgement.

The ATLAS PubMed abstract was recovered through indexed primary-source search, with
publisher highlights also available. ScienceDirect full text returned HTTP 403 and the
JGAR full-text endpoint was inaccessible; tables and supplements were not recovered.
Observation years are 2019–2023, publication year 2026. The study covers 6,036 E. faecium
isolates in 59 countries across six continents, but drug-specific tested denominators,
breakpoint version, and exact all-region resistant results remain unrecovered. Both
replacement references are therefore unresolved.

The full AGAR publisher PDF, methods, Table 2, and breakpoint reference were checked.
Australian bloodstream surveillance in 2024 reports one resistant E. faecium isolate among
593 tested for linezolid (rounded 0.2%), using EUCAST criteria with version 15.0 cited.
This is a restricted population, not a global transformation. The Klimik publisher
abstract and reference list were checked; its full Turkish article was not checked.
It reports tigecycline resistance of 1.2% and no detected linezolid resistance among
E. faecium from Hacettepe University Hospitals during 2010–2023. Those numbers remain
context only.

The Korean study's PMC endpoint returned a browser challenge; the full publisher PDF
was checked instead. It reports 118 recovered female-patient isolates from 2016–2019 and
experimental susceptibility distributions for four drugs. Its results are not a harmonised
global acquired-resistance prevalence series. No value is inferred for another drug.

## Version-bound inputs and reproduction

The v1 snapshot contains the original target/source/schema/manifest and all four original
wide/projection inputs. Snapshot bytes match the preserved live v1 files at review time.
The v2 review pins the snapshot manifest hash, and the builder verifies every snapshot
artifact before applying decisions. Expected old values and source identifiers guard each
decision against accidental replacement of newer evidence.

`build_resistance_targets_v2.py` reads the snapshot and review decisions, then applies the
existing eligibility rules to the current Rust potency and reachability projections. The
live legacy wide matrices cannot restore withdrawn values. V2 source, target, and revision
CSVs are generated companions. The manifest covers all version-bound inputs, the schema,
these outputs, and the current projections. The production loader verifies this contract
without relaxing v1 verification.

For active refresh, run from the project root:

```powershell
python -m amr_simulation_output_analysis.refresh_resistance_targets
```

Refresh stages the inputs and freshly exported Rust projections, validates the generated
rows and source registry, and publishes under `.resistance_targets_v2.update.lock`, with
the manifest last. For isolated development, `build_resistance_targets_v2(root=...)` builds
from a root containing the declared `INPUT_FILENAMES` under `data/`, plus the two current
projections. The focused tests confirm byte-identical regeneration and corruption rejection.

Historical v1 can be reproduced by copying `data/resistance_targets_v1_snapshot/` into
an isolated root as `data/`, then calling `build_resistance_targets_v1(root=isolated_root)`.
The target/source/manifest bytes reproduce the preserved v1 set. Do not refresh historical
inputs from the current projections when reproducing an old comparison.

Existing `calibration_summary_*.txt` files contain embedded historical target values.
They remain historical records. Selecting v2 overlays the current reference and preserves
the original run/reference information in the comparison audit; embedded targets cannot
override the selected set. Version-labelled comparison output belongs in
`paper_tables_resistance_targets_v2/`, separate from the archived `paper_tables/` bundle.
The selected original runs are `fff164402`, `fff311750`, `fff591413`, and `fff928228`.
Run measurements and the existing run-selection and mean/CI rules remain unchanged.

## References

Saunders J, et al. BASHH UK national guideline for doxycycline post-exposure prophylaxis for prevention of syphilis, 2025. doi:10.1177/09564624251352053. https://journals.sagepub.com/doi/10.1177/09564624251352053

Long GS, et al. Genomic Analysis of Doxycycline Resistance-Associated 16S rRNA Mutations in Treponema pallidum Subspecies pallidum. Emerging Infectious Diseases. 2026;32(2):242–245. doi:10.3201/eid3202.251060. https://wwwnc.cdc.gov/eid/article/32/2/25-1060_article

Liu Y, Zheng X, Wu B. ATLAS analysis of Enterococcus faecalis/faecium susceptibility, 2019–2023. Journal of Global Antimicrobial Resistance. 2026;46:254–263. doi:10.1016/j.jgar.2025.12.018. PMID:41490597. https://pubmed.ncbi.nlm.nih.gov/41490597/

Coombs G, et al. AGAR Australian Enterococcal Surveillance Outcome Program Bloodstream Infection Annual Report 2024. Communicable Diseases Intelligence. 2025;49. doi:10.33321/cdi.2025.49.053. PMID:41248469. https://ojs.cdi.cdc.gov.au/index.php/cdi/article/download/3420/4985/7267

Altun B, Hazirolan G, Gur D. Antibiotic Resistance of Enterococcus faecium and Enterococcus faecalis (2010–2023). Klimik Journal. 2025;38(1):13–18. doi:10.36519/kd.2025.5110. https://www.klimikdergisi.org/en/2025/03/21/antibiotic-resistance-of-enterococcus-faecium-and-enterococcus-faecalis-2010-2023/

Kang T, Choi YJ, Jang WJ. A Nationwide Investigation on the Distribution and Antimicrobial Susceptibility Profiling of Chlamydia trachomatis Collected from Patients in Korea. Infection & Chemotherapy. 2026;58(1):39–48. doi:10.3947/ic.2025.0081. https://icjournal.org/pdf/10.3947/ic.2025.0081
