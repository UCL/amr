# Observations and questions about the AMR simulation model

Prepared by an external analysis group, 2026-09-04.
This document describes behavior of **the model at the current repository tip
(`91c5c1aa9ac7b5ed1aa035e709b8c0559ec3d401`, 2026-09-02)** and, where noted, of the
older snapshot `ed0243b4` (2026-07-27) that we used for an earlier round of analysis.
Where a figure comes from our own second engine (an independent reimplementation of your
model's specification), we say so explicitly; everything else comes from unmodified
builds of your repository at the commits named, and every figure can be traced to a
checksummed output file (verification chain in the appendix).
We are raising **questions**, not making defect claims.

---

## 1. What we ran and observed

- **50 full-horizon runs** of the current tip: population 100,000, 35,040 passes (~96
  years), full output mode, 50 different RNG seeds (13 reference seeds and 37 expansion
  seeds, listed in the appendix), release build. Each run's complete daily output
  (33,203 columns, one row per day) was preserved and checksummed.
- **3 large-population runs** (~6.0M peak population each, same model, full horizon;
  file ids 42ab3663bf6e, 91683220c2a1, b78fa314e5ed) — provided by your team — which we
  used to test whether small counts that are zero at 100k are structural zeros or
  sampling floors.
- We also ran the same 50 seeds with **our second engine** (an independent
  reimplementation of your model's specification, built to compare behavior). Where a
  number below is from that second engine, it is labeled. The questions, however, are
  about **your model's** intended semantics.

## 2. Diagnostic-cascade accounting — confirmation of fixes, and one historical note

We previously reported an ordering issue in the diagnostic-cascade counters (an
unidentified bacterium inheriting another bacterium's targeted context, so
`resistance_testing_done` / `targeted_treatment_started` / `effective_targeted_treatment_started`
stages could be recorded out of order). We have verified in your repository history
that this is fixed at the current tip:

- `a3e463d` (2026-08-12, "Prevent post-death diagnostic cascade attribution") — stops
  post-death attribution; and
- `128d29c` + `2eef5b4` (2026-08-28) — both stage predicates identification-gated, with
  regression tests pinning effective ⇒ targeted ⇒ identification.

The verification rests on your commits' code and tests. For completeness, in our 50
current-tip runs the cascade collection is **enabled** — the
`diagnostic_cascade_collection_enabled` column reads 1 on every day of every run — and
the counters are active throughout the horizon (in one run — seed
3833834968269084703 — eligible max 56 with 29,562
nonzero days, identification max 25, resistance-testing max 13, targeted max 24,
effective max 13); the final-day values happen to be small or zero
(e.g. eligible 16, identification 1, the three downstream stages 0 on that day).

Historical note: our earlier comparisons used the older snapshot `ed0243b4` (2026-07-27),
which predates both fixes; in that snapshot we measured substantial entry-cohort churn
(19.1 entries per infected-day across 13 runs, from the per-day cascade columns; the
verification artifact is `CHURN.json`, sha256 cited in the appendix).

We also record your team's answer to our earlier questions, for completeness:
entry-cohort attribution (later stages written back to the episode's eligibility day) is
intended; each stage is counted at most once per episode; counters are bacterium-episode
level; treatment before AST is an intended sibling outcome; the counters never affect
simulation behavior.

## 3. Question 1 — at which instant is the hospital/community split intended to be read?

The per-bacterium `currently_infected_hospital_count` / `currently_infected_community_count`
columns are produced in the per-person counting pass, using the person's live
`hospital_status.is_hospitalized()` at that point in the timestep
(`src/simulation/simulation.rs:5874` and `:4916` at the tip).

**Observation (cross-engine).** When we compared the same seeds between your model and
our second engine — both implementing "current care location", but reading it at
different instants within the day — we found a small, persistent difference in the
split (all figures below are means over the 50 same-seed run pairs). Over the late horizon it is consistently +1 to +3 people per bacterium attributed
to the hospital side in our engine, growing slowly; earlier in the horizon it is noisy
(campylobacter dips as low as −0.22 around day 7,500; chlamydia stays below +1 for most
of the first half). The **totals** (hospital + community per bacterium, occupancy, total
infected) agree to noise throughout. This is the signature expected if care-location
events (admissions, discharges) that occur between the two reading instants are
attributed to different days.

**Question.** Is the counting-pass instant (after the day's care-location events have all
been applied) the intended semantics for these columns? Is that instant documented
anywhere? If a second implementation wants to reproduce these columns bit-exactly, what
is the canonical point in the day to read care location?

## 4. Question 2 — amplification of a rare resistant infection through the carriage profile pool

**Background (your model's machinery).** New carriage acquisitions can import a complete
resistance profile from the circulating profile cache (`mechanism_cache.sample_profile`,
`src/rules/mod.rs:5371` for carriage at the tip), and a per-mechanism local persistence
archive can supply an established genotype with probability
`min(virtual_mass / (n + virtual_mass), 0.10)` (`src/simulation/simulation.rs`,
`local_persistence_sampling_probability`). The archive law (first complete observed
genotype carrying each mechanism, never replaced) behaved as documented in our runs.

**Observation from your model's own runs.** Scanning the
`streptococcus_agalactiae_infected_with_mutation_pbp_mosaic` column across all 50 of our
current-tip runs of your model: exactly **two runs** contain an *S. agalactiae*
infection carrying `mutation_pbp_mosaic` — one person for days 24,360–24,361 (seed
11060706022568330210) and one person for days 16,542–16,546 (seed 7851651572952336106).
In the first of those runs, the β-lactam resistant-carriage columns
(`streptococcus_agalactiae_microbiome_r_positive_*`) stayed between 0 and 11 people
from their first nonzero day (~day 12,273; i.e. for the remaining ~22,767 days of the
horizon; measured from the establishing infection's last day, 24,361, the bounded
carriage window is ~10,679 days) — the amplification never took off, despite the same
archive and import machinery.

**Observation from our second engine (labeled).** In our reimplementation, one run
(seed 3833834968269084703) traced the full chain event-by-event: a single *S. agalactiae*
**infection** acquired `mutation_pbp_mosaic` on day 24,563 (in-infection mutation under
β-lactam selection; your config rate `3e-5`), the genotype entered the circulating cache
on day 24,564, carriage imports of the pbp_mosaic-bearing profile began day 24,565, and
by day 35,039 **623 people** carried β-lactam-resistant *S. agalactiae* (all 20
penicillin/cephalosporin `microbiome_r_positive` columns reading identical counts). A
second run of ours established late (day ~34,530; final count 10) and a third only
transiently (one day; final count 0). The establishment rate itself matched your model's
(2 of your runs, 3 of ours, out of 50 each); what differed was the amplification after
establishment. In our runs the establishing infections persisted 59 days (the amplifying
run) and 19 days (the late run); in yours, 2 and 5 days.

**Question.** Both the bounded outcome (your model: ≤11 carriers sustained for years) and
the large outcome (our engine: 623 carriers in ~10,500 days, accelerating at the end) are
consistent with the same cache and archive laws as far as we can tell; the difference
appears to ride on how long the establishing infection persists and keeps depositing
into the cache. Is a bounded outcome the intended behavior — i.e., is there a known
bound on how much a single archived genotype can amplify through carriage-acquisition
imports? Or is amplification of this kind expected at these rates, with small-population
runs simply rarely establishing? We would also be glad to know whether the
establishing-infection duration (2–5 days in your runs; 19–59 in ours) is a quantity you
have characterized.

## 5. Three cells that stayed zero at both population scales

Across the 50 current-tip runs at 100k **and** the 3 large-population runs (~6M, final
~4 years), the following were zero on every observed day. We report these as observations
only; they may simply be rarer than one event per ~6M people per ~1,460 days:

1. `campylobacter_jejuni_currently_on_drug_aztreonam_avibactam`
2. `chlamydia_trachomatis_infected_with_enzyme_cat`
3. `clostridioides_difficile_microbiome_r_positive_fidaxomicin`

(Caveat: the large-population outputs cover only days 33,580–35,039, so activity earlier
in the horizon would not be visible to us.)

## 6. What we are not claiming

- We are not claiming any defect in the current tip. The ordering/churn issues we reported
  earlier were fixed in your repository before the current tip, and we have verified the
  fixes in your code and tests.
- We are not asserting what the "right" answer is for either question above. Both concern
  semantics your code and comments implement precisely; our questions are about intent.
- We could not find a place in the repository documentation that states the counting
  instant for the hospital/community split (Question 1) or the expected amplification
  behavior of the persistence archive (Question 2). If such documentation exists and we
  missed it, we would be grateful for a pointer — that would answer both questions
  directly.
- Everything in this document is reproducible from unmodified builds of your repository
  at the cited commits (and, where labeled, from our second engine under the same
  configuration); our internal verification chain (run checksums and analysis scripts)
  is in the appendix and we are happy to share any of it.

---

## Appendix — our verification chain (internal, for reproducibility)

- **The 50 seeds**: 13 reference seeds 3809718317527778200, 11060706022568330210,
  2481154537346185106, 3022473303077927265, 18347766894618398142, 1743706249353741536,
  3076749006915123828, 1769517158426135365, 14771631115214586733, 14629093910719143219,
  7067765765970965996, 5654162616053096720, 8319410004210104027; plus 37 expansion seeds
  (complete list: 130946020741174916, 210944879901308425, 1579037585477334267,
  2394380033687139321, 3833834968269084703, 3890482914138087643, 4094287276153474633,
  4303517363679454187, 6053568347508581821, 6123623731046892286, 6187542720390436348,
  6243516156600009763, 6380790934850758466, 6977788455006284434, 7681091453701726579,
  7851651572952336106, 8004248509054568370, 8079806453050359259, 8392487890181735596,
  8939636767791906524, 9963118203015443267, 10981310032587039608, 11992120946073908750,
  12353354734928645595, 12838881076068827306, 13327053708039856915, 13439474077587650685,
  13772517095535012135, 14543018133750326352, 14677337251484143550, 14716506744925582465,
  16058959731774348343, 16211927700108261052, 16375040998042055813, 16618598027411416590,
  17806458247871106239, 17865820007111518327; run set
  `AMR-V1-TIP91C5C1A-FULLOUTPUT-20260903-R1`, per-run sha256 receipts; ledger over 29,325
  comparable output metrics derived from the 33,203 columns: `LEDGER-TIP54.csv`, sha256
  `9fc54774d666f8dd9a46b29a73e5bc53f08a4c9ed9e3b9b8efb23294cdeaa2ac`).
- **Large-population runs**: your files 42ab3663bf6e, 91683220c2a1, b78fa314e5ed (sha256
  1698396507…, cae996a109…, aaf34a835…; floor-test results `FLOOR_RESULTS.json`
  sha256 `bac4c56f24cc3a20621c0f480d859be8d32138ef21dfc75932960c8ba0243b05`).
- **Cascade fix verification**: `UPSTREAM_CHECK_R2.md`, sealed dir
  `/home/someone/scratch/amr-upstream-cascade-check-20260903-r1/` (review of your commits
  a3e463d, 128d29c, 2eef5b4). Frozen-snapshot churn measurement: `CHURN.json` in
  `/home/someone/scratch/amr-eq01-cascade-x4-20260903-r1/` (13-run, full-series; your model's
  19.137 entries per infected-day).
- **Hospital/community split trajectory (cross-engine)**: sealed dir
  `/home/someone/scratch/amr-eq01-f1f3-falsifiers-20260904-r1/` (per-day series from all
  50 pairs; independent re-derivation confirmed all 16 time points).
- **Second-engine amplification chain**: event-level telemetry replay, sealed dir
  `/home/someone/scratch/amr-eq01-f3-inflane-20260904-r1/` (verdict `a51920b6…`), plus
  the scan of your model's own `infected_with_mutation_pbp_mosaic` columns across all 100
  outputs (`pbp_scan/results.jsonl`, sha256 `6ec0922b…`).
