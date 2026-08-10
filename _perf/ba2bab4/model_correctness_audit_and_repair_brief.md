# AMR Model Correctness Audit and Repair Brief

## Purpose

This is a self-contained handoff for an agent auditing and repairing the AMR
model itself. It assumes no knowledge beyond the repository. The task is
correctness, scientific meaning, output meaning, and regression protection.

The evidence below was frozen from source commit
`54513895b1f31728f79256e1064a02b8309efadc` on 2026-07-21. The remote may have
advanced since that observation. Before changing anything, record the current
`git rev-parse HEAD`, confirm whether each cited behavior still exists, and
update the evidence boundary in the final report.

This brief consolidates the model-internal defects and concerns established by
the audit work completed through the stated evidence date. It is not a claim
that every model pathway has already received the same depth of correctness
audit.

This audit is deliberately honest about certainty:

- **confirmed active defect** means the implementation mechanism is wrong
  regardless of the unresolved scientific policy;
- **confirmed contradiction** means code, comments, names, or history cannot
  all describe the same behavior;
- **contract decision required** means the behavior is reachable and
  suspicious, but the correct scientific/output definition has not yet been
  established; and
- **closed regression item** means the audited source already contains the
  repair, so only verification or a regression guard may be needed.

Do not silently choose a contract merely because it is easiest to implement or
matches existing output. Where intent is unresolved, present the alternatives
and obtain a model/output-owner decision. Any output semantic change must be
documented and, where necessary, versioned.

## Executive issue list

| ID | Area | Status | Immediate action |
| --- | --- | --- | --- |
| `AMR-AUDIT-001` | microbiome species masks | **confirmed active defect** | replace literal bacterium ordinals with an explicit named policy, after resolving the intended mask per output family |
| `AMR-AUDIT-002` | clinical species exclusion | **confirmed contradiction; policy unresolved** | decide whether each field excludes H. pylori, Treponema, both, or neither; then align logic, names, comments, and tests |
| `AMR-AUDIT-003` | new-infection incidence | **contract decision required** | choose acquisition-event incidence or end-of-day surviving new infection and test same-day clearance explicitly |
| `AMR-AUDIT-004` | sepsis incidence | **contract decision required** | define person-level episode incidence and bacterium-level onset incidence separately, including same-day recovery and death |
| `AMR-AUDIT-005` | toxicity-stop incidence | **contract decision required** | choose all stop events or survivor-only stops; do not implement the distinction by incrementing and retracting a scalar |
| `AMR-AUDIT-006` | active-infection threshold | **contract decision required** | decide whether any positive residual load or only load above `INFECTION_EPS` qualifies a drug initiation as infection-associated |
| `AMR-AUDIT-007` | processing after death | **closed regression item** | verify the terminal return remains ordered immediately after death recording and keep a no-post-mortem-mutation/RNG guard |
| `AMR-AUDIT-008` | repository authority and stale metadata | **confirmed active documentation/tooling risk** | correct executable counts and quarantine or regenerate tracked configuration copies that use retired mechanism names |
| `AMR-AUDIT-009` | infection-death analysis fallback | **confirmed approximation** | require dedicated model-scope counters where possible and make approximate subtraction explicit and opt-in |
| `AMR-AUDIT-010` | resistance target provenance | **confirmed evidence limitation** | preserve corrected scoring logic but do not present placeholders or unweighted benchmarks as empirical cell estimates |
| `AMR-AUDIT-011` | activity output timing | **contract decision required** | decide whether one row may intentionally combine pre-HGT stored activity with post-HGT recomputed activity planes |

No population-level bias or magnitude is claimed for the open items. The
forcing paths establish that the behavior can occur; representative and
calibrated runs are still needed to measure its practical effect.

## AMR-AUDIT-001 — model behavior is coupled to bacterium slot 32

**Status:** confirmed active defect.

### Observed behavior

The audited bacterium list has 42 entries. Slot 32 is
`enterobacter_cloacae`, while `treponema_pallidum` is slot 35 and
`helicobacter_pylori` is slot 37
(`src/simulation/population.rs:566`). Nevertheless, three output-collection
loops skip `b_idx == 32` (`src/simulation/simulation.rs:5333`, `:5379`, and
`:5394`). The skip affects:

- per-bacterium microbiome presence and regional presence;
- minority/majority resistance classification;
- resistant microbiome presence;
- carriage-duration bins; and
- microbiome acquisitions and clearances split by drug exposure.

As a result, Enterobacter cloacae is suppressed from those cells solely
because of its current array position.

### Why this is definitely defective

Adding or reordering an unrelated bacterium changes reported biology without
changing any named scientific policy. Repository history reproduces the drift:

| Historical checkpoint | Organism in slot 32 | Relevant event |
| --- | --- | --- |
| `209e884` | H. pylori | first related ordinal exclusion, described as an H. pylori exclusion |
| `ed8b3912` | H. pylori | literal output exclusion introduced |
| `9233340` | Treponema | earlier list insertions changed the excluded organism |
| `394a1df6` | Treponema | some related sites became name-based, while three literals remained |
| `0db03990` | E. cloacae | further insertions changed the remaining literal again |
| `54513895` | E. cloacae | three literal skips still present |

The defect proven here is the ordinal coupling. The history alone does not
prove which named species, if any, each field should exclude.

### Required repair process

1. Inventory each affected field separately. Do not assume presence,
   acquisition, clearance, resistance, regional, and duration outputs all
   require the same mask.
2. Resolve the intended policy from the field definition and scientific
   design. The code currently supplies conflicting clues:
   `src/rules/mod.rs:5270` prevents H. pylori microbiome generation, while
   `is_microbiome_excluded` at `src/simulation/simulation.rs:1178` names
   Treponema.
3. Implement the accepted policy by bacterium identity or a typed capability
   property, never an array ordinal.
4. Document any historical-output discontinuity.

### Mandatory guards

- One forcing fixture with E. cloacae, Treponema, and H. pylori legs for every
  affected output family.
- A list-reordering test, or a direct test of a name/capability-based policy,
  proving that unrelated insertions cannot change the decision.
- A mutation probe showing that changing one intended species decision fails
  the test.

## AMR-AUDIT-002 — comments, names, and executable species policy disagree

**Status:** confirmed source contradiction; intended policy unresolved.

### Observed behavior

`is_microbiome_excluded` returns true only for
`treponema_pallidum` (`src/simulation/simulation.rs:1178`). The helper is used
by clinical summary logic whose comments repeatedly say that H. pylori is
being excluded and still describe it as slot 32
(`src/simulation/simulation.rs:5308`, `:5460`, and `:5511`). In the audited
list, H. pylori is actually slot 37.

This affects at least:

- diagnostic-cascade entry eligibility;
- the pre-rule microbiome/resistance scan trigger;
- `new_drug_initiations_count_infected`;
- `total_currently_infected`;
- `currently_infected_and_on_drug_count`;
- the greater-than-10-day and greater-than-21-day infection-duration counts;
  and
- `infected_on_drug_with_previous_failure`.

All six executable call sites must be audited, not only the three carrying
stale H. pylori comments: `src/simulation/simulation.rs:915`, `:4597`, `:5313`,
`:5461`, `:5512`, and `:5831`. The pre-rule scan has an especially subtle
effect: a person whose only carriage is Treponema can fail the outer
`has_any_microbiome` trigger, while adding another included carried species can
cause the inner scan to run and expose the Treponema state. Reporting should
not depend on an unrelated second species unless that dependency is an
explicit contract.

The rules layer separately sets `allows_microbiome = false` for H. pylori at
`src/rules/mod.rs:5270`. Therefore the current source contains at least three
different apparent policies: a Treponema-named helper, H. pylori comments, and
an H. pylori microbiome-generation restriction.

### Why a comment-only edit is not enough

The history of slot 32 suggests that some name-based Treponema exclusions may
have crystallized behavior after list insertions had already moved the
original ordinal away from H. pylori. That is strong evidence for a semantic
review, but not enough to declare that every affected field must exclude
H. pylori. Treponema may also need special treatment for some microbiome or
clinical fields.

### Required decision and guards

For each affected field, state explicitly whether H. pylori, Treponema, both,
or neither is included, and why. Then:

- give policy helpers names that describe the decision rather than a generic
  `is_microbiome_excluded` label;
- align comments and model documentation with executable behavior;
- test H. pylori and Treponema independently; and
- search the repository for stale ordinal/species comments and literal
  bacterium indices before closing the issue.

## AMR-AUDIT-003 — “new infection” outputs are attrited post-rule stocks

**Status:** contract decision required.

### Observed behavior

`newly_infected_count` and `new_active_infections_by_bacteria` are collected
from people who are alive and whose infection still satisfies
`level > INFECTION_EPS` after the day's rules. A record contributes when its
surviving `date_last_infected` equals the current day
(`src/simulation/simulation.rs:5459`, `:5533`, and `:5729`).

A successful acquisition sets `date_last_infected` and
`clearance_ready_day` to the current day (`src/rules/mod.rs:5893-5895`). The
same day's clearance path is therefore reachable and can reset the infection
and its date before collection. Eligibility/duration is evaluated at
`src/rules/mod.rs:6329-6333`; immune and drug-assisted forcing paths continue
at `:6344-6346` and `:6484-6487`; the final state reset is at `:6555-6572`.

The affected surface is larger than the two headline counters. Every counter
inside the `date_last_infected == t` branch inherits the same final-state
predicate: syndrome and region incidence; per-bacterium active incidence;
carrier/non-carrier and pre-acquisition-risk splits; any-resistance carrier
splits; age and hospital/community acquisition splits; serious-resistance
eligibility/status; new-infection resistance incidence; and the derived
past-year scalar. The repair agent must inventory the exact current fields
rather than changing only `newly_infected_count`.

### Minimal forcing case

1. Force a successful acquisition for one person and bacterium on day `t`.
2. Force clearance later on day `t`.
3. Observe that neither the scalar nor per-bacterium “new infection” output
   records the acquisition.

This is not necessarily wrong, but it means the outputs are “new active
infections still observable at the end of the day,” not pure acquisition-event
incidence. A countervailing path is also reachable: a persisting sub-threshold
record can be re-dated by acquisition and then counted. The field names and
comments can reasonably be read as event incidence, but source inspection
alone does not establish the net direction or magnitude of bias.

### Required decision

Choose and document one contract:

1. **Acquisition-event incidence:** count every successful acquisition once,
   even if it clears later that day. A person-level scalar should be explicitly
   deduplicated while bacterium detail remains per bacterium.
2. **End-of-day observed new infection:** retain the final-state predicate and
   rename/document the fields so users cannot treat them as event incidence.
3. **Both:** expose separate, unambiguous event and end-of-day fields.

Do not infer the choice from historical equality. Add same-day clearance,
multiple-bacterium acquisition, threshold-boundary, and re-acquisition forcing
tests. If event incidence is selected, use an explicit daily event journal or
rule-event result rather than reconstructing events from mutable final state.

## AMR-AUDIT-004 — sepsis outputs mix episode, onset, and final-state semantics

**Status:** contract decision required.

### Observed behavior

The scalar `new_sepsis_cases` is limited to one person-level case when no
active sepsis existed before the rules and at least one infection is still
septic, above `INFECTION_EPS`, and onset-dated today after the rules
(`src/simulation/simulation.rs:818-831` and `:4784-4791`). Per-bacterium
counts use final sepsis state, final infection level, and onset day in separate
living and died-today paths (`src/simulation/simulation.rs:5205-5213` and
`:5751-5760`).

Sepsis can start earlier in the day (`src/rules/mod.rs:2837`) and be cleared by
infection clearance before collection. Such an onset can disappear from the
incidence outputs. Same-day sepsis recovery is also possible when configuration
permits it, but the audited default `sepsis_minimum_duration_days = 1` prevents
that particular default-path case. Conversely, simultaneous onsets can produce
one scalar case and multiple bacterium cases, which may be correct but must be
explicit.

### Required decision

Define these independently:

- Does the scalar count a new person-level sepsis episode, a person with any
  onset that day, or an end-of-day observed case?
- Does bacterium detail count every onset, only onsets still active at the end
  of the day, or only bacteria assigned as the primary cause?
- How are simultaneous onsets, a new onset during an existing episode,
  same-day recovery, below-threshold infection, and death handled?

`new_sepsis_cases_by_bacteria` and `sepsis_onset_count_by_bacteria` are
incremented together in both the living and died-today branches. Decide and
document whether they are intentional aliases for compatibility, two fields
that should diverge under a defined case, or accidental duplicate semantics.

Required tests must force every case above. If onset incidence is intended,
record onset events at transition time and keep person-level episode
deduplication separate from bacterium-level multiplicity.

## AMR-AUDIT-005 — toxicity-stop count excludes people who die later that day

**Status:** contract decision required.

### Observed behavior

Sub-lethal toxicity can stop a drug and stamp
`toxicity_stopped_drug_day` (`src/rules/mod.rs:4808-4835`). Mortality is
evaluated later in the same person-day, using the already computed toxicity
hazard (`src/rules/mod.rs:4903-5157`). The output count is then collected only
inside the living-person branch (`src/simulation/simulation.rs:5176` and
`:5319-5323`).

Therefore a real toxicity stop followed by death on the same day is omitted
from `drug_stops_due_to_toxicity`. The path is reachable because stopping the
drug does not erase the toxicity death probability already computed for that
day.

### Required decision

- If the field means **drug discontinuation events**, count the stop even when
  the person later dies.
- If it means **toxicity stops among end-of-day survivors**, retain the filter
  but rename and document it explicitly.
- If both are useful, publish separate fields.

Also define multiplicity: the current collector iterates drugs, so its natural
unit is stamped drug courses rather than people. The current rule selects and
stops only the single worst active drug per person-day, even when several drugs
are active. Tests must prove that current one-stop behavior, then separately
define whether output multiplicity is courses or people if the stopping rule
ever changes.

Use a small daily stop journal or a rule-event result. Do not increment a
global count and later retract it: that makes correctness depend on stage
order and invites underflow or double-adjustment bugs. Force stop-only,
stop-then-death, death-without-stop, and multiple-drug cases.

## AMR-AUDIT-006 — infection-associated drug starts use a different threshold

**Status:** contract decision required.

### Observed behavior

When a person starts any drug, `new_drug_initiations_count_infected` uses
`level > 0.0` (`src/simulation/simulation.rs:5308-5314`). Most neighboring
clinical infection summaries use `level > INFECTION_EPS`. The same expression
also applies the disputed species helper from `AMR-AUDIT-002`.

Consequently a positive residual infection below `INFECTION_EPS` can classify
a drug initiation as infection-associated even though the person is not
counted as currently infected by the nearby active-infection fields.
`infected_on_drug_with_previous_failure` repeats the same `level > 0.0`
predicate at `src/simulation/simulation.rs:5831` and belongs in this threshold
decision as well as the species-policy audit.

### Required decision and tests

Decide whether the field means any residual modeled infection or a clinically
active infection. Test exact values `0`, just above `0`, exactly
`INFECTION_EPS`, and just above `INFECTION_EPS`, with H. pylori and Treponema
as separate species-policy legs. Align the field description, code predicate,
and downstream analysis. Do not change the threshold globally without an
inventory of every field that intentionally uses a different definition.

## AMR-AUDIT-007 — processing after death is already terminal

**Status:** closed regression item in the audited source.

An earlier implementation could record a death and continue through later
infection acquisition, progression, clearance, and attribution stages for the
same person-day. The audited source records death and returns before later
rules at `src/rules/mod.rs:5156-5192`.

Do not reintroduce post-mortem model transitions to preserve old output. Keep
and strengthen the existing
`newly_recorded_death_stops_remaining_person_day_rules` guard, rather than
creating a weaker duplicate. It must continue to prove:

- no later state mutation occurs;
- no later random draws are consumed;
- infection state at the moment of death remains available for attribution;
  and
- output collection includes the intended died-today fields without treating
  the person as a living stock.

## AMR-AUDIT-008 — tracked documentation and configuration copies are stale

**Status:** confirmed active documentation/tooling risk; runtime impact depends
on whether any external workflow consumes the stale files.

### Observed contradictions

- `README.md:3` claims 58 antibiotics and 35 resistance mechanisms. The
  executable lists contain 62 drugs and 46 mechanisms.
- `src/simulation/population.rs:1559` says per-bacterium arrays have
  `BACTERIA_COUNT = 39`; the executable list contains 42 bacteria.
- The tracked `src/config.rs.bak` contains retired mechanism names such as
  `global_porin_loss` and multiple `as_yet_unknown_*` families.
- Tracked files under `calibration_configs/` also contain retired
  `global_porin_loss` and `as_yet_unknown` keys.

Even when these files are not compiled, they can mislead contributors,
automated extraction, calibration tooling, or an agent asked to infer the
model's dimensions and parameter vocabulary.

### Required repair process

1. Establish and document the single executable authority for bacteria, drug,
   mechanism, and parameter inventories.
2. Generate user-facing counts from that authority, or add a test that fails
   when documented counts drift.
3. Determine whether each backup/calibration file is live input, reproducible
   output, or historical evidence.
4. Regenerate live inputs. Move historical evidence to a clearly labelled,
   source-stamped archive that loaders and search-based generators cannot
   mistake for active configuration. Remove redundant backups only through a
   reviewed, recoverable repository change.
5. Add a repository-wide retired-key test covering active and executable-
   adjacent configuration surfaces.

## AMR-AUDIT-009 — infection-death fallback is only an approximation

**Status:** confirmed approximation in an analysis compatibility path.

### Observed behavior

When the dedicated `deaths_sepsis_model_scope` and
`deaths_infection_non_sepsis_model_scope` columns exist, analysis uses them
directly. When they are absent, `amr_simulation_output_analysis/
calibration_summary.py:1107-1135` falls back to broad infection-death totals
minus the sum of per-bacterium deaths for excluded organisms.

That subtraction is not person-level attribution. With co-infections, an
excluded organism can be present in a death that still has an in-scope
contributor, and per-bacterium death detail can be multiplicative. The fallback
can therefore subtract the wrong number of people. Its source comment already
acknowledges concurrent excluded infections; the remaining risk is that a
consumer may receive the approximate value without treating it as approximate.

### Required repair and tests

- Prefer the dedicated model-scope counters and fail clearly when an analysis
  requires exact scope but they are unavailable.
- If old output must remain readable, require explicit approximate-mode opt-in
  or emit machine-readable metadata and a prominent warning.
- Force deaths with only included infection, only excluded infection, and
  included/excluded co-infection; include multiple bacteria on the same person.
- Never use this fallback as exact calibration or validation evidence.

## AMR-AUDIT-010 — resistance targets have unresolved evidence provenance

**Status:** confirmed evidence limitation. Several scoring errors are repaired,
but the target values do not all have empirical cell-level support.

### Corrections that must remain fixed

Source history records these analysis repairs:

- `51e7597`: use configured resistance weights, whose default is 4:1, instead
  of a hard-coded 3:1 ratio;
- `ce2c3e9`: stop reusing any-resistance setting anchors for serious
  resistance;
- `ff7ab93`: reject malformed target-matrix rows;
- `9e94bf6` and `9045fce`: version the target schema and make score eligibility
  explicit;
- `e11fd70` and `683ed95`: exclude phenotypes the model cannot represent and
  severity targets above the attainable model maximum; and
- `3daafed`: stop copying prevalence values into reserve-drug severity targets.

### Unresolved evidence limitation

`data/RESISTANCE_TARGETS.md` states that cell-level provenance has not been
recovered and `evidence_weight` is blank because evidence quality has not been
assessed. It also identifies 48 reserve-drug cells as coarse expert
placeholders and five positive-severity/zero-prevalence cells as rare-positive
structural priors. These can be useful modeling constraints, but they must not
be described as empirical cell estimates or silently receive evidence-derived
confidence.

### Required protection

- Preserve source and rationale identity per cell, and distinguish empirical
  estimates, expert placeholders, and structural priors in every report.
- Do not invent evidence weights. Add them only through a documented review.
- Regenerate potency, resistance reachability, attainable-severity, target,
  and eligibility artifacts whenever mechanism identity, host scope, drug
  route, effect magnitude, or potency changes.
- Require the long-form targets, source table, schema, and generated matrices
  to agree exactly, with hashes and negative malformed-row tests.
- Make the analysis report counts and total score weight by provenance class
  so placeholders cannot dominate invisibly.

## AMR-AUDIT-011 — one activity row mixes within-day observation points

**Status:** contract decision required; the mixed timing is observed behavior,
not yet classified as a scientific defect.

### Observed behavior

For an active infection/drug pair, stored `activity_r` is computed during the
infection rule path (`src/rules/mod.rs:6179-6183`). Later HGT can change the
infection's resistance state without recomputing that stored activity. At
final collection, `activity_r_sum_by_bacteria` reads the stored value, while
the maximum-possible and pure activity planes are recomputed from the later
infection, exposure, potency, penetration, and resistance state
(`src/simulation/simulation.rs:5506-5508` and `:5631-5647`).

There is a second asymmetry: an infection acquired after the activity stage
has no stored activity contribution for that day, but can contribute to the
recomputed final planes if it survives above the infection threshold and has
eligible exposure. Thus one output row can combine a pre-HGT stored numerator
with post-HGT or same-day-final denominators/pure metrics.

### Required decision and tests

Decide whether each field is intended to describe activity actually applied
during the day's progression, final counterfactual activity, or one coherent
snapshot. Then force:

- HGT that changes resistance after stored activity is computed;
- same-day acquisition after the activity stage;
- same-day clearance;
- exposure start/stop around the activity stage; and
- exact `INFECTION_EPS` boundary cases.

If mixed timing is retained intentionally, document it field by field so users
do not divide or compare planes as though they share one observation point.
If one coherent snapshot is selected, account explicitly for how changing the
timing affects dynamics versus reporting only.

## Additional correctness-sensitive surfaces

The following source-history findings are not additional proven live defects.
They are recent corrections or unusual contracts that deserve explicit
regression tests and should be re-audited if their code changes.

| Surface | Current understanding | Required protection |
| --- | --- | --- |
| day-zero event dates (`d4f9a95`) | Missing dates now use `MISSING_EVENT_DATE = -1` at `src/simulation/population.rs:68`, allowing day 0 to remain a real event date | force initialization, acquisition, clearance, and every day-zero output gate; reject code that again treats `0` as “missing” |
| person-level sepsis multiplicity (`9afdccc`) | the scalar was changed to count an incident person once; bacterium detail remains separately multiplicative | retain simultaneous-onset, existing-episode/new-bacterium, death-day, and same-day-clearance cases under `AMR-AUDIT-004` |
| infection-death scope (`4b398cd`) | dedicated model-scope death counters were added alongside broad cause counters | force excluded organisms, threshold boundaries, multi-infection attribution, and background/toxicity causes |
| resistance output care setting (`d3a53a7`) | selected hospital/community resistance fields now use current care location rather than infection-acquisition route | test a person who changes care setting after acquisition and assert each field's named context independently |
| HGT donor timing (`908338e`) | donors are intended to be snapshotted before daily transfer so a new recipient cannot retransmit on the same day | force a three-person transfer chain and prove next-day, not same-day, second-hop eligibility |
| resistance vocabulary, hosts, routes, and potency | recent changes added/removed mechanisms, tightened exact-host rules, and corrected drug spectra | fingerprint the ordered mechanism list, every host decision, route/effect matrices, and all 42×62 potency cells from executable state |

### Resistance-science correction inventory

Source history records the following concrete resistance errors as repaired in
the audited commit. Treat them as regression requirements, not as instructions
to reapply a patch blindly:

| Commit | Corrected failure mode | Minimum regression evidence |
| --- | --- | --- |
| `b1edef9` | incomplete or misspelled bacterium×mechanism emergence grid | exact 42×46 = 1,932-key inventory with missing, duplicate, and unknown-key rejection |
| `5217974` | host eligibility was not enforced consistently across emergence, imported profiles, inheritance, floors, HGT, and resistance projection | exhaustive 42×46 host-status matrix plus production-path tests at every consumption boundary |
| `e8f3e0c` | narrow-spectrum Gram-negative penicillinase identity and host mappings were conflated | reviewed exact-host list and distinct route/effect tests |
| `ee67cdd` | H. pylori tetracycline resistance used proxy mechanisms rather than its organism-specific 16S route | H. pylori-only 16S tetracycline mechanism and negative host controls |
| `ad5ce6c` | a generic porin-loss mechanism created an over-broad phenotype | retired identity/key rejection and exact named porin mechanisms only |
| `5cad4a4` | OmpK35/36 host/drug mapping was too broad | K. pneumoniae-only host test and complete reviewed drug-route mask |
| `df1b20e` | OprD host/drug mapping was too broad | P. aeruginosa-only host test and exact carbapenem routes |
| `43c18c7` | primary GyrA mutation omitted nalidixic-acid projection | positive nalidixic-acid route plus negative unrelated-drug controls |
| `56934e1` | environmental-floor routes existed that executable reachability analysis showed could never act | prove every retained floor has a live eligible route, or remove it |
| `bd88ab0` | negligible-potency drugs inflated the multidrug de-novo resistance penalty | inactive, below-threshold, exactly-threshold, and above-threshold drug cases |
| `908338e` | same-day HGT recipients could become donors and amplify transfer within one day | immutable pre-transfer donor snapshot and next-day publication test |
| `7518f5f` | a duplicate hospital cache guard forced retention of a resistant profile | empty/non-resistant cache cases proving no artificial resistant profile is manufactured |
| `8beadd2` | synthetic S. aureus lineage completion injected an unsupported MRSA mechanism package and random draws | absence test for the retired package and deterministic draw-cadence guard |
| `64d921a` | ratchet memory was global rather than region×bacterium×mechanism | two-region divergence, isolation, merge, and checkpoint round-trip tests |
| `ca9045c`, `1f00053` | redundant or saturated hospital amplification controls remained in configuration and machinery | retired-key rejection and proof that active hospital effects have one documented owner |
| `b7872fc` | cefiderocol siderophore-uptake resistance lacked a named active route | exact host eligibility, cefiderocol-only projection, HGT status, and mechanism identity tests |
| `9120092` | OmpK alone incorrectly supplied a ceftolozane-tazobactam resistance route | negative OmpK-only test and positive supported-route tests |
| `0fe8819`, `5451389` | ceftolozane-tazobactam and cefiderocol potency spectra contained incorrect organism cells | exhaustive 42×62 = 2,604-cell potency fingerprint, not a test of only the changed cells |

The local resistance persistence archive and bounded pseudo-reservoir added by
`daf5acb` and `b3d3143` are not classified here as biological defects or
confirmed repairs. They are new, high-impact design mechanisms. Validate
profile deduplication, virtual-mass bounds, region/care isolation, sampling
law and draw cadence, output counters, and checkpoint state before treating
their trajectories as scientifically established.

The resistance acceptance matrix should independently cover 42 bacteria, 46
ordered mechanisms, 62 ordered drugs, all 1,932 emergence entries, all host
decisions, boolean drug routes separately from effect magnitudes, and all
2,604 potency cells. At minimum, include exact-host anchors for OmpK35/36,
OprD, LiaFSR, BlaZ, the reviewed narrow-spectrum Gram-negative penicillinase
list, H. pylori 16S tetracycline resistance, and cefiderocol siderophore
uptake. Each fingerprint must specify ordering and encoding and must fail when
one source cell, name, bit, route, effect, or key is deliberately mutated.

## Cross-cutting repair rules

1. **Freeze the source identity.** Report the exact starting and final commits.
2. **Separate defects from policy choices.** A reachable surprising result is
   not proof of the correct alternative.
3. **Prefer explicit events for incidence.** Do not reconstruct event counts
   from mutable end-of-day state when later transitions can erase them.
4. **Prefer named capabilities to ordinals.** Species lists may grow or be
   reordered without changing biology.
5. **Specify multiplicity.** Every scalar and detail output must say whether it
   counts people, episodes, bacteria, infections, courses, or events.
6. **Specify observation time.** State whether each field observes transition
   time, immediately after a stage, end of person-day, or end of model-day.
7. **Protect output users.** If semantics change, update schema descriptions,
   analysis code, fixtures, and release notes together.
8. **Use mutation-sensitive tests.** A green test that only repeats the
   production formula is insufficient; perturb species, thresholds, event
   order, and multiplicity to prove the guard can fail.
9. **Measure impact only after correctness.** Quantify changed rows and
   population-level effects with fixed seeds and representative scenarios;
   do not invent a bias estimate from source inspection.

## Expected agent deliverable

The working agent should return one row per issue containing:

| Field | Required content |
| --- | --- |
| Source identity | starting commit, final commit, and whether cited behavior was still present |
| Adjudication | confirmed defect, intentional contract, ambiguous, disproved, or already repaired |
| Decision evidence | code, history, documentation, scientific rationale, and contradictory evidence |
| Implemented change | exact logic/comments/schema changed, or `none` if a decision is still required |
| Tests | forcing cases, mutation probe, commands, and complete results |
| Output consequence | fields changed, compatibility/versioning treatment, and measured impact if run |
| Residual risk | unresolved intent, untested scenario, or downstream consumer requiring review |

The agent must not report an issue as fixed merely because compilation passes
or a current-output fixture remains unchanged.
