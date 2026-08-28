# Person-level targeted therapy can advance an unidentified bacterium-specific episode — CONFIRMED\_DEFECT

## Validation anchor

Confirmed at official UCL/amr main commit `aa1c3012b8319e11ca85a1b5817d9aa0fbc1c70d`, verified as an actual upstream commit. The frozen full-run comparison core `ed0243b479699a7dec29477997980653acd98de4` is also affected. The treatment-specific path at `2c59d16c5b2f809bb37e7d876223462440ad17cc` already contains the correct ordering predicate.

## Mechanism

The diagnostic cascade is bacterium-specific, but targeted antibiotic context is stored per person and drug.

Relevant current-upstream source:

* `src/simulation/population.rs`: `AntibioticUseContext::Targeted`; `Individual::test\_identified\_infection`; `Individual::drug\_use\_context`; and the episode guards `diagnostic\_cascade\_bacterial\_identification\_recorded`, `diagnostic\_cascade\_targeted\_treatment\_recorded`, and `diagnostic\_cascade\_effective\_targeted\_treatment\_recorded`.
* `src/rules/mod.rs`: `start\_drug\_course` and `apply\_rules` establish and retain person/drug-level treatment context.
* `src/simulation/simulation.rs`: `targeted\_drug\_started\_today\_for\_bacterium`, `best\_active\_targeted\_antibiotic\_activity`, `record\_diagnostic\_cascade\_stage`, and `Simulation::run\_from` implement the relevant accounting.

`targeted\_drug\_started\_today\_for\_bacterium` correctly requires identification of the supplied bacterium and a targeted course started today. By contrast, `best\_active\_targeted\_antibiotic\_activity` accepts any active person-level targeted course and evaluates its activity against the supplied bacterium without requiring that bacterium's identification stage.

Minimal trigger and control flow:

1. One person has infections A and B.
2. A is identified, so an active drug course has `AntibioticUseContext::Targeted`.
3. B has an eligible, open cascade episode, but `test\_identified\_infection\[B]` and `diagnostic\_cascade\_bacterial\_identification\_recorded\[B]` are false.
4. The active targeted-context drug has sufficient stored activity against B. A prior-day start also makes the ordinary same-day targeted-start predicate false.
5. In `Simulation::run\_from`, B skips bacterial-identification accounting and the ordinary targeted-start branch. The effective-treatment branch nevertheless backfills targeted treatment through `record\_diagnostic\_cascade\_stage`, sets `diagnostic\_cascade\_targeted\_treatment\_recorded\[B]`, then records effective targeted treatment and sets its guard.

The episode can therefore record targeted and effective treatment without ever recording bacterial identification.

## Observable consequence

`TimeStepSummary::diagnostic\_cascade\_stage\_counts` and `diagnostic\_cascade\_stage\_counts\_by\_setting` feed `Simulation::export\_summary\_to\_csv`. The corresponding analysis in `amr\_simulation\_output\_analysis/make\_paper\_tables.py` uses `\_SF5\_STAGES`, `\_sf5\_reliability\_flag`, and `make\_supplementary\_figure\_s5\_diagnostic\_testing\_targeted\_treatment\_cascade`.

Affected output can overcount `diagnostic\_cascade\_targeted\_treatment\_started` and `diagnostic\_cascade\_effective\_targeted\_treatment\_started` for bacterium-specific episodes. Reports can show downstream counts above identification, percentages above 100% of the preceding applicable stage, or `non-monotonic cascade count`. Aggregation may hide non-monotonicity while retaining contaminated downstream numerators.

## Execution scope

The defective accounting executes in testing-enabled Full and all/default routes. It is dormant in the common-output-only route because testing summary content is disabled there. The treatment-specific full and diagnostic-window routes at `2c59d16c5b2f809bb37e7d876223462440ad17cc` share the corrected accounting path.

## Smallest safe repair

Reuse the proven treatment-specific ordering pattern:

1. Permit effective-treatment recording only when the current episode's bacterial-identification and targeted-treatment guards are already true, its effective-treatment guard is false, and activity meets the existing threshold.
2. Remove stage-3 backfilling from the effective-treatment branch. Targeted treatment must be recorded only by the existing identification-gated, same-day targeted-start path.

The existing `diagnostic\_effective\_targeted\_recordable` predicate at the corrected commit expresses this condition and can be ported narrowly.

## Focused regression test

Add a same-module deterministic unit test using two infections. With A identified and an active targeted-context drug effective against B, make B eligible for a new episode but leave it unidentified, and give the course a prior-day start. After one accounting step, B must have eligible count `1`, identification count `0`, targeted count `0`, effective count `0`, with all three downstream episode guards false.

Then exercise a fresh B episode where identification is recorded and a targeted course starts that day. Assert exactly one write per stage in this order: bacterial identification, targeted treatment, effective targeted treatment. The episode invariant is:

`effective recorded => targeted recorded => bacterial identification recorded`.

## Compatibility caution

Keep the repair output-accounting-only. Preserve rule/event ordering, random draws, treatment selection, antibiotic-activity calculation, and infection and drug state. Legitimate same-day identification followed by targeted and effective treatment must remain countable in its existing order.

