# B20 eligible stage can be recorded for a person who died earlier

Status: confirmed defect at the exact source identity below; correction pending

## Exact source identity

The model source snapshot is commit
`ed0243b479699a7dec29477997980653acd98de4`, tree
`5ea2967350b9c89c4bd2e705e7cae15c7b4ca7d6`. The affected file is
`src/simulation/simulation.rs`, blob
`ab6a407a3d5c22b7aa26321ef95baf919a9c984e`.

The relevant functions and caller are all in that file:

- `person_day_vital_status` defines whether a person is available for
  current-day event attribution;
- `diagnostic_cascade_entry_eligible` checks infection, symptom, delay and
  testing conditions but has no death-status condition; and
- the post-rule diagnostic-cascade loop computes the vital status and
  `active_now`, opens a closed episode from
  `diagnostic_cascade_entry_eligible`, records
  `DIAGNOSTIC_CASCADE_ELIGIBLE_IDX`, then resets an open episode when
  `!active_now || died_today`.

## Observed and expected behavior

Observed: a person who died before the current day can retain a positive,
symptomatic, old-enough infection record. The eligibility predicate returns
true for that stale state. The caller records a new eligible stage and then
resets the episode on the same day because the person is not active. Persistent
stale state permits the same person-bacterium pair to re-enter on later days.

Expected: a person who died before the current day is unavailable for any new
current-day event attribution. A person whose death is recorded on the current
day remains available at the established same-day observation point. Apply
that exact boundary before opening or recording a diagnostic-cascade episode;
do not change the behavior for a living person.

## Minimal reproduction

At the named commit, add a test beside the existing
`diagnostic_cascade_entry_eligible` tests in
`src/simulation/simulation.rs`:

```rust
#[test]
fn earlier_dead_stale_symptomatic_state_reaches_eligibility_predicate() {
    let mut rng = SmallRng::seed_from_u64(92);
    let mut individual =
        Individual::new(1, 40 * 365, "female".to_string(), &mut rng);
    let bacteria_idx = BACTERIA_LIST
        .iter()
        .position(|&name| name == "escherichia_coli")
        .expect("E. coli must be modelled");
    let param_cache = ParameterKeyCache::new();
    let availability_day = param_cache.bacteria_test_availability_day[bacteria_idx]
        .unwrap_or(param_cache.bacterial_testing_available_from_day as usize);
    let current_day = availability_day.max(param_cache.test_delay_days as usize);

    individual.level[bacteria_idx] = 1.0;
    individual.infection_has_caused_symptoms[bacteria_idx] = true;
    individual.date_last_infected[bacteria_idx] = 0;
    individual.date_of_death = Some(current_day - 1);

    assert!(
        !person_day_vital_status(&individual, current_day)
            .available_for_current_day_event_attribution
    );
    assert!(diagnostic_cascade_entry_eligible(
        &individual,
        bacteria_idx,
        current_day,
        &param_cache,
    ));
}
```

Run only the focused locked release test:

```text
cargo test --locked --release \
  earlier_dead_stale_symptomatic_state_reaches_eligibility_predicate \
  -- --nocapture
```

The two assertions pass before correction. The second assertion is a witness
to the missing vital-status boundary; it is not the desired post-correction
contract by itself. Source inspection then completes the production-path
reproduction: the closed-episode caller invokes that predicate without the
vital guard, records the eligible stage, and resets the episode through
`active_now`.

## Correction and regression guidance

Gate closed-episode entry with
`person_day_vital_status(...).available_for_current_day_event_attribution`
before opening the episode or recording `DIAGNOSTIC_CASCADE_ELIGIBLE_IDX`.
Reusing the caller's already computed vital status avoids inventing a second
death rule. An equivalent pure helper is acceptable if it has the same
observation point and same-day-death semantics.

Add mutation-sensitive checks for all of these cases:

1. living, otherwise eligible state opens once and records one eligible stage;
2. death on the current day preserves the established current-day attribution;
3. death before the current day opens nothing and records nothing;
4. persistent stale state after earlier death cannot re-enter on a later day;
5. existing infection-age and test-availability failures remain ineligible.

Prefer a caller-level or extracted-pure-gate test that observes both episode
state and eligible-stage count. Keep the focused predicate witness if it helps
document that death is owned by the caller rather than the infection/test
predicate.

## Impact boundary

Known: any eligible-stage event created by this exact earlier-dead stale-state
path is invalid current-day attribution.

Unknown: runtime prevalence, aggregate output magnitude, effects on later
diagnostic stages, and persistence at any source other than the exact commit
named above. Do not invalidate unrelated B20 contributions. After correction,
identify the new commit and tree, run the focused tests, and quantify impact
with a separately governed comparison if magnitude is needed.
