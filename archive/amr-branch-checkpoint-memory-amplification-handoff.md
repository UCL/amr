# AMR Branch Checkpoint Memory Amplification Handoff

Status: investigation handoff; no model fix has been implemented or approved

Date: 2026-08-07

Primary repository: `C:\Users\MAUG\Documents\DEV\amr-main`

Upstream: `https://github.com/UCL/amr`

## Purpose

This document gives another engineer or coding agent enough evidence to independently investigate and repair the AMR model's branch-checkpoint memory amplification. It deliberately does not put AMR behavior into the job-platform coordinator, runner, or WebUI.

The immediate incident was not evidence that a 10-million-person counterfactual simulation inherently requires approximately 690 GiB. The model normally held approximately 365-374 GiB during the baseline portion. Memory then rose to approximately 690 GiB over the final few minutes, at the point where the model was expected to capture counterfactual branch state. The source constructs the supposedly disk-backed checkpoint by cloning the full in-memory population before serializing it. That is an avoidable implementation cost.

## Executive conclusion

There are two separate facts:

1. Branch-enabled and full-history modes legitimately require more state than reduced calibration modes. The population, mechanism cache, retained summary rows, and branch checkpoint all have real costs.
2. The current branch-checkpoint implementation multiplies that cost unnecessarily. Disk-backed checkpointing first creates an owned snapshot by cloning the population, mechanism cache, and summary log. Later restoration can hold the completed baseline state, a deserialized checkpoint, and another clone used to restore the simulation. This is not required by the mathematical model or by counterfactual semantics.

The observed late increase of roughly 315 GiB is consistent with cloning a population already occupying hundreds of GiB. It is not consistent with ordinary gradual growth in the simulation.

## Incident identity

| Field | Value |
| --- | --- |
| Submission | `submission-189a61a8d13e469a` |
| Job | `job-65d29c50c38d4957` |
| Attempt | `bootstrap-attempt::submission-189a61a8d13e469a::f9d7e34d-d2b2-42c3-a133-5ae36078ebf9` |
| Run group | `run-group-4bdf1980191b4327` |
| Runner | `runner-node-vm240` |
| Project | `UCL/amr latest` |
| Source commit | `c281cbd95c24718791fd5de0bfcee630b5dc7149` |
| Population | 10,000,000 |
| Time steps | 35,040 |
| Mode | `partial_25_counterfactual` |
| Seed | `5965993106776908880` |
| CPU | 72 Rayon threads; CPUs 0-71 |
| NUMA memory nodes | 0-2, interleaved/balanced |
| Start | `2026-08-06T15:20:16.957Z` |
| Terminal time | `2026-08-07T04:24:54.420Z` |
| Result | exit 137; `systemd` result `oom-kill` |
| Last stdout progress | approximately timestep 33,500 |

Commit `c281cbd` changed the selected mode from `Full25Counterfactual` to `Partial25Counterfactual` in `src/main.rs`. A later upstream commit changed it back. Review the exact change before assuming the mode name still exists on current `main`:

`https://github.com/UCL/amr/commit/c281cbd95c24718791fd5de0bfcee630b5dc7149`

## Observed memory behavior

Durable runner observability captured approximately 47,053 samples over 47,050 seconds.

- Memory stayed mostly between approximately 365 and 374 GiB for nearly 13 hours.
- During the last approximately four minutes it rose through roughly 425, 506, 587, and 668 GiB.
- Recorded peak was approximately 690.05 GiB.
- The three NUMA nodes remained reasonably balanced at roughly 235.5, 241.7, and 246.6 GB near failure.
- CPU placement remained CPUs 0-71 and memory placement remained nodes 0-2.
- Guest `MemAvailable` reached zero.
- The kernel reported a guest-global OOM and selected `executable_amr` for termination.

The shape matters: this was a late, steep allocation event rather than baseline state slowly increasing throughout the run. NUMA imbalance did not cause the memory requirement.

Prior runs provide a useful comparison, but should not be treated as a formal benchmark because source revisions and machine contention must be controlled:

| Run | Approximate peak |
| --- | ---: |
| 10M Full | 365.8 GiB |
| 3M Full | 110.3 GiB |
| Failed 10M partial counterfactual | 690.0 GiB |

## Exact amplification path

### 1. The disk path clones before it checks the disk setting

Current `main` performs this sequence in `src/simulation/simulation.rs` around `run_from()`:

```rust
let snapshot = self.create_branch_snapshot();
if self.use_disk_branch_checkpoint {
    let path = self.persist_branch_snapshot_to_disk(&snapshot, step)?;
    branch_snapshot = Some(StoredBranchSnapshot::OnDisk(path));
} else {
    branch_snapshot = Some(StoredBranchSnapshot::InMemory(snapshot));
}
```

The decision to persist to disk occurs only after an owned `BranchSnapshot` already exists.

### 2. Snapshot construction clones all major state

Current `create_branch_snapshot()` is:

```rust
fn create_branch_snapshot(&self) -> BranchSnapshot {
    BranchSnapshot {
        population: self.population.clone(),
        mechanism_cache: self.mechanism_cache.clone(),
        summary_log: self.summary_log.clone(),
    }
}
```

For a 10-million-person simulation, `self.population.clone()` is an O(population) allocation. The cache and complete retained summary history are cloned as well. Disk-backed mode therefore needs enough memory for the active simulation and the complete owned snapshot before the first checkpoint byte can be written.

At the failed revision, this first clone is the most likely immediate OOM point because failure aligned with branch capture and the run did not produce normal completion evidence.

### 3. Materialization can duplicate the completed baseline

After the baseline finishes, an on-disk snapshot is deserialized into a new owned `BranchSnapshot` while the simulation still owns its current end-of-baseline population, cache, and summary log:

```rust
let data = self.load_branch_snapshot_from_disk(&path)?;
```

This can again hold two full generations of large state at once.

### 4. Branch restoration clones the snapshot again

Current `run_policy_branch()` restores with:

```rust
self.population = snapshot.population.clone();
self.mechanism_cache = snapshot.mechanism_cache.clone();
self.summary_log = snapshot.summary_log.clone();
```

The owned snapshot remains alive while another owned copy is created for `self`. Depending on allocation and drop order, restoration can temporarily retain the completed baseline state, the deserialized checkpoint, and the new restoration clone.

### 5. Multiple branches repeat restoration

The model loops through `branch_policy_adjustments` and calls `run_policy_branch(&branch_snapshot, ...)` for each policy. Reusing an immutable owned snapshot is logically convenient, but cloning it for every branch repeats O(population) allocation and copy work.

## Why the mode itself still costs more

Not every difference from a reduced calibration run is a bug.

Current UCL/amr `main` documents these broad semantics in `src/simulation/simulation.rs`:

- `None` retains every row, enables policy branches, and collects policy diagnostics.
- `Partial` retains every baseline row but skips policy branches.
- `FullMinimal` retains calibration-window rows and lean fields.
- `Full` retains calibration-window rows and the complete calibration fields.

`Partial` and `None` use `SummaryContentFlags::all()`, while reduced calibration modes discard disabled field groups. Full-history modes therefore legitimately retain more summary data. A branch-enabled mode also needs a recoverable branch-point state somewhere while the baseline continues.

What is not inherent is keeping that recoverable state as another complete in-memory population when disk checkpointing was selected, or cloning the checkpoint again during restoration. A sequential baseline-plus-branch model needs one active population plus durable branch state, not multiple simultaneous owned populations.

## Current-main relevance

The exact `Partial25Counterfactual` enum from the failed revision is no longer present on local current `main` (`093a50872ab4b97af3824ee91778617b9b26dbcd`). The problem has not disappeared:

- `src/main.rs` currently enables disk branch checkpointing for `CalibrationMode::None`.
- `run_from()` still calls `create_branch_snapshot()` before checking `use_disk_branch_checkpoint`.
- `create_branch_snapshot()` still clones population, cache, and summary log.
- `materialize_branch_snapshot()` still deserializes into an owned snapshot while active simulation state exists.
- `run_policy_branch()` still clones all snapshot state into the simulation.

The affected user-facing mode changed, but the defective ownership and serialization design remains in the branch-enabled execution path.

## Required investigation by the next agent

The next agent must not begin by changing code. It should first produce evidence for each stage below.

### A. Reproduce at a safe scale

Use fixed-seed 300k and 1M runs that cross the branch point. Capture:

- cgroup `memory.current`, `memory.peak`, and `memory.events`;
- process RSS/PSS;
- checkpoint start/end timestamps and file size;
- state immediately before and after `create_branch_snapshot`;
- state before and after deserialization;
- state before and after assignment in `run_policy_branch`;
- per-NUMA memory distribution;
- output hashes and model-local run identity.

Do not begin with another 10M run.

### B. Prove allocation ownership

Use an appropriate Rust/Linux allocation profiler or narrowly scoped instrumentation to quantify allocations attributable to:

- `Population::clone`;
- `MechanismCache::clone`;
- `summary_log.clone`;
- `bincode::serialize_into`;
- `bincode::deserialize_from`;
- assignment/restoration in `run_policy_branch`.

The investigation must distinguish anonymous heap growth from page cache used by checkpoint I/O.

### C. Confirm serialization behavior

Prove whether every nested type can be serialized through borrowed references without constructing an owned `BranchSnapshot`. Check custom `Serialize` implementations for hidden whole-collection buffers.

### D. Confirm branch lifecycle

Document exactly which baseline outputs must remain after baseline completion, which state can be dropped before restore, and whether each policy branch can reload the same checkpoint independently. Do not assume state can be discarded until output ownership and failure behavior are understood.

## Preferred implementation direction

The design should preserve semantics while changing ownership:

1. In disk mode, serialize a borrowed view of the active population, mechanism cache, and required summary state directly to a temporary checkpoint file. Do not call the owned cloning constructor.
2. Flush, validate, and atomically publish the checkpoint. Include format version, source/config identity, expected size or checksum, and branch step.
3. Continue the baseline with the original active state after serialization.
4. Before restoring a branch, finalize or move out required baseline outputs and release the no-longer-needed active population/cache/summary state.
5. Deserialize checkpoint state into owned values and move those values into `Simulation`; do not clone them into `self`.
6. For multiple branches, finalize and release the previous branch, then reload the durable checkpoint for the next branch. This trades controlled disk I/O for bounded memory.
7. Keep the in-memory checkpoint path only when explicitly selected and safe. It must not silently remain the large-run default.

This direction must be reviewed against the actual output and policy-branch lifecycle before implementation. It is not permission to replace the code mechanically.

## Failure handling requirements

The repaired model must handle these cases explicitly:

- insufficient checkpoint disk space before branch capture;
- write failure and partial checkpoint files;
- checksum or format mismatch;
- cancellation during write or restore;
- process restart with stale checkpoint files;
- multiple policy branches;
- no configured branch or a branch step outside the run horizon;
- fixed and generated seeds;
- cleanup after success and failure.

Checkpoint writes must use temporary-file plus atomic-rename semantics. Cleanup must never delete an unrelated path.

## Required correctness tests

For a controlled fixed seed and identical configuration, old and repaired implementations must produce equivalent baseline and branch results. At minimum verify:

- branch starts from the same population state;
- mechanism-cache state is identical;
- required summary history is identical;
- policy adjustment starts at the same timestep;
- final summary artifact hashes match when deterministic serialization/output ordering permits;
- differences, if any, are explained at field level rather than accepted as numerical noise.

Also cover zero, one, and multiple policy branches and both disk and intentionally in-memory checkpoint modes.

## Required performance and memory gates

The repair is not complete merely because a small run exits successfully.

1. Run the same fixed-seed branch case at 300k, 1M, and 3M population.
2. Measure steady-state memory and checkpoint/restore peak separately.
3. Demonstrate that checkpoint memory delta no longer scales as approximately one additional full population.
4. Demonstrate that restore does not temporarily retain multiple full populations.
5. Record checkpoint size, write/read throughput, branch transition latency, and CPU overhead.
6. Run 10M only after smaller scaling results predict it fits with explicit system and job headroom.
7. A 10M acceptance run must complete without guest-global OOM and stay within its declared memory envelope.

No arbitrary percentage should be used to declare success. The scaling evidence must show that the O(population) duplicate ownership has been removed.

## Component boundary

This is an AMR model implementation issue.

- Model checkpoint ownership, branch semantics, calibration modes, and model output correctness belong in UCL/amr.
- The AMR workload adapter may declare generic resource requirements and expose generic progress/evidence.
- The job-platform runner may enforce generic memory and CPU limits.
- The coordinator may admit against generic allocatable capacity and adapter-declared requirements.
- The coordinator, runner, and WebUI must not implement AMR checkpoint logic or interpret AMR mode names.

## Delivery requirements

- Work on a dedicated UCL/amr branch and submit a pull request.
- Do not modify UCL/amr `main` directly.
- Keep commits separated into measurement/test infrastructure, ownership repair, failure handling, and documentation/cleanup.
- Do not merge until fixed-seed parity, scaling evidence, Rust quality gates, and a controlled installed job-platform proof all pass.
- Do not claim the job-platform is the model fix; it can only contain and report a model overrun.

## Local source references

At local UCL/amr `main` commit `093a50872ab4b97af3824ee91778617b9b26dbcd`:

- `src/main.rs:394-400`: selects disk checkpointing for branch-enabled `None` mode.
- `src/simulation/simulation.rs:2891-2943`: disk serialization/deserialization.
- `src/simulation/simulation.rs:2956-2973`: owned snapshot is created before disk-path selection.
- `src/simulation/simulation.rs:6385-6434`: baseline and policy-branch lifecycle.
- `src/simulation/simulation.rs:6451-6456`: full state clone.
- `src/simulation/simulation.rs:6459-6468`: owned disk materialization.
- `src/simulation/simulation.rs:6472-6488`: restoration clones snapshot state.

Line numbers are evidence for the recorded commit only and must be refreshed after upstream changes.
