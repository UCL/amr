# Local population-3,000 comparison

This folder contains two optimized Windows executables built from the same workspace state:

- `bin/pre_refactor_pop3000.exe`: current committed baseline (`a3e463d`).
- `bin/post_refactor_pop3000.exe`: same baseline with only the two
  pre-acquisition `Vec<bool>` values replaced by `u64` bitmasks.

Both use population 3,000 and `CalibrationMode::Full`. The workspace launcher has been restored
to population 10,000,000; changing it later does not change these already-built executables.

Binary SHA-256 values:

- Pre-refactor: `96AA3341AA8AB322FD34E5A6D6B8957086BEF89E4548D594F2FE7C4D9DA8960F`
- Post-refactor: `AB5F8452FDAAD472ECE454312E976F70AABAA2DF5F68FA65775BE4AC5E4180AF`

Run one fixed-seed comparison from PowerShell:

```powershell
.\local_perf_comparison_3000\run_comparison.ps1 -Threads 8 -Repetitions 1
```

Replace `8` with the thread count you want to hold constant. The runner defaults to seed
`123456789`, creates separate timestamped output directories, records model and wall duration,
and compares the summary CSV files by SHA-256. It fails at the end if the summaries differ.

For timing after the first correctness run, use three alternating repetitions:

```powershell
.\local_perf_comparison_3000\run_comparison.ps1 -Threads 8 -Repetitions 3
```

Each batch writes `comparison_<timestamp>.csv` in this folder. Outputs and console logs are under
`runs/pre_refactor` and `runs/post_refactor`.
