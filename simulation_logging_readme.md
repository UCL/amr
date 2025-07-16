## Simulation Run Log

This file demonstrates the ongoing log functionality that tracks each simulation run.

### Features:
- **Timestamp**: UTC timestamp of when the simulation was run
- **Population Size**: Number of individuals in the simulation
- **Time Steps**: Number of time steps the simulation ran for
- **Duration**: Total execution time in seconds

### Log Format:
```csv
timestamp,population_size,time_steps,duration_seconds
```

### Example Usage:
Every time you run the simulation executable (`./target/release/executable_amr.exe`), a new line will be automatically appended to `simulation_run_log.csv` with the run details.

### Performance Tracking:
You can use this data to:
1. Monitor simulation performance over time
2. Compare execution times for different parameter settings
3. Track your experimentation history
4. Identify performance patterns or regressions

The log file will grow over time and never overwrite previous runs, giving you a complete history of your simulation experiments.
