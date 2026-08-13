[CmdletBinding()]
param(
    [ValidateRange(1, 512)]
    [int]$Threads = [Environment]::ProcessorCount,

    [ValidateRange(1, 20)]
    [int]$Repetitions = 1,

    [UInt64]$Seed = 123456789
)

$ErrorActionPreference = 'Stop'

$comparisonRoot = $PSScriptRoot
$runsRoot = Join-Path $comparisonRoot 'runs'
$executables = @{
    pre_refactor  = Join-Path $comparisonRoot 'bin\pre_refactor_pop3000.exe'
    post_refactor = Join-Path $comparisonRoot 'bin\post_refactor_pop3000.exe'
}

foreach ($executable in $executables.Values) {
    if (-not (Test-Path -LiteralPath $executable -PathType Leaf)) {
        throw "Missing comparison executable: $executable"
    }
}

$batchId = Get-Date -Format 'yyyyMMdd_HHmmss'
$previousSeed = $env:AMR_RNG_SEED
$previousThreads = $env:RAYON_NUM_THREADS
$results = [System.Collections.Generic.List[object]]::new()

try {
    $env:AMR_RNG_SEED = $Seed.ToString()
    $env:RAYON_NUM_THREADS = $Threads.ToString()

    Write-Host "Comparison batch: $batchId"
    Write-Host "Population: 3000"
    Write-Host "AMR_RNG_SEED: $env:AMR_RNG_SEED"
    Write-Host "RAYON_NUM_THREADS: $env:RAYON_NUM_THREADS"

    for ($repetition = 1; $repetition -le $Repetitions; $repetition++) {
        # Alternate execution order to reduce temperature/order bias across repeated runs.
        $order = if ($repetition % 2 -eq 1) {
            @('pre_refactor', 'post_refactor')
        } else {
            @('post_refactor', 'pre_refactor')
        }

        $runs = @{}
        foreach ($version in $order) {
            $runName = '{0}_rep{1:D2}' -f $batchId, $repetition
            $runDirectory = Join-Path (Join-Path $runsRoot $version) $runName
            New-Item -ItemType Directory -Path $runDirectory -ErrorAction Stop | Out-Null

            $logPath = Join-Path $runDirectory 'console.log'
            Write-Host ""
            Write-Host "Running $version repetition $repetition in $runDirectory"

            Push-Location $runDirectory
            try {
                $timer = [System.Diagnostics.Stopwatch]::StartNew()
                # Windows PowerShell 5 wraps native stderr as NativeCommandError records.
                # The model intentionally writes startup diagnostics to stderr, so tolerate
                # that stream here and use the native process exit code as the failure signal.
                $strictErrorActionPreference = $ErrorActionPreference
                $ErrorActionPreference = 'Continue'
                try {
                    & $executables[$version] 2>&1 |
                        ForEach-Object { $_.ToString() } |
                        Tee-Object -FilePath $logPath
                    $exitCode = $LASTEXITCODE
                } finally {
                    $ErrorActionPreference = $strictErrorActionPreference
                    $timer.Stop()
                }
            } finally {
                Pop-Location
            }

            if ($exitCode -ne 0) {
                throw "$version repetition $repetition exited with code $exitCode"
            }

            $outputDirectory = Join-Path $runDirectory 'amr_simulation_output_analysis_outputs'
            $summaryFiles = @(Get-ChildItem -LiteralPath $outputDirectory -Filter 'simulation_summary_*.csv' -File)
            if ($summaryFiles.Count -ne 1) {
                throw "Expected one summary CSV for $version repetition $repetition; found $($summaryFiles.Count)"
            }

            $runLogPath = Join-Path $runDirectory 'simulation_run_log.csv'
            $runLog = @(Import-Csv -LiteralPath $runLogPath)
            if ($runLog.Count -ne 1) {
                throw "Expected one simulation timing row for $version repetition $repetition"
            }

            $runs[$version] = [pscustomobject]@{
                Directory = $runDirectory
                Csv = $summaryFiles[0].FullName
                CsvHash = (Get-FileHash -Algorithm SHA256 -LiteralPath $summaryFiles[0].FullName).Hash
                ModelSeconds = [double]$runLog[0].duration_seconds
                WallSeconds = $timer.Elapsed.TotalSeconds
            }
        }

        $exactMatch = $runs.pre_refactor.CsvHash -eq $runs.post_refactor.CsvHash
        $results.Add([pscustomobject]@{
            repetition = $repetition
            seed = $Seed
            rayon_threads = $Threads
            pre_model_seconds = $runs.pre_refactor.ModelSeconds
            post_model_seconds = $runs.post_refactor.ModelSeconds
            pre_wall_seconds = [Math]::Round($runs.pre_refactor.WallSeconds, 3)
            post_wall_seconds = [Math]::Round($runs.post_refactor.WallSeconds, 3)
            pre_summary_sha256 = $runs.pre_refactor.CsvHash
            post_summary_sha256 = $runs.post_refactor.CsvHash
            exact_csv_match = $exactMatch
            pre_run_directory = $runs.pre_refactor.Directory
            post_run_directory = $runs.post_refactor.Directory
        })

        if ($exactMatch) {
            Write-Host "Exact summary match for repetition ${repetition}: $($runs.pre_refactor.CsvHash)" -ForegroundColor Green
        } else {
            Write-Warning "Summary mismatch for repetition $repetition"
        }
    }
} finally {
    $env:AMR_RNG_SEED = $previousSeed
    $env:RAYON_NUM_THREADS = $previousThreads
}

$reportPath = Join-Path $comparisonRoot "comparison_$batchId.csv"
$results | Export-Csv -LiteralPath $reportPath -NoTypeInformation
$results | Format-Table repetition, pre_model_seconds, post_model_seconds, exact_csv_match -AutoSize
Write-Host "Comparison report: $reportPath"

if (@($results | Where-Object { -not $_.exact_csv_match }).Count -gt 0) {
    throw 'At least one pre/post summary CSV pair differed'
}
