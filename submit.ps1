# submit.ps1
$ErrorActionPreference = 'Continue'

try {
    Write-Host "[0] Preflight SSH connectivity..." -ForegroundColor Cyan
    $testResult = ssh -o ConnectTimeout=10 -o BatchMode=yes dev@runner "echo ping" 2>&1
    if ($LASTEXITCODE -ne 0) {
        throw "[0] SSH preflight failed (exit $LASTEXITCODE). Output: $testResult"
    }
    Write-Host "[0] SSH OK" -ForegroundColor Cyan

    # 1. Generate job ID
    Write-Host "[1] Generating job ID..." -ForegroundColor Cyan
    $JID = (Get-Date -Format "yyyyMMddTHHmmss") + "-" + (Get-Random -Maximum 99999)
    Write-Host "[1] JID = $JID" -ForegroundColor Cyan

    # 2. Use cargo metadata to find exactly what we need
    Write-Host "[2] Running cargo metadata..." -ForegroundColor Cyan
    $metadataJson = cargo metadata --no-deps --format-version 1 2>&1
    if ($LASTEXITCODE -ne 0) {
        throw "[2] cargo metadata failed (exit $LASTEXITCODE). Output: $metadataJson"
    }
    $metadata = $metadataJson | ConvertFrom-Json

    # Get the workspace root (where top-level Cargo.toml lives)
    $workspaceRoot = $metadata.workspace_root
    Write-Host "[2] Workspace root: $workspaceRoot" -ForegroundColor Cyan

    # 3. Create remote dir - use same flat structure as tar version
    $remoteDir = "/home/dev/inbox/$JID"
    Write-Host "[3] Creating remote dir $remoteDir..." -ForegroundColor Cyan
    ssh -o ConnectTimeout=10 -o BatchMode=yes dev@runner "mkdir -p $remoteDir" 2>&1 | ForEach-Object {
        Write-Host "    remote-mkdir: $_"
    }
    if ($LASTEXITCODE -ne 0) { throw "[3] remote mkdir failed (exit $LASTEXITCODE)" }

    # 3b. Record the local source hash in the uploaded workspace for provenance
    Write-Host "[3b] Recording source hash..." -ForegroundColor Cyan
    $sourceHash = (git rev-parse HEAD 2>&1)
    if ($LASTEXITCODE -ne 0) {
        throw "[3b] git rev-parse HEAD failed (exit $LASTEXITCODE). Output: $sourceHash"
    }
    $sourceHash = $sourceHash.Trim()
    ssh -o ConnectTimeout=10 -o BatchMode=yes dev@runner "printf '%s\n' '$sourceHash' > '$remoteDir/source_hash.txt'" 2>&1 | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "[3b] remote source_hash.txt write failed (exit $LASTEXITCODE)" }

    # 4. Upload workspace-level files
    Write-Host "[4] Uploading workspace files..." -ForegroundColor Cyan

    # Cargo.toml (workspace manifest)
    $workspaceToml = Join-Path $workspaceRoot "Cargo.toml"
    scp -o ConnectTimeout=10 -o BatchMode=yes $workspaceToml "dev@runner:$remoteDir/Cargo.toml" 2>&1 | Out-Null
    if ($LASTEXITCODE -ne 0) { throw "[4] scp Cargo.toml failed (exit $LASTEXITCODE)" }

    # Cargo.lock if exists
    $lockFile = Join-Path $workspaceRoot "Cargo.lock"
    if (Test-Path $lockFile) {
        scp -o ConnectTimeout=10 -o BatchMode=yes $lockFile "dev@runner:$remoteDir/Cargo.lock" 2>&1 | Out-Null
        if ($LASTEXITCODE -ne 0) { throw "[4] scp Cargo.lock failed (exit $LASTEXITCODE)" }
    }

    # .cargo config if exists
    $cargoConfig = Join-Path $workspaceRoot ".cargo"
    if (Test-Path $cargoConfig) {
        scp -r -o ConnectTimeout=10 -o BatchMode=yes $cargoConfig "dev@runner:$remoteDir/" 2>&1 | Out-Null
        if ($LASTEXITCODE -ne 0) { throw "[4] scp .cargo failed (exit $LASTEXITCODE)" }
    }

    Write-Host "[4] Workspace files uploaded" -ForegroundColor Cyan

    # 5. Upload each package's necessary files
    Write-Host "[5] Uploading package sources..." -ForegroundColor Cyan
    foreach ($pkg in $metadata.packages) {
        $pkgManifest = $pkg.manifest_path
        $pkgDir = Split-Path $pkgManifest -Parent

        # Calculate relative path from workspace root
        $relPath = $pkgDir.Substring($workspaceRoot.Length).TrimStart('\', '/')

        if ([string]::IsNullOrEmpty($relPath)) {
            # Root package - files go directly into remote dir
            $targetDir = $remoteDir
        }
        else {
            # Nested package - preserve directory structure
            $relPathUnix = $relPath -replace '\\', '/'
            $targetDir = "$remoteDir/$relPathUnix"

            # Create nested dir on remote
            ssh -o ConnectTimeout=10 -o BatchMode=yes dev@runner "mkdir -p $targetDir" 2>&1 | Out-Null
            if ($LASTEXITCODE -ne 0) { throw "[5] remote mkdir $targetDir failed (exit $LASTEXITCODE)" }
        }

        Write-Host "    Uploading package: $($pkg.name)" -ForegroundColor Cyan

        # Upload package Cargo.toml (only if not workspace root)
        if (-not [string]::IsNullOrEmpty($relPath)) {
            scp -o ConnectTimeout=10 -o BatchMode=yes $pkgManifest "dev@runner:$targetDir/Cargo.toml" 2>&1 | Out-Null
            if ($LASTEXITCODE -ne 0) { throw "[5] scp package Cargo.toml failed (exit $LASTEXITCODE)" }
        }

        # Upload src/ if exists
        $srcDir = Join-Path $pkgDir "src"
        if (Test-Path $srcDir) {
            scp -r -o ConnectTimeout=10 -o BatchMode=yes $srcDir "dev@runner:$targetDir/" 2>&1 | Out-Null
            if ($LASTEXITCODE -ne 0) { throw "[5] scp src failed (exit $LASTEXITCODE)" }
        }

        # Upload build.rs if exists
        $buildRs = Join-Path $pkgDir "build.rs"
        if (Test-Path $buildRs) {
            scp -o ConnectTimeout=10 -o BatchMode=yes $buildRs "dev@runner:$targetDir/build.rs" 2>&1 | Out-Null
            if ($LASTEXITCODE -ne 0) { throw "[5] scp build.rs failed (exit $LASTEXITCODE)" }
        }
    }

    Write-Host "[5] All sources uploaded" -ForegroundColor Cyan

    # 6. Start systemd service
    Write-Host "[6] Starting rcq@$JID.service..." -ForegroundColor Cyan
    $remoteCmd = "bash -lc 'set -Eeuo pipefail; echo START; sudo -n systemctl start --no-block rcq@$JID.service; echo SYSTEMCTL_OK'"
    ssh -o ConnectTimeout=10 -o BatchMode=yes dev@runner $remoteCmd 2>&1 | ForEach-Object {
        Write-Host "    remote: $_"
    }
    if ($LASTEXITCODE -ne 0) { throw "[6] systemctl failed (exit $LASTEXITCODE)" }
    Write-Host "[6] Service started" -ForegroundColor Cyan

    Write-Host "Submitted JID: $JID" -ForegroundColor Green
    exit 0
}
catch {
    Write-Host "ERROR: $_" -ForegroundColor Red
    Write-Host "Common causes:" -ForegroundColor Yellow
    Write-Host "  - SSH/Tailscale auth failed (BatchMode=yes disables prompts)" -ForegroundColor Yellow
    Write-Host "  - cargo not in PATH or project not valid Cargo workspace" -ForegroundColor Yellow
    Write-Host "  - sudo on runner not configured for passwordless systemctl" -ForegroundColor Yellow
    exit 1
}
