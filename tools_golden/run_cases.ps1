# Golden-output harness for the WHAMP C translation.
#
# Runs a set of scripted terminal sessions against a WHAMP C executable
# (the original srcc build or the srccopt build) and stores/compares the
# captured stdout.
#
#   .\run_cases.ps1 -Exe ..\srcc\whamp.exe     -OutDir golden
#   .\run_cases.ps1 -Exe ..\srccopt\whamp.exe  -OutDir actual -CompareDir golden
#
# cases.json entries:
#   name    - case id (output file name)
#   model   - model file path relative to the repo root ("" = no -file)
#   args    - extra command line arguments
#   input   - text fed on stdin ('\n' escapes expanded)
#   engine  - $true when the executable is whamp_engine_test

param(
    [Parameter(Mandatory = $true)][string]$Exe,
    [Parameter(Mandatory = $true)][string]$OutDir,
    [string]$CompareDir = "",
    [string]$Cases = "$PSScriptRoot\cases.json"
)

$ErrorActionPreference = "Stop"
$repoRoot = Split-Path -Parent $PSScriptRoot
$Exe = (Resolve-Path $Exe).Path
$Cases = (Resolve-Path $Cases).Path
if ($CompareDir) { $CompareDir = (Resolve-Path $CompareDir).Path }
if (-not [System.IO.Path]::IsPathRooted($OutDir)) { $OutDir = Join-Path (Get-Location).Path $OutDir }
Push-Location $repoRoot
try {
    New-Item -ItemType Directory -Force -Path $OutDir | Out-Null
    $defs = Get-Content -Raw $Cases | ConvertFrom-Json
    $failures = @()
    $ran = 0

    $isEngineExe = [bool]($Exe -match 'engine')
    foreach ($c in $defs) {
        if ([bool]$c.engine -ne $isEngineExe) { continue }
        $argv = @()
        if ($c.model) { $argv += @("-file", $c.model) }
        if ($c.args) { $argv += ($c.args -split '\s+') | Where-Object { $_ } }

        $stdinText = $c.input -replace '\\n', "`n"
        $psi = New-Object System.Diagnostics.ProcessStartInfo
        $psi.FileName = (Resolve-Path $Exe).Path
        $psi.Arguments = ($argv | ForEach-Object { if ($_ -match '\s') { "`"$_`"" } else { $_ } }) -join ' '
        $psi.WorkingDirectory = $repoRoot
        $psi.RedirectStandardInput = $true
        $psi.RedirectStandardOutput = $true
        $psi.RedirectStandardError = $true
        $psi.UseShellExecute = $false

        $p = [System.Diagnostics.Process]::Start($psi)
        $p.StandardInput.Write($stdinText)
        $p.StandardInput.Close()
        # Read asynchronously so a chatty program cannot deadlock the pipe.
        $readTask = $p.StandardOutput.ReadToEndAsync()
        $errTask = $p.StandardError.ReadToEndAsync()
        if (-not $p.WaitForExit(60000)) {
            try { Stop-Process -Id $p.Id -Force } catch { }
            Write-Host "TIMEOUT $($c.name)"
            $out = ""; $err = "timeout"
        }
        else {
            $out = $readTask.Result
            $err = $errTask.Result
        }

        # Normalise: strip trailing blanks per line, drop trailing blank lines.
        $lines = ($out -split "`r?`n") | ForEach-Object { $_.TrimEnd() }
        while ($lines.Count -gt 0 -and $lines[-1] -eq "") {
            $lines = $lines[0..($lines.Count - 2)]
        }
        $normalised = ($lines -join "`n") + "`n"
        $ran++
        Set-Content -NoNewline -Encoding ascii -Path "$OutDir\$($c.name).txt" -Value $normalised
        if ($err) { Set-Content -NoNewline -Encoding ascii -Path "$OutDir\$($c.name).stderr.txt" -Value $err }

        if ($CompareDir) {
            $ref = Get-Content -Raw "$CompareDir\$($c.name).txt"
            if ($ref -ne $normalised) {
                $failures += $c.name
                Write-Host "FAIL $($c.name)"
            }
            else {
                Write-Host "ok   $($c.name)"
            }
        }
    }

    if ($CompareDir) {
        if ($failures.Count -gt 0) {
            Write-Host "`n$($failures.Count) case(s) differ: $($failures -join ', ')"
            exit 1
        }
        Write-Host "`nAll $ran cases match."
    }
    else {
        Write-Host "Wrote $ran golden files to $OutDir"
    }
}
finally {
    Pop-Location
}
