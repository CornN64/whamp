param(
    [Parameter(Mandatory = $true)][string]$Src,
    [Parameter(Mandatory = $true)][string]$Flags,
    [string]$ObjDir = "$env:TEMP\strictobj"
)
$ErrorActionPreference = 'Continue'
if (Test-Path $ObjDir) { Remove-Item -Recurse -Force $ObjDir }
New-Item -ItemType Directory -Force -Path $ObjDir | Out-Null

$sources = Get-ChildItem $Src -Filter *.c -File | Sort-Object Name
$all = @()
foreach ($f in $sources) {
    $out = "$ObjDir\$($f.BaseName).o"
    $args = @('-std=c11', '-O2') + ($Flags -split '\s+' | Where-Object { $_ -ne '' }) + @('-c', $f.FullName, '-o', $out)
    $res = & gcc @args 2>&1 | Out-String
    if ($res.Trim() -ne '') { $all += $res }
}
$all -join "`n"
$total = ($all | Measure-Object).Count
"===== SUMMARY $Src ====="
if ($all.Count -eq 0) { "CLEAN: no diagnostics from $($sources.Count) files" }
else {
    $joined = $all -join "`n"
    $warns = ([regex]::Matches($joined, 'warning:')).Count
    $errs  = ([regex]::Matches($joined, 'error:')).Count
    "files=$($sources.Count)  warning-lines=$warns  error-lines=$errs"
    "---- warning categories ----"
    $cats = [regex]::Matches($joined, '\[-W[a-z0-9+=-]+\]') | ForEach-Object { $_.Value }
    $cats | Group-Object | Sort-Object Count -Descending | ForEach-Object { "{0,6}  {1}" -f $_.Count, $_.Name }
}
