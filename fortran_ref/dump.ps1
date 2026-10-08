$psi = New-Object System.Diagnostics.ProcessStartInfo
$psi.FileName = "$PWD\probe13.exe"
$psi.WorkingDirectory = "$PWD"
$psi.RedirectStandardOutput = $true
$psi.UseShellExecute = $false
$p = [System.Diagnostics.Process]::Start($psi)
$out = $p.StandardOutput.ReadToEnd()
$p.WaitForExit()
[System.IO.File]::WriteAllText("$PWD\probe13_out.txt", $out)
$lines = $out -split "`r?`n"
$n = 0
foreach ($l in $lines) {
    $n++
    "=== L$n len=$($l.Length)"
    # print in chunks of 10 with index
    for ($i = 0; $i -lt $l.Length; $i += 10) {
        $chunk = $l.Substring($i, [Math]::Min(10, $l.Length - $i))
        "  [$i] $($chunk -replace ' ', '.')"
    }
}
