param([string]$InputText, [string]$Tag = "ET", [int]$Limit = 20000, [string]$Model = "Ex1", [string]$ExtraArgs = "")
function RunBounded([string]$exe, [string]$outFile, [string]$workDir) {
  $psi = New-Object System.Diagnostics.ProcessStartInfo
  $psi.FileName = $exe
  $psi.Arguments = "-file ..\Models\$Model $ExtraArgs"
  $psi.WorkingDirectory = $workDir
  $psi.RedirectStandardInput = $true
  $psi.RedirectStandardOutput = $true
  $psi.RedirectStandardError = $false
  $psi.UseShellExecute = $false
  $p = [System.Diagnostics.Process]::Start($psi)
  try { $p.StandardInput.Write($InputText); $p.StandardInput.Close() } catch {}
  $buf = New-Object byte[] 65536
  $ms = New-Object System.IO.MemoryStream
  $sw = [Diagnostics.Stopwatch]::StartNew()
  while ($ms.Length -lt $Limit -and $sw.Elapsed.TotalSeconds -lt 25) {
    $n = $p.StandardOutput.BaseStream.Read($buf, 0, [Math]::Min($buf.Length, $Limit - $ms.Length))
    if ($n -le 0) { break }
    $ms.Write($buf, 0, $n)
  }
  try { $p.Kill() } catch {}
  [IO.File]::WriteAllBytes($outFile, $ms.ToArray())
}
$fo = "c:\Users\WP\VSCODE\whamp\fortran_ref\${Tag}_f.txt"
$co = "c:\Users\WP\VSCODE\whamp\fortran_ref\${Tag}_c.txt"
RunBounded "c:\Users\WP\VSCODE\whamp\fortran_ref\engtest\whamp_engine_test_ref.exe" $fo "c:\Users\WP\VSCODE\whamp\fortran_ref\engtest"
RunBounded "c:\Users\WP\VSCODE\whamp\srcc\whamp_engine_test_c.exe" $co "c:\Users\WP\VSCODE\whamp\fortran_ref\engtest"
$fb = [IO.File]::ReadAllBytes($fo); $cb = [IO.File]::ReadAllBytes($co)
$n = [Math]::Min($fb.Length, $cb.Length)
$same = $true; $at = -1
for ($i = 0; $i -lt $n; $i++) { if ($fb[$i] -ne $cb[$i]) { $same = $false; $at = $i; break } }
if ($same -and $fb.Length -eq $cb.Length) { "${Tag}: IDENTICAL ($n bytes)"; return }
if ($same) { "${Tag}: PREFIX-MATCH first $n bytes (f=$($fb.Length) c=$($cb.Length))"; return }
"${Tag}: DIFFER at byte $at (f=$($fb.Length) c=$($cb.Length))"
$s = [Math]::Max(0, $at - 150)
"F: [" + ([Text.Encoding]::ASCII.GetString($fb[$s..[Math]::Min($fb.Length-1,$at+150)])) + "]"
"C: [" + ([Text.Encoding]::ASCII.GetString($cb[$s..[Math]::Min($cb.Length-1,$at+150)])) + "]"
