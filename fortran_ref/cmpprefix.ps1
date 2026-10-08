param([string]$InputText, [string]$Tag = "P", [int]$Limit = 6000, [string]$Model = "Ex1", [string]$ExtraArgs = "")
# Run both exes, read only the first $Limit bytes of stdout each, byte-compare.
function RunBounded([string]$exe, [string]$outFile) {
  $psi = New-Object System.Diagnostics.ProcessStartInfo
  $psi.FileName = $exe
  $psi.Arguments = "-file ..\Models\$Model $ExtraArgs"
  $psi.WorkingDirectory = "c:\Users\WP\VSCODE\whamp\src"
  $psi.RedirectStandardInput = $true
  $psi.RedirectStandardOutput = $true
  $psi.RedirectStandardError = $false
  $psi.UseShellExecute = $false
  $p = [System.Diagnostics.Process]::Start($psi)
  try { $p.StandardInput.Write($InputText); $p.StandardInput.Close() } catch {}
  $buf = New-Object byte[] 65536
  $ms = New-Object System.IO.MemoryStream
  $sw = [Diagnostics.Stopwatch]::StartNew()
  while ($ms.Length -lt $Limit -and $sw.Elapsed.TotalSeconds -lt 20) {
    $n = $p.StandardOutput.BaseStream.Read($buf, 0, [Math]::Min($buf.Length, $Limit - $ms.Length))
    if ($n -le 0) { break }
    $ms.Write($buf, 0, $n)
  }
  $bytes = $ms.ToArray()
  try { $p.Kill() } catch {}
  [IO.File]::WriteAllBytes($outFile, $bytes)
  return $bytes.Length
}
$fo = "c:\Users\WP\VSCODE\whamp\fortran_ref\${Tag}_f.txt"
$co = "c:\Users\WP\VSCODE\whamp\fortran_ref\${Tag}_c.txt"
$nf = RunBounded "c:\Users\WP\VSCODE\whamp\fortran_ref\whamp_ref.exe" $fo
$nc = RunBounded "c:\Users\WP\VSCODE\whamp\srcc\whamp_c.exe" $co
$fb = [IO.File]::ReadAllBytes($fo); $cb = [IO.File]::ReadAllBytes($co)
$n = [Math]::Min($fb.Length, $cb.Length)
$same = $true; $at = -1
for ($i = 0; $i -lt $n; $i++) { if ($fb[$i] -ne $cb[$i]) { $same = $false; $at = $i; break } }
if ($same -and $fb.Length -eq $cb.Length) { "${Tag}: IDENTICAL prefix ($n bytes)"; return }
if ($same) { "${Tag}: PREFIX-MATCH first $n bytes (f=$($fb.Length) c=$($cb.Length))"; return }
"${Tag}: DIFFER at byte $at (f=$($fb.Length) c=$($cb.Length))"
$s = [Math]::Max(0, $at - 120)
"F: [" + ([Text.Encoding]::ASCII.GetString($fb[$s..[Math]::Min($fb.Length-1,$at+120)])) + "]"
"C: [" + ([Text.Encoding]::ASCII.GetString($cb[$s..[Math]::Min($cb.Length-1,$at+120)])) + "]"
