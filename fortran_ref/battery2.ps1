param([string]$Only = "", [int]$Limit = 200000)
# Full differential battery: Fortran reference vs C translation.
# Uses bounded stdout reads so the known EOF infinite-loop cases stay small.
$ref = "c:\Users\WP\VSCODE\whamp\fortran_ref"
$S1 = "p0z.0022f.1`npzf`ns`n"
$S2 = "p0z.0022f.1`npzf`no`npzf/e`ns`n"

function RunBounded([string]$exe, [string]$outFile, [string]$inputText, [string]$extra, [string]$model, [int]$limit, [double]$secs) {
  $psi = New-Object System.Diagnostics.ProcessStartInfo
  $psi.FileName = $exe
  $psi.Arguments = "-file ..\Models\$model $extra"
  $psi.WorkingDirectory = "c:\Users\WP\VSCODE\whamp\src"
  $psi.RedirectStandardInput = $true
  $psi.RedirectStandardOutput = $true
  $psi.RedirectStandardError = $false
  $psi.UseShellExecute = $false
  $p = [System.Diagnostics.Process]::Start($psi)
  try { $p.StandardInput.Write($inputText); $p.StandardInput.Close() } catch {}
  $buf = New-Object byte[] 65536
  $ms = New-Object System.IO.MemoryStream
  $sw = [Diagnostics.Stopwatch]::StartNew()
  while ($ms.Length -lt $limit -and $sw.Elapsed.TotalSeconds -lt $secs) {
    $n = $p.StandardOutput.BaseStream.Read($buf, 0, [Math]::Min($buf.Length, $limit - $ms.Length))
    if ($n -le 0) { break }
    $ms.Write($buf, 0, $n)
  }
  $bytes = $ms.ToArray()
  try { $p.Kill() } catch {}
  [IO.File]::WriteAllBytes($outFile, $bytes)
  return $bytes.Length
}

$cases = New-Object System.Collections.ArrayList
function Add1($tag,$in,$model,$extra) { [void]$cases.Add([pscustomobject]@{Tag=$tag;In=$in;Model=$model;Extra=$extra}) }

Add1 "T1" $S1 "Ex1" ""
Add1 "T2" $S2 "Ex1" ""
Add1 "T3" $S1 "Ex1" "-debug"
Add1 "T4" $S1 "Ex1" "-debug -maxiterations 20"
foreach ($L in "b","d","e","f","g","h","l","m","n","o","p","r","s","t","u","v","x","y","z") {
  Add1 ("L_" + $L) ("p0z.0022f.1`npzf`n$L`ns`n") "Ex1" ""
}
foreach ($L in "b","d","e","f","g","h","l","m","n","o","p","r","s","t","u","v","x","y","z") {
  Add1 ("X2_" + $L) ("p0z.0022f.1`npzf`n$L`ns`n") "Ex2" ""
}
Add1 "R_help" "`n" "Ex1" ""
Add1 "R_bad1" "p0z.0022f.1`npzf`nqqq`ns`n" "Ex1" ""
Add1 "R_bad2" "p0z.0022f.1`npzf`n@#`$`ns`n" "Ex1" ""
Add1 "R_fonly" "f.05`nf`ns`n" "Ex1" ""
Add1 "R_comb" "p.1.2.01z.1.2.01f.1`npzf`ns`n" "Ex1" ""
Add1 "R_two" "p.1.2.01`nz.1.2.01`nf.1`npzf`ns`n" "Ex1" ""
Add1 "R_big" "p.001.5.005z.001.5.005f.001.5.005`npzf`ns`n" "Ex1" ""
Add1 "R_log" "p.01.10.3l`npzf`ns`n" "Ex1" ""
Add1 "R_par1" "N=1.e7 C=28.`np0z.0022f.1`npzf`ns`n" "Ex1" ""
Add1 "R_par2" "N=1.e9 C=28.`np0z.0022f.1`npzf`ns`n" "Ex1" ""
Add1 "R_par3" "N=1.e5 C=28.`np0z.0022f.1`npzf`ns`n" "Ex1" ""
Add1 "R_par4" "N=1.e7 C=1.`np0z.0022f.1`npzf`ns`n" "Ex1" ""
Add1 "R_par5" "N=1.e7 C=100.`np0z.0022f.1`npzf`ns`n" "Ex1" ""
Add1 "R_par6" "N=1.e7 C=28. B=5000.`np0z.0022f.1`npzf`ns`n" "Ex1" ""
Add1 "R_negf" "p0z.0022f.1`npzf`nf-.1`npzf`ns`n" "Ex1" ""
Add1 "X2_T1" $S1 "Ex2" ""
Add1 "X2_T2" $S2 "Ex2" ""
Add1 "X2_dbg" $S1 "Ex2" "-debug"

$results = @()
foreach ($c in $cases) {
  if ($Only -ne "" -and $c.Tag -notlike $Only) { continue }
  $fo = "$ref\$($c.Tag)_f.txt"; $co = "$ref\$($c.Tag)_c.txt"
  $nf = RunBounded "$ref\whamp_ref.exe" $fo $c.In $c.Extra $c.Model $Limit 30
  $nc = RunBounded "c:\Users\WP\VSCODE\whamp\srcc\whamp_c.exe" $co $c.In $c.Extra $c.Model $Limit 30
  $fb = [IO.File]::ReadAllBytes($fo); $cb = [IO.File]::ReadAllBytes($co)
  $n = [Math]::Min($fb.Length, $cb.Length)
  $same = $true; $at = -1
  for ($i=0; $i -lt $n; $i++) { if ($fb[$i] -ne $cb[$i]) { $same=$false; $at=$i; break } }
  $trunc = ""
  if ($fb.Length -ge $Limit -or $cb.Length -ge $Limit) { $trunc = " (truncated@${Limit})" }
  if ($same -and $fb.Length -eq $cb.Length) {
    "{0,-10} {1,-4} {2,8}  IDENTICAL{3}" -f $c.Tag,$c.Model,$fb.Length,$trunc
  } elseif ($same) {
    "{0,-10} {1,-4} {2,8}/{3,-8}  PREFIX-MATCH first $n bytes{4}" -f $c.Tag,$c.Model,$fb.Length,$cb.Length,$trunc
  } else {
    "{0,-10} {1,-4} {2,8}/{3,-8}  DIFFER at byte $at{4}" -f $c.Tag,$c.Model,$fb.Length,$cb.Length,$trunc
    $s = [Math]::Max(0,$at-100)
    $e = [Math]::Min([Math]::Max($fb.Length,$cb.Length)-1, $at+100)
    $sf = [Text.Encoding]::ASCII.GetString($fb[$s..[Math]::Min($e,$fb.Length-1)])
    $sc = [Text.Encoding]::ASCII.GetString($cb[$s..[Math]::Min($e,$cb.Length-1)])
    "    F: [$sf]"
    "    C: [$sc]"
  }
}
