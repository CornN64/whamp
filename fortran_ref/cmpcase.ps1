param([string]$InputText, [string]$ExtraArgs = "", [string]$Tag = "T")
$fexe = "c:\Users\WP\VSCODE\whamp\fortran_ref\whamp_ref.exe"
$cexe = "c:\Users\WP\VSCODE\whamp\srcc\whamp_c.exe"
$fo = "c:\Users\WP\VSCODE\whamp\fortran_ref\${Tag}_f.txt"
$co = "c:\Users\WP\VSCODE\whamp\fortran_ref\${Tag}_c.txt"
& "c:\Users\WP\VSCODE\whamp\fortran_ref\runcase.ps1" -Exe $fexe -WorkDir "c:\Users\WP\VSCODE\whamp\src" -OutFile $fo -InputText $InputText -ExtraArgs $ExtraArgs
& "c:\Users\WP\VSCODE\whamp\fortran_ref\runcase.ps1" -Exe $cexe -WorkDir "c:\Users\WP\VSCODE\whamp\src" -OutFile $co -InputText $InputText -ExtraArgs $ExtraArgs
$fb = [IO.File]::ReadAllBytes($fo)
$cb = [IO.File]::ReadAllBytes($co)
if ($fb.Length -eq $cb.Length) {
  $same = $true
  for ($i = 0; $i -lt $fb.Length; $i++) { if ($fb[$i] -ne $cb[$i]) { $same = $false; break } }
  if ($same) { "${Tag}: IDENTICAL ($($fb.Length) bytes)"; return }
}
"${Tag}: DIFFER (fortran=$($fb.Length) bytes, c=$($cb.Length) bytes)"
$fl = [Text.Encoding]::ASCII.GetString($fb) -split "`r`n"
$cl = [Text.Encoding]::ASCII.GetString($cb) -split "`r`n"
$n = [Math]::Max($fl.Count, $cl.Count)
$d = 0
for ($i = 0; $i -lt $n; $i++) {
  $a = ""; if ($i -lt $fl.Count) { $a = $fl[$i] }
  $b = ""; if ($i -lt $cl.Count) { $b = $cl[$i] }
  if ($a -cne $b) {
    "  L${i} F: [$a]"
    "  L${i} C: [$b]"
    $d++
    if ($d -ge 25) { "  ... (truncated)"; break }
  }
}
