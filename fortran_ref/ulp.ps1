param([string]$F, [string]$C)
# Compare numeric tokens between two outputs; report max relative difference and ULP spread.
$fl = [System.IO.File]::ReadAllLines($F)
$cl = [System.IO.File]::ReadAllLines($C)
$maxRel = 0.0
$maxWhere = ""
$nNum = 0
$nTok = 0
for ($i = 0; $i -lt [Math]::Max($fl.Count, $cl.Count); $i++) {
  $a = if ($i -lt $fl.Count) { $fl[$i] } else { "" }
  $b = if ($i -lt $cl.Count) { $cl[$i] } else { "" }
  if ($a -ceq $b) { continue }
  $ta = $a -split '\s+' | Where-Object { $_ -ne "" }
  $tb = $b -split '\s+' | Where-Object { $_ -ne "" }
  $nTok += [Math]::Max($ta.Count, $tb.Count)
  for ($k = 0; $k -lt [Math]::Min($ta.Count, $tb.Count); $k++) {
    $va = 0.0; $vb = 0.0
    $okA = [double]::TryParse($ta[$k], [ref]$va)
    $okB = [double]::TryParse($tb[$k], [ref]$vb)
    if (-not ($okA -and $okB)) { continue }
    $nNum++
    $den = [Math]::Max([Math]::Abs($va), [Math]::Abs($vb))
    if ($den -eq 0) { continue }
    $rel = [Math]::Abs($va - $vb) / $den
    if ($rel -gt $maxRel) {
      $maxRel = $rel
      $maxWhere = "L$($i+1) tok$k F=$($ta[$k]) C=$($tb[$k])"
    }
  }
}
"numeric tokens compared : $nNum / $nTok"
"max relative difference : $maxRel"
"max relative as ULP(52) : {0:N2}" -f ($maxRel * [Math]::Pow(2, 52))
"max at                  : $maxWhere"
