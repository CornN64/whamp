param([string]$F, [string]$C)
$fl = [System.IO.File]::ReadAllLines($F)
$cl = [System.IO.File]::ReadAllLines($C)
for ($i = 0; $i -lt [Math]::Max($fl.Count, $cl.Count); $i++) {
  $a = if ($i -lt $fl.Count) { $fl[$i] } else { "" }
  $b = if ($i -lt $cl.Count) { $cl[$i] } else { "" }
  if ($a -ceq $b) { continue }
  $ta = $a -split '\s+' | Where-Object { $_ -ne "" }
  $tb = $b -split '\s+' | Where-Object { $_ -ne "" }
  $lr = 0.0
  for ($k = 0; $k -lt [Math]::Min($ta.Count, $tb.Count); $k++) {
    $va = 0.0; $vb = 0.0
    if (-not ([double]::TryParse($ta[$k], [ref]$va) -and [double]::TryParse($tb[$k], [ref]$vb))) { continue }
    $den = [Math]::Max([Math]::Abs($va), [Math]::Abs($vb))
    if ($den -eq 0) { continue }
    $rel = [Math]::Abs($va - $vb) / $den
    if ($rel -gt $lr) { $lr = $rel }
  }
  "L{0,3}  maxRel={1,12:E4}  ~{2,8:N1} printed-digit-ULP" -f ($i + 1), $lr, ($lr * 1e8)
}
