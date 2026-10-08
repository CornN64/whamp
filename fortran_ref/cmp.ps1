$f = Get-Content f_out.txt
$c = Get-Content c_out.txt
"fortran lines: $($f.Count)  c lines: $($c.Count)"
$n = [Math]::Min($f.Count, $c.Count)
$bad = 0
for ($i = 0; $i -lt $n; $i++) {
    if ($f[$i] -cne $c[$i]) {
        $bad++
        if ($bad -le 30) {
            "line $($i+1):"
            "  F: [$($f[$i] -replace ' ', '.')]"
            "  C: [$($c[$i] -replace ' ', '.')]"
        }
    }
}
"total mismatched lines: $bad"
