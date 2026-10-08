$s = [Text.Encoding]::ASCII.GetString([IO.File]::ReadAllBytes("c:\Users\WP\VSCODE\whamp\fortran_ref\probe16_raw.bin"))
$l = $s -split "`r`n"
$l2 = $l[2]
"L2 len=" + $l2.Length
$idx = $l2.IndexOf("(")
$idx2 = $l2.IndexOf(")")
"index of open paren = $idx"
"index of close paren = $idx2"
"complex token len = " + ($idx2 - $idx + 1)
"after complex: [" + $l2.Substring($idx2 + 1) + "]"
$l3 = $l[3]
"L3 len=" + $l3.Length
$i3 = $l3.IndexOf("(")
"L3 idx of open paren = $i3"
"L3 before paren: [" + $l3.Substring(0, $i3) + "]"
$l0 = $l[0]
"L0 len=" + $l0.Length
"L0: [" + $l0 + "]"
