$env:PATH = "C:\MinGW\bin;$env:PATH"
$tmpC = Join-Path $env:TEMP "flagprobe.c"
$tmpO = Join-Path $env:TEMP "flagprobe.o"
Set-Content -NoNewline -Encoding ascii -Path $tmpC -Value 'int main(void){return 0;}'

$candidates = @(
 '-Wshadow','-Wcast-qual','-Wcast-align','-Wdouble-promotion','-Wformat=2',
 '-Wmissing-prototypes','-Wstrict-prototypes','-Wold-style-definition','-Wredundant-decls',
 '-Wundef','-Wunreachable-code','-Wswitch-enum','-Wimplicit-fallthrough','-Wlogical-op',
 '-Wshift-overflow','-Wduplicated-branches','-Wdangling-else','-Wnull-dereference',
 '-Wwrite-strings','-Wnested-externs','-Wpointer-arith','-Winit-self','-Wstrict-overflow=5',
 '-Wvla','-Wmissing-declarations','-Wconversion','-Wsign-conversion','-Wfloat-equal',
 '-Wpadded','-Wformat-overflow','-Wformat-truncation','-Warray-bounds=2',
 '-Wmaybe-uninitialized','-Walloc-zero','-Wduplicated-cond','-Wint-conversion',
 '-Wsizeof-pointer-memaccess','-Wstringop-overflow','-Wswitch-default','-Wunused-macros',
 '-Wvariadic-macros','-Wsystem-headers'
)
$ok = @(); $no = @()
foreach ($f in $candidates) {
    $res = & gcc $f -std=c11 -c $tmpC -o $tmpO 2>&1 | Out-String
    if ($res -match 'unrecognized|unknown option|invalid') { $no += $f } else { $ok += $f }
}
"===== SUPPORTED ($($ok.Count)) ====="
$ok -join ' '
""
"===== UNSUPPORTED ($($no.Count)) ====="
$no -join ' '
