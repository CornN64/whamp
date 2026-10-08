param([string]$Exe, [string]$WorkDir, [string]$OutFile, [string]$Input)
$psi = New-Object System.Diagnostics.ProcessStartInfo
$psi.FileName = $Exe
$psi.Arguments = "-file ..\Models\Ex1"
$psi.WorkingDirectory = $WorkDir
$psi.RedirectStandardInput = $true
$psi.RedirectStandardOutput = $true
$psi.RedirectStandardError = $true
$psi.UseShellExecute = $false
$p = [System.Diagnostics.Process]::Start($psi)
$p.StandardInput.Write($Input)
$p.StandardInput.Close()
$out = $p.StandardOutput.ReadToEnd()
$err = $p.StandardError.ReadToEnd()
$p.WaitForExit()
[IO.File]::WriteAllBytes($OutFile, [Text.Encoding]::UTF8.GetBytes($out))
if ($err) { "STDERR: $err" }
"ExitCode: $($p.ExitCode)"
