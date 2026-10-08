param([string]$Exe, [string]$WorkDir, [string]$OutFile, [string]$InputText, [string]$ExtraArgs = "", [int]$TimeoutSec = 25, [string]$Model = "Ex1")
$psi = New-Object System.Diagnostics.ProcessStartInfo
$psi.FileName = $Exe
$psi.Arguments = "-file ..\Models\$Model $ExtraArgs"
$psi.WorkingDirectory = $WorkDir
$psi.RedirectStandardInput = $true
$psi.RedirectStandardOutput = $true
$psi.RedirectStandardError = $true
$psi.UseShellExecute = $false
$p = [System.Diagnostics.Process]::Start($psi)
try { $p.StandardInput.Write($InputText); $p.StandardInput.Close() } catch {}
$readTask = $p.StandardOutput.ReadToEndAsync()
if (-not $p.WaitForExit($TimeoutSec * 1000)) {
  try { $p.Kill() } catch {}
  "TIMEOUT (killed after ${TimeoutSec}s)"
}
$out = $readTask.Result
"ExitCode: $($p.ExitCode)  bytes: $($out.Length)"
[IO.File]::WriteAllBytes($OutFile, [Text.Encoding]::UTF8.GetBytes($out))
