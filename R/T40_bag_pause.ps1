# Pause or resume the BAGofT run without losing the tests in flight: the runner and every R process it
# started are suspended (no CPU, memory kept) or resumed through the Windows NtSuspendProcess/NtResumeProcess
# calls. Used to lend the machine's cores to another simulation. A reboot while paused loses only the tests
# in flight; finished rows are on disk and the runner resumes from them.
#   powershell -File T40_bag_pause.ps1 -Action suspend|resume
param([ValidateSet("suspend", "resume")][string]$Action)
Add-Type -Namespace Win -Name Nt -MemberDefinition @'
[DllImport("ntdll.dll")] public static extern int NtSuspendProcess(IntPtr h);
[DllImport("ntdll.dll")] public static extern int NtResumeProcess(IntPtr h);
[DllImport("kernel32.dll")] public static extern IntPtr OpenProcess(int access, bool inherit, int pid);
[DllImport("kernel32.dll")] public static extern bool CloseHandle(IntPtr h);
'@
$procs = Get-CimInstance Win32_Process | Where-Object { $_.Name -eq "Rscript.exe" -and $_.CommandLine -match "T40_bag_(runner|one)\.R" }
$n = 0
foreach ($p in $procs) {
  $h = [Win.Nt]::OpenProcess(0x0800, $false, [int]$p.ProcessId)     # PROCESS_SUSPEND_RESUME
  if ($h -ne [IntPtr]::Zero) {
    if ($Action -eq "suspend") { [void][Win.Nt]::NtSuspendProcess($h) } else { [void][Win.Nt]::NtResumeProcess($h) }
    [void][Win.Nt]::CloseHandle($h); $n++
  }
}
"{0}: {1} BAGofT processes" -f $Action, $n
