# Queue after the T40 recheck (its process id is the argument), replacing the tail of T40_chain_bag.ps1:
#   1. the independent second block of the T42 p/n = 0.25 null cell (Amendment 6), 22 workers;
#   2. BAGofT, one R process per test (Amendment 3), 22 at a time.
param([int]$RecheckPid)
$r = "C:\Users\ebrah\.gemini\Projects\PDFs\arashi_proposal\gof-penalized-paper\R"
try { Wait-Process -Id $RecheckPid -ErrorAction Stop } catch { }
$env:T42_NW = "22"
Start-Process -FilePath "Rscript.exe" -ArgumentList "T42_intercept.R", "recheck" -WorkingDirectory $r `
  -RedirectStandardOutput "$r\..\data\T42_recheck.log" -RedirectStandardError "$r\..\data\T42_recheck.err" -WindowStyle Hidden -Wait
$env:T40_NW = "22"
Start-Process -FilePath "Rscript.exe" -ArgumentList "T40_bag_runner.R" -WorkingDirectory $r `
  -RedirectStandardOutput "$r\..\data\T40_bag.log" -RedirectStandardError "$r\..\data\T40_bag.err" -WindowStyle Hidden -Wait
