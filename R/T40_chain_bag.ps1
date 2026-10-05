# Queue after the T40 main run (its process id is the argument):
#   1. T42 (intercept cells, Amendment 4), 22 workers;
#   2. the independent second block of the p/n >= 0.25 null cells (Amendment 5), 22 workers;
#   3. BAGofT, one R process per test (Amendment 3), 22 at a time.
param([int]$MainPid)
$r = "C:\Users\ebrah\.gemini\Projects\PDFs\arashi_proposal\gof-penalized-paper\R"
try { Wait-Process -Id $MainPid -ErrorAction Stop } catch { }
$env:T42_NW = "22"
Start-Process -FilePath "Rscript.exe" -ArgumentList "T42_intercept.R", "main" -WorkingDirectory $r `
  -RedirectStandardOutput "$r\..\data\T42_intercept.log" -RedirectStandardError "$r\..\data\T42_intercept.err" -WindowStyle Hidden -Wait
$env:T40_NW = "22"
Start-Process -FilePath "Rscript.exe" -ArgumentList "T40_highdim_grid.R", "recheck" -WorkingDirectory $r `
  -RedirectStandardOutput "$r\..\data\T40_recheck.log" -RedirectStandardError "$r\..\data\T40_recheck.err" -WindowStyle Hidden -Wait
Start-Process -FilePath "Rscript.exe" -ArgumentList "T40_bag_runner.R" -WorkingDirectory $r `
  -RedirectStandardOutput "$r\..\data\T40_bag.log" -RedirectStandardError "$r\..\data\T40_bag.err" -WindowStyle Hidden -Wait
