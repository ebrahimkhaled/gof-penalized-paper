# Rewrites ../sessionInfo.txt with every package the archived scripts load attached, plus the Python
# versions used by the two post-processing scripts. Run from dev/:  Rscript record_session.R
pkgs <- c("glmnet", "TH.data", "ipred", "ggplot2", "patchwork", "scales", "ggrepel", "CompQuadForm",
          "ebrahim.gof", "GRPtests", "RPtests", "PLStests", "BAGofT", "processx", "parallel")
ok <- vapply(pkgs, function(p) suppressWarnings(suppressPackageStartupMessages(
  require(p, character.only = TRUE, quietly = TRUE))), logical(1))
py <- tryCatch(system2("python", c("-c", shQuote(paste(
  "import sys, numpy, pandas, scipy;",
  "print('Python', sys.version.split()[0], '| numpy', numpy.__version__,",
  "'| pandas', pandas.__version__, '| scipy', scipy.__version__)"))), stdout = TRUE),
  error = function(e) "Python versions not available")
out <- c(
  sprintf("Recorded on %s on the machine that ran the scripts, with every package the archived", Sys.Date()),
  "scripts load attached. The grid of Section 5.4 (T40-T45, BAGofT) ran on this machine from 2 to",
  "6 October 2026; earlier studies ran in August-September 2026 with the package versions then",
  "installed (the EDGE agreement check used ebrahim.gof 2.4.0; see README, Software).",
  "",
  "Packages loaded by the archived scripts, and whether each is installed here:",
  capture.output(print(ok)),
  "",
  "Python used by T40_validity_map_data.py and T40_si_tables.py:",
  py,
  "",
  capture.output(sessionInfo()))
writeLines(out, file.path("..", "sessionInfo.txt"))
cat("sessionInfo.txt written;", sum(ok), "of", length(ok), "packages attached\n")
