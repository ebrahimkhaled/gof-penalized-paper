# B2: promote GlaucomaMVF out of AMBER by re-deriving its FULL 12-cell grid from the
# saved run, so the supplement stops printing unregistered single-configuration numbers.
x <- readRDS("T26_app_b_results.rds")
g <- x$grid
cat("dataset(s):", paste(unique(g$dataset), collapse = ", "), "\n")
cat("rows:", nrow(g), " lambdas:", paste(unique(g$lambda), collapse = " / "), "\n\n")
print(g, digits = 4, row.names = FALSE)

cat("\n--- how often does each UNCORRECTED test reject at 5%? ---\n")
cat(sprintf("  decile : %d of %d   (max p = %.4f)\n", sum(g$naive_dec  <= 0.05), nrow(g), max(g$naive_dec)))
cat(sprintf("  EDGE   : %d of %d   (max p = %.4f)\n", sum(g$naive_edge <= 0.05), nrow(g), max(g$naive_edge)))
cat("\n--- and the CORRECTED test? ---\n")
cat(sprintf("  decile : %d of %d reject   range [%.3f, %.3f]\n",
            sum(g$corr_dec <= 0.05), nrow(g), min(g$corr_dec), max(g$corr_dec)))
cat(sprintf("  EDGE   : %d of %d reject   range [%.3f, %.3f]\n",
            sum(g$corr_edge <= 0.05), nrow(g), min(g$corr_edge), max(g$corr_edge)))
cat("\n--- seed stability at lambda.1se, G=10, NB=999 ---\n")
print(x$stab, digits = 3, row.names = FALSE)
cat(sprintf("  decile %.3f-%.3f ; EDGE %.3f-%.3f\n",
            min(x$stab$corr_dec), max(x$stab$corr_dec),
            min(x$stab$corr_edge), max(x$stab$corr_edge)))
cat("\nCELLS WHERE THE UNCORRECTED DECILE TEST DOES *NOT* REJECT:\n")
print(g[g$naive_dec > 0.05, c("lambda","G","naive_dec","corr_dec")], digits = 4, row.names = FALSE)
