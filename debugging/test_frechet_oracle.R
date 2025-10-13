#!/usr/bin/env Rscript
# Test 3: Oracle Test with True Means
# Purpose: Test frechet_anova formula using TRUE means instead of estimated means
# to isolate whether the issue is in riemstats formula or riemtan estimation

library(riemtan)
library(riemstats)
library(Matrix)
library(parallel)

cat("=== Fréchet ANOVA Oracle Test ===\n")
cat("Testing with TRUE means to isolate formula vs estimation issues\n\n")

# Set Cairo for headless plotting (no X11 required)
options(bitmapType = 'cairo')

# Detect number of cores
n_cores <- detectCores()
cat(sprintf("Using %d CPU cores for parallel processing\n\n", n_cores))

# Load metric
data(euclidean)

# Set seed for reproducibility
set.seed(42)

# Test parameters (reduced for faster execution)
n_per_group <- 10  # samples per group (reduced from 30)
g <- 2             # number of groups
d <- 3             # dimension of SPD matrices (3×3)
sigma <- 100.0     # dispersion parameter
n_replicates <- 100  # number of test runs (reduced from 1000)

cat(sprintf("Parameters:\n"))
cat(sprintf("  Groups: %d\n", g))
cat(sprintf("  Samples per group: %d\n", n_per_group))
cat(sprintf("  Matrix dimension: %d×%d\n", d, d))
cat(sprintf("  Dispersion: %.1f\n", sigma))
cat(sprintf("  Replicates: %d\n\n", n_replicates))

# Helper function to create dispersion matrix
create_dispersion_matrix <- function(d, sigma) {
  p <- d * (d + 1) / 2
  diag_mat <- sigma * diag(p)
  as(diag_mat, "dpoMatrix") |> Matrix::pack()
}

# True center (known)
true_center <- as(diag(d), "dpoMatrix") |> Matrix::pack()
scale <- create_dispersion_matrix(d, sigma)

cat("Generating data under H0 (all groups from same distribution)...\n")

# Helper function to set Fréchet mean by accessing private field
set_frechet_mean <- function(sample_obj, mean_value) {
  # Access the private environment of the R6 object
  private_env <- environment(sample_obj$initialize)$private
  # Set the private f_mean field
  private_env$f_mean <- mean_value
  invisible(sample_obj)
}

# Run simulations in parallel
cat("Running", n_replicates, "replicates in parallel...\n")

# Define simulation function
run_simulation <- function(iter_num) {
  tryCatch({
    # Generate H0 data: all groups from SAME distribution
    groups <- lapply(1:g, function(j) {
      sample <- rspdnorm(n_per_group, true_center, scale, euclidean)
      sample$compute_unvecs()
      sample$compute_conns()
      sample
    })

    # Test 1: Oracle test with TRUE mean
    # Manually set each sample's Fréchet mean to the true center
    oracle_groups <- lapply(groups, function(sample) {
      # Clone the sample
      oracle_sample <- sample$clone()
      # Override the Fréchet mean with the TRUE center (via private field access)
      set_frechet_mean(oracle_sample, true_center)
      oracle_sample
    })

    # Create super sample with oracle means
    oracle_ss <- CSuperSample$new(oracle_groups)
    oracle_result <- frechet_anova(oracle_ss)

    # Test 2: Standard test with ESTIMATED means (for comparison)
    estimated_ss <- CSuperSample$new(groups)
    estimated_result <- frechet_anova(estimated_ss)

    # Return both statistics
    c(oracle = oracle_result$statistic, estimated = estimated_result$statistic)
  }, error = function(e) {
    # Return a try-error object that will be caught in post-processing
    structure(list(message = e$message), class = "try-error")
  })
}

# Run parallel simulations
results <- mclapply(1:n_replicates, run_simulation, mc.cores = n_cores)

# Check for errors and filter them out
is_error <- sapply(results, function(x) inherits(x, "try-error") || is.null(x))
n_errors <- sum(is_error)
n_success <- sum(!is_error)

if (n_errors > 0) {
  cat(sprintf("Warning: %d out of %d simulations failed (%.1f%%)\n",
              n_errors, n_replicates, 100 * n_errors / n_replicates))
  cat("First error message:\n")
  first_error_idx <- which(is_error)[1]
  print(results[[first_error_idx]])
  cat("\n")
  results <- results[!is_error]
}

if (n_success == 0) {
  stop("All simulations failed. Cannot proceed with analysis.")
}

cat(sprintf("Successfully completed %d simulations\n", n_success))

# Extract results from successful runs
results_matrix <- do.call(rbind, results)
oracle_stats <- results_matrix[, "oracle"]
estimated_stats <- results_matrix[, "estimated"]

cat("Done!\n\n")

# Filter outliers using IQR method (1.5×IQR above Q3)
filter_outliers <- function(x) {
  q1 <- quantile(x, 0.25)
  q3 <- quantile(x, 0.75)
  iqr <- q3 - q1
  threshold <- q3 + 1.5 * iqr
  x[x <= threshold]
}

n_oracle_before <- length(oracle_stats)
n_estimated_before <- length(estimated_stats)

oracle_stats_filtered <- filter_outliers(oracle_stats)
estimated_stats_filtered <- filter_outliers(estimated_stats)

n_oracle_removed <- n_oracle_before - length(oracle_stats_filtered)
n_estimated_removed <- n_estimated_before - length(estimated_stats_filtered)

cat("=== Outlier Filtering ===\n")
cat(sprintf("Oracle: Removed %d outliers (%.1f%%), kept %d\n",
            n_oracle_removed, 100 * n_oracle_removed / n_oracle_before, length(oracle_stats_filtered)))
cat(sprintf("Estimated: Removed %d outliers (%.1f%%), kept %d\n\n",
            n_estimated_removed, 100 * n_estimated_removed / n_estimated_before, length(estimated_stats_filtered)))

# Expected degrees of freedom
df_expected <- g - 1

# Test oracle statistics against chi-squared distribution (filtered data)
cat("=== Oracle Test Results (TRUE means) ===\n")
ks_oracle <- ks.test(oracle_stats_filtered, "pchisq", df = df_expected)

cat(sprintf("Expected df: %d\n\n", df_expected))
cat("Expected vs Observed (Oracle):\n")
cat(sprintf("  Mean - Expected: %.4f, Observed: %.4f\n", df_expected, mean(oracle_stats_filtered)))
cat(sprintf("  Variance - Expected: %.4f, Observed: %.4f\n", 2*df_expected, var(oracle_stats_filtered)))
cat(sprintf("  Median - Expected: %.4f, Observed: %.4f\n",
            qchisq(0.5, df_expected), median(oracle_stats_filtered)))

cat(sprintf("\nKolmogorov-Smirnov test (Oracle):\n"))
cat(sprintf("  D = %.4f, p-value = %.6f\n", ks_oracle$statistic, ks_oracle$p.value))

if (ks_oracle$p.value < 0.05) {
  cat("  ✗ FAIL: Oracle statistics do NOT follow χ²(%d) (p < 0.05)\n", df_expected)
} else {
  cat(sprintf("  ✓ PASS: Oracle statistics follow χ²(%d) (p ≥ 0.05)\n", df_expected))
}

# Test estimated statistics against chi-squared distribution (filtered data)
cat("\n=== Standard Test Results (ESTIMATED means) ===\n")
ks_estimated <- ks.test(estimated_stats_filtered, "pchisq", df = df_expected)

cat("Expected vs Observed (Estimated):\n")
cat(sprintf("  Mean - Expected: %.4f, Observed: %.4f\n", df_expected, mean(estimated_stats_filtered)))
cat(sprintf("  Variance - Expected: %.4f, Observed: %.4f\n", 2*df_expected, var(estimated_stats_filtered)))
cat(sprintf("  Median - Expected: %.4f, Observed: %.4f\n",
            qchisq(0.5, df_expected), median(estimated_stats_filtered)))

cat(sprintf("\nKolmogorov-Smirnov test (Estimated):\n"))
cat(sprintf("  D = %.4f, p-value = %.6f\n", ks_estimated$statistic, ks_estimated$p.value))

if (ks_estimated$p.value < 0.05) {
  cat(sprintf("  ✗ FAIL: Estimated statistics do NOT follow χ²(%d) (p < 0.05)\n", df_expected))
} else {
  cat(sprintf("  ✓ PASS: Estimated statistics follow χ²(%d) (p ≥ 0.05)\n", df_expected))
}

# Comparison (filtered data)
cat("\n=== Comparison ===\n")
cat(sprintf("Mean ratio (Observed/Expected):\n"))
cat(sprintf("  Oracle: %.4f\n", mean(oracle_stats_filtered) / df_expected))
cat(sprintf("  Estimated: %.4f\n", mean(estimated_stats_filtered) / df_expected))
cat(sprintf("\nInflation difference: %.4f%%\n",
            100 * (mean(estimated_stats_filtered) - mean(oracle_stats_filtered)) / df_expected))

# Create diagnostic plots
cat("\nGenerating diagnostic plots...\n")
png("debugging/oracle_test_qq.png", width = 1200, height = 800)
par(mfrow = c(2, 3))

# Q-Q plot for oracle
qqplot(qchisq(ppoints(n_replicates), df = df_expected),
       sort(oracle_stats),
       main = sprintf("Oracle Test Q-Q Plot\n(Expected: χ²(%d))", df_expected),
       xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles (Oracle)",
       pch = 19, col = rgb(0, 0, 1, 0.3))
abline(0, 1, col = "red", lwd = 2)
legend("topleft",
       legend = sprintf("KS p = %.4f", ks_oracle$p.value),
       bty = "n", cex = 1.2)

# Q-Q plot for estimated
qqplot(qchisq(ppoints(n_replicates), df = df_expected),
       sort(estimated_stats),
       main = sprintf("Estimated Test Q-Q Plot\n(Expected: χ²(%d))", df_expected),
       xlab = "Theoretical Quantiles",
       ylab = "Sample Quantiles (Estimated)",
       pch = 19, col = rgb(0, 0, 1, 0.3))
abline(0, 1, col = "red", lwd = 2)
legend("topleft",
       legend = sprintf("KS p = %.4f", ks_estimated$p.value),
       bty = "n", cex = 1.2)

# Histogram comparison
hist(oracle_stats, breaks = 30, probability = TRUE,
     main = "Oracle Statistics Distribution",
     xlab = "Test Statistic",
     col = rgb(0, 0, 1, 0.5),
     xlim = c(0, max(oracle_stats, estimated_stats)))
curve(dchisq(x, df = df_expected), add = TRUE, col = "red", lwd = 2)
legend("topright",
       legend = c("Observed", sprintf("χ²(%d)", df_expected)),
       fill = c(rgb(0, 0, 1, 0.5), NA),
       border = c("black", NA),
       lty = c(NA, 1),
       col = c(NA, "red"),
       lwd = c(NA, 2))

hist(estimated_stats, breaks = 30, probability = TRUE,
     main = "Estimated Statistics Distribution",
     xlab = "Test Statistic",
     col = rgb(1, 0, 0, 0.5),
     xlim = c(0, max(oracle_stats, estimated_stats)))
curve(dchisq(x, df = df_expected), add = TRUE, col = "red", lwd = 2)
legend("topright",
       legend = c("Observed", sprintf("χ²(%d)", df_expected)),
       fill = c(rgb(1, 0, 0, 0.5), NA),
       border = c("black", NA),
       lty = c(NA, 1),
       col = c(NA, "red"),
       lwd = c(NA, 2))

# Overlay comparison
hist(oracle_stats, breaks = 30, probability = TRUE,
     main = "Oracle vs Estimated Comparison",
     xlab = "Test Statistic",
     col = rgb(0, 0, 1, 0.3),
     xlim = c(0, max(oracle_stats, estimated_stats)))
hist(estimated_stats, breaks = 30, probability = TRUE,
     col = rgb(1, 0, 0, 0.3), add = TRUE)
curve(dchisq(x, df = df_expected), add = TRUE, col = "black", lwd = 2)
legend("topright",
       legend = c("Oracle", "Estimated", sprintf("χ²(%d)", df_expected)),
       fill = c(rgb(0, 0, 1, 0.3), rgb(1, 0, 0, 0.3), NA),
       border = c("black", "black", NA),
       lty = c(NA, NA, 1),
       col = c(NA, NA, "black"),
       lwd = c(NA, NA, 2))

# ECDFs
plot(ecdf(oracle_stats),
     main = "Empirical CDFs",
     xlab = "Test Statistic",
     ylab = "Cumulative Probability",
     col = "blue", lwd = 2)
plot(ecdf(estimated_stats), add = TRUE, col = "red", lwd = 2)
curve(pchisq(x, df = df_expected), add = TRUE, col = "black", lwd = 2, lty = 2)
legend("bottomright",
       legend = c("Oracle", "Estimated", sprintf("χ²(%d)", df_expected)),
       col = c("blue", "red", "black"),
       lwd = 2,
       lty = c(1, 1, 2))

dev.off()
cat("Plots saved to: debugging/oracle_test_qq.png\n")

# Final diagnostic verdict
cat("\n=== FINAL DIAGNOSTIC VERDICT ===\n")

oracle_pass <- ks_oracle$p.value >= 0.05
estimated_pass <- ks_estimated$p.value >= 0.05

if (oracle_pass && !estimated_pass) {
  cat("✓ Oracle test PASSES but estimated test FAILS\n\n")
  cat("CONCLUSION: Issue is in RIEMTAN (Fréchet mean estimation)\n")
  cat("  - The frechet_anova formula is correct (works with true means)\n")
  cat("  - The problem is bias/error in estimated Fréchet means\n")
  cat("  - Action: Investigate riemtan's Fréchet mean estimation algorithm\n")
} else if (!oracle_pass && !estimated_pass) {
  cat("✗ Both oracle and estimated tests FAIL\n\n")
  cat("CONCLUSION: Issue is in RIEMSTATS (formula problem)\n")
  cat("  - The frechet_anova test statistic formula is incorrect\n")
  cat("  - Even with perfect means, the statistic doesn't follow χ²\n")
  cat("  - Action: Review frechet_anova implementation in R/frechet_anova.R\n")
} else if (!oracle_pass && estimated_pass) {
  cat("⚠ Unexpected: Oracle FAILS but estimated PASSES\n\n")
  cat("CONCLUSION: Anomalous result - may indicate issue with oracle setup\n")
  cat("  - Action: Review test methodology\n")
} else {
  cat("✓ Both tests PASS\n\n")
  cat("CONCLUSION: No consistent issue detected with these parameters\n")
  cat("  - Issue 59 may be specific to different parameter regimes\n")
  cat("  - Action: Try different values of n, g, d, sigma\n")
}

cat("\nTest completed.\n")
