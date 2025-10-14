#!/usr/bin/env Rscript
# Decomposition Analysis for Fréchet ANOVA
# Purpose: Verify asymptotic behavior - Term1 should → χ²(k-1), Term2 should → 0

library(riemtan)
library(riemstats)
library(Matrix)
library(parallel)

cat("=== Fréchet ANOVA Asymptotic Decomposition Test ===\n")
cat("Testing Term1 (Un) and Term2 (F_stat²) separately\n")
cat("Theory: Term1 → χ²(k-1), Term2 → 0 as n → ∞\n\n")

# Set Cairo for headless plotting (no X11 required)
options(bitmapType = 'cairo')

# Detect number of cores
n_cores <- detectCores()
cat(sprintf("Using %d CPU cores for parallel processing\n\n", n_cores))

# Load metric
data(euclidean)

# Set seed for reproducibility
set.seed(42)

# Test parameters
n_per_group <- 30  # samples per group
g <- 2             # number of groups
d <- 3             # dimension of SPD matrices (3×3)
sigma <- 0.01      # dispersion parameter
n_replicates <- 50 # number of test runs

cat(sprintf("Parameters:\n"))
cat(sprintf("  Groups: %d\n", g))
cat(sprintf("  Samples per group: %d\n", n_per_group))
cat(sprintf("  Matrix dimension: %d×%d\n", d, d))
cat(sprintf("  Dispersion: %.2f\n", sigma))
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

# Helper function to set Fréchet mean by accessing private field
set_frechet_mean <- function(sample_obj, mean_value) {
  private_env <- environment(sample_obj$initialize)$private
  private_env$f_mean <- mean_value
  invisible(sample_obj)
}

# MODIFIED frechet_anova that returns detailed decomposition
frechet_anova_detailed <- function(super_sample) {
  if (!inherits(super_sample, "CSuperSample")) {
    stop("Argument 'super_sample' must be an object of class 'CSuperSample'.")
  }

  if (length(super_sample$list_of_samples) < 2) {
    stop("CSuperSample must contain at least 2 groups for ANOVA analysis")
  }

  k <- super_sample$list_of_samples |> length()
  n <- super_sample$sample_size

  group_stats <- super_sample$list_of_samples |>
    purrr::map(\(sample) {
      sample$compute_dists()
      sample$distances
    }) |>
    purrr::map(\(vec_of_dists) {
      group_v <- mean(vec_of_dists^2)
      group_sig_2_raw <- mean(vec_of_dists^4) - mean(vec_of_dists^2)^2

      # Cap group_sig_2 to prevent numerical instability
      epsilon <- 0.1
      group_sig_2 <- max(group_sig_2_raw, epsilon * group_v^2)

      list(
        group_sig_2 = group_sig_2,
        group_sig_2_raw = group_sig_2_raw,
        group_v = group_v,
        group_sample_size = length(vec_of_dists)
      )
    })

  group_sig_2s <- group_stats |> purrr::map_dbl(\(stats) stats$group_sig_2)
  group_sig_2s_raw <- group_stats |> purrr::map_dbl(\(stats) stats$group_sig_2_raw)
  group_vs <- group_stats |> purrr::map_dbl(\(stats) stats$group_v)

  super_sample$gather()
  super_sample$full_sample$compute_dists()
  pooled_v <- mean(super_sample$full_sample$distances^2)

  group_rel_sizes <- group_stats |>
    purrr::map_dbl(\(stats) stats$group_sample_size) |>
    (\(x) x / super_sample$full_sample$sample_size)()

  F_stat <- pooled_v - sum(group_rel_sizes * group_vs)

  # Vectorized computation of Un
  idx_mat <- which(upper.tri(matrix(0, k, k)), arr.ind = TRUE)
  j_idx <- idx_mat[, 1]
  l_idx <- idx_mat[, 2]
  Un <- sum(
    (group_rel_sizes[j_idx] * group_rel_sizes[l_idx]) /
      (group_sig_2s[j_idx] * group_sig_2s[l_idx]) *
      (group_vs[j_idx] - group_vs[l_idx])^2
  )

  # Denominators
  denom1 <- sum(group_rel_sizes / group_sig_2s)
  denom2 <- sum((group_rel_sizes^2) * group_sig_2s)

  # Terms
  term1 <- (n * Un) / denom1
  term2 <- (n * (F_stat^2)) / denom2
  thestat <- term1 + term2

  # Compute p-value
  p_value <- stats::pchisq(thestat, df = (k - 1), lower.tail = FALSE)

  list(
    statistic = thestat,
    p_value = p_value,
    term1 = term1,
    term2 = term2,
    Un = Un,
    F_stat = F_stat,
    denom1 = denom1,
    denom2 = denom2,
    group_sig_2s = group_sig_2s,
    group_sig_2s_raw = group_sig_2s_raw,
    group_vs = group_vs,
    pooled_v = pooled_v,
    n_capped = sum(group_sig_2s != group_sig_2s_raw)
  )
}

cat("Generating data under H0 (all groups from same distribution)...\n")
cat("Running", n_replicates, "replicates in parallel...\n\n")

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

    # Oracle test with TRUE mean
    oracle_groups <- lapply(groups, function(sample) {
      oracle_sample <- sample$clone()
      set_frechet_mean(oracle_sample, true_center)
      oracle_sample
    })

    oracle_ss <- CSuperSample$new(oracle_groups)
    oracle_result <- frechet_anova_detailed(oracle_ss)

    # Standard test with ESTIMATED means
    estimated_ss <- CSuperSample$new(groups)
    estimated_result <- frechet_anova_detailed(estimated_ss)

    # Return detailed results
    list(
      oracle = oracle_result,
      estimated = estimated_result
    )
  }, error = function(e) {
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
  results <- results[!is_error]
}

if (n_success == 0) {
  stop("All simulations failed. Cannot proceed with analysis.")
}

cat(sprintf("Successfully completed %d simulations\n\n", n_success))

# Extract results
oracle_stats <- sapply(results, \(x) x$oracle$statistic)
oracle_term1 <- sapply(results, \(x) x$oracle$term1)
oracle_term2 <- sapply(results, \(x) x$oracle$term2)
oracle_Un <- sapply(results, \(x) x$oracle$Un)
oracle_F_stat <- sapply(results, \(x) x$oracle$F_stat)
oracle_n_capped <- sapply(results, \(x) x$oracle$n_capped)

estimated_stats <- sapply(results, \(x) x$estimated$statistic)
estimated_term1 <- sapply(results, \(x) x$estimated$term1)
estimated_term2 <- sapply(results, \(x) x$estimated$term2)
estimated_Un <- sapply(results, \(x) x$estimated$Un)
estimated_F_stat <- sapply(results, \(x) x$estimated$F_stat)
estimated_n_capped <- sapply(results, \(x) x$estimated$n_capped)

# Expected degrees of freedom
df_expected <- g - 1

cat("=== ASYMPTOTIC DECOMPOSITION ANALYSIS ===\n\n")

# ============= ORACLE RESULTS =============
cat("--- Oracle Test (TRUE means) ---\n\n")

cat("Term1 (Un-based, should → χ²(k-1)):\n")
cat(sprintf("  Mean: %.4f (Expected: %.4f)\n", mean(oracle_term1), df_expected))
cat(sprintf("  Variance: %.4f (Expected: %.4f)\n", var(oracle_term1), 2*df_expected))
cat(sprintf("  Median: %.4f (Expected: %.4f)\n", median(oracle_term1), qchisq(0.5, df_expected)))
cat(sprintf("  Range: [%.4f, %.4f]\n", min(oracle_term1), max(oracle_term1)))

ks_oracle_term1 <- ks.test(oracle_term1, "pchisq", df = df_expected)
cat(sprintf("  KS test: D = %.4f, p-value = %.6f ", ks_oracle_term1$statistic, ks_oracle_term1$p.value))
if (ks_oracle_term1$p.value >= 0.05) {
  cat("✓ PASS\n")
} else {
  cat("✗ FAIL\n")
}

cat("\nTerm2 (F_stat²-based, should → 0):\n")
cat(sprintf("  Mean: %.6f (Should be ≈ 0)\n", mean(oracle_term2)))
cat(sprintf("  Variance: %.6f\n", var(oracle_term2)))
cat(sprintf("  Median: %.6f\n", median(oracle_term2)))
cat(sprintf("  Range: [%.6f, %.6f]\n", min(oracle_term2), max(oracle_term2)))
cat(sprintf("  Max/Mean(Term1) ratio: %.4f%% (should be small)\n",
            100 * max(oracle_term2) / mean(oracle_term1)))

cat("\nTotal statistic:\n")
cat(sprintf("  Mean: %.4f (Expected: %.4f)\n", mean(oracle_stats), df_expected))
cat(sprintf("  Variance: %.4f (Expected: %.4f)\n", var(oracle_stats), 2*df_expected))

ks_oracle_total <- ks.test(oracle_stats, "pchisq", df = df_expected)
cat(sprintf("  KS test: D = %.4f, p-value = %.6f ", ks_oracle_total$statistic, ks_oracle_total$p.value))
if (ks_oracle_total$p.value >= 0.05) {
  cat("✓ PASS\n")
} else {
  cat("✗ FAIL\n")
}

cat("\nContribution analysis:\n")
cat(sprintf("  Mean Term1 / Total: %.2f%%\n", 100 * mean(oracle_term1) / mean(oracle_stats)))
cat(sprintf("  Mean Term2 / Total: %.2f%%\n", 100 * mean(oracle_term2) / mean(oracle_stats)))
cat(sprintf("  Variance capping: %.1f%% of simulations\n", 100 * mean(oracle_n_capped > 0)))

cat("\nIntermediate quantities:\n")
cat(sprintf("  Un - Mean: %.6f, Variance: %.6f\n", mean(oracle_Un), var(oracle_Un)))
cat(sprintf("  F_stat - Mean: %.6f, Variance: %.6f\n", mean(oracle_F_stat), var(oracle_F_stat)))

# ============= ESTIMATED RESULTS =============
cat("\n--- Estimated Test (ESTIMATED means) ---\n\n")

cat("Term1 (Un-based, should → χ²(k-1)):\n")
cat(sprintf("  Mean: %.4f (Expected: %.4f)\n", mean(estimated_term1), df_expected))
cat(sprintf("  Variance: %.4f (Expected: %.4f)\n", var(estimated_term1), 2*df_expected))
cat(sprintf("  Median: %.4f (Expected: %.4f)\n", median(estimated_term1), qchisq(0.5, df_expected)))
cat(sprintf("  Range: [%.4f, %.4f]\n", min(estimated_term1), max(estimated_term1)))

ks_estimated_term1 <- ks.test(estimated_term1, "pchisq", df = df_expected)
cat(sprintf("  KS test: D = %.4f, p-value = %.6f ", ks_estimated_term1$statistic, ks_estimated_term1$p.value))
if (ks_estimated_term1$p.value >= 0.05) {
  cat("✓ PASS\n")
} else {
  cat("✗ FAIL\n")
}

cat("\nTerm2 (F_stat²-based, should → 0):\n")
cat(sprintf("  Mean: %.6f (Should be ≈ 0)\n", mean(estimated_term2)))
cat(sprintf("  Variance: %.6f\n", var(estimated_term2)))
cat(sprintf("  Median: %.6f\n", median(estimated_term2)))
cat(sprintf("  Range: [%.6f, %.6f]\n", min(estimated_term2), max(estimated_term2)))
cat(sprintf("  Max/Mean(Term1) ratio: %.4f%% (should be small)\n",
            100 * max(estimated_term2) / mean(estimated_term1)))

cat("\nTotal statistic:\n")
cat(sprintf("  Mean: %.4f (Expected: %.4f)\n", mean(estimated_stats), df_expected))
cat(sprintf("  Variance: %.4f (Expected: %.4f)\n", var(estimated_stats), 2*df_expected))

ks_estimated_total <- ks.test(estimated_stats, "pchisq", df = df_expected)
cat(sprintf("  KS test: D = %.4f, p-value = %.6f ", ks_estimated_total$statistic, ks_estimated_total$p.value))
if (ks_estimated_total$p.value >= 0.05) {
  cat("✓ PASS\n")
} else {
  cat("✗ FAIL\n")
}

cat("\nContribution analysis:\n")
cat(sprintf("  Mean Term1 / Total: %.2f%%\n", 100 * mean(estimated_term1) / mean(estimated_stats)))
cat(sprintf("  Mean Term2 / Total: %.2f%%\n", 100 * mean(estimated_term2) / mean(estimated_stats)))
cat(sprintf("  Variance capping: %.1f%% of simulations\n", 100 * mean(estimated_n_capped > 0)))

cat("\nIntermediate quantities:\n")
cat(sprintf("  Un - Mean: %.6f, Variance: %.6f\n", mean(estimated_Un), var(estimated_Un)))
cat(sprintf("  F_stat - Mean: %.6f, Variance: %.6f\n", mean(estimated_F_stat), var(estimated_F_stat)))

# ============= DIAGNOSTIC PLOTS =============
cat("\nGenerating diagnostic plots...\n")
png("debugging/decomposition_diagnostics.png", width = 1600, height = 1200)
par(mfrow = c(3, 4), mar = c(4, 4, 3, 1))

# Oracle plots
# Term1 distribution
hist(oracle_term1, breaks = 30, probability = TRUE,
     main = "Oracle: Term1 Distribution",
     xlab = "Term1 (Un-based)",
     col = rgb(0, 0, 1, 0.5))
curve(dchisq(x, df = df_expected), add = TRUE, col = "red", lwd = 2)
legend("topright", legend = sprintf("KS p=%.3f", ks_oracle_term1$p.value), bty = "n")

# Term1 Q-Q plot
qqplot(qchisq(ppoints(n_success), df = df_expected), oracle_term1,
       main = "Oracle: Term1 Q-Q Plot",
       xlab = sprintf("χ²(%d) Quantiles", df_expected),
       ylab = "Term1 Quantiles",
       pch = 19, col = rgb(0, 0, 1, 0.5))
abline(0, 1, col = "red", lwd = 2)

# Term2 distribution
hist(oracle_term2, breaks = 30, probability = TRUE,
     main = "Oracle: Term2 Distribution",
     xlab = "Term2 (F_stat²-based)",
     col = rgb(1, 0, 0, 0.5))
abline(v = 0, col = "blue", lwd = 2, lty = 2)
legend("topright", legend = sprintf("Mean=%.4f", mean(oracle_term2)), bty = "n")

# Term1 vs Term2 scatter
plot(oracle_term1, oracle_term2,
     main = "Oracle: Term1 vs Term2",
     xlab = "Term1", ylab = "Term2",
     pch = 19, col = rgb(0, 0, 1, 0.3))
abline(h = 0, col = "red", lty = 2)

# Estimated plots
# Term1 distribution
hist(estimated_term1, breaks = 30, probability = TRUE,
     main = "Estimated: Term1 Distribution",
     xlab = "Term1 (Un-based)",
     col = rgb(0, 0, 1, 0.5))
curve(dchisq(x, df = df_expected), add = TRUE, col = "red", lwd = 2)
legend("topright", legend = sprintf("KS p=%.3f", ks_estimated_term1$p.value), bty = "n")

# Term1 Q-Q plot
qqplot(qchisq(ppoints(n_success), df = df_expected), estimated_term1,
       main = "Estimated: Term1 Q-Q Plot",
       xlab = sprintf("χ²(%d) Quantiles", df_expected),
       ylab = "Term1 Quantiles",
       pch = 19, col = rgb(0, 0, 1, 0.5))
abline(0, 1, col = "red", lwd = 2)

# Term2 distribution
hist(estimated_term2, breaks = 30, probability = TRUE,
     main = "Estimated: Term2 Distribution",
     xlab = "Term2 (F_stat²-based)",
     col = rgb(1, 0, 0, 0.5))
abline(v = 0, col = "blue", lwd = 2, lty = 2)
legend("topright", legend = sprintf("Mean=%.4f", mean(estimated_term2)), bty = "n")

# Term1 vs Term2 scatter
plot(estimated_term1, estimated_term2,
     main = "Estimated: Term1 vs Term2",
     xlab = "Term1", ylab = "Term2",
     pch = 19, col = rgb(0, 0, 1, 0.3))
abline(h = 0, col = "red", lty = 2)

# Comparison plots
# Term1 comparison
boxplot(list(Oracle = oracle_term1, Estimated = estimated_term1),
        main = "Term1 Comparison",
        ylab = "Term1 value",
        col = c(rgb(0, 0, 1, 0.3), rgb(1, 0, 0, 0.3)))
abline(h = df_expected, col = "darkgreen", lwd = 2, lty = 2)

# Term2 comparison
boxplot(list(Oracle = oracle_term2, Estimated = estimated_term2),
        main = "Term2 Comparison",
        ylab = "Term2 value",
        col = c(rgb(0, 0, 1, 0.3), rgb(1, 0, 0, 0.3)))
abline(h = 0, col = "darkgreen", lwd = 2, lty = 2)

# Contribution pie chart - Oracle
pie_oracle <- c(mean(oracle_term1), mean(oracle_term2))
pie(pie_oracle, labels = c(sprintf("Term1: %.2f%%", 100*pie_oracle[1]/sum(pie_oracle)),
                            sprintf("Term2: %.2f%%", 100*pie_oracle[2]/sum(pie_oracle))),
    main = "Oracle: Mean Contribution",
    col = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)))

# Contribution pie chart - Estimated
pie_estimated <- c(mean(estimated_term1), mean(estimated_term2))
pie(pie_estimated, labels = c(sprintf("Term1: %.2f%%", 100*pie_estimated[1]/sum(pie_estimated)),
                               sprintf("Term2: %.2f%%", 100*pie_estimated[2]/sum(pie_estimated))),
    main = "Estimated: Mean Contribution",
    col = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)))

dev.off()
cat("Plots saved to: debugging/decomposition_diagnostics.png\n\n")

# ============= DIAGNOSTIC VERDICT =============
cat("=== DIAGNOSTIC VERDICT ===\n\n")

# Check if Term1 follows chi-squared
term1_oracle_ok <- ks_oracle_term1$p.value >= 0.05
term1_estimated_ok <- ks_estimated_term1$p.value >= 0.05

# Check if Term2 is vanishing (mean should be << 1)
term2_oracle_vanishing <- mean(oracle_term2) < 0.1 * df_expected
term2_estimated_vanishing <- mean(estimated_term2) < 0.1 * df_expected

cat("Term1 (Un) Analysis:\n")
if (term1_oracle_ok && term1_estimated_ok) {
  cat("  ✓ Term1 follows χ²(k-1) for both oracle and estimated\n")
} else if (!term1_oracle_ok) {
  cat("  ✗ Term1 does NOT follow χ²(k-1) even with true means (FORMULA ISSUE)\n")
} else {
  cat("  ⚠ Term1 follows χ²(k-1) with oracle but not with estimation\n")
}

cat("\nTerm2 (F_stat²) Analysis:\n")
if (term2_oracle_vanishing && term2_estimated_vanishing) {
  cat("  ✓ Term2 is vanishing (→ 0) for both oracle and estimated\n")
} else if (!term2_oracle_vanishing) {
  cat("  ✗ Term2 is NOT vanishing even with true means (FORMULA ISSUE)\n")
  cat(sprintf("    Oracle Term2 mean = %.6f (should be ≈ 0)\n", mean(oracle_term2)))
} else {
  cat("  ⚠ Term2 vanishes with oracle but not with estimation\n")
  cat(sprintf("    Estimated Term2 mean = %.6f (should be ≈ 0)\n", mean(estimated_term2)))
}

cat("\nOverall Conclusion:\n")
if (term1_oracle_ok && term2_oracle_vanishing) {
  cat("  ✓ Asymptotic theory holds with oracle (true means)\n")
  if (!term1_estimated_ok || !term2_estimated_vanishing) {
    cat("  → Issue is in ESTIMATION, not formula\n")
  }
} else {
  cat("  ✗ Asymptotic theory FAILS even with oracle (true means)\n")
  cat("  → Issue is in FORMULA (riemstats implementation)\n")
  if (!term1_oracle_ok) {
    cat("    Problem: Term1 formula is incorrect\n")
  }
  if (!term2_oracle_vanishing) {
    cat("    Problem: Term2 is not vanishing as theory requires\n")
  }
}

cat("\nTest completed.\n")
