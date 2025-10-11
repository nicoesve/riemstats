#!/usr/bin/env Rscript
# Test 2: Distance Computation Verification
# Purpose: Verify that riemtan's compute_dists() produces accurate distances
# from points to their Fréchet mean

library(riemtan)
library(Matrix)
library(parallel)

cat("=== Distance Computation Verification Test ===\n\n")

# Detect number of cores
n_cores <- detectCores()
cat(sprintf("Using %d CPU cores for parallel processing\n\n", n_cores))

# Load AIRM metric
data(airm)

# Test parameters
d <- 3  # dimension of SPD matrices
n_samples <- 50  # number of samples per test
n_tests <- 100  # number of test runs
tolerance <- 0.05  # 5% relative error tolerance

# Helper function to set Fréchet mean by accessing private field
set_frechet_mean <- function(sample_obj, mean_value) {
  # Access the private environment of the R6 object
  private_env <- environment(sample_obj$initialize)$private
  # Set the private f_mean field
  private_env$f_mean <- mean_value
  invisible(sample_obj)
}

# Helper function: Create SPD matrix at specific distance from identity using AIRM
# Distance formula: d(I, A) = ||log(A)||_F
create_spd_at_distance <- function(dim, target_distance) {
  # Create a random symmetric matrix
  V <- matrix(rnorm(dim * dim), dim, dim)
  V <- (V + t(V)) / 2

  # Normalize to have Frobenius norm = target_distance
  V <- V * target_distance / norm(V, "F")

  # Exponentiate to get SPD matrix
  A <- expm::expm(V)

  # Convert to dppMatrix format
  A_pd <- as(A, "dpoMatrix")
  A_packed <- Matrix::pack(A_pd)

  return(A_packed)
}

# Test 1: Verify distances for samples at known distances from identity
cat("Test 1: Samples at known distances from identity\n")
cat("-----------------------------------------------\n")

test_distances <- c(0.5, 1.0, 2.0, 3.0)
errors <- numeric(length(test_distances))

for (i in seq_along(test_distances)) {
  target_dist <- test_distances[i]

  # Generate samples at this distance
  samples <- lapply(1:n_samples, function(j) {
    create_spd_at_distance(d, target_dist)
  })

  # Create CSample with identity as center
  identity_mat <- diag(d) |> as("dpoMatrix") |> Matrix::pack()
  sample_obj <- CSample$new(samples, metric_obj = airm)

  # Set Fréchet mean to identity (oracle)
  set_frechet_mean(sample_obj, identity_mat)

  # Compute distances
  sample_obj$compute_dists()

  # Get computed distances (these are squared distances)
  computed_dists_sq <- sample_obj$distances
  computed_dists <- sqrt(computed_dists_sq)

  # Calculate error
  mean_computed <- mean(computed_dists)
  relative_error <- abs(mean_computed - target_dist) / target_dist
  errors[i] <- relative_error

  cat(sprintf("  Target distance: %.2f\n", target_dist))
  cat(sprintf("  Mean computed: %.4f (SD: %.4f)\n", mean_computed, sd(computed_dists)))
  cat(sprintf("  Relative error: %.2f%%\n", relative_error * 100))

  if (relative_error > tolerance) {
    cat("  *** WARNING: Error exceeds tolerance! ***\n")
  } else {
    cat("  ✓ PASS\n")
  }
  cat("\n")
}

# Test 2: Distribution test - generate null data and check distance distribution
cat("\nTest 2: Null distribution distance verification\n")
cat("-----------------------------------------------\n")

sigma <- 100.0  # dispersion parameter
center <- diag(d) |> as("dpoMatrix") |> Matrix::pack()

# Create dispersion matrix
p <- d * (d + 1) / 2
scale <- sigma * diag(p) |> as("dpoMatrix") |> Matrix::pack()

# Load log_cholesky metric for generating data
data(log_cholesky)

# Generate multiple samples and compute distances in parallel
cat("Generating", n_tests, "samples in parallel...\n")

get_distances <- function(iter) {
  # Generate sample from Riemannian normal
  sample_list <- rspdnorm(n_samples, center, scale, log_cholesky)
  sample_list$compute_unvecs()
  sample_list$compute_conns()

  # Compute distances (this uses estimated Fréchet mean)
  sample_list$compute_dists()

  # Return squared distances
  sample_list$distances
}

all_distances_list <- mclapply(1:n_tests, get_distances, mc.cores = n_cores)
all_distances <- unlist(all_distances_list)

# Convert squared distances to distances
all_distances_unsquared <- sqrt(all_distances)

# Summary statistics
cat(sprintf("Number of samples: %d\n", length(all_distances)))
cat(sprintf("Mean squared distance: %.4f\n", mean(all_distances)))
cat(sprintf("SD squared distance: %.4f\n", sd(all_distances)))
cat(sprintf("Mean distance: %.4f\n", mean(all_distances_unsquared)))
cat(sprintf("SD distance: %.4f\n", sd(all_distances_unsquared)))

# For Riemannian normal, E[d²] is related to variance
# This is a consistency check, not a precise theoretical test
expected_mean_sq_dist <- sigma * p / (d * (d + 1))  # rough approximation
cat(sprintf("\nApproximate expected mean squared distance: %.4f\n", expected_mean_sq_dist))

# Test 3: Check consistency - same data, recompute distances
cat("\nTest 3: Consistency check\n")
cat("-------------------------\n")

# Generate one sample
test_sample <- rspdnorm(n_samples, center, scale, log_cholesky)
test_sample$compute_unvecs()
test_sample$compute_conns()

# Compute distances first time
test_sample$compute_dists()
dists1 <- test_sample$distances

# Compute distances second time
test_sample$compute_dists()
dists2 <- test_sample$distances

# Check if identical
max_diff <- max(abs(dists1 - dists2))
cat(sprintf("Max difference between recomputations: %.10e\n", max_diff))

if (max_diff < 1e-10) {
  cat("✓ PASS: Distance computation is consistent\n")
} else {
  cat("*** WARNING: Distance computation is not consistent! ***\n")
}

# Create Q-Q plot comparing theoretical vs computed for Test 1
cat("\n\nGenerating diagnostic plot...\n")
png("debugging/distance_verification_qq.png", width = 800, height = 600)
par(mfrow = c(2, 2))

# Plot for each test distance
for (i in seq_along(test_distances)) {
  target_dist <- test_distances[i]

  # Regenerate samples for plotting
  samples <- lapply(1:n_samples, function(j) {
    create_spd_at_distance(d, target_dist)
  })

  identity_mat <- diag(d) |> as("dpoMatrix") |> Matrix::pack()
  sample_obj <- CSample$new(samples, metric_obj = airm)
  set_frechet_mean(sample_obj, identity_mat)
  sample_obj$compute_dists()
  computed_dists <- sqrt(sample_obj$distances)

  # Simple diagnostic plot
  hist(computed_dists,
       main = sprintf("Distance Distribution\n(Target: %.1f)", target_dist),
       xlab = "Computed Distance",
       breaks = 20,
       col = "lightblue")
  abline(v = target_dist, col = "red", lwd = 2, lty = 2)
  abline(v = mean(computed_dists), col = "blue", lwd = 2)
  legend("topright",
         legend = c("Target", "Observed Mean"),
         col = c("red", "blue"),
         lty = c(2, 1), lwd = 2,
         cex = 0.8)
}

dev.off()
cat("Plot saved to: debugging/distance_verification_qq.png\n")

# Final verdict
cat("\n=== FINAL VERDICT ===\n")
max_error <- max(errors)
cat(sprintf("Maximum relative error: %.2f%%\n", max_error * 100))

if (max_error < tolerance) {
  cat("\n✓ PASS: Distance computation appears accurate\n")
  cat("Conclusion: riemtan distance computation is likely NOT the source of issue 59\n")
} else {
  cat("\n✗ FAIL: Distance computation has significant errors\n")
  cat("Conclusion: riemtan distance computation may be contributing to issue 59\n")
}

cat("\nTest completed.\n")
