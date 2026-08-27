#' Compute the Log Wilks' Lambda Statistic
#'
#' Calculates the log of Wilks' lambda statistic
#' for a given \code{super_sample} object.
#' This function ensures that the within-group
#' and total covariance matrices are computed,
#' then computes the difference of their log determinants.
#'
#' @param super_sample An object of class CSuperSample
#'
#' @return A numeric value representing the log Wilks' lambda statistic.
#'
#' @details
#' Wilks' lambda is a test statistic for the ANOVA test decribed in (to appear)
#' @export
log_wilks_lambda <- function(super_sample) {
  if (!inherits(super_sample, "CSuperSample")) {
    stop("Argument 'super_sample' must be an object of class 'CSuperSample'.")
  }

  # Check for minimum number of groups
  if (length(super_sample$list_of_samples) < 2) {
    stop("CSuperSample must contain at least 2 groups for ANOVA analysis")
  }

  # Compute within-group covariance matrix with error handling
  if (super_sample$Within |> is.null()) {
    tryCatch({
      super_sample$compute_W()
    }, error = function(e) {
      stop("Failed to compute within-group covariance matrix: ", e$message)
    })
  }

  # Compute total covariance matrix with error handling
  if (super_sample$Total |> is.null()) {
    tryCatch({
      super_sample$compute_T()
    }, error = function(e) {
      if (grepl("In index:", e$message)) {
        stop("Failed to compute total covariance matrix due to CSuperSample structure issues. ",
             "This may be caused by incompatible sample data or a bug in the riemtan package. ",
             "Try creating CSample objects first before combining into CSuperSample. ",
             "Original error: ", e$message)
      } else {
        stop("Failed to compute total covariance matrix: ", e$message)
      }
    })
  }

  # Check that matrices are valid
  if (is.null(super_sample$Within)) {
    stop("Within-group covariance matrix is NULL after computation")
  }
  
  if (is.null(super_sample$Total)) {
    stop("Total covariance matrix is NULL after computation")
  }

  # Compute determinants with error handling
  tryCatch({
    det_within <- Matrix::determinant(super_sample$Within)$modulus
    det_total <- Matrix::determinant(super_sample$Total)$modulus
    
    result <- det_within - det_total
    result
  }, error = function(e) {
    stop("Failed to compute log Wilks' lambda determinants: ", e$message)
  })
}

#' Compute Pillai's Trace Statistic
#'
#' Calculates Pillai's trace statistic for a given \code{super_sample} object.
#' This function ensures that the within-group
#' and total covariance matrices are computed,
#' then computes the sum of the eigenvalues of
#' the matrix (Total - Within) %*% solve(Total).
#'
#' @param super_sample An object of class CSuperSample
#'
#' @return A numeric value representing Pillai's trace statistic.
#'
#' @details
#' Pillai's trace is a test statistic for
#' the ANOVA test decribed in (to appear).
#' @export
pillais_trace <- function(super_sample) {
  if (!inherits(super_sample, "CSuperSample")) {
    stop("Argument 'super_sample' must be an object of class 'CSuperSample'.")
  }

  # Check for minimum number of groups
  if (length(super_sample$list_of_samples) < 2) {
    stop("CSuperSample must contain at least 2 groups for ANOVA analysis")
  }

  # Compute within-group covariance matrix with error handling
  if (super_sample$Within |> is.null()) {
    tryCatch({
      super_sample$compute_W()
    }, error = function(e) {
      stop("Failed to compute within-group covariance matrix: ", e$message)
    })
  }

  # Compute total covariance matrix with error handling
  if (super_sample$Total |> is.null()) {
    tryCatch({
      super_sample$compute_T()
    }, error = function(e) {
      if (grepl("In index:", e$message)) {
        stop("Failed to compute total covariance matrix due to CSuperSample structure issues. ",
             "This may be caused by incompatible sample data or a bug in the riemtan package. ",
             "Try creating CSample objects first before combining into CSuperSample. ",
             "Original error: ", e$message)
      } else {
        stop("Failed to compute total covariance matrix: ", e$message)
      }
    })
  }

  # Check that matrices are valid
  if (is.null(super_sample$Within)) {
    stop("Within-group covariance matrix is NULL after computation")
  }
  
  if (is.null(super_sample$Total)) {
    stop("Total covariance matrix is NULL after computation")
  }

  # Compute Pillai's trace with error handling
  tryCatch({
    result <- (
      (super_sample$Total - super_sample$Within) %*%
        solve(super_sample$Total)
    ) |>
      as.matrix() |>
      diag() |>
      sum()

    result
  }, error = function(e) {
    stop("Failed to compute Pillai's trace: ", e$message)
  })
}

#' Permutation Statistic for a Super Sample
#'
#' Randomly shuffles all observations across groups while preserving
#' group sample sizes, creating a new super sample under the null hypothesis
#' of no group differences, and computes a specified statistic on the
#' resulting permuted super sample.
#'
#' @param x An object of class \code{CSuperSample},
#' representing the original super sample.
#' @param stat_fun A function to compute a statistic
#' on the resulting \code{CSuperSample} object.
#'
#' @return The value returned by \code{stat_fun} when applied
#' to the permuted super sample.
#'
#' @details
#' This function performs a permutation test by:
#' 1. Extracting all data points from all groups
#' 2. Randomly shuffling the data
#' 3. Reassigning data to groups with the same sample sizes as the original
#' 4. Computing the test statistic on the permuted data
#'
#' This approach tests the null hypothesis that group labels are exchangeable,
#' which is natural for testing whether sub-populations differ.
#' @export
one_permutation <- function(x, stat_fun) {
  # Extract all connectomes (raw data) from all groups
  # Connectomes are the primary data format and don't require centering
  all_data <- x$list_of_samples |>
    purrr::map(\(s) s$connectomes) |>
    purrr::reduce(c)

  # Get original sample sizes
  sample_sizes <- x$list_of_samples |>
    purrr::map_int(\(s) s$sample_size)

  # Randomly permute all data
  shuffled_data <- sample(all_data)

  # Split shuffled data back into groups (preserving sample sizes)
  start_idx <- 1
  permuted_samples <- sample_sizes |>
    purrr::map(\(n) {
      end_idx <- start_idx + n - 1
      group_data <- shuffled_data[start_idx:end_idx]
      start_idx <<- end_idx + 1

      # Create new CSample with permuted connectomes
      riemtan::CSample$new(conns = group_data, metric_obj = x$riem_metric)
    })

  # Create new CSuperSample and compute statistic
  permuted_samples |>
    riemtan::CSuperSample$new() |>
    stat_fun()
}

# Orientation of the statistics shipped with the package: which tail of the
# permutation distribution is evidence against the null. log Wilks' Lambda =
# log|W| - log|T| is large-negative when the groups differ; Pillai's trace =
# tr((T - W) T^-1) is large-positive. riem_anova() reads these.
attr(log_wilks_lambda, "riemstats_tail") <- "left"
attr(pillais_trace, "riemstats_tail") <- "right"

#' Compute p-values using permutation test
#'
#' Computes a permutation-based p-value for a given super sample.
#' The statistic used for the permutation test can be specified
#' via the `stat_fun` argument.
#'
#' @param ss An object of class `CSuperSample`.
#' @param stat_fun A function to compute a statistic
#' on the `CSuperSample` object (default: `log_wilks_lambda`).
#' @param nperm The number of permutations to generate
#' for estimating the p-value (default: 1000).
#' @param p_method Character; which permutation p-value estimator to report.
#' `"plugin"` (the default) returns the proportion of permuted statistics
#' strictly more extreme than the observed one. `"exact"` returns
#' `(1 + k) / (nperm + 1)`, where `k` counts the permuted statistics *at least
#' as* extreme, so that the observed labelling is itself in the reference set.
#' Both are in common use and neither dominates the other; see the
#' *Choice of estimator* section below.
#' @param tail Character or `NULL`; which tail of the permutation distribution
#' counts as evidence against the null. `"left"` for statistics that are small
#' under the alternative (`log_wilks_lambda`), `"right"` for statistics that
#' are large under the alternative (`pillais_trace`). The default, `NULL`,
#' reads the `riemstats_tail` attribute that the statistics shipped with this
#' package carry, and falls back to `"left"` for a user-supplied statistic
#' that declares no orientation.
#'
#' @return numeric A permutation-based p-value in `[0, 1]`.
#'
#' @details
#' The function computes the statistic on the observed data and compares it to
#' the distribution of statistics computed on permuted samples, in whichever
#' tail `tail` selects.
#'
#' The permutation test:
#' 1. Computes the test statistic on the observed data
#' 2. Randomly shuffles group assignments while preserving sample sizes
#' 3. Recomputes the test statistic on each permuted dataset
#' 4. Locates the observed statistic in the resulting permutation
#'    distribution and converts its position into a p-value
#'
#' This approach is computationally efficient and does not require
#' parameter estimation or synthetic data generation.
#'
#' # Choice of estimator
#'
#' Write `B` for `nperm`. Let `k` be the number of permuted statistics
#' *strictly* more extreme than the observed one and `k2` the number *at least
#' as* extreme; the two differ only when a permuted statistic ties the observed
#' one. Two estimators of the permutation p-value are in common use, and this
#' function exposes both.
#'
#' `p_method = "plugin"` reports `k / B`. It is an *unbiased* estimator of the
#' underlying permutation p-value, but it is not level-preserving: it can
#' return exactly 0, which is not a valid p-value at any nominal level. Its
#' attainable size at nominal `alpha` is `(floor(alpha * B) + 1) / (B + 1)`,
#' which exceeds `alpha` unless `alpha * (B + 1)` happens to be an integer.
#' At `B = 49`, a nominal 0.05 test has size 0.060; at `B = 999` the
#' discrepancy is negligible.
#'
#' `p_method = "exact"` reports `(1 + k2) / (B + 1)`, treating the observed
#' labelling as one of the `B + 1` equally likely labellings under the null.
#' It is *valid* -- its size never exceeds the nominal level -- at the cost of a
#' small upward bias, and it can never return 0.
#'
#' `"plugin"` is the default because it is the construction this function has
#' always used: results published from earlier versions remain reproducible by
#' the default call. `"exact"` is recommended for new work, and the default is
#' expected to change in a future major version.
#' @export
riem_anova <- function(ss, stat_fun = log_wilks_lambda, nperm = 1000,
                       p_method = c("plugin", "exact"), tail = NULL) {
  if (!inherits(ss, "CSuperSample")) {
    stop("Argument 'ss' must be an object of class 'CSuperSample'.")
  }

  # Check for minimum number of groups
  if (length(ss$list_of_samples) < 2) {
    stop("CSuperSample must contain at least 2 groups for ANOVA analysis")
  }

  p_method <- match.arg(p_method)

  # Which tail counts as evidence against the null. The statistics shipped
  # with the package declare their own orientation; a user-supplied statistic
  # that declares none falls back to "left", the orientation this function
  # assumed unconditionally before the argument existed.
  if (is.null(tail)) {
    tail <- attr(stat_fun, "riemstats_tail")
    if (is.null(tail)) {
      tail <- "left"
    }
  }
  tail <- match.arg(tail, c("left", "right"))

  # Compute observed test statistic
  stat_val <- stat_fun(ss)

  # Generate permutation distribution
  perm_vals <- 1:nperm |>
    purrr::map_dbl(
      \(m) one_permutation(ss, stat_fun)
    )

  if (identical(p_method, "plugin")) {
    # #{strictly more extreme} / nperm. Preserved exactly as written before
    # p_method existed, so that the default call reproduces published results.
    if (identical(tail, "left")) {
      mean(stat_val > perm_vals)
    } else {
      mean(stat_val < perm_vals)
    }
  } else {
    # (1 + #{at least as extreme}) / (nperm + 1).
    n_extreme <- if (identical(tail, "left")) {
      sum(perm_vals <= stat_val)
    } else {
      sum(perm_vals >= stat_val)
    }
    (1 + n_extreme) / (nperm + 1)
  }
}
