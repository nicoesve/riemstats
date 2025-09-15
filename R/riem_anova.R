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

#' Bootstrap Statistic for a Super Sample
#'
#' For each subsample in a given super sample,
#' draws a normal subsample of the same size
#' using the provided parameters, collects them into a new super sample,
#' and computes a specified statistic on the resulting super sample.
#'
#' @param x An object of class \code{CSuperSample},
#' representing the original super sample.
#' @param hat_sigma The mean parameter for the normal distribution
#' used in resampling.
#' @param hat_gamma The covariance parameter for the normal distribution
#' used in resampling.
#' @param geom The geometry parameter to be passed to \code{riemtan::rspdnorm}.
#' @param stat_fun A function to compute a statistic
#' on the resulting \code{CSuperSample} object.
#'
#' @return The value returned by \code{stat_fun} when applied
#' to the bootstrapped super sample.
#'
#' @details
#' This function performs a parametric bootstrap by resampling
#' each subsample in \code{x} according to the specified parameters,
#' then aggregates the resampled data into a new
#' \code{CSuperSample} and computes the desired statistic.
#' @export
one_bootstrap <- function(x, hat_sigma, hat_gamma, geom, stat_fun) {
  x$list_of_samples |>
    purrr::map(
      \(s) riemtan::rspdnorm(s$sample_size, hat_sigma, hat_gamma, geom)
    ) |>
    riemtan::CSuperSample$new() |>
    stat_fun()
}

#' Compute p-values using bootstrap
#'
#' Computes a bootstrap-based p-value and for a given super sample.
#' The statistic used for the bootstrap can be specified
#' via the `stat_fun` argument.
#'
#' @param ss An object of class `CSuperSample`.
#' @param stat_fun A function to compute a statistic
#' on the `CSuperSample` object (default: `log_wilks_lambda`).
#' @param den The number of bootstrap samples to generate
#' for estimating the p-value.
#'
#' @return numeric A bootstrap p-value.
#'
#' @details
#' The function computes the statistic on the observed data
#' and compares it to the distribution
#' of statistics computed on bootstrapped samples.
#' @export
riem_anova <- function(ss, stat_fun = log_wilks_lambda, den = 5) {
  if (!inherits(ss, "CSuperSample")) {
    stop("Argument 'ss' must be an object of class 'CSuperSample'.")
  }

  # Check for minimum number of groups
  if (length(ss$list_of_samples) < 2) {
    stop("CSuperSample must contain at least 2 groups for ANOVA analysis")
  }

  ss$gather()
  ss$compute_fmean()

  # estimate the parameters (the dispersion is estimated by pooling)
  hat_sigma <- ss$frechet_mean
  # hat_gamma <- ss$list_of_samples |>
  #   purrr::walk(\(s) s$compute_sample_cov()) |>
  #   purrr::map(\(s) (s$sample_size - 1) * s$sample_cov) |>
  #   print() |>
  #   Reduce(`+`, x = _) |>
  #   print() |>
  #   (\(x) x / (ss$sample_size - 1))() |>
  #   methods::as("dpoMatrix") |>
  #   Matrix::pack()
  hat_gamma <- diag(ss$mfd_dim) |>
    methods::as("dpoMatrix") |>
    Matrix::pack()

  # reference statistic value
  stat_val <- stat_fun(ss)

  # bootstraping
  1:den |>
    purrr::map_dbl(
      \(m) one_bootstrap(ss, hat_sigma, hat_gamma, ss$riem_metric, stat_fun)
    ) |>
    (\(v) stat_val > v)() |>
    mean()
}
