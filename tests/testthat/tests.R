test_pd_mats <- list(
  Matrix::Matrix(c(2.0, 0.5, 0.5, 3.0), nrow = 2) |>
    Matrix::nearPD() |> _$mat |> Matrix::pack(),
  Matrix::Matrix(c(1.5, 0.3, 0.3, 2.5), nrow = 2) |>
    Matrix::nearPD() |> _$mat |> Matrix::pack()
)

test_that("OK", {
  TRUE |> expect_true()
})

test_that("ANOVA stats work", {
  library(riemtan)
  data("airm")
  sam1 <- test_pd_mats |>
    purrr::map(\(x) 2 * x) |>
    CSample$new(metric_obj = airm)
  sam2 <- test_pd_mats |> CSample$new(metric_obj = airm)
  ss <- list(sam1, sam2) |> CSuperSample$new()

  # ss$Log_Wilks_Lambda() |> (\(x) {
  ss |>
    log_wilks_lambda() |>
    (\(x){
      list(
        x |> is.null() |> expect_false(),
        x |> inherits("numeric") |> expect_true(),
        x |> expect_lt(0)
      )
    })()
  ss |>
    pillais_trace() |>
    (\(x) {
      list(
        x |> is.null() |> expect_false(),
        x |> inherits("numeric") |> expect_true(),
        x |> expect_gt(0)
      )
    })()
  ss |>
    frechet_anova() |>
    (\(x) {
      list(
        x |> is.null() |> expect_false(),
        x |> inherits("numeric") |> expect_true(),
        x |> (\(x) x >= 0)() |> expect_true(),
        x |> (\(x) x <= 1)() |> expect_true()
      )
    })()
  ss |>
    riem_anova() |>
    (\(x) {
      list(
        x |> is.null() |> expect_false(),
        x |> inherits("numeric") |> expect_true(),
        x |> (\(x) x >= 0)() |> expect_true(),
        x |> (\(x) x <= 1)() |> expect_true()
      )
    })()
})

test_that("log_wilks_lambda input validation", {
  expect_error(
    log_wilks_lambda("not_a_supersample"),
    "Argument 'super_sample' must be an object of class 'CSuperSample'"
  )
  expect_error(
    log_wilks_lambda(123),
    "Argument 'super_sample' must be an object of class 'CSuperSample'"
  )
  expect_error(
    log_wilks_lambda(list()),
    "Argument 'super_sample' must be an object of class 'CSuperSample'"
  )
  expect_error(
    log_wilks_lambda(matrix(1:4, nrow = 2)),
    "Argument 'super_sample' must be an object of class 'CSuperSample'"
  )
  expect_error(
    log_wilks_lambda(NULL),
    "Argument 'super_sample' must be an object of class 'CSuperSample'"
  )
})

test_that("pillais_trace input validation", {
  expect_error(
    pillais_trace("not_a_supersample"),
    "Argument 'super_sample' must be an object of class 'CSuperSample'"
  )
  expect_error(
    pillais_trace(123),
    "Argument 'super_sample' must be an object of class 'CSuperSample'"
  )
  expect_error(
    pillais_trace(list()),
    "Argument 'super_sample' must be an object of class 'CSuperSample'"
  )
  expect_error(
    pillais_trace(matrix(1:4, nrow = 2)),
    "Argument 'super_sample' must be an object of class 'CSuperSample'"
  )
  expect_error(
    pillais_trace(NULL),
    "Argument 'super_sample' must be an object of class 'CSuperSample'"
  )
})

test_that("frechet_anova input validation", {
  library(riemtan)
  
  # Test with invalid data types
  expect_error(frechet_anova("not_a_supersample"))
  expect_error(frechet_anova(123))
  expect_error(frechet_anova(list()))
  expect_error(frechet_anova(matrix(1:4, nrow = 2)))
  expect_error(frechet_anova(NULL))
  
  # Test with empty CSuperSample (should still error due to invalid structure)
  expect_error({
    empty_ss <- try(CSuperSample$new(list()), silent = TRUE)
    if (!inherits(empty_ss, "try-error")) {
      frechet_anova(empty_ss)
    }
  })
})

test_that("riem_anova input validation", {
  library(riemtan)
  
  # Test with invalid data types for super_sample
  expect_error(riem_anova("not_a_supersample"))
  expect_error(riem_anova(123))
  expect_error(riem_anova(list()))
  expect_error(riem_anova(matrix(1:4, nrow = 2)))
  expect_error(riem_anova(NULL))
})

test_that("empty and malformed CSuperSample objects", {
  library(riemtan)
  data("airm")
  
  # Test with CSuperSample containing empty CSample objects
  expect_error({
    empty_sample <- try({
      CSample$new(list(), metric_obj = airm)
    }, silent = TRUE)
    if (!inherits(empty_sample, "try-error")) {
      empty_ss <- CSuperSample$new(list(empty_sample))
      log_wilks_lambda(empty_ss)
    }
  })
  
  expect_error({
    empty_sample <- try({
      CSample$new(list(), metric_obj = airm)
    }, silent = TRUE)
    if (!inherits(empty_sample, "try-error")) {
      empty_ss <- CSuperSample$new(list(empty_sample))
      pillais_trace(empty_ss)
    }
  })
  
  expect_error({
    empty_sample <- try({
      CSample$new(list(), metric_obj = airm)
    }, silent = TRUE)
    if (!inherits(empty_sample, "try-error")) {
      empty_ss <- CSuperSample$new(list(empty_sample))
      frechet_anova(empty_ss)
    }
  })
  
  expect_error({
    empty_sample <- try({
      CSample$new(list(), metric_obj = airm)
    }, silent = TRUE)
    if (!inherits(empty_sample, "try-error")) {
      empty_ss <- CSuperSample$new(list(empty_sample))
      riem_anova(empty_ss)
    }
  })
  
  # Test with single sample (should error due to insufficient groups)
  single_sample <- test_pd_mats |> CSample$new(metric_obj = airm)
  single_ss <- CSuperSample$new(list(single_sample))
  
  expect_error(log_wilks_lambda(single_ss))
  expect_error(pillais_trace(single_ss))
  expect_error(frechet_anova(single_ss))
  expect_error(riem_anova(single_ss))
  
  # Test harmonization functions with malformed CSuperSample
  expect_error(combat_harmonization(single_ss))
  expect_error(rigid_harmonization(single_ss))
})

test_that("harmonization functions input validation", {
  library(riemtan)
  data("airm")
  
  # Test combat_harmonization with invalid inputs
  expect_error(
    combat_harmonization("not_a_supersample"),
    "subscript out of bounds|object of type 'character' is not subsettable"
  )
  expect_error(
    combat_harmonization(123),
    "subscript out of bounds|object of type 'double' is not subsettable"
  )
  expect_error(
    combat_harmonization(list()),
    "subscript out of bounds"
  )
  expect_error(
    combat_harmonization(matrix(1:4, nrow = 2)),
    "subscript out of bounds|object of type 'integer' is not subsettable"
  )
  expect_error(
    combat_harmonization(NULL),
    "subscript out of bounds|object of type 'NULL' is not subsettable"
  )
  
  # Test rigid_harmonization with invalid inputs
  expect_error(
    rigid_harmonization("not_a_supersample"),
    "subscript out of bounds|object of type 'character' is not subsettable"
  )
  expect_error(
    rigid_harmonization(123),
    "subscript out of bounds|object of type 'double' is not subsettable"
  )
  expect_error(
    rigid_harmonization(list()),
    "subscript out of bounds"
  )
  expect_error(
    rigid_harmonization(matrix(1:4, nrow = 2)),
    "subscript out of bounds|object of type 'integer' is not subsettable"
  )
  expect_error(
    rigid_harmonization(NULL),
    "subscript out of bounds|object of type 'NULL' is not subsettable"
  )
  
  # Test with empty CSuperSample (should error due to no samples)
  expect_error({
    empty_ss <- try(CSuperSample$new(list()), silent = TRUE)
    if (!inherits(empty_ss, "try-error")) {
      combat_harmonization(empty_ss)
    }
  })
  
  expect_error({
    empty_ss <- try(CSuperSample$new(list()), silent = TRUE)
    if (!inherits(empty_ss, "try-error")) {
      rigid_harmonization(empty_ss)
    }
  })
}

test_that("normalization input validation", {
  # Test with invalid matrix inputs
  expect_error(
    normalization("not_a_matrix"),
    "non-numeric argument|invalid 'type'|argument is not a matrix"
  )
  expect_error(
    normalization(123),
    "non-numeric argument|invalid 'type'|argument is not a matrix"
  )
  expect_error(
    normalization(NULL),
    "object 'si' not found|subscript out of bounds"
  )
  expect_error(
    normalization(list(1, 2, 3)),
    "non-numeric argument|invalid 'type'"
  )
  
  # Test with edge cases
  expect_error(normalization(matrix(numeric(0))))  # Empty matrix
  expect_error(normalization(matrix(NA, nrow = 2, ncol = 2)))  # Matrix with NAs
  expect_error(normalization(matrix(Inf, nrow = 2, ncol = 2)))  # Matrix with Inf
  
  # Test with valid matrices (should not error)
  expect_no_error(normalization(matrix(1:6, nrow = 2)))
  expect_no_error(normalization(matrix(rnorm(12), nrow = 3)))
  expect_no_error(normalization(matrix(1, nrow = 1, ncol = 5)))  # Single row
})

test_that("format_matr input validation", {
  # Test with invalid inputs
  expect_error(
    format_matr("not_a_matrix"),
    "cannot coerce|invalid 'type'"
  )
  expect_error(
    format_matr(123),
    "cannot coerce|invalid 'type'"
  )
  expect_error(
    format_matr(NULL),
    "cannot coerce|invalid 'type'"
  )
  expect_error(
    format_matr(list(1, 2, 3)),
    "cannot coerce|invalid 'type'"
  )
  
  # Test with edge cases
  expect_error(format_matr(matrix(numeric(0))))  # Empty matrix
  expect_error(format_matr(matrix(NA, nrow = 2, ncol = 2)))  # Matrix with NAs
  expect_error(format_matr(matrix(c(1, 2, 3, 4), nrow = 2)))  # Non-positive definite
  
  # Test with valid positive definite matrix
  valid_matrix <- matrix(c(2, 1, 1, 2), nrow = 2)
  expect_no_error(format_matr(valid_matrix))
})

test_that("ts2corr input validation", {
  # Test with invalid inputs
  expect_error(
    ts2corr("not_a_matrix"),
    "non-numeric argument|invalid 'type'|argument is not a matrix"
  )
  expect_error(
    ts2corr(123),
    "non-numeric argument|invalid 'type'|argument is not a matrix"
  )
  expect_error(
    ts2corr(NULL),
    "object 'ts' not found|subscript out of bounds"
  )
  expect_error(
    ts2corr(list(1, 2, 3)),
    "non-numeric argument|invalid 'type'"
  )
  
  # Test with edge cases
  expect_error(ts2corr(matrix(numeric(0))))  # Empty matrix
  expect_error(ts2corr(matrix(NA, nrow = 2, ncol = 2)))  # Matrix with NAs
  expect_error(ts2corr(matrix(Inf, nrow = 2, ncol = 2)))  # Matrix with Inf
  
  # Test with valid time series matrix
  valid_ts <- matrix(rnorm(20), nrow = 4, ncol = 5)
  expect_no_error(ts2corr(valid_ts))
  
  # Test with minimal valid matrix
  minimal_ts <- matrix(c(1, 2, 3, 4), nrow = 2, ncol = 2)
  expect_no_error(ts2corr(minimal_ts))
})
