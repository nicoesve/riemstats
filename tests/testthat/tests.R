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

test_that("harmonization functions input validation", {
  library(riemtan)
  
  # Test combat_harmonization with invalid inputs
  expect_error(combat_harmonization("not_a_supersample"))
  expect_error(combat_harmonization(123))
  expect_error(combat_harmonization(list()))
  expect_error(combat_harmonization(matrix(1:4, nrow = 2)))
  expect_error(combat_harmonization(NULL))
  
  # Test rigid_harmonization with invalid inputs
  expect_error(rigid_harmonization("not_a_supersample"))
  expect_error(rigid_harmonization(123))
  expect_error(rigid_harmonization(list()))
  expect_error(rigid_harmonization(matrix(1:4, nrow = 2)))
  expect_error(rigid_harmonization(NULL))
})

test_that("normalization input validation", {
  # Test with invalid matrix inputs
  expect_error(normalization("not_a_matrix"))
  expect_error(normalization(123))
  expect_error(normalization(NULL))
  
  # Test with valid matrices (should not error)
  expect_no_error(normalization(matrix(1:6, nrow = 2)))
  expect_no_error(normalization(matrix(rnorm(12), nrow = 3)))
})

test_that("format_matr input validation", {
  # Test with invalid inputs
  expect_error(format_matr("not_a_matrix"))
  expect_error(format_matr(123))
  expect_error(format_matr(NULL))
  
  # Test with valid positive definite matrix
  valid_matrix <- matrix(c(2, 1, 1, 2), nrow = 2)
  expect_no_error(format_matr(valid_matrix))
})

test_that("ts2corr input validation", {
  # Test with invalid inputs
  expect_error(ts2corr("not_a_matrix"))
  expect_error(ts2corr(123))
  expect_error(ts2corr(NULL))
  
  # Test with valid time series matrix
  valid_ts <- matrix(rnorm(20), nrow = 4, ncol = 5)
  expect_no_error(ts2corr(valid_ts))
})
