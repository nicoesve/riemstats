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
