testthat("OK", {
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

ss$Log_Wilks_Lambda() |> (\(x) {
    list(
      x |> is.null() |> expect_false(),
      x |> inherits("numeric") |> expect_true(),
      x |> expect_lt(0)
    )
  })()
  ss$Pillais_Trace() |> (\(x) {
    list(
      x |> is.null() |> expect_false(),
      x |> inherits("numeric") |> expect_true(),
      x |> expect_gt(0)
    )
  })()
})