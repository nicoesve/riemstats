library(riemtan)
library(riemstats)
library(Matrix)
library(purrr)

test_pd_mats <- list(
  Matrix::Matrix(c(2.0, 0.5, 0.5, 3.0), nrow = 2) |>
    Matrix::nearPD() |> _$mat |> Matrix::pack(),
  Matrix::Matrix(c(1.5, 0.3, 0.3, 2.5), nrow = 2) |>
    Matrix::nearPD() |> _$mat |> Matrix::pack()
)

data("airm")

sam1 <- test_pd_mats |>
  purrr::map(\(x) (2 * x) |> Matrix::unpack() |> as("dpoMatrix") |> Matrix::pack()) |>
  CSample$new(metric_obj = airm)

sam2 <- test_pd_mats |> CSample$new(metric_obj = airm)

ss <- list(sam1, sam2) |> CSuperSample$new()

cat("Testing log_wilks_lambda...\n")
result1 <- try(log_wilks_lambda(ss), silent = FALSE)
if (inherits(result1, "try-error")) {
  cat("ERROR in log_wilks_lambda\n")
  print(result1)
} else {
  cat("log_wilks_lambda succeeded:", result1, "\n")
}

cat("\nTesting pillais_trace...\n")
result2 <- try(pillais_trace(ss), silent = FALSE)
if (inherits(result2, "try-error")) {
  cat("ERROR in pillais_trace\n")
  print(result2)
} else {
  cat("pillais_trace succeeded:", result2, "\n")
}

cat("\nTesting frechet_anova...\n")
result3 <- try(frechet_anova(ss), silent = FALSE)
if (inherits(result3, "try-error")) {
  cat("ERROR in frechet_anova\n")
  print(result3)
} else {
  cat("frechet_anova succeeded:", result3, "\n")
}

cat("\nTesting riem_anova...\n")
result4 <- try(riem_anova(ss), silent = FALSE)
if (inherits(result4, "try-error")) {
  cat("ERROR in riem_anova\n")
  print(result4)
} else {
  cat("riem_anova succeeded:", result4, "\n")
}

cat("\nAll tests completed.\n")