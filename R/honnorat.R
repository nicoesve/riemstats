# normalization: Center and scale each row of a matrix
normalization <- function(si) {
  on <- matrix(1, nrow = 1, ncol = ncol(si))
  eps <- 1e-9
  row_means <- rowMeans(si)
  si <- si -
    matrix(
      row_means,
      nrow = nrow(si),
      ncol = ncol(si), byrow = FALSE
    ) %*% on
  row_norms <- sqrt(rowSums(si * si))
  aux_matr <- matrix(row_norms, nrow = nrow(si), ncol = ncol(si), byrow = FALSE)
  si / pmax(aux_matr %*% on, eps * on)
}

# ts2corr: Normalize time series and compute OAS-shrunk correlation matrix
ts2corr <- function(ts) {
  ts <- normalization(ts)
  c <- ts %*% t(ts)
  c |>
    CovTools::CovEst.2010OAS() |>
    _$S
}
