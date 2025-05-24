format_matr <- function(x) {
  x |>
    as("dpoMatrix") |>
    Matrix::pack()
}

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

harmonize_with_combat <- function(super_sample) {
  # applying ComBat
  harmonized_vector_images <- super_sample$list_of_samples |>
    purrr::imap(
      \(sample, idx) {
        sample$compute_tangents()
        sample$compute_vecs()
        data <- sample$vector_images
        batch <- rep(idx, data |> nrow())
        cbind(data, batch)
      }
    ) |>
    purrr::reduce(rbind) |>
    (\(m) list(m[, -ncol(m)], m[, ncol(m)]))() |>
    do.call(what = sva::ComBat, args = _)

  # Reconstructing
  batches <- harmonized_vector_images[, ncol(harmonized_vector_images)]
  vec_imgs <- harmonized_vector_images[
    , -ncol(harmonized_vector_images)
  ]

  harmonized_vector_images |>
    nrow() |>
    seq_len() |>
    split(batches) |>
    purrr::map(
      \(idx) {
        riemtan::CSample$new(
          vec_imgs = vec_imgs[idx, ],
          centered = FALSE,
          ref_pt = idx |> length() |> diag() |> format_matr(),
          metric_obj = super_sample$riem_metr
        )
      }
    ) |>
    riemtan::CSuperSample$new()
}
