#' Format a Matrix as a Packed dpoMatrix
#'
#' Converts a matrix to a packed symmetric positive definite matrix
#' (`dpoMatrix`) using the Matrix package.
#'
#' @param x A numeric matrix.
#' @return A packed `dpoMatrix` object.
#' @importFrom Matrix pack
#' @export
format_matr <- function(x) {
  x |>
    as("dpoMatrix") |>
    Matrix::pack()
}

#' Normalize Rows of a Matrix
#'
#' Centers and scales each row of a matrix to have zero mean and unit norm.
#'
#' @param si A numeric matrix.
#' @return A matrix with each row centered and scaled.
#' @export
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

#' Compute OAS-Shrunk Correlation Matrix from Time Series
#'
#' Normalizes a time series matrix and computes
#' its OAS-shrunk correlation matrix.
#'
#' @param ts A numeric matrix representing time series data.
#' @return A shrunk correlation matrix.
#' @importFrom CovTools CovEst.2010OAS
#' @export
ts2corr <- function(ts) {
  ts <- normalization(ts)
  c <- ts %*% t(ts)
  c |>
    CovTools::CovEst.2010OAS() |>
    _$S
}

#' Harmonize Vector Images Using ComBat
#'
#' Applies the ComBat batch effect correction to vector images
#' in a `CSuperSample` object and reconstructs harmonized samples.
#'
#' @param super_sample A `CSuperSample` object containing samples to harmonize.
#' @return A new `CSuperSample` object with harmonized vector images.
#' @importFrom sva ComBat
#' @export
combat_harmonization <- function(super_sample) {
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
          ref_pt = super_sample$matrix_size |> diag() |> format_matr(),
          metric_obj = super_sample$riem_metr
        )
      }
    ) |>
    riemtan::CSuperSample$new()
}


#' Harmonize Tangent Images Across Batches Using Rigid Correction
#'
#' Applies a rigid harmonization procedure to tangent images in a `CSuperSample` object.
#' First, batch means are subtracted from each sample's tangent images (batch correction),
#' then the overall mean is added back (global correction). The harmonized tangent images
#' are used to reconstruct new `CSample` objects, which are collected into a new `CSuperSample`.
#'
#' @param super_sample A `CSuperSample` object containing samples to harmonize.
#' @return A new `CSuperSample` object with harmonized tangent images.
#' @export
#'
#' @details
#' This function performs harmonization in two steps:
#' \enumerate{
#'   \item \strong{Batch Correction:} For each batch, the mean of its tangent images is subtracted from each tangent image in the batch.
#'   \item \strong{Global Correction:} The overall mean (across all batches) of tangent images is added back to each tangent image.
#' }
#' The harmonized tangent images are then used to reconstruct samples using the reference point and metric from the original `CSuperSample`.
#'
#' @examples
#' # Assuming `super_sample` is a CSuperSample object:
#' # harmonized <- rigid_harmonization(super_sample)
rigid_harmonization <- function(super_sample) {
  tangent_collection <- super_sample$list_of_samples |>
    purrr::map(
      \(sample) {
        sample$compute_tangents()
        sample$tangent_images
      }
    )

  batch_means <- tangent_collection |>
    purrr::map(
      \(l) {
        l |>
          purrr::reduce("+") |>
          (\(x) x / length(l))()
      }
    )

  overall_mean <- tangent_collection |>
    do.call(what = c, args = _) |>
    purrr::reduce("+") |>
    (\(x) x / super_sample$sample_size)()

  list(tangent_collection, batch_means) |>
    # Batch correction
    purrr::pmap(
      \(list_of_tangents, batch_mean){
        list_of_tangents |>
          purrr::map(
            \(tgt) tgt - batch_mean
          )
      }
    ) |>
    purrr::map(
      \(list_of_tangents) {
        aux_id <- super_sample$matrix_size |>
          diag() |>
          format_matr()

        # Global correction
        list_of_tangents |>
          purrr::map(
            \(tgt) tgt + overall_mean
          ) |>
          riemtan::CSample$new(
            tan_imgs = _,
            centered = FALSE,
            ref_pt = aux_id,
            metric_obj = super_sample$riem_metr
          )
      }
    ) |>
    riemtan::CSuperSample$new()
}
