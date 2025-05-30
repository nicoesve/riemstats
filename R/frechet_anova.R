fanova <- function(super_sample) {
  k <- super_sample$list_of_samples |> length()
  n <- super_sample$sample_size

  super_sample$list_of_samples |>
    purrr::map(
      \(sample) {
        sample$vector_images |>
          apply(X = _, MARGIN = 1, FUN = function(x) sum(x^2))
      }
    ) |>
    (\(x) {
      # TODO
    })()


  # get local information
  vec.nj <- rep(0, k)
  vec.Vj <- rep(0, k)
  vec.sig2j <- rep(0, k)
  for (j in 1:k) {
    distj <- as.vector(distvecs[[j]])
    nj <- length(distj)

    vec.nj[j] <- nj
    vec.Vj[j] <- sum(distj^2) / nj
    vec.sig2j[j] <- (sum(distj^4) / nj) - ((sum(distj^2) / nj)^2)
  }

  # get global information
  Vp <- sum(distall^2) / n
  lbdj <- vec.nj / n

  # compute statistics
  Fn <- Vp - sum(lbdj * vec.Vj)
  Un <- 0
  for (j in 1:(k - 1)) {
    for (l in (j + 1):k) {
      Un <- Un + ((lbdj[j] * lbdj[l]) / (vec.sig2j[j] * vec.sig2j[l])) * ((vec.Vj[j] - vec.Vj[l])^2)
    }
  }
  term1 <- (n * Un) / sum(vec.nj / vec.sig2j)
  term2 <- (n * (Fn^2)) / sum((vec.nj^2) * vec.sig2j)
  thestat <- term1 + term2

  # compute p-value
  pvalue <- stats::pchisq(thestat, df = (k - 1), lower.tail = FALSE)

  list(thestat, pvalue)
}
