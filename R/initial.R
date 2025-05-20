Log_Wilks_Lambda = function() {
      if (self$Within |> is.null()) self$compute_W()

      if (self$Total |> is.null()) self$compute_T()

      result <- Matrix::determinant(self$Within)$modulus -
        Matrix::determinant(self$Total)$modulus

      result
    }
    Pillais_Trace = function() {
      if (self$Within |> is.null()) self$compute_W()

      if (self$Total |> is.null()) self$compute_T()

      result <- ((self$Total - self$Within) %*% solve(self$Total)) |>
        as.matrix() |>
        diag() |>
        sum()

      result
    }