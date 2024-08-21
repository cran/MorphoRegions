#' Calculate PCO (principal co-ordinates analysis) based on SVD
#'
#' Calculates distance matrix from raw data, then conducts a PCO ordination using a
#' single value decomposition (SVD). This differs from other PCO functions which use [stats::cmdscale()] and rely on a
#' spectral decomposition.
#'
#' @param x a `regions_data` object; the output of a call to [process_measurements()].
#' @param metric string; the distance matrix calculation metric. Allowable options include those support by [cluster::daisy()], which are `"euclidean"`, `"manhattan"`, or `"gower"`. Default is `"gower"`. Abbreviations allowed.
#' @param scale `logical`; whether to scale the variables prior to including them in the PCO estimation. Default is `TRUE`, which is especially advisable when using the bootstrap to select the number of PCOs to use in downstream analyses. Passed to the `stand` argument of `cluster::daisy()`. Ignored if `metric = "gower"`.
#'
#' @return A `regions_pco` object, which contains eigenvectors in the `scores` component and eigenvalues in the `eigen.val` component. The original dataset is stored in the `data` attribute.
#'
#' @seealso
#' [plot.regions_pco()] for plotting PCO axes
#'
#' [cluster::daisy()], which is used to compute the distance matrix used in the calculation; [stats::cmdscale()] for a spectral decomposition-based implementation
#'
#' @example man/examples/example-svdPCO.R
#'

#' @export
svdPCO <- function(x, metric = "gower", scale = TRUE) {

  chk::chk_is(x, "regions_data")
  chk::chk_string(metric)
  metric <- tolower(metric)
  metric <- .match_arg(metric, eval(formals(cluster::daisy)[["metric"]]))
  chk::chk_flag(scale)

  dat <- .get_data_without_pos(x)

  if (metric != "gower" && !all(vapply(dat, is.numeric, logical(1L)))) {
    .wrn_immediate(sprintf("`metric = \"%s\" cannot be used when non-numeric measurement variables are present. Setting `metric` to `\"gower\"`", metric))
  }

  #Set distance metric
  dist <- cluster::daisy(dat, metric = metric, stand = scale)

  out <- .svdPCO_internal(dist)

  attr(out, "data") <- x
  attr(out, "metric") <- metric
  attr(out, "scale") <- scale
  attr(out, "specimen") <- factor(rep(seq_along(x), unlist(lapply(x, nrow))),
                                  labels = paste("Specimen", seq_along(x)),
                                  levels = seq_along(x))

  class(out) <- "regions_pco"

  out
}

.svdPCO_internal <- function(dist, val.only = FALSE) {

  ### PCO analysis based on SVD not eigen as in cmdscale()
  d2 <- as.matrix(dist^2)

  # Double centering code
  n <- nrow(d2)
  k <- ncol(d2)
  r_m <- matrix(rowMeans(d2, na.rm = TRUE), n, k, byrow = TRUE)
  c_m <- matrix(colMeans(d2, na.rm = TRUE), n, k, byrow = FALSE)
  m_m <- mean(d2, na.rm = TRUE)
  step3 <- -0.5 * (d2 - r_m - c_m + m_m)

  if (val.only) {
    #Primarily for use in PCOcutoff(); faster than full SVD
    eigen.val <- sqrt(abs(eigen(tcrossprod(step3), symmetric = TRUE,
                                only.values = TRUE)$values))
    return(list(eigen.val = eigen.val))
  }

  # Single value decomp
  step4 <- svd(step3)
  eigen.val <- step4$d
  step5 <- step4$v

  # Scale by root eigenval
  for (i in seq_len(ncol(step5))) {
    step5[, i] <- step5[, i] * sqrt(step4$d[i])
  }

  eigen.vect <- step5[, -ncol(step5), drop = FALSE]  #remove the last PCO which has no variance

  colnames(eigen.vect) <- paste0("PCO.", seq_len(ncol(eigen.vect)))

  list(scores = eigen.vect,
       eigen.val = eigen.val)
}

#' @exportS3Method print regions_pco
print.regions_pco <- function(x, digits = 3, ...) {
  cat("- Scores:\n")
  print(as.data.frame(x$scores, row.names = rownames(attr(x, "data")))[1:min(nrow(x$scores), 6),],
        digits = digits, ...)
  if (nrow(x$scores) > 6) {
    cat("(First 6 of", nrow(x$scores), "rows displayed.)\n")
  }
  cat("\n- Eigenvalues:\n")
  print(x$eigen.val, digits = digits, ...)
}