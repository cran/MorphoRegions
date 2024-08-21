#' Simulate regions data
#'
#' `simregions()` simulates vertebrae and PCOs that satisfy certain constraints.
#'
#' @param nvert `numeric`; the number of vertebrae for which to simulate data.
#' @param nregions `numeric`; the desired number of regions in the simulated data.
#' @param nvar `numeric`; the number of PCO axes to simulate. Default is 1.
#' @param r2 `numeric`; a vector containing the \eqn{R^2} of the true segmented regression model for each simulated PCO. If a single value is supplied, all PCOs will receive that value. Otherwise, one value should be supplied for each simulated PCO.
#' @param minvert `numeric`; the minimum number of vertebrae allowed in each region. Default is 3.
#' @param cont `logical`; whether to use models that are continuous (`TRUE`) or discontinuous (`FALSE`) at the breakpoints to generate the data. Default is `TRUE`.
#' @param sl.dif `numeric`; the minimum required difference in slopes between adjacent regions, expressed as a proportion of the maximal difference between allowable slopes. Must be between 0 and 1. See Details.
#' @param x a `regions_sim` object.
#' @param scores `numeric`; for which simulated PCO scores the simulated values should be plotted.
#' @param lines `logical`; whether to display the simulated regression lines on the plot. Default is `TRUE`.
#' @param \dots ignored.
#'
#' @returns `simregions()` returns a `regions_sim` object, which contains the vertebra indices in the `Xvar` entry, the PCO scores in the `Yvar` entry, the simulated breakpoints in the `BPs` entry, the simulated model coefficients in the `coefs` entry, and the simulated error standard deviation in the `ersd` entry. The attribute `"design"` contains the design matrix, which when multiplied by the coefficients and added to a random normal variate with standard deviation equal to the error standard deviation yields the observed PCO scores.
#'
#' `plot()` returns a `ggplot` object that can be manipulated using `ggplot2` syntax. The plot is similar to that produced by [plot.regions_pco()] and to that produced by [plotsegreg()] except that the displayed lines (if requested) are the true rather than fitted regression lines.
#'
#' @details
#' `simregions()` generates PCO scores for each requested vertebra such that certain conditions are met. The slopes for each region are drawn from a uniform distribution with limits of -.5 and .5. If a set of slopes contains two adjacent slopes that have a difference less than `sl.dif`, it is rejected and a new one is drawn. The scaling of the PCOs is determined by the slopes and the requested \eqn{R^2}. The PCOs will not necessarily be in order from most variable to least variable as they are in a traditional PCO analysis.
#'
#' Intercepts (the intercept of the first region when `cont = TRUE` and the intercept of all regions when `cont = FALSE`) are drawn from a uniform distribution with limits of \eqn{-n/4} and \eqn{n/4}, where \eqn{n} is the number of breakpoints, one less than `nregions`. Intercepts other than the first when `cont = TRUE` are determined by the slopes.
#'
#' The `cont`, `r2`, and `sl.dif` arguments control how easy it is for fitted segmented regression models to capture the true structure. When `cont = TRUE`, it can be harder to determine exactly where regions begin and end, especially if `sl.dif` is 0. When `r2` is high, there is little variation around the true line, so the fitted lines will be more precise and region boundaries clearer. When `sl.dif` is high, slopes of adjacent regions are different from each other, so it is easier to detect region boundaries. Setting `sl.dif` to between .5 and 1 ensures that the slopes in adjacent regions have different signs.
#'
#' @seealso
#' [calcregions()] for fitting segmented regression models to the simulated data; [calcmodel()] for fitting a single segmented regression model to the simulated data; [plotsegreg()] for plotting estimated regression lines.
#'
#' @example man/examples/example-simregions.R
#'

#' @export
simregions <- function(nvert, nregions, nvar = 1, r2 = .95,
                       minvert = 3, cont = TRUE, sl.dif = 0) {

  #Checks
  chk::chk_not_missing(nvert, "`nvert`")
  chk::chk_whole_number(nvert)
  chk::chk_gte(nvert, 4)

  chk::chk_not_missing(nregions, "`nregions`")
  chk::chk_count(nregions)

  chk::chk_count(nvar)
  chk::chk_gte(nvar, 1)

  chk::chk_numeric(r2)
  chk::chk_gt(r2, 0)
  chk::chk_lte(r2, 1)

  chk::chk_whole_number(minvert)
  chk::chk_gte(minvert, 2)

  chk::chk_flag(cont)

  chk::chk_number(sl.dif)
  chk::chk_gte(sl.dif, 0)
  chk::chk_lt(sl.dif, 1)

  chk::chk_lte(minvert * nregions, nvert)

  Xvar <- seq_len(nvert)

  nbp <- nregions - 1

  if (length(r2) == 1L) {
    r2 <- rep(r2, nvar)
  }

  chk::chk_equivalent(length(r2), nvar)

  #Draw possible combination of breakpoints
  bp <- .sample_combosR(Xvar, minvert, nbp)

  #Create design matrix
  x <- .design_matrix(Xvar, bp, cont)

  #Initialize matrix of coefficients
  coefs <- matrix(nrow = ncol(x), ncol = nvar)

  rownames(coefs) <- {
    if (cont) c("int1", paste0("slope", seq_len(nregions)))
    else paste0(rep(c("int", "slope"), each = nregions),
                rep(seq_len(nregions), 2))

  }
  colnames(coefs) <- paste0("PC.", seq_len(nvar))

  if (cont) {
    for (i in seq_len(nvar)) {
      # Rejection sampling for slopes
      repeat {
        if (sl.dif > .5) {
          test_slopes <- rep(0, nregions)
          test_slopes[1] <- (-1)^rbinom(1, 1, .5) * runif(1, min = sl.dif - .5, max = .5)

          for (j in seq_along(test_slopes)[-1]) {
            test_slopes[j] <- -sign(test_slopes[j - 1]) * runif(1, min = sl.dif - .5, max = .5)
          }
        }
        else {
          test_slopes <- runif(nregions, min = -.5, max = .5)
        }

        if (all(abs(diff(test_slopes)) >= sl.dif)) break
      }

      for (j in seq_along(test_slopes)[-1]) {
        test_slopes[j] <- test_slopes[j] - test_slopes[j - 1]
      }

      coefs[startsWith(rownames(coefs), "slope"), i] <- test_slopes

      coefs["int1", i] <- runif(1, -nbp/4, nbp/4)
    }
  }
  else {
    for (i in seq_len(nvar)) {
      # Rejection sampling for slopes
      repeat {
        if (sl.dif > .5) {
          test_slopes <- rep(0, nregions)
          test_slopes[1] <- (-1)^rbinom(1, 1, .5) * runif(1, min = sl.dif - .5, max = .5)

          for (j in seq_along(test_slopes)[-1]) {
            test_slopes[j] <- -sign(test_slopes[j - 1]) * runif(1, min = sl.dif - .5, max = .5)
          }
        }
        else {
          test_slopes <- runif(nregions, min = -.5, max = .5)
        }

        if (all(abs(diff(test_slopes)) >= sl.dif)) break
      }

      coefs[startsWith(rownames(coefs), "slope"), i] <- test_slopes

      ints <- runif(nregions, -nbp/4, nbp/4)
      coefs[startsWith(rownames(coefs), "int"), i] <- ints - test_slopes * c(0, bp)
    }
  }

  #Generate structural outcomes from coefs and design
  ypred <- x %*% coefs

  #Generate actual outcomes by adding error
  ersd <- numeric(nvar)
  y <- array(dim = dim(ypred), dimnames = dimnames(ypred))
  for (i in seq_len(nvar)) {
    #SD(error) required to make requested R2
    ersd[i] <- sqrt((1 - r2[i]) * var(ypred[,i]) / r2[i])

    y[,i] <- ypred[,i] + rnorm(nvert, 0, ersd[i])
  }

  out <- list(Xvar = Xvar, Yvar = y, BPs = bp, coefs = coefs, ersd = ersd)
  attr(out, "cont") <- cont
  attr(out, "design") <- x

  class(out) <- "regions_sim"
  out
}

#' @exportS3Method plot regions_sim
#' @rdname simregions
plot.regions_sim <- function(x, scores = 1, lines = TRUE, ...) {

  chk::chk_whole_numeric(scores)
  chk::chk_range(scores, c(1, ncol(x$Yvar)))
  chk::chk_flag(lines)

  Yvar <- x$Yvar[,scores, drop = FALSE]
  Xvar <- x$Xvar

  yhat <- attr(x, "design") %*% x$coef[,scores, drop = FALSE]

  .plotreg_internal(Xvar, Yvar, yhat, x$BPs, lines, scores = scores,
                    linescolor = "gray30")
}

#' @exportS3Method print regions_sim
print.regions_sim <- function(x, ...) {
  cat("A `regions_sim` object\n")
  cat(" - number of vertebrae:", length(x$Xvar), "\n")
  cat(" - number of regions:", length(x$BPs) + 1, "\n")
  cat(" - breakpoints:", paste(x$BPs, collapse = ", "), "\n")
  cat(" - model type:", if (attr(x, "cont")) "continuous" else "discontinuous", "\n")
  cat(" - number of PCO scores:", ncol(x$Yvar), "\n")
  cat("Use `plot()` to display the true lines and simulated data.")
}

# Sample combos without needing to store all combos
.sample_combosR <- function(vert, minvert, nbp) {

  n <- length(vert)

  #Choose combination at random
  maxcombos <- .ncombos(n, minvert, nbp)

  k <- sample.int(maxcombos, 1L)

  .kth_combo(k, vert, minvert, nbp)
}

# Find the kth combo; assumes k <= maxcombs
.kth_combo <- function(k, vert, minvert, nbp) {

  n <- length(vert)

  mins <- minvert * seq_len(nbp)
  maxes <- n - rev(mins)

  combo <- rep.int(0, nbp)

  a <- 0

  for (i in seq_len(nbp)) {
    if (i < nbp) {
      sums <- cumsum(.ncombos(n - mins[i]:maxes[i], minvert, nbp - i))
      b <- min(which(k <= a + sums))
      combo[i] <- mins[i] + b - 1
      a <- a + c(0, sums)[b]
      mins[i + 1] <- combo[i] + minvert
    }
    else {
      sums <- seq_len(maxes[i] - mins[i] + 1)
      b <- min(which(k <= a + sums))
      combo[i] <- mins[i] + b - 1
    }
  }

  combo[] <- vert[combo]

  setNames(combo, paste0("breakpoint", seq_along(combo)))
}