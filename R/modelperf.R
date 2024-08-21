#' Assess model performance
#'
#' `modelperf()` computes model performance statistics in the form of \eqn{R^2} measures for a given combination of breakpoints.
#'
#' @param x a `regions_pco` object, the output of a call to [svdPCO()], or a `regions_results_single` object, the output of a call to [calcmodel()].
#' @param scores `numeric`; the indices of the PCO scores for which fit statistics are to be computed.
#' @param modelsupport a `regions_modelsupport` object, the output of a call to [modelsupport()]. When `x` is a `regions_pco` object, either `modelsupport`, `criterion`, and `model` must be supplied or `bps` and `cont` must be supplied. See Details.
#' @param criterion string; the criterion to use to select the best model for which breakpoints are to be displayed when `modelsupport` is specified. Ignored otherwise. Allowable options include `"aic"` to use the AICc and `"bic"` to use the BIC. Abbreviations allowed. Default is `"aic"`. When `x` is a `regions_pco` object, either `modelsupport`, `criterion`, and `model` must be supplied or `bps` and `cont` must be supplied. See Details.
#' @param model `numeric`; for which model among the best as determined by `criterion` should fit statistics be computed. 1 is the best model, 2 the second best, etc. Default is 1. When `x` is a `regions_pco` object, either `modelsupport`, `criterion`, and `model` must be supplied or `bps` and `cont` must be supplied. See Details.
#' @param bps `numeric`; a vector of breakpoints for which model fit should be computed. When `x` is a `regions_pco` object, either `modelsupport`, `criterion`, and `model` must be supplied or `bps` and `cont` must be supplied. See Details.
#' @param cont `logical`; whether to fit a model that is continuous (`TRUE`) or discontinuous (`FALSE`) at the breakpoints supplied to `bps`. Default is `TRUE`. When `x` is a `regions_pco` object, either `modelsupport`, `criterion`, and `model` must be supplied or `bps` and `cont` must be supplied. See Details.
#' @param \dots ignored.
#'
#' @return A `regions_perf` object containing the breakpoints of the specified model, the univariate \eqn{R^2} and adjusted \eqn{R^2} statistics for each PCO score, and the multivariate \eqn{R^2} and adjusted \eqn{R^2} statistics.
#'
#' @details `modelperf()` operates on a single model identified by breakpoints and whether the model is continuous or not. When `x` is a `regions_pco` object, the model is selected either as the best model in the supplied `modelsupport` object (where "best" is determined by the arguments to `criterion` and `model`) or as specified by the user using the arguments to `bps` and `cont`. When `x` is a `regions_results_single` object, the breakpoints and model form are determined based on the supplied object.
#'
#' @seealso [modelsupport()] for assessing model support using information criteria; [calcmodel()] for fitting a single segmented regression model; [plotsegreg()] for plotting the results of a single segmented regression model.
#'
#' @example man/examples/example-modelperf.R
#'

#' @rdname modelperf
#' @export
modelperf <- function(x, ...) {
  UseMethod("modelperf")
}

#' @exportS3Method modelperf regions_pco
#' @rdname modelperf
modelperf.regions_pco <- function(x, scores,
                                  modelsupport = NULL, criterion = "aic", model = 1,
                                  bps = NULL, cont = TRUE, ...) {
  chk::chk_not_missing(scores, "`scores`")
  chk::chk_whole_numeric(scores)
  chk::chk_range(scores, c(1, ncol(x[["scores"]])))

  Xvar <- .get_pos(x)
  Yvar <- x[["scores"]][, scores, drop = FALSE]

  if (!is.null(bps)) {
    if (chk::vld_atomic(bps) && all(is.na(bps))) {
      bps <- NA_real_
      cont <- TRUE
    }
    else {
      chk::chk_numeric(bps)
      chk::chk_range(bps, range(Xvar))
      chk::chk_flag(cont)
    }

    BPs <- sort(.drop_na(bps))

    names(BPs) <- paste0("breakpoint", seq_along(BPs))
  }
  else if (!is.null(modelsupport)) {
    chk::chk_is(modelsupport, "regions_modelsupport")

    chk::chk_string(criterion)
    criterion <- tolower(criterion)
    criterion <- .match_arg(criterion, c("aic", "bic"))
    model_support_crit <- modelsupport[[switch(criterion, aic = "Model_support", bic = "Model_support_BIC")]]

    if (is.null(model)) {
      model <- 1L
    }
    else {
      chk::chk_whole_number(model)
      chk::chk_range(model, c(1, nrow(model_support_crit)))
    }
    cont <- attr(modelsupport, "cont")

    keep <- which(startsWith(names(model_support_crit), "breakpoint"))
    BPs <- .drop_na(unlist(model_support_crit[model, keep]))
  }
  else {
    chk::err("`bps` or `modelsupport` argument must be provided")
  }

  .modelperf_internal(Xvar, Yvar, BPs, cont)
}

#' @exportS3Method modelperf regions_sim
#' @rdname modelperf
modelperf.regions_sim <- function(x, scores = NULL,
                                  modelsupport = NULL, criterion = "aic", model = 1,
                                  bps = NULL, cont = TRUE, ...) {
  Xvar <- x$Xvar
  Yvar <- x$Yvar

  if (!is.null(scores)) {
    chk::chk_whole_numeric(scores)
    chk::chk_range(scores, c(1, ncol(Yvar)))

    Yvar <- Yvar[, scores, drop = FALSE]
  }

  if (!is.null(bps)) {
    if (chk::vld_atomic(bps) && all(is.na(bps))) {
      bps <- NA_real_
      cont <- TRUE
    }
    else {
      chk::chk_numeric(bps)
      chk::chk_range(bps, range(Xvar))
      chk::chk_flag(cont)
    }

    BPs <- sort(.drop_na(bps))

    names(BPs) <- paste0("breakpoint", seq_along(BPs))
  }
  else if (!is.null(modelsupport)) {
    chk::chk_is(modelsupport, "regions_modelsupport")

    chk::chk_string(criterion)
    criterion <- tolower(criterion)
    criterion <- .match_arg(criterion, c("aic", "bic"))
    model_support_crit <- modelsupport[[switch(criterion, aic = "Model_support", bic = "Model_support_BIC")]]

    if (is.null(model)) {
      model <- 1L
    }
    else {
      chk::chk_whole_number(model)
      chk::chk_range(model, c(1, nrow(model_support_crit)))
    }
    cont <- attr(modelsupport, "cont")

    keep <- which(startsWith(names(model_support_crit), "breakpoint"))
    BPs <- .drop_na(unlist(model_support_crit[model, keep]))
  }
  else {
    chk::err("`bps` or `modelsupport` argument must be provided")
  }

  .modelperf_internal(Xvar, Yvar, BPs, cont)
}

#' @exportS3Method modelperf regions_results_single
#' @rdname modelperf
modelperf.regions_results_single <- function(x, scores = NULL, ...) {
  Xvar <- attr(x, "pos")
  Yvar <- attr(x, "scores")
  BPs <- unlist(x$results[startsWith(names(x$results), "breakpoint")])
  cont <- attr(x, "cont")

  if (!is.null(scores)) {
    chk::chk_whole_numeric(scores)
    chk::chk_range(scores, c(1, ncol(Yvar)))

    Yvar <- Yvar[, scores, drop = FALSE]
  }

  .modelperf_internal(Xvar, Yvar, BPs, cont)
}

#' @exportS3Method print regions_perf
print.regions_perf <- function(x, digits = 3, ...) {
  x0 <- x
  cat("Breakpoints:", if (length(x[["BPs"]]) == 0) "(none)"
      else paste(x[["BPs"]], collapse = ", "), "\n")
  cat("\n- Univariate:\n")
  colnames(x$univariate) <- c("R\u00B2", "Adj. R\u00B2")
  print(round(x$univariate, digits))
  cat("\n- Multivariate:\n")
  names(x$multivariate) <- c("R\u00B2", "Adj. R\u00B2")
  print(round(x$multivariate, digits))

  invisible(x0)
}

.modelperf_internal <- function(Xvar, Yvar, BPs, cont) {

  #Calculate weights to ensure each vertebra counts equally
  vert_tab <- tabulate(Xvar)
  w <- 1/vert_tab[Xvar]

  totrsq <- do.call("rbind", lapply(seq_len(ncol(Yvar)), function(i) {
    x <- .design_matrix(Xvar, BPs, cont)
    fit <- .fast_lm(x = x, y = Yvar[,i], w = w)

    n <- length(Xvar)
    fitted <- Yvar[,i] - fit$residuals
    mss <- sum(w * (fitted - sum(w * fitted)/sum(w))^2)
    rss <- sum(w * fit$residuals^2)
    rdf <- n - fit$rank
    r.squared <- mss/(mss + rss)

    c(r2 =     r.squared,
      r2.adj = 1 - (1 - r.squared) * ((n - 1)/rdf),
      SSres =  rss,
      df =     fit$rank + rdf - 1,
      SStot =  mss + rss,
      dfe =    rdf)
  }))
  rownames(totrsq) <- paste0("PCO.", seq_len(ncol(Yvar)))

  # Calculate R2 for all PCOs multivariately:
  tot.rsq <- colSums(totrsq[, -(1:2), drop = FALSE])
  ord.tot.rsq <- unname(1 - tot.rsq["SSres"] / tot.rsq["SStot"])
  adj.tot.rsq <- unname(1 - (tot.rsq["SSres"] / tot.rsq["dfe"]) / (tot.rsq["SStot"] / tot.rsq["df"]))

  out <- list(BPs = BPs,
              univariate = totrsq[,1:2],
              multivariate = c(r2 = ord.tot.rsq, r2.adj = adj.tot.rsq))
  attr(out, "cont") <- cont

  class(out) <- "regions_perf"

  out
}
