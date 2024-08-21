#' Plot a segmented regression model
#'
#' `plotsegreg()` plots the fitted lines resulting from a segmented regression model.
#'
#' @inheritParams modelperf
#' @param scores `numeric`; the indices of the PCO scores for which the fitted lines should be plotted.
#' @param model `numeric`; for which model among the best as determined by `criterion` should fitted lines be plotted. 1 is the best model, 2 the second best, etc. Default is 1. When `x` is a `regions_pco` object, either `modelsupport`, `criterion`, and `model` must be supplied or `bps` and `cont` must be supplied. See Details.
#' @param bps `numeric`; a vector of breakpoints for which model fitted lines should be plotted. When `x` is a `regions_pco` object, either `modelsupport`, `criterion`, and `model` must be supplied or `bps` and `cont` must be supplied. See Details.
#' @param \dots ignored.
#'
#' @return A `ggplot` object that can be manipulated using *ggplot2* syntax.
#'
#' @details `plotsegreg()` operates on a single model identified by breakpoints and whether the model is continuous or not. When `x` is a `regions_pco` object, the model is selected either as the best model in the supplied `modelsupport` object (where "best" is determined by the arguments to `criterion` and `model`) or as specified by the user using the arguments to `bps` and `cont`. When `x` is a `regions_results_single` object, the breakpoints and model form are determined based on the supplied object.
#'
#' `plot()` is an alias for `plotsegreg()` for `regions_results_single` objects.
#'
#' @seealso [modelsupport()] for assessing model support using information criteria; [calcmodel()] for fitting a single segmented regression model; [modelperf()] for computing fit statistics for a single segmented regression model.
#'
#' @example man/examples/example-modelperf.R
#'

#' @export
plotsegreg <- function(x, scores, ...) {
  chk::chk_not_missing(scores, "`scores`")

  UseMethod("plotsegreg")
}

#' @exportS3Method plotsegreg regions_pco
#' @rdname plotsegreg
plotsegreg.regions_pco <- function(x, scores, modelsupport = NULL, criterion = "aic",
                                   model = 1, bps = NULL, cont = TRUE, ...) {
  chk::chk_not_missing(scores, "`scores`")
  chk::chk_whole_numeric(scores)
  chk::chk_range(scores, c(1, ncol(x[["scores"]])))
  scores <- sort(scores)

  Xvar <- .get_pos(x)
  Yvar <- x[["scores"]][, scores, drop = FALSE]

  #Calculate weights to ensure each vertebra counts equally
  vert_tab <- tabulate(Xvar)
  w <- 1/vert_tab[Xvar]

  if (!is.null(bps)) {
    if (!is.null(modelsupport)) {
      .wrn_immediate("`bps` specified; ignoring `modelsupport`")
    }

    chk::chk_numeric(bps)
    chk::chk_range(bps, range(Xvar))
    chk::chk_flag(cont)

    BPs <- sort(.drop_na(bps))
    nbp <- length(BPs)

    if (nbp > 0) {
      names(BPs) <- paste0("breakpoint", seq_len(nbp))
    }
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
    BPs <- NULL
    cont <- TRUE
  }

  fit <- .fast_lm(x = .design_matrix(Xvar, BPs, cont), y = Yvar, w = w)

  yhat <- Yvar - fit$residuals

  .plotreg_internal(Xvar, Yvar, yhat, BPs, lines = TRUE, scores = scores, weights = w)
}

#' @exportS3Method plotsegreg regions_sim
#' @rdname plotsegreg
plotsegreg.regions_sim <- function(x, scores, modelsupport = NULL, criterion = "aic",
                                   model = 1, bps = NULL, cont = TRUE, ...) {

  chk::chk_not_missing(scores, "`scores`")
  chk::chk_whole_numeric(scores)
  chk::chk_range(scores, c(1, ncol(x$Yvar)))

  Xvar <- x$Xvar
  Yvar <- x$Yvar[, scores, drop = FALSE]

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

  fit <- .fast_lm(x = .design_matrix(Xvar, BPs, cont), y = Yvar)

  yhat <- Yvar - fit$residuals

  .plotreg_internal(Xvar, Yvar, yhat, BPs, lines = TRUE, scores = scores)
}

#' @exportS3Method plotsegreg regions_results_single
#' @rdname plotsegreg
plotsegreg.regions_results_single <- function(x, scores, ...) {

  chk::chk_not_missing(scores, "`scores`")
  chk::chk_whole_numeric(scores)
  chk::chk_range(scores, c(1, ncol(attr(x, "scores"))))

  Xvar <- attr(x, "pos")
  Yvar <- attr(x, "scores")[, scores, drop = FALSE]

  #Calculate weights to ensure each vertebra counts equally
  vert_tab <- tabulate(Xvar)
  w <- 1/vert_tab[Xvar]

  BPs <- unlist(x$results[startsWith(names(x$results), "breakpoint")])
  cont <- attr(x, "cont")

  fit <- .fast_lm(x = .design_matrix(Xvar, BPs, cont), y = Yvar, w = w)

  yhat <- Yvar - fit$residuals

  .plotreg_internal(Xvar, Yvar, yhat, BPs, lines = TRUE, scores = scores, weights = w)
}

.plotreg_internal <- function(Xvar, Yvar, yhat = NULL, BPs = NULL, lines = TRUE, scores = 1,
                              linescolor = "darkturquoise", BPcolor = "coral", specimen = NULL,
                              weights = NULL) {

  plot_data <- data.frame(
    PCO = factor(rep(paste("PCO", scores), each = length(Xvar))),
    Xvar = rep(Xvar, ncol(Yvar))
  )

  use_specimen <- !is.null(specimen) && nlevels(specimen) > 1
  if (use_specimen) {
    plot_data$specimen <- rep(specimen, ncol(Yvar))
  }

  plot_data$weights <- {
    if (is.null(weights)) 1
    else rep(weights, ncol(Yvar))
  }


  #Flatten scores/predicted values
  plot_data$Yvar <- as.vector(Yvar)
  if (!is.null(yhat)) plot_data$yhat <- as.vector(yhat)

  p <- ggplot(plot_data,
              aes(x = .data$Xvar))

  if (use_specimen)
    p <- p + geom_point(aes(y = .data$Yvar,
                            color = .data$specimen,
                            size = .data$weights),
                        shape = "circle")
  else
    p <- p + geom_point(aes(y = .data$Yvar,
                            size = .data$weights),
                        shape = "circle",
                        color = "darkgray")

  if (length(BPs) > 0) {
    breakpoints <- BPs + .5
    p <- p + geom_vline(xintercept = breakpoints,
                        color = BPcolor,
                        linetype = "dashed")
  }
  else {
    breakpoints <- numeric(0)
  }

  if (!is.null(yhat) && lines) {
    p <- p + geom_line(aes(y = .data$yhat,
                           group = cut(.data$Xvar, c(-Inf, breakpoints, Inf))),
                       show.legend = FALSE, color = linescolor)
  }

  p <- p + scale_size_area(max_size = 2, guide = NULL) +
    facet_wrap(vars(.data$PCO), ncol = 1,
               scales = "free_y",
               strip.position = "left") +
    theme_bw() +
    theme(strip.placement = "outside",
          strip.background = element_blank(),
          strip.text = element_text(size = 11)) +
    labs(x = "Position",
         y = NULL,
         color = NULL)

  p
}
