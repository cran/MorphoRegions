#' Fit segmented regression models for all combinations of breakpoints
#'
#' `calcregions()` enumerates all possible combinations of breakpoints to fit multivariate segmented regression models. `addregions()` adds models with additional numbers of regions to the resulting output object. `ncombos()` computes an upper bound on the number of breakpoint combinations that will be tested.
#'
#' @param pco a `regions_pco` object; the output of a call to [svdPCO()].
#' @param scores `numeric`; the indices of the PCO scores to use as the outcomes in fitting the models (e.g., `1:4` to use the first four scores). Can also be the output of a call to [PCOselect()].
#' @param noregions `numeric`; for `calcregions()`, the maximum number of regions for which models are fit (e.g, 4 to request models with 1 to 4 regions); for `addregions()`, a vector containing the numbers of regions to add (e.g., 5:6 to request models with 5 and 6 regions); for `ncombos()`, a vector containing the numbers of regions to check.
#' @param minvert `numeric`; the minimum number of vertebrae allowed in each region. Default is 3.
#' @param cont `logical`; whether to fit models that are continuous (`TRUE`) or discontinuous (`FALSE`) at the breakpoints. Default is `TRUE`.
#' @param exhaus `logical`; whether to fit all possible models (`TRUE`) or use heuristics to reduce the number of models fit (`FALSE`). Default is `TRUE`. See Details. Setting to `FALSE` can reduce the size of the resulting object.
#' @param includebp an optional vector of vertebrae that must be included in any tested set of breakpoints, e.g., if it is known that two regions are divided at that vertebra. `includebp` does not have to obey the `minvert` rules, but a warning will be thrown if it doesn't.
#' @param omitbp an optional vector of vertebrae to be omitted from the list of possible breakpoints, e.g., if it is known that two adjacent vertebrae belong to the same region.
#' @param ncombos_file_trigger `numeric`; when the number of eligible combinations of breakpoints exceeds this number, the problem will be split into smaller problems, with the results of each stored in its own temporary file in the directory supplied to `temp_file_dir` before being re-read into memory. The primary purpose of this is to preserve memory when `exhaus = FALSE` by delegating storage of the results to disk instead of RAM.
#' @param temp_file_dir string; the directory where the temporary files will be saved (and then deleted) when the number of breakpoint combinations exceeds `ncombos_file_trigger`. Default is the directory produced by [tempdir()], but it is much safer to provide your own directory, which must already exist on your machine. See Details.
#' @param cl a cluster object created by [parallel::makeCluster()], an integer to indicate number of child-processes (integer values are ignored on Windows) for parallel evaluations, or `"future"` to use a future backend. `NULL` (the default) refers to sequential evaluation (no parallelization). See [pbapply::pbapply()] for details.
#' @param verbose `logical`; whether to print information about the fitting process, including a progress bar. Default is `TRUE`.
#' @param regions_results,object a `regions_results` object; the output of a call to `calcregions()` or `addregions()`.
#' @param \dots ignored.
#'
#' @returns A `regions_results` object with the following components:
#' * `results` - the results of the fitting process for each combination of breakpoints
#' * `stats` - statistics summarizing the fitting process. Use `summary()` to view this information in a clean format.
#'
#' `ncombos()` returns a numeric vector with the number of breakpoint combinations for each number of regions (which are stored as the names).
#'
#' @details `calcregions()` enumerates all possible combinations of breakpoints that satisfy the constraint imposed by `minvert` (i.e., that breakpoints need to be at least `minvert` vertebrae apart) and fits the segmented regression models implied by each combination. These are multivariate regression models with the PCO scores specified by `scores` as the outcomes. When `cont = TRUE`, these regression models are continuous; i.e., the regression lines for each region connect at the breakpoints. Otherwise, the models are discontinuous so that each region has its own intercept and slope. The models are fit using [.lm.fit()], which efficiently implements ordinary least squares regression.
#'
#' When `exhaus = FALSE`, heuristics are used to reduce the number of models to fit, which can be useful for keeping the size of the resulting object down by avoiding fitting models corresponding to breakpoint combinations that yield a poor fit to the data. Only breakpoint combinations that correspond to the breakpoints of the best fitting model with a smaller number of regions +/- 3 vertebrae are used, and only models that have an RSS smaller than half a standard deviation more the smallest RSS are kept.
#'
#' `addregions()` should be used on an existing `regions_results` object to add models with more regions. Internally, it works just the same as `calcregions()`.
#'
#' `ncomobs()` computes an upper bound on the number of possible breakpoint combinations. When `exhaus = FALSE` or `includebp` is specified, the actual number of combinations will be smaller than that produced by `ncombos()`.
#'
#' When the number of possible combinations of breakpoints for a given number of regions (as computed by `ncombos()`) is larger than `ncombos_file_trigger`, the problem will be split into smaller problems, with the results of each stored in temporary files that are deleted when the function completes. These temporary files will be stored in the directory supplied to `temp_file_dir`. By default, this is the temporary directory produced by [tempdir()]. However, this directory can be deleted by R at any time without warning, which will cause the function to crash, so it is a good idea to supply your own directory that will be preserved. You can use `ncombos()` to check to see if the number of breakpoint combinations exceeds `ncombos_file_trigger`.
#'
#' @seealso [calcmodel()] to fit a segmented regression model for a single set of breakpoints; [modelselect()] to select the best model for each number of regions based on RSS; [modelsupport()] to compute statistics the describe the support of the best models; [calcBPvar()] to compute the variability in the optimal breakpoints.
#'
#' @example man/examples/example-calcregions.R

#' @export
calcregions <- function(pco, scores, noregions, minvert = 3, cont = TRUE,
                        exhaus = TRUE, includebp = NULL, omitbp = NULL,
                        ncombos_file_trigger = 1e7, temp_file_dir = tempdir(TRUE),
                        cl = NULL, verbose = TRUE) {
  # Argument checks
  if (inherits(pco, "regions_pco")) {
    chk::chk_not_missing(scores, "`scores`")
    chk::chk_whole_numeric(scores)
    chk::chk_range(scores, c(1, ncol(pco[["scores"]])))

    Xvar <- .get_pos(pco)
    Yvar <- pco[["scores"]][, scores, drop = FALSE]

    eligible_vertebrae <- .get_eligible_vertebrae(pco)
  }
  else if (inherits(pco, "regions_sim")) {
    if (missing(scores)) {
      scores <- seq_len(ncol(pco[["Yvar"]]))
    }
    else {
      chk::chk_whole_numeric(scores)
      chk::chk_range(scores, c(1, ncol(pco[["Yvar"]])))
    }

    Xvar <- pco[["Xvar"]]
    Yvar <- pco[["Yvar"]][, scores, drop = FALSE]

    eligible_vertebrae <- sort(unique(Xvar))
  }
  else {
    chk::err("`pco` must be a `regions_pco` or `regions_sim` object")
  }

  chk::chk_not_missing(noregions, "`noregions`")
  chk::chk_count(noregions)
  chk::chk_gte(noregions, 1)

  chk::chk_count(minvert)
  chk::chk_gte(minvert, 2)
  chk::chk_flag(cont)
  chk::chk_flag(exhaus)
  chk::chk_flag(verbose)

  if (!is.null(includebp)) {
    chk::chk_whole_numeric(includebp)

    if (!all(includebp %in% eligible_vertebrae)) {
      chk::err("all breakpoints specified in `includebp` must correspond to measured vertebrae")
    }

    if (length(includebp) > 1) {
      includebp <- sort(includebp)

      if (sum(eligible_vertebrae < min(includebp)) < 1 ||
          sum(eligible_vertebrae > max(includebp)) < 2) {
        chk::err("breakpoints specified in `includebp` must be be at least 2 measured vertebrae from the ends")
      }

      for (i in seq_along(includebp)[-1]) {
        if (sum(eligible_vertebrae > includebp[i - 1] &
                eligible_vertebrae <= includebp[i]) < 2) {
          chk::err("breakpoints specified in `includebp` must have at least 2 measured vertebrae between them")
        }
      }

      if (length(includebp) >= noregions) {
        chk::err("the length of `includebp` must be less than `noregions`")
      }

      if (any(diff(c(eligible_vertebrae[1] - 1, includebp, max(eligible_vertebrae))) < minvert)) {
        .wrn_immediate("the breakpoints specified in `includebp` do not obey the `minvert` crtierion")
      }
    }
  }

  eligible_noregions <- {
    if (is.null(includebp)) seq_len(noregions)
    else seq(length(includebp) + 1, noregions)
  }

  if (!is.null(omitbp)) {
    chk::chk_whole_numeric(omitbp)
    if (any(omitbp %in% includebp)) {
      chk::err("breakpoints specified in `omitbp` cannot also be present in `includebp`")
    }
    eligible_vertebrae <- eligible_vertebrae[!eligible_vertebrae %in% omitbp]
  }

  # Re-order data
  o <- order(Xvar)
  Xvar <- Xvar[o]
  Yvar <- Yvar[o, , drop = FALSE]

  ncombs <- vapply(eligible_noregions, function(i) {
    .ncombos(length(eligible_vertebrae), minvert, i - 1)
  }, numeric(1L))

  if (any(ncombs == 0)) {
    max.regions <- max(which(ncombs > 0))
    .wrn_immediate(sprintf("models for %s vertebrae cannot be fit with more than %s region%%s and `minvert` of %s",
                           length(eligible_vertebrae), max.regions, minvert), n = max.regions)
    ncombs <- ncombs[eligible_noregions <= max.regions]
    eligible_noregions <- eligible_noregions[eligible_noregions <= max.regions]
  }

  if (any(eligible_noregions > 3)) {
    chk::chk_count(ncombos_file_trigger)
    chk::chk_gt(ncombos_file_trigger, 0)
    if (any(ncombs[eligible_noregions > 3] > ncombos_file_trigger)) {
      if (!"temp_file_dir" %in% names(match.call())[-1] ||
          identical(temp_file_dir, tempdir())) {
        fn <- "calcregions"
        # .wrn_immediate(sprintf("`%s()` will use temporary files stored in the directory supplied to `temp_file_dir`. Because this was unspecified or left as its default, you may see the error \"cannot write to connection\" after a long run of this function. See `help(\"%s\")` for solutions and how to avoid this problem", fn, fn))
      }
      else {
        chk::chk_string(temp_file_dir)
      }
    }
  }

  opb <- pbapply::pboptions(type = if (verbose) "timer" else "none")
  on.exit(pbapply::pboptions(opb))

  if (inherits(cl, "cluster")) {
    parallel::clusterExport(cl,
                            unclass(utils::lsf.str(envir = asNamespace(utils::packageName()),
                                                   all.names = TRUE)),
                            envir = getNamespace(utils::packageName()))
  }

  # Prep output object
  out <- list(results = NULL, stats = NULL)

  attr(out, "scores") <- Yvar
  attr(out, "cont") <- cont
  attr(out, "minvert") <- minvert
  attr(out, "includebp") <- includebp
  attr(out, "omitbp") <- omitbp
  attr(out, "pos") <- Xvar
  attr(out, "eligible_vertebrae") <- eligible_vertebrae

  class(out) <- "regions_results"

  for (noregion in eligible_noregions) {
    out <- .addregions_internal(out, noregion, exhaus = exhaus,
                                ncombos_file_trigger = ncombos_file_trigger,
                                temp_file_dir = temp_file_dir,
                                cl = cl, verbose = verbose)
  }

  out
}

#' @export
#' @rdname calcregions
addregions <- function(regions_results, noregions, exhaus = TRUE,
                       ncombos_file_trigger = 1e7, temp_file_dir = tempdir(TRUE),
                       cl = NULL, verbose = TRUE) {
  # Argument checks
  chk::chk_is(regions_results, "regions_results")

  chk::chk_not_missing(noregions, "`noregions`")
  chk::chk_whole_numeric(noregions)
  chk::chk_gte(noregions, 1)
  noregions <- sort(noregions)

  chk::chk_flag(exhaus)
  chk::chk_flag(verbose)

  opb <- pbapply::pboptions(type = if (verbose) "timer" else "none")
  on.exit(pbapply::pboptions(opb))

  eligible_vertebrae <- attr(regions_results, "eligible_vertebrae")

  minvert <- attr(regions_results, "minvert")

  ncombs <- vapply(noregions, function(i) {
    .ncombos(length(eligible_vertebrae), minvert, i - 1)
  }, numeric(1L))

  if (any(ncombs == 0)) {
    max.regions <- max(noregions[ncombs > 0])
    .wrn_immediate(sprintf("models for %s vertebrae cannot be fit with more than %s region%%s and `minvert` of %s",
                           length(eligible_vertebrae), max.regions, minvert), n = max.regions)
    ncombs <- ncombs[noregions <= max.regions]
    noregions <- noregions[noregions <= max.regions]
  }

  if (any(noregions > 3)) {
    chk::chk_count(ncombos_file_trigger)
    chk::chk_gt(ncombos_file_trigger, 0)
    if (any(ncombs[noregions > 3] > ncombos_file_trigger)) {
      if (!"temp_file_dir" %in% names(match.call())[-1] ||
          identical(temp_file_dir, tempdir())) {
        fn <- "addregions"
        # .wrn_immediate(sprintf("`%s()` will use temporary files stored in the directory supplied to `temp_file_dir`. Because this was unspecified or left as its default, you may see the error \"cannot write to connection\" after a long run of this function. See `help(\"%s\")` for solutions and how to avoid this problem", fn, fn))
      }
      else {
        chk::chk_string(temp_file_dir)
      }
    }
  }

  if (inherits(cl, "cluster")) {
    parallel::clusterExport(cl,
                            unclass(utils::lsf.str(envir = asNamespace(utils::packageName()),
                                                   all.names = TRUE)),
                            envir = getNamespace(utils::packageName()))
  }

  for (n in noregions) {
    regions_results <- .addregions_internal(regions_results, noregions = n,
                                            exhaus = exhaus,
                                            ncombos_file_trigger = ncombos_file_trigger,
                                            temp_file_dir = temp_file_dir,
                                            cl = cl, verbose = verbose)
  }

  regions_results
}

#' @exportS3Method print regions_results
print.regions_results <- function(x, ...) {
  cat("A `regions_results` object\n")
  cat(" - number of PCOs used:", ncol(attr(x, "scores")), "\n")
  cat(" - number of regions:", paste(sort(unique(x$stats[["Nregions"]])), collapse = ", "), "\n")
  cat(" - model type:", if (attr(x, "cont")) "continuous" else "discontinuous", "\n")
  cat(" - min vertebrae per region:", attr(x, "minvert"), "\n")
  if (!is.null(attr(x, "omitbp"))) {
    cat(" - omitted breakpoints:", paste(attr(x, "omitbp"), collapse = ", "), "\n")
  }
  cat(" - total models saved:", nrow(x$results), "\n")
  cat("Use `summary()` to examine summaries of the fitting process.\n")

  invisible(x)
}

#' @exportS3Method summary regions_results
#' @rdname calcregions
summary.regions_results <- function(object, ...) {
  stats <- object$stats[1:6]
  names(stats) <- c("Regions", "Possible", "Tested", "Saved",
                    "Comp. method", "Saving method")

  class(stats) <- c("summary.regions_results", class(stats))

  stats
}

#' @exportS3Method print summary.regions_results
print.summary.regions_results <- function(x, ...) {
  print.data.frame(x, row.names = FALSE)

  invisible(x)
}

#' @export
#' @rdname calcregions
ncombos <- function(pco, noregions, minvert = 3, includebp = NULL, omitbp = NULL) {
  # Argument checks
  if (inherits(pco, "regions_pco")) {
    eligible_vertebrae <- .get_eligible_vertebrae(pco)
  }
  else if (inherits(pco, "regions_sim")) {
    Xvar <- pco[["Xvar"]]

    eligible_vertebrae <- sort(unique(Xvar))
  }
  else {
    chk::err("`pco` must be a `regions_pco` or `regions_sim` object")
  }

  chk::chk_not_missing(noregions, "`noregions`")
  chk::chk_whole_numeric(noregions)
  chk::chk_gte(noregions, 1)
  noregions <- sort(noregions)

  chk::chk_count(minvert)
  chk::chk_gte(minvert, 2)

  if (!is.null(includebp)) {
    chk::chk_whole_numeric(includebp)

    if (!all(includebp %in% eligible_vertebrae)) {
      chk::err("all breakpoints specified in `includebp` must correspond to measured vertebrae")
    }

    if (length(includebp) > 1) {
      includebp <- sort(includebp)

      if (sum(eligible_vertebrae < min(includebp)) < 1 ||
          sum(eligible_vertebrae > max(includebp)) < 2) {
        chk::err("breakpoints specified in `includebp` must be be at least 2 measured vertebrae from the ends")
      }

      for (i in seq_along(includebp)[-1]) {
        if (sum(eligible_vertebrae > includebp[i - 1] &
                eligible_vertebrae <= includebp[i]) < 2) {
          chk::err("breakpoints specified in `includebp` must have at least 2 measured vertebrae between them")
        }
      }

      if (length(includebp) >= noregions) {
        chk::err("the length of `includebp` must be less than `noregions`")
      }

      if (any(diff(c(eligible_vertebrae[1] - 1, includebp, max(eligible_vertebrae))) < minvert)) {
        .wrn_immediate("the breakpoints specified in `includebp` do not obey the `minvert` crtierion")
      }
    }
  }

  eligible_noregions <- {
    if (is.null(includebp)) noregions
    else noregions[noregions >= length(includebp)]
  }

  if (!is.null(omitbp)) {
    chk::chk_whole_numeric(omitbp)
    if (any(omitbp %in% includebp)) {
      chk::err("breakpoints specified in `omitbp` cannot also be present in `includebp`")
    }
    eligible_vertebrae <- eligible_vertebrae[!eligible_vertebrae %in% omitbp]
  }

  ncombs <- vapply(eligible_noregions, function(i) {
    .ncombos(length(eligible_vertebrae), minvert, i - 1)
  }, numeric(1L))

  setNames(ncombs, eligible_noregions)
}

# Does the actual work of fitting the segmented regression models
.addregions_internal <- function(regions_results, noregions, exhaus = FALSE,
                                 ncombos_file_trigger = 1e6,
                                 temp_file_dir = tempdir(TRUE),
                                 cl = NULL, verbose = TRUE) {

  minvert <- attr(regions_results, "minvert")

  cont <- attr(regions_results, "cont")

  Xvar <- attr(regions_results, "pos")
  Yvar <- attr(regions_results, "scores")
  eligible_vertebrae <- attr(regions_results, "eligible_vertebrae")

  #Calculate weights to ensure each vertebra counts equally
  vert_tab <- tabulate(Xvar)
  w <- 1/vert_tab[Xvar]

  includebp <- attr(regions_results, "includebp")

  nbp <- noregions - 1

  noPC <- ncol(Yvar)

  colhead <- c("regions", paste0("breakpoint", seq_len(max(1, nbp))),
               "sumRSS", paste("RSS", 1:noPC, sep = "."))

  # Set verbosity of progress bar
  if (!verbose) {
    opb <- pbapply::pboptions(type = "none")
    on.exit(pbapply::pboptions(opb))
  }

  # Function to fit regression based on BPs
  fregions <- function(BPs, Xvar, Yvar, w) {
    x <- .design_matrix(Xvar, BPs, cont)

    lines <- .fast_lm(x = x, y = Yvar, w = w)
    RSS <- sum(w * lines$residuals^2)

    rsq <- {
      if (noPC > 1) .colSums(w * lines$residuals^2,
                             m = length(Xvar), n = ncol(Yvar))
      else RSS
    }

    c((length(BPs) + 1), BPs, RSS, rsq)
  }

  if (nbp == 0) {
    # 1 Regions (0 bp)

    if (!is.null(includebp)) return(regions_results)

    if (verbose) {
      cat("Fitting model with 1 region (0 breakpoints)...\n")
    }

    lines <- .fast_lm(x = cbind(1, Xvar), y = Yvar, w = w)

    RSS <- sum(w * lines$residuals^2)

    rsq <- {
      if (noPC > 1) .colSums(w * lines$residuals^2,
                             m = length(Xvar), n = ncol(Yvar))
      else RSS
    }

    res <- as.data.frame(rbind(c(nbp + 1, NA_real_, RSS, rsq)))

    names(res) <- colhead

    bpkeep <- NA_character_

    stats <- data.frame(Nregions = nbp + 1L,
                        Nmodel_possible = 1L,
                        Nmodel_tested = 1L,
                        Nmodel_saved = 1L,
                        Comp_method = "Exhaustive",
                        Saving_method = "All")

    best_bps <- matrix(bpkeep, nrow = 1, dimnames = list(NULL, "Best_BPs1"))

    if (verbose) {
      cat("Done. 1 model tested, 1 model saved.\n")
    }
  }
  else if (nbp == 1) {
    # 2 Regions (1 bp)

    # Get all possible combinations for 1 breakpoint:
    goodcomb <- .combos(eligible_vertebrae, minvert, nbp, includebp = includebp)

    if (nrow(goodcomb) == 0) return(regions_results)

    if (verbose) {
      cat("Fitting models with 2 regions (1 breakpoint)...\n")
    }

    # Run fitting on good combinations:
    res <- do.call("rbind", pbapply::pbapply(goodcomb, 1, function(bps) {
      fregions(bps, Xvar = Xvar, Yvar = Yvar, w = w)
    },
    simplify = FALSE, cl = cl))

    res <- setNames(as.data.frame(res), colhead)

    bpkeep <- {
      NA_character_
      # if (!exhaus) paste(goodcomb[,1], collapse = "|")
      # else NA_character_
    }

    stats <- data.frame(Nregions = nbp + 1L,
                        Nmodel_possible = nrow(res),
                        Nmodel_tested = nrow(res),
                        Nmodel_saved = nrow(res),
                        Comp_method = "Exhaustive",
                        Saving_method = "All")

    best_bps <- matrix(bpkeep, nrow = 1, dimnames = list(NULL, "Best_BPs1"))

    if (verbose) {
      cat(sprintf("Done. %s models tested, %s models saved.\n",
                  nrow(res), nrow(res)))
    }
  }
  else {
    # 3 and more Regions (2+ bp)

    # Define BPs to include if non-exhaustive search:
    if (!exhaus) {

      stats <- regions_results$stats
      max_nregion_so_far <- max(stats$Nregions[stats$Nregions < noregions])
      best_bps <- unlist(stats[stats$Nregions == max_nregion_so_far,
                               startsWith(names(stats), "Best_BPs")])
      bpkeep <- .drop_na(best_bps)

      if (length(bpkeep) == 0 && max_nregion_so_far > 2) {
        # Find best_bps
        res0 <- regions_results$res[regions_results$res$regions == max_nregion_so_far,, drop = FALSE]
        bestBPs <- drop(as.matrix(res0)[which.min(res0$sumRSS), startsWith(names(res0), "breakpoint")])

        # Add 3 BPs left and right to best BPs
        bestBPs <- lapply(bestBPs, function(b) {
          bpp <- match(b, eligible_vertebrae) + seq(-3, 3)

          # Convert position of bp to keep to actual value of bp
          eligible_vertebrae[bpp[bpp %in% seq_along(eligible_vertebrae)]]
        })

        bpkeep <- vapply(seq_along(bestBPs), function(i) {
          paste(
            sort(.drop_na(unique(c(res0[, i + 1], bestBPs[[i]]), nmax = length(eligible_vertebrae)))),
            collapse = "|"
          )
        }, character(1L))

        regions_results$stats[stats$Nregions == max_nregion_so_far,
                              startsWith(names(stats), "Best_BPs")] <- bpkeep
      }
    }

    # Number of combinations
    ncombi <- .ncombos(length(eligible_vertebrae), minvert, nbp)

    # Max length for an R object to not exceed 2Gb in size (Memory safety net)
    # lim <- 5e7

    # Find nrow max for given number of cols (bps) to avoid exceeding 2Gb objects
    # max_rows <- ceiling(lim/nbp)

    if (ncombi <= ncombos_file_trigger) {

      goodcomb <- .combos(eligible_vertebrae, minvert, nbp, includebp = includebp)

      if (nrow(goodcomb) == 0) return(regions_results)

      nmodel_possible <- nrow(goodcomb)

      if (!exhaus && length(bpkeep) > 0) {
        # Keep only probable combinations (non-exhaustive search)
        for (m in strsplit(bpkeep, "|", fixed = TRUE)) {
          mn <- as.numeric(m)
          goodcomb <- goodcomb[.any_mat_in(goodcomb, mn),, drop = FALSE]
        }
      }

      if (nrow(goodcomb) == 0) return(regions_results)

      if (verbose) {
        cat(sprintf("Fitting models with %s regions (%s breakpoints)...\n",
                    nbp + 1, nbp))
      }

      # Run fitting on good combinations:
      res <- do.call("rbind", pbapply::pbapply(goodcomb, 1, function(bps) {
        fregions(bps, Xvar = Xvar, Yvar = Yvar, w = w)
      },
      simplify = FALSE, cl = cl))

      res <- setNames(as.data.frame(res), colhead)

      nmodel_tested <- nrow(res)

      if (!exhaus) {
        # If non-exhaustive search, keep only models with sumRSS <= min(sumRSS)+(0.5*sd(sumRSS))
        cutoff <- min(res$sumRSS) + sd(res$sumRSS) / 2
        if (!is.na(cutoff)) {
          # Subsample res only if more than 1 option
          res <- res[res$sumRSS <= cutoff,, drop = FALSE]
        }
      }
    }
    else {

      #Check directory
      if (!dir.exists(temp_file_dir)) {
        dir.create(dir.exists)
      }

      # Max number of rows per iteration; keeps sizes small and
      # improves progress bar
      max_rows <- ceiling(ncombos_file_trigger/100)

      # Required number of splits of combos
      aa <- ceiling(ncombi / max_rows)

      # File names to store intermediate res objects
      files <- character(aa)
      # files <- vapply(1:aa, function(i) tempfile("regions_", fileext = ".Rds"),
      #                 character(1L))
      on.exit(unlink(files))

      # Initialize statistics
      nmodel_possible <- nmodel_tested <- 0
      last_combo <- NULL

      n_total <- 0
      if (!exhaus) {
        min_RSS <- Inf
        Xbar_total <- M2_total <- 0
      }

      for (i in 1:aa) {

        # Only returns next max_rows rows of combos, using last_combo as
        # starting place
        goodcomb <- .combos(eligible_vertebrae, minvert, nbp, ncombos = max_rows,
                            last_combo = last_combo,
                            includebp = includebp)

        if (i == 1) {
          if (nrow(goodcomb) == 0) {
            return(regions_results)
          }

          if (verbose) {
            cat(sprintf("Fitting models with %s regions (%s breakpoints)...\n",
                        nbp + 1, nbp))
          }

          # Initialize progress bar (pb). Won't activate if `verbose = FALSE`
          pb <- pbapply::startpb(0, ncombi)
        }

        last_combo <- goodcomb[nrow(goodcomb),]
        nmodel_possible <- nmodel_possible + nrow(goodcomb)

        if (!exhaus) {
          for (m in strsplit(bpkeep, "|", fixed = TRUE)) {
            mn <- as.numeric(m)
            goodcomb <- goodcomb[.any_mat_in(goodcomb, mn),, drop = FALSE]
          }
        }

        if (nrow(goodcomb) == 0) {
          pbapply::setpb(pb, min(i * max_rows, ncombi))
          next
        }

        res <- do.call("rbind", .lapply_selector(seq_len(nrow(goodcomb)), function(g) {
          fregions(goodcomb[g,], Xvar = Xvar, Yvar = Yvar, w = w)
        }, cl = cl))

        res <- setNames(as.data.frame(res), colhead)
        nmodel_tested <- nmodel_tested + nrow(res)

        # Save intermediate object to disk
        files[i] <- tempfile("regions_", tmpdir = temp_file_dir, fileext = ".Rds")
        saveRDS(res, files[i])

        n_old <- n_total #old n
        n_i <- nrow(res) #n of new subset
        n_total <- n_old + n_i #total n so far

        if (!exhaus) {
          #Accumulate values for cutoff
          Xbar_old <- Xbar_total #old mean
          Xbar_i <- sum(res$sumRSS)/n_i #mean of new subset
          Xbar_total <- (n_old * Xbar_old + n_i * Xbar_i)/(n_old + n_i) #new mean

          delta <- Xbar_i - Xbar_old

          M2_old <- M2_total #old M2
          M2_i <- sum((res$sumRSS - Xbar_i)^2) #M2 of new subset

          M2_total <- M2_old + M2_i + delta^2 * n_old * n_i / (n_old + n_i) #new M2

          min_RSS <- min(min_RSS, res$sumRSS) #new min
        }

        pbapply::setpb(pb, nmodel_possible)
      }

      pbapply::setpb(pb, ncombi)
      pbapply::closepb(pb)

      if (!exhaus) {
        sd_RSS <- sqrt(M2_total/(nmodel_tested - 1))
        cutoff <- min_RSS + sd_RSS / 2
      }

      if (verbose) {
        cat("Reading files...\n")
      }

      files <- files[files != character(1L)]

      res <- do.call("rbind", .lapply_selector(files, function(file) {
        tryCatch({res_tmp <- readRDS(file)},
                 error = function(err) {
                   e <- conditionMessage(err)
                   if (!grepl("cannot open connection", e)) {
                     chk::err(e, tidy = FALSE)
                   }
                   if (!dir.exists(temp_file_dir)) {
                     chk::err("the directory containing the stored files has been deleted. See `help(\"calcregions\")` for more information. You will likely have to re-run the function")
                   }
                   chk::err("at least one of the stored files cannot be found, suggesting that something has happened to the directory in which the files were stored. See `help(\"calcregions\")` for more information. You will likely have to re-run the function")
                 })

        if (exhaus) return(res_tmp)

        # If non-exhaustive search, keep only models with sumRSS <= min(sumRSS)+(0.5*sd(sumRSS))
        res_tmp[res_tmp$sumRSS <= cutoff,, drop = FALSE]

      }, cl = cl))

      unlink(files)
    }

    nmodel_saved <- nrow(res)

    # Create string with info on best BPs to use for next region:	# BPs kept are
    # all  > cutoff (min+0.5*sd) and BPs of best models +/- 3 vertebrae
    if (!exhaus) {
      bestBPs <- drop(res[which.min(res$sumRSS), startsWith(names(res), "breakpoint")])

      # Add 3 BPs left and right to best BPs
      bestBPs <- lapply(bestBPs, function(x) {
        #Convert bp value to position, add range
        bpp <- match(x, eligible_vertebrae) + seq(-3, 3)
        # Convert position of bp to keep to actual value of bp
        eligible_vertebrae[bpp[bpp > 0]]
      })

      bpkeep <- vapply(seq_along(bestBPs), function(i) {
        paste(
          sort(.drop_na(unique(c(res[, i + 1], bestBPs[[i]]), nmax = length(eligible_vertebrae)))),
          collapse = "|"
        )
      }, character(1L))

      stats <- data.frame(Nregions = noregions,
                          Nmodel_possible = nmodel_possible,
                          Nmodel_tested = nmodel_tested,
                          Nmodel_saved = nmodel_saved,
                          Comp_method = "Non-exhaus",
                          Saving_method = "SD/2")
    }
    else {
      bpkeep <- rep(NA_character_, nbp)

      stats <- data.frame(Nregions = noregions,
                          Nmodel_possible = nmodel_possible,
                          Nmodel_tested = nmodel_tested,
                          Nmodel_saved = nmodel_saved,
                          Comp_method = "Exhaustive",
                          Saving_method = "All")
    }

    best_bps <- matrix(bpkeep, nrow = 1, dimnames = list(NULL, paste0("Best_BPs", 1:nbp)))

    if (verbose) {
      cat(sprintf("Done. %s models tested, %s models saved.\n",
                  nmodel_tested, nmodel_saved))
    }
  }

  rownames(res) <- NULL

  stats <- cbind(stats, as.data.frame(best_bps))

  #Remove analysis with same number of regions
  if (any(regions_results$results$regions == noregions)) {
    regions_results$results <- regions_results$results[regions_results$results$regions != noregions,, drop = FALSE]
  }
  if (any(regions_results$stats$Nregions == noregions)) {
    regions_results$stats <- regions_results$stats[regions_results$stats$Nregions != noregions,, drop = FALSE]
  }

  #Combine res and stats with exist components of input object
  regions_results$results <- .rbind_larger(regions_results$results, res)
  regions_results$stats <- .rbind_larger(regions_results$stats, stats)

  regions_results$results <- regions_results$results[order(regions_results$results$regions),, drop = FALSE]
  regions_results$stats <- regions_results$stats[order(regions_results$stats$Nregions),, drop = FALSE]

  rownames(regions_results$results) <- NULL
  rownames(regions_results$stats) <- NULL

  regions_results
}

# Returns number of combos possible; vectorized across nvert
.ncombos <- function(nvert, minvert, nbp) {
  if (nbp == 0) return(rep.int(1, length(nvert)))

  # Number of possibilities for each BP
  p <- nvert - minvert * (nbp + 1) + 1

  out <- rep(0, length(nvert))

  if (all(p <= 0)) return(out)

  # Number of combinations - this is magic, kind of, but it works
  # prod(p - 1 + 1:nbp)/prod(1:nbp)
  out[p > 0] <- choose(p[p > 0] - 1 + nbp, nbp)

  out
}

# Returns all combos satisfying rules by dispatching to .combosR_slow()
# or .combosR_quick() depending on rules
.combos <- function(vert, minvert, nbp, ncombos = Inf, last_combo = NULL,
                    includebp = NULL) {

  maxcombos <- .ncombos(length(vert), minvert, nbp)
  if (!is.null(last_combo) || ncombos < maxcombos || maxcombos == 0 ||
      (!is.null(includebp) && any(diff(includebp) < minvert))) {
    return(.combosR_slow(vert, minvert, nbp, ncombos, last_combo,
                         includebp))
  }

  out <- .combosR_quick(vert, minvert, nbp)

  if (!is.null(includebp)) {
    out <- out[.all_mat_in(out, includebp),, drop = FALSE]
  }

  out
}

# Returns matrix of all combos, uses counter and allows other options
.combosR_slow <- function(vert, minvert, nbp, ncombos = Inf, last_combo = NULL,
                          includebp = NULL) {

  if (!is.null(includebp) && length(includebp) == nbp) {
    combos <- matrix(includebp, nrow = 1, ncol = nbp)

    colnames(combos) <- paste0("breakpoint", seq_len(nbp))
    rownames(combos) <- NULL

    return(combos)
  }

  n <- length(vert)

  y <- seq_len(n)

  if (!is.null(includebp)) {
    includebp <- match(includebp, vert)
    minvert_override <- min(includebp, diff(includebp), minvert)
  }
  else {
    minvert_override <- minvert
  }

  mins <- minvert_override * seq_len(nbp)
  maxes <- n - rev(mins)

  maxcombos <- .ncombos(n, minvert_override, nbp)

  ncombos <- min(ncombos, maxcombos)
  combos <- matrix(NA_integer_, nrow = ncombos, ncol = max(nbp, 1))

  colnames(combos) <- paste0("breakpoint", seq_len(max(nbp, 1)))
  rownames(combos) <- NULL

  if (ncombos == 0 || nbp == 0) return(combos)

  reset0 <- rep(FALSE, nbp)

  if (is.null(last_combo)) {
    counters <- mins
  }
  else {
    counters <- match(last_combo, vert)

    #Increment combo before inserting into combos matrix
    reset <- reset0
    counters[nbp] <- counters[nbp] + 1

    i <- nbp
    while (i > 1 && counters[i] > maxes[i]) {
      reset[i] <- TRUE
      counters[i - 1] <- counters[i - 1] + 1
      i <- i - 1
    }
    if (any(reset)) {
      for (i in which(reset)) {
        counters[i] <- counters[i - 1] + minvert_override
      }
    }
  }

  k <- 1

  while (k <= nrow(combos) && counters[1] <= maxes[1]) {
    reset <- reset0

    if (is.null(includebp) || all(includebp %in% counters)) {
      ok <- TRUE

      if (minvert_override != minvert) {
        # If any violations of minvert aren't in includebp, move on
        for (i in which(diff(y[counters]) < minvert)) {
          if (!all(counters[c(i, i + 1)] %in% includebp)) {
            ok <- FALSE
            break
          }
        }
      }

      if (ok) {
        combos[k, ] <- counters
        k <- k + 1
      }
    }

    counters[nbp] <- counters[nbp] + 1
    i <- nbp
    while (i > 1 && counters[i] > maxes[i]) {
      reset[i] <- TRUE
      counters[i - 1] <- counters[i - 1] + 1
      i <- i - 1
    }

    if (any(reset)) {
      for (i in which(reset)) {
        counters[i] <- counters[i - 1] + minvert_override
      }
    }
  }

  combos <- combos[seq_len(k - 1),, drop = FALSE]
  combos[] <- vert[combos]

  combos
}

#Bare bones, uses recursive solution from https://stackoverflow.com/a/76975508/6348551
#(fastest)
.combosR_quick <- function(vert, minvert, nbp) {

  y <- minvert:(length(vert) - minvert)
  vn <- y[length(y)]

  f0 <- function(v, k, m) {
    if (k == 1) return(matrix(v))

    d <- Recall(v, k - 1, m)

    u <- (m * (k - 1)):(vn - m)

    lst <- lapply(u, function(i) {
      p <- (i + m):vn
      dd <- d[d[, k - 1] == i, , drop = FALSE]
      cbind(dd[rep(1:nrow(dd), each = length(p)), , drop = FALSE], p)
    })

    unname(do.call("rbind", lst))
  }

  combos <- f0(y, nbp, minvert)

  colnames(combos) <- paste0("breakpoint", seq_len(nbp))

  combos <- combos[do.call("order", as.data.frame(combos)),, drop = FALSE]

  combos[] <- vert[combos]

  combos
}

.design_matrix <- function(Xvar, BPs, cont) {
  if (length(BPs) == 0) {
    return(cbind(1, Xvar))
  }

  bpmat <- matrix(BPs,
                  ncol = length(BPs),
                  nrow = length(Xvar),
                  byrow = TRUE)

  if (cont) {
    # Continuous fit
    cbind(1, Xvar, .pmax2(Xvar - bpmat, 0))
  }
  else {
    # Discontinuous fit
    m <- cbind(TRUE, Xvar > bpmat) & cbind(Xvar <= bpmat, TRUE)

    cbind(1 * m, Xvar * m)
  }
}
