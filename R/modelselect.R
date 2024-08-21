#' Select the best models
#'
#' `modelselect()` narrows down the search for the best model by identifying the best model for each number of regions as determined by its residual sums of squares (RSS).
#'
#' @param results a `regions_results` object; the output of a call to `calcregions()` or `addregions()`.
#' @param scores `numeric`; a vector corresponding to the indices of the PCOs the \eqn{R^2} of which will be used to determine the best model for each number of regions. If `NULL`, the default, all PCOs used included in the fitting will be used.
#'
#' @returns A `regions_modelselect` object, which contains information about the best models for each number of regions extracted from `results`.
#'
#' @seealso [modelsupport()] for computing statistics that describe the support of each model using information criteria; [modelperf()] for computing fit statistics for selected models.
#'
#' @example man/examples/example-modelselect.R
#'

#' @export
modelselect <- function(results, scores = NULL) {

  chk::chk_is(results, "regions_results")

  regiondata <- results$results

  noregions <- sort(unique(results$stats$Nregions))

  if (!is.null(scores)) {
    chk::chk_whole_numeric(scores)
    chk::chk_range(scores, c(1, sum(startsWith(names(regiondata), "RSS."))))

    keep.pcos <- paste0("RSS.", sort(scores))
    regiondata$sumRSS <- rowSums(regiondata[keep.pcos])

    regiondata[startsWith(names(regiondata), "RSS.") & !names(regiondata) %in% keep.pcos] <- NULL
  }

  models <- do.call("rbind", lapply(noregions, function(i) {
    allmodels <- regiondata[regiondata$regions == i,, drop = FALSE]	#select only models with correct region no
    allmodels[which.min(allmodels$sumRSS),, drop = FALSE]
  }))

  attr(models, "cont") <- attr(results, "cont")
  attr(models, "nvert") <- length(attr(results, "pos"))

  class(models) <- c("regions_modelselect", class(models))

  models
}

#' @exportS3Method print regions_modelselect
print.regions_modelselect <- function(x, digits = 3, ...) {
  x0 <- x
  class(x) <- setdiff(class(x), "regions_modelselect")
  for (i in which(startsWith(names(x), "breakpoint"))) {
    x[[i]] <- ifelse(is.na(x[[i]]), ".",
                     format(x[[i]], justify = "right"))
    names(x) <- sub("breakpoint", "BP ", names(x), fixed = TRUE)
  }

  for (i in which(grepl("RSS", names(x), fixed = TRUE))) {
    x[[i]] <- round(x[[i]], digits)
  }

  names(x)[1] <- "Regions"

  print(x, row.names = FALSE)

  invisible(x0)
}
