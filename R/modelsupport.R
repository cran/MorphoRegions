#' Evaluate model support
#'
#' `modelsupport()` computes measures of the relative support of each of the best models identified by [modelselect()] to facilitate selecting the optimal number and position of regions. These measures are in the form of information criteria (AICc and BIC).
#'
#' @param models a `regions_modelselect` object; the output of a call to [modelselect()].
#'
#' @returns A `regions_modelsupport` object, which contains the best model for each number of regions as determined by the AICc and BIC. The computed statistics are `AICc`/`BIC`--the value of the information criterion (IC) for each model, `deltaAIC`/`deltaBIC`--the difference between the IC for the corresponding model and that of the model with the lowest IC value, `model_lik`--the likelihood ratio of the model against the model with the lowest IC value, and `Ak_weight`/`BIC_weight`--the Akaike weights for each model used to compute the region score. The region score is a weighted average of the numbers of regions, weighted by the Akaike weights to represent the variability around the optimal number of regions.
#'
#' @seealso [modelselect()], [calcregions()], [calcBPvar()], [modelperf()], [plotsegreg()]
#'
#' @example man/examples/example-modelselect.R
#'

#' @export
modelsupport <- function(models) {

  chk::chk_is(models, "regions_modelselect")

  cont <- attr(models, "cont")
  nvert <- attr(models, "nvert")

  nPC <- sum(startsWith(names(models), "RSS."))

  models <- models[!startsWith(names(models), "RSS.")]

  AICc <- numeric(nrow(models))
  BIC <- numeric(nrow(models))

  for (i in seq_len(nrow(models))) {
    probs <- .AICcalc(models$sumRSS[i], nPC, nvert,
                     models$regions[i], cont) #Calculate AIC score

    AICc[i] <- probs["AICc"]
    BIC[i] <- probs["BIC"]
  }

  # AICc:
  deltaAIC <- AICc - min(AICc)	#Calculate AIC difference
  model_lik <- exp(-deltaAIC / 2)	#Likelihood of the model

  Ak_weight <- model_lik/sum(model_lik)	#Akaike Weights

  AIC_models <- cbind(models, AICc, deltaAIC, model_lik, Ak_weight)
  AIC_models <- AIC_models[order(AIC_models$AICc),, drop = FALSE] #Sort so best at top

  Regions_score <- sum(AIC_models$regions * AIC_models$Ak_weight)

  # BIC:
  deltaBIC <- BIC - min(BIC)	#Calculate AIC difference
  model_lik <- exp(-deltaBIC / 2)	#Likelihood of the model

  BIC_weight <- model_lik/sum(model_lik)	#BIC Weights

  BIC_models <- cbind(models, BIC, deltaBIC, model_lik, BIC_weight)
  BIC_models <- BIC_models[order(BIC_models$BIC),, drop = FALSE] #Sort so best at top

  Regions_scoreB <- sum(BIC_models$regions * BIC_models$BIC_weight)

  out <- list(Model_support = AIC_models,
              Region_score = Regions_score,
              Model_support_BIC = BIC_models,
              Region_score_BIC = Regions_scoreB)

  attr(out, "cont") <- cont

  class(out) <- "regions_modelsupport"

  out
}

#' @exportS3Method print regions_modelsupport
print.regions_modelsupport <- function(x, digits = 3, ...) {
  x0 <- x
  for (m in c("Model_support", "Model_support_BIC")) {
    for (i in which(startsWith(names(x[[m]]), "breakpoint"))) {
      x[[m]][[i]] <- ifelse(is.na(x[[m]][[i]]), ".",
                       format(x[[m]][[i]], justify = "right"))
      names(x[[m]]) <- sub("breakpoint", "BP ", names(x[[m]]), fixed = TRUE)
    }

    for (i in setdiff(names(x[[m]]),
                      c("regions", names(x[[m]])[startsWith(names(x[[m]]), "BP ")]))) {
      x[[m]][[i]] <- round(x[[m]][[i]], digits)
    }

    names(x[[m]])[1] <- "Regions"
  }

  cat("- Model support (AICc)\n")
  print.data.frame(x[["Model_support"]], row.names = FALSE)
  cat("Region score:", round(x[["Region_score"]], 2), "\n\n")

  cat("- Model support (BIC)\n")
  print.data.frame(x[["Model_support_BIC"]], row.names = FALSE)
  cat("Region score:", round(x[["Region_score_BIC"]], 2), "\n")

  invisible(x0)
}
