data("alligator")

alligator_data <- process_measurements(alligator,
                                       pos = "Vertebra")

# Compute PCOs
alligator_PCO <- svdPCO(alligator_data)

# Evaluate model performance (R2) given supplied
# breakpoints for a continuous model
modelperf(alligator_PCO, scores = 1:3,
          bps = c(7, 15, 20), cont = TRUE)

plotsegreg(alligator_PCO, scores = 1:3,
           bps = c(7, 15, 20), cont = TRUE)
## See also `?calcmodel` for use with a single model

# Fit segmented regression models for 1 to 7 regions
# using PCOs 1 to 4 and a continuous model with a
# non-exhaustive search
regionresults <- calcregions(alligator_PCO,
                             scores = 1:4,
                             noregions = 7,
                             minvert = 3,
                             cont = TRUE,
                             exhaus = FALSE,
                             verbose = FALSE)

regionresults

# For each number of regions, identify best
# model based on minimizing RSS
bestresults <- modelselect(regionresults)

# Evaluate support for each model and rank
supp <- modelsupport(bestresults)

# Evaluate model performance (R2) for best model
# as chosen by BIC
modelperf(alligator_PCO, scores = 1:4,
          modelsupport = supp,
          criterion = "bic", model = 1)

# Plot that model for the first PCO score
plotsegreg(alligator_PCO, scores = 1:4,
           modelsupport = supp,
           criterion = "bic", model = 1)

## See `?simregions` for use with simulated data
