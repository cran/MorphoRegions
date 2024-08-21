data("alligator")

alligator_data <- process_measurements(alligator,
                                       pos = "Vertebra")

# Compute PCOs
alligator_PCO <- svdPCO(alligator_data)

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
bestresults

# Evaluate support for each model and rank models
supp <- modelsupport(bestresults)
supp

# 5 regions best based on AICc; 6 regions based on BIC
