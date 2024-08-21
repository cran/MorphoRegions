data("alligator")

alligator_data <- process_measurements(alligator,
                                       pos = "Vertebra")

# Compute PCOs
alligator_PCO <- svdPCO(alligator_data)

# Calculate a single segmented regression model
# using first 2 PCOs and a discontinuous model
regionsmodel <- calcmodel(alligator_PCO,
                          scores = 1:3,
                          bps = c(8, 16, 21),
                          cont = FALSE)

regionsmodel

#Evaluate performance (R2) on that model
modelperf(regionsmodel)

#Plot model results:
plotsegreg(regionsmodel, scores = 1:3)
