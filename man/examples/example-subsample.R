data("alligator")

alligator_data <- process_measurements(alligator,
                                       pos = "Vertebra")

# Compute PCOs
alligator_PCO <- svdPCO(alligator_data)

# Plot vertebrae before subsampling
plotvertmap(alligator_PCO, dropNA = FALSE,
            text = TRUE)

# Subsample data after estimating PCOs; subsample down
# to 15 vertebrae
alligator_PCO_sample <- subsample(alligator_PCO,
                                  sample = 15)

# Plot vertebrae after subsampling; gray vertebrae
# have been dropped
plotvertmap(alligator_PCO_sample, dropNA = FALSE,
            text = TRUE)
