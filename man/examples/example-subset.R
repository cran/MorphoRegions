data("alligator")

alligator_data <- process_measurements(alligator,
                                       pos = "Vertebra")

# Compute PCOs
alligator_PCO <- svdPCO(alligator_data)

# Plot vertebrae before subsetting
plotvertmap(alligator_PCO, text = TRUE)

# Subset data after estimating PCOs
alligator_PCO_subset <- subset(alligator_PCO,
                               Vertebra < 20)

# Plot vertebrae after subsetting; vertebrae outside
# subset are dropped
plotvertmap(alligator_PCO_subset, text = TRUE)

# Subset data after estimating PCOs with `drop = FALSE`
alligator_PCO_subset <- subset(alligator_PCO,
                               Vertebra < 20,
                               drop = FALSE)

# Plot vertebrae after subsetting; vertebrae outside
# subset remain and are treated as missing
plotvertmap(alligator_PCO_subset, text = TRUE)
