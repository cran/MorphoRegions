data("alligator")

alligator_data <- process_measurements(alligator,
                                       pos = "Vertebra")

# Compute PCOs
alligator_PCO <- svdPCO(alligator_data)

# Compute PCO loadings
loadings <- PCOload(alligator_PCO, scores = 1:4)
loadings

# Plot loadings
plot(loadings)
