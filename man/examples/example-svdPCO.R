data("alligator")

alligator_data <- process_measurements(alligator,
                                       pos = "Vertebra")

# Compute PCOs
alligator_PCO <- svdPCO(alligator_data,
                        metric = "gower")

alligator_PCO

# Plot PCOs against vertebra index
plot(alligator_PCO, pco_y = 1:2)

# Plot PCOs against each other
plot(alligator_PCO, pco_y = 1, pco_x = 2)
