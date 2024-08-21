# Process dataset; vertebra index in "Vertebra" column
data("alligator")

alligator_data <- process_measurements(alligator,
                                       pos = "Vertebra")

# Process multiple datasets; vertebra index in first column
data("porpoise")

porpoise_data <- process_measurements(list(porpoise1,
                                           porpoise2,
                                           porpoise3),
                                      pos = 1)
