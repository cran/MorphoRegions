data("alligator")

alligator_data <- process_measurements(alligator,
                                       pos = "Vertebra")

# Compute PCOs
alligator_PCO <- svdPCO(alligator_data)

# Fit segmented regression models for 1 to 7 regions
# using PCOs 1 to 4 and a continuous model with a
# exhaustive search
regionresults <- calcregions(alligator_PCO,
                             scores = 1:4,
                             noregions = 7,
                             minvert = 3,
                             cont = TRUE,
                             exhaus = TRUE,
                             verbose = FALSE)

# Compute Akaike-weighted location and SD of optimal
# breakpoints using top 10% of models with 4 regions
calcBPvar(regionresults, noregions = 4,
          pct = .1, criterion = "aic")
