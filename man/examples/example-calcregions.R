data("alligator")

alligator_data <- process_measurements(alligator,
                                       pos = "Vertebra")

# Compute PCOs
alligator_PCO <- svdPCO(alligator_data)

# Fit segmented regression models for 1 to 5 regions
# using PCOs 1 to 4 and a continuous model with a
# non-exhaustive search
regionresults <- calcregions(alligator_PCO,
                             scores = 1:4,
                             noregions = 5,
                             minvert = 3,
                             cont = TRUE,
                             exhaus = FALSE,
                             verbose = FALSE)

regionresults

# View model fitting summary
summary(regionresults)

# Add additional regions to existing results,
# exhaustive search this time
regionresults <- addregions(regionresults,
                            noregions = 6:7,
                            exhaus = TRUE,
                            verbose = FALSE)

regionresults

summary(regionresults)

# Fit segmented regression models for 1 to 5 regions
# using PCOs 1 to 4 and a discontinuous model with a
# exhaustive search, excluding breakpoints at vertebrae
# 10 and 15
regionresults <- calcregions(alligator_PCO,
                             scores = 1:4,
                             noregions = 5,
                             minvert = 3,
                             cont = FALSE,
                             omitbp = c(10, 15),
                             verbose = FALSE)

regionresults

summary(regionresults)

# Compute the number of breakpoint combinations for given
# specification using `ncombos()`; if any number exceeds
# the value supplied to `ncombos_file_trigger`, results
# will temporary be stored in files before being read in and
# deleted.
ncombos(alligator_PCO,
        noregions = 1:8,
        minvert = 3)
