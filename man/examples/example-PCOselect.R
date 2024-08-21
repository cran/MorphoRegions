data("alligator")

alligator_data <- process_measurements(alligator,
                                       pos = "Vertebra")

# Compute PCOs
alligator_PCO <- svdPCO(alligator_data)

# Select which PCOs to use
## Manually (first 4 PCOs)
(PCOs <- PCOselect(alligator_PCO, "manual", scores = 4))

## Using variance cutoff: PCOs that explain 5% or more
## of total PCO variance
(PCOs <- PCOselect(alligator_PCO, "variance", cutoff = .05))

## Using bootstrapping with 50 reps (more reps should
## be used in practice; default is fine)
(PCOs <- PCOselect(alligator_PCO, "boot", nreps = 50))

plot(PCOs) #plot true eigenvalues against bootstrapped

## Using PCOs that optimize region score:
regionresults <- calcregions(alligator_PCO, scores = 1:21, noregions = 7,
                             minvert = 3, cont = TRUE, exhaus = TRUE,
                             verbose = FALSE)

(PCOs <- PCOselect(alligator_PCO, "max",
                   results = regionresults,
                   criterion = "bic"))

plot(PCOs)

summary(PCOs)
