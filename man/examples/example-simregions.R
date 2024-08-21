# Simulate 40 vertebra, 4 regions (3 breakpoints), 3 PCOs,
# true model R2 of .9, continuous
set.seed(11)
sim <- simregions(nvert = 30, nregions = 4, nvar = 3, r2 = .95,
                  minvert = 3, cont = TRUE)

sim

# Plot the true data-generating lines and breakpoints
plot(sim, scores = 1:3)

# Run segmented regression models on simulated data,
# up to 6 regions
simresults <- calcregions(sim, scores = 1:3, noregions = 6,
                          minvert = 3, cont = TRUE,
                          verbose = FALSE)

summary(simresults)

# Select best model for each number of regions
(simmodels <- modelselect(simresults))

# Evaluate support for each model and rank models
(simsupp <- modelsupport(simmodels))
# AICc supports 3-4 regions

# Evaluate model performance of best model
modelperf(sim, modelsupport = simsupp,
          criterion = "aic", model = 1)
# Second best model (3 regions) does quite well, too
modelperf(sim, modelsupport = simsupp,
          criterion = "aic", model = 2)

#Plot best model fit
plotsegreg(sim, scores = 1:3,
           modelsupport = simsupp,
           criterion = "aic", model = 1)

# Calculate variability of estimate breakpoints for
# 3-region model; high uncertainty for breakpoints
# 1 and 2. Note weighted value for breakpoint 2
# differs from that of best model
bpvar <- calcBPvar(simresults, noregions = 4,
                   pct = .05, criterion = "aic")
bpvar

# Plot estimated vertebral map with variability
plotvertmap(sim, modelsupport = simsupp, model = 1,
            criterion = "aic", text = TRUE)

# True map; pretty close
plotvertmap(sim, bps = c(3, 7, 24),
            text = TRUE)
