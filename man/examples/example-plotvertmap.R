data("alligator")

alligator_data <- process_measurements(alligator,
                                       pos = "Vertebra")

# Compute PCOs
alligator_PCO <- svdPCO(alligator_data)

# Plot vertebral map with specified breakpoints
plotvertmap(alligator_PCO,
            type = "percent",
            name = "Alligator",
            bps = c(8, 15, 19),
            text = TRUE)

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

# For each number of regions, identify best
# model based on minimizing RSS
bestresults <- modelselect(regionresults)

# Evaluate support for each model and rank models
supp <- modelsupport(bestresults)

# Plot vertebral map with breakpoints corresponding to
# best segmented regression model as determined by
# AICc
plotvertmap(alligator_PCO,
            type = "percent",
            name = "Alligator",
            modelsupport = supp,
            model = 1,
            criterion = "aic",
            text = TRUE)

# Plot vertebral map with breakpoints corresponding to
# best segmented regression model as determined by
# AICc, using centrum length to size vertebrae
plotvertmap(alligator_PCO,
            name = "Alligator",
            modelsupport = supp,
            model = 1,
            criterion = "aic",
            centraL = "CL",
            text = TRUE)

# Compute Akaike-weighted location and SD of optimal
# breakpoints using top 10% of models with 4 regions
bpvar <- calcBPvar(regionresults, noregions = 5,
                   pct = .1, criterion = "aic")

#Using weighted BPs and SDs from calcBPvar()
plotvertmap(alligator_PCO, name = "Dolphin",
            bpvar = bpvar,
            text = TRUE)
