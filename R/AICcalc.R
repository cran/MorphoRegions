.AICcalc <- function(RSS, noPC, nvert, noregions, cont) {
  n <- noPC * nvert 	# No of variables used
  var <- RSS/n 	# Variance calculated ML way
  k <- { # Calculating the number of parameters (k) being fitted in the model
    if (cont) {
      # For a continuous fit, a slope is estimated for each region and each PCO,
      # only the intercept of the fist segment is estimated for each PCO since
      # other intercepts are forced by the previous slope and breakpoint
      # position, + a noregions-1 number of breakpoint is estimated
      (noregions * noPC) + noPC + noregions - 1
    }
    else {
      # For discontinuous fit, a slope and an intercept (2 params) are evaluated
      # for region and for each PCO, + a noregions-1 number of breakpoint is estimated
      (2 * noregions * noPC) + noregions - 1
    }
  }

  if (n < k + 2) {
    chk::err("the ratio of variables to parameters is too small. Reduce the number of regions or increase the number of variables")
  }

  AIC <- n * log(var) + (2 * k)
  corr <- (2 * k * (k + 1)) / (n - k - 1)	# Correct for number of parameters and small sample
  AICc <- AIC + corr 	# Calculate AICc
  BIC <- n * log(var) + log(n) * k

  c(AICc = AICc, BIC = BIC)
}
