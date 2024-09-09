library(pwr)

# Define parameters
n_levels <- 2 # Number of levels for each factor
n_subjects <- 52 # Estimated number of subjects
alpha <- 0.05
power <- 0.8

# Perform power analysis
f2 <- pwr.f2.test(u = (n_levels - 1)^2,  # degrees of freedom for numerator
                  v = (n_subjects - 1),  # degrees of freedom for denominator
                  sig.level = alpha,
                  power = power)$f2

# Convert f2 to partial eta squared
partial_eta_squared <- f2 / (1 + f2)
