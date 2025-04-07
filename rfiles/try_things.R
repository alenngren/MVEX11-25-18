# Load necessary libraries
library(spatstat)
library(spatstat.core)
library(spatstat.geom)

# Generate a random spatial point pattern in a unit square
set.seed(123)
pp <- rpoispp(lambda = 100)  # Poisson process with intensity 100

# Compute Ripley's K-function with simulation envelopes
K_env <- envelope(pp, Lest, nsim = 999)  # 99 simulations for confidence envelope

# Plot the envelope
plot(K_env, main = "Envelope Plot for Ripley's K-Function")
