# skatta parametrar

# Ladda paketet
library(spatstat)

# Skapa ett punktmönster (t.ex. slumpmässigt inom en kvadrat)
pp <- rpoispp(lambda = 100)  # 100 punkter per enhet area

# Skatta intensiteten (lambda)
lambda_hat <- intensity(pp)
print(lambda_hat)


# Simulera en inhomogen Poisson-process
lambda_function <- function(x, y) { 100 * (x + y) }
pp_inhom <- rpoispp(lambda_function)

# Skatta intensiteten icke-parametriskt
intensity_map <- density(pp_inhom, sigma = 0.05)
plot(intensity_map)


# Anpassa en Poisson-modell (homogen)
ppm_model <- ppm(pp, ~1)
summary(ppm_model)

# Anpassa en inhomogen modell med covariater
ppm_inhom <- ppm(pp, ~x + y)
summary(ppm_inhom)

# Anpassa en Strauss-modell
ppm_gibbs <- ppm(pp, ~1, Strauss(r=0.05))
summary(ppm_gibbs)
