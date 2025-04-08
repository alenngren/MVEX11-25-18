# load data from csv 
# install.packages("readr")
library(readr)
iterpath <- '/Users/alexander/Chalmers/MVEX11-25-18/dataframes/unif-iter10.csv'
iterframe <- read_csv(iterpath, col_names = c("x", "y"))
plot(iterframe)

filename <- "/Users/alexander/Chalmers/MVEX11-25-18/ERIKA.rds"
ERIKA <- readRDS(filename)
x <- ERIKA$VES13$small$x
y <- ERIKA$VES13$small$y
X <- ppp(x=x, y=y, window=square(40))



library(spatstat)
X_iter <-  ppp(x=iterframe$x, y=iterframe$y, window=square(40))
# 
clarkevans.test(X_iter, alternative="clustered")

E <- envelope(X, Kest, nsim=9999, verbose=TRUE)
E_iter <- envelope(X_iter, Kest, nsim=9)

library(fields)
jet_colors <- tim.colors(5)

plot(E, sqrt(./pi) -r  ~ r, main = 'L ~ med iterativ', legend = FALSE) 
lines(E_iter$r, sqrt(E_iter$obs/pi) -E_iter$r , col = jet_colors[1], lty = 2, lwd = 2)

legend("topleft", legend = c("data", "E_iter"), 
       col = c(1, jet_colors[1]), lty = c(1, 2))

legend("topleft", legend = c("data", 'E_iter'), 
       col = c(1, jet_colors[1]), lty = c(1, 2))

text(9, sqrt(max(E$obs)/pi) - 10.2, "data", col = "black", font = 1, cex = 1)
text(9, sqrt(max(E_iter$obs)/pi) - 10.2, "E_iter", col = jet_colors[1], font = 1, cex = 1)


E <- envelope(X, Lest, nsim=999, verbose=FALSE) # L function, same no difference 
plot(E-E$r, legend=FALSE)


# Plot with more steps in iterative 
iter0 <- '/Users/alexander/Chalmers/MVEX11-25-18/dataframes/unif-iter0-aut.csv'
iter5 <- '/Users/alexander/Chalmers/MVEX11-25-18/dataframes/unif-iter5-aut.csv'
iter10 <- '/Users/alexander/Chalmers/MVEX11-25-18/dataframes/unif-iter10-aut.csv'
iterframe0 <- read_csv(iter0, col_names = c("x", "y"))
iterframe5 <- read_csv(iter5, col_names = c("x", "y"))
iterframe10 <- read_csv(iter10, col_names = c("x", "y"))

plot(iterframe0, win=square(40))
plot(iterframe5, win=square(40))
plot(iterframe10, win=square(40))

X_iter0 <-  ppp(x=iterframe0$x, y=iterframe0$y, window=square(40))
X_iter5 <-  ppp(x=iterframe5$x, y=iterframe5$y, window=square(40))
X_iter10 <-  ppp(x=iterframe10$x, y=iterframe10$y, window=square(40))
# 
clarkevans.test(X_iter0, alternative="clustered")
clarkevans.test(X_iter5, alternative="clustered")
clarkevans.test(X_iter10, alternative="clustered")

E_0 <- envelope(X_iter0, Kest, nsim=999, verbose=TRUE)
E_5 <- envelope(X_iter5, Kest, nsim=999, verbose=TRUE)
E_10 <- envelope(X_iter10, Kest, nsim=999, verbose=TRUE)

E_iter <- envelope(X_iter, Kest, nsim=9)


plot(E, sqrt(./pi) -r  ~ r, main = 'L ~ iterativ steg 1,6,11', legend = FALSE) 
lines(E_0$r, sqrt(E_0$obs/pi) -E_0$r , col = jet_colors[3], lty = 1, lwd = 2) 
lines(E_5$r, sqrt(E_5$obs/pi) -E_5$r , col = jet_colors[2], lty = 2, lwd = 2) 
lines(E_10$r, sqrt(E_10$obs/pi) -E_10$r , col = jet_colors[1], lty = 2, lwd = 2) 
lines(E_half$r, sqrt(E_half$obs/pi) -E_half$r, col = jet_colors[4], lty = 4, lwd = 2)

text(9, sqrt(max(E$obs)/pi) - 10.1, "data", col = "black", font = 1, cex = 1)
text(9, sqrt(max(E_0$obs)/pi) - 9.8, "E_1", col = jet_colors[3], font = 1, cex = 1)
text(9, sqrt(max(E_5$obs)/pi) - 10.4, "E_6", col = jet_colors[2], font = 1, cex = 1)
text(9, sqrt(max(E_10$obs)/pi) - 10.2, "E_11", col = jet_colors[1], font = 1, cex = 1)

text(9, sqrt(max(E_half$obs)/pi) - 10.2, "E_half", col = jet_colors[4], font = 1, cex = 1)


E <- envelope(X, Kest, nsim=99)

plot(E, sqrt(./pi) -r  ~ r, main = 'L ~ iterativ steg 1,6,11', legend = FALSE) 

# Summary statistics 









# doesnt quite work how I want it 
env <- envelope(X, fun = Kest, nsim = 99, rank = 1, global = FALSE)

# Plot envelope
plot(env, main = "Envelope of CSR around K-function")


# Simulate CSR envelopes (e.g., 9 simulations)
nsim <- 9
csr_sims <- replicate(nsim, rpoispp(lambda = intensity(X), win = Window(X)), simplify = FALSE)

# Plot observed + CSR simulations
par(mfrow = c(2, 5), mar = c(1, 1, 2, 1))  # layout: 2 rows, 5 columns

# Plot observed pattern
plot(X, main = "Observed", pch = 20, cols = "red")

# Plot CSR simulations
for (i in 1:nsim) {
  plot(csr_sims[[i]], main = paste("CSR", i), pch = 20)
}



# Load required package
library(spatstat)

# Assume X is your observed spatial point pattern (ppp object)
# Here's an example (replace this with your actual pattern)
set.seed(123)
X <- rpoispp(lambda = 100, win = owin(c(0, 1), c(0, 1)))  # Example pattern

# Number of CSR simulations to generate
nsim <- 9

# Simulate CSR patterns (same intensity and window as X)
csr_sims <- replicate(nsim, rpoispp(lambda = intensity(X), win = Window(X)), simplify = FALSE)

# Plot observed pattern and CSR simulations
par(mfrow = c(2, 5), mar = c(1, 1, 2, 1))  # layout: 2x5 grid

# Plot the observed pattern
plot(X, main = "Observed", pch = 20, cols = "red")

# Plot each CSR simulation
for (i in 1:nsim) {
  plot(csr_sims[[i]], main = paste("CSR", i), pch = 20)
}




# Load required package
library(spatstat)

# Assume X is your observed point pattern (ppp object)
# Here's an example (replace this with your own pattern)
set.seed(123)
# X <- rpoispp(100)  # or use your own: X <- your_ppp_object

# Generate envelope using Kest (CSR simulation envelope)
env <- envelope(X, fun = Kest, nsim = 99, rank = 1, global = FALSE)

# Plot the envelope
par(mfrow = c(1, 1))  # la
env <- envelope(X, fun = Kest, nsim = 99, global = TRUE)
plot(env, main = "CSR Envelope Around K-function", legendargs = list(cex = 0.8))





# Load required package
library(spatstat)

# Your observed point pattern X (replace with your own ppp object)
set.seed(123)
X <- rpoispp(100)  # or your own: X <- your_ppp_object

# Compute observed K function
K_obs <- Kest(X)

# Generate envelope from CSR, but exclude observed pattern from envelope band
env <- envelope(X, fun = Kest, nsim = 99, savefuns = TRUE, simulate = expression(rpoispp(lambda = intensity(X), win = Window(X))), 
                global = FALSE, include = FALSE)  # include = FALSE excludes observed from envelope stats

en
# Plot the CSR envelope (as a shaded area)
plot(env, . - theo ~ r, shade = c("hi", "lo"), main = "Observed K with CSR Envelope")

# Overlay the observed K function
lines(K_obs$r, K_obs$iso - K_obs$theo, col = "red", lwd = 2)
legend("topright", legend = c("Observed", "CSR Envelope"), col = c("red", "grey"), lwd = c(2, NA), fill = c(NA, "lightgray"))



# gpt fix 
# Load spatstat
library(spatstat)

# Create observed point pattern X (replace with your own)
set.seed(123)
X <- rpoispp(100)

# Compute observed K function (with 'theo')
K_obs <- Kest(X)

# Generate CSR envelope without including the observed pattern in the stats
env <- envelope(X, 
                fun = Kest, 
                nsim = 99, 
                simulate = expression(rpoispp(lambda = intensity(X), win = Window(X))), 
                savefuns = TRUE, 
                savepatterns = TRUE,
                global = FALSE,
                include = FALSE)

# Just to be safe: Add theo back manually if it's missing
if (!"theo" %in% names(env)) {
  env$theo <- K_obs$theo
}
# theo
# Plot envelope with shading
plot(env, . - env$obs ~ r, shade = c("hi", "lo"), main = "Observed K with CSR Envelope")

# Overlay observed K-function minus theo
lines(K_obs$r, K_obs$iso - K_obs$theo, col = "red", lwd = 2)

legend("topright", legend = c("Observed", "CSR Envelope"), 
       col = c("red", "lightgray"), lwd = c(2, NA), pch = c(NA, 15), pt.cex = 2)


# Third try 

library(spatstat)

# Observed point pattern
set.seed(123)
X <- rpoispp(100)

# Compute observed K-function
K_obs <- Kest(X)

# Simulate CSR envelopes
env <- envelope(X,
                fun = Kest,
                nsim = 99,
                simulate = expression(rpoispp(lambda = intensity(X), win = Window(X))),
                savefuns = TRUE,
                savepatterns = TRUE,
                global = FALSE,
                include = FALSE)

# Extract simulation results
sim_curves <- as.data.frame(attr(env, "simfuns"))

# simfuns includes columns: r, sim1, sim2, ..., simN
r_vals <- sim_curves$r
sim_matrix <- as.matrix(sim_curves[, -1])  # Remove r column

# Get pointwise min and max across simulations
env_lo <- apply(sim_matrix, 1, min)
env_hi <- apply(sim_matrix, 1, max)

# Plot observed K with envelope around it
plot(K_obs$r, K_obs$iso, type = "l", col = "red", lwd = 2,
     ylim = range(c(env_lo, env_hi, K_obs$iso)),
     main = "CSR Envelope Around Observed K-function",
     xlab = "r", ylab = "K(r)")

# Add shaded CSR envelope around the observed line
polygon(c(r_vals, rev(r_vals)),
        c(env_hi, rev(env_lo)),
        col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)

# Redraw observed line on top
lines(K_obs$r, K_obs$iso, col = "red", lwd = 2)

legend("topleft", legend = c("Observed", "CSR Envelope"),
       col = c("red", rgb(0.8, 0.8, 0.8)), lwd = c(2, NA), pch = c(NA, 15), pt.cex = 2)



# 4 th try
X <- ppp(x=x, y=y, window=square(40))
library(spatstat)

# Create your point pattern
set.seed(123)
X <- rpoispp(100)

# Compute observed K-function
K_obs <- Kest(X)

# Simulate envelope (excluding observed from stats)
env <- envelope(X,
                fun = Kest,
                nsim = 99,
                simulate = expression(rpoispp(lambda = intensity(X), win = Window(X))),
                savefuns = TRUE,
                include = FALSE)

# Extract simulated K functions
sim_funs <- as.data.frame(attr(env, "simfuns"))
r_vals <- sim_funs$r
sim_vals <- sim_funs[, -1]  # remove r column

# Subtract observed K from each simulation to get deviation from observed
obs_vals <- K_obs$iso
sim_deviation <- sweep(sim_vals, 1, obs_vals, FUN = "-")

# Compute min and max deviations from the observed line
lo_diff <- apply(sim_deviation, 1, min)
hi_diff <- apply(sim_deviation, 1, max)

# Plot observed K-function > START HERE! WORKS 
plot(K_obs$r, K_obs$iso, type = "l", col = "red", lwd = 2,
     #ylim = range(c(K_obs$iso + hi_diff, K_obs$iso + lo_diff)),
     main = "CSR Envelope Around Observed K-function",
     xlab = "r", ylab = "K(r)")

enve <- envelope(X, Kest, nsim=99)

lines(enve$r, enve$hi-enve$theo+enve$obs) # add the line we want 
lines(enve$r, enve$lo-enve$theo+enve$obs)

# Add shaded CSR deviation envelope around the observed line
polygon(c(r_vals, rev(r_vals)),
        c(K_obs$iso + hi_diff, rev(K_obs$iso + lo_diff)),
        col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)

# this works as envelope
polygon(c(r_vals, rev(r_vals)),
        c( enve$hi-enve$theo+enve$obs, rev( enve$lo-enve$theo+enve$obs )),
        col = rgb(0.8, 0.8, 0.8, 0.9), border = NA)

# Re-draw observed line on top
lines(K_obs$r, K_obs$iso, col = "red", lwd = 1)

legend("topleft", legend = c("Observed", "CSR Envelope"),
       col = c("red", rgb(0.8, 0.8, 0.8)), lwd = c(2, NA), pch = c(NA, 15), pt.cex = 2)

# --------- THIS IS PRETTY PLOT ------------ # 
# NOW DO IT FOR L FUNCTION
# Plot observed K-function > START HERE! WORKS 
plot(K_obs$r, sqrt(K_obs$iso/pi)-K_obs$r, type = "l", col = "red", lwd = 1,
     #ylim = range(c(K_obs$iso + hi_diff, K_obs$iso + lo_diff)),
     main = "CSR Envelope Around Observed K-function",
     xlab = "r", ylab = "K(r)")

enve <- envelope(X, Lest, nsim=99)

lines(enve$r, (enve$hi-enve$theo+enve$obs)-enve$r) # add the line we want 
lines(enve$r, (enve$lo-enve$theo+enve$obs)-enve$r)

# Add shaded CSR deviation envelope around the observed line
polygon(c(r_vals, rev(r_vals)),
        c(K_obs$iso + hi_diff, rev(K_obs$iso + lo_diff)),
        col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)

# this works as envelope
polygon(c(r_vals, rev(r_vals)),
        c( (enve$hi-enve$theo+enve$obs)-enve$r, rev( (enve$lo-enve$theo+enve$obs)-enve$r )),
        col = rgb(0.8, 0.8, 0.8, 0.9), border = NA)

# Re-draw observed line on top
lines(K_obs$r, sqrt(K_obs$iso/pi)-K_obs$r, col = "red", lwd = 1)

# Add other data models 
lines(E_0$r, sqrt(E_0$obs/pi) -E_0$r , col = jet_colors[3], lty = 1, lwd = 2) 
lines(E_5$r, sqrt(E_5$obs/pi) -E_5$r , col = jet_colors[2], lty = 2, lwd = 2) 
lines(E_10$r, sqrt(E_10$obs/pi) -E_10$r , col = jet_colors[1], lty = 2, lwd = 2) 
lines(E_half$r, sqrt(E_half$obs/pi) -E_half$r, col = jet_colors[4], lty = 4, lwd = 2)

text(9, sqrt(max(E$obs)/pi) - 10.1, "data", col = "black", font = 1, cex = 1)
text(9, sqrt(max(E_0$obs)/pi) - 9.8, "E_1", col = jet_colors[3], font = 1, cex = 1)
text(9, sqrt(max(E_5$obs)/pi) - 10.4, "E_6", col = jet_colors[2], font = 1, cex = 1)
text(9, sqrt(max(E_10$obs)/pi) - 10.2, "E_11", col = jet_colors[1], font = 1, cex = 1)

text(9, sqrt(max(E_half$obs)/pi) - 10.2, "E_half", col = jet_colors[4], font = 1, cex = 1)

legend("topleft", legend = c("Observed", "CSR Envelope"),
       col = c("red", rgb(0.8, 0.8, 0.8)), lwd = c(2, NA), pch = c(NA, 15), pt.cex = 2)





# 5th try
X <- ppp(x=x, y=y, window=square(40))
library(spatstat)

# Your observed point pattern
# Replace with your real pattern!
# X <- your_ppp_object
set.seed(123)
X <- rpoispp(100)

# Compute observed K-function
K_obs <- Kest(X, correction = "iso")

# Simulate CSR and compute K-function for each
nsim <- 99
sim_K <- simulate(expression(rpoispp(lambda = intensity(X), win = Window(X))), nsim = nsim, envir = environment())

# Compute K-function for each simulation
K_list <- lapply(sim_K, function(sim) Kest(sim, correction = "iso"))

# Extract the r values (they are all the same)
r_vals <- K_obs$r

# Build matrix of iso K(r) values from all simulations
K_matrix <- sapply(K_list, function(K) K$iso)

# Compute pointwise quantiles
lower <- apply(K_matrix, 1, quantile, probs = 0.025)
upper <- apply(K_matrix, 1, quantile, probs = 0.975)

# Plot observed K-function
plot(K_obs$r, K_obs$iso, type = "l", col = "red", lwd = 2,
     ylim = range(c(lower, upper, K_obs$iso)),
     main = "CSR Envelope (Simulated) Around Observed K-function",
     xlab = "r", ylab = "K(r)")

# Add shaded CSR envelope
polygon(c(r_vals, rev(r_vals)),
        c(upper, rev(lower)),
        col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)

# Re-draw observed line
lines(K_obs$r, K_obs$iso, col = "red", lwd = 2)

legend("topleft", legend = c("Observed", "CSR envelope (2.5%–97.5%)"),
       col = c("red", rgb(0.8, 0.8, 0.8)), lwd = c(2, NA), pch = c(NA, 15), pt.cex = 2)


# 6th try

library(spatstat)

# Replace with your actual point pattern
# X <- your_ppp_object
set.seed(123)
X <- rpoispp(100)

# Observed K-function
K_obs <- Kest(X, correction = "iso")

# Simulation settings
nsim <- 99
lambda_hat <- intensity(X)
win <- Window(X)

# Simulate nsim CSR patterns and compute K-functions
K_sims <- lapply(1:nsim, function(i) {
  sim <- rpoispp(lambda = lambda_hat, win = win)
  Kest(sim, correction = "iso")
})

# Matrix of simulated iso K(r) values
K_matrix <- sapply(K_sims, function(K) K$iso)

# r values (same for all)
r_vals <- K_obs$r

# Compute 2.5% and 97.5% quantiles across simulations at each r
lo <- apply(K_matrix, 1, quantile, probs = 0.025)
hi <- apply(K_matrix, 1, quantile, probs = 0.975)

# Plot observed K with simulated envelope around it
plot(r_vals, K_obs$iso, type = "l", col = "red", lwd = 2,
     ylim = range(c(lo, hi, K_obs$iso)),
     main = "CSR Envelope Around Observed K-function",
     xlab = "r", ylab = "K(r)")

# Add shaded envelope around observed data
polygon(c(r_vals, rev(r_vals)),
        c(hi, rev(lo)),
        col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)

# Redraw observed K
lines(r_vals, K_obs$iso, col = "red", lwd = 2)

legend("topleft", legend = c("Observed", "CSR Envelope"),
       col = c("red", rgb(0.8, 0.8, 0.8)), lwd = c(2, NA), pch = c(NA, 15), pt.cex = 2)

# 7th try 

# install.packages("ggplot2")
# install.packages("dplyr")


library(spatstat)
library(ggplot2)
library(dplyr)

# Example point pattern (replace with your real one)
set.seed(123)
X <- rpoispp(100)

# Observed K-function
K_obs <- Kest(X, correction = "iso")

# Simulate CSR patterns and compute Kest
nsim <- 99
lambda_hat <- intensity(X)
win <- Window(X)

K_sims <- lapply(1:nsim, function(i) {
  sim <- rpoispp(lambda_hat, win = win)
  Kest(sim, correction = "iso")
})

# Stack all simulated iso K(r) into a data frame
sim_df <- do.call(cbind, lapply(K_sims, function(k) k$iso))
r_vals <- K_obs$r

# Compute envelope (2.5% and 97.5% quantiles) at each r
envelope_df <- data.frame(
  r = r_vals,
  lower = apply(sim_df, 1, quantile, probs = 0.025),
  upper = apply(sim_df, 1, quantile, probs = 0.975),
  observed = K_obs$iso
)

# Plot with ggplot2
ggplot(envelope_df, aes(x = r)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightgray") +
  geom_line(aes(y = observed), color = "red", size = 1) +
  labs(
    title = "CSR Envelope Around Observed K-function",
    x = "r",
    y = "K(r)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# 8th try

library(spatstat)
library(ggplot2)
library(dplyr)

# Replace with your actual pattern
X <- ppp(x=x, y=y, window=square(40))
set.seed(123)
X <- rpoispp(100)

# Observed K-function
K_obs <- Kest(X, correction = "iso")
r_vals <- K_obs$r
K_obs_vals <- K_obs$iso

# Simulate CSR and compute deviations from observed
nsim <- 99
lambda_hat <- intensity(X)
win <- Window(X)

# Compute simulated K-functions
K_sims <- lapply(1:nsim, function(i) {
  sim <- rpoispp(lambda = lambda_hat, win = win)
  Kest(sim, correction = "iso")$iso
})

# Stack simulations into matrix
sim_matrix <- do.call(cbind, K_sims)

# Compute deviations from observed at each r
deviation_matrix <- sweep(sim_matrix, 1, K_obs_vals, FUN = "-")

# Compute 2.5% and 97.5% quantiles of deviation
lo_deviation <- apply(deviation_matrix, 1, quantile, probs = 0.025)
hi_deviation <- apply(deviation_matrix, 1, quantile, probs = 0.975)

# Envelope centered on observed K
envelope_df <- data.frame(
  r = r_vals,
  observed = K_obs_vals,
  lower = K_obs_vals + lo_deviation,
  upper = K_obs_vals + hi_deviation
)

# Plot using ggplot2
ggplot(envelope_df, aes(x = r)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightgray") +
  geom_line(aes(y = observed), color = "red", size = 1.2) +
  labs(
    title = "CSR Envelope Centered Around Observed K-function",
    x = "r",
    y = "K(r)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# 9 th try 

library(spatstat)
library(ggplot2)

# Your actual point pattern here
set.seed(123)
X <- rpoispp(100)

# Observed K-function
K_obs <- Kest(X, correction = "iso")
r_vals <- K_obs$r
K_obs_vals <- K_obs$iso

# Simulate CSR patterns and compute K-functions
nsim <- 99
lambda_hat <- intensity(X)
win <- Window(X)

K_sims <- lapply(1:nsim, function(i) {
  sim <- rpoispp(lambda_hat, win)
  Kest(sim, correction = "iso")$iso
})

# Convert to matrix (rows = r, cols = simulations)
K_matrix <- do.call(cbind, K_sims)

# Compute deviation from observed at each r
# BUT: store only the range of variation from simulations
residuals <- sweep(K_matrix, 1, K_obs_vals, "-")

# Compute envelope of residuals (quantiles)
res_lo <- apply(residuals, 1, quantile, probs = 0.025)
res_hi <- apply(residuals, 1, quantile, probs = 0.975)

# Construct envelope around the observed K

lines(enve$r, enve$hi-enve$theo+enve$obs)
lines(enve$r, enve$lo-enve$theo+enve$obs)

enve <- envelope(X, Kest, nsim=99)

envelope_df <- data.frame( # original 
  r = r_vals,
  observed = K_obs_vals,
  lower = K_obs_vals + res_lo,
  upper = K_obs_vals + res_hi
)

envelope_df <- data.frame(
  r = r_vals,
  observed = K_obs_vals,
  lower = enve$lo-enve$theo+enve$obs-10, # change value
  upper = enve$hi-enve$theo+enve$obs+10
)

# GGPlot
ggplot(envelope_df, aes(x = r)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightgray") +
  geom_line(aes(y = observed), color = "red", size = 1) +
  labs(
    title = "CSR Envelope HUGGING Observed K(r)",
    x = "r",
    y = "K(r)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.position = "none"
  )

# 10 th try 

library(spatstat)
library(ggplot2)

# Use your own point pattern here
set.seed(123)
X <- rpoispp(100)

# Observed K
K_obs <- Kest(X, correction = "iso")
r_vals <- K_obs$r
K_obs_vals <- K_obs$iso

# Simulate CSR and extract K-functions
nsim <- 99
win <- Window(X)
lambda_hat <- intensity(X)

# Simulate and extract only K(r) at the same r
K_sims <- lapply(1:nsim, function(i) {
  sim <- rpoispp(lambda = lambda_hat, win = win)
  Kest(sim, correction = "iso")$iso
})

# Matrix: each column = simulation, rows = r
K_sim_matrix <- do.call(cbind, K_sims)

# Now — calculate deviations of CSR from **itself**
# Not subtracting observed!
# Compute how much CSR varies at each r
csr_sd <- apply(K_sim_matrix, 1, sd)
csr_q_lo <- apply(K_sim_matrix, 1, quantile, probs = 0.025)
csr_q_hi <- apply(K_sim_matrix, 1, quantile, probs = 0.975)
csr_iqr <- data.frame(lo = csr_q_lo, hi = csr_q_hi)

# Now — use K_obs as the CENTER, but apply CSR variation to build band
envelope_df <- data.frame(
  r = r_vals,
  observed = K_obs_vals,
  lower = K_obs_vals - (K_obs_vals - csr_q_lo),
  upper = K_obs_vals + (csr_q_hi - K_obs_vals)
)

# Plot: CSR envelope HUGGING observed data at every r
ggplot(envelope_df, aes(x = r)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "gray80") +
  geom_line(aes(y = observed), color = "red", size = 1.2) +
  labs(
    title = "CSR Envelope Stuck to Observed K(r)",
    x = "r",
    y = "K(r)"
  ) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold")
  )