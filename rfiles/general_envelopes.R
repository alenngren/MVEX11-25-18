# Envelopes for other functions! and make it dynamic ---- 
setwd("/Users/alexander/Chalmers/MVEX11-25-18/fig/conf")
csv_dir <- "/Users/alexander/Chalmers/MVEX11-25-18/dataframes/LGCP_punkter_2"
# csv_files <- list.files(path = csv_dir, pattern = "\\.csv$", full.names = TRUE)
csv_files <- list.files(path = csv_dir, pattern = "\\.csv$", full.names = TRUE)
rx <- Lest(X, correction="iso")$r # ves13 r values 
# rx
# Initialize a list to hold the Lr values
L_iso_values <- vector("list", length(csv_files))

# Loop through each file
for (i in seq_along(csv_files)) {
  df <- read.csv(csv_files[i], header = FALSE, skip = 1)
  colnames(df) <- c("x", "y")
  
  # Set window (adjust as needed)
  pp <- ppp(x = df$x, y = df$y, window = square(40))
  
  # Compute L-function
  L <- Lest(pp, correction = "iso")$iso - rx
  
  # Store only the 'iso' values (as a row)
  L_iso_values[[i]] <- L
}

# Convert list of vectors to data frame
L_matrix <- do.call(rbind, L_iso_values)
L_matrix <- t(L_matrix)
L_df35 <- as.data.frame(L_matrix)
# rownames(L_df) <- paste0("rep_", seq_along(csv_files))
# View(L_df35) # check if it works 


# r <- seq(0, 10, length.out = 100) # not necessarY
# Beräkna 2.5:e och 97.5:e percentilen för varje punkt
lower <- apply(L_df35, 1, quantile, probs = 0.025)
upper <- apply(L_df35, 1, quantile, probs = 0.975)
mean_sim <- apply(L_df35, 1, mean)

l_X <- Lest(X, correction = "iso") # ves13 
# par(mfrow=c(1,1)) # standard plot 

{ # new new mats plot 
  # png("conf-hierarkisk1-col.png", width=694, height= 482)
  plot(rx, l_X$iso-rx, type = "l", col = "black", lwd = 2, lty=1,
       ylim = range(0, 8),
       ylab = expression(L(r)-r), xlab = expression(italic(r)),
       cex.main = 2, cex.lab = 1.2, cex.axis = 1.5,
       main = "Konfidenshölje för Hierarkisk Poissonprocess")
  
  coli <- jet_colors[10]
  rgb_vals <- col2rgb(coli)
  polygon(c(rx, rev(rx)), c(upper, rev(lower)),
          col = rgb(rgb_vals[1]/255,rgb_vals[2]/255,rgb_vals[3]/255, 0.1), border = NA)
  #col = rgb(0.8, 0.8, 1, 0.5), border = NA)
  
  # Lägg till medelvärdeskurvan ovanpå igen (för tydlighet)
  lines(rx, mean_sim, col = jet_colors[10], lwd = 2, lty=2)
  # lines(rx, Lest(X, Lest, correction='isotropic')$iso - E$r, lwd=2)
  lines(rx, l_X$theo-rx, lwd=2, lty=2, col='red')
  
  legend("topleft", legend = c('data', 'Hierarkisk','Poisson'),
         col = c('black',jet_colors[10], 'red'), lty = c(1,2,2), lwd=c(2,1,1)) 
  # dev.off()
}

# make it with gpt ----

generate_envelope_plot <- function(csv_dir, statfun = "Lest", X, 
                                   window_size = 40, correction = "iso") {
  csv_files <- list.files(path = csv_dir, pattern = "\\.csv$", full.names = TRUE)
  
  if (length(csv_files) == 0) stop("No CSV files found in directory.")
  
  # Evaluate the correct function
  fun <- match.fun(statfun)
  
  # Reference r values from observed data
  rx <- fun(X, correction = correction)$r
  obs_stat <- fun(X, correction = correction)[[correction]] - rx
  
  # Initialize list for simulated values
  stat_list <- vector("list", length(csv_files))
  
  for (i in seq_along(csv_files)) {
    df <- read.csv(csv_files[i], header = FALSE, skip = 1)
    colnames(df) <- c("x", "y")
    pp <- ppp(x = df$x, y = df$y, window = square(window_size))
    sim_stat <- fun(pp, correction = correction)[[correction]] - rx
    stat_list[[i]] <- sim_stat
  }
  
  # Convert to matrix and data frame
  stat_matrix <- do.call(rbind, stat_list)
  stat_matrix <- t(stat_matrix)
  stat_df <- as.data.frame(stat_matrix)
  
  # Calculate envelope
  lower <- apply(stat_df, 1, quantile, probs = 0.025)
  upper <- apply(stat_df, 1, quantile, probs = 0.975)
  mean_sim <- rowMeans(stat_df)
  theo_line <- fun(X, correction = correction)$theo - rx
  
  # Plot
  jet_colors <- colorRampPalette(brewer.pal(9, "Set1"))(10)
  plot(rx, obs_stat, type = "l", col = "black", lwd = 2, lty = 1,
       ylim = range(c(lower, upper, obs_stat)), 
       ylab = expression(S(r)-r), xlab = expression(italic(r)),
       main = paste("Envelope for", statfun, "function"),
       cex.lab = 1.2, cex.axis = 1.2)
  
  rgb_vals <- col2rgb(jet_colors[2])
  polygon(c(rx, rev(rx)), c(upper, rev(lower)),
          col = rgb(rgb_vals[1]/255, rgb_vals[2]/255, rgb_vals[3]/255, 0.1), border = NA)
  
  lines(rx, mean_sim, col = jet_colors[2], lwd = 2, lty = 2)
  lines(rx, theo_line, col = "red", lwd = 2, lty = 2)
  
  legend("topleft", legend = c("Observed", "Mean Simulated", "Poisson"),
         col = c("black", jet_colors[2], "red"), lty = c(1, 2, 2), lwd = 2)
}


# New version ----
generate_envelope_plot <- function(csv_dir, statfun = "Lest", X, 
                                   window_size = 40, correction = "iso",
                                   main_title = "Envelope Plot",
                                   model_name = "Mean Model Simulation",
                                   line_color = "blue") {
  # Required packages
  # if (!requireNamespace("spatstat.geom", quietly = TRUE) ||
  #     !requireNamespace("spatstat.core", quietly = TRUE)) {
  #   stop("Please install the spatstat.geom and spatstat.core packages.")
  # }
  
  # Get list of simulation CSVs
  csv_files <- list.files(path = csv_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) stop("No CSV files found in the directory.")
  
  # Select function
  fun <- match.fun(statfun)
  
  # Reference r from observed pattern
  rx <- fun(X, correction = correction)$r
  obs_stat <- fun(X, correction = correction)[[correction]] - rx
  
  # Loop over all files and compute summary function
  stat_list <- vector("list", length(csv_files))
  for (i in seq_along(csv_files)) {
    df <- read.csv(csv_files[i], header = FALSE, skip = 1)
    colnames(df) <- c("x", "y")
    pp <- spatstat.geom::ppp(x = df$x, y = df$y, window = spatstat.geom::square(window_size))
    stat_list[[i]] <- fun(pp, correction = correction)[[correction]] - rx
  }
  
  # Combine results
  stat_matrix <- do.call(rbind, stat_list)
  stat_matrix <- t(stat_matrix)
  
  # Compute percentiles and mean
  lower <- apply(stat_matrix, 1, quantile, probs = 0.025)
  upper <- apply(stat_matrix, 1, quantile, probs = 0.975)
  mean_sim <- rowMeans(stat_matrix)
  theo_line <- fun(X, correction = correction)$theo - rx
  
  # Convert color to RGB for transparency
  rgb_vals <- col2rgb(line_color)
  fill_color <- rgb(rgb_vals[1]/255, rgb_vals[2]/255, rgb_vals[3]/255, 0.1)
  
  # Plot
  plot(rx, obs_stat, type = "l", col = "black", lwd = 2, lty = 1,
       ylim = range(c(lower, upper, obs_stat)), 
       ylab = expression(L(r) - r), xlab = expression(italic(r)),
       main = main_title,
       cex.lab = 1.2, cex.axis = 1.2)
  
  # Envelope
  polygon(c(rx, rev(rx)), c(upper, rev(lower)),
          col = fill_color, border = NA)
  
  # Mean line and theoretical line
  lines(rx, mean_sim, col = line_color, lwd = 2, lty = 2)
  lines(rx, theo_line, col = "red", lwd = 2, lty = 2)
  
  # Legend
  legend("topleft", legend = c("Data", model_name, "Poisson"),
         col = c("black", line_color, "red"), lty = c(1, 2, 2), lwd = 2,
         bty = "n", cex = 1.2)
}

# Observed data
# plot(X)
# csv_dir <- "/Users/alexander/Chalmers/MVEX11-25-18/dataframes/LGCP_punkter_2"

# obs_df <- read.csv("your_observed.csv", header = FALSE, skip = 1)
#colnames(obs_df) <- c("x", "y")
#X_obs <- spatstat.geom::ppp(x = obs_df$x, y = obs_df$y, window = spatstat.geom::square(40))

## Generate plot, works for L very nicely! ----
generate_envelope_plot(
  csv_dir = csv_dir,
  statfun = "Lest",
  correction = 'iso',
  X = X,
  main_title = "Konfidenshölje för Hierarkisk Poissonprocess",
  model_name = "Hierarkisk",
  line_color = jet_colors[10]
)


# for G envelope ---- 
generate_envelope_plot <- function(csv_dir, statfun = "Lest", X, 
                                   window_size = 40, correction = "iso",
                                   main_title = "Envelope Plot",
                                   model_name = "Mean Model Simulation",
                                   line_color = "blue",
                                   legend_pos = "topleft") {
  # Load CSV files
  csv_files <- list.files(path = csv_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) stop("No CSV files found in the directory.")
  
  # Match summary function
  fun <- match.fun(statfun)
  
  # Get observed summary function
  obs_fun <- fun(X, correction = correction)
  rx <- obs_fun$r
  
  # Decide whether to subtract r (only for Lest)
  subtract_r <- statfun == "Lest"
  
  obs_stat <- obs_fun[[correction]]
  if (subtract_r) obs_stat <- obs_stat - rx
  
  # Compute simulated summary functions
  stat_list <- vector("list", length(csv_files))
  for (i in seq_along(csv_files)) {
    df <- read.csv(csv_files[i], header = FALSE, skip = 1)
    colnames(df) <- c("x", "y")
    pp <- spatstat.geom::ppp(x = df$x, y = df$y, window = spatstat.geom::square(window_size))
    
    sim_val <- fun(pp, correction = correction)[[correction]]
    if (subtract_r) sim_val <- sim_val - rx
    stat_list[[i]] <- sim_val
  }
  
  # Combine and summarize
  stat_matrix <- t(do.call(rbind, stat_list))
  lower <- apply(stat_matrix, 1, quantile, probs = 0.025)
  upper <- apply(stat_matrix, 1, quantile, probs = 0.975)
  mean_sim <- rowMeans(stat_matrix)
  
  theo_line <- obs_fun$theo
  if (subtract_r) theo_line <- theo_line - rx
  
  # Prepare color
  rgb_vals <- col2rgb(line_color)
  fill_color <- rgb(rgb_vals[1]/255, rgb_vals[2]/255, rgb_vals[3]/255, 0.1)
  
  # Determine axis label
  ylab_expr <- switch(statfun,
                      "Lest" = expression(L(r) - r),
                      "Gest" = expression(G(r)),
                      "Fest" = expression(F(r)),
                      "Jest" = expression(J(r)),
                      expression(S(r)))
  
  # Plot
  plot(rx, obs_stat, type = "l", col = "black", lwd = 2, lty = 1,
       ylim = range(c(lower, upper, obs_stat)), 
       ylab = ylab_expr, xlab = expression(italic(r)),
       main = main_title,
       cex.lab = 1.2, cex.axis = 1.2)
  
  # Envelope polygon
  polygon(c(rx, rev(rx)), c(upper, rev(lower)),
          col = fill_color, border = NA)
  
  # Mean simulation and theoretical
  lines(rx, mean_sim, col = line_color, lwd = 2, lty = 2)
  lines(rx, theo_line, col = "red", lwd = 2, lty = 2)
  
  # Legend
  legend(legend_pos, legend = c("Data", model_name, "Poisson"),
         col = c("black", line_color, "red"), lty = c(1, 2, 2), lwd = 2,
         bty = "n", cex = 1.2)
}

# Updated function !
generate_envelope_plot <- function(csv_dir, statfun = "Lest", X, 
                                   window_size = 40, correction = "iso",
                                   main_title = "Envelope Plot",
                                   model_name = "Mean Model Simulation",
                                   line_color = "blue",
                                   legend_pos = "topleft",
                                   xlim = NULL) {
  # Load CSV files
  csv_files <- list.files(path = csv_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) stop("No CSV files found in the directory.")
  
  # Match summary function
  fun <- match.fun(statfun)
  
  # Observed pattern
  obs_fun <- fun(X, correction = correction)
  rx <- obs_fun$r
  subtract_r <- statfun == "Lest"
  
  obs_stat <- obs_fun[[correction]]
  if (subtract_r) obs_stat <- obs_stat - rx
  
  # Simulated values
  stat_list <- list()
  max_r_sim <- numeric(length(csv_files))  # for dynamic xlim
  
  for (i in seq_along(csv_files)) {
    df <- read.csv(csv_files[i], header = FALSE, skip = 1)
    colnames(df) <- c("x", "y")
    pp <- spatstat.geom::ppp(x = df$x, y = df$y, window = spatstat.geom::square(window_size))
    
    sim_fun <- fun(pp, correction = correction)
    r_vals <- sim_fun$r
    max_r_sim[i] <- max(r_vals)
    
    sim_vals <- sim_fun[[correction]]
    if (subtract_r) sim_vals <- sim_vals - r_vals
    
    stat_list[[i]] <- sim_vals
  }
  
  # Use common rx range (based on observed pattern)
  stat_matrix <- t(do.call(rbind, stat_list))
  lower <- apply(stat_matrix, 1, quantile, probs = 0.025)
  upper <- apply(stat_matrix, 1, quantile, probs = 0.975)
  mean_sim <- rowMeans(stat_matrix)
  
  theo_line <- obs_fun$theo
  if (subtract_r) theo_line <- theo_line - rx
  
  # Auto xlim from simulation if not provided
  if (is.null(xlim)) {
    max_r <- min(max(rx), min(max_r_sim))  # ensures consistency
    xlim <- c(0, max_r)
  }
  
  # Prepare color
  rgb_vals <- col2rgb(line_color)
  fill_color <- rgb(rgb_vals[1]/255, rgb_vals[2]/255, rgb_vals[3]/255, 0.1)
  
  # Y-label auto expression
  ylab_expr <- switch(statfun,
                      "Lest" = expression(L(r) - r),
                      "Gest" = expression(G(r)),
                      "Fest" = expression(F(r)),
                      "Jest" = expression(J(r)),
                      expression(S(r)))
  
  # Plot
  plot(rx, obs_stat, type = "l", col = "black", lwd = 2, lty = 1,
       ylim = range(c(lower, upper, obs_stat)), xlim = xlim,
       ylab = ylab_expr, xlab = expression(italic(r)),
       main = main_title,
       cex.lab = 1.2, cex.axis = 1.2)
  
  # Envelope polygon
  polygon(c(rx, rev(rx)), c(upper, rev(lower)),
          col = fill_color, border = NA)
  
  # Mean simulation and theoretical lines
  lines(rx, mean_sim, col = line_color, lwd = 2, lty = 2)
  lines(rx, theo_line, col = "red", lwd = 2, lty = 2)
  
  # Legend
  legend(legend_pos, legend = c("Data", model_name, "Poisson"),
         col = c("black", line_color, "red"), lty = c(1, 2, 2), lwd = 2,
         bty = "n", cex = 1.2)
}

# Usage 
# X_obs <- ppp(x = df_obs$x, y = df_obs$y, window = square(40))

# generate_envelope_plot(
#   csv_dir = "path/to/your/csvs",
#   statfun = "Gest",
#   X = X_obs,
#   main_title = "G-function Envelope (Cluster Model)",
#   model_name = "Cluster Simulation",
#   line_color = "#E69F00",
#   legend_pos = "bottomright"
# )

## Works for G! Nice! ----
generate_envelope_plot(
  csv_dir = csv_dir,
  statfun = "Gest",
  correction = 'rs',
  X = X,
  main_title = "Konfidenshölje för Hierarkisk Poissonprocess",
  model_name = "Hierarkisk",
  line_color = jet_colors[10],
  legend_pos = "bottomright",
  # xlim = c(0, 5)
)

plot(Gest(X,correction='rs'))

# Works for F?
generate_envelope_plot <- function(csv_dir, statfun = "Lest", X, 
                                   window_size = 40, correction = "iso",
                                   main_title = "Envelope Plot",
                                   model_name = "Mean Model Simulation",
                                   line_color = "blue",
                                   legend_pos = "topleft",
                                   xlim = NULL) {
  # Load CSV files
  csv_files <- list.files(path = csv_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) stop("No CSV files found in the directory.")
  
  # Match summary function
  fun <- match.fun(statfun)
  
  # Observed pattern
  obs_fun <- fun(X, correction = correction)
  rx <- obs_fun$r
  subtract_r <- statfun == "Lest"
  
  obs_stat <- obs_fun[[correction]]
  if (subtract_r) obs_stat <- obs_stat - rx
  
  # Simulated values
  stat_list <- list()
  max_r_sim <- numeric(length(csv_files))
  
  for (i in seq_along(csv_files)) {
    df <- read.csv(csv_files[i], header = FALSE, skip = 1)
    colnames(df) <- c("x", "y")
    pp <- spatstat.geom::ppp(x = df$x, y = df$y, window = spatstat.geom::square(window_size))
    
    sim_fun <- fun(pp, correction = correction)
    r_vals <- sim_fun$r
    max_r_sim[i] <- max(r_vals)
    
    sim_vals <- sim_fun[[correction]]
    if (subtract_r) sim_vals <- sim_vals - r_vals
    
    # Interpolate to match observed r if needed
    if (length(r_vals) != length(rx) || any(abs(r_vals - rx) > 1e-6)) {
      interp_vals <- approx(r_vals, sim_vals, xout = rx, rule = 2)$y
    } else {
      interp_vals <- sim_vals
    }
    
    stat_list[[i]] <- interp_vals
  }
  
  # Use common rx
  stat_matrix <- t(do.call(rbind, stat_list))
  lower <- apply(stat_matrix, 1, quantile, probs = 0.025)
  upper <- apply(stat_matrix, 1, quantile, probs = 0.975)
  mean_sim <- rowMeans(stat_matrix)
  
  theo_line <- obs_fun$theo
  if (subtract_r) theo_line <- theo_line - rx
  
  # Auto xlim
  if (is.null(xlim)) {
    max_r <- min(max(rx), min(max_r_sim))
    xlim <- c(0, max_r)
  }
  
  # Prepare color
  rgb_vals <- col2rgb(line_color)
  fill_color <- rgb(rgb_vals[1]/255, rgb_vals[2]/255, rgb_vals[3]/255, 0.1)
  
  # Y-label
  ylab_expr <- switch(statfun,
                      "Lest" = expression(L(r) - r),
                      "Gest" = expression(G(r)),
                      "Fest" = expression(F(r)),
                      "Jest" = expression(J(r)),
                      expression(S(r)))
  
  # Plot
  plot(rx, obs_stat, type = "l", col = "black", lwd = 2, lty = 1,
       ylim = range(c(lower, upper, obs_stat)), xlim = xlim,
       ylab = ylab_expr, xlab = expression(italic(r)),
       main = main_title,
       cex.lab = 1.2, cex.axis = 1.2)
  
  polygon(c(rx, rev(rx)), c(upper, rev(lower)),
          col = fill_color, border = NA)
  
  lines(rx, mean_sim, col = line_color, lwd = 2, lty = 2)
  lines(rx, theo_line, col = "red", lwd = 2, lty = 2)
  
  legend(legend_pos, legend = c("Data", model_name, "Poisson"),
         col = c("black", line_color, "red"), lty = c(1, 2, 2), lwd = 2,
         bty = "n", cex = 1.2)
}
## generate_envelope_plot: Updated again for ylim ----
generate_envelope_plot <- function(csv_dir, statfun = "Lest", X, 
                                   window_size = 40, correction = "iso",
                                   main_title = "Envelope Plot",
                                   model_name = "Mean Model Simulation",
                                   line_color = "blue",
                                   legend_pos = "topleft",
                                   xlim = NULL) {
  # Load CSV files
  csv_files <- list.files(path = csv_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) stop("No CSV files found in the directory.")
  
  # Match summary function
  fun <- match.fun(statfun)
  
  # Observed pattern
  obs_fun <- fun(X, correction = correction)
  rx <- obs_fun$r
  subtract_r <- statfun == "Lest"
  
  obs_stat <- obs_fun[[correction]]
  if (subtract_r) obs_stat <- obs_stat - rx
  
  # Simulated values
  stat_list <- list()
  max_r_sim <- numeric(length(csv_files))
  
  for (i in seq_along(csv_files)) {
    df <- read.csv(csv_files[i], header = FALSE, skip = 1)
    colnames(df) <- c("x", "y")
    pp <- spatstat.geom::ppp(x = df$x, y = df$y, window = spatstat.geom::square(window_size))
    
    sim_fun <- fun(pp, correction = correction)
    r_vals <- sim_fun$r
    max_r_sim[i] <- max(r_vals)
    
    sim_vals <- sim_fun[[correction]]
    if (subtract_r) sim_vals <- sim_vals - r_vals
    
    # Interpolate to match observed r if needed
    if (length(r_vals) != length(rx) || any(abs(r_vals - rx) > 1e-6)) {
      interp_vals <- approx(r_vals, sim_vals, xout = rx, rule = 2)$y
    } else {
      interp_vals <- sim_vals
    }
    
    stat_list[[i]] <- interp_vals
  }
  
  # Use common rx
  stat_matrix <- t(do.call(rbind, stat_list))
  lower <- apply(stat_matrix, 1, quantile, probs = 0.025)
  upper <- apply(stat_matrix, 1, quantile, probs = 0.975)
  mean_sim <- rowMeans(stat_matrix)
  
  theo_line <- obs_fun$theo
  if (subtract_r) theo_line <- theo_line - rx
  
  # Auto xlim
  if (is.null(xlim)) {
    max_r <- min(max(rx), min(max_r_sim))
    xlim <- c(0, max_r)
  }
  
  # Y-label
  ylab_expr <- switch(statfun,
                      "Lest" = expression(L(r) - r),
                      "Gest" = expression(G(r)),
                      "Fest" = expression(F(r)),
                      "Jest" = expression(J(r)),
                      expression(S(r)))
  
  # Force ylim for cumulative functions
  ylim_forced <- if (statfun %in% c("Gest", "Fest", "Jest")) c(0, 1) else range(c(lower, upper, obs_stat))
  
  # Prepare color
  rgb_vals <- col2rgb(line_color)
  fill_color <- rgb(rgb_vals[1]/255, rgb_vals[2]/255, rgb_vals[3]/255, 0.1)
  
  # Plot
  plot(rx, obs_stat, type = "l", col = "black", lwd = 2, lty = 1,
       ylim = ylim_forced, xlim = xlim,
       ylab = ylab_expr, xlab = expression(italic(r)),
       main = main_title,
       cex.main = 3, cex.lab = 2.5, cex.axis = 1.4,) # 2.5 lab, axis
       # cex.lab = 1.2, cex.axis = 1.2)
  
  polygon(c(rx, rev(rx)), c(upper, rev(lower)),
          col = fill_color, border = NA)
  
  lines(rx, mean_sim, col = line_color, lwd = 2, lty = 2)
  lines(rx, theo_line, col = "red", lwd = 2, lty = 2)
  
  legend(legend_pos, legend = c("Data", model_name, "Poisson"),
         col = c("black", line_color, "red"), lty = c(1, 2, 2), lwd = 2,
         bty = "n", cex = 2.5)
}

# Works for Fest nice! 
generate_envelope_plot(
  csv_dir = csv_dir,
  statfun = "Fest",
  correction = 'rs',
  X = X,
  main_title = "Konfidenshölje för Hierarkisk Poissonprocess",
  model_name = "Hierarkisk",
  line_color = jet_colors[10],
  legend_pos = "bottomright",
  xlim=c(0, 2.5)
)

plot(Fest(X,correction='rs'))

# Works for Jest nice! 
generate_envelope_plot(
  csv_dir = csv_dir,
  statfun = "Jest",
  correction = 'rs',
  X = X,
  main_title = "Konfidenshölje för Hierarkisk Poissonprocess",
  model_name = "Hierarkisk",
  line_color = jet_colors[10],
  legend_pos = "topright",
  xlim=c(0, 2.5)
)

# Alla samlade med senaste funktionen ----
setwd("/Users/alexander/Chalmers/MVEX11-25-18/fig/conf")
## hierarkisk ----
{png("conf-g-hierarkisk-col.png", width=694, height= 482)
generate_envelope_plot(
  csv_dir = csv_dir,
  statfun = "Gest",
  correction = 'rs',
  X = X,
  main_title = "Konfidenshölje för Hierarkisk Poissonprocess",
  model_name = "Hierarkisk",
  line_color = jet_colors[10],
  legend_pos = "bottomright",
  # xlim = c(0, 5)
)
dev.off()}

{png("conf-f-hierarkisk-col.png", width=694, height= 482)
  generate_envelope_plot(
    csv_dir = csv_dir,
    statfun = "Fest",
    correction = 'rs',
    X = X,
    main_title = "Konfidenshölje för Hierarkisk Poissonprocess",
    model_name = "Hierarkisk",
    line_color = jet_colors[10],
    legend_pos = "bottomright",
    xlim = c(0, 2.45) # leave off if possible
  )
  dev.off()
}

{png("conf-j-hierarkisk-col.png", width=694, height= 482)
  generate_envelope_plot(
    csv_dir = csv_dir,
    statfun = "Jest",
    correction = 'rs',
    X = X,
    main_title = "Konfidenshölje för Hierarkisk Poissonprocess",
    model_name = "Hierarkisk",
    line_color = jet_colors[10],
    legend_pos = "topright",
    xlim = c(0, 2.45) # leave off if possible
  )
  dev.off()
}

## Intensitetsstyrd Thomasprocess (mats final) ----
{ # run all 
  csv_dir <- "/Users/alexander/Chalmers/MVEX11-25-18/dataframes/punkterz_1" # Intensitets Poisson
  inthom_name <- "Intensitetsstyrd Thomasprocess"
  inthom_legend <- "Intensitets"
{png("conf-l-intensitets-col.png", width=694, height= 482)
  par(mar = c(5, 6.5, 4, 2) + 0.1)
  generate_envelope_plot(
    csv_dir = csv_dir,
    statfun = "Lest",
    correction = 'iso',
    X = X,
    main_title = inthom_name,
    model_name = inthom_legend,
    line_color = jet_colors[13],
    legend_pos = "bottomright",
    # xlim = c(0, 5)
  )
  dev.off()}

{png("conf-g-intensitets-col.png", width=694, height= 482)
  par(mar = c(5, 6.5, 4, 2) + 0.1)
  generate_envelope_plot(
    csv_dir = csv_dir,
    statfun = "Gest",
    correction = 'rs',
    X = X,
    main_title = inthom_name,
    model_name = inthom_legend,
    line_color = jet_colors[13],
    legend_pos = "bottomright",
    # xlim = c(0, 5)
  )
  dev.off()}

{png("conf-f-intensitets-col.png", width=694, height= 482)
  par(mar = c(5, 6.5, 4, 2) + 0.1)
  generate_envelope_plot(
    csv_dir = csv_dir,
    statfun = "Fest",
    correction = 'rs',
    X = X,
    main_title = inthom_name,
    model_name = inthom_legend,
    line_color = jet_colors[13],
    legend_pos = "bottomright",
    xlim = c(0, 2.45) # leave off if possible
  )
  dev.off()
}

{png("conf-j-intensitets-col.png", width=694, height= 482)
  par(mar = c(5, 6.5, 4, 2) + 0.1)
  generate_envelope_plot(
    csv_dir = csv_dir,
    statfun = "Jest",
    correction = 'rs',
    X = X,
    main_title = inthom_name,
    model_name = inthom_legend,
    line_color = jet_colors[13],
    legend_pos = "topright",
    # xlim = c(0, 2.45) # leave off if possible
  )
  dev.off()}
}

## 

## Iterativ Thomasprocess ----
{ # Run all 
  csv_dir <- "/Users/alexander/Chalmers/MVEX11-25-18/dataframes/path" # Iterativ Thomas
  model_main <- "Iterativ Thomasprocess"
  model_name <- "Iterativ"
{png("conf-l-iterativt-col.png", width=694, height= 482)
  par(mar = c(5, 6.5, 4, 2) + 0.1)
  generate_envelope_plot(
    csv_dir = csv_dir,
    statfun = "Lest",
    correction = 'iso',
    X = X,
    main_title = model_main,
    model_name = model_name,
    line_color = jet_colors[2],
    legend_pos = "bottomright",
    # xlim = c(0, 5)
  )
  dev.off()}
{png("conf-g-iterativt-col.png", width=694, height= 482)
  par(mar = c(5, 6.5, 4, 2) + 0.1)
  generate_envelope_plot(
    csv_dir = csv_dir,
    statfun = "Gest",
    correction = 'rs',
    X = X,
    main_title = model_main,
    model_name = model_name,
    line_color = jet_colors[2],
    legend_pos = "bottomright",
    # xlim = c(0, 5)
  )
  dev.off()}

{png("conf-f-iterativt-col.png", width=694, height= 482)
  par(mar = c(5, 6.5, 4, 2) + 0.1)
  generate_envelope_plot(
    csv_dir = csv_dir,
    statfun = "Fest",
    correction = 'rs',
    X = X,
    main_title = model_main,
    model_name = model_name,
    line_color = jet_colors[2],
    legend_pos = "bottomright",
    xlim = c(0, 2.45) # leave off if possible
  )
  dev.off()
}

{png("conf-j-iterativt-col.png", width=694, height= 482)
  par(mar = c(5, 6.5, 4, 2) + 0.1)
  generate_envelope_plot(
    csv_dir = csv_dir,
    statfun = "Jest",
    correction = 'rs',
    X = X,
    main_title = model_main,
    model_name = model_name,
    line_color = jet_colors[2],
    legend_pos = "topright",
    # xlim = c(0, 2.45) # leave off if possible
  )
  dev.off()}
}

# par(mar = c(5, 6.5, 4, 2) + 0.1)
## Herarkisk_3 Poisson, last mats  ----
# Last mats ---- run all pictures 
{ # Run all
  #csv_dir <- "/Users/alexander/Chalmers/MVEX11-25-18/dataframes/Herarkisk_3" # Iterativ Thomas
  csv_dir <- "/Users/alexander/Chalmers/MVEX11-25-18/dataframes/final_final_herarkisk"
  hier3 <- "Hierarkisk Poissonprocess"
  
  {png("conf-l-hierarkisk3-col.png", width=694, height= 482)
  par(mar = c(5, 6.5, 4, 2) + 0.1)
  generate_envelope_plot(
    csv_dir = csv_dir,
    statfun = "Lest",
    correction = 'iso',
    X = X,
    main_title = hier3,
    model_name = "Hierarkisk",
    line_color = jet_colors[10],
    legend_pos = "topleft",
    # xlim = c(0, 2.45) # leave off if possible
  )
  dev.off()
  }

{png("conf-g-hierarkisk3-col.png", width=694, height= 482)
  par(mar = c(5, 6.5, 4, 2) + 0.1)
  generate_envelope_plot(
    csv_dir = csv_dir,
    statfun = "Gest",
    correction = 'rs',
    X = X,
    main_title = hier3,
    model_name = "Hierarkisk",
    line_color = jet_colors[10],
    legend_pos = "bottomright",
    # xlim = c(0, 5)
  )
  dev.off()}

{png("conf-f-hierarkisk3-col.png", width=694, height= 482)
  par(mar = c(5, 6.5, 4, 2) + 0.1) # Shows the ylabel
  generate_envelope_plot(
    csv_dir = csv_dir,
    statfun = "Fest",
    correction = 'rs',
    X = X,
    main_title = hier3,
    model_name = "Hierarkisk",
    line_color = jet_colors[10],
    legend_pos = "bottomright",
    xlim = c(0, 2.45) # leave off if possible
  )
  dev.off()
}

# old_par <- par(no.readonly = TRUE)  # save current par settings
# par(mar = c(5, 6.5, 4, 2) + 0.1)
{png("conf-j-hierarkisk3-col.png", width=694, height= 482)
  par(mar = c(5, 6.5, 4, 2) + 0.1)
  generate_envelope_plot(
    csv_dir = csv_dir,
    statfun = "Jest",
    correction = 'rs',
    X = X,
    main_title = hier3,
    model_name = "Hierarkisk",
    line_color = jet_colors[10],
    legend_pos = "topright",
    xlim = c(0, 2.45) # leave off if possible
  )
  dev.off()
  }
}






# New generate_enve... see how it works for all data ----
generate_envelope_plot <- function(csv_dir, statfun = "Lest", X, 
                                   window_size = 40, correction = "iso",
                                   main_title = "Envelope Plot",
                                   model_name = "Mean Model Simulation",
                                   line_color = "blue",
                                   legend_pos = "topleft",
                                   xlim = NULL) {
  # Load CSV files
  csv_files <- list.files(path = csv_dir, pattern = "\\.csv$", full.names = TRUE)
  if (length(csv_files) == 0) stop("No CSV files found in the directory.")
  
  # Match summary function
  fun <- match.fun(statfun)
  
  # Observed pattern
  obs_fun <- fun(X, correction = correction)
  rx <- obs_fun$r
  subtract_r <- statfun == "Lest"
  
  obs_stat <- obs_fun[[correction]]
  if (subtract_r) obs_stat <- obs_stat - rx
  
  # Simulated values
  stat_list <- list()
  max_r_sim <- numeric(length(csv_files))
  
  for (i in seq_along(csv_files)) {
    df <- read.csv(csv_files[i], header = FALSE, skip = 1)
    colnames(df) <- c("x", "y")
    pp <- spatstat.geom::ppp(x = df$x, y = df$y, window = spatstat.geom::square(window_size))
    
    sim_fun <- fun(pp, correction = correction)
    r_vals <- sim_fun$r
    max_r_sim[i] <- max(r_vals)
    
    sim_vals <- sim_fun[[correction]]
    if (subtract_r) sim_vals <- sim_vals - r_vals
    
    # Interpolate to match observed r
    if (length(r_vals) != length(rx) || any(abs(r_vals - rx) > 1e-6)) {
      interp_vals <- approx(r_vals, sim_vals, xout = rx, rule = 2)$y
    } else {
      interp_vals <- sim_vals
    }
    
    stat_list[[i]] <- interp_vals
  }
  
  # Use common rx
  stat_matrix <- t(do.call(rbind, stat_list))
  lower <- apply(stat_matrix, 1, quantile, probs = 0.025)
  upper <- apply(stat_matrix, 1, quantile, probs = 0.975)
  mean_sim <- rowMeans(stat_matrix)
  
  theo_line <- obs_fun$theo
  if (subtract_r) theo_line <- theo_line - rx
  
  # Auto xlim
  if (is.null(xlim)) {
    max_r <- min(max(rx), min(max_r_sim))
    xlim <- c(0, max_r)
  }
  
  # Labels
  ylab_expr <- switch(statfun,
                      "Lest" = expression(L(r) - r),
                      "Gest" = expression(G(r)),
                      "Fest" = expression(F(r)),
                      "Jest" = expression(J(r)),
                      expression(S(r)))
  
  # Determine ylim and whether to manually control y-axis
  force_ylim_01 <- statfun %in% c("Gest", "Fest", "Jest")
  ylim_forced <- if (force_ylim_01) c(0, 1) else range(c(lower, upper, obs_stat))
  
  # Color conversion
  rgb_vals <- col2rgb(line_color)
  fill_color <- rgb(rgb_vals[1]/255, rgb_vals[2]/255, rgb_vals[3]/255, 0.1)
  
  # Plot
  plot(rx, obs_stat, type = "l", col = "black", lwd = 2, lty = 1,
       ylim = ylim_forced, xlim = xlim,
       ylab = ylab_expr, xlab = expression(italic(r)),
       main = main_title,
       cex.main = 3, cex.lab = 2.5, cex.axis = 2,
       yaxt = if (force_ylim_01) "n" else "s")
  
  # Force ticks at 0 and 1 for relevant functions
  if (force_ylim_01) axis(2, at = seq(0, 1, by = 0.2), labels = TRUE, las = 1, cex.axis = 2)
  
  # Envelope
  polygon(c(rx, rev(rx)), c(upper, rev(lower)),
          col = fill_color, border = NA)
  
  lines(rx, mean_sim, col = line_color, lwd = 2, lty = 2)
  lines(rx, theo_line, col = "red", lwd = 2, lty = 2)
  
  legend(legend_pos, legend = c("Data", model_name, "Poisson"),
         col = c("black", line_color, "red"), lty = c(1, 2, 2), lwd = 2,
         bty = "n", cex = 2.5)
}
# old_par2 <- par(no.readonly = TRUE)
# par_gfj <- par(mar = c(5, 6.5, 4, 2) + 0.1, mgp = c(3.5, 1, 0))
## Intensitets Poisson 
{ # run all 
  csv_dir <- "/Users/alexander/Chalmers/MVEX11-25-18/dataframes/punkterz_1" # Intensitets Poisson
  inthom_name <- "Intensitetsstyrd Thomasprocess"
  inthom_legend <- "Intensitets"
  {png("conf-l-intensitets-col.png", width=694, height= 482)
    par(mar = c(5, 6.5, 4, 2) + 0.1)
    generate_envelope_plot(
      csv_dir = csv_dir,
      statfun = "Lest",
      correction = 'iso',
      X = X,
      main_title = inthom_name,
      model_name = inthom_legend,
      line_color = jet_colors[13],
      legend_pos = "bottomright",
      # xlim = c(0, 5)
    )
    dev.off()}
  
  {png("conf-g-intensitets-col.png", width=694, height= 482)
    par(mar = c(5, 6.5, 4, 2) + 0.1)
    generate_envelope_plot(
      csv_dir = csv_dir,
      statfun = "Gest",
      correction = 'rs',
      X = X,
      main_title = inthom_name,
      model_name = inthom_legend,
      line_color = jet_colors[13],
      legend_pos = "bottomright",
      # xlim = c(0, 5)
    )
    dev.off()}
  
  {png("conf-f-intensitets-col.png", width=694, height= 482)
    par(mar = c(5, 6.5, 4, 2) + 0.1)
    generate_envelope_plot(
      csv_dir = csv_dir,
      statfun = "Fest",
      correction = 'rs',
      X = X,
      main_title = inthom_name,
      model_name = inthom_legend,
      line_color = jet_colors[13],
      legend_pos = "bottomright",
      xlim = c(0, 2.45) # leave off if possible
    )
    dev.off()
  }
  
  {png("conf-j-intensitets-col.png", width=694, height= 482)
    par(mar = c(5, 6.5, 4, 2) + 0.1)
    generate_envelope_plot(
      csv_dir = csv_dir,
      statfun = "Jest",
      correction = 'rs',
      X = X,
      main_title = inthom_name,
      model_name = inthom_legend,
      line_color = jet_colors[13],
      legend_pos = "topright",
      # xlim = c(0, 2.45) # leave off if possible
    )
    dev.off()}
}

## Iterativ Thomasprocess ----
{ # Run all 
  csv_dir <- "/Users/alexander/Chalmers/MVEX11-25-18/dataframes/path" # Iterativ Thomas
  model_main <- "Iterativ Thomasprocess"
  model_name <- "Iterativ"
  {png("conf-l-iterativt-col.png", width=694, height= 482)
    par(mar = c(5, 6.5, 4, 2) + 0.1)
    generate_envelope_plot(
      csv_dir = csv_dir,
      statfun = "Lest",
      correction = 'iso',
      X = X,
      main_title = model_main,
      model_name = model_name,
      line_color = jet_colors[2],
      legend_pos = "bottomright",
      # xlim = c(0, 5)
    )
    dev.off()}
  {png("conf-g-iterativt-col.png", width=694, height= 482)
    par(mar = c(5, 6.5, 4, 2) + 0.1)
    generate_envelope_plot(
      csv_dir = csv_dir,
      statfun = "Gest",
      correction = 'rs',
      X = X,
      main_title = model_main,
      model_name = model_name,
      line_color = jet_colors[2],
      legend_pos = "bottomright",
      xlim = c(0, 1.5)
    )
    dev.off()}
  
  {png("conf-f-iterativt-col.png", width=694, height= 482)
    par(mar = c(5, 6.5, 4, 2) + 0.1)
    generate_envelope_plot(
      csv_dir = csv_dir,
      statfun = "Fest",
      correction = 'rs',
      X = X,
      main_title = model_main,
      model_name = model_name,
      line_color = jet_colors[2],
      legend_pos = "bottomright",
      xlim = c(0, 2.45) # leave off if possible
    )
    dev.off()
  }
  
  {png("conf-j-iterativt-col.png", width=694, height= 482)
    par(mar = c(5, 6.5, 4, 2) + 0.1)
    generate_envelope_plot(
      csv_dir = csv_dir,
      statfun = "Jest",
      correction = 'rs',
      X = X,
      main_title = model_main,
      model_name = model_name,
      line_color = jet_colors[2],
      legend_pos = "topright",
      # xlim = c(0, 2.45) # leave off if possible
    )
    dev.off()}
}


## Herarkisk_3 Poisson, last mats  ----
# Last mats ---- run all pictures 
{ # Run all
  #csv_dir <- "/Users/alexander/Chalmers/MVEX11-25-18/dataframes/Herarkisk_3" # Iterativ Thomas
  csv_dir <- "/Users/alexander/Chalmers/MVEX11-25-18/dataframes/final_final_herarkisk"
  hier3 <- "Hierarkisk Poissonprocess"
  
  {png("conf-l-hierarkisk3-col.png", width=694, height= 482)
    par(mar = c(5, 6.5, 4, 2) + 0.1)
    generate_envelope_plot(
      csv_dir = csv_dir,
      statfun = "Lest",
      correction = 'iso',
      X = X,
      main_title = hier3,
      model_name = "Hierarkisk",
      line_color = jet_colors[10],
      legend_pos = "topleft",
      # xlim = c(0, 2.45) # leave off if possible
    )
    dev.off()
  }
  
  {png("conf-g-hierarkisk3-col.png", width=694, height= 482)
    par(mar = c(5, 6.5, 4, 2) + 0.1)
    generate_envelope_plot(
      csv_dir = csv_dir,
      statfun = "Gest",
      correction = 'rs',
      X = X,
      main_title = hier3,
      model_name = "Hierarkisk",
      line_color = jet_colors[10],
      legend_pos = "bottomright",
      # xlim = c(0, 5)
    )
    dev.off()}
  
  {png("conf-f-hierarkisk3-col.png", width=694, height= 482)
    par(mar = c(5, 6.5, 4, 2) + 0.1) # Shows the ylabel
    generate_envelope_plot(
      csv_dir = csv_dir,
      statfun = "Fest",
      correction = 'rs',
      X = X,
      main_title = hier3,
      model_name = "Hierarkisk",
      line_color = jet_colors[10],
      legend_pos = "bottomright",
      xlim = c(0, 2.45) # leave off if possible
    )
    dev.off()
  }
  
  # old_par <- par(no.readonly = TRUE)  # save current par settings
  # par(mar = c(5, 6.5, 4, 2) + 0.1)
  {png("conf-j-hierarkisk3-col.png", width=694, height= 482)
    par(mar = c(5, 6.5, 4, 2) + 0.1)
    generate_envelope_plot(
      csv_dir = csv_dir,
      statfun = "Jest",
      correction = 'rs',
      X = X,
      main_title = hier3,
      model_name = "Hierarkisk",
      line_color = jet_colors[10],
      legend_pos = "topright",
      xlim = c(0, 2.45) # leave off if possible
    )
    dev.off()
  }
}


# Try things out ----- 
## Hierarkisk Poissonprocess ----
{ # Run all
  #csv_dir <- "/Users/alexander/Chalmers/MVEX11-25-18/dataframes/Herarkisk_3" # Iterativ Thomas
  csv_dir <- "/Users/alexander/Chalmers/MVEX11-25-18/dataframes/final_final_herarkisk"
  hier3 <- "Hierarkisk Poissonprocess"
  
  {
    png("conf-g-hierarkisk3-col1.png", width=694, height= 482)
    # png("conf-g-hierarkisk3-col1.png", width=726, height= 457)
    par(mar = c(5, 6.5, 4, 2) + 0.1, mgp = c(3.5, 1, 0))
    generate_envelope_plot(
      csv_dir = csv_dir,
      statfun = "Gest",
      correction = 'rs',
      X = X,
      main_title = hier3,
      model_name = "Hierarkisk",
      line_color = jet_colors[10],
      legend_pos = "bottomright"
    )
    dev.off()
  }

}
