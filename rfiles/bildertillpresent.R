pp <- runifpoint(100, square(40), nsim=10)

plot(pp[[6]], main='')


# 2D grid from 0 to 10 in both x and y directions, with 10 steps each
x <- seq(0, 10, length.out = 10)
y <- seq(0, 10, length.out = 10)

ppsquare <- ppp(x,y, square(10))
plot(ppsquare)
grid_points <- expand.grid(x = x, y = y)
grid_p <- ppp(grid_points$x, grid_points$y,square(10))
plot(grid_p, main='')
plot(grid_points, pch = 16, asp = 1)


pp <- rpoispp(100, square(40))

plot(pp, main='')


X
plot(Kest(X, correction = 'Ripley'),main='')

plot(Lest(X, correction = 'Ripley'), main = '')
plot(Fest(X, correction = 'rs'), main = '')
plot(Gest(X, correction = 'rs'), main = '')
plot(Jest(X, correction = 'rs'), main = '')

# Calculate Lest
L <- Lest(X, correction = 'Ripley')

# Linearize by subtracting r
L$iso <- L$iso - L$r
L$trans <- L$trans - L$r
L$Ripley <- L$Ripley - L$r


setwd('/Users/alexander/Library/CloudStorage/OneDrive-Personal/__kandidatfigpresentaiton')
# Plot the linearized function
{
  png("L-data.png", width = 673, height = 517)
  par(mar = c(5, 4.5, 4, 2) + 0.1, mgp = c(2.5, 1, 0))
  plot(L$r, L$iso, main = '', type = 'l', lwd = 2, 
       ylab = expression(italic("L(r) - r")),
       xlab = expression(italic("r")))
  lines(L$r, L$theo - L$r, col = 'red', lwd = 2, lty = 2)
  
  legend("topleft",
         legend = c('Data', 'Poisson'),
         col = c('black', 'red'),
         lty = c(1, 2),
         lwd = c(2, 2),
         bty = "n", cex = 1.3)
  dev.off()
}
# ---- Kest ----
{
  K <- Kest(X, correction = "Ripley")
  
  png("K-data.png", width = 673, height = 517)
  par(mar = c(5, 4.5, 4, 2) + 0.1, mgp = c(2.5, 1, 0))
  plot(K$r, K$iso, type = 'l', lwd = 2,
       ylab = expression(italic("K(r)")),
       xlab = expression(italic("r")),
       main = '')
  lines(K$r, K$theo, col = 'red', lwd = 2, lty = 2)
  legend("topleft",
         legend = c("Data", "Poisson"),
         col = c("black", "red"),
         lty = c(1, 2),
         lwd = c(2, 2),
         bty = "n", cex = 1.3)
  dev.off()
}
# ---- Gest ----
{
  G <- Gest(X, correction = "rs")
  
  png("G-data.png", width = 673, height = 517)
  par(mar = c(5, 4.5, 4, 2) + 0.1, mgp = c(2.5, 1, 0))
  plot(G$r, G$rs, type = 'l', lwd = 2,
       ylab = expression(italic("G(r)")),
       xlab = expression(italic("r")),
       main = '')
  lines(G$r, G$theo, col = 'red', lwd = 2, lty = 2)
  legend("topleft",
         legend = c("Data", "Poisson"),
         col = c("black", "red"),
         lty = c(1, 2),
         lwd = c(2, 2),
         bty = "n", cex = 1.3)
  dev.off()
}

# ---- Fest ----
{
  F <- Fest(X, correction = "rs")
  
  png("F-data.png", width = 673, height = 517)
  par(mar = c(5, 4.5, 4, 2) + 0.1, mgp = c(2.5, 1, 0))
  plot(F$r, F$rs, type = 'l', lwd = 2,
       ylab = expression(italic("F(r)")),
       xlab = expression(italic("r")),
       ylim=c(0,1),
       main = '')
  lines(F$r, F$theo, col = 'red', lwd = 2, lty = 2)
  legend("topleft",
         legend = c("Data", "Poisson"),
         col = c("black", "red"),
         lty = c(1, 2),
         lwd = c(2, 2),
         bty = "n", cex = 1.3)
  dev.off()
}
# ---- Jest ----
{
  J <- Jest(X, correction = "rs")
  
  png("J-data.png", width = 673, height = 517)
  par(mar = c(5, 4.5, 4, 2) + 0.1, mgp = c(2.5, 1, 0))
  plot(J$r, J$rs, type = 'l', lwd = 2,
       ylab = expression(italic("J(r)")),
       xlab = expression(italic("r")),
       main = '')
  lines(J$r, J$theo, col = 'red', lwd = 2, lty = 2)
  legend("topright",
         legend = c("Data", "Poisson"),
         col = c("black", "red"),
         lty = c(1, 2),
         lwd = c(2, 2),
         bty = "n", cex = 1.3)
  dev.off()
}


{
  png("thom-halvvags.png", width = 673, height = 517)
  # halfway_thomas <- ppp(offspring_x, offspring_y, square(40))
  plot(halfway_thomas, main='') # main="Hybrid Thomas Process (Averaged Parameters)"
  dev.off()
}

X_ordinary_thomas <- simulate(fitT,nsim=9)
simulated_ppp <- X_ordinary_thomas[[7]] # [[1]]
{
  png("thom-thomass.png", width = 673, height = 517)
  plot(simulated_ppp, main='')
  dev.off()
}



{
  png("thom-omodifierad.png", width = 673, height = 517)
  plot(thomas_process, main='')
  dev.off()
}


# Fix the rest -----
plot_summary_stat <- function(point_patterns, pattern_names, colors,
                              method = "Kest", xmax = 30,
                              ylim = NULL, main_title = NULL,
                              legend_position = "topright",
                              cex_main = 2.5, cex_lab = 2, cex_axis = 1.8,
                              cex_legend = 1.5,
                              envelope_first = FALSE, nsim = 99) {
  
  stopifnot(length(point_patterns) == length(pattern_names),
            length(colors) == length(pattern_names))
  
  method_function <- match.fun(method)
  y_component <- if (method %in% c("Fest", "Gest", "Jest")) "rs" else "iso"
  
  if (envelope_first) {
    env <- envelope(point_patterns[[1]], method_function, nsim = nsim, verbose = TRUE)
    
    # Transform to L(r) - r if Lest
    if (method == "Lest") {
      env$obs <- env$obs - env$r
      env$theo <- env$theo - env$r
      env$lo <- env$lo - env$r
      env$hi <- env$hi - env$r
    }
    
    ylab_expr <- switch(method,
                        "Lest" = expression(L(r) - r),
                        "Kest" = expression(K(r)),
                        "Fest" = expression(F(r)),
                        "Gest" = expression(G(r)),
                        "Jest" = expression(J(r)),
                        bquote(.(as.name(method))(r)))
    
    if (is.null(main_title)) {
      method_label <- switch(method,
                             "Kest" = "K-funktion",
                             "Lest" = "L-funktion",
                             "Fest" = "F-funktion",
                             "Gest" = "G-funktion",
                             "Jest" = "J-funktion",
                             method)
      main_title <- paste(method_label, "för processer")
    }
    
    plot(env, . ~ r, main = main_title, legend = FALSE,
         xlim = c(0, xmax), ylim = ylim,
         cex.main = cex_main, cex.lab = cex_lab, cex.axis = cex_axis,
         xlab = expression(italic(r)), ylab = ylab_expr,
         col = colors[1], lwd = 2)
    
    # Optional: add other patterns
    if (length(point_patterns) > 1) {
      for (i in 2:length(point_patterns)) {
        stat <- method_function(point_patterns[[i]], correction = y_component)
        if (method == "Lest") stat[[y_component]] <- stat[[y_component]] - stat$r
        lines(stat$r, stat[[y_component]], col = colors[i], lwd = 2)
      }
    }
    
    # Add Poisson theoretical line
    lines(env$r, env$theo, lty = 2, col = "red", lwd = 2)
    
    legend(legend_position,
           legend = c(pattern_names, "Poisson"),
           col = c(colors, "red"),
           lty = c(rep(1, length(pattern_names)), 2),
           lwd = c(rep(2.5, length(pattern_names)), 2.5),
           bty = "n", cex = cex_legend)
    
  } else {
    # Standard (non-envelope) plotting
    results <- lapply(point_patterns, function(pp) method_function(pp, correction = y_component))
    
    if (method == "Lest") {
      results <- lapply(results, function(res) {
        res[[y_component]] <- res[[y_component]] - res$r
        res
      })
    }
    
    if (is.null(main_title)) {
      method_label <- switch(method,
                             "Kest" = "K-funktion",
                             "Lest" = "L-funktion",
                             "Fest" = "F-funktion",
                             "Gest" = "G-funktion",
                             "Jest" = "J-funktion",
                             method)
      main_title <- paste(method_label, "för processer")
    }
    
    ylab_expr <- switch(method,
                        "Lest" = expression(L(r) - r),
                        "Kest" = expression(K(r)),
                        "Fest" = expression(F(r)),
                        "Gest" = expression(G(r)),
                        "Jest" = expression(J(r)),
                        bquote(.(as.name(method))(r)))
    
    r_vals <- results[[1]]$r
    y_vals <- results[[1]][[y_component]]
    
    plot(r_vals, y_vals, type = 'l',
         main = main_title,
         xlab = expression(italic(r)), ylab = ylab_expr,
         xlim = c(0, xmax), ylim = ylim,
         col = colors[1], lwd = 3,
         cex.main = cex_main, cex.lab = cex_lab, cex.axis = cex_axis)
    
    for (i in 2:length(results)) {
      lines(results[[i]]$r, results[[i]][[y_component]],
            col = colors[i], lwd = 2)
    }
    
    theo_line <- if (method == "Lest") results[[1]]$theo - results[[1]]$r else results[[1]]$theo
    lines(r_vals, theo_line, lwd = 2.5, lty = 2, col = 'red')
    
    legend(legend_position,
           legend = c(pattern_names, 'Poisson'),
           col = c(colors, 'red'),
           lty = c(rep(1, length(pattern_names)), 2),
           lwd = c(rep(2.5, length(pattern_names)), 2.5),
           bty = "n", cex = cex_legend)
  }
}


{
  png("thom-all.png", width = 650, height = 457)
  par(mar = c(5, 6.5, 4, 2) + 0.1, mgp = c(2.5, 1, 0))
  # x_j <- Jest(X, correction = 'rs')
  xmax_j <- max(x_j$r)
  plot_summary_stat(
    point_patterns = list(X, simulated_ppp, halfway_thomas, thomas_process),
    pattern_names = c("Data", "Thomas", "Halvvägs", "Omodifierad"),
    colors = c("black", jet_colors[6], jet_colors[7], jet_colors[8]),
    method = "Lest",
    xmax = 10, # 10
    # ylim = c(0, 1),
    legend_position = "topleft",
    cex_legend = 1.3, 
    cex_main = 0.01,
  )
  dev.off()
}
# # L -data 
# {
#   png("L-data.png", width = 650, height = 457)
#   par(mar = c(5, 6.5, 4, 2) + 0.1, mgp = c(2.5, 1, 0))
#   # x_j <- Jest(X, correction = 'rs')
#   xmax_j <- max(x_j$r)
#   plot_summary_stat(
#     point_patterns = list(X),
#     pattern_names = c("Data", NULL),
#     colors = c("black", NULL),
#     method = "Lest",
#     xmax = 10, # 10
#     # ylim = c(0, 1),
#     legend_position = "topleft",
#     cex_legend = 1.3, 
#     cex_main = 0.01,
#   )
#   # lines(r_vals, theo_line, lwd = 2.5, lty = 2, col = 'red')
#   abline(h = 0, lwd=2.5, lty=2, col="red")
#   dev.off()
# }




# all l ---- 
{
  png("L-all-really.png", width = 650, height = 457)
  par(mar = c(5, 6.5, 4, 2) + 0.1, mgp = c(2.5, 1, 0))
  # x_j <- Jest(X, correction = 'rs')
  xmax_j <- max(x_j$r)
  plot_summary_stat(
    point_patterns = list(X, simulated_ppp, halfway_thomas, thomas_process, X_iter10, pp_hier, pp_intp),
    pattern_names = c("Data", "Thomas", "Halvvägs", "Omodifierad", "Iterativ", "Hierarkisk", "Intensitetsstyrd"),
    colors = c("black", jet_colors[6], jet_colors[7], jet_colors[8], jet_colors[2], jet_colors[10], jet_colors[13]),
    method = "Lest",
    xmax = 10, # 10
    ylim = c(0, 5),
    legend_position = "topleft",
    cex_legend = 1.1, 
    cex_main = 0.01,
  )
  dev.off()
}


{
  png("jall.png", width = 650, height = 457)
  par(mar = c(5, 6.5, 4, 2) + 0.1, mgp = c(2.5, 1, 0))
  x_j <- Jest(X, correction = 'rs')
  xmax_j <- max(x_j$r)
  plot_summary_stat(
    point_patterns = list(X, simulated_ppp, halfway_thomas, X_iter10, pp_hier, pp_intp),
    pattern_names = c("Data", "Thomas", "Halvvägs", "Iterativ", "Hierarkisk", "Intensitetsstyrd"),
    colors = c("black", jet_colors[6], jet_colors[7], jet_colors[2], jet_colors[10], jet_colors[13]),
    method = "Jest",
    xmax = xmax_j, # 10
    ylim = c(0, 1),
    legend_position = "topright",
    cex_main = 2,
  )
  dev.off()
}



Vinf <- rSSI(0.07)

V100 <- rSSI(0.07, 100)

X <- runifpoint(100)
Y <- rSSI(0.05,100,x.init=X) # Y consists of X together with
# 42 added points.
plot(Y, main="rSSI")
plot(X,add=TRUE,chars=20,cols="red")

L_matern <- Lest(Y,correction='Ripley')


plot(L_matern)

## inhomogeneous
Z <- rSSI(0.07, 50, f=function(x,y){x})
plot(Z)


# matern 2 ----
X <- rMaternII(50, 0.05)
Y <- rMaternII(20, 0.05, stationary=FALSE)
# matern_x <- X # store the value 
plot(matern_x, main='Matern 2, radie = 0.05') # matern 2
plot(Y)
