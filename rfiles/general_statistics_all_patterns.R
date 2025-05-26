setwd("/Users/alexander/Chalmers/MVEX11-25-18/fig/all")


plot_summary_stat <- function(point_patterns, pattern_names, colors,
                              method = "Kest", xmax = 30,
                              ylim = NULL, main_title = NULL,
                              legend_position = "topleft") {
  stopifnot(length(point_patterns) == length(pattern_names),
            length(colors) == length(pattern_names))
  
  method_function <- switch(method,
                            "Kest" = Kest,
                            "Lest" = Lest,
                            "Gest" = Gest,
                            "Fest" = Fest,
                            "Jest" = Jest,
                            stop("Unsupported method"))
  
  y_component <- if (method %in% c("Fest", "Gest", "Jest")) "rs" else "iso"
  print(y_component)
  results <- lapply(point_patterns, function(pp) method_function(pp, correction = y_component))
  
  # For L-function, subtract r from iso
  if (method == "Lest") {
    results <- lapply(results, function(res) {
      res$iso <- res$iso - res$r
      res
    })
  }
  
  # Determine which component to plot (iso or rs)
  # y_component <- if (method %in% c("Fest", "Gest", "Jest")) "rs" else "iso"
  
  if (is.null(main_title)) {
    main_title <- switch(method,
                         "Kest" = "K-funktion för processer",
                         "Lest" = "L-funktion för processer",
                         "Fest" = "F-funktion för processer",
                         "Gest" = "G-funktion för processer",
                         "Jest" = "J-funktion för processer",
                         paste(method, "för processer"))
  }
  ylab_expr <- if (method == "Lest") expression(L(r) - r) else bquote(.(as.name(method))(r))
  
  plot(results[[1]]$r, results[[1]][[y_component]], type = 'l',
       main = main_title, xlab = expression(r), ylab = ylab_expr,
       cex.main = 2, cex.lab = 1.35, cex.axis = 1.5,
       xlim = c(0, xmax), ylim = ylim,
       col = colors[1], lwd = 2)
  
  for (i in 2:length(results)) {
    lines(results[[i]]$r[results[[i]]$r <= xmax],
          results[[i]][[y_component]][results[[i]]$r <= xmax],
          type = 'l', col = colors[i], lwd = 1)
  }
  
  # Adjust theoretical line for Lest
  theo_line <- if (method == "Lest") results[[1]]$theo - results[[1]]$r else results[[1]]$theo
  
  lines(results[[1]]$r, theo_line, type = 'l', lwd = 2, lty = 2, col = 'red')
  
  legend(legend_position,
         legend = c(pattern_names, 'Poisson'),
         col = c(colors, 'red'),
         lty = c(rep(1, length(pattern_names)), 2),
         lwd = c(2, rep(2, length(pattern_names) - 1), 2),
         bty = "n", cex = 1.2)  # <-- Adjust legend size here if needed
} 


plot_summary_stat <- function(point_patterns, pattern_names, colors,
                              method = "Kest", xmax = 30,
                              ylim = NULL, main_title = NULL,
                              legend_position = "topleft") {
  stopifnot(length(point_patterns) == length(pattern_names),
            length(colors) == length(pattern_names))
  
  # Select method function and y-component
  method_function <- match.fun(method)
  y_component <- if (method %in% c("Fest", "Gest", "Jest")) "rs" else "iso"
  
  # Apply summary function
  results <- lapply(point_patterns, function(pp) method_function(pp, correction = y_component))
  
  # If L-function, subtract r from component
  if (method == "Lest") {
    results <- lapply(results, function(res) {
      res[[y_component]] <- res[[y_component]] - res$r
      res
    })
  }
  
  # Default main title
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
  
  # Y-axis label
  ylab_expr <- if (method == "Lest") expression(L(r) - r) else bquote(.(as.name(method))(r))
  
  # Plot the first line
  r_vals <- results[[1]]$r
  y_vals <- results[[1]][[y_component]]
  plot(r_vals, y_vals, type = 'l',
       main = main_title, xlab = expression(r), ylab = ylab_expr,
       cex.main = 2, cex.lab = 1.35, cex.axis = 1.5,
       xlim = c(0, xmax), ylim = ylim,
       col = colors[1], lwd = 2)
  
  # Overlay additional patterns
  for (i in 2:length(results)) {
    lines(results[[i]]$r, results[[i]][[y_component]],
          col = colors[i], lwd = 1)
  }
  
  # Theoretical line
  theo_line <- if (method == "Lest") results[[1]]$theo - results[[1]]$r else results[[1]]$theo
  lines(r_vals, theo_line, lwd = 2, lty = 2, col = 'red')
  
  # Legend
  legend(legend_position,
         legend = c(pattern_names, 'Poisson'),
         col = c(colors, 'red'),
         lty = c(rep(1, length(pattern_names)), 2),
         lwd = c(rep(2, length(pattern_names)), 2),
         bty = "n", cex = 1.5)
}


plot_summary_stat <- function(point_patterns, pattern_names, colors,
                              method = "Kest", xmax = 30,
                              ylim = NULL, main_title = NULL,
                              legend_position = "topright") {
  stopifnot(length(point_patterns) == length(pattern_names),
            length(colors) == length(pattern_names))
  
  method_function <- match.fun(method)
  y_component <- if (method %in% c("Fest", "Gest", "Jest")) "rs" else "iso"
  
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
  
  #ylab_expr <- if (method == "Lest") expression(L(r) - r) else bquote(.(as.name(method))(r))
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
       cex.main = 2.5, cex.lab = 2, cex.axis = 1.8)
  
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
         bty = "n", cex = 1.5)
}

plot_summary_stat <- function(point_patterns, pattern_names, colors,
                              method = "Kest", xmax = 30,
                              ylim = NULL, main_title = NULL,
                              legend_position = "topright",
                              cex_main = 2.5) {
  stopifnot(length(point_patterns) == length(pattern_names),
            length(colors) == length(pattern_names))
  
  method_function <- match.fun(method)
  y_component <- if (method %in% c("Fest", "Gest", "Jest")) "rs" else "iso"
  
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
       cex.main = cex_main, cex.lab = 2, cex.axis = 1.8)
  
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
         bty = "n", cex = 1.5)
}


# End function ----
# All for kall, lall, gall, fall, jall ----
{
  {
    png("kall.png", width = 650, height = 457)
    par(mar = c(5, 6.5, 4, 2) + 0.1, mgp = c(2.5, 1, 0))
    plot_summary_stat(
      point_patterns = list(X, simulated_ppp, halfway_thomas, X_iter10, pp_hier, pp_intp),
      pattern_names = c("Data", "Thomas", "Halvvägs", "Iterativ", "Hierarkisk", "Intensitets"),
      colors = c("black", jet_colors[6], jet_colors[7], jet_colors[2], jet_colors[10], jet_colors[13]),
      method = "Kest",
      xmax = 10,
      ylim = c(0, 660),
      legend_position = "topleft",
      cex_main = 2.5
    )
    dev.off()
  }
  
  {
    # png("lall.png", width = 650, height = 457)
    par(mar = c(5, 6.5, 4, 2) + 0.1, mgp = c(2.5, 1, 0))
    plot_summary_stat(
      point_patterns = list(X, simulated_ppp, halfway_thomas, X_iter10, pp_hier, pp_intp, thomas_process),
      pattern_names = c("Data", "Thomas", "Halvvägs", "Iterativ", "Hierarkisk", "Intensitets", "Omodifierad"),
      colors = c("black", jet_colors[6], jet_colors[7], jet_colors[2], jet_colors[10], jet_colors[13]),
      method = "Lest",
      xmax = 10,
      ylim = c(0, 6),
      legend_position = "topleft",
      cex_main = 2.5
    )
    dev.off()
  }
  
  
  {
    png("gall.png", width = 650, height = 457)
    par(mar = c(5, 6.5, 4, 2) + 0.1, mgp = c(2.5, 1, 0))
    x_g <- Gest(X, correction = 'rs')
    xmax_g <- max(x_g$r)
    
    plot_summary_stat(
      point_patterns = list(X, simulated_ppp, halfway_thomas, X_iter10, pp_hier, pp_intp),
      pattern_names = c("Data", "Thomas", "Halvvägs", "Iterativ", "Hierarkisk", "Intensitets"),
      colors = c("black", jet_colors[6], jet_colors[7], jet_colors[2], jet_colors[10], jet_colors[13]),
      method = "Gest",
      xmax = xmax_g, # 10
      ylim = c(0, 1),
      legend_position = "bottomright",
      cex_main = 2
    )
    
    dev.off()
  }
  
  {
    png("fall.png", width = 650, height = 457)
    xmax_f <- max(x_f$r)
    par(mar = c(5, 6.5, 4, 2) + 0.1, mgp = c(2.5, 1, 0))
    plot_summary_stat(
      point_patterns = list(X, simulated_ppp, halfway_thomas, X_iter10, pp_hier, pp_intp),
      pattern_names = c("Data", "Thomas", "Halvvägs", "Iterativ", "Hierarkisk", "Intensitets"),
      colors = c("black", jet_colors[6], jet_colors[7], jet_colors[2], jet_colors[10], jet_colors[13]),
      method = "Fest",
      xmax = xmax_f, # 10
      ylim = c(0, 1),
      legend_position = "topleft",
      cex_main = 2
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

}




# Small fixes with plots! ----
plot_summary_stat <- function(point_patterns, pattern_names, colors,
                              method = "Kest", xmax = 30,
                              ylim = NULL, main_title = NULL,
                              legend_position = "topright",
                              cex_main = 2.5, cex_lab = 2, cex_axis = 1.8,
                              cex_legend = 1.5) {
  stopifnot(length(point_patterns) == length(pattern_names),
            length(colors) == length(pattern_names))
  
  method_function <- match.fun(method)
  y_component <- if (method %in% c("Fest", "Gest", "Jest")) "rs" else "iso"
  
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

# Iterativ -----
library(readr)
library(spatstat)
iter0 <- '/Users/alexander/Chalmers/MVEX11-25-18/dataframes/unif-iter0-aut.csv'
iter5 <- '/Users/alexander/Chalmers/MVEX11-25-18/dataframes/unif-iter5-aut.csv'
iter10 <- '/Users/alexander/Chalmers/MVEX11-25-18/dataframes/unif-iter10-aut.csv'
iterframe0 <- read_csv(iter0, col_names = c("x", "y"))
iterframe5 <- read_csv(iter5, col_names = c("x", "y"))
iterframe10 <- read_csv(iter10, col_names = c("x", "y"))

X_iter0 <-  ppp(x=iterframe0$x, y=iterframe0$y, window=square(40))
X_iter5 <-  ppp(x=iterframe5$x, y=iterframe5$y, window=square(40))
X_iter10 <-  ppp(x=iterframe10$x, y=iterframe10$y, window=square(40))

{
  png("Iterativ_all.png", width = 650, height = 457)
  par(mar = c(5, 6.5, 4, 2) + 0.1, mgp = c(2.5, 1, 0))
  x_l <- Jest(X, correction = 'rs')
  xmax_l <- max(x_l$r)
  plot_summary_stat(
    point_patterns = list(X, X_iter0, X_iter5, X_iter10),
    pattern_names = c("Data", "Iterativ 1", "Iterativ 6", "Iterativ 11"),
    colors = c("black", jet_colors[2], jet_colors[3], jet_colors[4]),
    method = "Lest",
    xmax = 10, # 10
    # ylim = c(0, 1),
    legend_position = "topleft",
    cex_main = 1.5,
    cex_lab = 1.5, cex_axis = 1.2, cex_legend = 1.3,
    main_title = "Iterativ Thomasprocess"
  )
  dev.off()
}



# Crazy change for envelope -----
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
    env <- envelope(point_patterns[[1]], method_function, nsim = nsim, verbose = TRUE, correction = y_component)
    
    if (method == "Lest") {
      env[[y_component]] <- env[[y_component]] - env$r
      env$obs <- env$obs - env$r
      env$theo <- env$theo - env$r
      env$lo <- env$lo - env$r
      env$hi <- env$hi - env$r
    }
    
    plot(env, .[[y_component]] ~ r,
         main = main_title %||% paste(method, "envelope for", pattern_names[1]),
         xlim = c(0, xmax), ylim = ylim,
         legend = FALSE,
         cex.main = cex_main, cex.lab = cex_lab, cex.axis = cex_axis,
         col = colors[1], lwd = 3,
         xlab = expression(italic(r)),
         ylab = switch(method,
                       "Lest" = expression(L(r) - r),
                       "Kest" = expression(K(r)),
                       "Fest" = expression(F(r)),
                       "Gest" = expression(G(r)),
                       "Jest" = expression(J(r)),
                       bquote(.(as.name(method))(r))))
  } else {
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

## Run the code ----
{
  png("Iterativ_all.png", width = 650, height = 457)
  par(mar = c(5, 6.5, 4, 2) + 0.1, mgp = c(2.5, 1, 0))
  x_l <- Jest(X, correction = 'rs')
  xmax_l <- max(x_l$r)
  plot_summary_stat(
    point_patterns = list(X, X_iter0, X_iter5, X_iter10),
    pattern_names = c("Data", "Iterativ 1", "Iterativ 6", "Iterativ 11"),
    colors = c("black", jet_colors[2], jet_colors[3], jet_colors[4]),
    method = "Lest",
    xmax = 10, # 10
    # ylim = c(0, 1),
    legend_position = "topleft",
    cex_main = 1.5,
    cex_lab = 1.5, cex_axis = 1.2, cex_legend = 1.3,
    main_title = "Iterativ Thomasprocess",
    envelope_first = TRUE
  )
  dev.off()
}


plot(X_k, sqrt(./pi) -r  ~ r, type='l', main='L-funktion för processer', legend=FALSE, xlim=c(0,xmax), ylim=c(0,8), lwd=2, ylab=expression(L(r) -r))
