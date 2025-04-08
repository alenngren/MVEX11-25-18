# Load data

library(readr)
library(spatstat)
iterpath <- '/Users/alexander/Chalmers/MVEX11-25-18/dataframes/unif-iter10.csv'
iterframe <- read_csv(iterpath, col_names = c("x", "y"))
plot(iterframe)

filename <- "/Users/alexander/Chalmers/MVEX11-25-18/ERIKA.rds"
ERIKA <- readRDS(filename)
x <- ERIKA$VES13$small$x
y <- ERIKA$VES13$small$y
X <- ppp(x=x, y=y, window=square(40))

# Plot observed K-function > L-function START HERE! WORKS 
K_obs <- Kest(X)
enve <- envelope(X, Lest, nsim=99)
#setEPS()
#postscript("CSR-L-func.eps")
png("CSR-L-func.png", width = 1000, height = 600)
plot(K_obs$r, sqrt(K_obs$iso/pi)-K_obs$r, type = "l", col = "red", lwd = 1,
     # ylim = range(c(K_obs$iso + hi_diff, K_obs$iso + lo_diff)),
     main = "CSR Envelope Around Observed L-function",
     #xlab = "r", ylab = "K(r)"
     )



#lines(enve$r, (enve$hi-enve$theo+enve$obs)-enve$r) # add the line we want 
#lines(enve$r, (enve$lo-enve$theo+enve$obs)-enve$r)

# Add shaded CSR deviation envelope around the observed line
#polygon(c(r_vals, rev(r_vals)),
#        c(K_obs$iso + hi_diff, rev(K_obs$iso + lo_diff)),
#        col = rgb(0.8, 0.8, 0.8, 0.5), border = NA)

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
       col = c("red", rgb(0.8, 0.8, 0.8)), lwd = c(2, NA),) 

# setwd("/Users/alexander/Chalmers/MVEX11-25-18/fig")




dev.off()
       
       