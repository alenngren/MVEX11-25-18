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
