library(spatstat)
library(mgcv)

library(spatstat)
# install.packages("feather")
# install.packages("spatstat")


filename <- "/Users/alexander/Chalmers/MVEX11-25-18/ERIKA.rds"

ERIKA <- readRDS(filename)

x <- ERIKA$VES13$small$x
y <- ERIKA$VES13$small$y



events <- data.frame(x = x, y = y) # alex fix 


# events <- read.csv("C:/Users/matti/Downloads/VES13_small.csv")
events$pt <- 1      
events$wt <- 1e-6

xrange <- range(events$x)
yrange <- range(events$y)
area   <- diff(xrange) * diff(yrange)
W      <- owin(xrange, yrange)  

x_big <- ERIKA$VES13$large$x
y_big <- ERIKA$VES13$large$y
dbh_big <- ERIKA$VES13$large$dbh
bigTrees <- data.frame(x=x_big, y=y_big, dbh=dbh_big) # my version 
bigTrees$dbh
# bigTrees <- read.csv("C:/Users/matti/Downloads/VES13_large.csv")

print(names(bigTrees))

t <- 5  

computeCovariate <- function(sx, sy, bigTrees, t) {
  dx <- sx - bigTrees$x
  dy <- sy - bigTrees$y
  d  <- sqrt(dx^2 + dy^2)
  sum( bigTrees$dbh * exp(- d / t) )
}

q <- 10000
quad <- data.frame(
  x  = runif(q, min = xrange[1], max = xrange[2]),
  y  = runif(q, min = yrange[1], max = yrange[2]),
  pt = 0,            
  wt = area / q      
)

dat <- rbind(events, quad)

dat$C <- mapply(
  computeCovariate, 
  dat$x, 
  dat$y,
  MoreArgs = list(bigTrees = bigTrees, theta = t)
)

fit <- gam(
  pt / wt ~ C + s(x, y, bs = "gp", k = 200),
  data    = dat,
  family  = poisson,
  weights = wt,
  method  = "REML"
)


grid_x <- seq(xrange[1], xrange[2], length.out = 200)
grid_y <- seq(yrange[1], yrange[2], length.out = 200)
grid   <- expand.grid(x = grid_x, y = grid_y)
grid$pt <- 0  
grid$wt <- 1  

grid$C <- mapply(
  computeCovariate,
  grid$x,
  grid$y,
  MoreArgs = list(bigTrees = bigTrees, theta = t)
)

grid$intensity <- predict(fit, newdata = grid, type = "response")
im_intensity <- as.im(
  data.frame(x = grid$x, y = grid$y, value = grid$intensity),
  W = W
)
plot(im_intensity, zlim = c(0, 6), main = "Fält med kovariater")
points(events$x, events$y, pch = 20, col = "black")

summary(fit)
AIC(fit)
