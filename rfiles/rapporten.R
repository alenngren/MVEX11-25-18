# Load data
setwd("/Users/alexander/Chalmers/MVEX11-25-18/fig")
library(readr)
library(spatstat)
iterpath <- '/Users/alexander/Chalmers/MVEX11-25-18/dataframes/unif-iter10.csv'
path_csv <- "/Users/alexander/Chalmers/MVEX11-25-18/python-data/ves13-thomasprocess-modified.csv"


iterframe <- read_csv(iterpath, col_names = c("x", "y"))
plot(iterframe)


# Uniform plot -----------
# Number of points
n <- 886
# Define the window (bounding box)
win <- owin(xrange=c(0, 40), yrange=c(0, 40))
# Generate uniformly distributed points
pp_uniform <- runifpoint(n, win=win)
# Plot the result
plot(pp_uniform, main='880 punkter ~ uniform')


filename <- "/Users/alexander/Chalmers/MVEX11-25-18/ERIKA.rds"
ERIKA <- readRDS(filename)
x <- ERIKA$VES13$small$x
y <- ERIKA$VES13$small$y
x_p <- ERIKA$VES13$large$x
y_p <- ERIKA$VES13$large$y
X <- ppp(x=x, y=y, window=square(40))
X_p <- ppp(x=x_p, y=y_p, window=square(40))

# Data VES13 postscript("LE-halv-thomas.eps", width = 6, height = 3.7)
setEPS()
postscript("ves13.eps", width = 6, height = 6) # square
par(mar=c(0,0,0,0))
plot(X, main="") # all the children points 
dev.off()

setEPS()
postscript("ves13_p.eps", width = 6, height = 6) # square
plot(X_p, main="") # parents points 
dev.off()

par(
  mfrow=c(1,1),     # One plot per window (default is no grid)
  mar=c(5,4,4,2),   # Margins: bottom, left, top, right (default values)
  oma=c(0,0,0,0),   # Outer margins: top, right, bottom, left (default is none)
  cex=1,            # Character expansion (default text size)
  las=0,            # Axis label orientation (default is horizontal)
  bty="o",          # Box type (default is "o" - open box)
  xpd=FALSE,        # No clipping to the plot region
  tcl=-0.5,         # Tick mark length (default)
  pch=1,            # Default plotting character (empty circle)
  lwd=1,            # Line width (default is 1)
  col="black"       # Default color (black)
)

setEPS()
postscript("ves13_double.eps") # square , width = 6, height = 3.7
png("ves13_double.png", width = 800, height = 640)
par(mfrow= c(1,2))
par(mfrow=c(1,2), mar=c(0,0,2,0)) # , oma=c(0,0,2,0)
plot(X_p, main="Föräldrarpunkter") # parents points 
plot(X, main="Dotterpunkter ") # all the children points 
dev.off()

par(mfrow=c(1, 2), mar=c(4, 4, 2, 1), asp=1) # mar adjusts the margins

# First plot
# plot(x_p, y_p, col='black', main='Plot 1: y = x^2', xlab='X', ylab='Y')
plot(x_p, y_p, col='black', main='Föräldrarpunkter', xaxt = "n", yaxt = "n",
     xlab = "", ylab = "")

# plot(X_p, main='Plot 1: y = x^2', xaxt = "n", yaxt = "n")

# Second plot
plot(x, y, col='black', main='Dotterpunkter', xaxt = "n", yaxt = "n",
     xlab = "", ylab = "")

# plot(0:40, 0:40, main="", type="n")
# text(x=20, y=5, "VES13", cex=8, col=rgb(0.9, 0.9, 0.9, 0.3))  # big, light, transparent
# points(X)

all_plot <- allstats(X, envelope=TRUE, nsim=99)
plot(all_plot,legend=FALSE, main="")
plot(all_plot,legend=TRUE, main="")

# Plot förklara uniform
setEPS()
postscript("motivera_kluster.eps") # square , width = 6, height = 3.7
#par(mfrow=c(1,3))
par(mfrow=c(1,3),            # 1 row, 3 columns
    mar=c(2, 2, 3, 1),       # margins: bottom, left, top, right (smaller inner margins)
    oma=c(0, 0, 0, 0))       # outer 
plot(X_p, main="Föräldrarpunkter")
plot(X, main="Dotterpunkter") 
plot(pp_uniform, main='Uniform fördelning')
dev.off()





# Plot observed K-function > L-function START HERE! WORKS 
K_obs <- Kest(X)

# Poisson 
en <- envelope(X, nsim=99)
plot(en, main="VES13 - Envelope of K function with CSR")

plot( rThomas(intensity(X), 22, 12))

fitT <- kppm(X ~ 1, "Thomas")
X_ordinary_thomas <- simulate(fitT,nsim=1)
simulated_ppp <- X_ordinary_thomas[[1]] # happy with this 
plot(simulated_ppp, main='thomas') # nsim = 9 # plot Thomas 
points_matrix <- coords(simulated_ppp)
print(points_matrix)
X_thomas <- simulate(fitT,nsim=1)
plot(simulated_ppp, main='thomas') # nsim = 9
fitL <- kppm(X ~ 1, "LGCP")
plot(simulate(fitL,nsim=9), 'thomas')
fitT # uniform intensity = l

# Poisson omodifierad
window <- owin(c(0, 40), c(0, 40)) #creating the window
lambda_p_hat <- 72 / (40*40) #intensity of the parent process
mu_hat <- 882 / 74  #mean number of offspring

#calculating standard deviation
offspring_x <- ERIKA$VES13$small$x
offspring_y <- ERIKA$VES13$small$y
parent_x <- ERIKA$VES13$large$x
parent_y <- ERIKA$VES13$large$y
distances <- sqrt((offspring_x - parent_x)^2 + (offspring_y - parent_y)^2)
sigma_hat <- sd(distances)

#Simulate parent points
num_parents <- rpois(1, lambda_p_hat * area.owin(window))  # Poisson-distributed parent count
parents <- runifpoint(num_parents, win = window)  # Generate parent points

#Simulating offspring around parents
offspring_x <- c()
offspring_y <- c()

for (i in 1:num_parents) {
  num_offspring <- rpois(1, mu_hat)  # Poisson-distributed number of offspring
  if (num_offspring > 0) {
    # Generate offspring positions using normal distribution around parent
    x_offset <- rnorm(num_offspring, mean = 0, sd = sigma_hat)
    y_offset <- rnorm(num_offspring, mean = 0, sd = sigma_hat)
    
    new_x <- parents$x[i] + x_offset
    new_y <- parents$y[i] + y_offset
    
    # Keep points inside window
    valid <- inside.owin(new_x, new_y, window)
    offspring_x <- c(offspring_x, new_x[valid])
    offspring_y <- c(offspring_y, new_y[valid])
  }
}

# Creating a ppp for the process and plotting
thomas_process <- ppp(offspring_x, offspring_y, window=window)  # thomas unmodief 
plot(thomas_process, main="Simulated Unmodified Thomas Process")
points(parents, col="red", pch=3)  # Parent points in red

K_thomas <- Kest(thomas_process, correction="Ripley")
plot(K_thomas, sqrt(./pi) -r ~ r, ylab="L(r) - r ", main="L-r, un-modified Thomas")

fitUT <- kppm(thomas_process ~ 1, "Thomas")
fitUT


# Halfway Thomas -----------
lambda_avg <- 2 / (1/lambda_p_hat + 1/0.003943241)  # Harmonic mean
mu_avg <- sqrt(mu_hat * 139.7962)  # Geometric mean
sigma_avg <- sqrt(sigma_hat * 2.097576813)  # Geometric mean, skattad för halvvägs

num_parents_avg <- rpois(1, lambda_avg * area.owin(window)) # skattad för halvvägs
parents_avg <- runifpoint(num_parents_avg, win = window) #Simulate parents using average intensity

offspring_x <- c()
offspring_y <- c()

for (i in 1:num_parents_avg) {
  num_offspring <- rpois(1, mu_avg)
  if (num_offspring > 0) {
    x_offset <- rnorm(num_offspring, mean = 0, sd = sigma_avg)
    y_offset <- rnorm(num_offspring, mean = 0, sd = sigma_avg)
    
    new_x <- parents_avg$x[i] + x_offset
    new_y <- parents_avg$y[i] + y_offset
    
    # Keep points inside window
    valid <- inside.owin(new_x, new_y, window)
    offspring_x <- c(offspring_x, new_x[valid])
    offspring_y <- c(offspring_y, new_y[valid])
  }
}

#Create point pattern object and plot
halfway_thomas <- ppp(offspring_x, offspring_y, window=window) # lisa mellan model
plot(halfway_thomas, main="Hybrid Thomas Process (Averaged Parameters)")
# points(parents_avg, col="red", pch=3)  # Parent points in red


# Iterativ med DATA 
iter0 <- '/Users/alexander/Chalmers/MVEX11-25-18/dataframes/unif-iter0-aut.csv'
iter5 <- '/Users/alexander/Chalmers/MVEX11-25-18/dataframes/unif-iter5-aut.csv'
iter10 <- '/Users/alexander/Chalmers/MVEX11-25-18/dataframes/unif-iter10-aut.csv'
iterframe0 <- read_csv(iter0, col_names = c("x", "y"))
iterframe5 <- read_csv(iter5, col_names = c("x", "y"))
iterframe10 <- read_csv(iter10, col_names = c("x", "y"))

X_iter0 <-  ppp(x=iterframe0$x, y=iterframe0$y, window=square(40))
X_iter5 <-  ppp(x=iterframe5$x, y=iterframe5$y, window=square(40))
X_iter10 <-  ppp(x=iterframe10$x, y=iterframe10$y, window=square(40))

X_iter <-  ppp(x=iterframe$x, y=iterframe$y, window=square(40))

E_0 <- envelope(X_iter0, Kest, nsim=999, verbose=TRUE)
E_5 <- envelope(X_iter5, Kest, nsim=999, verbose=TRUE)
E_10 <- envelope(X_iter10, Kest, nsim=999, verbose=TRUE)

E_iter <- envelope(X_iter, Kest, nsim=9)

E <- envelope(X, Kest, nsim=99, verbose=TRUE)

envelope(X, Kest, )

E_iter <- envelope(X_iter, Kest, nsim=9)

K_csv <- Kest(X_csv, correction="Ripley") # CSV data 

plot(E, sqrt(./pi) -r  ~ r, main = 'L ~ iterativ steg', legend = FALSE) 
lines(E_0$r, sqrt(E_0$obs/pi) -E_0$r , col = jet_colors[3], lty = 1, lwd = 2) 
lines(E_5$r, sqrt(E_5$obs/pi) -E_5$r , col = jet_colors[2], lty = 2, lwd = 2) 
lines(E_10$r, sqrt(E_10$obs/pi) -E_10$r , col = jet_colors[1], lty = 2, lwd = 2) 
lines(E_half$r, sqrt(E_half$obs/pi) -E_half$r, col = jet_colors[4], lty = 4, lwd = 2)

enve <- envelope(X, Lest, nsim=99)
#setEPS()
#postscript("CSR-L-func.eps")
# png("CSR-L-func.png", width = 1000, height = 600)
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
plot(E, sqrt(./pi) -r  ~ r, main = 'L ~ iterativ steg', legend = FALSE, cex.main = 1) 
lines(E_0$r, sqrt(E_0$obs/pi) -E_0$r , col = jet_colors[3], lty = 1, lwd = 2) 
lines(E_5$r, sqrt(E_5$obs/pi) -E_5$r , col = jet_colors[2], lty = 2, lwd = 2) 
lines(E_10$r, sqrt(E_10$obs/pi) -E_10$r , col = jet_colors[1], lty = 2, lwd = 2) 
# lines(E_half$r, sqrt(E_half$obs/pi) -E_half$r, col = jet_colors[4], lty = 4, lwd = 2)


text(9, sqrt(max(E$obs)/pi) - 10.1, "VES13", col = "black", font = 1, cex = 1)
text(9, sqrt(max(E_0$obs)/pi) - 9.8, "Iterativ 1", col = jet_colors[3], font = 1, cex = 1)
text(9, sqrt(max(E_5$obs)/pi) - 10.4, "Iterativ 6", col = jet_colors[2], font = 1, cex = 1)
text(9, sqrt(max(E_10$obs)/pi) - 10.2, "Iterativ 11", col = jet_colors[1], font = 1, cex = 1)

legend("topleft", legend = c("VES13", "Iterativ 1", "Iterativ 6", "Iterativ 11"),
       col = c('black', jet_colors[3],jet_colors[2],jet_colors[1]), lty = c(1,2,2,2 )) 

# text(9, sqrt(max(E_half$obs)/pi) - 10.2, "E_half", col = jet_colors[4], font = 1, cex = 1)

#legend("topleft", legend = c("Observed", "CSR Envelope"),
       col = c("red", rgb(0.8, 0.8, 0.8)), lwd = c(2, NA),) 

# setwd("/Users/alexander/Chalmers/MVEX11-25-18/fig")
# dev.off()


# K funktion plot 




# Clark Evans statistics 
clarkevans.test(X_iter, alternative="clustered")
clarkevans.test(X_iter0, alternative="clustered")
clarkevans.test(X_iter5, alternative="clustered")
clarkevans.test(X_iter10, alternative="clustered")
       

# LGPC modell (mats) -----
window <- owin(range(x), range(y))

ppp_data <- ppp(x, y, window = window)
densitet <- density(ppp_data, sigma = 0.6) 
kordinater <- as.data.frame(densitet)
colnames(kordinater) <- c("x", "y", "ca")

kordinater$ca <- pmax(kordinater$ca, 0)
lgcp <- gam(ca ~ s(x, y, k = 50), family = poisson, data = kordinater)
kordinater$vänte_lam <- exp(predict(lgcp, newdata = kordinater))
skal <- 0.6 
kordinater$skalad_lamda1 <- skal * pmax(kordinater$vänte_lam - 0.05, 0.1)
matris_skalad <- matrix(kordinater$skalad_lamda1, 
                        nrow = dim(densitet$v)[1], 
                        ncol = dim(densitet$v)[2], 
                        byrow = FALSE)
skalad_lamda2 <- im(matris_skalad, xcol = densitet$xcol, yrow = densitet$yrow)
LGCP <- rpoispp(lambda = skalad_lamda2, win = window) # X patter



# Mats metod för kluster - ej data baserad -----
win <- owin(xrange = c(0, 40), yrange = c(0, 40))
m <- as.im(function(x, y) {4 - 1.5 * ((x / 40) - 0.5)^2 + 2 * ((y / 40) - 0.5)^2}, W = win)
X_rlgcp <- rLGCP("gauss", m, var = 0.6, scale = 10, win = win, saveLambda = TRUE)
lambda_field <- (attr(X_rlgcp, "Lambda") - min(attr(X_rlgcp, "Lambda"))) / (max(attr(X_rlgcp, "Lambda")) - min(attr(X_rlgcp, "Lambda"))) * 2.5
plot(lambda_field, main = ":D")
win <- Window(lambda_field)
n_f <- 160
f_punkt <- data.frame(
  x = runif(n_f, min = win$xrange[1], max = win$xrange[2]),
  y = runif(n_f, min = win$yrange[1], max = win$yrange[2])
)
densitet <- as.function(lambda_field)


f_punkt$intensity <- densitet(f_punkt$x, f_punkt$y)
hist(f_punkt$intensity, main = "Histogram predicted intensiteter",
     xlab = "Predicted intensiteter", col = "lightblue", breaks = 20)

t <- 0.3 # se slides
A <- 15  

f_punkt$lambda <- ifelse(f_punkt$intensity <= t,
                         0,
                         A * (f_punkt$intensity - t))
f_punkt$n_children <- rpois(n_f, lambda = f_punkt$lambda)
hist(f_punkt$lambda, main = "Histogram lambda",
     xlab = "lambda(A*(Intensitet - t))", col = "lightgreen", breaks = 20)

barn_std <- 0.8 # varians barn

barn <- do.call(rbind, lapply(1:n_f, function(i) {
  nchild <- f_punkt$n_children[i]
  if (nchild > 0) {
    data.frame(
      x = rnorm(nchild, mean = f_punkt$x[i], sd = barn_std),
      y = rnorm(nchild, mean = f_punkt$y[i], sd = barn_std),
      parent = i
    )
  } else {
    NULL
  }
}))

plot(lambda_field, main  = "")
points(f_punkt$x, f_punkt$y, col = "black", pch = 1)
if (!is.null(barn)) {
  points(barn$x, barn$y, col = "black", pch = 1)
}

x2 <- range(f_punkt$x, if (!is.null(barn)) barn$x)
y2 <- range(f_punkt$y, if (!is.null(barn)) barn$y)
plot(NA, xlim = x2, ylim = y2, type = "n", axes = FALSE, xlab = "", ylab = "", frame.plot = FALSE)
points(f_punkt$x, f_punkt$y, col = "black", pch = 1)
points(barn$x, barn$y, col = "black", pch = 1)

Xlambda <- ppp(x=barn$x, y=barn$y, window=square(40))
Xlambda

# Alla point patterns X ppp 
thomas_process
halfway_thomas
simulated_ppp
X_iter10
Xlambda
LGCP

# F function, empty space statistics -----
Fest(X, ..., eps, r=NULL, breaks=NULL,
     correction=c("rs", "km", "cs"),
     domain=NULL)

plot(Fest(X))
X_f <- Fest(X, correction = 'rs')
Xthomas_f <- Fest(simulated_ppp, correction = 'rs') 
Xhalf_f <- Fest(halfway_thomas, correction = 'rs')
Xiter_f <- Fest(X_iter10, correction = 'rs')
Xlgcp_f <- Fest(LGCP)



par(
  mfrow=c(1,1),     # One plot per window (default is no grid)
  mar=c(5,4,4,2),   # Margins: bottom, left, top, right (default values)
  oma=c(0,0,0,0),   # Outer margins: top, right, bottom, left (default is none)
  cex=1,            # Character expansion (default text size)
  las=0,            # Axis label orientation (default is horizontal)
  bty="o",          # Box type (default is "o" - open box)
  xpd=FALSE,        # No clipping to the plot region
  tcl=-0.5,         # Tick mark length (default)
  pch=1,            # Default plotting character (empty circle)
  lwd=1,            # Line width (default is 1)
  col="black"       # Default color (black)
)
library(fields)
jet_colors <- tim.colors(5)
plot(X_f, main = " F-function", legend = FALSE, , col=jet_colors[1])
lines(Xthomas_f, col=jet_colors[2])
lines(Xhalf_f, col=jet_colors[3])
lines(Xiter_f, col=jet_colors[4])
lines(Xlgcp_f, col=jet_colors[5])

text(9, max(X_f), "data", col = "black", font = 1, cex = 1)
text(9, max(Xthomas_f), "E_1", col = jet_colors[3], font = 1, cex = 1)
text(9, max(Xhalf_f), "E_6", col = jet_colors[2], font = 1, cex = 1)
text(9, max(Xiter_f), "E_11", col = jet_colors[1], font = 1, cex = 1)
text(9, max(Xlgcp_f), "E_half", col = jet_colors[4], font = 1, cex = 1)

legend("bottomright", legend = c("data", "border", "Thomas", "Halv-thomas", "Iterativ", "lgcp"),
       col = c(jet_colors[1], 'red', jet_colors[2], jet_colors[3], jet_colors[4], jet_colors[5]), 
       lwd = c(1,1,1,1,1,1),
       lty=c(1,2,1,1,1))



fitData <- kppm(X ~ 1, "Thomas")
fit <- ppm(X ~ 1, "Thomas") # doesnt work 

en_data_t <- envelope(fitData,Lest)

lgcp_ppp <- as.ppp(lgcp)
k_lgcp <- Kest(lgcp)

plot(en_data_t, main = 'Data and iterativ')
lines(E_10$r, sqrt(E_10$obs/pi), col = jet_colors[2], lty = 1, lwd = 1)
lines(K_obs$r, sqrt(K_obs$iso/pi), col = "red", lwd = 1) # - K_obs$r

thomas_process
halfway_thomas
simulated_ppp
X_iter10
Xlambda
LGCP

# fit envelopes
fit_tp <- kppm(thomas_process ~ 1, "Thomas")
fit_thom <- kppm(simulated_ppp ~ 1, "Thomas")
fit_ht <- kppm(halfway_thomas ~ 1, "Thomas")
fit_iter10 <- kppm(X_iter10 ~ 1, "Thomas")
fit_lam <- kppm(Xlambda ~ 1, "Thomas")
fit_lgcp <- kppm(LGCP ~ 1, "Thomas")

# fit <- ppm(X ~ 1, "Thomas") # doesnt work 

env_tp <- envelope(fit_tp, Lest, nsim=99)
env_thom <- envelope(fit_tp, Lest, nsim=99)
env_ht <-  envelope(fit_ht, Lest, nsim=99)
env_iter10 <- envelope(fit_iter10, Lest, nsim=99)
env_lam <- envelope(fit_lam, Lest, nsim=99)
env_lgcp <- envelope(fit_lgcp, Lest, nsim=99)


# setwd("/Users/alexander/Chalmers/MVEX11-25-18/fig")
# png("Omodifierad.png")
setEPS()
postscript("LE-Omodifierad.eps", width = 6, height = 3.7)
plot(env_tp - E$r, main = 'Omodifierad Thomas', legend=FALSE, col = "blue", ylab=expression(L(r) - r))
lines(K_obs$r, sqrt(K_obs$iso/pi)- E$r, col = "black", lwd = 1)  # data under lest 
# lines(E_10$r, sqrt(E_10$obs/pi)- E$r, col = jet_colors[2], lty = 1, lwd = 1)
dev.off()

legend("topleft", legend = c("data", "Omodifierad"),
       col = c('black', 'blue'), 
       lwd = c(1,1),
       lty=c(1,1))
dev.off()

setEPS()
postscript("LE-vanlig-thomas.eps", width = 6, height = 3.7)
plot(env_thom - E$r, main = 'Vanlig Thomas', legend=FALSE, col = "blue",
     ylab=expression(L(r) - r))
lines(K_obs$r, sqrt(K_obs$iso/pi)- E$r, col = "black", lwd = 1)  # data under lest 
# lines(E_10$r, sqrt(E_10$obs/pi)- E$r, col = jet_colors[2], lty = 1, lwd = 1)

legend("topleft", legend = c("data", "Vanlig"),
       col = c('black', 'blue'), 
       lwd = c(1,1),
       lty=c(1,1))
dev.off()


setEPS()
postscript("LE-halv-thomas.eps", width = 6, height = 3.7)
plot(env_ht - E$r, main = 'Halv Thomas', legend=FALSE, col = "blue",
     ylab=expression(L(r) - r))
lines(K_obs$r, sqrt(K_obs$iso/pi)- E$r, col = "black", lwd = 1)  # data under lest 
# lines(E_10$r, sqrt(E_10$obs/pi)- E$r, col = jet_colors[2], lty = 1, lwd = 1)

legend("topleft", legend = c("data", "Halv"),
       col = c('black', 'blue'), 
       lwd = c(1,1),
       lty=c(1,1))
dev.off()


setEPS()
postscript("LE-iterativ10.eps", width = 6, height = 3.7)
plot(env_iter10 - E$r, main = 'Envelop of iterativ, blå = data ', legend=FALSE, col = "blue",
     ylab=expression(L(r) - r))
lines(K_obs$r, sqrt(K_obs$iso/pi)- E$r, col = "black", lwd = 1)  # data under lest 
# lines(E_10$r, sqrt(E_10$obs/pi)- E$r, col = jet_colors[2], lty = 1, lwd = 1)

legend("topleft", legend = c("data", "Iterativ"),
       col = c('black', 'blue'), 
       lwd = c(1,1),
       lty=c(1,1))
dev.off()

setEPS()
postscript("LE-lambda.eps", width = 6, height = 3.7)
plot(env_lam - E$r, main = 'lambda latent fält', legend=FALSE, col = "blue",
     ylab=expression(L(r) - r))
lines(K_obs$r, sqrt(K_obs$iso/pi)- E$r, col = "black", lwd = 1)  # data under lest 
# lines(E_10$r, sqrt(E_10$obs/pi)- E$r, col = jet_colors[2], lty = 1, lwd = 1)

legend("topleft", legend = c("data", "latent"),
       col = c('black', 'blue'), 
       lwd = c(1,1),
       lty=c(1,1))
dev.off()

setEPS()
postscript("LE-lgcp.eps", width = 6, height = 3.7)
plot(env_lgcp - E$r, main = 'lgcp', legend=FALSE, col = "blue",
     ylab=expression(L(r) - r))
lines(K_obs$r, sqrt(K_obs$iso/pi)- E$r, col = "black", lwd = 1)  # data under lest 
# lines(E_10$r, sqrt(E_10$obs/pi)- E$r, col = jet_colors[2], lty = 1, lwd = 1)

legend("topleft", legend = c("data", "LGCP"),
       col = c('black', 'blue'), 
       lwd = c(1,1),
       lty=c(1,1))
dev.off()









# Plot different Thomas processes ----
par(mfrow=c(1,4))

plot(thomas_process, main="Unmodified Thomas Process")
plot(simulated_ppp, main='Thomas') # nsim = 9 # plot Thomas 
plot(halfway_thomas, main="Halvvägs Thomas")
plot(X, main = 'VES13')

# works for thomasprocesses
png("widescreen_plot.png", width = 2400, height = 600, res = 150)

par(mfrow = c(1, 4), mar = c(2, 2, 2, 1), oma = c(0, 0, 2, 0))
par(mfrow = c(1, 4), mar = c(0, 0, 2, 0), oma = c(0, 0, 2, 0), cex.main = 2.5)

plot(thomas_process, main = "Omodifierad")
plot(simulated_ppp, main = "Thomas")
plot(halfway_thomas, main = "Halvvägs")
plot(X, main = "VES13")

dev.off()


# Plot L for different Thomasprocesses 
E_unthomas <- envelope(thomas_process, Kest, nsim=99, verbose=TRUE)
E_thomas <- envelope(simulated_ppp, Kest, nsim=99, verbose=TRUE)
E_halfway <- envelope(halfway_thomas, Kest, nsim=99, verbose=TRUE)

plot(E, sqrt(./pi) -r  ~ r, main = 'L ~ Olika Thomasprocesser', legend = FALSE, ylim = c(0, 6), cex.main = 1) 
lines(E_unthomas$r, sqrt(E_unthomas$obs/pi) -E_unthomas$r , col = jet_colors[3], lty = 1, lwd = 2) 
lines(E_thomas$r, sqrt(E_thomas$obs/pi) -E_thomas$r , col = jet_colors[2], lty = 2, lwd = 2) 
lines(E_halfway$r, sqrt(E_halfway$obs/pi) -E_halfway$r , col = jet_colors[1], lty = 2, lwd = 2) 
# lines(E_half$r, sqrt(E_half$obs/pi) -E_half$r, col = jet_colors[4], lty = 4, lwd = 2)

text(9, sqrt(max(E$obs)/pi) - 10.1, "VES13", col = "black", font = 1, cex = 1)
text(9, sqrt(max(E_unthomas$obs)/pi) - 9.8, "Omodifierad", col = jet_colors[3], font = 1, cex = 1)
text(9, sqrt(max(E_thomas$obs)/pi) - 10.2, "Thomas", col = jet_colors[2], font = 1, cex = 1)
text(9, sqrt(max(E_halfway$obs)/pi) - 10.0, "Halvvägs", col = jet_colors[1], font = 1, cex = 1)

legend("topleft", legend = c("VES13", "Omodifierad", "'Thomas", "Halvvägs", "Poisson"),
       col = c('black', jet_colors[3],jet_colors[2],jet_colors[1],'red'), lty = c(1,2,2,2,3)) 



# K funktion plot 
plot(E, (.)   ~ r, main = 'K ~ Olika Thomasprocesser', legend = FALSE, cex.main = 1) # ylim = c(0, 6) 
lines(E_unthomas$r, (E_unthomas$obs) , col = jet_colors[3], lty = 1, lwd = 2) 
lines(E_thomas$r, (E_thomas$obs) , col = jet_colors[2], lty = 2, lwd = 2) 
lines(E_halfway$r, (E_halfway$obs) , col = jet_colors[1], lty = 2, lwd = 2) 
lines(E$r, E$theo, col='red', lty=3)
# lines(E_half$r, sqrt(E_half$obs/pi) -E_half$r, col = jet_colors[4], lty = 4, lwd = 2)

text(9, (max(E$obs)) - 10.1, "VES13", col = "black", font = 1, cex = 1)
text(9, (max(E_unthomas$obs)) - 9.8, "Omodifierad", col = jet_colors[3], font = 1, cex = 1)
text(9, (max(E_thomas$obs)) - 10.2, "Thomas", col = jet_colors[2], font = 1, cex = 1)
text(9, (max(E_halfway$obs)) - 10.0, "Halvvägs", col = jet_colors[1], font = 1, cex = 1)
text(9, (max(E$theo)) - 110, "Poisson", col = 'red', font = 1, cex = 1)

legend("topleft", legend = c("VES13", "Omodifierad", "'Thomas", "Halvvägs", "Poisson"),
       col = c('black', jet_colors[3],jet_colors[2],jet_colors[1],'red'), lty = c(1,2,2,2,3)) 






# K och L jämförelse 
par(mfrow=c(1,2))
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1) + 0.1, oma = c(0, 0, 2, 0))
plot(E, (.)   ~ r, main = 'K ~ Olika Thomasprocesser', legend = FALSE, cex.main = 1,
     ylab = expression(K(r) )) # ylim = c(0, 6) 
lines(E_unthomas$r, (E_unthomas$obs) , col = jet_colors[3], lty = 1, lwd = 2) 
lines(E_thomas$r, (E_thomas$obs) , col = jet_colors[2], lty = 2, lwd = 2) 
lines(E_halfway$r, (E_halfway$obs) , col = jet_colors[1], lty = 2, lwd = 2) 
lines(E$r, E$theo, col='red', lty=3)
# lines(E_half$r, sqrt(E_half$obs/pi) -E_half$r, col = jet_colors[4], lty = 4, lwd = 2)

text(9, (max(E$obs)) - 10.1, "VES13", col = "black", font = 1, cex = 1)
text(9, (max(E_unthomas$obs)) - 9.8, "Omodifierad", col = jet_colors[3], font = 1, cex = 1)
text(9, (max(E_thomas$obs)) - 10.2, "Thomas", col = jet_colors[2], font = 1, cex = 1)
text(9, (max(E_halfway$obs)) - 10.0, "Halvvägs", col = jet_colors[1], font = 1, cex = 1)
text(9, (max(E$theo)) - 110, "Poisson", col = 'red', font = 1, cex = 1)

legend("topleft", legend = c("VES13", "Omodifierad", "'Thomas", "Halvvägs", "Poisson"),
       col = c('black', jet_colors[3],jet_colors[2],jet_colors[1],'red'), lty = c(1,2,2,2,3),
       cex = 0.7) 


plot(E, sqrt(./pi) -r  ~ r, main = 'L ~ Olika Thomasprocesser', legend = FALSE, ylim = c(0, 6), cex.main = 1,
     ylab = expression(L(r) -r ))
lines(E_unthomas$r, sqrt(E_unthomas$obs/pi) -E_unthomas$r , col = jet_colors[3], lty = 1, lwd = 2) 
lines(E_thomas$r, sqrt(E_thomas$obs/pi) -E_thomas$r , col = jet_colors[2], lty = 2, lwd = 2) 
lines(E_halfway$r, sqrt(E_halfway$obs/pi) -E_halfway$r , col = jet_colors[1], lty = 2, lwd = 2) 
# lines(E_half$r, sqrt(E_half$obs/pi) -E_half$r, col = jet_colors[4], lty = 4, lwd = 2)

text(9, sqrt(max(E$obs)/pi) - 10.1, "VES13", col = "black", font = 1, cex = 1)
text(9, sqrt(max(E_unthomas$obs)/pi) - 9.8, "Omodifierad", col = jet_colors[3], font = 1, cex = 1)
text(9, sqrt(max(E_thomas$obs)/pi) - 10.2, "Thomas", col = jet_colors[2], font = 1, cex = 1)
text(9, sqrt(max(E_halfway$obs)/pi) - 10.0, "Halvvägs", col = jet_colors[1], font = 1, cex = 1)

legend("topleft", legend = c("VES13", "Omodifierad", "'Thomas", "Halvvägs", "Poisson"),
       col = c('black', jet_colors[3],jet_colors[2],jet_colors[1],'red'), lty = c(1,2,2,2,3),
       cex=0.7) 



# Konfidenshölje för Thomasprocessen ----
simulations1 <- replicate(100, {
  # repeat the block above and return simulated_ppp
  Lest(simulate(kppm(X~1, "Thomas"), nsim=1)[[1]], Lest, correction='iso')$iso - E$r
}, simplify = TRUE)

# Beräkna 2.5:e och 97.5:e percentilen för varje punkt
lower <- apply(simulations1, 1, quantile, probs = 0.025)
upper <- apply(simulations1, 1, quantile, probs = 0.975)
mean_sim <- apply(simulations1, 1, mean)

# Plot mean curve
plot(E$r, mean_sim, type = "l", col = "blue", lwd = 2, lty=2,
     ylim = range(lower, upper),
     ylab = expression(L(r)-r), xlab = expression(italic(r)),
     main = "Konfidenshölje för Thomasprocessen")

polygon(c(E$r, rev(E$r)), c(upper, rev(lower)),
        col = rgb(0.8, 0.8, 1, 0.5), border = NA)

# Lägg till medelvärdeskurvan ovanpå igen (för tydlighet)
lines(E$r, mean_sim, col = "blue", lwd = 2, lty=2)
lines(E$r, Lest(X, Lest, correction='isotropic')$iso - E$r, lwd=2)

legend("topleft", legend = c('VES13-data', 'Thomasprocessen'),
       col = c('black','blue'), lty = c(1,2 ), lwd=c(2,2)) 


simulations_lambda <- replicate(100, {
  # repeat the block above and return simulated_ppp
  
}, simplify = TRUE)




# Lambda konfidenshölje ----
par(mfrow = c(2, 3), mar = c(2, 2, 2, 1), oma = c(0, 0, 4, 0))  # oma reserves space for outer margin title

results <- replicate(100, {
  X_rlgcp <- rLGCP("gauss", m, var = 0.6, scale = 10, win = win, saveLambda = TRUE)
  
  lambda_raw <- attr(X_rlgcp, "Lambda")
  lambda_field <- (lambda_raw - min(lambda_raw)) / (max(lambda_raw) - min(lambda_raw)) * 2.5
  # added 
  win <- Window(lambda_field)
  n_f <- 160
  f_punkt <- data.frame(
    x = runif(n_f, min = win$xrange[1], max = win$xrange[2]),
    y = runif(n_f, min = win$yrange[1], max = win$yrange[2])
  )
  densitet <- as.function(lambda_field)
  
  
  f_punkt$intensity <- densitet(f_punkt$x, f_punkt$y)
  # hist(f_punkt$intensity, main = "Histogram predicted intensiteter",
  #      xlab = "Predicted intensiteter", col = "lightblue", breaks = 20)
  
  t <- 0.3 # se slides
  A <- 15  
  
  f_punkt$lambda <- ifelse(f_punkt$intensity <= t,
                           0,
                           A * (f_punkt$intensity - t))
  f_punkt$n_children <- rpois(n_f, lambda = f_punkt$lambda)
  # hist(f_punkt$lambda, main = "Histogram lambda",
  #      xlab = "lambda(A*(Intensitet - t))", col = "lightgreen", breaks = 20)
  
  barn_std <- 0.8 # varians barn
  
  barn <- do.call(rbind, lapply(1:n_f, function(i) {
    nchild <- f_punkt$n_children[i]
    if (nchild > 0) {
      data.frame(
        x = rnorm(nchild, mean = f_punkt$x[i], sd = barn_std),
        y = rnorm(nchild, mean = f_punkt$y[i], sd = barn_std),
        parent = i
      )
    } else {
      NULL
    }
  }))
  
  # plot(lambda_field, main  = "")
  # points(f_punkt$x, f_punkt$y, col = "black", pch = 1)
  # if (!is.null(barn)) {
  #   points(barn$x, barn$y, col = "black", pch = 1)
  # }
  
  x2 <- range(f_punkt$x, if (!is.null(barn)) barn$x)
  y2 <- range(f_punkt$y, if (!is.null(barn)) barn$y)
  plot(NA, xlim = x2, ylim = y2, type = "n", axes = FALSE, xlab = "", ylab = "", frame.plot = FALSE)
  points(f_punkt$x, f_punkt$y, col = "black", pch = 1)
  points(barn$x, barn$y, col = "black", pch = 1)
  
  Xlambda <- ppp(x=barn$x, y=barn$y, window=square(40))
  
  savelest <- Lest(Xlambda, Lest, correction='iso')$iso - E$r
  # savelest
  # savelest_matrix <- do.call(cbind, savelest_list)
  # savelest_matrix
  savelest
  # list(Xlambda = savelest)
}, simplify = TRUE)

mtext("Simulerade lambdafält punktprocesser (6 varianter)", 
      outer = TRUE, cex = 1.5)

par(
  mfrow=c(1,1),     # One plot per window (default is no grid)
  mar=c(5,4,4,2),   # Margins: bottom, left, top, right (default values)
  oma=c(0,0,0,0),   # Outer margins: top, right, bottom, left (default is none)
  cex=1,            # Character expansion (default text size)
  las=0,            # Axis label orientation (default is horizontal)
  bty="o",          # Box type (default is "o" - open box)
  xpd=FALSE,        # No clipping to the plot region
  tcl=-0.5,         # Tick mark length (default)
  pch=1,            # Default plotting character (empty circle)
  lwd=1,            # Line width (default is 1)
  col="black"       # Default color (black)
)

# Beräkna 2.5:e och 97.5:e percentilen för varje punkt
lower_lam <- apply(results, 1, quantile, probs = 0.025)
upper_lam <- apply(results, 1, quantile, probs = 0.975)
mean_sim_lam <- apply(results, 1, mean)

# Plot mean curve
plot(E$r, mean_sim_lam, type = "l", col = "blue", lwd = 2, lty=2,
     ylim = range(lower_lam, upper_lam),
     ylab = expression(L(r)-r), xlab = expression(italic(r)),
     main = "Konfidenshölje för lambdafält")

polygon(c(E$r, rev(E$r)), c(upper_lam, rev(lower_lam)),
        col = rgb(0.8, 0.8, 1, 0.5), border = NA)

# Lägg till medelvärdeskurvan ovanpå igen (för tydlighet)
lines(E$r, mean_sim_lam, col = "blue", lwd = 2, lty=2)
lines(E$r, Lest(X, Lest, correction='isotropic')$iso - E$r, lwd=2)

legend("topleft", legend = c('VES13-data', 'Lambdafält'),
       col = c('black','blue'), lty = c(1,2 ), lwd=c(2,2)) 



# maybe later 
plot(env_iter, main = 'Envelop of iterativ, blå = data ', legend=FALSE) # 
lines(E_10$r, sqrt(E_10$obs/pi), col = jet_colors[2], lty = 1, lwd = 1)
lines(K_obs$r, sqrt(K_obs$iso/pi), col = "black", lwd = 1)  # data under lest 

legend("bottomright", legend = c("data", "Iterativ"),
       col = c('blue', 'black'), 
       lwd = c(1,1),
       lty=c(1,1))

env_iter




