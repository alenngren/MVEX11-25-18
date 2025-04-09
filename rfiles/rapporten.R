# Load data

library(readr)
library(spatstat)
iterpath <- '/Users/alexander/Chalmers/MVEX11-25-18/dataframes/unif-iter10.csv'
path_csv <- "/Users/alexander/Chalmers/MVEX11-25-18/python-data/ves13-thomasprocess-modified.csv"


iterframe <- read_csv(iterpath, col_names = c("x", "y"))
plot(iterframe)

filename <- "/Users/alexander/Chalmers/MVEX11-25-18/ERIKA.rds"
ERIKA <- readRDS(filename)
x <- ERIKA$VES13$small$x
y <- ERIKA$VES13$small$y
X <- ppp(x=x, y=y, window=square(40))

# Plot observed K-function > L-function START HERE! WORKS 
K_obs <- Kest(X)

# Poisson 
en <- envelope(X, nsim=99)
plot(en, main="VES13 - Envelope of K function with CSR")

plot( rThomas(intensity(X), 22, 12))

fitT <- kppm(X ~ 1, "Thomas")
X_ordinary_thomas <- simulate(fitT,nsim=1)
simulated_ppp <- X_ordinary_thomas[[1]] # happy with this 
points_matrix <- coords(simulated_ppp)
print(points_matrix)
plot(simulate(fitT,nsim=1)) # nsim = 9
fitL <- kppm(X ~ 1, "LGCP")
plot(simulate(fitL,nsim=9))

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
sigma_avg <- sqrt(sigma_hat * 2.097576813)  # Geometric mean

num_parents_avg <- rpois(1, lambda_avg * area.owin(window))
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

plot(E, sqrt(./pi) -r  ~ r, main = 'L ~ iterativ steg 1,6,11', legend = FALSE) 
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
# dev.off()






# Clark Evans statistics 
clarkevans.test(X_iter, alternative="clustered")
clarkevans.test(X_iter0, alternative="clustered")
clarkevans.test(X_iter5, alternative="clustered")
clarkevans.test(X_iter10, alternative="clustered")
       

# LGPC modell 
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











# F function, empty space statistics

Fest(X, ..., eps, r=NULL, breaks=NULL,
     correction=c("rs", "km", "cs"),
     domain=NULL)

plot(Fest(X))
X_f <- Fest(X, correction = 'rs')
Xthomas_f <- Fest(simulated_ppp, correction = 'rs') 
Xhalf_f <- Fest(halfway_thomas, correction = 'rs')
Xiter_f <- Fest(X_iter10, correction = 'rs')
Xlgcp_f <- Fest(LGCP)


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



# fit new 
fit_iter <- kppm(X_iter10 ~ 1, "Thomas")
fit <- ppm(X ~ 1, "Thomas") # doesnt work 

env_iter <- envelope(fit_iter, Lest, nsim=99)

lgcp_ppp <- as.ppp(lgcp)
k_lgcp <- Kest(lgcp)

plot(env_iter - E$r, main = 'Envelop of iterativ, blå = data ', legend=FALSE)
lines(E_10$r, sqrt(E_10$obs/pi)- E$r, col = jet_colors[2], lty = 1, lwd = 1)
lines(K_obs$r, sqrt(K_obs$iso/pi)- E$r, col = "black", lwd = 1)  # data under lest 

legend("bottomright", legend = c("data", "Iterativ"),
       col = c('blue', 'black'), 
       lwd = c(1,1),
       lty=c(1,1))

plot(env_iter, main = 'Envelop of iterativ, blå = data ', legend=FALSE) # 
lines(E_10$r, sqrt(E_10$obs/pi), col = jet_colors[2], lty = 1, lwd = 1)
lines(K_obs$r, sqrt(K_obs$iso/pi), col = "black", lwd = 1)  # data under lest 

legend("bottomright", legend = c("data", "Iterativ"),
       col = c('blue', 'black'), 
       lwd = c(1,1),
       lty=c(1,1))

env_iter


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