# Load necessary libraries
library(spatstat)
library(spatstat.core)
library(spatstat.geom)
library(mgcv)

# Generate a random spatial point pattern in a unit square
set.seed(123)
pp <- rpoispp(lambda = 100)  # Poisson process with intensity 100

# Compute Ripley's K-function with simulation envelopes
K_env <- envelope(pp, Lest, nsim = 999)  # 99 simulations for confidence envelope

# Plot the envelope
plot(K_env, main = "Envelope Plot for Ripley's K-Function")




# shade 

# mock data: random walk plus a sinus curve.
# two envelopes for added contrast.
install.packages("reshape2")
library(reshape2)
tt=10*sin(c(1:100)/(3*pi))
rr=apply(matrix(rnorm(5000),100,50),2,cumsum) +tt
rr2=apply(matrix(rnorm(5000),100,50),2,cumsum)/1.5 +tt

# stuff data into a dataframe and melt it.
df=data.frame(c(1:100),cbind(rr,rr2) )
names(df)=c("step",paste("ser",c(1:100),sep=""))
dfm=melt(df,id.vars = 1)

# ensemble average
ensemble_av=data.frame(step=df[,1],ensav=apply(df[,-1],1,mean))
ensemble_av$variable=as.factor("Mean")


ggplot(dfm,aes(step,value,group=variable))+
  stat_binhex(alpha=0.2) + geom_line(alpha=0.2) +
  geom_line(data=ensemble_av,aes(step,ensav,size=2))+
  theme(legend.position="none") 

qs = data.frame(
  do.call(
    rbind,
    tapply(
      dfm$value, dfm$step, function(i){quantile(i)})),
  t=1:100)

head(qs)

ggplot() + 
  geom_ribbon(data=qs, aes(x=t, ymin=X0., ymax=X100.),fill="gray30", alpha=0.2) +
  geom_ribbon(data=qs, aes(x=t, ymin=X25., ymax=X75.),fill="gray30", alpha=0.2)

# Different 3
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

enve <- envelope(X, Kest, nsim=99)

plot(enve$r, enve$obs,  type = "l", lty = 1)
# plot(x, y, type = "l", lty = 1)
upper <- enve$hi-enve$theo+enve$obs
lines(enve$r, enve$hi-enve$theo+enve$obs)
lines(enve$r, enve$lo-enve$theo+enve$obs)

ggplot() + 
  geom_ribbon(data=qs, aes(x=t, ymin=X0., ymax=X100.),fill="gray30", alpha=0.2) +
  geom_ribbon(data=qs, aes(x=t, ymin=X25., ymax=X75.),fill="gray30", alpha=0.2)

x <- K_obs$r
f_x <- K_obs$iso
plot(x, f_x, xlim=range(x), ylim=range(f_x), xlab="x", ylab="y", 
     main = "noise-less data",pch=16)
lines(x[order(x)], f_x[order(x)], xlim=range(x), ylim=range(f_x), pch=16)



# Create some sample data
x <- 1:10
y1 <- x^2
y2 <- sqrt(x)

# Set up the layout for two plots in a single row
par(mfrow=c(1, 2), mar=c(4, 4, 2, 1)) # mar adjusts the margins

# First plot
plot(x, y1, type='o', col='blue', main='Plot 1: y = x^2', xlab='X', ylab='Y')

# Second plot
plot(x, y2, type='o', col='red', main='Plot 2: y = sqrt(x)', xlab='X', ylab='Y')


# Set up a tight 1x4 plot layout
par(mfrow = c(1, 4),           # 1 row, 4 columns
    mar = c(2, 2, 2, 1),       # margins: bottom, left, top, right
    oma = c(0, 0, 2, 0))       # outer margins for overall title if needed

# Plot the 4 panels
plot(thomas_process, main = "Unmodified", cex.main=0.9)
plot(simulated_ppp, main = "Thomas", cex.main=0.9)
plot(halfway_thomas, main = "Halvvägs", cex.main=0.9)
plot(X, main = "VES13", cex.main=0.9)

# Optional: Add an overall title
# mtext("Thomas Process Plots", outer = TRUE, cex = 1.5)


setEPS()
postscript("motivera_kluster.eps", width=6, height=3.7)

par(mfrow=c(1,3),
    mar=c(2, 2, 2, 1),     # reduced top margin further
    oma=c(0, 0, 0, 0))    

plot(X_p, main="", axes=FALSE)
box()
mtext("Föräldrarpunkter", side=3, line=0.5, cex=0.9)

plot(X, main="", axes=FALSE)
box()
mtext("Dotterpunkter", side=3, line=0.5, cex=0.9)

plot(pp_uniform, main="", axes=FALSE)
box()
mtext("Uniformfördelning", side=3, line=0.5, cex=0.9)

dev.off()


setEPS()
postscript("motivera_kluster.eps", width=6, height=3.7)

par(mfrow=c(1,3), 
    mar=c(1, 1, 2, 1),  # very small inner margins
    oma=c(0, 0, 0, 0))  # no outer margin

# 1. Föräldrarpunkter
plot.new()  # start new empty plot
plot.window(xlim=c(0,1), ylim=c(0,1))  # define plotting area
points(X_p, pch=16, cex=0.5)  # plot the points
mtext("Föräldrarpunkter", side=3, line=0.2, cex=0.9)  # tight title

# 2. Dotterpunkter
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1))
points(X, pch=16, cex=0.5)
mtext("Dotterpunkter", side=3, line=0.2, cex=0.9)

# 3. Uniformfördelning
plot.new()
plot.window(xlim=c(0,1), ylim=c(0,1))
points(pp_uniform, pch=16, cex=0.5)
mtext("Uniformfördelning", side=3, line=0.2, cex=0.9)

dev.off()



setEPS()
postscript("motivera_kluster.eps", width=6, height=3.7)

par(mfrow=c(1,3), 
    mar=c(1, 1, 2, 1),  
    oma=c(0, 0, 0, 0))  

# 1. Föräldrarpunkter
plot.new()  
plot.window(xlim=c(0,40), ylim=c(0,40))  
points(X_p$x, X_p$y, pch=16, cex=0.5)  
mtext("Föräldrarpunkter", side=3, line=0.2, cex=0.9)

# 2. Dotterpunkter
plot.new()
plot.window(xlim=c(0,40), ylim=c(0,40))
points(X$x, X$y, pch=16, cex=0.5)
mtext("Dotterpunkter", side=3, line=0.2, cex=0.9)

# 3. Uniformfördelning
plot.new()
plot.window(xlim=c(0,40), ylim=c(0,40))
points(pp_uniform$x, pp_uniform$y, pch=16, cex=0.5)
mtext("Uniformfördelning", side=3, line=0.2, cex=0.9)

dev.off()



setEPS()
postscript("motivera_kluster.eps", width=6, height=3.7)

par(mfrow=c(1,3), 
    mar=c(1, 1, 2, 1),  
    oma=c(0, 0, 0, 0))  

# 1. Föräldrarpunkter
plot.new()  
plot.window(xlim=c(0,40), ylim=c(0,40), asp=1)  # asp=1 ensures a square plot
points(X_p$x, X_p$y, pch=16, cex=0.5)  
mtext("Föräldrarpunkter", side=3, line=0.2, cex=0.9)

# 2. Dotterpunkter
plot.new()
plot.window(xlim=c(0,40), ylim=c(0,40), asp=1)
points(X$x, X$y, pch=16, cex=0.5)
mtext("Dotterpunkter", side=3, line=0.2, cex=0.9)

# 3. Uniformfördelning
plot.new()
plot.window(xlim=c(0,40), ylim=c(0,40), asp=1)
points(pp_uniform$x, pp_uniform$y, pch=16, cex=0.5)
mtext("Uniformfördelning", side=3, line=0.2, cex=0.9)

dev.off()


setEPS()
postscript("motivera_kluster.eps", width=6, height=3.7)

par(mfrow=c(1,3), 
    mar=c(1, 1, 2, 1),  
    oma=c(0, 0, 0, 0))  

# 1. Föräldrarpunkter
plot.new()
usr <- c(0, 40, 0, 40)
plot.window(xlim=usr[1:2], ylim=usr[3:4])
rect(0, 0, 40, 40)  # optional frame
points(X_p$x, X_p$y, pch=16, cex=0.5)
mtext("Föräldrarpunkter", side=3, line=0.5, cex=0.9)

# 2. Dotterpunkter
plot.new()
plot.window(xlim=usr[1:2], ylim=usr[3:4])
rect(0, 0, 40, 40)
points(X$x, X$y, pch=16, cex=0.5)
mtext("Dotterpunkter", side=3, line=0.5, cex=0.9)

# 3. Uniformfördelning
plot.new()
plot.window(xlim=usr[1:2], ylim=usr[3:4])
rect(0, 0, 40, 40)
points(pp_uniform$x, pp_uniform$y, pch=16, cex=0.5)
mtext("Uniformfördelning", side=3, line=0.5, cex=0.9)

dev.off()



setEPS()
postscript("motivera_kluster.eps", width=6, height=3.7)

par(mfrow=c(1,3), 
    mar=c(1, 1, 1, 1),  # very tight margins
    oma=c(0, 0, 0, 0))  

# 1. Föräldrarpunkter
plot.new()
usr <- c(0, 40, 0, 40)
plot.window(xlim=usr[1:2], ylim=usr[3:4])
rect(0, 0, 40, 40)
points(X_p$x, X_p$y, pch=16, cex=0.5)
text(20, 40.5, "Föräldrarpunkter", cex=0.9, adj=0.5)  # centered, just above the box

# 2. Dotterpunkter
plot.new()
plot.window(xlim=usr[1:2], ylim=usr[3:4])
rect(0, 0, 40, 40)
points(X$x, X$y, pch=16, cex=0.5)
text(20, 40.5, "Dotterpunkter", cex=0.9, adj=0.5)

# 3. Uniformfördelning
plot.new()
plot.window(xlim=usr[1:2], ylim=usr[3:4])
rect(0, 0, 40, 40)
points(pp_uniform$x, pp_uniform$y, pch=16, cex=0.5)
text(20, 40.5, "Uniformfördelning", cex=0.9, adj=0.5)

dev.off()




set.seed(123)

# 100 punkter, 10 simuleringar
x <- 1:100
simulations <- replicate(10, sin(x / 10) + rnorm(100, 0, 0.2))

# Beräkna 2.5:e och 97.5:e percentilen för varje rad
lower <- apply(simulations, 1, quantile, probs = 0.025)
upper <- apply(simulations, 1, quantile, probs = 0.975)
mean_sim <- apply(simulations, 1, mean)

library(ggplot2)

# Gör ett data.frame för plotting
df <- data.frame(
  x = x,
  mean = mean_sim,
  lower = lower,
  upper = upper
)

# Plot
ggplot(df, aes(x = x)) +
  geom_ribbon(aes(ymin = lower, ymax = upper), fill = "lightblue", alpha = 0.5) +
  geom_line(aes(y = mean), color = "blue") +
  labs(title = "Konfidenshölje från simuleringar",
       y = "Simulerat värde")


library(spatstat)

# Simulate a random point pattern in a unit square
pp <- rpoispp(100)
# Calculate the envelope with 10 simulations under CSR
K_env <- envelope(pp, Kest, nsim=100)

# Plot the envelope
plot(K_env, main="Konfidenshölje för K-funktionen")
K_env <- envelope(pp, Kest, nsim=10, rank=1) # rank=1 for approx 95% if nsim=10

L_env <- envelope(pp, Lest, nsim=10)
plot(L_env, main="Konfidenshölje för L-funktionen")


# Bästa kod för att göra konfidensintervall! 
set.seed(123)

# 100 punkter, 10 simuleringar
x <- 1:100
simulations <- replicate(100, sin(x / 10) + rnorm(100, 0, 0.2))

# Beräkna 2.5:e och 97.5:e percentilen för varje punkt
lower <- apply(simulations, 1, quantile, probs = 0.025)
upper <- apply(simulations, 1, quantile, probs = 0.975)
mean_sim <- apply(simulations, 1, mean)

# Plot mean curve
plot(x, mean_sim, type = "l", col = "blue", lwd = 2,
     ylim = range(lower, upper),
     ylab = "Simulerat värde", xlab = "x",
     main = "Konfidenshölje med plot()")

# Lägg till höljet som ett polygon-område
polygon(c(x, rev(x)), c(upper, rev(lower)),
        col = rgb(0.8, 0.8, 1, 0.5), border = NA)

# Lägg till medelvärdeskurvan ovanpå igen (för tydlighet)
lines(x, mean_sim, col = "blue", lwd = 2)


plot( rThomas(0.003953809, 2.089521688, 140.0548, win=square(40))) # most accurate one 
plot( rThomas(0.55375, 2.089521688, 140.0548, win=square(40))) 
plot( rThomas(intensity(X), sqrt(4.609785), 10, win = square(40)))
plot( rThomas(intensity(X), win = square(40)))

fitT <- kppm(X ~ 1, "Thomas")
X_ordinary_thomas <- simulate(fitT,nsim=1)
simulated_ppp <- X_ordinary_thomas[[1]] # happy with this 
plot(simulated_ppp, main='thomas') # nsim = 9 # plot Thomas 

fit <- Lest(simulate(kppm(X~1, "Thomas"), nsim=1)[[1]])
Lest(fit)


# Simulate point pattern with replicate to create envelope ----
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



# Tror jag fick samma kod som lisa
library(spatstat)

# Fitted model
fitT <- kppm(X ~ 1, "Thomas")

# Extract parameters
lambda_parent <- fitT$poissonIntensity  # parent intensity (note: may need 'fitT$intensity' depending on spatstat version)
sigma <- fitT$clustersigma
mu <- fitT$clusterpar

# Good code snipper 
simulations <- replicate(100, {
  # repeat the block above and return simulated_ppp
}, simplify = FALSE)

# Simulation window
W <- as.owin(X)

# Number of parent points
n_parents <- rpois(1, lambda_parent * area.owin(W))

# Parent locations
parents_x <- runif(n_parents, min = W$xrange[1], max = W$xrange[2])
parents_y <- runif(n_parents, min = W$yrange[1], max = W$yrange[2])

# Offspring locations
offspring_x <- c()
offspring_y <- c()

for (i in 1:n_parents) {
  n_offspring <- rpois(1, mu)
  offspring_x <- c(offspring_x, rnorm(n_offspring, mean = parents_x[i], sd = sigma))
  offspring_y <- c(offspring_y, rnorm(n_offspring, mean = parents_y[i], sd = sigma))
}

# Keep only offspring inside window
inside <- inside.owin(offspring_x, offspring_y, W)
offspring_x <- offspring_x[inside]
offspring_y <- offspring_y[inside]

# Convert to point pattern object
simulated_ppp <- ppp(offspring_x, offspring_y, window = W)

# Plot
plot(simulated_ppp)





# experiment to get what want 
results <- replicate(5, {
  X_rlgcp <- rLGCP("gauss", m, var = 0.6, scale = 10, win = win, saveLambda = TRUE)
  
  lambda_raw <- attr(X_rlgcp, "Lambda")
  lambda_field <- (lambda_raw - min(lambda_raw)) / (max(lambda_raw) - min(lambda_raw)) * 2.5
  
  list(X_rlgcp = X_rlgcp, lambda_field = lambda_field)
}, simplify = FALSE)

# Set up a 2x3 plotting window (for 5 plots)
par(mfrow = c(2, 3), mar = c(2, 2, 2, 1))  # 2 rows, 3 columns, small margins

# Loop over the 5 simulations and plot them
for (i in 1:5) {
  plot(results[[i]]$X_rlgcp, 
       main = paste("Iteration", i),
       cols = "black", pch = 20)
}

# Leave the 6th plot empty
plot.new()


# Set up a 2x3 plotting window
par(mfrow = c(2, 3), mar = c(2, 2, 2, 1))  # 2 rows, 3 columns, small margins

# Loop over the 5 simulations and plot them
for (i in 1:5) {
  plot(results[[i]]$X_rlgcp$x, results[[i]]$X_rlgcp$y, 
       main = paste("Iteration", i),
       cols = "black", pch = 20,
       xlim = c(0, 40), ylim = c(0, 40))
}

# Leave the 6th plot empty
plot.new()



# The one that works, start -----
par(mfrow = c(2, 3), mar = c(2, 2, 2, 1))  # 2 rows, 3 columns, small margins
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

results1 <- t(results)

# Beräkna 2.5:e och 97.5:e percentilen för varje punkt
lower_lam <- apply(results, 1, quantile, probs = 0.025)
upper_lam <- apply(results, 1, quantile, probs = 0.975)
mean_sim_lam <- apply(results, 1, mean)

# Plot mean curve
par(mfrow=c(1,1))
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





window <- owin(range(x), range(y))

ppp_data <- ppp(x, y, window = square(40)) # found you exclude 6 points 
densitet <- density(ppp_data, sigma = 0.6) 
kordinater <- as.data.frame(densitet)
colnames(kordinater) <- c("x", "y", "ca")

kordinater$ca <- pmax(kordinater$ca, 0)
# Careful to run this one! (takes long time) -----
results_lgcp <- replicate(100,{
  lgcp <- gam(ca ~ s(x, y, k = 50), family = poisson, data = kordinater)
  kordinater$vänte_lam <- exp(predict(lgcp, newdata = kordinater))
  skal <- 0.6 
  kordinater$skalad_lamda1 <- skal * pmax(kordinater$vänte_lam - 0.05, 0.1)
  matris_skalad <- matrix(kordinater$skalad_lamda1, 
                          nrow = dim(densitet$v)[1], 
                          ncol = dim(densitet$v)[2], 
                          byrow = FALSE)
  skalad_lamda2 <- im(matris_skalad, xcol = densitet$xcol, yrow = densitet$yrow)
  LGCP <- rpoispp(lambda = skalad_lamda2, win = window) # X pattern
  # LGCP # stores x,y and n 
  X_lgcpsim <- ppp(x=LGCP$x, y=LGCP$y, window=square(40))
  Lest(X_lgcpsim, Lest, correction='iso')$iso - E$r

})

# X_lgcpsim <- ppp(x=LGCP$x, y=LGCP$y, window=square(40))
# Lest(X_lgcpsim, Lest, correction='iso')$iso - E$r

# Beräkna 2.5:e och 97.5:e percentilen för varje punkt
lower_lgcp <- apply(results_lgcp, 1, quantile, probs = 0.025)
upper_lgcp <- apply(results_lgcp, 1, quantile, probs = 0.975)
mean_sim_lgcp <- apply(results_lgcp, 1, mean)

# Plot mean curve for LGCP ----
# par(mfrow=c(1,1))
# eventuellt plotta några bredvid varandra, men fint såhär, eller 2,3 plots
plot(E$r, mean_sim_lgcp, type = "l", col = "blue", lwd = 2, lty=2,
     ylim = range(lower_lgcp, upper_lgcp),
     ylab = expression(L(r)-r), xlab = expression(italic(r)),
     main = "Konfidenshölje för LGCP process")

polygon(c(E$r, rev(E$r)), c(upper_lgcp, rev(lower_lgcp)),
        col = rgb(0.8, 0.8, 1, 0.5), border = NA)

# Lägg till medelvärdeskurvan ovanpå igen (för tydlighet)
lines(E$r, mean_sim_lgcp, col = "blue", lwd = 2, lty=2)
lines(E$r, Lest(X, Lest, correction='isotropic')$iso - E$r, lwd=2)

legend("topleft", legend = c('VES13-data', 'LGCP process'),
       col = c('black','blue'), lty = c(1,2 ), lwd=c(2,2)) 


# G funktion konfidensintervall ----
simulations_f_thom <- replicate(10, {
  # repeat the block above and return simulated_ppp
  fs <- Fest(simulate(kppm(X~1, "Thomas"), nsim=1)[[1]], correction='rs')$rs
  fs
}, simplify = TRUE)

# Beräkna 2.5:e och 97.5:e percentilen för varje punkt
lower_f_thom <- apply(simulations_f_thom, 1, quantile, probs = 0.025)
upper_f_thom <- apply(simulations_f_thom, 1, quantile, probs = 0.975)
mean_sim_f_thom <- apply(simulations_f_thom, 1, mean)

# Plot mean curve
plot(E$r, mean_sim_f_thom, type = "l", col = "blue", lwd = 2, lty=2,
     ylim = range(lower_f_thom, upper_f_thom),
     ylab = expression(G(r)), xlab = expression(italic(r)),
     main = "Konfidenshölje för Thomasprocessen")

polygon(c(E$r, rev(E$r)), c(upper_g_thom, rev(lower_f_thom)),
        col = rgb(0.8, 0.8, 1, 0.5), border = NA)

# Lägg till medelvärdeskurvan ovanpå igen (för tydlighet)
lines(E$r, mean_sim_g_thom, col = "blue", lwd = 2, lty=2)
lines(E$r, Gest(X, Gest, correction='isotropic')$iso - E$r, lwd=2)

legend("topleft", legend = c('VES13-data', 'Thomasprocessen'),
       col = c('black','blue'), lty = c(1,2 ), lwd=c(2,2)) 

fs <- Fest(X, correction='rs')$rs
plot(fs, type='l')

n=20
fs_probs <- rep(0.5, n)
fs <- Fest(X, correction='rs', probs = fs_probs)$rs
fs_thom <- Fest(simulate(kppm(X~1, "Thomas"), nsim=1)[[1]], correction='rs', probs = fs_probs)$rs
length(fs)
length(fs_thom)
