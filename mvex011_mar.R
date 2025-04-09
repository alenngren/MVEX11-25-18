# library(feather)
library(spatstat)
# install.packages("feather")
install.packages("spatstat")

filename <- "/Users/alexander/Chalmers/MVEX11-25-18/ERIKA.rds"
path <- "/Users/alexander/Chalmers/MVEX11-25-18/python-data/ves13-thomasprocess-modified.feather"
path_csv <- "/Users/alexander/Chalmers/MVEX11-25-18/python-data/ves13-thomasprocess-modified.csv"
file_path <- "/Users/alexander/Chalmers/MVEX11-25-18/python-data/"

ERIKA <- readRDS(filename)
modelves13 <- read.csv(path_csv)
x_c <- modelves13$x
y_c <- modelves13$y

x_evo03 <- ERIKA$EVO03$small$x
y_evo03 <- ERIKA$EVO03$small$y
X_evo03 <- ppp(x=x_evo03, y=y_evo03, window=square(40))

x <- ERIKA$VES13$small$x
y <- ERIKA$VES13$small$y
X <- ppp(x=x, y=y, window=square(40))


plot(X, main='data')
# plot(K, main="K function for VES13")
plot(K, sqrt(./pi) -r ~ r, ylab="L(r) - r ", main="L(r) - r ~ VES13")
# fitT <- kppm(V13~1, "Thomas")

X_csv <- ppp(x=x_c, y=y_c, window=square(40))
plot(X_csv, main="VES13 model 1")
K_csv <- Kest(X_csv, correction="Ripley")
plot(K_csv, sqrt(./pi) -r ~ r, ylab="L(r) - r ", main="L(r) - r ~ VES13 Model 1")

# från lisa 
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

#Creating a ppp for the process and plotting
thomas_process <- ppp(offspring_x, offspring_y, window=window) 
plot(thomas_process, main="Simulated Unmodified Thomas Process")
points(parents, col="red", pch=3)  # Parent points in red

K_thomas <- Kest(thomas_process, correction="Ripley")
plot(K_thomas, sqrt(./pi) -r ~ r, ylab="L(r) - r ", main="L-r, un-modified Thomas")

fitUT <- kppm(thomas_process ~ 1, "Thomas")
fitUT

par(mfrow=c(1,1))

#par(mfrow=c(3,1))
plot(K_csv, sqrt(./pi) -r ~ r, ylab="L(r) - r ", main="L(r) - r ~ VES13 Model 1")
plot(K, sqrt(./pi) -r ~ r, ylab="L(r) - r ", main="L(r) - r ~ VES13")
plot(K_thomas, sqrt(./pi) -r ~ r, ylab="L(r) - r ", main="L-r, un-modified Thomas")

K_csv
saveRDS(K_csv, paste(file_path,"K_csv.rds"))
write.csv(K_thomas, paste(file_path,"K_thomas.csv"))
write.csv(K_csv, paste(file_path,"K_csv.csv"))
write.csv(K, paste(file_path,"K.csv"))

write.csv(ERIKA, paste(file_path,"ERIKA.csv"))



## sätt ihop ny bild 
#Modified-Thomas Test : Halfway point between modified and normal Thomas
#####
#define new parameters
lambda_avg <- 2 / (1/lambda_p_hat + 1/0.003943241)  # Harmonic mean
mu_avg <- sqrt(mu_hat * 139.7962)  # Geometric mean
sigma_avg <- sqrt(sigma_hat * 2.097576813)  # Geometric mean

num_parents_avg <- rpois(1, lambda_avg * area.owin(window))
parents_avg <- runifpoint(num_parents_avg, win = window) #Simulate parents using average intensity

#Generate offspring
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
halfway_thomas <- ppp(offspring_x, offspring_y, window=window)
plot(halfway_thomas, main="Hybrid Thomas Process (Averaged Parameters)")
points(parents_avg, col="red", pch=3)  # Parent points in red

K_halfthomas <- Kest(halfway_thomas, correction="Ripley")
plot(K_halfthomas, sqrt(./pi) -r ~ r, ylab="L(r) - r ", main="L-r, halfway Thomas")


# Combine the plot
plot(K, sqrt(./pi) -r ~ r, ylab="L(r) - r ", main="L -r ~ VES13, Thomas(half), Thomas(mod), Thomas(un) ")
plot(K_csv, sqrt(./pi) -r ~ r, ylab="L(r) - r ", main="L(r) - r ~ VES13 Model 1", add=TRUE, col='blue')
plot(K_thomas, sqrt(./pi) -r ~ r, ylab="L(r) - r ", main="L-r, un-modified Thomas", add=TRUE, col='red')
plot(K_halfthomas, sqrt(./pi) -r ~ r, ylab="L(r) - r ", main="L-r, halfway Thomas", add=TRUE, col='green')

# save ppp half thomas  
data_file <- "/Users/alexander/Downloads/"

write.csv(halfway_thomas, paste(data_file,"half_thomas_points.csv"))
## Plot för k function istället

plot(K, main="K function for VES13")
plot(K_csv, main="K function for VES13", add=TRUE, col='blue')
plot(K_thomas, main="K function for VES13", add=TRUE, col='red')
plot(K_halfthomas, main="K function for VES13", add=TRUE, col='green')
# plot(Kest(X))
# plot(Kest(halfway_thomas))
swp <- rescale(X)
Kvb <- varblock(swp, Kest, nx=3, ny=3)
Kloh <- lohboot(swp, Kest)
xyV <- plot(Kvb, limitsonly=TRUE)
xlim <- xyV$xlim
xyL <- plot(Kloh, limitsonly=TRUE, xlim=xlim)
ylim <- range(xyV$ylim, xyL$ylim)

newplot(12, 0.95)
setmargins(0.5+c(3,3,0,1))

par(mfrow=c(1,2))
plot(Kvb, main="", legend=FALSE, xlim=xlim, ylim=ylim,lty=c(1,2,1,1))
plot(Kloh, main="", legend=FALSE, xlim=xlim, ylim=ylim)
par(mfrow=c(1,1))

swp <- rescale(X)
Kloh <- lohboot(Kest(X))

Kloh_half <- lohboot(Kest(halfway_thomas))

Kloh_thomas <- lohboot(Kest(thomas_process))

Kloh_csv <- lohboot(Kest(X_csv))

plot(Kloh)
plot(Kloh_half, add=TRUE)
plot(Kloh_csv, add=TRUE)
plot(Kloh_thomas, add=TRUE)
X
thomas_process
halfway_thomas          
X_csv

swp <- rescale(X)

E <- envelope(X, Kest, nsim=39, verbose=FALSE)
E_thomas <- envelope(thomas_process, Kest, nsim=39)
E_half <- envelope(halfway_thomas, Kest, nsim=39)
E_csv <- envelope(X_csv, Kest, nsim=39)

plot(E, main='Envelope (CSR) Kest ~ data, half thomas, mod thomas, thomas')
plot(E_thomas, add=TRUE)
plot(E_half, add=TRUE)
plot(E_csv, add=TRUE)

E <- envelope(X, Lest, nsim=99, verbose=FALSE)
E_thomas <- envelope(thomas_process, Lest, nsim=99)
E_half <- envelope(halfway_thomas, Lest, nsim=99)
E_csv <- envelope(X_csv, Lest, nsim=99)

plot(E, main='Envelope Lest ~ data, half thomas, mod thomas, thomas')
plot(E_thomas, add=TRUE)
plot(E_half, add=TRUE)
plot(E_csv, add=TRUE)

ce_x <- clarkevans.test(X, alternative="clustered")
ce_tp <- clarkevans.test(thomas_process, alternative="clustered")
ce_csv <- clarkevans.test(X_csv, alternative="clustered")
ce_half <- clarkevans.test(halfway_thomas, alternative="clustered")
ce_x$statistic

ce <- c(ce_x, ce_tp, ce_csv, ce_half)
ce$data.name
# type = c('VES13', 'Thomasproces', 'Modifierad', 'Halfway')

var1 <- 1:3
var2 <- 1:3



df <- data.frame(
  score = c(ce_x$statistic, ce_tp$statistic, ce_csv$statistic, ce_half$statistic),
  type = c('VES13', 'Thomasproces', 'Modifierad Thomas', 'Halfway Thomas')
  
)
                


ce_bar <- barplot(height=df$score, names=type, 
        col=rgb(0.8,0.1,0.1,0.6),
        xlab="punktmönster", 
        ylab="R värde", 
        main="Clark-Evans test", 
        ylim=c(0,1)
)
text(ce_bar, df$score, pos=1 , paste("R= ", round(df$score,2), sep="") ,cex=1.2) 

install.packages("spatstat.core")
library(spatstat.core)


fit(halfway_thomas, ~1, Thom). # doesnt work 
envelope(k)



k = Kest(X, nsim=1000)

envelope(k)

E <- envelope(X, Lest, nsim=1000)
E$lo
E$hi
plot(E)


# Mats plot 
install.packages("mgcv")
library(mgcv)
library(spatstat)
# data <- read.csv("C:/Users/matti/Downloads/VES13_small.csv")
filename <- "/Users/alexander/Chalmers/MVEX11-25-18/ERIKA.rds"
ERIKA <- readRDS(filename)
x <- ERIKA$VES13$small$x
y <- ERIKA$VES13$small$y


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

# this is the ppp pattern. Plot L med mats data 
plot(LGCP)

E_lgcp <- envelope(LGCP, Kest, nsim=999, verbose=TRUE)
E_lgcp_evo013 <- envelope(LGCP, Kest, nsim=999, verbose=TRUE)

E <- envelope(X, Kest, nsim=999, verbose=TRUE)
plot(E, sqrt(./pi) -r  ~ r, main = 'L ~ med iterativ', legend = FALSE) 
lines(E_lgcp$r, sqrt(E_lgcp$obs/pi) -E_lgcp$r , col = jet_colors[1], lty = 2, lwd = 2)


par(mfrow=c(1,3))
plot(ppp_data, main = "Original Data")
plot(densitet, main = "Latent skattat fält")
plot(LGCP, main = "LGCP från skattad fält")

all_ves13 <- allstats(ppp_data)
all_lgcp <- allstats(LGCP)
plot(all_ves13)
plot(all_lgcp)


E <- envelope(ppp_data, Kest, correction="border", nsim=999)
# E <- envelope(X, Lest, nsim=99, verbose=FALSE)
E_thomas <- envelope(thomas_process, Kest, nsim=9999)
E_half <- envelope(halfway_thomas, Kest, nsim=9999)
E_csv <- envelope(X_csv, Kest, nsim=9999)


plot(E, sqrt(./pi) ~ r, main = 'K function envelope')
plot(E_thomas, sqrt(./pi) ~ r, add=TRUE)
plot(E_half, sqrt(./pi) ~ r, add=TRUE)
plot(E_csv, sqrt(./pi) ~ r, add=TRUE)

E <- envelope(X, Kest, nsim=9999, verbose=FALSE)
E_thomas <- envelope(thomas_process, Kest, nsim=9999)
E_half <- envelope(halfway_thomas, Kest, nsim=9999)
E_csv <- envelope(X_csv, Kest, nsim=9999)

#install.packages("fields")  # Install if not already installed
library(fields)

jet_colors <- tim.colors(5)  # Generate 100 jet-like colors
image(matrix(1:100, 10, 10), col = jet_colors)  # Example visualization

plot(E, sqrt(./pi) -r  ~ r, main = 'L funktion', legend = FALSE) 
lines(E_thomas$r, sqrt(E_thomas$obs/pi) -E_thomas$r , col = jet_colors[1], lty = 2, lwd = 2)  # Theoretical function
lines(E_half$r, sqrt(E_half$obs/pi) -E_half$r, col = jet_colors[2], lty = 4, lwd = 2)  # Theoretical function
lines(E_csv$r, sqrt(E_csv$obs/pi) -E_csv$r, col = jet_colors[4], lty = 2, lwd = 2)  # Theoretical function
legend("topleft", legend = c("data", "theo", 'E_Thomas', 'E_half', 'E_csv'), 
       col = c(1, 2, jet_colors[1], jet_colors[2], jet_colors[4]), lty = c(1, 2, 2,2,2))

text(9, sqrt(max(E$obs)/pi) - 10.2, "data", col = "black", font = 1, cex = 1)
text(9, sqrt(max(E_thomas$obs)/pi) - 10.2, "E_thomas", col = jet_colors[1], font = 1, cex = 1)
text(9, sqrt(max(E_half$obs)/pi) - 10.2, "E_half", col = jet_colors[2], font = 1, cex = 1)
text(9, sqrt(max(E_csv$obs)/pi) - 10.2, "E_csv", col = jet_colors[4], font = 1, cex = 1)
# updated the colors and plots 



