# filename <- file.choose()
filename <- "/Users/alexander/Chalmers/MVEX11-25-18/ERIKA.rds"
ERIKA <- readRDS(filename)


# install.packages("spatstat")
library(spatstat)

x <- ERIKA$VES13$small$x
y <- ERIKA$VES13$small$y
X <- ppp(x=x, y=y, window=square(40))
X
X <- unique(X)
X
plot(X, ylab="m", xlab="m", main="VES13 data seedlings")

K <- Kest(X, correction="Ripley")
# K <- Kest(cells, correction="isotropic")
plot(K, main="K function for VES13")
plot(K, sqrt(./pi) -r ~ r, ylab="L(r) - r ", main="L(r) - r ~ VES13")


online <- interactive()
Nsim <- if(online) 19 else 3
en <- envelope(X, nsim=19)
plot(en, main="VES13 - Envelope of K function with CSR")

data(ERIKA$VES13$large)
ERIKA$VES13$small = unclass(ERIKA$VES13$large)

write.table(ERIKA$VES13$large, file="/Users/alexander/Chalmers/MVEX11-25-18/ERIKA_np.txt", row.names=F, sep=",")


ERIKA$VES13$small$x 
ERIKA$VES13$small$y
ERIKA$VES13$large

sum(y == y)
length(which(dbh==dbh))

K$theo

plot( rThomas(intensity(X), 22, 12))


fitT <- kppm(X ~ 1, "Thomas")
plot(simulate(fitT,nsim=9))
fitL <- kppm(X ~ 1, "LGCP")
plot(simulate(fitL,nsim=9))

fitT # uniform intensity = lambda, kappa = intensity of parents points, scale = std, 
# phi = styrka 

  

## Plot för k function istället 

plot(K, main="K function for VES13")
plot(K, main="K function for VES13")
plot(K, main="K function for VES13")
plot(K, main="K function for VES13")


ves13 <- allstats(X, dataname='VES13')
plot(ves13)

halfthomas <- allstats(halfway_thomas)
plot(halfthomas)
