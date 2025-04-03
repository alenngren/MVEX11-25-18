# setwd("~/Documents/School/Kandidatarbete")
library("spatstat")
#vignette('getstart')
#ERIKA <- readRDS("ERIKA.rds")

filename <- "/Users/alexander/Chalmers/MVEX11-25-18/ERIKA.rds"
ERIKA <- readRDS(filename)

#extracting VES13 data from whole ERIKA data
xs

#plotting VES13 data
plot(V13_x,V13_y)

# Extract unique x and y values to get rid of duplicates
unique_points <- unique(data.frame(x = V13_x, y = V13_y))
V13_xU <- unique_points$x
V13_yU <- unique_points$y

#Making data into a point pattern
window <- owin(xrange = c(0, 40), yrange = c(0, 40))
V13 <- ppp(x = V13_xU, y = V13_yU, window = window)

#####
#Clark Evans Test
#####
#Clark-Evans test for CSR (based on shortest distances)
clarkevans(V13)
clarkevans.test(V13, correction = "donnelly", alternative = "clustered")
# result: R = 0.45624, p-value < 2.2e-16

#####
# Thomas Test
#####
#Thomas process under assumption of homogenuity
fitT <- kppm(V13~1, "Thomas") # ~1 is the predictor in a model for log intensity
fitT # a fitted model of class kppm
#kappa       scale 
# 0.003943241 2.097576813 
plot(simulate(fitT,nsim=6), main = "Simulated Thomas Process") # Den är för klustrad! 


#####
# Modified-Thomas Test
#####
