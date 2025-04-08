# Load necessary libraries
library(spatstat)
library(spatstat.core)
library(spatstat.geom)

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
