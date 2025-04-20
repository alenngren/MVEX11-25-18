library(spatstat)
win <- owin(xrange = c(0, 40), yrange = c(0, 40))
m <- as.im(function(x, y) {4 - 1.5 * ((x / 40) - 0.5)^2 + 2 * ((y / 40) - 0.5)^2}, W = win)
X <- rLGCP("gauss", m, var = 0.6, scale = 10, win = win, saveLambda = TRUE)
lambda_field <- (attr(X, "Lambda") - min(attr(X, "Lambda"))) / (max(attr(X, "Lambda")) - min(attr(X, "Lambda"))) * 2.5
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

