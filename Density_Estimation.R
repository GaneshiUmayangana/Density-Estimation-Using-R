##############1.	How to decide the number of bins

n <- 200
sig <- sample(c(10, 60), n, replace=T)
x <- rnorm(n, 0, sig)


# rounding the constant in Scott’s rule
# and using sample standard deviation to estimate sigma
h <- 3.5 * sd(x) * n^(-1/3)


# number of classes is determined by the range and h
m <- min(x)
M <- max(x)
nclass <- ceiling((M - m) / h)
scott.breaks <- m + h * 0:nclass


par(mfrow=c(2,2))
hist(x, breaks = "sturges", freq = FALSE, main = "Sturges",col = "red")
hist(x, breaks = "scott", freq = FALSE,main = "Scott",col = "yellow")
hist(x, breaks = "FD", freq = FALSE,
     main = "Freedman-Diaconis",col = "green")
hist(x, breaks = scott.breaks, freq = FALSE,
     main = "Scott rule from basic formula",col = "pink")
par(mfrow=c(1,1))




################Frequency polygon estimation

##############
# Load the mtcars dataset
data(mtcars)

# Display the first few rows of the dataset
head(mtcars)

# Use the mtcars dataset
w <- mtcars$mpg
n <- length(w)

# freq poly bin width using normal ref rule
h <- 2.15 * sqrt(var(w)) * n^(-1/5)

# calculate the sequence of breaks and histogram
br <- pretty(w, diff(range(w)) / h)
brplus <- c(min(br)-h, max(br+h))
histg <- hist(w, breaks = br, freq = FALSE,
              main = "Frequency Polygon Density Estimate of MPG in mtcars",
              xlim = brplus, col = "skyblue", xlab = "Mpg")

vx <- histg$mids # density est at vertices of polygon
vy <- histg$density

delta <- diff(vx)[1] # h after pretty is applied
k <- length(vx)
vx <- vx + delta # the bins on the ends
vx <- c(vx[1] - 2 * delta, vx[1] - delta, vx)
vy <- c(0, vy, 0)

# add the polygon to the histogram
polygon(vx, vy)




##################### ASH density estimation

library(MASS)
data(mtcars)

mpg_data <- mtcars$mpg
n <- length(mpg_data)
m <- 20
a <- min(mpg_data) - .5
b <- max(mpg_data) + .5
h <- 1.5  # Adjusted bandwidth for demonstration purposes
delta <- h / m

# Get the bin counts on the delta-width mesh.
br <- seq(a - delta * m, b + 2 * delta * m, delta)
histg <- hist(mpg_data, breaks = br, plot = FALSE)
nk <- histg$counts

K <- abs((1 - m):(m - 1))
fhat <- function(x) {
  # Locate the leftmost interval containing x
  i <- max(which(x > br))
  k <- (i - m + 1):(i + m - 1)
  # Get the 2m-1 bin counts centered at x
  vk <- nk[k]
  sum((1 - K / m) * vk) / (n * h)  # f.hat
}

# Density can be computed at any points in the range of data
z <- as.matrix(seq(a, b + h, .1))
f.ash <- apply(z, 1, fhat)  # Density estimates at midpoints

# Plot ASH density estimate over histogram
br2 <- seq(a, b + h, h)
hist(mpg_data, breaks = br2, freq = FALSE, main = "",
     ylim = c(0, max(f.ash)), xlab = "Miles Per Gallon (mpg)",col = "skyblue")
lines(z, f.ash, col = "red", lwd = 2)




##################Kernel density estimation

# Load the mtcars dataset
data(mtcars)

# Extract the mpg variable
mpg_data <- mtcars$mpg

# Create a histogram
hist(mpg_data, main = "Kernel Density Estimate of MPG in mtcars",
     col = "skyblue", xlab = "Miles Per Gallon (mpg)", ylab = "Frequency", probability = TRUE)

# Add Kernel Density Estimate
lines(density(mpg_data), col = "blue", lwd = 2)


###### Kernel density with different bandwiths

# Load the mtcars dataset
data(mtcars)

# Extract the mpg variable
mpg_data <- mtcars$mpg

# Set up a layout for multiple plots
par(mfrow = c(2, 2))

# Generate kernel density estimates with different bandwidths
bandwidths <- c(0.25, 0.4, 0.6, 1)

for (bw in bandwidths) {
  kde <- density(mpg_data, bw = bw, kernel = "gaussian")
  
  # Plot Kernel Density Estimate
  plot(kde, main = "", col = "blue", lwd = 2, xlab = "Miles Per Gallon (mpg)", ylab = "Density")
  
  # Add Legend
  legend("topright", legend = paste( bw), col = "red", lwd = 2, bty = "n")
}

# Add a single main title for all plots
mtext("Kernel Density Estimation (Gaussian) with Different Bandwidths", line = -1, outer = TRUE)

# Reset the layout to a single plot
par(mfrow = c(1, 1))

#######################################



# Create an empty plot for the first density curve
kde <- density(mpg_data, bw = 0.25, kernel = "gaussian")
plot(kde, col = "skyblue", lwd = 2, xlab = "Miles Per Gallon (mpg)", ylab = "Density",
     main = "Kernel Density Estimation (Gaussian) with Different Bandwidths")

# Add subsequent density curves with different colors
bandwidths <- c(1, 3, 10)
colors <- c("red", "green", "blue")

for (i in seq_along(bandwidths)) {
  kde <- density(mpg_data, bw = bandwidths[i], kernel = "gaussian")
  lines(kde, col = colors[i], lwd = 2)
}

# Add a legend
legend("topright", legend = c("Bandwidth = 0.25", "Bandwidth = 1", "Bandwidth = 3", "Bandwidth = 10"),
       col = c("skyblue", "red", "green", "blue"), lwd = 2, bty = "n")


########################

# Load the mtcars dataset
data(mtcars)

# Extract the mpg variable
mpg_data <- mtcars$mpg

# Create an empty plot for the first density curve
kde_gaussian <- density(mpg_data, kernel = "gaussian")
plot(kde_gaussian, col = "blue", lwd = 2, xlab = "Miles Per Gallon (mpg)", ylab = "Density",
     main = "Kernel Density Estimation - Different Kernels")

# Add subsequent density curves with different colors
kernels <- c("epanechnikov", "rectangular", "triangular", "biweight")
colors <- c("red", "green", "purple", "orange")

for (i in seq_along(kernels)) {
  kde <- density(mpg_data, kernel = kernels[i])
  lines(kde, col = colors[i], lwd = 2)
}

# Add a legend
legend("topright", legend = c("Gaussian", "Epanechnikov", "Rectangular", "Triangular", "Biweight"),
       col = c("blue", "red", "green", "purple", "orange"), lwd = 2, bty = "n")

# Reset the layout to a single plot
par(mfrow = c(1, 1))

####################
# Load the mtcars dataset
data(mtcars)

# Extract the mpg variable
mpg_data <- mtcars$mpg

# Set up a layout for multiple plots
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))

# Specify kernel types
kernels <- c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight")
colors <- c("blue", "red", "green", "purple", "orange")

# Generate density plots for each kernel
for (i in seq_along(kernels)) {
  kde <- density(mpg_data, kernel = kernels[i])
  
  # Plot Density Estimate
  plot(kde, col = colors[i], lwd = 2, xlab = "Miles Per Gallon (mpg)", ylab = "Density",
       main = paste("Kernel Density Estimation -", kernels[i]))
}

# Reset the layout to a single plot
par(mfrow = c(1, 1))


###########Boundary Kernels
# Exponential density

# Exponential density
set.seed(123)  # for reproducibility
x <- rexp(1000, 1)

# Plot estimated density
plot(density(x), xlim = c(-1, 6), ylim = c(0, 1), col = "blue", main = "Exponential density", xlab = bw, ylab = "Density")

# Add a vertical line at x = 0
abline(v = 0)

# Add the true density to compare
y <- seq(0, 6, 0.01)
lines(y, dexp(y, 1), lty = 2, col = "red")

# Add legend
legend("topright", legend = c("Estimated Density", "True Density"), col = c("blue", "red"), lty = c(1, 2), bty = "n")


#####################3
###Reflected boundary

set.seed(123)  # for reproducibility
x <- rexp(1000, 1)

# Combine with mirrored values
xx <- c(x, -x)

# Estimate density
g <- density(xx, bw = bw.nrd0(x))
a <- seq(0, 6, .01)
ghat <- approx(g$x, g$y, xout = a)
fhat <- 2 * ghat$y  # Density estimate along a

# Plot the estimated density and true density
bw <- paste("Bandwidth = ", round(g$bw, 5))
plot(a, fhat, type="l", xlim=c(-1, 6), ylim=c(0, 1), col = "blue",
     main = "Exponential density with Reflected boundary ", xlab = bw, ylab = "Density")
abline(v = 0, col = "black")

# Add the true density for comparison
y <- seq(0, 6, 0.01)
lines(y, dexp(y, 1), lty = 2, col = "red")

# Add legend
legend("topright", legend = c("Estimated Density", "True Density"), col = c("blue", "red"), lty = c(1, 2), bty = "n")



###########Bivariate Frequency
bin2d <-
  function(x, breaks1 = "Sturges", breaks2 = "Sturges"){
    # Data matrix x is n by 2
    # breaks1, breaks2: any valid breaks for hist function
    # using same defaults as hist
    histg1 <- hist(x[,1], breaks = breaks1, plot = FALSE)
    histg2 <- hist(x[,2], breaks = breaks2, plot = FALSE)
    brx <- histg1$breaks
    bry <- histg2$breaks
    # bin frequencies
    freq <- table(cut(x[,1], brx), cut(x[,2], bry))
    return(list(call = match.call(), freq = freq,
                breaks1 = brx, breaks2 = bry,
                mids1 = histg1$mids, mids2 = histg2$mids))
  }
bin2d(iris[1:50,1:2])

bin2d(x = iris[1:50, 1:2])


#
#generate standard bivariate normal random sample
n <- 2000; d <- 2
x <- matrix(rnorm(n*d), n, d)
# compute the frequency table and density estimates
# using bin2d function from the previous example
b <- bin2d(x)
h1 <- diff(b$breaks1)
h2 <- diff(b$breaks2)
# matrix h contains the areas of the bins in b
h <- outer(h1, h2, "*")
Z <- b$freq / (n * h) # the density estimate
persp(x=b$mids1, y=b$mids2, z=Z, shade=TRUE,
      xlab="X", ylab="Y", main="",
      theta=45, phi=30, ltheta=60)