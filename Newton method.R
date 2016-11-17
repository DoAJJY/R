########################## Newton's method ###########################
##### Copyright (C) 2016.11.16 Jaeyoung Jung all rights reserved. ####

rm(list=ls(all=TRUE))
gc(reset=TRUE)

library(mosaic)


### Define function and derivatives
# Original function
a <- 0
b <- 4

z <- makeFun(2*x^3 + 3*x^2 - 12*x + 5 ~ x)

curve(z, a, b)

# Tolerence
epsilon <- 0.001

# First order derivative
dz <- D(z(x) ~ x)

# Second order derivative
dzdz <- D(dz(x) ~ x)


### Newton method
# Initial values(condition of problem 4.7)
x_0 <- 0

Niter <- 1000 # Sufficiently large number

# Define newton method as a function
newton <- function(f, tol, x_0, N) {
    i <- 1
    x_1 <- x_0
    p <- numeric(N)
    
    while (i <= N) {
        x_1 <- x_0 - dz(x_0) / dzdz(x_0)
        p[i] <- x_1
        i <- i+1
        if (abs(x_1 - x_0) < tol) break
        x_0 <- x_1
        print(paste("Niter = ", i, ", ", "x^(n+1) = ", x_0, sep = ""))
    }
    
    return(p[1:(i-1)])
}

# Optimization
newton(z, epsilon, x_0, Niter)
