####################### False position method ########################
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


### Newton method
# Initial values(condition of problem 4.7)
x_0 <- 0
x_1 <- 0.5

Niter <- 1000 # Sufficiently large number

# Define newton method as a function
false_position <- function(f, tol, x_0, x_1, N) {
    # First order derivative
    df <- D(f(x) ~ x)
    
    i <- 1
    p <- numeric(N)
    
    while (i <= N) {
        x_2 <- x_1 - df(x_1) / ((df(x_0) - df(x_1)) / (x_0 - x_1))
        p[i] <- x_2
        i <- i+1
        if (abs(x_2 - x_1) < tol) break
        x_0 <- x_1
        x_1 <- x_2
        print(paste("Niter = ", i, ", ", "x^(n+1) = ", x_1, sep = ""))
    }
    
    return(p[1:(i-1)])
}

# Optimization
false_position(z, epsilon, x_0, x_1, Niter)
