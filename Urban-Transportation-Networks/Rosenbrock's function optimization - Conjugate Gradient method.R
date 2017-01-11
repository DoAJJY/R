# Optimize Rosenbrock's function by using Conjugate Gradient method and inexact line search #
################# Copyright (C) 2016.12.15 Jaeyoung Jung all rights reserved. ###############
library(mosaic)
library(limSolve)

rm(list=ls(all=TRUE))
gc(reset=TRUE)


### Define function
# Target function - Rosenbrock's function
f <- function(x) {
  y <- (1 - x[1])^2 + 100 * (x[2] - x[1]^2)^2 # Generally, a = 1 and b = 100
  return(y)
}

# First order derivative
grad <- function(x) {
  y <- c(400 * x[1]^3 - 400 * x[1] * x[2] + 2 * x[1] - 2, 200 * (x[2] - x[1]^2))
  return(y)
}

# Contour graph


# Inexact line search - Wolfe Powell step size strategy
c1 <- 10^(-4)
c2 <- 0.1 # 0 < c1 < c2 < 1

wolfe_alp <- function(x, g, p) {
    s <- 1/2
    l <- 0
    alp <- 1
    while (TRUE) {
        if (f(x+alp*p) - f(x) <= c1*alp*g%*%p) break # Wolfe condition  / strong ->  & (grad(x+alp*p) %*% p >= c2*grad(x) %*% p)
        l <- l + 1
        alp <- alp * s^l
    }
    return(alp)
}


# Inexact line search - Armijo step size strategy
armijo_alp <- function(p, x) {
  l <- 0
  alp <- 1
  nu_1 <- 0.66
  nu_2 <- 0.99
  while (TRUE) {
    if (f(x+alp*p) <= f(x) + c1*alp*grad(x)%*%p) break
    alp <- runif(1, nu_1*alp, nu_2*alp)
  }
  return(alp)
}



### Initial value
x <- c(-1, 2) # Set random initial point x^1 in feasible region

Ni <- 0 # Set initial number of iteration
MaxIter <- 200
hist_alp <- NULL
hist_dir <- NULL
hist_beta <- NULL
hist_x <- as.data.frame(matrix(rep(0, MaxIter * (length(x)+1)), ncol = 3))
hist_x[1, ] <- c(x, f(x))
# Start optimization
epsilon <- 10^(-3)

while (TRUE) {
  if (Ni == 0) {
    # Step 1 : Determine step size alpha
    g <- grad(x)
    p <- -g
    alp <- wolfe_alp(x, g, p)
    #alp <- armijo_alp(p, x)
    
    # Step 2 : Determin next point x
    x_new <- x + alp * p
    
    # Step 3 : Calculate beta
    beta <- sum(grad(x_new)^2) / (-p %*% g)

    # Step 4 : Determine the next direction p
    p_new <- -grad(x_new) + beta * p
    
    hist_alp <- c(hist_alp, alp)
    hist_dir <- rbind(hist_dir, p, p_new)
    hist_beta <- c(hist_beta, beta)
    hist_x[Ni+1, ] <- c(x_new, f(x_new))
    
    # Print the information
    print(paste("Niter = ", Ni, ", ", Ni+1, "th X is ", "[", format(x_new[1], nsmall = 6), ", ", 
                format(x_new[2], nsmall = 6), "]", ", ", "Convergence is ",
                format((grad(x_new) %*% grad(x_new)) /(grad(x) %*% grad(x)), nsmall = 6), ", ",
                "Step size is ", format(alp, nsmall = 6), sep = ""))
    
  } else {
    # Step 1 : Determine step size alpha
    g <- grad(x)
    p <- p
    alp <- wolfe_alp(x, g, p)
    #alp <- armijo_alp(p, x)
    
    # Step 2 : Determin next point x
    x_new <- x + alp * p
    
    # Step 3 : Calculate beta
    beta <- sum(grad(x_new)^2) / (-p %*% g)
    
    # Step 4 : Determine the next direction p
    p_new <- -grad(x_new) + beta * p
    
    hist_alp <- c(hist_alp, alp)
    hist_dir <- rbind(hist_dir, p_new)
    hist_beta <- c(hist_beta, beta)
    hist_x[Ni+1, ] <- c(x_new, f(x_new))
    
    # Print the information
    print(paste("Niter = ", Ni, ", ", Ni+1, "th X is ", "[", format(x_new[1], nsmall = 6), ", ", 
                format(x_new[2], nsmall = 6), "]", ", ", "Convergence is ",
                format((grad(x_new) %*% grad(x_new)) /(grad(x) %*% grad(x)), nsmall = 6), ", ",
                "Step size is ", format(alp, nsmall = 6), sep = ""))
  }

  # Convergence test
  if ((grad(x_new) %*% grad(x_new)) /(grad(x) %*% grad(x)) <= epsilon) break
  
  # Go to the next iteration
  Ni <- Ni + 1
  x <- x_new
  p <- p_new
  #Sys.sleep(1.5)
}