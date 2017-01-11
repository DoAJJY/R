#### Optimize Rosenbrock's function by using Quasi-Newton method and inexact line search ####
################# Copyright (C) 2016.12.15 Jaeyoung Jung all rights reserved. ###############
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

# Inexact line search - Wolfe Powell step size strategy
c1 <- 10^(-4)
c2 <- 0.9 # 0 < c1 < c2 < 1

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
armijo_alp <- function(x, g, p) {
  l <- 0
  alp <- 1
  nu_1 <- 0.66
  nu_2 <- 0.99
  while (TRUE) {
    if (f(x+alp*p) <= f(x) + c1*alp*g%*%p) break
    alp <- runif(1, nu_1*alp, nu_2*alp)
  }
  return(alp)
}

#########################################################################################################################
######################################################## DFP method #####################################################
#########################################################################################################################


### Quasi-Newton method
# Decide initial point
x0 <- c(-1, -1)

# Initial Hessian matrix identity matrix
H0 <- matrix(rep(0, length(x0)*length(x0)), nrow = 2)
diag(H0) <- 1

# Optimization
Ni <- 0
epsilon <- 10^-4

while (TRUE) {
  # For initial iteration
  if (Ni == 0) {
    # Step 1 : Calculate the gradient
    g <- grad(x0)
    
    # Step 2 : Calculate the direction
    p <- -H0 %*% g
    
    # Step 3 : Inexact line search and find the alpha
    #alp <- wolfe_alp(x0, g, p)
    alp <- armijo_alp(x0, g, p)
    
    # Step 4 : Find the next point
    x <- x0 + alp * p
    
    # Step 5 : DFP method - Update the next iteration the H matrix
    delta_x <- x - x0
    delta_grad <- grad(x) - grad(x0)
    term_U_1 <- (delta_x %*% t(delta_x)) / matrix(rep((t(delta_x) %*% delta_grad), 4), nrow = 2)
    term_U_2 <- (H0 %*% delta_grad) %*% (t(H0 %*% delta_grad)) / matrix(rep(t(delta_grad) %*% H0 %*% delta_grad, 4), nrow = 2)
    U <- term_U_1 - term_U_2
    H <- H0 + U
    
    # Print the information
    print(paste("Niter = ", Ni, ", ", Ni+1, "th X is ", "[", format(x[1], nsmall = 6), ", ", 
                format(x[2], nsmall = 6), "]", ", ", "Convergence is ",
                format(sqrt(sum((x-x0)^2)), nsmall = 6), ", ",
                "Step size is ", format(alp, nsmall = 6), sep = ""))
    
  } else { # For iteration number > 0
    # Step 1 : Calculate the gradient
    g <- grad(x0)
    
    # Step 2 : Calculate the direction
    p <- -H0 %*% g
    
    # Step 3 : Inexact line search amd find the H matrix
    #alp <- wolfe_alp(x0, g, p)
    alp <- armijo_alp(x0, g, p)
    
    # Step 4 : Find the next point
    x <- x0 + alp * p
    
    # Step 5 : DFP method - Update the next iteration the H matrix
    delta_x <- x - x0
    delta_grad <- grad(x) - grad(x0)
    term_U_1 <- (delta_x %*% t(delta_x)) / matrix(rep((t(delta_x) %*% delta_grad), 4), nrow = 2)
    term_U_2 <- (H0 %*% delta_grad) %*% (t(H0 %*% delta_grad)) / matrix(rep(t(delta_grad) %*% H0 %*% delta_grad, 4), nrow = 2)
    U <- term_U_1 - term_U_2
    H <- H0 + U
    
    # Print the information
    print(paste("Niter = ", Ni, ", ", Ni+1, "th X is ", "[", format(x[1], nsmall = 6), ", ", 
                format(x[2], nsmall = 6), "]", ", ", "Convergence is ",
                format(sqrt(sum((x-x0)^2)), nsmall = 6), ", ",
                "Step size is ", format(alp, nsmall = 6), sep = ""))
  }
  
  if ((sqrt(sum((x-x0)^2)) < epsilon) && (sum(grad(x)^2) < epsilon)) break # Convergence
  
  Ni <- Ni + 1
  x0 <- x
  H0 <- H
}

#########################################################################################################################
######################################################## BFGS method ####################################################
#########################################################################################################################

## Remove lists except functions
rm(list = setdiff(ls(), lsf.str()))

# Inexact line search - Wolfe Powell step size strategy
c1 <- 10^(-4)
c2 <- 0.9 # 0 < c1 < c2 < 1

# Decide initial point
x0 <- c(-1, -1)

# Initial Hessian matrix identity matrix
H0 <- matrix(rep(0, length(x0)*length(x0)), nrow = 2)
diag(H0) <- 1

# Optimization
Ni <- 0
epsilon <- 10^-4

while (TRUE) {
  # For initial iteration
  if (Ni == 0) {
    # Step 1 : Calculate the gradient
    g <- grad(x0)
    
    # Step 2 : Calculate the direction
    p <- -H0 %*% g
    
    # Step 3 : Inexact line search and find the alpha
    #alp <- wolfe_alp(x0, g, p)
    alp <- armijo_alp(x0, g, p)
    
    # Step 4 : Find the next point
    x <- x0 + alp * p
    
    # Step 5 : DFP method - Update the next iteration the H matrix
    delta_x <- x - x0
    delta_grad <- grad(x) - grad(x0)
    term_U_1 <- rep((1 + (t(delta_grad) %*% H0 %*% delta_grad) / (t(delta_x) %*% delta_grad)), 4) *
      delta_x %*% t(delta_x) / rep((t(delta_x) %*% delta_grad), 4)
    term_U_2 <- (H0 %*% delta_grad %*% t(delta_x) + t(H0 %*% delta_grad %*% t(delta_x))) / 
      rep(t(delta_x) %*% delta_grad, 4)
    U <- term_U_1 - term_U_2
    H <- H0 + U
    
    # Print the information
    print(paste("Niter = ", Ni, ", ", Ni+1, "th X is ", "[", format(x[1], nsmall = 6), ", ", 
                format(x[2], nsmall = 6), "]", ", ", "Convergence is ",
                format(sqrt(sum((x-x0)^2)), nsmall = 6), ", ",
                "Step size is ", format(alp, nsmall = 6), sep = ""))
    
  } else { # For iteration number > 0
    # Step 1 : Calculate the gradient
    g <- grad(x0)
    
    # Step 2 : Calculate the direction
    p <- -H0 %*% g
    
    # Step 3 : Inexact line search amd find the H matrix
    alp <- armijo_alp(x0, g, p)
    
    # Step 4 : Find the next point
    x <- x0 + alp * p
    
    # Step 5 : DFP method - Update the next iteration the H matrix
    delta_x <- x - x0
    delta_grad <- grad(x) - grad(x0)
    term_U_1 <- rep((1 + (t(delta_grad) %*% H0 %*% delta_grad) / (t(delta_x) %*% delta_grad)), 4) *
      delta_x %*% t(delta_x) / rep((t(delta_x) %*% delta_grad), 4)
    term_U_2 <- (H0 %*% delta_grad %*% t(delta_x) + t(H0 %*% delta_grad %*% t(delta_x))) / 
      rep(t(delta_x) %*% delta_grad, 4)
    U <- term_U_1 - term_U_2
    H <- H0 + U
    
    # Print the information
    print(paste("Niter = ", Ni, ", ", Ni+1, "th X is ", "[", format(x[1], nsmall = 6), ", ", 
                format(x[2], nsmall = 6), "]", ", ", "Convergence is ",
                format(sqrt(sum((x-x0)^2)), nsmall = 6), ", ",
                "Step size is ", format(alp, nsmall = 6), sep = ""))
  }
  
  if ((sqrt(sum((x-x0)^2)) < epsilon) && (sum(grad(x)^2) < epsilon)) break # Convergence
  
  Ni <- Ni + 1
  x0 <- x
  H0 <- H
}