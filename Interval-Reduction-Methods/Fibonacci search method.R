#################### Fibonacci line search method ####################
##### Copyright (C) 2016.11.16 Jaeyoung Jung all rights reserved. ####

rm(list=ls(all=TRUE))
gc(reset=TRUE)


### Target function : Problem 4.3
z <- function(x) {
  y = 3 * x^2 - 4 * x + 1
  return(y)
}

# Draw function graph
curve(z, 0, 2)


### Initial value|
# Interval
a_0 <- 0
b_0 <- 2
e <- 0.01

# Convergence ratio
epsilon <- 0.001

## Make fibonacci number by using the predetermined number of iterations
# Determine # of iteration meets this condition
Fn <- (b_0 - a_0) / epsilon

Niter <- 10000 # Sufficient number
Nfib <- rep(0, Niter)
Nfib[1] <- 1 # Because the start number of R is 1, not 0 like python => Same with F_0 in book
Nfib[2] <- 1

Ni <- 3
while (Ni < Niter) {
    Nfib[Ni] <- Nfib[Ni-1] + Nfib[Ni-2]
    
    if (Nfib[Ni] > Fn) {
        Maxiter <- Ni
        print(paste("Niter = ", Ni, ", ", "Fibonacci number is ", Nfib[Ni], sep=""))
        break
    } else {
        Ni <- Ni + 1
    }
}

a <- rep(0, Maxiter)
b <- rep(0, Maxiter)
a[1] <- a_0
b[1] <- b_0

r <- rep(0, Maxiter)
x_L <- rep(0, Maxiter)
x_R <- rep(0, Maxiter)

x_L[1] <- a[1] + (1 - Nfib[Maxiter-1] / Nfib[Maxiter]) * (b[1] - a[1])
x_R[1] <- a[1] + (Nfib[Maxiter-1] / Nfib[Maxiter]) * (b[1] - a[1])

### Optimization
for (i in 1:(Maxiter-2)) {
  
    if (b[i] - a[i] > epsilon) {
      
        if (z(x_L[i]) < z(x_R[i])) {
          a[i+1] <- a[i]
          b[i+1] <- x_R[i]
          x_R[i+1] <- x_L[i]
          x_L[i+1] <- a[i+1] + (1 - Nfib[Maxiter - i - 1]/Nfib[Maxiter - i] - e) * (b[i+1] - a[i+1])
          print(paste("Subinterval is ", "[", a[i], ",", b[i], "]", ", ", "Niter = ", i-1, ", ","Convergence interval = ", (b[i] - a[i]), sep = ""))
        } else {
          a[i+1] <- x_L[i]
          b[i+1] <- b[i]
          x_L[i+1] <- x_R[i]
          x_R[i+1] <- a[i+1] + (Nfib[Maxiter - i - 1]/Nfib[Maxiter - i] + e) * (b[i+1] - a[i+1])
          print(paste("Subinterval is ", "[", a[i], ",", b[i], "]", ", ", "Niter = ", i-1, ", ","Convergence interval = ", (b[i] - a[i]), sep = ""))
        }
    } else {
      print(paste("Subinterval is ", "[", a[i], ",", b[i], "]", ", ", "Niter = ", i-1, ", ","Convergence interval = ", (b[i] - a[i]), sep = ""))
      break()
    }
}
