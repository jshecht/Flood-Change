# Nonlinear equation solver test

dslnex <- function(x) {
  y <- numeric(2)
  y[1] <- x[1]^2 + x[2]^2 - 2
  y[2] <- exp(x[1]-1) + x[2]^3 - 2
  y
}

xstart <- c(2,0.5)
fstart <- dslnex(xstart)
xstart 
fstart 