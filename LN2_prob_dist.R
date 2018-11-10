# Probability distributions

zz <- seq(0,200,0.1)
plot(zz,dlnorm(zz,4,0.5),typ="l")
lines(zz,dlnorm(zz,4,0.75),typ="l",lty=2)
lines(zz,dlnorm(zz,4.5,0.5),typ="l",col="red")
lines(zz,dlnorm(zz,4.5,0.75),typ="l",lty=2,col="red")

# Plot 1: Change in mean

# Plot 2: Change in mean & Cv

# Plot 3: No change in upper tail