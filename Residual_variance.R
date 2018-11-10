# Normal dist properties

# Generate 1,000,000 homoscedastic residuals
E_e_homo <- rnorm(100000,0,1)
mean(abs(E_e_homo))
E_e_homo_2 <- E_e_homo^2
mean(E_e_homo_2)

# Generate a trend to put into the residuals
t_res <- seq(0.500001,1.5,0.000001)
sqrt_t_res <- sqrt(t_res)

# Generate residuals using the trend
E_e <- rnorm(1000000,0,sqrt_t_res)
# Mean absolute residual
mean(abs(E_e))
# Compute residual variance
E_e_2 <- E_e^2
# Compute mean of residual variance
mean(E_e_2)


# Generate trend in mean
aa = 10
bb = seq(0.000005,5,0.000005)
yy = aa + bb + E_e_homo
var(yy)

yy = aa + bb + E_e
var(yy)
