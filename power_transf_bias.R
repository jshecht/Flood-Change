# This sheet examines the effect of a power transformation on the dependent variable of a regression
# equation. I only examine the performance of this procedure on the second-stage of the regression 
# modeling approach. To do this, I assume that the dependent variable is an observed quantity and that 
# "beta" expresses a true value of the unit increase in the dependent variable yy in response to a 
# unit change in an independent variable xx. This way, we can pinpoint the extent to which the power 
# transformation biases estimates of the value of yy at a maxmium xx value without the confoudning 
# effects of the estimatation accuracy of yy or beta.

# One million observations are simulated to avoid accuracy issues related to sampling error. 

# Computation of residual variance using n-2 degrees of freedom ignored due to large sample size
# Sample variance computed instead.

# Wilson (1990) shows that lognormal quantiles are relatively unbiased. 

# Create output vectors
cond_max <- vector(length=1000)
cond_mean_trans <- vector(length=1000)
cond_max_trans <- vector(length=1000)
cond_max_trans_hh <- vector(length=1000)
xx_yy_mean <- vector(length=1000)

for (j in 1:1000){

#Create a normally distributed variable
max_xx <- 10   # Maximum value of xx
mu <- 0
beta <- 0.01
xx <- runif(10000,0,max_xx)
sigma <- rnorm(10000,0,(0.50 - 0.01 * xx))
yy <- mu + sigma  

# Dependent variable 
# Theoretical value = mu + beta * mean(xx) = 0.75

# Square normally distributed dependent variable
# Equivalent to squaring RHS
yy_2 <- yy^2
  
# Compute sample mean of yy^2
mean(yy_2) # ~ 0.203

# Adjusting for transformation bias
(mu + beta * mean(xx))^2 + sigma^2 # ~ 4.50

# Compute max value of y^2 - THE TRUTH
(mu + beta * max_xx)^2 + sigma^2 # ~4.25

# Check the alternative formula for THE TRUTH
mean(yy_2) + cov(xx,yy_2)/var(xx) * (max(xx)-mean(xx)) # ~4.0972

plot(xx[1:100],yy[1:100])

# Compute modified Anscombe transformation
yy_2_3 <- (yy^2)^(1/3)

lm_test_1 <- lm(yy_2_3 ~ xx) 
hh0 <- as.numeric(lm_test_1$coefficients[1]) 
hh1 <- as.numeric(lm_test_1$coefficients[2]) 
res_hh <- residuals(lm_test_1)

# Estimate maximum value of Y without transformation bias correction
cond_max[j] = (hh0 + hh1*max_xx)^3 

# Estimate mean value of Y with transformation bias correction
cond_mean_trans[j] = (hh0+hh1*mean(xx))^3 + 3*var(res_hh)*(hh0+hh1*mean(xx)) + mean(res_hh^3)

# Next estimate component beta(xx,yy^2) * (xx_max - xx|E[yy^2])

# Estimate xx|E[yy^2]

# Solve the following equation: 

# E[Y^2] = (hh0 + hh1*xx)^3 + 3*var(res_hh)*(hh0 + hh1*xx) + mean((res_hh)^3)
# Let Z = hh0 + hh1*xx
# E[Y^2] = Z^3 + 3*var(res_hh)*Z + mean((res_hh)^3)
# Rewrite as cubic root equation
# Z^3 + 3*var(res_hh)*Z -(E[Y^2] - mean((res_hh)^3)) = 0


# Solve cubic equation
QQ <- var(res_hh)
RR <- 0.5*(cond_mean_trans[j]-mean(res_hh^3))
SS <- (RR + sqrt(QQ^3 + RR^2))^(1/3)
TT <- sign(RR - sqrt(QQ^3 + RR^2))*abs((RR - sqrt(QQ^3 + RR^2)))^(1/3)

# The real root of Z is S + T since QQ^3 + RR^2 > 0 

# xx|E[Y^2] = (Z - hh0)/hh1

xx_yy_mean[j] = ((SS + TT)-hh0)/hh1  

# See if the estimate of Y^2|X(max) is close to THE TRUTH
cond_max_trans[j] = cond_mean_trans[j] + cov(xx,yy^2)/var(xx) * (max(xx) - xx_yy_mean[j])

# Compare with method without beta (xx,yy^2)
cond_max_trans_hh[j] = (hh0+hh1*max(xx))^3 + 3*var(res_hh)*(hh0+hh1*max(xx)) + mean(res_hh^3)

cond_max_trans_hh_hat[j] = (hh0+hh1*max(xx))^3 + 3*var(res_hh)*(hh0+hh1*max(xx)) + mean(res_hh^3)

}

mean(yy)
mean(yy_2)
mu + beta * max_xx
(mu + beta * max_xx)^2
mean(cond_max)
mean(cond_mean_trans)
mean(xx_yy_mean)
mean(cond_max_trans)
mean(cond_max_trans_hh)



