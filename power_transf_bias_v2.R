# This sheet examines the effect of a power transformation on the dependent variable of a regression
# equation. I only examine the performance of this procedure on the second-stage of the regression 
# modeling approach.  

# Create output vectors
cond_mean_no_transf <- vector(mode="numeric",length=1000)
cond_mean_transf <- vector(mode="numeric",length=1000)
cond_max_no_transf <- vector(mode="numeric",length=1000)
cond_max_transf <- vector(mode="numeric",length=1000)

for (j in 1:1){

#Create a normally distributed variable
max_xx <- 10   # Maximum value of xx
mu <- 0
beta <- 0.01
xx <- seq(0.000001,max_xx,0.000001)
yy <- 0.50 + 0.01*xx + rnorm(1000000,0,0.1)

ppcc.test(yy)
cv_yy <- sd(yy)/mean(yy)
cor(xx,yy)

mean(yy^3)

# We know that yy|(xx = 10) is 0.60 (0.05 + 0.01*10)
# Therefore yy^3

# How well can we simulate it? 

# Now fit a model to the normal variable
lm.yy_xx <- lm(yy~xx)
bb0 <- as.numeric(lm.yy_xx$coefficients[1])
bb1 <- as.numeric(lm.yy_xx$coefficients[2])
ww <- residuals(lm.yy_xx)

# Estimate mean response without accounting for transformation bias 
cond_mean_no_transf[j] <- (bb0 + bb1*mean(xx))^3 
# Estimate mean response without accounting for transformation bias 
cond_mean_transf[j] <- (bb0 + bb1*mean(xx))^3 + 3*var(ww)*(bb0 + bb1*mean(xx)) 

# Estimate conditional mean at xx = 10 without accounting for transformation bias
cond_max_no_transf[j] <- (bb0 + bb1*max(xx))^3  
# Estimate conditional mean xx = 10 without accounting for transformation bias
cond_max_transf[j] <- (bb0 + bb1*max(xx))^3 + 3*var(ww)*(bb0 + bb1*max(xx)) 

} 



