# Run a Monte Carlo experiment with linear and exponential models

nn <- 50
tt <- runif(1000000,1,nn)

aa <- 5
bb <- 0.02
ss <- 0.02
ee <- rnorm(1000000,0,0.5+ss*sqrt(tt))

# Compute "true" values of y
yy <- aa + bb*tt + ee
sd(yy)/mean(yy)
# Compute Q99 at tt = nn
exp(aa+bb*nn+qnorm(0.99,0,1)*(0.5+ss*sqrt(nn))) #4131

yy_tt.lm <- lm(yy~tt)
bb0 <- as.numeric(yy_tt.lm$coefficients[1])
bb1 <- as.numeric(yy_tt.lm$coefficients[2])
mean_nn <- bb0 + bb1*nn
res_yy_tt <- residuals(yy_tt.lm)
res_yy_tt_2 <- res_yy_tt^2
res_tt.lm <- lm(res_yy_tt_2 ~ tt) 
res2_tt <- residuals(res_tt.lm)
cc0 <- as.numeric(res_tt.lm$coefficients[1])
cc1 <- as.numeric(res_tt.lm$coefficients[2])
# Estimate residual variance at nn
res_var_nn <- cc0 + cc1*nn 

# Compute quantile based on time conditional residual variance
exp(mean_nn+qnorm(0.99,0,1)*sqrt(res_var_nn)) #3910

# Compute yy_tt 
res_yy_tt_2_3 <- (res_yy_tt^2)^(1/3)
res_tt_2_3_A.lm <- lm(res_yy_tt_2_3 ~ tt) 
cc0 <- as.numeric(res_tt_2_3_A.lm$coefficients[1])
cc1 <- as.numeric(res_tt_2_3_A.lm$coefficients[2])
# Estimate residual variance at nn
res_var_nn_2_3_A <- (cc0 + cc1*nn + 0.861*var(residuals(res_tt_2_3_A.lm))^(3/2))^3
# Compute quantile based on time conditional residual variance
exp(mean_nn+qnorm(0.99,0,1)*sqrt(res_var_nn_2_3_A)) #2184

res_yy_tt_2_3 <- (res_yy_tt^2)^(1/3)
tt_1_3 <- tt^(1/3)
res_tt_2_3_B.lm <- lm(res_yy_tt_2_3 ~ tt_1_3) 
cc0 <- as.numeric(res_tt_2_3_B.lm$coefficients[1])
cc1 <- as.numeric(res_tt_2_3_B.lm$coefficients[2])
# Estimate residual variance at nn
res_var_nn_2_3_B <- (cc0 + cc1*nn^(1/3))^3
# Compute quantile based on time conditional residual variance
exp(mean_nn+qnorm(0.99,0,1)*sqrt(res_var_nn_2_3_B)) #1955

# Now try exponential model
res_tt_exp.lm <- lm(res_yy_tt_2 ~ exp(tt))
cc0_exp <- as.numeric(res_tt_exp.lm$coefficients[1])
cc1_exp <- as.numeric(res_tt_exp.lm$coefficients[2])
res_tt_exp_res <- residuals(res_tt_exp.lm)
# Estimate residual variance at nn
res_var_nn_exp <- cc0_exp+cc1_exp*exp(nn)

# Compute quantile estimate based on time conditional variance
exp(mean_nn+qnorm(0.99,0,1)*sqrt(res_var_nn_exp)) #6259


# Try our model: new 1-parameter with optimization

cc1 = 1
# estimate residual variance 
int <- (length(tt)-1)/length(tt)*var(res_yy_tt)
# find optimal cc1 value
cc1_min_RSS <- function(cc1) {
  sum((res_yy_tt^2 - (int + cc1*(tt-mean(tt))))^2)
}

result <- optimize(cc1_min_RSS,interval=c(-1,1))

cc1_est <- result$minimum 
RSS <- result$objective
totSS <- sum((res_yy_tt - mean(res_yy_tt))^2)
rsquared_est <- 1 - RSS/totSS

# Compute conditional variance and quantile
cond_var <- int + cc1_est*(tt-mean(tt))
cond_var_at_n <- int + cc1_est*(nn-mean(tt))
Q_99_at_n <- exp(mean_nn+qnorm(0.99,0,1)*sqrt(cond_var_at_n)) #3940


# Try our model: Old 2-parameter 

nu_2_mod2a <- var(res_yy_tt)/mean(tt)  
res_var_nn_mod2a <- nu_2_mod2a * nn
exp(mean_nn+qnorm(0.99,0,1)*sqrt(res_var_nn_mod2a)) #4898

nu_2_mod2b <- var(res_yy_tt)/(mean(tt)^2+var(tt))
res_var_nn_mod2b <- nu_2_mod2b * nn^2
exp(mean_nn+qnorm(0.99,0,1)*sqrt(res_var_nn_mod2b)) #8581

# New derived model

# Estimate variance regression parameters using method of moments
rho_tt_ee2 <- cor(tt,res_yy_tt_2)
var_ee2 <- var(res_yy_tt_2)
cc0 <- var(res_yy_tt)- (rho_tt_ee2 * var_ee2 * mean(tt)) / sd(tt)
cc1 <- rho_tt_ee2 * var_ee2 / sd(tt)
res_ee2_tt <- cc0 + cc1 * tt
res_var_nn_1param <- cc0 + cc1 * nn
plot(tt[1:1000],res_yy_tt_2[1:1000],col="blue")
lines(tt[1:1000],res_ee2_tt[1:1000])
exp(mean_nn+qnorm(0.99,0,1)*sqrt(res_var_nn_1param))

# Transform the new model
rho_tt_1_3_ee_2_3 <- cor(tt_1_3,res_yy_tt_2_3)
var_ee_2_3 <- var(res_yy_tt_2_3)

res_var_nn_1param_1_3 <- cc0_1_3 + cc1_1_3 * nn^(1/3)
plot(tt_1_3[1:1000],res_yy_tt_2_3[1:1000],col="blue")
lines(tt_1_3[1:1000],res_ee_2_3_tt[1:1000])
exp(mean_nn+qnorm(0.99,0,1)*sqrt((res_var_nn_1param_1_3)^3))

cc0_1_3 <- sign(cc0) * abs(var(res_yy_tt) - (rho_tt_ee2 * var_ee2 * mean(tt)) / sd(tt)) ^ (1/3)
cc1_1_3 <- sign(cc1) * abs(var(res_yy_tt) * var_ee2 / sd(tt)) ^ (1/3)

# Try GLM with gamma
glm_Gamma <- glm(res_yy_tt_2 ~ tt,family=Gamma)
Gamma_b0 <- as.numeric(glm_Gamma$coefficients[1])
Gamma_b1 <- as.numeric(glm_Gamma$coefficients[2])
fitted_Gamma <- fitted(glm_Gamma)
res_var_nn_gamma <- Gamma_b0 + Gamma_b1*nn
# anscombe.residuals(glm_Gamma)
exp(mean_nn+qnorm(0.99,0,1)*sqrt(res_var_nn_gamma)) #3351

# Try GLM with gamma and log link function 
glm_Gamma_log <- glm(res_yy_tt_2 ~ tt,family=Gamma(link="log"))
Gamma_log_b0 <- as.numeric(glm_Gamma_log$coefficients[1])
Gamma_log_b1 <- as.numeric(glm_Gamma_log$coefficients[2])
fitted_Gamma_log <- fitted(glm_Gamma_log)
res_var_nn_gamma_log <- exp(Gamma_log_b0 + Gamma_log_b1*nn) # NOT SURE ABOUT THIS STEP
# anscombe.residuals(glm_Gamma)
exp(mean_nn+qnorm(0.99,0,1)*sqrt(res_var_nn_gamma_log)) #4473

# Try GLM with gamma and inverse link function
glm_Gamma_inv <- glm(res_yy_tt_2 ~ tt,family=Gamma(link="inverse"))
Gamma_inv_b0 <- as.numeric(glm_Gamma_inv$coefficients[1])
Gamma_inv_b1 <- as.numeric(glm_Gamma_inv$coefficients[2])
fitted_Gamma_inv <- fitted(glm_Gamma_inv)
res_var_nn_gamma_inv <- Gamma_inv_b0 + Gamma_inv_b1*nn
# anscombe.residuals(glm_Gamma)
exp(mean_nn+qnorm(0.99,0,1)*sqrt(res_var_nn_gamma_inv)) #3350
