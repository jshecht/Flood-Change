# MONTE CARLO ANALYSIS

# Number of runs
n_runs = 100;
n_models = 7;

# Years of data
t <- seq(1,76,1)

# Establish the "true" 99th percentile (100-year flood) for each year
qpeak_99 <- array(NA,dim=c(76,1))

# Create array for sampled peak flows
qpeak_smpl <- array(NA,dim=c(76,1))

# Establish arrays for storing Monte Carlo results
ts_results_all <- array(NA,dim=c(n_runs,length(qpeak_99),n_models)) #Time series of different models
pbias_results_all <- array(NA,dim=c(n_runs,n_models)) #Mean bias of different models (MAKE RELATIVE)
rmse_results_all <- array(NA,dim=c(n_runs,n_models)) #RMSE of different models (MAKE RELATIVE)
res_norm_all <- array(NA,dim=c(n_runs,4)) #Residual normality, Anscombe-transformed residual normality

# Statistical parameters
mean_0 = 600
sd_0 = 150

b_mean = 2
b_sd = 0.5

# Try different values of d parameter
#d <- c(1/3,1/2,2/3)

for (m in 1:n_runs){

# Compute "true" values   
for (i in 1:76){
  qpeak_99[i]=qnorm(0.99,mean_0 + b_mean*(i-1),sd_0 + b_sd*(i-1))
}

plot(qpeak_99)

# Now sample values from this same distribution
for (i in 1:76){
  qpeak_smpl[i]=qnorm(runif(1,0,1),mean_0 + b_mean*(i-1),sd_0 + b_sd*(i-1))
  if(qpeak_smpl[i]<0){
    qpeak_smpl[i]=0
  }
}

plot(qpeak_smpl)
ppcc_qpeak_smpl <- ppcc.test(qpeak_smpl)
pval_ppcc_qpeak_smpl <- as.numeric(ppcc_qpeak_smpl[2])

# Conditional mean model for both methods
ols_reg <- lm(qpeak_smpl~t)
ols_fitted <- fitted(ols_reg)
ols_summary <- summary(ols_reg)
ols_coeff <- as.matrix(coefficients(ols_reg))
int_mean_est <- ols_coeff[1,]
b_mean_est <- ols_coeff[2,]
ols_cor <- cor(t,qpeak_smpl)

# Compute residuals
residuals_fit_ols <- residuals(ols_reg)

# Test residual normality
ppcc_res <- ppcc.test(residuals_fit_ols)
pval_ppcc_res <- as.numeric(ppcc_res[2])
if (pval_ppcc_res < 0.05) {
  res_norm = FALSE
  #stop("Residuals not normally distributed")
} else {
  res_norm = TRUE
}

# Absolute residuals do not yield a normal distribution

# Simplified Anscombe residual transformation from Davididan (2009 class notes)
residuals_fit_ols_2_3 <- (residuals(ols_reg)^2)^(1/3)

# Residuals using formula from McCullagh and Nelder (1989)
residuals_fit_ols_ansc2 <- 3*((qpeak_smpl^2)^(1/3)-(ols_fitted^2)^(1/3))/(ols_fitted^2)^(1/3)


# Check for normality
ppcc_res_2_3 <- ppcc.test(residuals_fit_ols_2_3)
pval_ppcc_res_2_3 <- as.numeric(ppcc_res_2_3[2])
if (pval_ppcc_res_2_3 < 0.05) {
  #stop("Anscombe-transformed residuals not normally distributed")
  res_2_3_norm = FALSE
} else {
  res_2_3_norm = TRUE
}

ppcc_res_2_3_ansc2 <- ppcc.test(residuals_fit_ols_ansc2)
pval_ppcc_res_2_3_ansc2 <- as.numeric(ppcc_res_2_3_ansc2[2])
if (pval_ppcc_res_2_3_ansc2 < 0.05) {
  #stop("Anscombe-transformed residuals not normally distributed")
  res_2_3_ansc2_norm = FALSE
} else {
  res_2_3_ansc2_norm = TRUE
}

# Time transformation
#t[d] <- (t^2)^(1/2*d)
t_1_3 <- t^(1/3)
t_2_3 <- t^(2/3)

# Plot time transformation vs. residual transformation
plot(t_1_3,residuals_fit_ols_2_3)
# time becomes non-uniform when power-transformed

# Run second regression of transformed residuals on transformed time
ols2_residuals_reg_t_1_3 <- lm(residuals_fit_ols_2_3~t_1_3)
summary(ols2_residuals_reg_t_1_3)
ols2_residuals_t_1_3_coeff <- as.matrix(coefficients(ols2_residuals_reg_t_1_3))
int_residual_t_1_3_est <- ols2_residuals_t_1_3_coeff[1,]
b_residual_t_1_3_est <- ols2_residuals_t_1_3_coeff[2,]
ols2_residuals_fitted_t_1_3 <- fitted(ols2_residuals_reg_t_1_3)
lines(ols2_residuals_fitted_t_1_3)

ols2_residuals_reg_t_2_3 <- lm(residuals_fit_ols_2_3~t_2_3)
summary(ols2_residuals_reg_t_2_3)
ols2_residuals_t_2_3_coeff <- as.matrix(coefficients(ols2_residuals_reg_t_2_3))
int_residual_t_2_3_est <- ols2_residuals_t_2_3_coeff[1,]
b_residual_t_2_3_est <- ols2_residuals_t_2_3_coeff[2,]
ols2_residuals_fitted_t_2_3 <- fitted(ols2_residuals_reg_t_2_3)
lines(ols2_residuals_fitted_t_2_3)

ols2_residuals_reg_t <- lm(residuals_fit_ols_2_3~t)
summary(ols2_residuals_reg_t)
ols2_residuals_t_coeff <- as.matrix(coefficients(ols2_residuals_reg_t))
int_residual_t_est <- ols2_residuals_t_coeff[1,]
b_residual_t_est <- ols2_residuals_t_coeff[2,]
ols2_residuals_fitted_t <- fitted(ols2_residuals_reg_t)
lines(ols2_residuals_fitted_t)

# Transform fitted Anscombe residuals into estimates of residual variance
ols2_res_var_t_1_3 <- ols2_residuals_fitted_t_1_3^3
lines(ols2_res_var_t_1_3)

ols2_res_var_t_2_3 <- ols2_residuals_fitted_t_2_3^3
lines(ols2_res_var_t_2_3)

ols2_res_var_t <- ols2_residuals_fitted_t^3

# Compute total standard deviation
ols2_stdev_t_1_3  <- sqrt(ols2_res_var_t_1_3 + (b_mean_est^2*var(t)))
ols2_stdev_t_2_3  <- sqrt(ols2_res_var_t_2_3 + (b_mean_est^2*var(t)))
ols2_stdev_t  <- sqrt(ols2_res_var_t + (b_mean_est^2*var(t)))

# Estimate quantiles
q99_est_ols2_t_1_3 <- qnorm(0.99,ols_fitted,ols2_stdev_t_1_3)
q99_est_ols2_t_2_3 <- qnorm(0.99,ols_fitted,ols2_stdev_t_2_3)
q99_est_ols2_t <- qnorm(0.99,ols_fitted,ols2_stdev_t)

# Derived variance model for t^(1/3), r^(2/3)
nu_2_t_1_3 <- (b_mean_est^2*(length(t)-1)^2)/(6*(length(t)+1))*(1-ols_cor^2)/ols_cor^2
res_var_t_1_3 <- nu_2_t_1_3*t
deriv_stdev_t_1_3 <- sqrt(res_var_t_1_3 + b_mean_est^2*var(t))
q99_est_deriv_t_1_3 <- qnorm(0.99,ols_fitted,deriv_stdev_t_1_3)
q99_est_deriv_t_1_3_sd_0 <- qnorm(0.99,ols_fitted+2.33*sd_0,deriv_stdev_t_1_3)

# Derived variance model for t^(2/3), r^(2/3)
nu2_t_2_3 <- (b_mean_est^2*(length(t)-1)^2)/(4*(length(t)^2+length(t)+1))*(1-ols_cor^2)/ols_cor^2
res_var_t_2_3 <- nu2_t_2_3*t^2
deriv_stdev_t_2_3 <- sqrt(res_var_t_2_3 + b_mean_est^2*var(t))
q99_est_deriv_t_2_3 <- qnorm(0.99,ols_fitted,deriv_stdev_t_2_3)
q99_est_deriv_t_2_3_sd_0 <- qnorm(0.99,ols_fitted+2.33*sd_0,deriv_stdev_t_2_3)

c = 0.2
# Derived variance model for t^(1/3), r^(2/3), c > 0
# UPDATE FORMULA 
c_nu2_t_1_3 <- b_mean_est^2*(length(t)-1)^2/(4*(3*c^2+length(t)^2+(3*c+1)*(length(t)+1)))*(1-ols_cor^2)/ols_cor^2
res_var_c_t_1_3 <- c_nu2_t_1_3*t^2
deriv_stdev_c_t_1_3 <- sqrt(res_var_c_t_1_3 + b_mean_est^2*var(t))
q99_est_deriv_c_t_1_3 <- qnorm(0.99,ols_fitted,deriv_stdev_c_t_1_3)

# Derived variance model for t^(1/3), r^(2/3), c > 0
c_nu2_t_2_3 <- b_mean_est^2*(length(t)-1)^2/(4*(3*c^2+length(t)^2+(3*c+1)*(length(t)+1)))*(1-ols_cor^2)/ols_cor^2
res_var_c_t_2_3 <- c_nu2_t_2_3*t^2
deriv_stdev_c_t_2_3 <- sqrt(res_var_c_t_2_3 + b_mean_est^2*var(t))
q99_est_deriv_c_t_2_3 <- qnorm(0.99,ols_fitted,deriv_stdev_c_t_2_3)

# Plot results 
plot.new()
plot(qpeak_99,typ="l",lwd=2,col="black",ylim=c(0,max(qpeak_99*2)),xlab="Year in Record",ylab="Value")
lines(q99_est_ols2_t_1_3,col="red")
lines(q99_est_ols2_t_2_3,col="purple")
lines(q99_est_ols2_t,col="pink")
lines(q99_est_deriv_t_1_3,col="green")
lines(q99_est_deriv_t_2_3,lty=2,col="blue")
lines(q99_est_deriv_c_t_1_3,col="green")
lines(q99_est_deriv_c_t_2_3,lty=2,col="blue")
leg.txt <- c("Truth",
             "2-step OLS t^1/3",
             "2-step OLS t^2/3",
             "2-step OLS t",
             "Deriv t^1/3",
             "Deriv t^2/3",
             "Deriv t^1/3 w c",
             "Deriv t^2/3 w c")

#legend("bottom",leg.txt)

# Compile results for all estimates
ts_results <- array(NA,dim=c(length(qpeak_99),7))
ts_results = c(q99_est_ols2_t_1_3,
               q99_est_ols2_t_2_3,
               q99_est_ols2_t,
               q99_est_deriv_t_1_3,
               q99_est_deriv_t_2_3,
               q99_est_deriv_c_t_1_3,
               q99_est_deriv_c_t_2_3)

# Compute errors
pbias_results <- array(NA,dim=7)
pbias_q99_est_ols2_t_1_3 = mean(abs(q99_est_ols2_t_1_3-qpeak_99)/qpeak_99)
pbias_q99_est_ols2_t_2_3 = mean(abs(q99_est_ols2_t_2_3-qpeak_99)/qpeak_99)
pbias_q99_est_ols2_t = mean(abs(q99_est_ols2_t-qpeak_99)/qpeak_99)
pbias_q99_est_deriv_t_1_3 = mean(abs(q99_est_deriv_t_1_3-qpeak_99)/qpeak_99)
pbias_q99_est_deriv_t_2_3 = mean(abs(q99_est_deriv_t_2_3-qpeak_99)/qpeak_99)
pbias_q99_est_deriv_c_t_1_3 = mean(abs(q99_est_deriv_c_t_1_3-qpeak_99)/qpeak_99)
pbias_q99_est_deriv_c_t_2_3 = mean(abs(q99_est_deriv_c_t_2_3-qpeak_99)/qpeak_99)

pbias_results = c(pbias_q99_est_ols2_t_1_3,
                 pbias_q99_est_ols2_t_2_3,
                 pbias_q99_est_ols2_t,
                 pbias_q99_est_deriv_t_1_3,
                 pbias_q99_est_deriv_t_2_3,
                 pbias_q99_est_deriv_c_t_1_3,
                 pbias_q99_est_deriv_c_t_2_3)

# Compute RMSE 
rmse_results <- array(NA,dim=7)
rmse_q99_est_ols2_t_1_3 = sqrt(mean((q99_est_ols2_t_1_3-qpeak_99)^2))
rmse_q99_est_ols2_t_2_3 = sqrt(mean((q99_est_ols2_t_2_3-qpeak_99)^2))
rmse_q99_est_ols2_t = sqrt(mean((q99_est_ols2_t-qpeak_99)^2))
rmse_q99_est_deriv_t_1_3 = sqrt(mean((q99_est_deriv_t_1_3-qpeak_99)^2))
rmse_q99_est_deriv_t_2_3 = sqrt(mean((q99_est_deriv_t_2_3-qpeak_99)^2))
rmse_q99_est_deriv_c_t_1_3 = sqrt(mean((q99_est_deriv_c_t_1_3-qpeak_99)^2))
rmse_q99_est_deriv_c_t_2_3 = sqrt(mean((q99_est_deriv_c_t_2_3-qpeak_99)^2))

# CHANGE TO RELATIVE RMSE!!!
rmse_results = c(rmse_q99_est_ols2_t_1_3,
                 rmse_q99_est_ols2_t_2_3,
                 rmse_q99_est_ols2_t,
                 rmse_q99_est_deriv_t_1_3,
                 rmse_q99_est_deriv_t_2_3,
                 rmse_q99_est_deriv_c_t_1_3,
                 rmse_q99_est_deriv_c_t_2_3)

ts_results_all[m,,] = ts_results
# Add: Comparison for last year in time series only

# Boxplots

boxplot_labels <- c("2-step OLS t^1/3",
                    "2-step OLS t^2/3",
                    "2-step OLS t",
                    "Deriv t^1/3",
                    "Deriv t^2/3",
                    "Deriv t^1/3 with c",
                    "Deriv t^2/3 with c")

pbias_results_all[m,] = pbias_results


rmse_results_all[m,] = rmse_results


res_norm_all[m,] = c(pval_ppcc_qpeak_smpl,pval_ppcc_res,pval_ppcc_res_2_3,pval_ppcc_res_2_3_ansc2)

# Add relative RMSE, slope difference
}

# Analyze data
boxplot(pbias_results_all,names=boxplot_labels,ylab="Absolute bias (%)")
boxplot(rmse_results_all,names=boxplot_labels,ylab="RMSE")
boxplot(res_norm_all,names=c("Sampled Data","Residuals","Residuals^(2/3)","SAS Formula"),main="PPCC tests for normality",ylab="p-value (Null hyp. = Normal)")

# Residual normality
plot(res_norm_all[,1],res_norm_all[,2],xlab="Residual Normal PPCC Significance",ylab="Residual^(2/3) Normal PPCC Significance",xlim=c(0,1),ylim=c(0,1))

