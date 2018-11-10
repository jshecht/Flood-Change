# Must have NWIS_Qpeak.R open while running this sheet

# EXTRACT RECORDS ---------------------------------------------------------------------

# Create vectors and lists to store output

station_group <- rbind(output_lin_mplus_sig05_vminus_sig05_clean,
                   output_lin_mplus_sig05_vplus_sig05_clean)


# Data generated in NWIS_Qpeak

# Conditional mean model

wy_order <- list()
wy_order_2 <- list()
wy_order_mean <- vector(length=nrow(station_group))
wy_order_max <- vector(length=nrow(station_group))

skew_peak_flow <- vector(length=nrow(station_group))
cv_peak_flow <- vector(length=nrow(station_group))
skew_ln_peak_flow <- vector(length=nrow(station_group))

b0_reg1 <- vector(length=nrow(station_group))
b1_reg1 <- vector(length=nrow(station_group))
rho_reg1 <- vector(length=nrow(station_group))
pval_b0_reg1 <- vector(length=nrow(station_group))
pval_b1_reg1 <- vector(length=nrow(station_group))
rsquared_reg1 <- vector(length=nrow(station_group))
res_dof_corr <- vector(length=nrow(station_group))
adj_r2_corr <- vector(length=nrow(station_group))
res_reg1 <- list()
res_var <- vector(length=nrow(station_group))
cond_med_at_n <- vector(length=nrow(station_group))
cond_mean_at_n <- vector(length=nrow(station_group)) # Accounts for transformation bias
Q_99_at_n_med_only <- vector(length=nrow(station_group))
Q_99_at_n_no_trend <- vector(length=nrow(station_group))

# Type II error analysis 
delta_b1_reg1_true <- vector(length=nrow(station_group))
tt_b1_reg1 <- vector(length=nrow(station_group))
t2_error_b1_reg1 <- vector(length=nrow(station_group))

# Two-stage least squares - Model A
b0_reg2_0a <- vector(length=nrow(station_group))
b1_reg2_0a <- vector(length=nrow(station_group))
pval_b0_reg2_0a <- vector(length=nrow(station_group))
pval_b1_reg2_0a <- vector(length=nrow(station_group))
adj_rsquared_reg2_0a <- vector(length=nrow(station_group))
res_reg2_0a <- list()
cond_var_at_n_2sls_0a <- vector(length=nrow(station_group)) 
Q_99_at_n_2sls_0a <- vector(length=nrow(station_group)) 

# Two-stage least squares - Model B
b0_reg2_0b <- vector(length=nrow(station_group))
b1_reg2_0b <- vector(length=nrow(station_group))
pval_b0_reg2_0b <- vector(length=nrow(station_group))
pval_b1_reg2_0b <- vector(length=nrow(station_group))
adj_rsquared_reg2_0b <- vector(length=nrow(station_group))
res_reg2_0b <- list()
cond_var_at_n_2sls_0b <- vector(length=nrow(station_group)) 
Q_99_at_n_2sls_0b <- vector(length=nrow(station_group)) 

# Two-stage least squares without Anscombe residuals - Model A
b0_reg2_1a <- vector(length=nrow(station_group))
b1_reg2_1a <- vector(length=nrow(station_group))
pval_b0_reg2_1a <- vector(length=nrow(station_group))
pval_b1_reg2_1a <- vector(length=nrow(station_group))
adj_rsquared_reg2_1a <- vector(length=nrow(station_group))
res_reg2_1a <- list()
cond_var_at_n_2sls_1a <- vector(length=nrow(station_group)) 
Q_99_at_n_2sls_1a <- vector(length=nrow(station_group))

# Two-stage least squares without Anscombe residuals - Model B
b0_reg2_1b <- vector(length=nrow(station_group))
b1_reg2_1b <- vector(length=nrow(station_group))
pval_b0_reg2_1b <- vector(length=nrow(station_group))
pval_b1_reg2_1b <- vector(length=nrow(station_group))
adj_rsquared_reg2_1b <- vector(length=nrow(station_group))
res_reg2_1b <- list()
cond_var_at_n_2sls_1b <- vector(length=nrow(station_group))
Q_99_at_n_2sls_1b <- vector(length=nrow(station_group))


# Two-parameter derived model - Model A 
cond_var_at_n_2a <- vector(length=nrow(station_group))
Q_99_at_n_2a <- vector(length=nrow(station_group))

# Two-parameter derived model - Model B
cond_var_at_n_2b <- vector(length=nrow(station_group)) 
Q_99_at_n_2b <- vector(length=nrow(station_group)) 

# Three-parameter model - Model A
cc1_3a <- vector(length=nrow(station_group))
rsquared_3a <- vector(length=nrow(station_group))
cond_var_at_n_3a <- vector(length=nrow(station_group))
Q_99_at_n_3a <- vector(length=nrow(station_group))

# Three-parameter model - Model C
cc1_3c <- vector(length=nrow(station_group))
rsquared_3c <- vector(length=nrow(station_group))
cond_var_at_n_3c <- vector(length=nrow(station_group))
Q_99_at_n_3c <- vector(length=nrow(station_group))

# Four-parameter - Model A
res_mult1_4a <- list()
cond_var_4a <- list()
cond_var_at_n_4a <- vector(length=nrow(station_group))
Q_99_at_n_4a <- vector(length=nrow(station_group))

# Model 5A (Method of moments, linear trend in variance)
cc0_5a <- vector(length=nrow(station_group))
cc1_5a <- vector(length=nrow(station_group))
rsquared_5a <- vector(length=nrow(station_group))
adj_rsquared_5a <- vector(length=nrow(station_group))
rmse_5a <- vector(length=nrow(station_group))
rrmse_5a <- vector(length=nrow(station_group))
mape_5a <- vector(length=nrow(station_group))
cond_var_5a <- list()
cond_var_at_n_5a <- vector(length=nrow(station_group))
Q_99_at_n_5a <- vector(length=nrow(station_group))
t_cc1_5a <- vector(length=nrow(station_group))
p_cc0_5a <- vector(length=nrow(station_group))
p_cc1_5a <- vector(length=nrow(station_group)) 

zp_RI_if_stnry_5a <- vector(length=nrow(station_group)) 
RI_if_stnry_5a <- vector(length=nrow(station_group)) 

delta_cc1_5a_true <- vector(length=nrow(station_group))
tt_cc1_5a <- vector(length=nrow(station_group))
t2_error_cc1_5a <- vector(length=nrow(station_group))

pval_ppcc_res_reg2_5a <- vector(length=nrow(station_group))
pval_dw_res_reg2_5a <- vector(length=nrow(station_group))
p_dd1_5a <- vector(length=nrow(station_group))

# Model 5B (Method of moments, quadratic trend in variance)
cc0_5b <- vector(length=nrow(station_group))
cc1_5b <- vector(length=nrow(station_group))
rsquared_5b <- vector(length=nrow(station_group))
adj_rsquared_5b <- vector(length=nrow(station_group))
rmse_5b <- vector(length=nrow(station_group))
rrmse_5b <- vector(length=nrow(station_group))

cond_var_5b <- list()
cond_var_at_n_5b <- vector(length=nrow(station_group))
Q_99_at_n_5b <- vector(length=nrow(station_group))
p_cc0_5b <- vector(length=nrow(station_group))
p_cc1_5b <- vector(length=nrow(station_group)) 

pval_ppcc_res_reg2_5b <- vector(length=nrow(station_group))
pval_dw_res_reg2_5b <- vector(length=nrow(station_group))
p_dd1_5b <- vector(length=nrow(station_group))

# Model 5C (Method of moments, log trend in variance)
cc0_5c<- vector(length=nrow(station_group))
cc1_5c <- vector(length=nrow(station_group))
rsquared_5c <- vector(length=nrow(station_group))
adj_rsquared_5c <- vector(length=nrow(station_group))
rmse_5c <- vector(length=nrow(station_group))
rrmse_5c <- vector(length=nrow(station_group))

cond_var_5c <- list()
cond_var_at_n_5c <- vector(length=nrow(station_group))

Q_99_at_n_5c <- vector(length=nrow(station_group))
p_cc0_5c <- vector(length=nrow(station_group))
p_cc1_5c <- vector(length=nrow(station_group)) 

pval_ppcc_res_reg2_5c <- vector(length=nrow(station_group))
pval_dw_res_reg2_5c <- vector(length=nrow(station_group))
p_dd1_5c <- vector(length=nrow(station_group))

# Model 5D (Method of moments, logarithmic trend in variance)
cc0_5d<- vector(length=nrow(station_group))
cc1_5d <- vector(length=nrow(station_group))
rsquared_5d <- vector(length=nrow(station_group))
adj_rsquared_5d <- vector(length=nrow(station_group))
rmse_5d <- vector(length=nrow(station_group))
rrmse_5d <- vector(length=nrow(station_group))

cond_var_5d <- list()
cond_var_at_n_5d <- vector(length=nrow(station_group))

Q_99_at_n_5d <- vector(length=nrow(station_group))
p_cc0_5d <- vector(length=nrow(station_group))
p_cc1_5d <- vector(length=nrow(station_group)) 

pval_ppcc_res_reg2_5d <- vector(length=nrow(station_group))
pval_dw_res_reg2_5d <- vector(length=nrow(station_group))
p_dd1_5d <- vector(length=nrow(station_group))

# Model 5E (Method of moments, log-transformed squared residuals)
cc0_5e<- vector(length=nrow(station_group))
cc1_5e <- vector(length=nrow(station_group))
rsquared_5e <- vector(length=nrow(station_group))
adj_rsquared_5e <- vector(length=nrow(station_group))
rmse_5e <- vector(length=nrow(station_group))
rrmse_5e <- vector(length=nrow(station_group))

cond_var_5e <- list()
cond_var_at_n_5e <- vector(length=nrow(station_group))
Q_99_at_n_5e <- vector(length=nrow(station_group))
p_cc0_5e <- vector(length=nrow(station_group))
p_cc1_5e <- vector(length=nrow(station_group)) 

pval_ppcc_res_reg2_5e <- vector(length=nrow(station_group))
pval_dw_res_reg2_5e <- vector(length=nrow(station_group))
p_dd1_5e <- vector(length=nrow(station_group))

# Model 5F (Method of moments, direct estimation s.d.)
cc0_5f<- vector(length=nrow(station_group))
cc1_5f <- vector(length=nrow(station_group))
rsquared_5f <- vector(length=nrow(station_group))
adj_rsquared_5f <- vector(length=nrow(station_group))
rmse_5f <- vector(length=nrow(station_group))
rrmse_5f <- vector(length=nrow(station_group))

cond_var_5f <- list()
cond_var_at_n_5f <- vector(length=nrow(station_group))
Q_99_at_n_5f <- vector(length=nrow(station_group))
p_cc0_5f <- vector(length=nrow(station_group))
p_cc1_5f <- vector(length=nrow(station_group)) 

pval_ppcc_res_reg2_5f <- vector(length=nrow(station_group))
pval_dw_res_reg2_5f <- vector(length=nrow(station_group))
p_dd1_5f <- vector(length=nrow(station_group))

# Model 6A (IWLS, linear trend in variance)
b0_reg1_iwls_6a <- matrix(nrow=nrow(station_group),ncol=10)
b1_reg1_iwls_6a <- matrix(nrow=nrow(station_group),ncol=10)
cc0_6a <- matrix(nrow=nrow(station_group),ncol=10)
cc1_6a <- matrix(nrow=nrow(station_group),ncol=10)

wt_6a_list <- list()
res_reg1_iwls_6a_list <- list()
res_reg1_iwls_2_3_6a_list <- list()
cond_var_6a_list <- list()

rsquared_6a <- vector(length=nrow(station_group))
adj_rsquared_6a <- vector(length=nrow(station_group))
rmse_6a <- vector(length=nrow(station_group))
rrmse_6a <- vector(length=nrow(station_group))

cond_var_at_n_6a <- vector(length=nrow(station_group))
Q_99_at_n_6a <- vector(length=nrow(station_group))
p_cc1_6a <- vector(length=nrow(station_group)) 

# Model 6B (IWLS, quadratic trend in variance)
b0_reg1_iwls_6b <- matrix(nrow=nrow(station_group),ncol=10)
b1_reg1_iwls_6b <- matrix(nrow=nrow(station_group),ncol=10)
cc0_6b <- matrix(nrow=nrow(station_group),ncol=10)
cc1_6b <- matrix(nrow=nrow(station_group),ncol=10)

wt_6b_list <- list()
res_reg1_iwls_6b_list <- list()
res_reg1_iwls_2_3_6b_list <- list()
cond_var_6b_list <- list()

rsquared_6b <- vector(length=nrow(station_group))
adj_rsquared_6b <- vector(length=nrow(station_group))
rmse_6b <- vector(length=nrow(station_group))
rrmse_6b <- vector(length=nrow(station_group))

cond_var_at_n_6b <- vector(length=nrow(station_group))
Q_99_at_n_6b <- vector(length=nrow(station_group))
p_cc1_6b <- vector(length=nrow(station_group)) 

# Model 6c (IWLS, exponential trend in variance)
b0_reg1_iwls_6c <- matrix(nrow=nrow(station_group),ncol=10)
b1_reg1_iwls_6c <- matrix(nrow=nrow(station_group),ncol=10)
cc0_6c <- matrix(nrow=nrow(station_group),ncol=10)
cc1_6c <- matrix(nrow=nrow(station_group),ncol=10)

wt_6c_list <- list()
res_reg1_iwls_6c_list <- list()
res_reg1_iwls_2_3_6c_list <- list()
cond_var_6c_list <- list()

rsquared_6c <- vector(length=nrow(station_group))
adj_rsquared_6c <- vector(length=nrow(station_group))
rmse_6c <- vector(length=nrow(station_group))
rrmse_6c <- vector(length=nrow(station_group))

cond_var_at_n_6c <- vector(length=nrow(station_group))
Q_99_at_n_6c <- vector(length=nrow(station_group))
p_cc1_6c <- vector(length=nrow(station_group))

# Model 6d (IWLS, logarithmic trend in variance)
b0_reg1_iwls_6d <- matrix(nrow=nrow(station_group),ncol=10)
b1_reg1_iwls_6d <- matrix(nrow=nrow(station_group),ncol=10)
cc0_6d <- matrix(nrow=nrow(station_group),ncol=10)
cc1_6d <- matrix(nrow=nrow(station_group),ncol=10)

wt_6d_list <- list()
res_reg1_iwls_6d_list <- list()
res_reg1_iwls_2_3_6d_list <- list()
cond_var_6d_list <- list()

rsquared_6d <- vector(length=nrow(station_group))
adj_rsquared_6d <- vector(length=nrow(station_group))
rmse_6d <- vector(length=nrow(station_group))
rrmse_6d <- vector(length=nrow(station_group))

cond_var_at_n_6d <- vector(length=nrow(station_group))
Q_99_at_n_6d <- vector(length=nrow(station_group))
p_cc1_6d <- vector(length=nrow(station_group))

# GLM - Model A 
b0_glm_cond_var_a <- vector(length=nrow(station_group))
b1_glm_cond_var_a <- vector(length=nrow(station_group)) 
pseudo_rsquared_glm_a <- vector(length=nrow(station_group)) 
cond_var_at_n_glm_a <- vector(length=nrow(station_group)) 
Q_99_at_n_glm_a <- vector(length=nrow(station_group))

# Gamma GLM model without IWLS
b0_glm_cond_var_gamma <- vector(length=nrow(station_group))
b1_glm_cond_var_gamma <- vector(length=nrow(station_group)) 
pseudo_rsquared_glm_gamma <- vector(length=nrow(station_group)) 
cond_var_at_n_glm_gamma <- vector(length=nrow(station_group)) 
Q_99_at_n_glm_gamma <- vector(length=nrow(station_group))

# Gamma GLM model with IWLS
b0_reg1_iwls_glm_gamma <- matrix(nrow=nrow(station_group),ncol=10)
b1_reg1_iwls_glm_gamma <- matrix(nrow=nrow(station_group),ncol=10)
cc0_iwls_glm_gamma <- matrix(nrow=nrow(station_group),ncol=10)
cc1_iwls_glm_gamma <- matrix(nrow=nrow(station_group),ncol=10)
res_reg1_iwls_glm_gamma <- matrix(nrow=nrow(station_group),ncol=10)
res_reg1_2_iwls_glm_gamma <- matrix(nrow=nrow(station_group),ncol=10)

wt_iwls_glm_gamma_list <- list()
res_reg1_iwls_glm_gamma_list <- list()
res_reg1_2_iwls_glm_gamma_list <- list()
cvar_iwls_glm_gamma_list <- list()

rsquared_iwls_glm_gamma <- vector(length=nrow(station_group))
adj_rsquared_iwls_glm_gamma <- vector(length=nrow(station_group))
rmse_iwls_glm_gamma <- vector(length=nrow(station_group))
rrmse_iwls_glm_gamma <- vector(length=nrow(station_group))

cond_var_at_n_iwls_glm_gamma <- vector(length=nrow(station_group))
Q_99_at_n_iwls_glm_gamma <- vector(length=nrow(station_group))
p_cc1_iwls_glm_gamma <- vector(length=nrow(station_group)) 

# Gamma GLM model with IWLS - quadratic
b0_reg1_iwls_glm_gamma_7b <- matrix(nrow=nrow(station_group),ncol=10)
b1_reg1_iwls_glm_gamma_7b <- matrix(nrow=nrow(station_group),ncol=10)
cc0_iwls_glm_gamma_7b <- matrix(nrow=nrow(station_group),ncol=10)
cc1_iwls_glm_gamma_7b <- matrix(nrow=nrow(station_group),ncol=10)
res_reg1_iwls_glm_gamma_7b <- matrix(nrow=nrow(station_group),ncol=10)
res_reg1_2_iwls_glm_gamma_7b <- matrix(nrow=nrow(station_group),ncol=10)

wt_iwls_glm_gamma_7b_list <- list()
res_reg1_iwls_glm_gamma_7b_list <- list()
res_reg1_2_iwls_glm_gamma_7b_list <- list()
cvar_iwls_glm_gamma_7b_list <- list()

rsquared_iwls_glm_gamma_7b <- vector(length=nrow(station_group))
adj_rsquared_iwls_glm_gamma_7b <- vector(length=nrow(station_group))
rmse_iwls_glm_gamma_7b <- vector(length=nrow(station_group))
rrmse_iwls_glm_gamma_7b <- vector(length=nrow(station_group))

cond_var_at_n_iwls_glm_gamma_7b <- vector(length=nrow(station_group))
Q_99_at_n_iwls_glm_gamma_7b <- vector(length=nrow(station_group))
p_cc1_iwls_glm_gamma_7b <- vector(length=nrow(station_group)) 

# Gamma GLM model with IWLS - logarithmic
b0_reg1_iwls_glm_gamma_7d <- matrix(nrow=nrow(station_group),ncol=10)
b1_reg1_iwls_glm_gamma_7d <- matrix(nrow=nrow(station_group),ncol=10)
cc0_iwls_glm_gamma_7d <- matrix(nrow=nrow(station_group),ncol=10)
cc1_iwls_glm_gamma_7d <- matrix(nrow=nrow(station_group),ncol=10)
res_reg1_iwls_glm_gamma_7d <- matrix(nrow=nrow(station_group),ncol=10)
res_reg1_2_iwls_glm_gamma_7d <- matrix(nrow=nrow(station_group),ncol=10)

wt_iwls_glm_gamma_7d_list <- list()
res_reg1_iwls_glm_gamma_7d_list <- list()
res_reg1_2_iwls_glm_gamma_7d_list <- list()
cvar_iwls_glm_gamma_7d_list <- list()

rsquared_iwls_glm_gamma_7d <- vector(length=nrow(station_group))
adj_rsquared_iwls_glm_gamma_7d <- vector(length=nrow(station_group))
rmse_iwls_glm_gamma_7d <- vector(length=nrow(station_group))
rrmse_iwls_glm_gamma_7d <- vector(length=nrow(station_group))

cond_var_at_n_iwls_glm_gamma_7d <- vector(length=nrow(station_group))
Q_99_at_n_iwls_glm_gamma_7d <- vector(length=nrow(station_group))
p_cc1_iwls_glm_gamma_7d <- vector(length=nrow(station_group)) 

# Gamma GLM model with IWLS - standard deviation
b0_reg1_iwls_glm_gamma_7f <- matrix(nrow=nrow(station_group),ncol=10)
b1_reg1_iwls_glm_gamma_7f <- matrix(nrow=nrow(station_group),ncol=10)
cc0_iwls_glm_gamma_7f <- matrix(nrow=nrow(station_group),ncol=10)
cc1_iwls_glm_gamma_7f <- matrix(nrow=nrow(station_group),ncol=10)
res_reg1_iwls_glm_gamma_7f <- matrix(nrow=nrow(station_group),ncol=10)
res_reg1_abs_iwls_glm_gamma_7f <- matrix(nrow=nrow(station_group),ncol=10)

wt_iwls_glm_gamma_7f_list <- list()
res_reg1_iwls_glm_gamma_7f_list <- list()
res_reg1_abs_iwls_glm_gamma_7f_list <- list()
cvar_iwls_glm_gamma_7f_list <- list()

rsquared_iwls_glm_gamma_7f <- vector(length=nrow(station_group))
adj_rsquared_iwls_glm_gamma_7f <- vector(length=nrow(station_group))
rmse_iwls_glm_gamma_7f <- vector(length=nrow(station_group))
rrmse_iwls_glm_gamma_7f <- vector(length=nrow(station_group))

cond_var_at_n_iwls_glm_gamma_7f <- vector(length=nrow(station_group))
Q_99_at_n_iwls_glm_gamma_7f <- vector(length=nrow(station_group))
p_cc1_iwls_glm_gamma_7f <- vector(length=nrow(station_group)) 

# MAIN LOOP STARTS 
for (k in 1:nrow(station_group)){  
wy <- wy_station[[station_group[k,1]]]
wy_order[[k]] = wy - min(wy) + 1
wy_order_2[[k]] = (wy - min(wy) + 1)^2
wy_order_mean[k] = mean(wy_order[[k]])
wy_order_max[k] = max(wy_order[[k]])
wyear <- wy_order[[k]] # Add in water year stuff
peak_flow <- peaks[[station_group[k,1]]]/35.3146662127
skew_peak_flow[k] = skewness(peak_flow)
cv_peak_flow[k] = sd(peak_flow)/mean(peak_flow)
skew_ln_peak_flow[k] = skewness(log(peak_flow+0.0001))
flood_reg1.lm = lm(log(peak_flow+0.0001) ~ wy_order[[k]]) # Add constant to avoid problems with ephemeral streams
b0_reg1[k] = summary(flood_reg1.lm)$coefficients[1,1]
b1_reg1[k] = summary(flood_reg1.lm)$coefficients[2,1]
rho_reg1[k] = cor(wy_order[[k]],log(peak_flow+0.01))
pval_b0_reg1[k] = summary(flood_reg1.lm)$coefficients[1,4]
pval_b1_reg1[k] = summary(flood_reg1.lm)$coefficients[2,4]
rsquared_reg1[k] = as.numeric(summary(flood_reg1.lm)$r.squared)
res_reg1[[k]] = residuals(flood_reg1.lm)
res_dof_corr[k] = length(wy_order[[k]])/(length(wy_order[[k]])-2)
adj_r2_corr[k] = (length(wy_order[[k]]) - 1)/(length(wy_order[[k]]) - 2)
res_var[k] = res_dof_corr[k] * var(res_reg1[[k]])
cond_med_at_n[k] = b0_reg1[k] + b1_reg1[k] * max(wy_order[[k]]) 
cond_mean_at_n[k] = cond_med_at_n[k] + 0.5 * res_var[k] #Assuming homoscedasticity
Q_99_at_n_med_only[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1) * sqrt(res_var[k]))
Q_99_at_n_no_trend[k] = exp(mean(log(peak_flow + 0.01)) + qnorm(0.99,0,1) * sd(log(peak_flow+0.01)))  

# View conditional mean model
plot(wy,log(peak_flow+0.01),xlab="Water year",ylab= "ln(Annual peak flow, cfs)") 
title(c(site_names[station_group[k,1]],
      as.character(site_id[station_group[k,1]]),
      "DA",
      as.numeric(as.character(site_info[station_group[k,1],6]))),cex.main=0.8) 
text(min(wy)+5,max(log(peak_flow+0.01))-0.5,bquote(~R^2 ==. (round(rsquared_reg1[k],3))),cex=0.7) 

# View transformed residuals assuming wy_order~res
plot(wy,(res_reg1[[k]]^2)^(1/3),xlab="Water year in record (t)",ylab= "Residual Variance ^ (2/3)")
title(c(site_names[station_group[k,1]],
     as.character(site_id[station_group[k,1]]),
     "DA",
     as.numeric(as.character(site_info[station_group[k,1],6]))),cex.main=0.8)
text(min(wy)+5,max((res_reg1[[k]]^2)^(1/3))-0.08,bquote(~R^2 ==.(round(rsquared_5a[k],3))),cex=0.7)

# Type II errors using methods from Vogel et al. (2013) and Rosner et al. (2014)
delta_b1_reg1_true[k] = 1/(sqrt(1/cor(wy_order[[k]],log(peak_flow+0.0001))^2-1))
tt_b1_reg1[k] = abs(qt(1-pval_b1_reg1[k],length(wy_order[[k]])-2))
t2_error_b1_reg1[k] = pt(tt_b1_reg1[k] - delta_b1_reg1_true[k]*sqrt(length(wy_order[[k]])),length(wy_order[[k]])-2)


# View transformed residuals assuming wy_order[[k]]~ln(res^2)
#log_res_reg1_2 <- log(res_reg1[[k]]^2)
#plot(wy_order[[k]],log_res_reg1_2,xlab="Water year in record (t)",ylab= "Residual Variance",pch=4)
#title(c(site_info[station_group[k,1]],
#        as.character(site_info[station_group[k,1]]),
#        "DA",
#        as.numeric(as.character(site_info[station_group[k,1],6]))),cex.main=0.8)

# Fit conditional variance model A using two-stage least squares with an Anscombe transformation
res_reg1_2_3 <- (res_reg1[[k]]^2)^(1/3)
wy_order_1_3 <- wy_order[[k]]^(1/3)
flood_reg2_0a.lm = lm(res_reg1_2_3 ~ wy_order_1_3) # Add constant to avoid problems with ephemeral streams
b0_reg2_0a[k] = summary(flood_reg2_0a.lm)$coefficients[1,1] 
b1_reg2_0a[k] = summary(flood_reg2_0a.lm)$coefficients[2,1] 
pval_b0_reg2_0a[k] = summary(flood_reg2_0a.lm)$coefficients[1,4] 
pval_b1_reg2_0a[k] = summary(flood_reg2_0a.lm)$coefficients[2,4] 
adj_rsquared_reg2_0a[k] = as.numeric(summary(flood_reg2_0a.lm)$adj.r.squared) 
res_reg2_0a[[k]] = residuals(flood_reg2_0a.lm)

# View transformed residuals
#plot(wy_order_1_3,(res_reg1[[k]]^2)^(1/3),xlab="Water year in record (t)",ylab= "Residual Variance",pch=3)
#title(c(site_names[station_group[k,1]],
#        as.character(site_id[station_group[k,1]]),
#        "DA",
#        as.numeric(as.character(site_info[station_group[k,1],6]))),cex.main=0.8)

# Compute residual variance in last year of record
cond_var_at_n_2sls_0a[k] = (b0_reg2_0a[k] + b1_reg2_0a[k] * max(wy_order_1_3))^3 
Q_99_at_n_2sls_0a[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_2sls_0a[k])) 

# Fit conditional variance model B using two-stage least squares with Anscombe transformation
res_reg1_2_3 <- (res_reg1[[k]]^2)^(1/3)
wy_order_2_3 <- wy_order[[k]]^(2/3)
flood_reg2_0b.lm = lm(res_reg1_2_3 ~ wy_order_2_3) # Add constant to avoid problems with ephemeral streams
b0_reg2_0b[k] = summary(flood_reg2_0b.lm)$coefficients[1,1]
b1_reg2_0b[k] = summary(flood_reg2_0b.lm)$coefficients[2,1]
pval_b0_reg2_0b[k] = summary(flood_reg2_0b.lm)$coefficients[1,4]
pval_b1_reg2_0b[k] = summary(flood_reg2_0b.lm)$coefficients[2,4]
adj_rsquared_reg2_0b[k] = as.numeric(summary(flood_reg2_0b.lm)$adj.r.squared)
res_reg2_0b[[k]] = residuals(flood_reg2_0b.lm)

# Compute residual variance in last year of record
cond_var_at_n_2sls_0b[k] = (b0_reg2_0b[k] + b1_reg2_0b[k] * max(wy_order_2_3))^3
Q_99_at_n_2sls_0b[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_2sls_0b[k]))

# Fit conditional variance Model A with two-stage least squares without an Anscombe transformation
res_reg1_2 <- (res_reg1[[k]]^2)
flood_reg2_1a.lm <- lm(res_reg1_2 ~ wy_order[[k]])
b0_reg2_1a[k] = summary(flood_reg2_1a.lm)$coefficients[1,1]
b1_reg2_1a[k] = summary(flood_reg2_1a.lm)$coefficients[2,1]
pval_b0_reg2_1a[k] = summary(flood_reg2_1a.lm)$coefficients[1,4]
pval_b1_reg2_1a[k] = summary(flood_reg2_1a.lm)$coefficients[2,4]
adj_rsquared_reg2_1a[k] = as.numeric(summary(flood_reg2_1a.lm)$adj.r.squared)
res_reg2_1a[[k]] = residuals(flood_reg2_1a.lm)

cond_var_at_n_2sls_1a[k] = b0_reg2_1a[k] + b1_reg2_1a[k] * max(wy_order[[k]])
Q_99_at_n_2sls_1a[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_2sls_1a[k]))

#plot(res_reg1[[k]]^2)
#lines(b0_reg2_1a[k]+b1_reg2_1a[k]*seq(1,length(res_reg1[[k]]),1))
#13,18,22,29,33,38

# Fit quadratic conditional variance model with two-stage least squares without an Anscombe transformation
res_reg1_2 <- (res_reg1[[k]]^2)
flood_reg2_1b.lm <- lm(res_reg1_2 ~ wy_order_2[[k]])
b0_reg2_1b[k] = summary(flood_reg2_1b.lm)$coefficients[1,1]
b1_reg2_1b[k] = summary(flood_reg2_1b.lm)$coefficients[2,1]
pval_b0_reg2_1b[k] = summary(flood_reg2_1b.lm)$coefficients[1,4]
pval_b1_reg2_1b[k] = summary(flood_reg2_1b.lm)$coefficients[2,4]
adj_rsquared_reg2_1b[k] = as.numeric(summary(flood_reg2_1b.lm)$adj.r.squared)
res_reg2_1b[[k]] = residuals(flood_reg2_1b.lm)

cond_var_at_n_2sls_1b[k] = b0_reg2_1b[k] + b1_reg2_1b[k] * max(wy_order_2[[k]])
Q_99_at_n_2sls_1b[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_2sls_1b[k]))


# Fit Model A using 2-parameter derived equation
var_res_reg1 <- length(wy_order[[k]])/(length(wy_order[[k]])-2)*var(as.numeric(res_reg1[[k]]))
mu_t <- (length(wy_order[[k]])+1)/2
cond_var_reg2_2a <- var_res_reg1/mu_t*wy_order[[k]]
cond_var_at_n_2a[k] = max(cond_var_reg2_2a)
Q_99_at_n_2a[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_2a[k]))

# Fit Model B using 2-parameter derived equation
var_res_reg1 <- length(wy_order[[k]])/(length(wy_order[[k]])-2)*var(as.numeric(res_reg1[[k]]))
mu_t_2 <- ((length(wy_order[[k]])+1)/2)^2
sigma_t_2 <- (length(wy_order[[k]])-1)^2/12
cond_var_reg2_2b <- var_res_reg1/(mu_t_2+sigma_t_2)*wy_order[[k]]^2
cond_var_at_n_2b[k] = max(cond_var_reg2_2b)
Q_99_at_n_2b[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_2b[k]))
# NOT DONE USING ANSCOMBE TRANSFORMATION

# Fit Model 4A - Multiplicative model
res_mult1_4a[[k]] = as.numeric((res_reg1[[k]] + b0_reg1[k] + b1_reg1[k] * wy_order[[k]])/(b0_reg1[k] + b1_reg1[k] * wy_order[[k]]))
cond_var_4a[[k]] = as.numeric((b0_reg1[k] + b1_reg1[k] * wy_order[[k]])^2*var(res_mult1_4a[[k]]))
cond_var_at_n_4a[k] = (b0_reg1[k] + b1_reg1[k] * max(wy_order[[k]]))^2*var(res_mult1_4a[[k]])
Q_99_at_n_4a[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_4a[k])) 

# Fit Model 5A 
rho_t_res_reg2 <- cor(wy_order[[k]],res_reg1_2_3)
sd_res_reg1_2_3 <- sd(res_reg1_2_3)
cc0_5a[k] = mean(res_reg1_2_3) - (rho_t_res_reg2 * sd_res_reg1_2_3  * mean(wy_order[[k]])) / sd(wy_order[[k]])
cc1_5a[k] = rho_t_res_reg2 * sd_res_reg1_2_3 / sd(wy_order[[k]])
res_reg2_5a_fit <- (cc0_5a[k] + cc1_5a[k] * wy_order[[k]]) 
res_reg2_5a_var <- res_dof_corr[k]*var(res_reg1_2_3-res_reg2_5a_fit)

cond_var_5a[[k]] <- (cc0_5a[k] + cc1_5a[k] * wy_order[[k]])^3 + 3*res_reg2_5a_var*(cc0_5a[k] + cc1_5a[k] * wy_order[[k]])
#cond_var_at_n_5a[k] = res_dof_corr[k]*(cc0_5a[k] + cc1_5a[k] * wy_order_max[k])^3 
cond_var_at_n_5a[k] = res_dof_corr[k]*(cc0_5a[k] + cc1_5a[k] * wy_order_max[k])^3 + 3*res_reg2_5a_var*(cc0_5a[k] + cc1_5a[k] * wy_order_max[k]) + mean((res_reg1_2_3-res_reg2_5a_fit)^3)

#plot(wy_order,res_reg1_2_3,col="blue")
#lines(wy_order,res_reg2_5a)
 
Q_99_at_n_5a[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5a[k])) 

rsquared_5a[k] = 1 - sum((res_reg1_2_3 - res_reg2_5a_fit)^2)/sum((res_reg1_2_3 - mean(res_reg1_2_3))^2)
adj_rsquared_5a[k] = 1 - adj_r2_corr[k]*sum((res_reg1_2_3 - res_reg2_5a_fit)^2)/sum((res_reg1_2_3 - mean(res_reg1_2_3))^2)
rmse_5a[k] = sqrt(1/length(wy_order[[k]])*sum((res_reg1_2_3^3 - cond_var_at_n_5a[[k]])^2))
rrmse_5a[k] = sqrt(1/length(wy_order[[k]])*sum(((res_reg1_2_3^3 - cond_var_at_n_5a[[k]])/cond_var_at_n_5a[[k]])^2))
mape_5a[k] = 100 * 1/length(wy_order[[k]])*sum(abs(res_reg1_2_3^3 - cond_var_at_n_5a[[k]])/cond_var_at_n_5a[[k]])
se_cc0_5a <- sqrt(sum((res_reg1_2_3 - res_reg2_5a_fit)^2)/(length(wy_order[[k]]) - 2)*(1/length(wy_order[[k]])+mean(wy_order[[k]])^2/sum((wy_order[[k]]-mean(wy_order[[k]]))^2))) 
se_cc1_5a <- sqrt(sum((res_reg1_2_3 - res_reg2_5a_fit)^2)/((length(wy_order[[k]]) - 2)*sum((wy_order[[k]]-mean(wy_order[[k]]))^2)))
t_cc0_5a <- cc0_5a[k]/se_cc0_5a 
t_cc1_5a[k] = cc1_5a[k]/se_cc1_5a
p_cc0_5a[k] = 2*(pt(-abs(t_cc0_5a), df=length(wy_order[[k]]) - 1))
p_cc1_5a[k] = 2*(pt(-abs(t_cc1_5a[k]), df=length(wy_order[[k]]) - 1)) 
#p_cc1_5a[k] = 1 - 2*pt(t_cc1_5a, df=length(wy_order[[k]]) - 1))

# Compute return period difference
# COmpute stationary 100-year flood
zp_RI_if_stnry_5a[k] = (log(Q_99_at_n_no_trend[k])-cond_med_at_n[k])/sqrt(cond_var_at_n_5a[k])
RI_if_stnry_5a[k] = 1/(1-pnorm(zp_RI_if_stnry_5a[k]))

# Type II errors using methods from Vogel et al. (2013) and Rosner et al. (2014)
delta_cc1_5a_true[k] = 1/(sqrt(1/cor(wy_order[[k]],res_reg1_2_3)^2-1)) 
tt_cc1_5a[k] = qt(1-p_cc1_5a[k],length(wy_order[[k]])-2)
t2_error_cc1_5a[k] = pt(tt_cc1_5a[k] - delta_cc1_5a_true[k]*sqrt(length(wy_order[[k]])),length(wy_order[[k]])-2)

# Test residual adequacy 
res_reg2_5a <- as.numeric(res_reg1_2_3 - res_reg2_5a_fit)
ppcc_res_reg2 <- ppcc.test(res_reg2_5a)
pval_ppcc_res_reg2_5a[k] = as.numeric(ppcc_res_reg2[2])
dw_res_reg2 <- dwtest(res_reg2_5a ~ wy_order[[k]])
pval_dw_res_reg2_5a[k] = as.numeric(dw_res_reg2[2]) 

# Test for heteroscedasticity of residuals of second regression
res_reg2_2_3_5a <- ((res_reg1_2_3 - res_reg2_5a_fit)^2)^(1/3)
res_reg2_2_3_5a.lm <- lm(res_reg2_2_3_5a ~ wy_order[[k]])
p_dd1_5a[k] = summary(res_reg2_2_3_5a.lm)$coefficients[2,4]

# View transformed residuals assuming wy_order~res
#plot(wy_order,(res_reg1[[k]]^2)^(1/3),xlab="Water year in record (t)",ylab= "Residual Variance ^ (2/3)")
#title(c(site_names[station_group[k,1]],
 #      as.character(site_id[station_group[k,1]]),
 #      "DA",
 #      as.numeric(as.character(site_info[station_group[k,1],6]))),cex.main=0.8) 
#text(max(wy_order)-10,max((res_reg1[[k]]^2)^(1/3))-0.08,bquote(~R^2 ==. (round(rsquared_5a[k],3))),cex=0.7)
# max((res_reg1[[k]]^2)^(1/3))-0.08   


# Fit Model 5B 
rho_t_res_reg2 <- cor(wy_order_2[[k]],res_reg1_2_3) 
sd_res_reg1_2_3 <- sd(res_reg1_2_3)
cc0_5b[k] = mean(res_reg1_2_3) - (rho_t_res_reg2 * sd_res_reg1_2_3  * mean(wy_order_2[[k]]) / sd(wy_order_2[[k]]))
cc1_5b[k] = rho_t_res_reg2 * sd_res_reg1_2_3 / sd(wy_order_2[[k]])
res_reg2_5b_fit <- cc0_5b[k] + cc1_5b[k] * wy_order_2[[k]]
res_reg2_5b_var <- res_dof_corr[k]*var(res_reg1_2_3-res_reg2_5b_fit)

# Test alternative conditional max estimation method
res_reg2_5b_trans_mean <- res_dof_corr[k]*(cc0_5b[k] + cc1_5b[k] * mean(wy_order_2[[k]]))^3 + 3*res_reg2_5b_var*(cc0_5b[k] + cc1_5b[k] * mean(wy_order_2[[k]]))

qq <- var(res_reg1_2_3-res_reg2_5b_fit) #0.0315
rr <- 0.5*(res_reg2_5b_trans_mean - mean((res_reg1_2_3-res_reg2_5b_fit)^3)) #0.5*(0.08086177 - 0.00227) = 0.0382209
ss <- (rr + sqrt(qq^3+rr^2))^(1/3)
tt <- sign(rr - sqrt(qq^3+rr^2))*abs((rr - sqrt(qq^3+rr^2)))^(1/3)
ss_tt <- (ss+tt) 

res_reg2_5b_trans_mean_t <- sqrt(((ss + tt) - cc0_5b[k])/cc1_5b[k]) 
  
res_reg2_5b_trans_cond <- cov(wy_order_2[[k]],res_reg1_2)/var(wy_order_2[[k]]) * (max(wy_order_2[[k]]) - res_reg2_5b_trans_mean_t^2)
#cond_var_at_n_5b[k] = res_dof_corr[k]*(res_reg2_5b_trans_mean + res_reg2_5b_trans_cond)

cond_var_5b[[k]] <- (cc0_5b[k] + cc1_5b[k] * wy_order_2[[k]])^3 + 3*res_reg2_5b_var*(cc0_5b[k] + cc1_5b[k] * wy_order_2[[k]]) + mean((res_reg1_2_3-res_reg2_5b_fit)^3) 
cond_var_at_n_5b[k] = res_dof_corr[k]*((cc0_5b[k] + cc1_5b[k] * max(wy_order_2[[k]]))^3 + 3*res_reg2_5b_var*(cc0_5b[k] + cc1_5b[k] * max(wy_order_2[[k]])) + res_dof_corr[k]*mean((res_reg1_2_3-res_reg2_5b_fit)^3))

#plot(wy_order,res_reg1_2_3,col="blue")
#lines(wy_order,res_reg2_5b)
Q_99_at_n_5b[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5b[k])) 
rsquared_5b[k] = 1-sum((res_reg1_2_3 - res_reg2_5b_fit)^2)/sum((res_reg1_2_3 - mean(res_reg1_2_3))^2) 
adj_rsquared_5b[k] = 1- adj_r2_corr[k]*sum((res_reg1_2_3 - res_reg2_5b_fit)^2)/sum((res_reg1_2_3 - mean(res_reg1_2_3))^2) 
rmse_5b[k] = sqrt(1/length(wy_order[[k]])*sum((res_reg1_2_3^3 - cond_var_5b[[k]])^2))
rrmse_5b[k] = sqrt(1/length(wy_order[[k]])*sum(((res_reg1_2_3^3 - cond_var_5b[[k]])/cond_var_5b[[k]])^2))

se_cc0_5b <- sqrt(sum((res_reg1_2_3 - res_reg2_5b_fit)^2)/(length(wy_order_2[[k]]) - 2)*(1/length(wy_order_2[[k]])+mean(wy_order_2[[k]])^2/sum((wy_order_2[[k]]-mean(wy_order_2[[k]]))^2)))
se_cc1_5b <- sqrt(sum((res_reg1_2_3 - res_reg2_5b_fit)^2)/((length(wy_order_2[[k]]) - 2)*sum((wy_order_2[[k]]-mean(wy_order_2[[k]]))^2)))
t_cc0_5b <- cc0_5b[k]/se_cc0_5b 
t_cc1_5b <- cc1_5b[k]/se_cc1_5b 
p_cc0_5b[k] = 2*(pt(-abs(t_cc0_5b), df=length(wy_order[[k]]^2) - 1))
p_cc1_5b[k] = 2*(pt(-abs(t_cc1_5b), df=length(wy_order[[k]]^2) - 1)) 

# Test residual adequacy 
res_reg2_5b <- as.numeric(res_reg1_2_3 - res_reg2_5b_fit)
ppcc_res_reg2 <- ppcc.test(res_reg2_5b)
pval_ppcc_res_reg2_5b[k] = as.numeric(ppcc_res_reg2[2])
dw_res_reg2 <- dwtest(res_reg2_5b ~ wy_order_2[[k]])
pval_dw_res_reg2_5b[k] = as.numeric(dw_res_reg2[2]) 

# Test for heteroscedasticity of residuals of second regression
res_reg2_2_3_5b <- ((res_reg1_2_3 - res_reg2_5b_fit)^2)^(1/3)
res_reg2_2_3_5b.lm <- lm(res_reg2_2_3_5b ~ wy_order_2[[k]])
p_dd1_5b[k] = summary(res_reg2_2_3_5b.lm)$coefficients[2,4]

# Fit Model 5C (exponential)
exp_wy_order_nrmlz <- exp(wy_order[[k]]/max(wy_order[[k]]))
rho_t_res_reg2 <- cor(exp_wy_order_nrmlz,res_reg1_2_3)
sd_res_reg1_2_3 <- sd(res_reg1_2_3)
cc0_5c[k] = mean(res_reg1_2_3) - (rho_t_res_reg2 * sd_res_reg1_2_3  * mean(exp_wy_order_nrmlz) / sd(exp_wy_order_nrmlz))
cc1_5c[k] = rho_t_res_reg2 * sd_res_reg1_2_3 / sd(exp_wy_order_nrmlz)
res_reg2_5c_fit <- cc0_5c[k] + cc1_5c[k] * exp_wy_order_nrmlz
res_reg2_5c_var <- res_dof_corr[k] * var(res_reg1_2_3-res_reg2_5c_fit) 

cond_var_5c[[k]] = (cc0_5c[k] + cc1_5c[k] * exp_wy_order_nrmlz)^3 + 3*res_reg2_5c_var*(cc0_5c[k] + cc1_5c[k] * exp_wy_order_nrmlz) + mean((res_reg1_2_3-res_reg2_5c_fit)^3)
cond_var_at_n_5c[k] = res_dof_corr[k]*((cc0_5c[k] + cc1_5c[k] * max(exp_wy_order_nrmlz))^3 + 3*res_reg2_5c_var*(cc0_5c[k] + cc1_5c[k] * max(exp_wy_order_nrmlz)) + mean((res_reg1_2_3-res_reg2_5c_fit)^3))

Q_99_at_n_5c[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5c[k])) 
rsquared_5c[k] = 1-sum((res_reg1_2_3 - res_reg2_5c_fit)^2)/sum((res_reg1_2_3 - mean(res_reg1_2_3))^2)  
adj_rsquared_5c[k] = 1 - adj_r2_corr[k] * sum((res_reg1_2_3 - res_reg2_5c_fit)^2)/sum((res_reg1_2_3 - mean(res_reg1_2_3))^2)  
rmse_5c[k] = sqrt(1/length(wy_order[[k]])*sum((res_reg1_2_3^3 - cond_var_5c[[k]])^2))
rrmse_5c[k] = sqrt(1/length(wy_order[[k]])*sum(((res_reg1_2_3^3 - cond_var_5c[[k]])/cond_var_5c[[k]])^2))

se_cc0_5c <- sqrt(sum((res_reg1_2_3 - res_reg2_5c_fit)^2)/(length(wy_order[[k]]) - 2)*(1/length(wy_order[[k]])+mean(exp_wy_order_nrmlz)^2/sum((exp_wy_order_nrmlz-mean(exp_wy_order_nrmlz))^2)))
se_cc1_5c <- sqrt(sum((res_reg1_2_3 - res_reg2_5c_fit)^2)/((length(wy_order[[k]]) - 2)*sum((exp_wy_order_nrmlz-mean(exp_wy_order_nrmlz))^2)))
t_cc0_5c <- cc0_5c[k]/se_cc0_5c
t_cc1_5c <- cc1_5c[k]/se_cc1_5c
p_cc0_5c[k] = 2*(pt(-abs(t_cc0_5c), df=length(wy_order[[k]]) - 1))
p_cc1_5c[k] = 2*(pt(-abs(t_cc1_5c), df=length(wy_order[[k]]) - 1)) 

# Test residual adequacy 
res_reg2_5c <- as.numeric(res_reg1_2_3 - res_reg2_5c_fit)
ppcc_res_reg2 <- ppcc.test(res_reg2_5c)
pval_ppcc_res_reg2_5c[k] = as.numeric(ppcc_res_reg2[2])
dw_res_reg2 <- dwtest(res_reg2_5c ~ wy_order[[k]])
pval_dw_res_reg2_5c[k] = as.numeric(dw_res_reg2[2]) 

# Test for heteroscedasticity of residuals of second regression
res_reg2_2_3_5c <- ((res_reg1_2_3 - res_reg2_5c_fit)^2)^(1/3)
res_reg2_2_3_5c.lm <- lm(res_reg2_2_3_5c ~ wy_order[[k]])
p_dd1_5c[k] = summary(res_reg2_2_3_5c.lm)$coefficients[2,4]

# Fit Model 5d (logarithmic)
rho_t_res_reg2 <- cor(log(wy_order[[k]]),res_reg1_2_3)
sd_res_reg1_2_3 <- sd(res_reg1_2_3)
cc0_5d[k] = mean(res_reg1_2_3) - (rho_t_res_reg2 * sd_res_reg1_2_3  * mean(log(wy_order[[k]])) / sd(log(wy_order[[k]])))
cc1_5d[k] = rho_t_res_reg2 * sd_res_reg1_2_3 / sd(log(wy_order[[k]]))
res_reg2_5d_fit <- cc0_5d[k] + cc1_5d[k] * log(wy_order[[k]])
res_reg2_5d_var <- length(wy_order[[k]])/(length(wy_order[[k]])-2)*var(res_reg1_2_3-res_reg2_5d_fit)

cond_var_5d[[k]] = (cc0_5d[k] + cc1_5d[k] * log(wy_order[[k]]))^3 + res_reg2_5d_var*(3*cc0_5d[k] + 3*cc1_5d[k] * log(wy_order[[k]])) + mean((res_reg1_2_3-res_reg2_5d_fit)^3)
cond_var_at_n_5d[k] = res_dof_corr[k]*((cc0_5d[k] + cc1_5d[k] * max(log(wy_order[[k]])))^3 + res_reg2_5d_var*(3*cc0_5d[k] + 3*cc1_5d[k] * max(log(wy_order[[k]]))) + mean((res_reg1_2_3-res_reg2_5d_fit)^3))
  #plot(wy_order,res_reg1_2_3,col="blue")
#lines(wy_order,res_reg2_5d)
Q_99_at_n_5d[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5d[k]))  
rsquared_5d[k] = 1-sum((res_reg1_2_3 - res_reg2_5d_fit)^2)/sum((res_reg1_2_3 - mean(res_reg1_2_3))^2)
adj_rsquared_5d[k] = 1 - adj_r2_corr[k] * sum((res_reg1_2_3 - res_reg2_5d_fit)^2)/sum((res_reg1_2_3 - mean(res_reg1_2_3))^2)
rmse_5d[k] = sqrt(1/length(wy_order[[k]])*sum((res_reg1_2_3 - res_reg2_5d_fit)^2))
rrmse_5d[k] = sqrt(1/length(wy_order[[k]])*sum(((res_reg1_2_3 - res_reg2_5d_fit)/res_reg1_2_3)^2))

se_cc0_5d <- sqrt(sum((res_reg1_2_3 - res_reg2_5d_fit)^2)/(length(wy_order[[k]]) - 2)*(1/length(wy_order[[k]])+mean(log(wy_order[[k]]))^2/sum((log(wy_order[[k]])-mean(log(wy_order[[k]])))^2))) 
se_cc1_5d <- sqrt(sum((res_reg1_2_3 - res_reg2_5d_fit)^2)/((length(wy_order[[k]]) - 2)*sum((log(wy_order[[k]])-mean(log(wy_order[[k]])))^2)))
t_cc0_5d <- cc0_5d[k]/se_cc0_5d
t_cc1_5d <- cc1_5d[k]/se_cc1_5d
p_cc0_5d[k] = 2*(pt(-abs(t_cc0_5d), df=length(wy_order[[k]]) - 1))
p_cc1_5d[k] = 2*(pt(-abs(t_cc1_5d), df=length(wy_order[[k]]) - 1)) 

# Test residual adequacy 
res_reg2_5d <- as.numeric(res_reg1_2_3 - res_reg2_5d_fit)
ppcc_res_reg2 <- ppcc.test(res_reg2_5d)
pval_ppcc_res_reg2_5d[k] = as.numeric(ppcc_res_reg2[2])
dw_res_reg2 <- dwtest(res_reg2_5d ~ wy_order[[k]])
pval_dw_res_reg2_5d[k] = as.numeric(dw_res_reg2[2]) 

# Test for heteroscedasticity of residuals of second regression
res_reg2_2_3_5d <- ((res_reg1_2_3 - res_reg2_5d_fit)^2)^(1/3)
res_reg2_2_3_5d.lm <- lm(res_reg2_2_3_5d ~ wy_order[[k]])
p_dd1_5d[k] = summary(res_reg2_2_3_5d.lm)$coefficients[2,4]

# Model 5E (Method of moments, log-transformed squared residuals)
log_res_reg1_2 <- log(res_reg1[[k]]^2)
rho_t_res_reg2 <- cor(wy_order[[k]],log_res_reg1_2)
sd_log_res_reg1_2 <- sd(log_res_reg1_2)
cc0_5e[k] = mean(log_res_reg1_2) - (rho_t_res_reg2 * sd_log_res_reg1_2  * mean(wy_order[[k]]) / sd(wy_order[[k]]))
cc1_5e[k] = rho_t_res_reg2 * sd_log_res_reg1_2 / sd(wy_order[[k]])
res_reg2_5e_fit <- cc0_5e[k] + cc1_5e[k] * wy_order[[k]]
res_reg2_5e_var <- var(log_res_reg1_2-res_reg2_5e_fit)

cond_var_5e[[k]] <- exp((cc0_5e[k] + cc1_5e[k] * wy_order[[k]] + 0.5*res_reg2_5e_var))
cond_var_at_n_5e[k] = exp(res_dof_corr[k]*(cc0_5e[k] + cc1_5e[k] * max(wy_order[[k]])) + 0.5*res_reg2_5e_var)
#plot(wy_order,log_res_reg1_2,col="blue")
#lines(wy_order,log_res_reg2_5e) 
Q_99_at_n_5e[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5e[k])) 
rsquared_5e[k] = 1-sum((log_res_reg1_2 - res_reg2_5e_fit)^2)/sum((log_res_reg1_2 - mean(log_res_reg1_2))^2) 
adj_rsquared_5e[k] = 1-adj_r2_corr[k]*sum((log_res_reg1_2 - res_reg2_5e_fit)^2)/sum((log_res_reg1_2 - mean(log_res_reg1_2))^2) 
rmse_5e[k] = sqrt(1/length(wy_order[[k]])*sum((res_reg1_2_3^3 - cond_var_5e[[k]])^2)) 
rrmse_5e[k] = sqrt(1/length(wy_order[[k]])*sum(((res_reg1_2_3^3 - cond_var_5e[[k]])/cond_var_5e[[k]])^2))

se_cc0_5e <- sqrt(sum((log_res_reg1_2 - res_reg2_5e_fit)^2)/(length(wy_order[[k]]) - 2)*(1/length(wy_order[[k]])+mean(wy_order[[k]])^2/sum((wy_order[[k]]-mean(wy_order[[k]]))^2))) 
se_cc1_5e <- sqrt(sum((log_res_reg1_2 - res_reg2_5e_fit)^2)/((length(wy_order[[k]]) - 2)*sum((wy_order[[k]]-mean(wy_order[[k]]))^2))) 
t_cc0_5e <- cc0_5e[k]/se_cc0_5e 
t_cc1_5e <- cc1_5e[k]/se_cc1_5e
p_cc0_5e[k] = 2*(pt(-abs(t_cc0_5e), df=length(wy_order[[k]]) - 1)) 
p_cc1_5e[k] = 2*(pt(-abs(t_cc1_5e), df=length(wy_order[[k]]) - 1)) 

# Test residual adequacy 
res_reg2_5e <- as.numeric(log_res_reg1_2 - res_reg2_5e_fit) 
ppcc_res_reg2 <- ppcc.test(res_reg2_5e) 
pval_ppcc_res_reg2_5e[k] = as.numeric(ppcc_res_reg2[2]) 
dw_res_reg2 <- dwtest(res_reg2_5e ~ wy_order[[k]]) 
pval_dw_res_reg2_5e[k] = as.numeric(dw_res_reg2[2]) 

# Test for heteroscedasticity of residuals of second regression
log_res_reg2_5e <- ((log_res_reg1_2 - res_reg2_5e_fit)^2)^(1/3)
log_res_reg2_5e.lm <- lm(log_res_reg2_5e ~ wy_order[[k]])
p_dd1_5e[k] = summary(log_res_reg2_5e.lm)$coefficients[2,4]

# Model 5f: Direct estimation of standard deviation
res_reg1_1_3 <- (res_reg1[[k]]^2)^(1/6)
rho_t_res_reg2 <- cor(wy_order[[k]],res_reg1_1_3)
sd_res_reg1_1_3  <- sd(res_reg1_1_3)
cc0_5f[k] = mean(res_reg1_1_3) - (rho_t_res_reg2 * sd_res_reg1_1_3  * mean(wy_order[[k]]) / sd(wy_order[[k]]))
cc1_5f[k] = rho_t_res_reg2 * sd_res_reg1_1_3 / sd(wy_order[[k]])
res_reg2_5f_fit <- cc0_5f[k] + cc1_5f[k] * wy_order[[k]]
res_reg2_5f_var <- var(res_reg1_1_3-res_reg2_5f_fit)

cond_var_5f[[k]] = (cc0_5f[k] + cc1_5f[k] * wy_order[[k]])^6  
cond_var_at_n_5f[k] = res_dof_corr[k]*(cc0_5f[k] + cc1_5f[k] * max(wy_order[[k]]))^6  
                           
Q_99_at_n_5f[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5f[k])) 
rsquared_5f[k] = 1-sum((res_reg1_1_3 - res_reg2_5f_fit)^2)/sum((res_reg1_1_3 - mean(res_reg1_1_3))^2) 
adj_rsquared_5f[k] = 1-adj_r2_corr[k]*sum((res_reg1_1_3 - res_reg2_5f_fit)^2)/sum((res_reg1_1_3 - mean(res_reg1_1_3))^2) 
rmse_5f[k] = sqrt(1/length(wy_order[[k]])*sum((res_reg1_1_3 - cond_var_5f[[k]])^2))
rrmse_5f[k] = sqrt(1/length(wy_order[[k]])*sum(((res_reg1_1_3 - cond_var_5f[[k]])/cond_var_5f[[k]])^2))

se_cc0_5f <- sqrt(sum((res_reg1_1_3 - res_reg2_5f_fit)^2)/(length(wy_order[[k]]) - 2)*(1/length(wy_order[[k]])+mean(wy_order[[k]])^2/sum((wy_order[[k]]-mean(wy_order[[k]]))^2))) 
se_cc1_5f <- sqrt(sum((res_reg1_1_3 - res_reg2_5f_fit)^2)/((length(wy_order[[k]]) - 2)*sum((wy_order[[k]]-mean(wy_order[[k]]))^2))) 
t_cc0_5f <- cc0_5f[k]/se_cc0_5f 
t_cc1_5f <- cc1_5f[k]/se_cc1_5f
p_cc0_5f[k] = 2*(pt(-abs(t_cc0_5f), df=length(wy_order[[k]]) - 1)) 
p_cc1_5f[k] = 2*(pt(-abs(t_cc1_5f), df=length(wy_order[[k]]) - 1)) 

# Test residual adequacy 
res_reg2_5f <- as.numeric(res_reg1_1_3 - res_reg2_5f_fit) 
ppcc_res_reg2 <- ppcc.test(res_reg2_5f) 
pval_ppcc_res_reg2_5f[k] = as.numeric(ppcc_res_reg2[2]) 
dw_res_reg2 <- dwtest(res_reg2_5f ~ wy_order[[k]]) 
pval_dw_res_reg2_5f[k] = as.numeric(dw_res_reg2[2]) 

# Test for heteroscedasticity of residuals of second regression
res_reg1_1_3_5f <- ((res_reg1_1_3 - res_reg2_5f_fit)^2)^(1/3)
res_reg1_1_3_5f.lm <- lm(res_reg1_1_3_5f ~ wy_order[[k]])
p_dd1_5f[k] = summary(res_reg1_1_3_5f.lm)$coefficients[2,4]


# MODEL 6: ITERATIVE RE(WEIGHTED) LEAST SQUARES --------------------------------------------------

# Fit Model 6A 

# Establish output arrays for each record (length-dependent)
wt_6a <- array(NA,dim=c(length(wy_order[[k]]),10))
res_reg1_iwls_6a <- array(NA,dim=c(length(wy_order[[k]]),10))
res_reg1_iwls_2_3_6a <- array(NA,dim=c(length(wy_order[[k]]),10))
cond_var_6a <- array(NA,dim=c(length(wy_order[[k]]),10))

for (i in 1:10){
  if(i == 1){
    wt_6a[,i] = rep(1,length(wy_order[[k]])) #Make a list
  } else {
    wt_6a[,i] = 1/(cond_var_6a[,i-1]) #Make a list
  }
  flood_reg1_iwls.lm <- lm(log(peak_flow + 0.01) ~ wy_order[[k]], weights=wt_6a[,i]) # Make a list
  b0_reg1_iwls_6a[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[1])
  b1_reg1_iwls_6a[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[2])
  res_reg1_iwls_6a[,i] = as.numeric(residuals(flood_reg1_iwls.lm)) # Make a list
  res_reg1_iwls_2_3_6a[,i] = (as.numeric(residuals(flood_reg1_iwls.lm)^2))^(1/3)# Make a list
  res_reg1_iwls_2_3_6a.lm <- lm(res_reg1_iwls_2_3_6a[,i] ~ wy_order[[k]]) # Make a list
  cc0_6a[k,i] = as.numeric((res_reg1_iwls_2_3_6a.lm$coefficients[1]))
  cc1_6a[k,i] = as.numeric((res_reg1_iwls_2_3_6a.lm$coefficients[2]))
  cond_var_6a[,i] = as.numeric(fitted(res_reg1_iwls_2_3_6a.lm)) # Make this a list
  for (l in 1:length(peak_flow))
  if (cond_var_6a[l,i] < 0.001){
    cond_var_6a[l,i] <- 0.001
  } 
 
} 

wt_6a_list[[k]] = wt_6a[,10]  
res_reg1_iwls_6a_list[[k]] = res_reg1_iwls_6a[,10]  
res_reg1_iwls_2_3_6a_list[[k]] = res_reg1_iwls_2_3_6a[,10]  
cond_var_6a_list[[k]] = cond_var_6a[,10]  
res_reg2_6a_var <- res_dof_corr[k]*var(res_reg1_iwls_2_3_6a[,10]-cond_var_6a[,10]) 

cond_var_at_n_6a[k] = res_dof_corr[k]*((cc0_6a[k,10] + cc1_6a[k,10] * max(wy_order[[k]]))^3 + 3*res_reg2_6a_var*(cc0_6a[k,10] + cc1_6a[k,10] * max(wy_order[[k]])) + mean((res_reg1_2_3-cond_var_6a[,10])^3))
#plot(wy_order,res_reg1_2_3,col="blue")
#lines(wy_order,res_reg2_5a)
Q_99_at_n_6a[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_6a[k]))
rsquared_6a[k] = 1-sum((res_reg1_iwls_2_3_6a[,max(i)] - cond_var_6a[,max(i)])^2)/sum((res_reg1_iwls_2_3_6a[,max(i)] - mean(res_reg1_iwls_2_3_6a[,max(i)]))^2)
adj_rsquared_6a[k] = 1-adj_r2_corr[k]*sum((res_reg1_iwls_2_3_6a[,max(i)] - cond_var_6a[,max(i)])^2)/sum((res_reg1_iwls_2_3_6a[,max(i)] - mean(res_reg1_iwls_2_3_6a[,max(i)]))^2)
rmse_6a[k] = sqrt(1/length(wy_order[[k]])*sum((res_reg1_iwls_2_3_6a[,max(i)] - cond_var_6a[,max(i)])^2))
rrmse_6a[k] = sqrt(1/length(wy_order[[k]])*sum(((res_reg1_iwls_2_3_6a[,max(i)] - cond_var_6a[,max(i)])/res_reg1_iwls_2_3_6a[,max(i)])^2))


se_cc1_6a <- sqrt(sum((res_reg1_iwls_2_3_6a[,max(i)] - cond_var_6a[,max(i)])^2)/((length(wy_order[[k]]) - 2)*sum((wy_order[[k]]-mean(wy_order[[k]]))^2)))
t_cc1_6a <- cc1_6a[k]/se_cc1_6a
p_cc1_6a[k] = 2*(pt(-abs(t_cc1_6a), df=length(wy_order[[k]]) - 1)) 

# Fit model 6B (Quadratic with IWLS)

# Establish output arrays for each record (length-dependent)
wt_6b <- array(NA,dim=c(length(wy_order[[k]]),10))
res_reg1_iwls_6b <- array(NA,dim=c(length(wy_order[[k]]),10))
res_reg1_iwls_2_3_6b <- array(NA,dim=c(length(wy_order[[k]]),10))
cond_var_6b <- array(NA,dim=c(length(wy_order[[k]]),10))

for (i in 1:10){
  if(i == 1){
    wt_6b[,i] = rep(1,length(wy_order[[k]])) #Make a list
  } else {
    wt_6b[,i] = 1/(cond_var_6b[,i-1]) #Make a list
  }
  flood_reg1_iwls.lm <- lm(log(peak_flow + 0.01) ~ wy_order[[k]], weights=wt_6b[,i]) # Make a list
  b0_reg1_iwls_6b[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[1])
  b1_reg1_iwls_6b[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[2])
  res_reg1_iwls_6b[,i] = as.numeric(residuals(flood_reg1_iwls.lm)) # Make a list
  res_reg1_iwls_2_3_6b[,i] = (as.numeric(residuals(flood_reg1_iwls.lm)^2))^(1/3)# Make a list
  res_reg1_iwls_2_3_6b.lm <- lm(res_reg1_iwls_2_3_6b[,i] ~ wy_order_2[[k]]) # Make a list
  cc0_6b[k,i] = as.numeric((res_reg1_iwls_2_3_6b.lm$coefficients[1]))
  cc1_6b[k,i] = as.numeric((res_reg1_iwls_2_3_6b.lm$coefficients[2]))
  cond_var_6b[,i] = as.numeric(fitted(res_reg1_iwls_2_3_6b.lm)) # Make this a list
  for (l in 1:length(peak_flow))
    if (cond_var_6b[l,i] < 0.001){
      cond_var_6b[l,i] <- 0.001
    } 
  
} 

wt_6b_list[[k]] = wt_6b[,10]
res_reg1_iwls_6b_list[[k]] = res_reg1_iwls_6b[,10]
res_reg1_iwls_2_3_6b_list[[k]] = res_reg1_iwls_2_3_6b[,10]
cond_var_6b_list[[k]] = cond_var_6b[,10]
res_reg2_6b_var <- res_dof_corr[k]*var(res_reg1_iwls_2_3_6b[,10]-cond_var_6b[,10])

cond_var_at_n_6b[k] = res_dof_corr[k]*((cc0_6b[k,10] + cc1_6b[k,10] * max(wy_order_2[[k]]))^3 + 3*res_reg2_6b_var*(cc0_6b[k,10] + cc1_6b[k,10] * max(wy_order_2[[k]])) + mean((res_reg1_iwls_2_3_6b[,10]-cond_var_6b[,10])^3))

#plot(wy_order,res_reg1_2_3,col="blue")
#lines(wy_order,res_reg2_5a)
Q_99_at_n_6b[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_6b[k]))
rsquared_6b[k] = 1-sum((res_reg1_iwls_2_3_6b[,max(i)] - cond_var_6b[,max(i)])^2)/sum((res_reg1_iwls_2_3_6b[,max(i)] - mean(res_reg1_iwls_2_3_6b[,max(i)]))^2)
adj_rsquared_6b[k] = 1-adj_r2_corr[k]*sum((res_reg1_iwls_2_3_6b[,max(i)] - cond_var_6b[,max(i)])^2)/sum((res_reg1_iwls_2_3_6b[,max(i)] - mean(res_reg1_iwls_2_3_6b[,max(i)]))^2)
rmse_6b[k] = sqrt(1/length(wy_order[[k]])*sum((res_reg1_iwls_2_3_6b[,max(i)] - cond_var_6b[,max(i)])^2))
rrmse_6b[k] = sqrt(1/length(wy_order[[k]])*sum(((res_reg1_iwls_2_3_6b[,max(i)] - cond_var_6b[,max(i)])/res_reg1_iwls_2_3_6b[,max(i)])^2))

se_cc1_6b <- sqrt(sum((res_reg1_iwls_2_3_6b[,max(i)] - cond_var_6b[,max(i)])^2)/((length(wy_order[[k]]) - 2)*sum((wy_order[[k]]^2-mean(wy_order[[k]]^2))^2)))
t_cc1_6b <- cc1_6b[k]/se_cc1_6b
p_cc1_6b[k] = 2*(pt(-abs(t_cc1_6b), df=length(wy_order[[k]]) - 1)) 

# Fit model 6c (Exponential with IWLS)

# Establish output arrays for each record (length-dependent)
wt_6c <- array(NA,dim=c(length(wy_order[[k]]),10))
res_reg1_iwls_6c <- array(NA,dim=c(length(wy_order[[k]]),10))
res_reg1_iwls_2_3_6c <- array(NA,dim=c(length(wy_order[[k]]),10))
cond_var_6c <- array(NA,dim=c(length(wy_order[[k]]),10))

for (i in 1:10){
  if(i == 1){
    wt_6c[,i] = rep(1,length(wy_order[[k]])) #Make a list
  } else {
    wt_6c[,i] = 1/(cond_var_6c[,i-1]) #Make a list
  }
  flood_reg1_iwls.lm <- lm(log(peak_flow + 0.01) ~ wy_order[[k]], weights=wt_6c[,i]) # Make a list
  b0_reg1_iwls_6c[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[1])
  b1_reg1_iwls_6c[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[2])
  res_reg1_iwls_6c[,i] = as.numeric(residuals(flood_reg1_iwls.lm)) # Make a list
  res_reg1_iwls_2_3_6c[,i] = (as.numeric(residuals(flood_reg1_iwls.lm)^2))^(1/3)# Make a list
  res_reg1_iwls_2_3_6c.lm <- lm(res_reg1_iwls_2_3_6c[,i] ~ exp(wy_order[[k]])) # Make a list
  cc0_6c[k,i] = as.numeric((res_reg1_iwls_2_3_6c.lm$coefficients[1])) 
  cc1_6c[k,i] = as.numeric((res_reg1_iwls_2_3_6c.lm$coefficients[2])) 
  cond_var_6c[,i] = as.numeric(fitted(res_reg1_iwls_2_3_6c.lm))  # Make this a list
  for (l in 1:length(peak_flow))
    if (cond_var_6c[l,i] < 0.001){
      cond_var_6c[l,i] <- 0.001
    } 
  
} 

wt_6c_list[[k]] = wt_6c[,10]
res_reg1_iwls_6c_list[[k]] = res_reg1_iwls_6c[,10]
res_reg1_iwls_2_3_6c_list[[k]] = res_reg1_iwls_2_3_6c[,10]
cond_var_6c_list[[k]] = cond_var_6c[,10]
cond_var_at_n_6c[k] = (exp(cc0_6c[k,10] + cc1_6c[k,10] * max(wy_order[[k]])))^3
#cond_var_at_n_6c[k] = (cc0_6c[k,10] + cc1_6c[k,10] * max(exp(wy_order[[k]])))^3 + 
#  var(res_reg1_iwls_2_3_6c[,10]-cond_var_6c[,10])*(3*cc0_6c[k,10] + 3*cc1_6c[k,10] * max(exp(wy_order[[k]]))) + mean((res_reg1_iwls_2_3_6c[,10]-cond_var_6c[,10])^3)
#plot(wy_order[[k]],res_reg1_2_3,col="blue")
#lines(wy_order[[k]],res_reg2_5a)
Q_99_at_n_6c[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_6c[k]))
rsquared_6c[k] = 1-sum((res_reg1_iwls_2_3_6c[,max(i)] - cond_var_6c[,max(i)])^2)/sum((res_reg1_iwls_2_3_6c[,max(i)] - mean(res_reg1_iwls_2_3_6c[,max(i)]))^2)
adj_rsquared_6c[k] = 1-adj_r2_corr[k]*sum((res_reg1_iwls_2_3_6c[,max(i)] - cond_var_6c[,max(i)])^2)/sum((res_reg1_iwls_2_3_6c[,max(i)] - mean(res_reg1_iwls_2_3_6c[,max(i)]))^2)
rmse_6c[k] = sqrt(1/length(wy_order[[k]])*sum((res_reg1_iwls_2_3_6c[,max(i)] - cond_var_6c[,max(i)])^2))
rrmse_6c[k] = sqrt(1/length(wy_order[[k]])*sum(((res_reg1_iwls_2_3_6c[,max(i)] - cond_var_6c[,max(i)])/res_reg1_iwls_2_3_6a[,max(i)])^2))

se_cc1_6c <- sqrt(sum((res_reg1_iwls_2_3_6c[,max(i)] - cond_var_6c[,max(i)])^2)/((length(wy_order[[k]]) - 2)*sum((exp(wy_order[[k]])-mean(exp(wy_order[[k]])))^2)))
t_cc1_6c <- cc1_6c[k]/se_cc1_6c
p_cc1_6c[k] = 2*(pt(-abs(t_cc1_6c), df=length(wy_order[[k]]) - 1)) 


# Fit model 6d (Logarithmic with IWLS)

# Establish output arrays for each record (length-dependent)
wt_6d <- array(NA,dim=c(length(wy_order[[k]]),10))
res_reg1_iwls_6d <- array(NA,dim=c(length(wy_order[[k]]),10))
res_reg1_iwls_2_3_6d <- array(NA,dim=c(length(wy_order[[k]]),10))
cond_var_6d <- array(NA,dim=c(length(wy_order[[k]]),10))

for (i in 1:10){
  if(i == 1){
    wt_6d[,i] = rep(1,length(wy_order[[k]])) #Make a list
  } else {
    wt_6d[,i] = 1/(cond_var_6d[,i-1]) #Make a list
  }
  flood_reg1_iwls.lm <- lm(log(peak_flow + 0.01) ~ wy_order[[k]], weights=wt_6d[,i]) # Make a list
  b0_reg1_iwls_6d[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[1])
  b1_reg1_iwls_6d[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[2])
  res_reg1_iwls_6d[,i] = as.numeric(residuals(flood_reg1_iwls.lm)) # Make a list
  res_reg1_iwls_2_3_6d[,i] = (as.numeric(residuals(flood_reg1_iwls.lm)^2))^(1/3)# Make a list
  res_reg1_iwls_2_3_6d.lm <- lm(res_reg1_iwls_2_3_6d[,i] ~ log(wy_order[[k]])) # Make a list
  cc0_6d[k,i] = as.numeric((res_reg1_iwls_2_3_6d.lm$coefficients[1])) 
  cc1_6d[k,i] = as.numeric((res_reg1_iwls_2_3_6d.lm$coefficients[2])) 
  cond_var_6d[,i] = as.numeric(fitted(res_reg1_iwls_2_3_6d.lm))  # Make this a list
  for (l in 1:length(peak_flow))
    if (cond_var_6d[l,i] < 0.001){
      cond_var_6d[l,i] <- 0.001
    } 
  
} 

wt_6d_list[[k]] = wt_6d[,10]
res_reg1_iwls_6d_list[[k]] = res_reg1_iwls_6d[,10]
res_reg1_iwls_2_3_6d_list[[k]] = res_reg1_iwls_2_3_6d[,10]
cond_var_6d_list[[k]] = cond_var_6d[,10]
res_reg2_6d_var <- res_dof_corr[k]*var(res_reg1_iwls_2_3_6d[,10]-cond_var_6d[,10])

cond_var_at_n_6d[k] = res_dof_corr[k]*((cc0_6d[k,10] + cc1_6d[k,10] * max(log(wy_order[[k]])))^3 + res_reg2_6d_var*(3*cc0_6d[k,10] + 3*cc1_6d[k,10] * max(log(wy_order[[k]]))) + mean((res_reg1_iwls_2_3_6d[,10]-cond_var_6d[,10])^3))
#plot(wy_order,res_reg1_2_3,col="blue") 
#lines(wy_order,res_reg2_5a)
Q_99_at_n_6d[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_6d[k]))
rsquared_6d[k] = 1-sum((res_reg1_iwls_2_3_6d[,max(i)] - cond_var_6d[,max(i)])^2)/sum((res_reg1_iwls_2_3_6d[,max(i)] - mean(res_reg1_iwls_2_3_6d[,max(i)]))^2)
adj_rsquared_6d[k] = 1-adj_r2_corr[k]*sum((res_reg1_iwls_2_3_6d[,max(i)] - cond_var_6d[,max(i)])^2)/sum((res_reg1_iwls_2_3_6d[,max(i)] - mean(res_reg1_iwls_2_3_6d[,max(i)]))^2)
rmse_6d[k] = sqrt(1/length(wy_order[[k]])*sum((res_reg1_iwls_2_3_6d[,max(i)] - cond_var_6d[,max(i)])^2))
rrmse_6d[k] = sqrt(1/length(wy_order[[k]])*sum(((res_reg1_iwls_2_3_6d[,max(i)] - cond_var_6d[,max(i)])/res_reg1_iwls_2_3_6d[,max(i)])^2))

se_cc1_6d <- sqrt(sum((res_reg1_iwls_2_3_6d[,max(i)] - cond_var_6d[,max(i)])^2)/((length(wy_order[[k]]) - 2)*sum((log(wy_order[[k]])-mean(log(wy_order[[k]])))^2)))
t_cc1_6d <- cc1_6d[k]/se_cc1_6d 
p_cc1_6d[k] = 2*(pt(-abs(t_cc1_6d), df=length(wy_order[[k]]) - 1))  

# Fit Model A using GLM 
glm.cond_var_a <- glm2(formula = res_reg1_2_3 ~ wy_order_1_3, family=gaussian) 
b0_glm_cond_var_a[k] = as.numeric(glm.cond_var_a$coefficients[1]) 
b1_glm_cond_var_a[k] = as.numeric(glm.cond_var_a$coefficients[2]) 
res_dev <- as.numeric(glm.cond_var_a$deviance) 
null_dev <- as.numeric(glm.cond_var_a$null.deviance) 
pseudo_rsquared_glm_a[k] = 1 - res_dev/null_dev 
cond_var_at_n_glm_a[k] = (b0_glm_cond_var_a[k] + b1_glm_cond_var_a[k] * max(wy_order_1_3))^3
Q_99_at_n_glm_a[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_glm_a[k]))

# Fit model using GLM with gamma regression without IWLS
res_reg1_2 <- as.numeric((res_reg1[[k]])^2)
glm.cond_var_gamma <- glm2(formula = res_reg1_2 ~ wy_order[[k]], family=Gamma(link=log))
#Initial conditions: start=c(mean(log(res_reg1_2)),0)
b0_glm_cond_var_gamma[k] = as.numeric(glm.cond_var_gamma$coefficients[1])
b1_glm_cond_var_gamma[k] = as.numeric(glm.cond_var_gamma$coefficients[2])
res_dev <- as.numeric(glm.cond_var_gamma$deviance)
null_dev <- as.numeric(glm.cond_var_gamma$null.deviance) # intercept-only model
pseudo_rsquared_glm_gamma[k] = 1 - res_dev/null_dev
cond_var_at_n_glm_gamma[k] = exp(b0_glm_cond_var_gamma[k] + b1_glm_cond_var_gamma[k] * max(wy_order[[k]])) #Transformation bias?
Q_99_at_n_glm_gamma[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_glm_gamma[k]))

# Fit model using GLM with gamma regression with IWLS
wt_iwls_glm_gamma <- array(NA,dim=c(length(wy_order[[k]]),10))
res_reg1_iwls_glm_gamma <- array(NA,dim=c(length(wy_order[[k]]),10))
res_reg1_2_iwls_glm_gamma <- array(NA,dim=c(length(wy_order[[k]]),10))
cvar_iwls_glm_gamma <- array(NA,dim=c(length(wy_order[[k]]),10))

for (i in 1:10){
  if(i == 1){
    wt_iwls_glm_gamma[,i] = rep(1,length(wy_order[[k]])) #Make a list
  } else {
    wt_iwls_glm_gamma[,i] = 1/(cvar_iwls_glm_gamma[,i-1]) #Make a list
  }
  flood_reg1_iwls.lm <- lm(log(peak_flow + 0.01) ~ wy_order[[k]], weights=wt_iwls_glm_gamma[,i]) # Make a list
  b0_reg1_iwls_glm_gamma[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[1])
  b1_reg1_iwls_glm_gamma[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[2])
  res_reg1_iwls_glm_gamma[,i] = as.numeric(residuals(flood_reg1_iwls.lm)) # Make a list
  res_reg1_2_iwls_glm_gamma[,i] = as.numeric((residuals(flood_reg1_iwls.lm))^2)# Make a list
  res_reg1_2_iwls_glm_gamma.lm <- glm2(res_reg1_2_iwls_glm_gamma[,i] ~ wy_order[[k]],family=Gamma(link=log)) # Make a list
  cc0_iwls_glm_gamma[k,i] = as.numeric((res_reg1_2_iwls_glm_gamma.lm$coefficients[1]))
  cc1_iwls_glm_gamma[k,i] = as.numeric((res_reg1_2_iwls_glm_gamma.lm$coefficients[2]))
  cvar_iwls_glm_gamma[,i] = as.numeric(fitted(res_reg1_2_iwls_glm_gamma.lm)) # Make this a list
  for (l in 1:length(peak_flow))
    if (cvar_iwls_glm_gamma[l,i] < 0.001){
      cvar_iwls_glm_gamma[l,i] <- 0.001
    } 
  
} 

wt_iwls_glm_gamma_list[[k]] = wt_iwls_glm_gamma[,10]  
res_reg1_iwls_glm_gamma_list[[k]] = res_reg1_iwls_glm_gamma[,10]  
res_reg1_2_iwls_glm_gamma_list[[k]] = res_reg1_2_iwls_glm_gamma[,10]  
cvar_iwls_glm_gamma_list[[k]] = cvar_iwls_glm_gamma[,10]  
res_reg2_glm_gamma_var <- res_dof_corr[k]*var(res_reg1_2_iwls_glm_gamma[,10]-cvar_iwls_glm_gamma[,10]) 

cond_var_at_n_iwls_glm_gamma[k] = res_dof_corr[k]*exp((cc0_iwls_glm_gamma[k,10] + cc1_iwls_glm_gamma[k,10] * max(wy_order[[k]]))) 

Q_99_at_n_iwls_glm_gamma[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_iwls_glm_gamma[k])) 
rsquared_iwls_glm_gamma[k] = 1-sum((res_reg1_2_iwls_glm_gamma[,max(i)] - cvar_iwls_glm_gamma[,max(i)])^2)/sum((res_reg1_2_iwls_glm_gamma[,max(i)] - mean(res_reg1_2_iwls_glm_gamma[,max(i)]))^2) 
adj_rsquared_iwls_glm_gamma[k] = 1 -adj_r2_corr[k]*sum((res_reg1_2_iwls_glm_gamma[,max(i)] - cvar_iwls_glm_gamma[,max(i)])^2)/sum((res_reg1_2_iwls_glm_gamma[,max(i)] - mean(res_reg1_2_iwls_glm_gamma[,max(i)]))^2) 
rmse_iwls_glm_gamma[k] = sqrt(1/length(wy_order[[k]])*sum((res_reg1_2_iwls_glm_gamma[,max(i)] - cvar_iwls_glm_gamma[,max(i)])^2))
rrmse_iwls_glm_gamma[k] = sqrt(1/length(wy_order[[k]])*sum(((res_reg1_2_iwls_glm_gamma[,max(i)] - cvar_iwls_glm_gamma[,max(i)])/cvar_iwls_glm_gamma[,max(i)])^2))

se_cc1_iwls_glm_gamma <- sqrt(sum((res_reg1_2_iwls_glm_gamma[,max(i)] - cvar_iwls_glm_gamma[,max(i)])^2)/((length(wy_order[[k]]) - 2)*sum((wy_order[[k]]-mean(wy_order[[k]]))^2))) 
t_cc1_iwls_glm_gamma <- cc1_iwls_glm_gamma[k]/se_cc1_iwls_glm_gamma  
p_cc1_iwls_glm_gamma[k] = 2*(pt(-abs(t_cc1_iwls_glm_gamma), df=length(wy_order[[k]]) - 1)) 

# Fit model using GLM with gamma regression with IWLS and quadratic model (7b)
wt_iwls_glm_gamma_7b <- array(NA,dim=c(length(wy_order[[k]]),10))
res_reg1_iwls_glm_gamma_7b <- array(NA,dim=c(length(wy_order[[k]]),10))
res_reg1_2_iwls_glm_gamma_7b <- array(NA,dim=c(length(wy_order[[k]]),10))
cvar_iwls_glm_gamma_7b <- array(NA,dim=c(length(wy_order[[k]]),10))

for (i in 1:10){
  if(i == 1){
    wt_iwls_glm_gamma_7b[,i] = rep(1,length(wy_order_2[[k]])) #Make a list
  } else {
    wt_iwls_glm_gamma_7b[,i] = 1/(cvar_iwls_glm_gamma_7b[,i-1]) #Make a list
  }
  flood_reg1_iwls_7b.lm <- lm(log(peak_flow + 0.01) ~ wy_order[[k]], weights=wt_iwls_glm_gamma_7b[,i]) # Make a list
  b0_reg1_iwls_glm_gamma_7b[k,i] = as.numeric((flood_reg1_iwls_7b.lm)$coefficients[1])
  b1_reg1_iwls_glm_gamma_7b[k,i] = as.numeric((flood_reg1_iwls_7b.lm)$coefficients[2])
  res_reg1_iwls_glm_gamma_7b[,i] = as.numeric(residuals(flood_reg1_iwls_7b.lm)) # Make a list
  res_reg1_2_iwls_glm_gamma_7b[,i] = as.numeric((residuals(flood_reg1_iwls_7b.lm))^2)# Make a list
  res_reg1_2_iwls_glm_gamma_7b.lm <- glm2(res_reg1_2_iwls_glm_gamma_7b[,i] ~ wy_order_2[[k]],family=Gamma(link=log)) # Make a list
  cc0_iwls_glm_gamma_7b[k,i] = as.numeric((res_reg1_2_iwls_glm_gamma_7b.lm$coefficients[1]))
  cc1_iwls_glm_gamma_7b[k,i] = as.numeric((res_reg1_2_iwls_glm_gamma_7b.lm$coefficients[2]))
  cvar_iwls_glm_gamma_7b[,i] = as.numeric(fitted(res_reg1_2_iwls_glm_gamma_7b.lm)) # Make this a list
  for (l in 1:length(peak_flow))
    if (cvar_iwls_glm_gamma_7b[l,i] < 0.001){
      cvar_iwls_glm_gamma_7b[l,i] <- 0.001
    } 

} 

wt_iwls_glm_gamma_7b_list[[k]] = wt_iwls_glm_gamma_7b[,10]  
res_reg1_iwls_glm_gamma_7b_list[[k]] = res_reg1_iwls_glm_gamma_7b[,10]  
res_reg1_2_iwls_glm_gamma_7b_list[[k]] = res_reg1_2_iwls_glm_gamma_7b[,10]  
cvar_iwls_glm_gamma_7b_list[[k]] = cvar_iwls_glm_gamma_7b[,10]  
res_reg2_glm_gamma_7b_var <- res_dof_corr[k]*var(res_reg1_2_iwls_glm_gamma_7b[,10]-cvar_iwls_glm_gamma_7b[,10]) 

cond_var_at_n_iwls_glm_gamma_7b[k] = res_dof_corr[k]*exp((cc0_iwls_glm_gamma_7b[k,10] + cc1_iwls_glm_gamma_7b[k,10] * max(wy_order_2[[k]]))) 

Q_99_at_n_iwls_glm_gamma_7b[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_iwls_glm_gamma_7b[k])) 
rsquared_iwls_glm_gamma_7b[k] = 1-sum((res_reg1_2_iwls_glm_gamma_7b[,max(i)] - cvar_iwls_glm_gamma_7b[,max(i)])^2)/sum((res_reg1_2_iwls_glm_gamma_7b[,max(i)] - mean(res_reg1_2_iwls_glm_gamma_7b[,max(i)]))^2) 
adj_rsquared_iwls_glm_gamma_7b[k] = 1 -adj_r2_corr[k]*sum((res_reg1_2_iwls_glm_gamma_7b[,max(i)] - cvar_iwls_glm_gamma_7b[,max(i)])^2)/sum((res_reg1_2_iwls_glm_gamma_7b[,max(i)] - mean(res_reg1_2_iwls_glm_gamma_7b[,max(i)]))^2) 
rmse_iwls_glm_gamma_7b[k] = sqrt(1/length(wy_order[[k]])*sum((res_reg1_2_iwls_glm_gamma_7b[,max(i)] - cvar_iwls_glm_gamma_7b[,max(i)])^2))
rrmse_iwls_glm_gamma_7b[k] = sqrt(1/length(wy_order[[k]])*sum(((res_reg1_2_iwls_glm_gamma_7b[,max(i)] - cvar_iwls_glm_gamma_7b[,max(i)])/cvar_iwls_glm_gamma_7b[,max(i)])^2))

se_cc1_iwls_glm_gamma_7b <- sqrt(sum((res_reg1_2_iwls_glm_gamma_7b[,max(i)] - cvar_iwls_glm_gamma_7b[,max(i)])^2)/((length(wy_order_2[[k]]) - 2)*sum((wy_order_2[[k]]-mean(wy_order_2[[k]]))^2))) 
t_cc1_iwls_glm_gamma_7b <- cc1_iwls_glm_gamma_7b[k]/se_cc1_iwls_glm_gamma_7b  
p_cc1_iwls_glm_gamma_7b[k] = 2*(pt(-abs(t_cc1_iwls_glm_gamma_7b), df=length(wy_order_2[[k]]) - 1)) 

# Fit model using GLM with gamma regression with IWLS and logarithmic model
wt_iwls_glm_gamma_7d <- array(NA,dim=c(length(log(wy_order[[k]])),10))
res_reg1_iwls_glm_gamma_7d <- array(NA,dim=c(length(log(wy_order[[k]])),10))
res_reg1_2_iwls_glm_gamma_7d <- array(NA,dim=c(length(log(wy_order[[k]])),10))
cvar_iwls_glm_gamma_7d <- array(NA,dim=c(length(log(wy_order[[k]])),10))

for (i in 1:10){
  if(i == 1){
    wt_iwls_glm_gamma_7d[,i] = rep(1,length(log(wy_order[[k]]))) #Make a list
  } else {
    wt_iwls_glm_gamma_7d[,i] = 1/(cvar_iwls_glm_gamma_7d[,i-1]) #Make a list
  }
  flood_reg1_iwls_7d.lm <- lm(log(peak_flow + 0.01) ~ wy_order[[k]], weights=wt_iwls_glm_gamma_7d[,i]) # Make a list
  b0_reg1_iwls_glm_gamma_7d[k,i] = as.numeric((flood_reg1_iwls_7d.lm)$coefficients[1])
  b1_reg1_iwls_glm_gamma_7d[k,i] = as.numeric((flood_reg1_iwls_7d.lm)$coefficients[2])
  res_reg1_iwls_glm_gamma_7d[,i] = as.numeric(residuals(flood_reg1_iwls_7d.lm)) # Make a list
  res_reg1_2_iwls_glm_gamma_7d[,i] = as.numeric((residuals(flood_reg1_iwls_7d.lm))^2)# Make a list
  res_reg1_2_iwls_glm_gamma_7d.lm <- glm2(res_reg1_2_iwls_glm_gamma_7d[,i] ~ log(wy_order[[k]]),family=Gamma(link=log)) # Make a list
  cc0_iwls_glm_gamma_7d[k,i] = as.numeric((res_reg1_2_iwls_glm_gamma_7d.lm$coefficients[1]))
  cc1_iwls_glm_gamma_7d[k,i] = as.numeric((res_reg1_2_iwls_glm_gamma_7d.lm$coefficients[2]))
  cvar_iwls_glm_gamma_7d[,i] = as.numeric(fitted(res_reg1_2_iwls_glm_gamma_7d.lm)) # Make this a list
  for (l in 1:length(peak_flow))
    if (cvar_iwls_glm_gamma_7d[l,i] < 0.001){
      cvar_iwls_glm_gamma_7d[l,i] <- 0.001
    } 
  
  
} 

wt_iwls_glm_gamma_7d_list[[k]] = wt_iwls_glm_gamma_7d[,10]  
res_reg1_iwls_glm_gamma_7d_list[[k]] = res_reg1_iwls_glm_gamma_7d[,10]  
res_reg1_2_iwls_glm_gamma_7d_list[[k]] = res_reg1_2_iwls_glm_gamma_7d[,10]  
cvar_iwls_glm_gamma_7d_list[[k]] = cvar_iwls_glm_gamma_7d[,10]  
res_reg2_glm_gamma_7d_var <- res_dof_corr[k]*var(res_reg1_2_iwls_glm_gamma_7d[,10]-cvar_iwls_glm_gamma_7d[,10]) 

cond_var_at_n_iwls_glm_gamma_7d[k] = res_dof_corr[k]*exp((cc0_iwls_glm_gamma_7d[k,10] + cc1_iwls_glm_gamma_7d[k,10] * max(log(wy_order[[k]])))) 

Q_99_at_n_iwls_glm_gamma_7d[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_iwls_glm_gamma_7d[k])) 
rsquared_iwls_glm_gamma_7d[k] = 1-sum((res_reg1_2_iwls_glm_gamma_7d[,max(i)] - cvar_iwls_glm_gamma_7d[,max(i)])^2)/sum((res_reg1_2_iwls_glm_gamma_7d[,max(i)] - mean(res_reg1_2_iwls_glm_gamma_7d[,max(i)]))^2) 
adj_rsquared_iwls_glm_gamma_7d[k] = 1 -adj_r2_corr[k]*sum((res_reg1_2_iwls_glm_gamma_7d[,max(i)] - cvar_iwls_glm_gamma_7d[,max(i)])^2)/sum((res_reg1_2_iwls_glm_gamma_7d[,max(i)] - mean(res_reg1_2_iwls_glm_gamma_7d[,max(i)]))^2) 
rmse_iwls_glm_gamma_7d[k] = sqrt(1/length(wy_order[[k]])*sum((res_reg1_2_iwls_glm_gamma_7d[,max(i)] - cvar_iwls_glm_gamma_7d[,max(i)])^2)) 
rrmse_iwls_glm_gamma_7d[k] = sqrt(1/length(wy_order[[k]])*sum(((res_reg1_2_iwls_glm_gamma_7d[,max(i)] - cvar_iwls_glm_gamma_7d[,max(i)])/cvar_iwls_glm_gamma_7d[,max(i)])^2))

se_cc1_iwls_glm_gamma_7d <- sqrt(sum((res_reg1_2_iwls_glm_gamma_7d[,max(i)] - cvar_iwls_glm_gamma_7d[,max(i)])^2)/((length(log(wy_order[[k]])) - 2)*sum((log(wy_order[[k]])-mean(log(wy_order[[k]])))^2))) 
t_cc1_iwls_glm_gamma_7d <- cc1_iwls_glm_gamma_7d[k]/se_cc1_iwls_glm_gamma_7d  
p_cc1_iwls_glm_gamma_7d[k] = 2*(pt(-abs(t_cc1_iwls_glm_gamma_7d), df=length(log(wy_order[[k]])) - 1)) 

# Fit model using GLM with gamma regression with IWLS and "standard deviation" model
wt_iwls_glm_gamma_7f <- array(NA,dim=c(length(wy_order[[k]]),10))
res_reg1_iwls_glm_gamma_7f <- array(NA,dim=c(length(wy_order[[k]]),10))
res_reg1_abs_iwls_glm_gamma_7f <- array(NA,dim=c(length(wy_order[[k]]),10))
cvar_iwls_glm_gamma_7f <- array(NA,dim=c(length(wy_order[[k]]),10))

for (i in 1:10){
  if(i == 1){
    wt_iwls_glm_gamma_7f[,i] = rep(1,length(wy_order[[k]])) #Make a list
  } else {
    wt_iwls_glm_gamma_7f[,i] = 1/(cvar_iwls_glm_gamma_7f[,i-1]) #Make a list
  }
  flood_reg1_iwls_7f.lm <- lm(log(peak_flow + 0.01) ~ wy_order[[k]], weights=wt_iwls_glm_gamma_7f[,i]) # Make a list
  b0_reg1_iwls_glm_gamma_7f[k,i] = as.numeric((flood_reg1_iwls_7f.lm)$coefficients[1])
  b1_reg1_iwls_glm_gamma_7f[k,i] = as.numeric((flood_reg1_iwls_7f.lm)$coefficients[2])
  res_reg1_iwls_glm_gamma_7f[,i] = as.numeric(residuals(flood_reg1_iwls_7f.lm)) # Make a list
  res_reg1_abs_iwls_glm_gamma_7f[,i] = abs(as.numeric(residuals(flood_reg1_iwls_7f.lm))) # Make a list
  res_reg1_abs_iwls_glm_gamma_7f.lm <- glm2(res_reg1_abs_iwls_glm_gamma_7f[,i] ~ wy_order[[k]],family=Gamma(link=log)) # Make a list
  cc0_iwls_glm_gamma_7f[k,i] = as.numeric((res_reg1_abs_iwls_glm_gamma_7f.lm$coefficients[1])) 
  cc1_iwls_glm_gamma_7f[k,i] = as.numeric((res_reg1_abs_iwls_glm_gamma_7f.lm$coefficients[2])) 
  cvar_iwls_glm_gamma_7f[,i] = as.numeric(fitted(res_reg1_abs_iwls_glm_gamma_7f.lm)^2) # Make this a list
  for (l in 1:length(peak_flow))
    if (cvar_iwls_glm_gamma_7f[l,i] < 0.001){
      cvar_iwls_glm_gamma_7f[l,i] <- 0.001
    } 
  
  
} 

wt_iwls_glm_gamma_7f_list[[k]] = wt_iwls_glm_gamma_7f[,10]  
res_reg1_iwls_glm_gamma_7f_list[[k]] = res_reg1_iwls_glm_gamma_7f[,10]  
res_reg1_abs_iwls_glm_gamma_7f_list[[k]] = res_reg1_abs_iwls_glm_gamma_7f[,10] 
cvar_iwls_glm_gamma_7f_list[[k]] = cvar_iwls_glm_gamma_7f[,10]  
res_reg2_glm_gamma_7f_var <- res_dof_corr[k]*var(res_reg1_abs_iwls_glm_gamma_7f[,10] - sqrt(cvar_iwls_glm_gamma_7f[,10]))

cond_var_at_n_iwls_glm_gamma_7f[k] = res_dof_corr[k]*exp(cc0_iwls_glm_gamma_7f[k,10] + cc1_iwls_glm_gamma_7f[k,10] * max(wy_order[[k]]))^2

Q_99_at_n_iwls_glm_gamma_7f[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_iwls_glm_gamma_7f[k])) 
rsquared_iwls_glm_gamma_7f[k] = 1-sum((res_reg1_abs_iwls_glm_gamma_7f[,max(i)] - sqrt(cvar_iwls_glm_gamma_7f[,10])))/sum((res_reg1_abs_iwls_glm_gamma_7f[,max(i)] - mean(res_reg1_abs_iwls_glm_gamma_7f[,max(i)]))^2) 
adj_rsquared_iwls_glm_gamma_7f[k] = 1 -adj_r2_corr[k]*sum((res_reg1_abs_iwls_glm_gamma_7f[,max(i)] - sqrt(cvar_iwls_glm_gamma_7f[,max(i)]))^2)/sum((res_reg1_abs_iwls_glm_gamma_7f[,max(i)] - mean(res_reg1_abs_iwls_glm_gamma_7f[,max(i)]))^2) 
rmse_iwls_glm_gamma_7f[k] = sqrt(1/length(wy_order[[k]])*sum((res_reg1_abs_iwls_glm_gamma_7f[,max(i)] - sqrt(cvar_iwls_glm_gamma_7f[,max(i)]))^2)) 
rrmse_iwls_glm_gamma_7f[k] = sqrt(1/length(wy_order[[k]])*sum(((res_reg1_abs_iwls_glm_gamma_7f[,max(i)] - sqrt(cvar_iwls_glm_gamma_7f[,max(i)]))/sqrt(cvar_iwls_glm_gamma_7f[,max(i)]))^2))

se_cc1_iwls_glm_gamma_7f <- sqrt(sum((res_reg1_abs_iwls_glm_gamma_7f[,max(i)] - sqrt(cvar_iwls_glm_gamma_7f[,max(i)]))^2)/((length(wy_order[[k]]) - 2)*sum((wy_order[[k]]-mean(wy_order[[k]]))^2))) 
t_cc1_iwls_glm_gamma_7f <- cc1_iwls_glm_gamma_7f[k]/se_cc1_iwls_glm_gamma_7f  
p_cc1_iwls_glm_gamma_7f[k] = 2*(pt(-abs(t_cc1_iwls_glm_gamma_7f), df=length(wy_order[[k]]) - 1)) 

# MAIN LOOP ENDS
}

# Compare quantiles
Q99_compar <- c(Q_99_at_n_2sls_0a,
              Q_99_at_n_2sls_0b,
              Q_99_at_n_2sls_1a,
              Q_99_at_n_2sls_1b,
              Q_99_at_n_2a,
              Q_99_at_n_2b,
              Q_99_at_n_5a,
              Q_99_at_n_5b,
              Q_99_at_n_5c,
              Q_99_at_n_5d,
              Q_99_at_n_5e,
              Q_99_at_n_6a,
              Q_99_at_n_6b,
              Q_99_at_n_6c,
              Q_99_at_n_6d,
              Q_99_at_n_glm_gamma,
              Q_99_at_n_iwls_glm_gamma)

# COMPARE R^2(adj) of different models for all stations
par(mfrow=c(1,3))

boxplot(adj_rsquared_5a,
        adj_rsquared_iwls_glm_gamma,
        names=c("OLS","GLM-IWLS"),
        ylab="Adjusted R^2",
        ylim=c(-0.1,0.7),
        main="Linear")

boxplot(adj_rsquared_5b,
        adj_rsquared_iwls_glm_gamma_7b,
        names=c("OLS","GLM-IWLS"),
        ylim=c(-0.1,0.7),
        main="Quadratic")

boxplot(adj_rsquared_5d,
        adj_rsquared_iwls_glm_gamma_7d,
        names=c("OLS","GLM-IWLS"),
        ylim=c(-0.1,0.7),
        main="Logarithmic")

par(mfrow=c(1,1))

# Compute range of estimates for each site for model 5
Q_99_at_n_5 <- matrix(data=c(Q_99_at_n_5a,Q_99_at_n_5b,Q_99_at_n_5d),nrow=3,ncol=length(Q_99_at_n_5a),byrow=TRUE)
apply(Q_99_at_n_5,2,min)/apply(Q_99_at_n_5,2,max)

# Compare R^2(adj) for each site for model 5
adj_rsquared_5 <- matrix(c(adj_rsquared_5a,adj_rsquared_5b,adj_rsquared_5d),nrow=3,ncol=length(Q_99_at_n_5a),byrow=TRUE)
apply(adj_rsquared_5,2,max)-apply(adj_rsquared_5,2,min)

# Which model of model 5 is best? 
order(adj_rsquared_5[,2],decreasing=T)[1]
adj_rsquared_5_best <- apply(adj_rsquared_5,2,function(x)order(x,decreasing=T)[1])

# COmpare OLS and GLM-IWLS R2adj for best models

par(mfrow=c(1,3))

boxplot(adj_rsquared_5a[adj_rsquared_5_best==1],
        adj_rsquared_iwls_glm_gamma[adj_rsquared_5_best==1],
        names=c("OLS","GLM-IWLS"),
        ylab="Adjusted R^2",
        ylim=c(-0.1,0.6),
        main="Linear")

boxplot(adj_rsquared_5b[adj_rsquared_5_best==2],
        adj_rsquared_iwls_glm_gamma_7b[adj_rsquared_5_best==2],
        names=c("OLS","GLM-IWLS"),
        ylim=c(-0.1,0.6),
        main="Quadratic")

boxplot(adj_rsquared_5d[adj_rsquared_5_best==3],
        adj_rsquared_iwls_glm_gamma_7d[adj_rsquared_5_best==3],
        names=c("OLS","GLM-IWLS"),
        ylim=c(-0.1,0.6),
        main="Logarithmic")

par(mfrow=c(1,1))

# Compare linear model quantiles with alternatives
par(mfrow=c(1,2))
plot(Q_99_at_n_5a,Q_99_at_n_5b,log="xy",xlab="100-year flood (Linear)",ylab="100-year flood (Quadratic)",cex.axis=0.85)
lines(seq(0,250000,100),seq(0,250000,100))
plot(Q_99_at_n_5a,Q_99_at_n_5d,log="xy",xlab="100-year flood (Linear)",ylab="100-year flood (Logarithmic)",cex.axis=0.85)
lines(seq(0,250000,100),seq(0,250000,100))
par(mfrow=c(1,1))

# Mean differences
mean((Q_99_at_n_5b - Q_99_at_n_5a)/Q_99_at_n_5a)
mean((Q_99_at_n_5d - Q_99_at_n_5a)/Q_99_at_n_5a)


# Compare OLS-IWLS vs. OLS in terms of bias (%) for each model structure
boxplot((Q_99_at_n_6a-Q_99_at_n_5a)/Q_99_at_n_5a,
        (Q_99_at_n_6b-Q_99_at_n_5b)/Q_99_at_n_5b,
        (Q_99_at_n_6d-Q_99_at_n_5d)/Q_99_at_n_5d,
        names=c("Linear","Quadratic","Logarithmic"),ylab="Bias(%)")

# Compare with stationary estimates and trends in the mean only----------------------

# Make comparison for linear model with histogram
par(mfrow=c(1,2),oma = c(0, 0, 3.5, 0))
hist(Q_99_at_n_5a/Q_99_at_n_no_trend,breaks=seq(1,3,0.2),main="Difference with \n stationary estimate",ylim=c(0,20),xlab="Ratio",cex.main=0.85,cex.axis=0.8)
hist(Q_99_at_n_5a/Q_99_at_n_med_only,breaks=seq(1,3,0.2),main="Difference with \n conditional mean estimate",ylim=c(0,20),xlab="Ratio",cex.main=0.85,cex.axis=0.8)
mtext("Current 100-Year Flood Estimates", outer = TRUE, cex = 1.3)
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

# Make comparison for linear model with 1:1 plot
par(mfrow=c(1,2),mar = c(5,6,3,1.5)+0.1)
plot(Q_99_at_n_no_trend,Q_99_at_n_5a,xlim=c(1000,50000),log="xy",xlab="Stationary Q100",ylab="Nonstationary Q100 \n (Trend in mean and Cv)")  
lines(seq(0,400000,100),seq(0,400000,100)) 
plot(Q_99_at_n_med_only,Q_99_at_n_5a,xlim=c(1000,50000),log="xy",xlab="Q100 for trend in mean only",ylab="Nonstationary Q100 \n (Trend in mean and Cv)") 
lines(seq(0,400000,100),seq(0,400000,100))  
mtext("Increasing trend in mean, \n Decreasing trend in Cv",outer = TRUE, cex = 1.3)
par(mfrow=c(1,1),mar = c(5,4,4,2)+0.1)

# Make comparison for best estimates

# Create vector with best quantile estimates
Q_99_at_n_5_best <- c(Q_99_at_n_5a[adj_rsquared_5_best==1],
                      Q_99_at_n_5b[adj_rsquared_5_best==2],
                      Q_99_at_n_5d[adj_rsquared_5_best==3])

# QUICK AND DIRTY SOLUTION TO ALIGN WITH ORDERING OF Q_99_at_n_5_best 
Q_99_at_n_no_trend_best <- c(Q_99_at_n_no_trend[adj_rsquared_5_best==1],
                             Q_99_at_n_no_trend[adj_rsquared_5_best==2],
                             Q_99_at_n_no_trend[adj_rsquared_5_best==3])

Q_99_at_n_med_only_best <- c(Q_99_at_n_med_only[adj_rsquared_5_best==1], 
                             Q_99_at_n_med_only[adj_rsquared_5_best==2],
                             Q_99_at_n_med_only[adj_rsquared_5_best==3])


# Make comparison for all models
par(mfrow=c(1,2),oma = c(0, 0, 3.5, 0))
hist(Q_99_at_n_5_best/Q_99_at_n_no_trend_best,breaks=seq(1,3,0.2),main="Ratio with \n stationary estimate",ylim=c(0,20),xlab="Ratio",col="gray",cex.main=0.85,cex.axis=0.8)
hist(Q_99_at_n_5_best/Q_99_at_n_med_only_best,breaks=seq(1,3,0.2),main="Ratio with \n conditional mean estimate",ylim=c(0,20),xlab="Ratio",col="gray",cex.main=0.85,cex.axis=0.8)
mtext("Effect of modeling \n increasing trends in variance (N = 36)", outer = TRUE, cex = 1.3)
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

# Compare OLS and IWLS-GLM

par(mfrow=c(1,3))
#plot(Q_99_at_n_5a,Q_99_at_n_iwls_glm_gamma,xlab="OLS Linear",ylab="GLM-IWLS Linear")
plot(Q_99_at_n_5a,Q_99_at_n_iwls_glm_gamma,
     xlab="OLS",ylab="GLM-IWLS",log="xy",xlim=c(100,1000000),ylim=c(100,1000000))
lines(seq(100,1000000,100),seq(100,1000000,100))
title("Linear model") 

plot(Q_99_at_n_5b,Q_99_at_n_iwls_glm_gamma_7b,
     xlab="OLS",ylab="GLM-IWLS",log="xy",xlim=c(100,1000000),ylim=c(100,1000000))
lines(seq(100,1000000,100),seq(100,1000000,100))
title("Quadratic model") 

plot(Q_99_at_n_5d,Q_99_at_n_iwls_glm_gamma_7d,
     xlab="OLS",ylab="GLM-IWLS",log="xy",xlim=c(100,1000000),ylim=c(100,1000000))
lines(seq(100,1000000,100),seq(100,1000000,100))
title("Logarithmic model") 

par(mfrow=c(1,1))

# Median bias
median((Q_99_at_n_5a - Q_99_at_n_iwls_glm_gamma)/Q_99_at_n_iwls_glm_gamma) 
# 0.07
median((Q_99_at_n_5b - Q_99_at_n_iwls_glm_gamma_7b)/Q_99_at_n_iwls_glm_gamma_7b)
# 0.10
median((Q_99_at_n_5d - Q_99_at_n_iwls_glm_gamma_7d)/Q_99_at_n_iwls_glm_gamma_7d)
# 0.06

# Boxplots of bias: OLS vs. GLM-IWLS for all stations
op <- par(mar=c(5, 6, 4, 2) + 0.1)
boxplot((Q_99_at_n_5a - Q_99_at_n_iwls_glm_gamma)/Q_99_at_n_iwls_glm_gamma*100,
        (Q_99_at_n_5b - Q_99_at_n_iwls_glm_gamma_7b)/Q_99_at_n_iwls_glm_gamma_7b*100,
        (Q_99_at_n_5d - Q_99_at_n_iwls_glm_gamma_7d)/Q_99_at_n_iwls_glm_gamma_7d*100,
        names=c("Linear","Quadratic","Logarithmic"),
        ylab="Percent difference \n (OLS - GLM_IWLS)/GLM_IWLS",cex.lab=0.9,cex.axis=0.8)
abline(h=0,lty=2,col="gray")

IQR((Q_99_at_n_5a - Q_99_at_n_iwls_glm_gamma)/Q_99_at_n_iwls_glm_gamma*100)

# Boxplots of bias: OLS vs. GLM-IWLS for stations where models perform best
op <- par(mar=c(5, 6, 4, 2) + 0.1)
boxplot((Q_99_at_n_5a[adj_rsquared_5_best==1] - Q_99_at_n_iwls_glm_gamma[adj_rsquared_5_best==1])/Q_99_at_n_iwls_glm_gamma[adj_rsquared_5_best==1]*100,
        (Q_99_at_n_5b[adj_rsquared_5_best==2] - Q_99_at_n_iwls_glm_gamma_7b[adj_rsquared_5_best==2])/Q_99_at_n_iwls_glm_gamma_7b[adj_rsquared_5_best==2]*100,
        (Q_99_at_n_5d[adj_rsquared_5_best==3] - Q_99_at_n_iwls_glm_gamma_7d[adj_rsquared_5_best==3])/Q_99_at_n_iwls_glm_gamma_7d[adj_rsquared_5_best==3]*100,
        names=c("Linear \n (N = 14)","Quadratic \n (N = 15)","Logarithmic \n (N = 22)"),
        ylab="Percent difference \n (OLS - GLM_IWLS)/GLM_IWLS",cex.lab=0.9,cex.axis=0.8)
abline(h=0,lty=2,col="gray")

median((Q_99_at_n_5a[adj_rsquared_5_best==1] - Q_99_at_n_iwls_glm_gamma[adj_rsquared_5_best==1])/Q_99_at_n_iwls_glm_gamma[adj_rsquared_5_best==1]*100)
median((Q_99_at_n_5b[adj_rsquared_5_best==2] - Q_99_at_n_iwls_glm_gamma_7b[adj_rsquared_5_best==2])/Q_99_at_n_iwls_glm_gamma_7b[adj_rsquared_5_best==2]*100)
median((Q_99_at_n_5d[adj_rsquared_5_best==3] - Q_99_at_n_iwls_glm_gamma_7d[adj_rsquared_5_best==3])/Q_99_at_n_iwls_glm_gamma_7d[adj_rsquared_5_best==3]*100)


# Boxplots of bias: OLS-IWLS vs. GLM-IWLS
boxplot((Q_99_at_n_6a - Q_99_at_n_iwls_glm_gamma)/Q_99_at_n_iwls_glm_gamma*100,
        (Q_99_at_n_6b - Q_99_at_n_iwls_glm_gamma_7b)/Q_99_at_n_iwls_glm_gamma_7b*100,
        (Q_99_at_n_6d - Q_99_at_n_iwls_glm_gamma_7d)/Q_99_at_n_iwls_glm_gamma_7d*100,
        names=c("Linear","Quadratic","Logarithmic"),
        ylab="Percent bias \n (OLS - GLM-IWLS)/GLM-IWLS")
abline(h=0,lty=2,col="gray")

plot(Q_99_at_n_5b,Q_99_at_n_iwls_glm_gamma_7b)
plot(Q_99_at_n_5b,Q_99_at_n_iwls_glm_gamma_7b,xlim=c(0,20000),ylim=c(0,20000),
     xlab="OLS Estimate",
     ylab="IWLS-GLM Estimate") 
lines(seq(0,20000,1000),seq(0,20000,1000))
title("Quadratic model: OLS vs. IWLS-GLM")

plot(Q_99_at_n_5d,Q_99_at_n_iwls_glm_gamma_7d)
plot(Q_99_at_n_5d,Q_99_at_n_iwls_glm_gamma_7d,xlim=c(0,20000),ylim=c(0,20000))
lines(seq(0,20000,1000),seq(0,20000,1000))
title("Logarithmic model: OLS vs. IWLS-GLM")

# Compare RMSE for OLS and IWLS-GLM
par(mfrow=c(1,3))
boxplot(rrmse_5a,rrmse_iwls_glm_gamma,
        names=c("OLS","GLM-IWLS"),main="Linear",
        ylim=c(0,10000),ylab = "Relative RMSE")
boxplot(rrmse_5b,rrmse_iwls_glm_gamma_7b,
        names=c("OLS","GLM-IWLS"),main="Quadratic",
        ylim=c(0,10000),ylab = "Relative RMSE")
boxplot(rrmse_5d,rrmse_iwls_glm_gamma_7d,
        names=c("OLS","GLM-IWLS"),main="Logarithmic",
        ylim=c(0,10000),ylab = "Relative RMSE")
par(mfrow=c(1,1)) 

# Compare RMSE for OLS and IWLS-OLS
par(mfrow=c(1,3))
boxplot(rrmse_5a,rrmse_6a,names=c("OLS","IWLS"),main="Linear")
boxplot(rrmse_5b,rrmse_6b,names=c("OLS","IWLS"),main="Quadratic")
boxplot(rrmse_5d,rrmse_6d,names=c("OLS","IWLS"),main="Logarithmic")
par(mfrow=c(1,1))

# Compare quantiles for OLS and IWLS-OLS
par(mfrow=c(1,3))
boxplot((Q_99_at_n_6a - Q_99_at_n_5a)/Q_99_at_n_5a*100,
        ylim=c(-20,5),
        ylab="Percent Difference (OLS-IWLS vs. OLS)",
        main="Linear")
boxplot((Q_99_at_n_6b - Q_99_at_n_5b)/Q_99_at_n_5b*100,
        ylim=c(-20,5),
        main="Quadratic")
boxplot((Q_99_at_n_6d - Q_99_at_n_5d)/Q_99_at_n_5d*100,
        ylim=c(-20,5),
        main="Logarithmic")
par(mfrow=c(1,1))


# Explore stations without LN2 distribution
station_group_no_ln2 <- station_group[which(station_group$adq_ppcc_ln2==0),]

# Quantile estimate comparison
par(mfrow=c(1,2))
plot(Q_99_at_n_med_only,Q_99_at_n_5a_no_trans_adj,log="xy",xlim=c(100,100000),ylim=c(100,100000),
xlab="Trend in mean only",ylab="Trend in mean and Cv",main="No transf bias \n adjustment",cex.main=0.9)
lines(seq(100,100000,100),seq(100,100000,100))

plot(Q_99_at_n_med_only,Q_99_at_n_5a,log="xy",xlim=c(100,100000),ylim=c(100,100000),
     xlab="Trend in mean only",ylab="Trend in mean and Cv",main="With transf bias \n adjustment",cex.main=0.9)
lines(seq(100,100000,100),seq(100,100000,100))
par(mfrow=c(1,1))

# Compare histograms to evaluate impact of transformation bias 
par(mfrow=c(1,2))
par(mar=c(5.1,5.1,4.1,2.1),mgp=c(4,1,0))
hist(Q_99_at_n_5a_no_trans_adj/Q_99_at_n_med_only,breaks=seq(0.8,1.8,0.1),xlab="100Y flood magnification \n (Trend in mean and Cv/ \n Trend in mean only)",ylim=c(0,10),
     cex.lab=0.8,
     main="No transf bias \n adjustment",
     cex.main=0.9,
     cex.axis=0.85)
hist(Q_99_at_n_5a/Q_99_at_n_med_only,breaks=seq(0.8,1.8,0.1),xlab="100Y flood magnification \n (Trend in mean and Cv/ \n Trend in mean only)",ylim=c(0,10),
     cex.lab=0.8,
     main="With transf bias \n adjustment",
     cex.main=0.9,
     cex.axis=0.85)
par(mfrow=c(1,1))
par(mar=c(5.1,4.1,4.1,2.1),mgp=c(3,1,0))

# Compute flood magnification of mean and Cv
d_mean_stnry <- (Q_99_at_n_med_only-Q_99_at_n_no_trend)/Q_99_at_n_no_trend*100
d_cv_d_mean <- (Q_99_at_n_5a - Q_99_at_n_med_only)/Q_99_at_n_no_trend*100
d_net <- d_mean_stnry + d_cv_d_mean
plot(d_mean_stnry,d_cv_d_mean,
     xlim=c(0,200),ylim=c(-200,200),
     xlab="Percent change in 100-year flood due to trend in mean",
     ylab="Percent change in 100-year flood due to trend in Cv") 
points((Q_99_at_n_med_only-Q_99_at_n_no_trend)/Q_99_at_n_no_trend*100,(Q_99_at_n_5a-Q_99_at_n_med_only)/Q_99_at_n_no_trend*100)
abline(h=0,lty=2)
d_mean_stnry <- (Q_99_at_n_med_only-Q_99_at_n_no_trend)/Q_99_at_n_no_trend*100
d_cv_d_mean <- (Q_99_at_n_5a - Q_99_at_n_med_only)/Q_99_at_n_no_trend*100
d_net <- d_mean_stnry + d_cv_d_mean

x <- seq(0,200,25)
y <- seq(-200,200,25)
z <- matrix(nrow=length(x),ncol=length(y))
for (i in 1:length(x)){
  for (j in 1:length(y)){
    z[i,j] = x[i] + y[j]
  }  
}

contour(x,y,z,        
          levels=c(-100,-50,0,50,100,150,200,250,300,350,400),
          labels = NULL,
          xlim = range(x, finite = TRUE),
          ylim = range(y, finite = TRUE),
          zlim = range(z, finite = TRUE),
          xlab="Corr. trend in mean",ylab="Corr. trend in Cv",main="Linear Model (Cv = 0.5)",cex.main=0.95,
          labcex = 0.75, drawlabels = TRUE, method = "flattest", 
          col = gray.colors(5), lty = par("lty"), lwd = par("lwd"),
          add = FALSE) 
#abline(a = 0, b = -1)


# Compare R2 of +mean,-Cv and + mean,+Cv
boxplot(rsquared_5a_vminus,rsquared_5a,
        names=c(expression(paste("+",mu,", -",Cv)), 
                expression(paste("+",mu,", +",Cv))),
        ylab=expression(R^2))
        
