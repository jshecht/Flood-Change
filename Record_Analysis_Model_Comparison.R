# Must have NWIS_Qpeak.R open while running this sheet

# EXTRACT RECORDS ---------------------------------------------------------------------

# Create vectors and lists to store output

station_group_urb <- rbind(output_lin_mplus_sig05_vminus_sig05_urb,
                            output_lin_mplus_sig05_vplus_sig05_urb)


# Data generated in NWIS_Qpeak

# Conditional mean model
wy_order_urb  <- list()
wy_order_2_urb  <- list()
#wy_order_mean_urb  <- vector(length=nrow(station_group_urb))
#wy_order_max_urb <- vector(length=nrow(station_group_urb))

bhc_ts_TIA_vec_station <- list()
bhc_ts_TIA_vec_station_2 <- list()
bhc_ts_TIA_vec_station_1_3 <- list()

skew_peak_flow_urb  <- vector(length=nrow(station_group_urb))
cv_peak_flow_urb  <- vector(length=nrow(station_group_urb))
skew_ln_peak_flow_urb  <- vector(length=nrow(station_group_urb))

b0_reg1_urb  <- vector(length=nrow(station_group_urb))
b1_reg1_urb  <- vector(length=nrow(station_group_urb))
rho_reg1_urb  <- vector(length=nrow(station_group_urb))
pval_b0_reg1_urb  <- vector(length=nrow(station_group_urb))
pval_b1_reg1_urb  <- vector(length=nrow(station_group_urb))
rsquared_reg1_urb  <- vector(length=nrow(station_group_urb))
# FIX: Look into res_dof_corr
res_dof_corr_urb  <- vector(length=nrow(station_group_urb))
adj_r2_corr_urb  <- vector(length=nrow(station_group_urb))
res_reg1_urb <- list()
res_var_urb  <- vector(length=nrow(station_group_urb))
cond_med_at_n_urb  <- vector(length=nrow(station_group_urb))
cond_mean_at_n_urb <- vector(length=nrow(station_group_urb)) # Accounts for transformation bias
Q_99_at_n_med_only_urb<- vector(length=nrow(station_group_urb))
Q_99_at_n_no_trend_urb <- vector(length=nrow(station_group_urb))

# Type II error analysis 
delta_b1_reg1_true_urb <- vector(length=nrow(station_group_urb))
tt_b1_reg1_urb <- vector(length=nrow(station_group_urb))
t2_error_b1_reg1_urb <- vector(length=nrow(station_group_urb))

# Two-stage least squares - Model A
b0_reg2_0a_urb <- vector(length=nrow(station_group_urb))
b1_reg2_0a_urb <- vector(length=nrow(station_group_urb))
pval_b0_reg2_0a_urb <- vector(length=nrow(station_group_urb))
pval_b1_reg2_0a_urb <- vector(length=nrow(station_group_urb))
adj_rsquared_reg2_0a_urb <- vector(length=nrow(station_group_urb))
res_reg2_0a_urb <- list()
cond_var_at_n_2sls_0a_urb <- vector(length=nrow(station_group_urb)) 
Q_99_at_n_2sls_0a_urb <- vector(length=nrow(station_group_urb)) 

# Two-stage least squares - Model B
b0_reg2_0b_urb <- vector(length=nrow(station_group_urb))
b1_reg2_0b_urb <- vector(length=nrow(station_group_urb))
pval_b0_reg2_0b_urb <- vector(length=nrow(station_group_urb))
pval_b1_reg2_0b_urb <- vector(length=nrow(station_group_urb))
adj_rsquared_reg2_0b_urb <- vector(length=nrow(station_group_urb))
res_reg2_0b_urb <- list()
cond_var_at_n_2sls_0b_urb <- vector(length=nrow(station_group_urb)) 
Q_99_at_n_2sls_0b_urb <- vector(length=nrow(station_group_urb)) 

# Two-stage least squares without Anscombe residuals - Model A
b0_reg2_1a_urb <- vector(length=nrow(station_group_urb))
b1_reg2_1a_urb <- vector(length=nrow(station_group_urb))
pval_b0_reg2_1a_urb <- vector(length=nrow(station_group_urb))
pval_b1_reg2_1a_urb <- vector(length=nrow(station_group_urb))
adj_rsquared_reg2_1a_urb <- vector(length=nrow(station_group_urb))
res_reg2_1a_urb <- list()
cond_var_at_n_2sls_1a_urb <- vector(length=nrow(station_group_urb)) 
Q_99_at_n_2sls_1a_urb <- vector(length=nrow(station_group_urb))

# Two-stage least squares without Anscombe residuals - Model B
b0_reg2_1b_urb <- vector(length=nrow(station_group_urb))
b1_reg2_1b_urb <- vector(length=nrow(station_group_urb))
pval_b0_reg2_1b_urb <- vector(length=nrow(station_group_urb))
pval_b1_reg2_1b_urb <- vector(length=nrow(station_group_urb))
adj_rsquared_reg2_1b_urb <- vector(length=nrow(station_group_urb))
res_reg2_1b_urb <- list()
cond_var_at_n_2sls_1b_urb <- vector(length=nrow(station_group_urb))
Q_99_at_n_2sls_1b_urb <- vector(length=nrow(station_group_urb))


# Two-parameter derived model - Model A 
cond_var_at_n_2a_urb <- vector(length=nrow(station_group_urb))
Q_99_at_n_2a_urb <- vector(length=nrow(station_group_urb))

# Two-parameter derived model - Model B
cond_var_at_n_2b_urb <- vector(length=nrow(station_group_urb)) 
Q_99_at_n_2b_urb <- vector(length=nrow(station_group_urb)) 

# Three-parameter model - Model A
cc1_3a_urb <- vector(length=nrow(station_group_urb))
rsquared_3a_urb <- vector(length=nrow(station_group_urb))
cond_var_at_n_3a_urb <- vector(length=nrow(station_group_urb))
Q_99_at_n_3a_urb <- vector(length=nrow(station_group_urb))

# Three-parameter model - Model C
cc1_3c_urb <- vector(length=nrow(station_group_urb))
rsquared_3c_urb <- vector(length=nrow(station_group_urb))
cond_var_at_n_3c_urb <- vector(length=nrow(station_group_urb))
Q_99_at_n_3c_urb <- vector(length=nrow(station_group_urb))

# Four-parameter - Model A
res_mult1_4a_urb <- list()
cond_var_4a_urb <- list()
cond_var_at_n_4a_urb <- vector(length=nrow(station_group_urb))
Q_99_at_n_4a_urb <- vector(length=nrow(station_group_urb))

# Model 5A (Method of moments, linear trend in variance)
cc0_5a_urb <- vector(length=nrow(station_group_urb))
cc1_5a_urb <- vector(length=nrow(station_group_urb))
rsquared_5a_urb <- vector(length=nrow(station_group_urb))
adj_rsquared_5a_urb <- vector(length=nrow(station_group_urb))
rmse_5a_urb <- vector(length=nrow(station_group_urb))
rrmse_5a_urb <- vector(length=nrow(station_group_urb))
mape_5a_urb <- vector(length=nrow(station_group_urb))
cond_var_5a_urb <- list()
cond_var_at_n_5a_urb <- vector(length=nrow(station_group_urb))
Q_99_at_n_5a_urb <- vector(length=nrow(station_group_urb))
t_cc1_5a_urb <- vector(length=nrow(station_group_urb))
p_cc0_5a_urb <- vector(length=nrow(station_group_urb))
p_cc1_5a_urb <- vector(length=nrow(station_group_urb)) 

zp_RI_if_stnry_5a_urb <- vector(length=nrow(station_group_urb)) 
RI_if_stnry_5a_urb <- vector(length=nrow(station_group_urb)) 

delta_cc1_5a_true_urb <- vector(length=nrow(station_group_urb))
tt_cc1_5a_urb <- vector(length=nrow(station_group_urb))
t2_error_cc1_5a_urb <- vector(length=nrow(station_group_urb))

pval_ppcc_res_reg2_5a_urb <- vector(length=nrow(station_group_urb))
pval_dw_res_reg2_5a_urb <- vector(length=nrow(station_group_urb))
p_dd1_5a_urb <- vector(length=nrow(station_group_urb))

# Model 5B (Method of moments, quadratic trend in variance)
cc0_5b_urb <- vector(length=nrow(station_group_urb))
cc1_5b_urb <- vector(length=nrow(station_group_urb))
rsquared_5b_urb <- vector(length=nrow(station_group_urb))
adj_rsquared_5b_urb <- vector(length=nrow(station_group_urb))
rmse_5b_urb <- vector(length=nrow(station_group_urb))
rrmse_5b_urb <- vector(length=nrow(station_group_urb))

cond_var_5b_urb <- list()
cond_var_at_n_5b_urb <- vector(length=nrow(station_group_urb))
Q_99_at_n_5b_urb <- vector(length=nrow(station_group_urb))
p_cc0_5b_urb <- vector(length=nrow(station_group_urb))
p_cc1_5b_urb <- vector(length=nrow(station_group_urb)) 

pval_ppcc_res_reg2_5b_urb <- vector(length=nrow(station_group_urb))
pval_dw_res_reg2_5b_urb <- vector(length=nrow(station_group_urb))
p_dd1_5b_urb <- vector(length=nrow(station_group_urb))

# Model 5C (Method of moments, log trend in variance)
cc0_5c_urb <- vector(length=nrow(station_group_urb))
cc1_5c_urb <- vector(length=nrow(station_group_urb))
rsquared_5c_urb <- vector(length=nrow(station_group_urb))
adj_rsquared_5c_urb <- vector(length=nrow(station_group_urb))
rmse_5c_urb <- vector(length=nrow(station_group_urb))
rrmse_5c_urb <- vector(length=nrow(station_group_urb))

cond_var_5c_urb <- list()
cond_var_at_n_5c_urb <- vector(length=nrow(station_group_urb))

Q_99_at_n_5c_urb <- vector(length=nrow(station_group_urb))
p_cc0_5_urbc <- vector(length=nrow(station_group_urb))
p_cc1_5c_urb <- vector(length=nrow(station_group_urb)) 

pval_ppcc_res_reg2_5c_urb <- vector(length=nrow(station_group_urb))
pval_dw_res_reg2_5c_urb <- vector(length=nrow(station_group_urb))
p_dd1_5c_urb <- vector(length=nrow(station_group_urb))

# Model 5D (Method of moments, logarithmic trend in variance)
cc0_5d_urb <- vector(length=nrow(station_group_urb))
cc1_5d_urb <- vector(length=nrow(station_group_urb))
rsquared_5d_urb <- vector(length=nrow(station_group_urb))
adj_rsquared_5d_urb <- vector(length=nrow(station_group_urb))
rmse_5d_urb <- vector(length=nrow(station_group_urb))
rrmse_5d_urb <- vector(length=nrow(station_group_urb))

cond_var_5d_urb <- list()
cond_var_at_n_5d_urb <- vector(length=nrow(station_group_urb))

Q_99_at_n_5d_urb <- vector(length=nrow(station_group_urb))
p_cc0_5d_urb <- vector(length=nrow(station_group_urb))
p_cc1_5d_urb <- vector(length=nrow(station_group_urb)) 

pval_ppcc_res_reg2_5d_urb <- vector(length=nrow(station_group_urb))
pval_dw_res_reg2_5d_urb <- vector(length=nrow(station_group_urb))
p_dd1_5d_urb <- vector(length=nrow(station_group_urb))

# Model 5E (Method of moments, log-transformed squared residuals)
cc0_5e_urb <- vector(length=nrow(station_group_urb))
cc1_5e_urb <- vector(length=nrow(station_group_urb))
rsquared_5e_urb <- vector(length=nrow(station_group_urb))
adj_rsquared_5e_urb <- vector(length=nrow(station_group_urb))
rmse_5e_urb <- vector(length=nrow(station_group_urb))
rrmse_5e_urb <- vector(length=nrow(station_group_urb))

cond_var_5e_urb <- list()
cond_var_at_n_5e_urb <- vector(length=nrow(station_group_urb))
Q_99_at_n_5e_urb <- vector(length=nrow(station_group_urb))
p_cc0_5e_urb <- vector(length=nrow(station_group_urb))
p_cc1_5e_urb <- vector(length=nrow(station_group_urb)) 

pval_ppcc_res_reg2_5e_urb <- vector(length=nrow(station_group_urb))
pval_dw_res_reg2_5e_urb <- vector(length=nrow(station_group_urb))
p_dd1_5e_urb <- vector(length=nrow(station_group_urb))

# Model 5F (Method of moments, direct estimation s.d.)
cc0_5f_urb <- vector(length=nrow(station_group_urb))
cc1_5f_urb <- vector(length=nrow(station_group_urb))
rsquared_5f_urb <- vector(length=nrow(station_group_urb))
adj_rsquared_5f_urb <- vector(length=nrow(station_group_urb))
rmse_5f_urb <- vector(length=nrow(station_group_urb))
rrmse_5f_urb <- vector(length=nrow(station_group_urb))

cond_var_5f_urb <- list()
cond_var_at_n_5f_urb <- vector(length=nrow(station_group_urb))
Q_99_at_n_5f_urb <- vector(length=nrow(station_group_urb))
p_cc0_5f_urb <- vector(length=nrow(station_group_urb))
p_cc1_5f_urb <- vector(length=nrow(station_group_urb)) 

pval_ppcc_res_reg2_5f_urb <- vector(length=nrow(station_group_urb))
pval_dw_res_reg2_5f_urb <- vector(length=nrow(station_group_urb))
p_dd1_5f_urb <- vector(length=nrow(station_group_urb))

# Model 6A (IWLS, linear trend in variance)
b0_reg1_iwls_6a_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
b1_reg1_iwls_6a_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
cc0_6a_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
cc1_6a_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)

wt_6a_list_urb <- list()
res_reg1_iwls_6a_list_urb <- list()
res_reg1_iwls_2_3_6a_list_urb <- list()
cond_var_6a_list_urb <- list()

rsquared_6a_urb <- vector(length=nrow(station_group_urb))
adj_rsquared_6a_urb <- vector(length=nrow(station_group_urb))
rmse_6a_urb <- vector(length=nrow(station_group_urb))
rrmse_6a_urb <- vector(length=nrow(station_group_urb))

cond_var_at_n_6a_urb <- vector(length=nrow(station_group_urb))
Q_99_at_n_6a_urb <- vector(length=nrow(station_group_urb))
p_cc1_6a_urb <- vector(length=nrow(station_group_urb)) 

# Model 6B (IWLS, quadratic trend in variance)
b0_reg1_iwls_6b_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
b1_reg1_iwls_6b_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
cc0_6b_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
cc1_6b_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)

wt_6b_list_urb <- list()
res_reg1_iwls_6b_list_urb <- list()
res_reg1_iwls_2_3_6b_list_urb <- list()
cond_var_6b_list_urb <- list()

rsquared_6b_urb <- vector(length=nrow(station_group_urb))
adj_rsquared_6b_urb <- vector(length=nrow(station_group_urb))
rmse_6b_urb <- vector(length=nrow(station_group_urb))
rrmse_6b_urb <- vector(length=nrow(station_group_urb))

cond_var_at_n_6b_urb <- vector(length=nrow(station_group_urb))
Q_99_at_n_6b_urb <- vector(length=nrow(station_group_urb))
p_cc1_6b_urb <- vector(length=nrow(station_group_urb)) 

# Model 6c (IWLS, exponential trend in variance)
b0_reg1_iwls_6c_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
b1_reg1_iwls_6c_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
cc0_6c_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
cc1_6c_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)

wt_6c_list_urb <- list()
res_reg1_iwls_6c_list_urb <- list()
res_reg1_iwls_2_3_6c_list_urb <- list()
cond_var_6c_list_urb <- list()

rsquared_6c_urb <- vector(length=nrow(station_group_urb))
adj_rsquared_6c_urb <- vector(length=nrow(station_group_urb))
rmse_6c_urb <- vector(length=nrow(station_group_urb))
rrmse_6c_urb <- vector(length=nrow(station_group_urb))

cond_var_at_n_6c_urb <- vector(length=nrow(station_group_urb))
Q_99_at_n_6c_urb <- vector(length=nrow(station_group_urb))
p_cc1_6c_urb <- vector(length=nrow(station_group_urb))

# Model 6d (IWLS, logarithmic trend in variance)
b0_reg1_iwls_6d_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
b1_reg1_iwls_6d_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
cc0_6d_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
cc1_6d_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)

wt_6d_list_urb <- list()
res_reg1_iwls_6d_list_urb <- list()
res_reg1_iwls_2_3_6d_list_urb <- list()
cond_var_6d_list_urb <- list()

rsquared_6d_urb <- vector(length=nrow(station_group_urb))
adj_rsquared_6d_urb <- vector(length=nrow(station_group_urb))
rmse_6d_urb <- vector(length=nrow(station_group_urb))
rrmse_6d_urb <- vector(length=nrow(station_group_urb))

cond_var_at_n_6d_urb <- vector(length=nrow(station_group_urb))
Q_99_at_n_6d_urb <- vector(length=nrow(station_group_urb))
p_cc1_6d_urb <- vector(length=nrow(station_group_urb))

# GLM - Model A 
b0_glm_cond_var_a_urb <- vector(length=nrow(station_group_urb))
b1_glm_cond_var_a_urb <- vector(length=nrow(station_group_urb)) 
pseudo_rsquared_glm_a_urb <- vector(length=nrow(station_group_urb)) 
cond_var_at_n_glm_a_urb <- vector(length=nrow(station_group_urb)) 
Q_99_at_n_glm_a_urb <- vector(length=nrow(station_group_urb))

# Gamma GLM model without IWLS
b0_glm_cond_var_gamma_urb <- vector(length=nrow(station_group_urb))
b1_glm_cond_var_gamma_urb <- vector(length=nrow(station_group_urb)) 
pseudo_rsquared_glm_gamma_urb <- vector(length=nrow(station_group_urb)) 
cond_var_at_n_glm_gamma_urb <- vector(length=nrow(station_group_urb)) 
Q_99_at_n_glm_gamma_urb <- vector(length=nrow(station_group_urb))

# Gamma GLM model with IWLS
b0_reg1_iwls_glm_gamma_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
b1_reg1_iwls_glm_gamma_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
cc0_iwls_glm_gamma_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
cc1_iwls_glm_gamma_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
res_reg1_iwls_glm_gamma_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
res_reg1_2_iwls_glm_gamma_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)

wt_iwls_glm_gamma_list_urb <- list()
res_reg1_iwls_glm_gamma_list_urb <- list()
res_reg1_2_iwls_glm_gamma_list_urb <- list()
cvar_iwls_glm_gamma_list_urb <- list()

rsquared_iwls_glm_gamma_urb <- vector(length=nrow(station_group_urb))
adj_rsquared_iwls_glm_gamma_urb <- vector(length=nrow(station_group_urb))
rmse_iwls_glm_gamma_urb <- vector(length=nrow(station_group_urb))
rrmse_iwls_glm_gamma_urb <- vector(length=nrow(station_group_urb))

cond_var_at_n_iwls_glm_gamma_urb <- vector(length=nrow(station_group_urb))
Q_99_at_n_iwls_glm_gamma_urb <- vector(length=nrow(station_group_urb))
p_cc1_iwls_glm_gamma_urb <- vector(length=nrow(station_group_urb)) 

# Gamma GLM model with IWLS - quadratic
b0_reg1_iwls_glm_gamma_7b_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
b1_reg1_iwls_glm_gamma_7b_urb<- matrix(nrow=nrow(station_group_urb),ncol=10)
cc0_iwls_glm_gamma_7b_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
cc1_iwls_glm_gamma_7b_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
res_reg1_iwls_glm_gamma_7b_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
res_reg1_2_iwls_glm_gamma_7b_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)

wt_iwls_glm_gamma_7b_list_urb <- list()
res_reg1_iwls_glm_gamma_7b_list_urb <- list()
res_reg1_2_iwls_glm_gamma_7b_list_urb <- list()
cvar_iwls_glm_gamma_7b_list_urb <- list()

rsquared_iwls_glm_gamma_7b_urb <- vector(length=nrow(station_group_urb))
adj_rsquared_iwls_glm_gamma_7b_urb <- vector(length=nrow(station_group_urb))
rmse_iwls_glm_gamma_7b_urb <- vector(length=nrow(station_group_urb))
rrmse_iwls_glm_gamma_7b_urb <- vector(length=nrow(station_group_urb))

cond_var_at_n_iwls_glm_gamma_7b_urb <- vector(length=nrow(station_group_urb))
Q_99_at_n_iwls_glm_gamma_7b_urb <- vector(length=nrow(station_group_urb))
p_cc1_iwls_glm_gamma_7b_urb <- vector(length=nrow(station_group_urb)) 

# Gamma GLM model with IWLS - logarithmic
b0_reg1_iwls_glm_gamma_7d_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
b1_reg1_iwls_glm_gamma_7d_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
cc0_iwls_glm_gamma_7d_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
cc1_iwls_glm_gamma_7d_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
res_reg1_iwls_glm_gamma_7d_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
res_reg1_2_iwls_glm_gamma_7d_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)

wt_iwls_glm_gamma_7d_list_urb <- list()
res_reg1_iwls_glm_gamma_7d_list_urb <- list()
res_reg1_2_iwls_glm_gamma_7d_list_urb <- list()
cvar_iwls_glm_gamma_7d_list_urb <- list()

rsquared_iwls_glm_gamma_7d_urb <- vector(length=nrow(station_group_urb))
adj_rsquared_iwls_glm_gamma_7d_urb <- vector(length=nrow(station_group_urb))
rmse_iwls_glm_gamma_7d_urb <- vector(length=nrow(station_group_urb))
rrmse_iwls_glm_gamma_7d_urb <- vector(length=nrow(station_group_urb))

cond_var_at_n_iwls_glm_gamma_7d_urb <- vector(length=nrow(station_group_urb))
Q_99_at_n_iwls_glm_gamma_7d_urb <- vector(length=nrow(station_group_urb))
p_cc1_iwls_glm_gamma_7d_urb <- vector(length=nrow(station_group_urb)) 

# Gamma GLM model with IWLS - standard deviation
b0_reg1_iwls_glm_gamma_7f_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
b1_reg1_iwls_glm_gamma_7f_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
cc0_iwls_glm_gamma_7f_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
cc1_iwls_glm_gamma_7f_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
res_reg1_iwls_glm_gamma_7f_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)
res_reg1_abs_iwls_glm_gamma_7f_urb <- matrix(nrow=nrow(station_group_urb),ncol=10)

wt_iwls_glm_gamma_7f_list_urb <- list()
res_reg1_iwls_glm_gamma_7f_list_urb <- list()
res_reg1_abs_iwls_glm_gamma_7f_list_urb <- list()
cvar_iwls_glm_gamma_7f_list_urb <- list()

rsquared_iwls_glm_gamma_7f_urb <- vector(length=nrow(station_group_urb))
adj_rsquared_iwls_glm_gamma_7f_urb <- vector(length=nrow(station_group_urb))
rmse_iwls_glm_gamma_7f_urb <- vector(length=nrow(station_group_urb))
rrmse_iwls_glm_gamma_7f_urb <- vector(length=nrow(station_group_urb))

cond_var_at_n_iwls_glm_gamma_7f_urb <- vector(length=nrow(station_group_urb))
Q_99_at_n_iwls_glm_gamma_7f_urb <- vector(length=nrow(station_group_urb))
p_cc1_iwls_glm_gamma_7f_urb <- vector(length=nrow(station_group_urb)) 

# MAIN LOOP STARTS 

# Loop over stations with +mu,-Cv or +mu,+Cv
for (k in 1:nrow(station_group_urb)){  
  
  # FIX: Import main_list_4 as a table
  # Extract peak flows
  peak_va <- as.numeric(main_list_4[[station_group_urb[k,1]]][,11])/35.3146662127 # Convert from cfs to m^3/s
  
  # Determine water years with peak flow
  wy <- as.numeric(main_list_4[[station_group_urb[k,1]]][,9]) 
  
  # Remove peak flows earlier than 1940 (in cfs)
  wy_pre1940 <- vector(mode="logical",length(wy))
  for(i in 1:length(wy)){
    if(wy[i] < 1940){
      wy_pre1940[i] = TRUE
    } else {
      wy_pre1940[i] = FALSE
    }
  }
  
  # Remove peak flows prior to 1940
  peak_va_urb <- peak_va[wy_pre1940==FALSE]
  
  # Store in subset of urban peak flows for modeled stations (subset of trends)
  peaks_urb[[k]] <- peak_va_urb
  
  # FIX: Remove this unnecessary variable. Kept here due to prior use of peak_flow_urb
  # Extract vector just to use for model
  peak_flow_urb <- peaks_urb[[k]]

  # Create vector of TIA to match water years
  # Determines how many years between 1940 (start of housing data) and 2016 have annual peak flows
  wy_in_1940_2016 <- vector(mode="logical",length=ncol(bhc_ts_ann_mat))
  for(i in 1:77){
    if((1939 + i) %in% wy == TRUE){
      wy_in_1940_2016[i] = TRUE
    } else {
      wy_in_1940_2016[i] = FALSE
    }
  }
  
  # Determine water year order for urban equation, set missing values to zero, filter them out
  wy_order_1940_2016 <- wy_in_1940_2016*seq(1,77,1)
  wy_order_urb[[k]] = wy_order_1940_2016[wy_order_1940_2016>0]
  
  # Find TIA values for years with peak flows
  bhc_ts_TIA_vec <- as.numeric(df_bhc_ts_urb_final_ann[station_group_urb[[k,1]],6:82][wy_in_1940_2016 == TRUE])
  
  # Assign output variable to list, and then compute values raised to 2 and 1/3 powers
  bhc_ts_TIA_vec_station[[k]] = bhc_ts_TIA_vec
  bhc_ts_TIA_vec_station_2[[k]] = bhc_ts_TIA_vec^2
  bhc_ts_TIA_vec_station_1_3[[k]] = bhc_ts_TIA_vec^(1/3)
  
  # COMPUTE LINEAR MODEL OF THE CONDITIONAL MEAN OF LOGS AND BASIC STATS
  
  skew_peak_flow_urb[k] = skewness(peak_flow_urb)
  cv_peak_flow_urb[k] = sd(peak_flow_urb)/mean(peak_flow_urb)
  skew_ln_peak_flow_urb[k] = skewness(log(peak_flow_urb+0.0001))
  flood_reg1_urb.lm = lm(log(peak_flow_urb+0.0001) ~  bhc_ts_TIA_vec_station[[k]]) # Add constant to avoid problems with ephemeral streams
  b0_reg1_urb[k] = summary(flood_reg1_urb.lm)$coefficients[1,1]
  b1_reg1_urb[k] = summary(flood_reg1_urb.lm)$coefficients[2,1]
  rho_reg1_urb[k] = cor(bhc_ts_TIA_vec_station[[k]],log(peak_flow_urb+0.01))
  pval_b0_reg1_urb[k] = summary(flood_reg1_urb.lm)$coefficients[1,4]
  pval_b1_reg1_urb[k] = summary(flood_reg1_urb.lm)$coefficients[2,4]
  rsquared_reg1_urb[k] = as.numeric(summary(flood_reg1_urb.lm)$r.squared)
  res_reg1_urb[[k]] = residuals(flood_reg1_urb.lm)
  res_dof_corr_urb[k] = length(bhc_ts_TIA_vec_station[[k]])/(length(bhc_ts_TIA_vec_station[[k]])-2)
  adj_r2_corr_urb[k] = (length(bhc_ts_TIA_vec_station[[k]]) - 1)/(length(bhc_ts_TIA_vec_station[[k]]) - 2)
  res_var_urb[k] = res_dof_corr_urb[k] * var(res_reg1_urb[[k]])
  cond_med_at_n_urb[k] = b0_reg1_urb[k] + b1_reg1_urb[k] * max(bhc_ts_TIA_vec_station[[k]]) 
  cond_mean_at_n_urb[k] = cond_med_at_n_urb[k] + 0.5 * res_var_urb[k] #Assuming homoscedasticity
  Q_99_at_n_med_only_urb[k] = exp(cond_med_at_n_urb[k] + qnorm(0.99,0,1) * sqrt(res_var_urb[k]))
  Q_99_at_n_no_trend_urb[k] = exp(mean(log(peak_flow_urb + 0.01)) + qnorm(0.99,0,1) * sd(log(peak_flow_urb+0.01)))  
  
  # View conditional mean model 
  plot(bhc_ts_TIA_vec_station[[k]],log(peak_flow_urb+0.01),xlab="TIA (%)",ylab= "ln(Annual peak flow, cfs)") 
  title(c(as.character(station_group_urb$bhc_ts_urb_final.STANAME[k]),
          as.character(station_group_urb$bhc_ts_urb_final.STAID_TEXT[k]),
          "DA",
          as.numeric(station_group_urb$bhc_ts_urb_final.DRAIN_SQKM[k])),
        cex.main=0.8) 
  text(min(bhc_ts_TIA_vec_station[[k]])+1,max(log(peak_flow_urb+0.01))-0.5,bquote(~R^2 ==. (round(rsquared_reg1_urb[k],3))),cex=0.7) 
  
  # View transformed residuals assuming wy_order~res
  plot(bhc_ts_TIA_vec_station[[k]],(res_reg1_urb[[k]]^2)^(1/3),xlab="TIA (%)",ylab= "Residual Variance ^ (2/3)")
  title(c(as.character(station_group_urb$bhc_ts_urb_final.STANAME[k]),
          as.character(station_group_urb$bhc_ts_urb_final.STAID_TEXT[k]),
          "DA",
          as.numeric(station_group_urb$bhc_ts_urb_final.DRAIN_SQKM[k])),
        cex.main=0.8)
  text(min(bhc_ts_TIA_vec_station[[k]])+1,max((res_reg1_urb[[k]]^2)^(1/3))-0.08,bquote(~R^2 ==.(round(rsquared_5a_urb[k],3))),cex=0.7)
  
  # Type II errors using methods from Vogel et al. (2013) and Rosner et al. (2014)
  delta_b1_reg1_true_urb[k] = 1/(sqrt(1/cor(bhc_ts_TIA_vec_station[[k]],log(peak_flow_urb+0.0001))^2-1))
  tt_b1_reg1_urb[k] = abs(qt(1-pval_b1_reg1_urb[k],length(bhc_ts_TIA_vec_station[[k]])-2))
  t2_error_b1_reg1_urb[k] = pt(tt_b1_reg1_urb[k] - delta_b1_reg1_true_urb[k]*sqrt(length(bhc_ts_TIA_vec_station[[k]])),length(bhc_ts_TIA_vec_station[[k]])-2)
  
  
  # View transformed residuals assuming wy_order[[k]]~ln(res^2)
  #log_res_reg1_2 <- log(res_reg1[[k]]^2)
  #plot(wy_order[[k]],log_res_reg1_2,xlab="Water year in record (t)",ylab= "Residual Variance",pch=4)
  #title(c(site_info[station_group_urb[k,1]],
  #        as.character(site_info[station_group_urb[k,1]]),
  #        "DA",
  #        as.numeric(as.character(site_info[station_group_urb[k,1],6]))),cex.main=0.8)
  
  # Fit conditional variance model A using two-stage least squares with an Anscombe transformation
  res_reg1_2_3_urb <- (res_reg1_urb[[k]]^2)^(1/3)
  wy_order_1_3_urb <- wy_order_urb[[k]]^(1/3)
  
  # FIX REPLACE: Models 0-4
 
  # Fit Model 5A 
  rho_t_res_reg2_urb <- cor(bhc_ts_TIA_vec_station[[k]],res_reg1_2_3_urb)
  sd_res_reg1_2_3_urb <- sd(res_reg1_2_3_urb)
  cc0_5a_urb[k] = mean(res_reg1_2_3_urb) - (rho_t_res_reg2_urb * sd_res_reg1_2_3_urb  * mean(bhc_ts_TIA_vec_station[[k]])) / sd(bhc_ts_TIA_vec_station[[k]])
  cc1_5a_urb[k] = rho_t_res_reg2_urb * sd_res_reg1_2_3_urb / sd(bhc_ts_TIA_vec_station[[k]])
  res_reg2_5a_fit_urb <- (cc0_5a_urb[k] + cc1_5a_urb[k] * bhc_ts_TIA_vec_station[[k]])
  res_reg2_5a_var_urb <- res_dof_corr_urb[k]*var(res_reg1_2_3_urb-res_reg2_5a_fit_urb)
  
  cond_var_5a_urb[[k]] <- (cc0_5a_urb[k] + cc1_5a_urb[k] * bhc_ts_TIA_vec_station[[k]])^3 + 3*res_reg2_5a_var_urb*(cc0_5a_urb[k] + cc1_5a_urb[k] * bhc_ts_TIA_vec_station[[k]])
  #cond_var_at_n_5a[k] = res_dof_corr[k]*(cc0_5a[k] + cc1_5a[k] * wy_order_max[k])^3 
  cond_var_at_n_5a_urb[k] = res_dof_corr_urb[k]*(cc0_5a_urb[k] + cc1_5a_urb[k] * max(bhc_ts_TIA_vec_station[[k]]))^3 + 3*res_reg2_5a_var_urb*(cc0_5a_urb[k] + cc1_5a_urb[k] * max(bhc_ts_TIA_vec_station[[k]])) + mean((res_reg1_2_3_urb-res_reg2_5a_fit_urb)^3)
  
  #plot(wy_order,res_reg1_2_3,col="blue")
  #lines(wy_order,res_reg2_5a)
  
  Q_99_at_n_5a_urb[k] = exp(cond_med_at_n_urb[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5a_urb[k])) 
  
  rsquared_5a_urb[k] = 1 - sum((res_reg1_2_3_urb - res_reg2_5a_fit_urb)^2)/sum((res_reg1_2_3_urb - mean(res_reg1_2_3_urb))^2)
  adj_rsquared_5a_urb[k] = 1 - adj_r2_corr_urb[k]*sum((res_reg1_2_3_urb - res_reg2_5a_fit_urb)^2)/sum((res_reg1_2_3_urb - mean(res_reg1_2_3_urb))^2)
  rmse_5a_urb[k] = sqrt(1/length(bhc_ts_TIA_vec_station[[k]])*sum((res_reg1_2_3_urb^3 - cond_var_at_n_5a_urb[[k]])^2))
  rrmse_5a_urb[k] = sqrt(1/length(bhc_ts_TIA_vec_station[[k]])*sum(((res_reg1_2_3_urb^3 - cond_var_at_n_5a_urb[[k]])/cond_var_at_n_5a_urb[[k]])^2))
  mape_5a_urb[k] = 100 * 1/length(bhc_ts_TIA_vec_station[[k]])*sum(abs(res_reg1_2_3_urb^3 - cond_var_at_n_5a_urb[[k]])/cond_var_at_n_5a_urb[[k]])
  se_cc0_5a_urb <- sqrt(sum((res_reg1_2_3_urb - res_reg2_5a_fit_urb)^2)/(length(bhc_ts_TIA_vec_station[[k]]) - 2)*(1/length(bhc_ts_TIA_vec_station[[k]])+mean(bhc_ts_TIA_vec_station[[k]])^2/sum((bhc_ts_TIA_vec_station[[k]]-mean(bhc_ts_TIA_vec_station[[k]]))^2))) 
  se_cc1_5a_urb <- sqrt(sum((res_reg1_2_3_urb - res_reg2_5a_fit_urb)^2)/((length(bhc_ts_TIA_vec_station[[k]]) - 2)*sum((bhc_ts_TIA_vec_station[[k]]-mean(bhc_ts_TIA_vec_station[[k]]))^2)))
  t_cc0_5a_urb <- cc0_5a_urb[k]/se_cc0_5a_urb 
  t_cc1_5a_urb[k] = cc1_5a_urb[k]/se_cc1_5a_urb
  p_cc0_5a_urb[k] = 2*(pt(-abs(t_cc0_5a_urb), df=length(bhc_ts_TIA_vec_station[[k]]) - 1))
  p_cc1_5a_urb[k] = 2*(pt(-abs(t_cc1_5a_urb[k]), df=length(bhc_ts_TIA_vec_station[[k]]) - 1)) 
  #p_cc1_5a[k] = 1 - 2*pt(t_cc1_5a, df=length(wy_order[[k]]) - 1))
  
  # Compute return period difference
  # COmpute stationary 100-year flood
  zp_RI_if_stnry_5a_urb[k] = (log(Q_99_at_n_no_trend_urb[k])-cond_med_at_n_urb[k])/sqrt(cond_var_at_n_5a_urb[k])
  RI_if_stnry_5a_urb[k] = 1/(1-pnorm(zp_RI_if_stnry_5a_urb[k]))
  
  # Type II errors using methods from Vogel et al. (2013) and Rosner et al. (2014)
  delta_cc1_5a_true_urb[k] = 1/(sqrt(1/cor(bhc_ts_TIA_vec_station[[k]],res_reg1_2_3_urb)^2-1)) 
  tt_cc1_5a_urb[k] = qt(1-p_cc1_5a_urb[k],length(bhc_ts_TIA_vec_station[[k]])-2)
  t2_error_cc1_5a_urb[k] = pt(tt_cc1_5a_urb[k] - delta_cc1_5a_true_urb[k]*sqrt(length(bhc_ts_TIA_vec_station[[k]])),length(bhc_ts_TIA_vec_station[[k]])-2)
  
  # Test residual adequacy 
  res_reg2_5a_urb <- as.numeric(res_reg1_2_3_urb - res_reg2_5a_fit_urb)
  ppcc_res_reg2_urb <- ppcc.test(res_reg2_5a_urb)
  pval_ppcc_res_reg2_5a_urb[k] = as.numeric(ppcc_res_reg2_urb[2])
  dw_res_reg2_urb <- dwtest(res_reg2_5a_urb ~ bhc_ts_TIA_vec_station[[k]])
  pval_dw_res_reg2_5a_urb[k] = as.numeric(dw_res_reg2_urb$p.value) 
  
  # Test for heteroscedasticity of residuals of second regression
  res_reg2_2_3_5a_urb <- ((res_reg1_2_3_urb - res_reg2_5a_fit_urb)^2)^(1/3)
  res_reg2_2_3_5a_urb.lm <- lm(res_reg2_2_3_5a_urb ~ bhc_ts_TIA_vec_station[[k]])
  p_dd1_5a_urb[k] = summary(res_reg2_2_3_5a_urb.lm)$coefficients[2,4]
  
  # View transformed residuals assuming wy_order~res
  #plot(wy_order,(res_reg1[[k]]^2)^(1/3),xlab="Water year in record (t)",ylab= "Residual Variance ^ (2/3)")
  #title(c(site_names[station_group_urb[k,1]],
  #      as.character(site_id[station_group_urb[k,1]]),
  #      "DA",
  #      as.numeric(as.character(site_info[station_group_urb[k,1],6]))),cex.main=0.8) 
  #text(max(wy_order)-10,max((res_reg1[[k]]^2)^(1/3))-0.08,bquote(~R^2 ==. (round(rsquared_5a[k],3))),cex=0.7)
  # max((res_reg1[[k]]^2)^(1/3))-0.08   
  
  
  # Fit Model 5B 
  rho_t_res_reg2_urb <- cor(bhc_ts_TIA_vec_station_2[[k]],res_reg1_2_3_urb) 
  sd_res_reg1_2_3_urb <- sd(res_reg1_2_3_urb)
  cc0_5b_urb[k] = mean(res_reg1_2_3_urb) - (rho_t_res_reg2_urb * sd_res_reg1_2_3_urb * mean(bhc_ts_TIA_vec_station_2[[k]]) / sd(bhc_ts_TIA_vec_station_2[[k]]))
  cc1_5b_urb[k] = rho_t_res_reg2_urb * sd_res_reg1_2_3_urb / sd(bhc_ts_TIA_vec_station_2[[k]])
  res_reg2_5b_fit_urb <- cc0_5b_urb[k] + cc1_5b_urb[k] * bhc_ts_TIA_vec_station_2[[k]]
  res_reg2_5b_var_urb <- res_dof_corr_urb[k]*var(res_reg1_2_3_urb - res_reg2_5b_fit_urb)
  
  # Test alternative conditional max estimation method
  res_reg2_5b_trans_mean_urb <- res_dof_corr_urb[k]*(cc0_5b_urb[k] + cc1_5b_urb[k] * mean(bhc_ts_TIA_vec_station_2[[k]]))^3 + 3*res_reg2_5b_var_urb*(cc0_5b_urb[k] + cc1_5b_urb[k] * mean(bhc_ts_TIA_vec_station_2[[k]]))
  
  #qq_urb <- var(res_reg1_2_3_urb - res_reg2_5b_fit_urb)
  #rr_urb <- 0.5*(res_reg2_5b_trans_mean_urb - mean((res_reg1_2_3_urb-res_reg2_5b_fit_urb)^3)) #0.5*(0.08086177 - 0.00227) = 0.0382209
  #ss_urb <- (rr_urb + sqrt(qq_urb^3+rr_urb^2))^(1/3)
  #tt_urb <- sign(rr_urb - sqrt(qq_urb^3+rr_urb^2))*abs((rr_urb - sqrt(qq_urb^3+rr_urb^2)))^(1/3)
  #ss_tt_urb <- (ss_urb+tt_urb) 
  
  #res_reg2_5b_trans_mean_urb <- sqrt(((ss_urb  + tt_urb ) - cc0_5b_urb[k])/cc1_5b_urb[k]) 
  
  #res_reg2_5b_trans_cond_urb <- cov(wy_order_2_urb[[k]],res_reg1_2_urb)/var(wy_order_2_urb[[k]]) * (max(wy_order_2_urb[[k]]) - res_reg2_5b_trans_mean_urb^2)
  #cond_var_at_n_5b[k] = res_dof_corr[k]*(res_reg2_5b_trans_mean + res_reg2_5b_trans_cond)
  
  cond_var_5b_urb[[k]] <- (cc0_5b_urb[k] + cc1_5b_urb[k] * bhc_ts_TIA_vec_station_2[[k]])^3 + 3*res_reg2_5b_var_urb*(cc0_5b_urb[k] + cc1_5b_urb[k] * bhc_ts_TIA_vec_station_2[[k]]) + mean((res_reg1_2_3_urb - res_reg2_5b_fit_urb)^3) 
  cond_var_at_n_5b_urb[k] = res_dof_corr_urb[k]*((cc0_5b_urb[k] + cc1_5b_urb[k] * max(bhc_ts_TIA_vec_station_2[[k]]))^3 + 3*res_reg2_5b_var_urb*(cc0_5b_urb[k] + cc1_5b_urb[k] * max(bhc_ts_TIA_vec_station_2[[k]])) + res_dof_corr_urb[k]*mean((res_reg1_2_3_urb - res_reg2_5b_fit_urb)^3))
  
  #plot(wy_order,res_reg1_2_3,col="blue")
  #lines(wy_order,res_reg2_5b)
  Q_99_at_n_5b_urb[k] = exp(cond_med_at_n_urb[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5b_urb[k])) 
  rsquared_5b_urb[k] = 1-sum((res_reg1_2_3_urb - res_reg2_5b_fit_urb)^2)/sum((res_reg1_2_3_urb - mean(res_reg1_2_3_urb))^2) 
  adj_rsquared_5b_urb[k] = 1- adj_r2_corr_urb[k]*sum((res_reg1_2_3_urb - res_reg2_5b_fit_urb)^2)/sum((res_reg1_2_3_urb - mean(res_reg1_2_3_urb))^2) 
  rmse_5b_urb[k] = sqrt(1/length(bhc_ts_TIA_vec_station_2[[k]])*sum((res_reg1_2_3_urb^3 - cond_var_5b_urb[[k]])^2))
  rrmse_5b_urb[k] = sqrt(1/length(bhc_ts_TIA_vec_station_2[[k]])*sum(((res_reg1_2_3_urb^3 - cond_var_5b_urb[[k]])/cond_var_5b_urb[[k]])^2))
  
  se_cc0_5b_urb <- sqrt(sum((res_reg1_2_3_urb - res_reg2_5b_fit_urb)^2)/(length(bhc_ts_TIA_vec_station_2[[k]]) - 2)*(1/length(bhc_ts_TIA_vec_station_2[[k]])+mean(bhc_ts_TIA_vec_station_2[[k]])^2/sum((bhc_ts_TIA_vec_station_2[[k]]-mean(bhc_ts_TIA_vec_station_2[[k]]))^2)))
  se_cc1_5b_urb <- sqrt(sum((res_reg1_2_3_urb - res_reg2_5b_fit_urb)^2)/((length(bhc_ts_TIA_vec_station_2[[k]]) - 2)*sum((bhc_ts_TIA_vec_station_2[[k]]-mean(bhc_ts_TIA_vec_station_2[[k]]))^2)))
  t_cc0_5b_urb <- cc0_5b_urb[k]/se_cc0_5b_urb 
  t_cc1_5b_urb <- cc1_5b_urb[k]/se_cc1_5b_urb
  p_cc0_5b_urb[k] = 2*(pt(-abs(t_cc0_5b_urb), df=length(bhc_ts_TIA_vec_station_2[[k]]) - 1))
  p_cc1_5b_urb[k] = 2*(pt(-abs(t_cc1_5b_urb), df=length(bhc_ts_TIA_vec_station_2[[k]]) - 1)) 
  
  # Test residual adequacy 
  res_reg2_5b_urb <- as.numeric(res_reg1_2_3_urb - res_reg2_5b_fit_urb)
  ppcc_res_reg2_urb <- ppcc.test(res_reg2_5b_urb)
  pval_ppcc_res_reg2_5b_urb[k] = as.numeric(ppcc_res_reg2_urb[2])
  dw_res_reg2_urb <- dwtest(res_reg2_5b_urb ~ bhc_ts_TIA_vec_station_2[[k]])
  pval_dw_res_reg2_5b_urb[k] = as.numeric(dw_res_reg2_urb$p.value) 
  
  # Test for heteroscedasticity of residuals of second regression
  res_reg2_2_3_5b_urb <- ((res_reg1_2_3_urb - res_reg2_5b_fit_urb)^2)^(1/3)
  res_reg2_2_3_5b_urb.lm <- lm(res_reg2_2_3_5b_urb ~ bhc_ts_TIA_vec_station_2[[k]])
  p_dd1_5b_urb[k] = summary(res_reg2_2_3_5b_urb.lm)$coefficients[2,4]
  
  # Fit Model 5C (exponential)
  '
  exp_wy_order_nrmlz_time <- exp(wy_order_time[[k]]/max(wy_order_time[[k]]))
  rho_t_res_reg2_time <- cor(exp_wy_order_nrmlz_time,res_reg1_2_3_time)
  sd_res_reg1_2_3_time <- sd(res_reg1_2_3_time)
  cc0_5c_time[k] = mean(res_reg1_2_3_time) - (rho_t_res_reg2_time * sd_res_reg1_2_3_time  * mean(exp_wy_order_nrmlz_time) / sd(exp_wy_order_nrmlz_time))
  cc1_5c_time[k] = rho_t_res_reg2_time * sd_res_reg1_2_3_time / sd(exp_wy_order_nrmlz_time)
  res_reg2_5c_fit_time <- cc0_5c_time[k] + cc1_5c_time[k] * exp_wy_order_nrmlz_time
  res_reg2_5c_var_time <- res_dof_corr_time[k] * var(res_reg1_2_3_time-res_reg2_5c_fit_time) 
  
  cond_var_5c_time[[k]] = (cc0_5c_time[k] + cc1_5c_time[k] * exp_wy_order_nrmlz_time)^3 + 3*res_reg2_5c_var_time*(cc0_5c_time[k] + cc1_5c_time[k] * exp_wy_order_nrmlz_time) + mean((res_reg1_2_3_time-res_reg2_5c_fit_time)^3)
  cond_var_at_n_5c_time[k] = res_dof_corr_time[k]*((cc0_5c_time[k] + cc1_5c_time[k] * max(exp_wy_order_nrmlz_time))^3 + 3*res_reg2_5c_var_time*(cc0_5c_time[k] + cc1_5c_time[k] * max(exp_wy_order_nrmlz_time)) + mean((res_reg1_2_3_time-res_reg2_5c_fit_time)^3))
  
  Q_99_at_n_5c_time[k] = exp(cond_med_at_n_time[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5c_time[k])) 
  rsquared_5c_time[k] = 1-sum((res_reg1_2_3_time - res_reg2_5c_fit_time)^2)/sum((res_reg1_2_3_time - mean(res_reg1_2_3_time))^2)  
  adj_rsquared_5c_time[k] = 1 - adj_r2_corr_time[k] * sum((res_reg1_2_3_time - res_reg2_5c_fit_time)^2)/sum((res_reg1_2_3_time - mean(res_reg1_2_3_time))^2)  
  rmse_5c_time[k] = sqrt(1/length(wy_order_time[[k]])*sum((res_reg1_2_3_time^3 - cond_var_5c_time[[k]])^2))
  rrmse_5c_time[k] = sqrt(1/length(wy_order_time[[k]])*sum(((res_reg1_2_3_time^3 - cond_var_5c_time[[k]])/cond_var_5c_time[[k]])^2))
  
  se_cc0_5c_time <- sqrt(sum((res_reg1_2_3_time - res_reg2_5c_fit_time)^2)/(length(wy_order_time[[k]]) - 2)*(1/length(wy_order_time[[k]])+mean(exp_wy_order_nrmlz_time)^2/sum((exp_wy_order_nrmlz_time-mean(exp_wy_order_nrmlz_time))^2)))
  se_cc1_5c_time <- sqrt(sum((res_reg1_2_3_time - res_reg2_5c_fit_time)^2)/((length(wy_order_time[[k]]) - 2)*sum((exp_wy_order_nrmlz_time-mean(exp_wy_order_nrmlz_time))^2)))
  t_cc0_5c_time <- cc0_5c_time[k]/se_cc0_5c_time
  t_cc1_5c_time <- cc1_5c_time[k]/se_cc1_5c_time
  #p_cc0_5c_time[k] = 2*(pt(-abs(t_cc0_5c_time), df=length(wy_order_time[[k]]) - 1))
  p_cc1_5c_time[k] = 2*(pt(-abs(t_cc1_5c_time), df=length(wy_order_time[[k]]) - 1)) 
  
  # Test residual adequacy 
  res_reg2_5c_time <- as.numeric(res_reg1_2_3_time - res_reg2_5c_fit_time)
  ppcc_res_reg2_time <- ppcc.test(res_reg2_5c_time)
  pval_ppcc_res_reg2_5c_time[k] = as.numeric(ppcc_res_reg2_time[2])
  dw_res_reg2_time <- dwtest(res_reg2_5c_time ~ wy_order_time[[k]])
  pval_dw_res_reg2_5c_time[k] = as.numeric(dw_res_reg2_time$p.value) 
  
  # Test for heteroscedasticity of residuals of second regression
  res_reg2_2_3_5c_time <- ((res_reg1_2_3_time - res_reg2_5c_fit_time)^2)^(1/3)
  res_reg2_2_3_5c_time.lm <- lm(res_reg2_2_3_5c_time ~ wy_order_time[[k]])
  p_dd1_5c_time[k] = summary(res_reg2_2_3_5c_time.lm)$coefficients[2,4]
  '
  
  # Fit Model 5d (logarithmic)
  rho_t_res_reg2_urb <- cor(log(bhc_ts_TIA_vec_station[[k]]),res_reg1_2_3_urb)
  sd_res_reg1_2_3_urb <- sd(res_reg1_2_3_urb)
  cc0_5d_urb[k] = mean(res_reg1_2_3_urb) - (rho_t_res_reg2_urb * sd_res_reg1_2_3_urb  * mean(log(bhc_ts_TIA_vec_station[[k]])) / sd(log(bhc_ts_TIA_vec_station[[k]])))
  cc1_5d_urb[k] = rho_t_res_reg2_urb * sd_res_reg1_2_3_urb / sd(log(bhc_ts_TIA_vec_station[[k]]))
  res_reg2_5d_fit_urb <- cc0_5d_urb[k] + cc1_5d_urb[k] * log(bhc_ts_TIA_vec_station[[k]])
  res_reg2_5d_var_urb <- length(bhc_ts_TIA_vec_station[[k]])/(length(bhc_ts_TIA_vec_station[[k]])-2)*var(res_reg1_2_3_urb-res_reg2_5d_fit_urb)
  
  cond_var_5d_urb[[k]] = (cc0_5d_urb[k] + cc1_5d_urb[k] * log(bhc_ts_TIA_vec_station[[k]]))^3 + res_reg2_5d_var_urb*(3*cc0_5d_urb[k] + 3*cc1_5d_urb[k] * log(bhc_ts_TIA_vec_station[[k]])) + mean((res_reg1_2_3_urb-res_reg2_5d_fit_urb)^3)
  cond_var_at_n_5d_urb[k] = res_dof_corr_urb[k]*((cc0_5d_urb[k] + cc1_5d_urb[k] * max(log(bhc_ts_TIA_vec_station[[k]])))^3 + res_reg2_5d_var_urb*(3*cc0_5d_urb[k] + 3*cc1_5d_urb[k] * max(log(bhc_ts_TIA_vec_station[[k]]))) + mean((res_reg1_2_3_urb-res_reg2_5d_fit_urb)^3))
  #plot(wy_order,res_reg1_2_3,col="blue")
  #lines(wy_order,res_reg2_5d)
  Q_99_at_n_5d_urb[k] = exp(cond_med_at_n_urb[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5d_urb[k]))  
  rsquared_5d_urb[k] = 1-sum((res_reg1_2_3_urb - res_reg2_5d_fit_urb)^2)/sum((res_reg1_2_3_urb - mean(res_reg1_2_3_urb))^2)
  adj_rsquared_5d_urb[k] = 1 - adj_r2_corr_urb[k] * sum((res_reg1_2_3_urb - res_reg2_5d_fit_urb)^2)/sum((res_reg1_2_3_urb - mean(res_reg1_2_3_urb))^2)
  rmse_5d_urb[k] = sqrt(1/length(bhc_ts_TIA_vec_station[[k]])*sum((res_reg1_2_3_urb - res_reg2_5d_fit_urb)^2))
  rrmse_5d_urb[k] = sqrt(1/length(bhc_ts_TIA_vec_station[[k]])*sum(((res_reg1_2_3_urb - res_reg2_5d_fit_urb)/res_reg1_2_3_urb)^2))
  
  se_cc0_5d_urb <- sqrt(sum((res_reg1_2_3_urb - res_reg2_5d_fit_urb)^2)/(length(bhc_ts_TIA_vec_station[[k]]) - 2)*(1/length(bhc_ts_TIA_vec_station[[k]])+mean(log(bhc_ts_TIA_vec_station[[k]]))^2/sum((log(bhc_ts_TIA_vec_station[[k]])-mean(log(bhc_ts_TIA_vec_station[[k]])))^2))) 
  se_cc1_5d_urb <- sqrt(sum((res_reg1_2_3_urb - res_reg2_5d_fit_urb)^2)/((length(bhc_ts_TIA_vec_station[[k]]) - 2)*sum((log(bhc_ts_TIA_vec_station[[k]])-mean(log(bhc_ts_TIA_vec_station[[k]])))^2)))
  t_cc0_5d_urb <- cc0_5d_urb[k]/se_cc0_5d_urb
  t_cc1_5d_urb <- cc1_5d_urb[k]/se_cc1_5d_urb
  p_cc0_5d_urb[k] = 2*(pt(-abs(t_cc0_5d_urb), df=length(bhc_ts_TIA_vec_station[[k]]) - 1))
  p_cc1_5d_urb[k] = 2*(pt(-abs(t_cc1_5d_urb), df=length(bhc_ts_TIA_vec_station[[k]]) - 1)) 
  
  # Test residual adequacy 
  res_reg2_5d_urb <- as.numeric(res_reg1_2_3_urb - res_reg2_5d_fit_urb)
  ppcc_res_reg2_urb <- ppcc.test(res_reg2_5d_urb)
  pval_ppcc_res_reg2_5d_urb[k] = as.numeric(ppcc_res_reg2_urb[2])
  dw_res_reg2_urb<- dwtest(res_reg2_5d_urb ~ bhc_ts_TIA_vec_station[[k]])
  pval_dw_res_reg2_5d_urb[k] = as.numeric(dw_res_reg2_urb$p.value) 
  
  # Test for heteroscedasticity of residuals of second regression
  res_reg2_2_3_5d_urb <- ((res_reg1_2_3_urb - res_reg2_5d_fit_urb)^2)^(1/3)
  res_reg2_2_3_5d_urb.lm <- lm(res_reg2_2_3_5d_urb ~ bhc_ts_TIA_vec_station[[k]])
  p_dd1_5d_urb[k] = summary(res_reg2_2_3_5d_urb.lm)$coefficients[2,4]
  
  '
  # Model 5E (Method of moments, log-transformed squared residuals)
  log_res_reg1_2_time <- log(res_reg1_time[[k]]^2)
  rho_t_res_reg2_time <- cor(wy_order_time[[k]],log_res_reg1_2_time)
  sd_log_res_reg1_2_time <- sd(log_res_reg1_2_time)
  cc0_5e_time[k] = mean(log_res_reg1_2_time) - (rho_t_res_reg2_time * sd_log_res_reg1_2_time  * mean(wy_order_time[[k]]) / sd(wy_order_time[[k]]))
  cc1_5e_time[k] = rho_t_res_reg2_time * sd_log_res_reg1_2_time / sd(wy_order_time[[k]])
  res_reg2_5e_fit_time <- cc0_5e_time[k] + cc1_5e_time[k] * wy_order_time[[k]]
  res_reg2_5e_var_time <- var(log_res_reg1_2_time - res_reg2_5e_fit_time)
  
  cond_var_5e_time[[k]] <- exp((cc0_5e_time[k] + cc1_5e_time[k] * wy_order_time[[k]] + 0.5*res_reg2_5e_var_time))
  cond_var_at_n_5e_time[k] = exp(res_dof_corr_time[k]*(cc0_5e_time[k] + cc1_5e_time[k] * max(wy_order_time[[k]])) + 0.5*res_reg2_5e_var_time)
  #plot(wy_order,log_res_reg1_2,col="blue")
  #lines(wy_order,log_res_reg2_5e) 
  Q_99_at_n_5e_time[k] = exp(cond_med_at_n_time[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5e_time[k])) 
  rsquared_5e_time[k] = 1-sum((log_res_reg1_2_time - res_reg2_5e_fit_time)^2)/sum((log_res_reg1_2_time - mean(log_res_reg1_2_time))^2) 
  adj_rsquared_5e_time[k] = 1-adj_r2_corr_time[k]*sum((log_res_reg1_2_time - res_reg2_5e_fit_time)^2)/sum((log_res_reg1_2_time - mean(log_res_reg1_2_time))^2) 
  rmse_5e_time[k] = sqrt(1/length(wy_order_time[[k]])*sum((res_reg1_2_3_time^3 - cond_var_5e_time[[k]])^2)) 
  rrmse_5e_time[k] = sqrt(1/length(wy_order_time[[k]])*sum(((res_reg1_2_3_time^3 - cond_var_5e_time[[k]])/cond_var_5e_time[[k]])^2))
  
  se_cc0_5e_time <- sqrt(sum((log_res_reg1_2_time - res_reg2_5e_fit_time)^2)/(length(wy_order_time[[k]]) - 2)*(1/length(wy_order_time[[k]])+mean(wy_order_time[[k]])^2/sum((wy_order[[k]]-mean(wy_order_time[[k]]))^2))) 
  se_cc1_5e_time <- sqrt(sum((log_res_reg1_2_time - res_reg2_5e_fit_time)^2)/((length(wy_order_time[[k]]) - 2)*sum((wy_order_time[[k]]-mean(wy_order_time[[k]]))^2))) 
  t_cc0_5e_time <- cc0_5e_time[k]/se_cc0_5e_time
  t_cc1_5e_time <- cc1_5e_time[k]/se_cc1_5e_time
  p_cc0_5e_time[k] = 2*(pt(-abs(t_cc0_5e_time), df=length(wy_order_time[[k]]) - 1)) 
  p_cc1_5e_time[k] = 2*(pt(-abs(t_cc1_5e_time), df=length(wy_order_time[[k]]) - 1)) 
  
  # Test residual adequacy 
  res_reg2_5e_time <- as.numeric(log_res_reg1_2_time - res_reg2_5e_fit_time) 
  ppcc_res_reg2_time <- ppcc.test(res_reg2_5e_time) 
  pval_ppcc_res_reg2_5e_time[k] = as.numeric(ppcc_res_reg2_time[2]) 
  dw_res_reg2_time <- dwtest(res_reg2_5e_time ~ wy_order_time[[k]]) 
  pval_dw_res_reg2_5e_time[k] = as.numeric(dw_res_reg2_time$p.value) 
  
  # Test for heteroscedasticity of residuals of second regression
  log_res_reg2_5e_time <- ((log_res_reg1_2_time - res_reg2_5e_fit_time)^2)^(1/3)
  log_res_reg2_5e_time.lm <- lm(log_res_reg2_5e_time ~ wy_order_time[[k]])
  p_dd1_5e_time[k] = summary(log_res_reg2_5e_time.lm)$coefficients[2,4]
  
  # Model 5f: Direct estimation of standard deviation
  res_reg1_1_3_time <- (res_reg1_time[[k]]^2)^(1/6)
  rho_t_res_reg2_time <- cor(wy_order_time[[k]],res_reg1_1_3_time)
  sd_res_reg1_1_3_time  <- sd(res_reg1_1_3_time)
  cc0_5f_time[k] = mean(res_reg1_1_3_time) - (rho_t_res_reg2_time * sd_res_reg1_1_3_time  * mean(wy_order_time[[k]]) / sd(wy_order_time[[k]]))
  cc1_5f_time[k] = rho_t_res_reg2_time * sd_res_reg1_1_3_time / sd(wy_order_time[[k]])
  res_reg2_5f_fit_time <- cc0_5f_time[k] + cc1_5f_time[k] * wy_order_time[[k]]
  res_reg2_5f_var_time <- var(res_reg1_1_3_time-res_reg2_5f_fit_time)
  
  cond_var_5f_time[[k]] = (cc0_5f_time[k] + cc1_5f_time[k] * wy_order_time[[k]])^6  
  cond_var_at_n_5f_time[k] = res_dof_corr_time[k]*(cc0_5f_time[k] + cc1_5f_time[k] * max(wy_order_time[[k]]))^6  
  
  Q_99_at_n_5f_time[k] = exp(cond_med_at_n_time[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5f_time[k])) 
  rsquared_5f_time[k] = 1-sum((res_reg1_1_3_time - res_reg2_5f_fit_time)^2)/sum((res_reg1_1_3_time - mean(res_reg1_1_3_time))^2) 
  adj_rsquared_5f_time[k] = 1-adj_r2_corr_time[k]*sum((res_reg1_1_3_time - res_reg2_5f_fit_time)^2)/sum((res_reg1_1_3_time - mean(res_reg1_1_3_time))^2) 
  rmse_5f_time[k] = sqrt(1/length(wy_order_time[[k]])*sum((res_reg1_1_3_time - cond_var_5f_time[[k]])^2))
  rrmse_5f_time[k] = sqrt(1/length(wy_order_time[[k]])*sum(((res_reg1_1_3_time - cond_var_5f_time[[k]])/cond_var_5f_time[[k]])^2))
  
  se_cc0_5f_time <- sqrt(sum((res_reg1_1_3_time - res_reg2_5f_fit_time)^2)/(length(wy_order_time[[k]]) - 2)*(1/length(wy_order_time[[k]])+mean(wy_order_time[[k]])^2/sum((wy_order_time[[k]]-mean(wy_order_time[[k]]))^2))) 
  se_cc1_5f_time <- sqrt(sum((res_reg1_1_3_time - res_reg2_5f_fit_time)^2)/((length(wy_order_time[[k]]) - 2)*sum((wy_order_time[[k]]-mean(wy_order_time[[k]]))^2))) 
  t_cc0_5f_time <- cc0_5f_time[k]/se_cc0_5f_time 
  t_cc1_5f_time <- cc1_5f_time[k]/se_cc1_5f_time
  p_cc0_5f_time[k] = 2*(pt(-abs(t_cc0_5f_time), df=length(wy_order_time[[k]]) - 1)) 
  p_cc1_5f_time[k] = 2*(pt(-abs(t_cc1_5f_time), df=length(wy_order_time[[k]]) - 1)) 
  
  # Test residual adequacy 
  res_reg2_5f_time <- as.numeric(res_reg1_1_3_time - res_reg2_5f_fit_time) 
  ppcc_res_reg2_time <- ppcc.test(res_reg2_5f_time) 
  pval_ppcc_res_reg2_5f_time[k] = as.numeric(ppcc_res_reg2_time[2]) 
  dw_res_reg2_time <- dwtest(res_reg2_5f_time ~ wy_order_time[[k]]) 
  pval_dw_res_reg2_5f_time[k] = as.numeric(dw_res_reg2_time$p.value) 
  
  # Test for heteroscedasticity of residuals of second regression
  res_reg1_1_3_5f_time <- ((res_reg1_1_3_time - res_reg2_5f_fit_time)^2)^(1/3)
  res_reg1_1_3_5f_time.lm <- lm(res_reg1_1_3_5f_time ~ wy_order_time[[k]])
  p_dd1_5f_time[k] = summary(res_reg1_1_3_5f_time.lm)$coefficients[2,4]
  '
  
  # MODEL 6: ITERATIVE RE(WEIGHTED) LEAST SQUARES --------------------------------------------------
  
  # Fit Model 6A 
  
  # Establish output arrays for each record (length-dependent)
  wt_6a_urb <- array(NA,dim=c(length(bhc_ts_TIA_vec_station[[k]]),10))
  res_reg1_iwls_6a_urb <- array(NA,dim=c(length(bhc_ts_TIA_vec_station[[k]]),10))
  res_reg1_iwls_2_3_6a_urb <- array(NA,dim=c(length(bhc_ts_TIA_vec_station[[k]]),10))
  cond_var_6a_urb <- array(NA,dim=c(length(bhc_ts_TIA_vec_station[[k]]),10))
  
  for (i in 1:10){
    if(i == 1){
      wt_6a_urb[,i] = rep(1,length(bhc_ts_TIA_vec_station[[k]])) #Make a list
    } else {
      wt_6a_urb[,i] = 1/(cond_var_6a_urb[,i-1]) #Make a list
    }
    flood_reg1_iwls_urb.lm <- lm(log(peak_flow_urb + 0.01) ~ bhc_ts_TIA_vec_station_2[[k]], weights=wt_6a_urb[,i]) # Make a list
    b0_reg1_iwls_6a_urb[k,i] = as.numeric((flood_reg1_iwls_urb.lm)$coefficients[1])
    b1_reg1_iwls_6a_urb[k,i] = as.numeric((flood_reg1_iwls_urb.lm)$coefficients[2])
    res_reg1_iwls_6a_urb[,i] = as.numeric(residuals(flood_reg1_iwls_urb.lm)) # Make a list
    res_reg1_iwls_2_3_6a_urb[,i] = (as.numeric(residuals(flood_reg1_iwls_urb.lm)^2))^(1/3)# Make a list
    res_reg1_iwls_2_3_6a_urb.lm <- lm(res_reg1_iwls_2_3_6a_urb[,i] ~ bhc_ts_TIA_vec_station[[k]]) # Make a list
    cc0_6a_urb[k,i] = as.numeric((res_reg1_iwls_2_3_6a_urb.lm$coefficients[1]))
    cc1_6a_urb[k,i] = as.numeric((res_reg1_iwls_2_3_6a_urb.lm$coefficients[2]))
    cond_var_6a_urb[,i] = as.numeric(fitted(res_reg1_iwls_2_3_6a_urb.lm)) # Make this a list
    for (l in 1:length(peak_flow_urb))
      if (cond_var_6a_urb[l,i] < 0.001){
        cond_var_6a_urb[l,i] <- 0.001
      } 
    
  } 
  
  wt_6a_list_urb[[k]] = wt_6a_urb[,10]  
  res_reg1_iwls_6a_list_urb[[k]] = res_reg1_iwls_6a_urb[,10]  
  res_reg1_iwls_2_3_6a_list_urb[[k]] = res_reg1_iwls_2_3_6a_urb[,10]  
  cond_var_6a_list_urb[[k]] = cond_var_6a_urb[,10]  
  res_reg2_6a_var_urb <- res_dof_corr_urb[k]*var(res_reg1_iwls_2_3_6a_urb[,10]-cond_var_6a_urb[,10]) 
  
  cond_var_at_n_6a_urb[k] = res_dof_corr_urb[k]*((cc0_6a_urb[k,10] + cc1_6a_urb[k,10] * max(bhc_ts_TIA_vec_station[[k]]))^3 + 3*res_reg2_6a_var_urb*(cc0_6a_urb[k,10] + cc1_6a_urb[k,10] * max(bhc_ts_TIA_vec_station[[k]])) + mean((res_reg1_2_3_urb-cond_var_6a_urb[,10])^3))
  #plot(wy_order,res_reg1_2_3,col="blue")
  #lines(wy_order,res_reg2_5a)
  Q_99_at_n_6a_urb[k] = exp(cond_med_at_n_urb[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_6a_urb[k]))
  rsquared_6a_urb[k] = 1-sum((res_reg1_iwls_2_3_6a_urb[,max(i)] - cond_var_6a_urb[,max(i)])^2)/sum((res_reg1_iwls_2_3_6a_urb[,max(i)] - mean(res_reg1_iwls_2_3_6a_urb[,max(i)]))^2)
  adj_rsquared_6a_urb[k] = 1-adj_r2_corr_urb[k]*sum((res_reg1_iwls_2_3_6a_urb[,max(i)] - cond_var_6a_urb[,max(i)])^2)/sum((res_reg1_iwls_2_3_6a_urb[,max(i)] - mean(res_reg1_iwls_2_3_6a_urb[,max(i)]))^2)
  rmse_6a_urb[k] = sqrt(1/length(bhc_ts_TIA_vec_station[[k]])*sum((res_reg1_iwls_2_3_6a_urb[,max(i)] - cond_var_6a_urb[,max(i)])^2))
  rrmse_6a_urb[k] = sqrt(1/length(bhc_ts_TIA_vec_station[[k]])*sum(((res_reg1_iwls_2_3_6a_urb[,max(i)] - cond_var_6a_urb[,max(i)])/res_reg1_iwls_2_3_6a_urb[,max(i)])^2))
  
  se_cc1_6a_urb <- sqrt(sum((res_reg1_iwls_2_3_6a_urb[,max(i)] - cond_var_6a_urb[,max(i)])^2)/((length(bhc_ts_TIA_vec_station[[k]]) - 2)*sum((bhc_ts_TIA_vec_station[[k]] - mean(bhc_ts_TIA_vec_station[[k]]))^2)))
  t_cc1_6a_urb <- cc1_6a_urb[k]/se_cc1_6a_urb
  p_cc1_6a_urb[k] = 2*(pt(-abs(t_cc1_6a_urb), df=length(bhc_ts_TIA_vec_station[[k]]) - 1)) 
  
  # Fit model 6B (Quadratic with IWLS)
  
  # Establish output arrays for each record (length-dependent)
  wt_6b_urb <- array(NA,dim=c(length(bhc_ts_TIA_vec_station_2[[k]]),10))
  res_reg1_iwls_6b_urb <- array(NA,dim=c(length(bhc_ts_TIA_vec_station_2[[k]]),10))
  res_reg1_iwls_2_3_6b_urb <- array(NA,dim=c(length(bhc_ts_TIA_vec_station_2[[k]]),10))
  cond_var_6b_urb <- array(NA,dim=c(length(bhc_ts_TIA_vec_station_2[[k]]),10))
  
  for (i in 1:10){
    if(i == 1){
      wt_6b_urb[,i] = rep(1,length(bhc_ts_TIA_vec_station_2[[k]])) #Make a list
    } else {
      wt_6b_urb[,i] = 1/(cond_var_6b_urb[,i-1]) #Make a list
    }
    flood_reg1_iwls_urb.lm <- lm(log(peak_flow_urb + 0.01) ~ bhc_ts_TIA_vec_station_2[[k]], weights=wt_6b_urb[,i]) # Make a list
    b0_reg1_iwls_6b_urb[k,i] = as.numeric((flood_reg1_iwls_urb.lm)$coefficients[1])
    b1_reg1_iwls_6b_urb[k,i] = as.numeric((flood_reg1_iwls_urb.lm)$coefficients[2])
    res_reg1_iwls_6b_urb[,i] = as.numeric(residuals(flood_reg1_iwls_urb.lm)) # Make a list
    res_reg1_iwls_2_3_6b_urb[,i] = (as.numeric(residuals(flood_reg1_iwls_urb.lm)^2))^(1/3)# Make a list
    res_reg1_iwls_2_3_6b_urb.lm <- lm(res_reg1_iwls_2_3_6b_urb[,i] ~ bhc_ts_TIA_vec_station_2[[k]]) # Make a list
    cc0_6b_urb[k,i] = as.numeric((res_reg1_iwls_2_3_6b_urb.lm$coefficients[1]))
    cc1_6b_urb[k,i] = as.numeric((res_reg1_iwls_2_3_6b_urb.lm$coefficients[2]))
    cond_var_6b_urb[,i] = as.numeric(fitted(res_reg1_iwls_2_3_6b_urb.lm)) # Make this a list
    for (l in 1:length(peak_flow_urb))
      if (cond_var_6b_urb[l,i] < 0.001){
        cond_var_6b_urb[l,i] <- 0.001
      } 
    
  } 
  
  wt_6b_list_urb[[k]] = wt_6b_urb[,10]
  res_reg1_iwls_6b_list_urb[[k]] = res_reg1_iwls_6b_urb[,10]
  res_reg1_iwls_2_3_6b_list_urb[[k]] = res_reg1_iwls_2_3_6b_urb[,10]
  cond_var_6b_list_urb[[k]] = cond_var_6b_urb[,10]
  res_reg2_6b_var_urb <- res_dof_corr_urb[k]*var(res_reg1_iwls_2_3_6b_urb[,10]-cond_var_6b_urb[,10])
  
  cond_var_at_n_6b_urb[k] = res_dof_corr_urb[k]*((cc0_6b_urb[k,10] + cc1_6b_urb[k,10] * max(bhc_ts_TIA_vec_station_2[[k]]))^3 + 3*res_reg2_6b_var_urb*(cc0_6b_urb[k,10] + cc1_6b_urb[k,10] * max(bhc_ts_TIA_vec_station_2[[k]])) + mean((res_reg1_iwls_2_3_6b_urb[,10]-cond_var_6b_urb[,10])^3))
  
  #plot(wy_order,res_reg1_2_3,col="blue")
  #lines(wy_order,res_reg2_5a)
  Q_99_at_n_6b_urb[k] = exp(cond_med_at_n_urb[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_6b_urb[k]))
  rsquared_6b_urb[k] = 1-sum((res_reg1_iwls_2_3_6b_urb[,max(i)] - cond_var_6b_urb[,max(i)])^2)/sum((res_reg1_iwls_2_3_6b_urb[,max(i)] - mean(res_reg1_iwls_2_3_6b_urb[,max(i)]))^2)
  adj_rsquared_6b_urb[k] = 1-adj_r2_corr_urb[k]*sum((res_reg1_iwls_2_3_6b_urb[,max(i)] - cond_var_6b_urb[,max(i)])^2)/sum((res_reg1_iwls_2_3_6b_urb[,max(i)] - mean(res_reg1_iwls_2_3_6b_urb[,max(i)]))^2)
  rmse_6b_urb[k] = sqrt(1/length(bhc_ts_TIA_vec_station_2[[k]])*sum((res_reg1_iwls_2_3_6b_urb[,max(i)] - cond_var_6b_urb[,max(i)])^2))
  rrmse_6b_urb[k] = sqrt(1/length(bhc_ts_TIA_vec_station_2[[k]])*sum(((res_reg1_iwls_2_3_6b_urb[,max(i)] - cond_var_6b_urb[,max(i)])/res_reg1_iwls_2_3_6b_urb[,max(i)])^2))
  
  se_cc1_6b_urb <- sqrt(sum((res_reg1_iwls_2_3_6b_urb[,max(i)] - cond_var_6b_urb[,max(i)])^2)/((length(bhc_ts_TIA_vec_station_2[[k]]) - 2)*sum((bhc_ts_TIA_vec_station_2[[k]]^2-mean(bhc_ts_TIA_vec_station_2[[k]]^2))^2)))
  t_cc1_6b_urb <- cc1_6b_urb[k]/se_cc1_6b_urb
  p_cc1_6b_urb[k] = 2*(pt(-abs(t_cc1_6b_urb), df=length(bhc_ts_TIA_vec_station_2[[k]]) - 1)) 
  
  # Fit model 6c (Exponential with IWLS)
  
  '
  # Establish output arrays for each record (length-dependent)
  wt_6c_urb <- array(NA,dim=c(length(wy_order_urb[[k]]),10))
  res_reg1_iwls_6c_urb <- array(NA,dim=c(length(wy_order_urb[[k]]),10))
  res_reg1_iwls_2_3_6c_urb <- array(NA,dim=c(length(wy_order_urb[[k]]),10))
  cond_var_6c_urb <- array(NA,dim=c(length(wy_order_urb[[k]]),10))
  
  for (i in 1:10){
    if(i == 1){
      wt_6c_urb[,i] = rep(1,length(wy_order_urb[[k]])) #Make a list
    } else {
      wt_6c_urb[,i] = 1/(cond_var_6c_urb[,i-1]) #Make a list
    }
    flood_reg1_iwls_urb.lm <- lm(log(peak_flow_urb + 0.01) ~ wy_order_urb[[k]], weights=wt_6c_urb[,i]) # Make a list
    b0_reg1_iwls_6c_urb[k,i] = as.numeric((flood_reg1_iwls_urb.lm)$coefficients[1])
    b1_reg1_iwls_6c_urb[k,i] = as.numeric((flood_reg1_iwls_urb.lm)$coefficients[2])
    res_reg1_iwls_6c_urb[,i] = as.numeric(residuals(flood_reg1_iwls_urb.lm)) # Make a list
    res_reg1_iwls_2_3_6c_urb[,i] = (as.numeric(residuals(flood_reg1_iwls_urb.lm)^2))^(1/3)# Make a list
    res_reg1_iwls_2_3_6c_urb.lm <- lm(res_reg1_iwls_2_3_6c_urb[,i] ~ exp(wy_order_urb[[k]])) # Make a list
    cc0_6c_urb[k,i] = as.numeric((res_reg1_iwls_2_3_6c_urb.lm$coefficients[1])) 
    cc1_6c_urb[k,i] = as.numeric((res_reg1_iwls_2_3_6c_urb.lm$coefficients[2])) 
    cond_var_6c_urb[,i] = as.numeric(fitted(res_reg1_iwls_2_3_6c_urb.lm))  # Make this a list
    for (l in 1:length(peak_flow_urb))
      if (cond_var_6c_urb[l,i] < 0.001){
        cond_var_6c_urb[l,i] <- 0.001
      } 
    
  } 
  
  wt_6c_list_urb[[k]] = wt_6c_urb[,10]
  res_reg1_iwls_6c_list_urb[[k]] = res_reg1_iwls_6c_urb[,10]
  res_reg1_iwls_2_3_6c_list_urb[[k]] = res_reg1_iwls_2_3_6c_urb[,10]
  cond_var_6c_list_urb[[k]] = cond_var_6c_urb[,10]
  cond_var_at_n_6c_urb[k] = (exp(cc0_6c_urb[k,10] + cc1_6c_urb[k,10] * max(wy_order_urb[[k]])))^3
  #cond_var_at_n_6c[k] = (cc0_6c[k,10] + cc1_6c[k,10] * max(exp(wy_order[[k]])))^3 + 
  #  var(res_reg1_iwls_2_3_6c[,10]-cond_var_6c[,10])*(3*cc0_6c[k,10] + 3*cc1_6c[k,10] * max(exp(wy_order[[k]]))) + mean((res_reg1_iwls_2_3_6c[,10]-cond_var_6c[,10])^3)
  #plot(wy_order[[k]],res_reg1_2_3,col="blue")
  #lines(wy_order[[k]],res_reg2_5a)
  Q_99_at_n_6c_urb[k] = exp(cond_med_at_n_urb[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_6c_urb[k]))
  rsquared_6c_urb[k] = 1-sum((res_reg1_iwls_2_3_6c_urb[,max(i)] - cond_var_6c_urb[,max(i)])^2)/sum((res_reg1_iwls_2_3_6c_urb[,max(i)] - mean(res_reg1_iwls_2_3_6c_urb[,max(i)]))^2)
  adj_rsquared_6c_urb[k] = 1-adj_r2_corr_urb[k]*sum((res_reg1_iwls_2_3_6c_urb[,max(i)] - cond_var_6c_urb[,max(i)])^2)/sum((res_reg1_iwls_2_3_6c_urb[,max(i)] - mean(res_reg1_iwls_2_3_6c_urb[,max(i)]))^2)
  rmse_6c_urb[k] = sqrt(1/length(wy_order_urb[[k]])*sum((res_reg1_iwls_2_3_6c_urb[,max(i)] - cond_var_6c_urb[,max(i)])^2))
  rrmse_6c_urb[k] = sqrt(1/length(wy_order_urb[[k]])*sum(((res_reg1_iwls_2_3_6c_urb[,max(i)] - cond_var_6c_urb[,max(i)])/res_reg1_iwls_2_3_6a_urb[,max(i)])^2))
  
  se_cc1_6c_urb <- sqrt(sum((res_reg1_iwls_2_3_6c_urb[,max(i)] - cond_var_6c_urb[,max(i)])^2)/((length(wy_order_urb[[k]]) - 2)*sum((exp(wy_order_urb[[k]])-mean(exp(wy_order_urb[[k]])))^2)))
  t_cc1_6c_urb <- cc1_6c_urb[k]/se_cc1_6c_urb
  p_cc1_6c_urb[k] = 2*(pt(-abs(t_cc1_6c_urb), df=length(wy_order_urb[[k]]) - 1)) 
  '
  
  # Fit model 6d (Logarithmic with IWLS)
  
  # Establish output arrays for each record (length-dependent)
  wt_6d_urb <- array(NA,dim=c(length(bhc_ts_TIA_vec_station[[k]]),10))
  res_reg1_iwls_6d_urb <- array(NA,dim=c(length(bhc_ts_TIA_vec_station[[k]]),10))
  res_reg1_iwls_2_3_6d_urb <- array(NA,dim=c(length(bhc_ts_TIA_vec_station[[k]]),10))
  cond_var_6d_urb <- array(NA,dim=c(length(bhc_ts_TIA_vec_station[[k]]),10))
  
  for (i in 1:10){
    if(i == 1){
      wt_6d_urb[,i] = rep(1,length(bhc_ts_TIA_vec_station[[k]])) #Make a list
    } else {
      wt_6d_urb[,i] = 1/(cond_var_6d_urb[,i-1]) #Make a list
    }
    flood_reg1_iwls_urb.lm <- lm(log(peak_flow_urb + 0.01) ~ bhc_ts_TIA_vec_station[[k]], weights=wt_6d_urb[,i]) # Make a list
    b0_reg1_iwls_6d_urb[k,i] = as.numeric((flood_reg1_iwls_urb.lm)$coefficients[1])
    b1_reg1_iwls_6d_urb[k,i] = as.numeric((flood_reg1_iwls_urb.lm)$coefficients[2])
    res_reg1_iwls_6d_urb[,i] = as.numeric(residuals(flood_reg1_iwls_urb.lm)) # Make a list
    res_reg1_iwls_2_3_6d_urb[,i] = (as.numeric(residuals(flood_reg1_iwls_urb.lm)^2))^(1/3)# Make a list
    res_reg1_iwls_2_3_6d_urb.lm <- lm(res_reg1_iwls_2_3_6d_urb[,i] ~ log(bhc_ts_TIA_vec_station[[k]])) # Make a list
    cc0_6d_urb[k,i] = as.numeric((res_reg1_iwls_2_3_6d_urb.lm$coefficients[1])) 
    cc1_6d_urb[k,i] = as.numeric((res_reg1_iwls_2_3_6d_urb.lm$coefficients[2])) 
    cond_var_6d_urb[,i] = as.numeric(fitted(res_reg1_iwls_2_3_6d_urb.lm))  # Make this a list
    for (l in 1:length(peak_flow_urb))
      if (cond_var_6d_urb[l,i] < 0.001){
        cond_var_6d_urb[l,i] <- 0.001
      } 
    
  } 
  
  wt_6d_list_urb[[k]] = wt_6d_urb[,10]
  res_reg1_iwls_6d_list_urb[[k]] = res_reg1_iwls_6d_urb[,10]
  res_reg1_iwls_2_3_6d_list_urb[[k]] = res_reg1_iwls_2_3_6d_urb[,10]
  cond_var_6d_list_urb[[k]] = cond_var_6d_urb[,10]
  res_reg2_6d_var_urb <- res_dof_corr_urb[k]*var(res_reg1_iwls_2_3_6d_urb[,10]-cond_var_6d_urb[,10])
  
  cond_var_at_n_6d_urb[k] = res_dof_corr_urb[k]*((cc0_6d_urb[k,10] + cc1_6d_urb[k,10] * max(log(bhc_ts_TIA_vec_station[[k]])))^3 + res_reg2_6d_var_urb*(3*cc0_6d_urb[k,10] + 3*cc1_6d_urb[k,10] * max(log(bhc_ts_TIA_vec_station[[k]]))) + mean((res_reg1_iwls_2_3_6d_urb[,10]-cond_var_6d_urb[,10])^3))
  #plot(wy_order,res_reg1_2_3,col="blue") 
  #lines(wy_order,res_reg2_5a)
  Q_99_at_n_6d_urb[k] = exp(cond_med_at_n_urb[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_6d_urb[k]))
  rsquared_6d_urb[k] = 1-sum((res_reg1_iwls_2_3_6d_urb[,max(i)] - cond_var_6d_urb[,max(i)])^2)/sum((res_reg1_iwls_2_3_6d_urb[,max(i)] - mean(res_reg1_iwls_2_3_6d_urb[,max(i)]))^2)
  adj_rsquared_6d_urb[k] = 1-adj_r2_corr_urb[k]*sum((res_reg1_iwls_2_3_6d_urb[,max(i)] - cond_var_6d_urb[,max(i)])^2)/sum((res_reg1_iwls_2_3_6d_urb[,max(i)] - mean(res_reg1_iwls_2_3_6d_urb[,max(i)]))^2)
  rmse_6d_urb[k] = sqrt(1/length(bhc_ts_TIA_vec_station[[k]])*sum((res_reg1_iwls_2_3_6d_urb[,max(i)] - cond_var_6d_urb[,max(i)])^2))
  rrmse_6d_urb[k] = sqrt(1/length(bhc_ts_TIA_vec_station[[k]])*sum(((res_reg1_iwls_2_3_6d_urb[,max(i)] - cond_var_6d_urb[,max(i)])/res_reg1_iwls_2_3_6d_urb[,max(i)])^2))
  
  se_cc1_6d_urb <- sqrt(sum((res_reg1_iwls_2_3_6d_urb[,max(i)] - cond_var_6d_urb[,max(i)])^2)/((length(bhc_ts_TIA_vec_station[[k]]) - 2)*sum((log(bhc_ts_TIA_vec_station[[k]])-mean(log(bhc_ts_TIA_vec_station[[k]])))^2)))
  t_cc1_6d_urb <- cc1_6d_urb[k]/se_cc1_6d_urb 
  p_cc1_6d_urb[k] = 2*(pt(-abs(t_cc1_6d_urb), df=length(bhc_ts_TIA_vec_station[[k]]) - 1))  
  
  # Fit Model A using GLM 
  glm.cond_var_a_urb <- glm2(formula = res_reg1_2_3_urb ~ bhc_ts_TIA_vec_station_1_3[[k]], family=gaussian) 
  b0_glm_cond_var_a_urb[k] = as.numeric(glm.cond_var_a_urb$coefficients[1]) 
  b1_glm_cond_var_a_urb[k] = as.numeric(glm.cond_var_a_urb$coefficients[2]) 
  res_dev_urb <- as.numeric(glm.cond_var_a_urb$deviance) 
  null_dev_urb <- as.numeric(glm.cond_var_a_urb$null.deviance) 
  pseudo_rsquared_glm_a_urb[k] = 1 - res_dev_urb/null_dev_urb
  cond_var_at_n_glm_a_urb[k] = (b0_glm_cond_var_a_urb[k] + b1_glm_cond_var_a_urb[k] * max(bhc_ts_TIA_vec_station_1_3[[k]]))^3
  Q_99_at_n_glm_a_urb[k] = exp(cond_med_at_n_urb[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_glm_a_urb[k]))
  
  # Fit model using GLM with gamma regression without IWLS
  res_reg1_2_urb <- as.numeric((res_reg1_urb[[k]])^2)
  glm.cond_var_gamma_urb <- glm2(formula = res_reg1_2_urb ~ bhc_ts_TIA_vec_station[[k]], family=Gamma(link=log))
  #Initial conditions: start=c(mean(log(res_reg1_2)),0)
  b0_glm_cond_var_gamma_urb[k] = as.numeric(glm.cond_var_gamma_urb$coefficients[1])
  b1_glm_cond_var_gamma_urb[k] = as.numeric(glm.cond_var_gamma_urb$coefficients[2])
  res_dev_urb <- as.numeric(glm.cond_var_gamma_urb$deviance)
  null_dev_urb <- as.numeric(glm.cond_var_gamma_urb$null.deviance) # intercept-only model
  pseudo_rsquared_glm_gamma_urb[k] = 1 - res_dev_urb/null_dev_urb
  cond_var_at_n_glm_gamma_urb[k] = exp(b0_glm_cond_var_gamma_urb[k] + b1_glm_cond_var_gamma_urb[k] * max(bhc_ts_TIA_vec_station[[k]])) #Transformation bias?
  Q_99_at_n_glm_gamma_urb[k] = exp(cond_med_at_n_urb[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_glm_gamma_urb[k]))
  
  # Fit model using GLM with gamma regression with IWLS
  wt_iwls_glm_gamma_urb <- array(NA,dim=c(length(bhc_ts_TIA_vec_station[[k]]),10))
  res_reg1_iwls_glm_gamma_urb <- array(NA,dim=c(length(bhc_ts_TIA_vec_station[[k]]),10))
  res_reg1_2_iwls_glm_gamma_urb <- array(NA,dim=c(length(bhc_ts_TIA_vec_station[[k]]),10))
  cvar_iwls_glm_gamma_urb <- array(NA,dim=c(length(bhc_ts_TIA_vec_station[[k]]),10))
  
  for (i in 1:10){
    if(i == 1){
      wt_iwls_glm_gamma_urb[,i] = rep(1,length(bhc_ts_TIA_vec_station[[k]])) #Make a list
    } else {
      wt_iwls_glm_gamma_urb[,i] = 1/(cvar_iwls_glm_gamma_urb[,i-1]) #Make a list
    }
    flood_reg1_iwls_urb.lm <- lm(log(peak_flow_urb + 0.01) ~ bhc_ts_TIA_vec_station[[k]], weights=wt_iwls_glm_gamma_urb[,i]) # Make a list
    b0_reg1_iwls_glm_gamma_urb[k,i] = as.numeric((flood_reg1_iwls_urb.lm)$coefficients[1])
    b1_reg1_iwls_glm_gamma_urb[k,i] = as.numeric((flood_reg1_iwls_urb.lm)$coefficients[2])
    res_reg1_iwls_glm_gamma_urb[,i] = as.numeric(residuals(flood_reg1_iwls_urb.lm)) # Make a list
    res_reg1_2_iwls_glm_gamma_urb[,i] = as.numeric((residuals(flood_reg1_iwls_urb.lm))^2)# Make a list
    res_reg1_2_iwls_glm_gamma_urb.lm <- glm2(res_reg1_2_iwls_glm_gamma_urb[,i] ~ bhc_ts_TIA_vec_station[[k]],family=Gamma(link=log)) # Make a list
    cc0_iwls_glm_gamma_urb[k,i] = as.numeric((res_reg1_2_iwls_glm_gamma_urb.lm$coefficients[1]))
    cc1_iwls_glm_gamma_urb[k,i] = as.numeric((res_reg1_2_iwls_glm_gamma_urb.lm$coefficients[2]))
    cvar_iwls_glm_gamma_urb[,i] = as.numeric(fitted(res_reg1_2_iwls_glm_gamma_urb.lm)) # Make this a list
    for (l in 1:length(peak_flow_urb))
      if (cvar_iwls_glm_gamma_urb[l,i] < 0.001){
        cvar_iwls_glm_gamma_urb[l,i] <- 0.001
      } 
    
  } 
  
  wt_iwls_glm_gamma_list_urb[[k]] = wt_iwls_glm_gamma_urb[,10]  
  res_reg1_iwls_glm_gamma_list_urb[[k]] = res_reg1_iwls_glm_gamma_urb[,10]  
  res_reg1_2_iwls_glm_gamma_list_urb[[k]] = res_reg1_2_iwls_glm_gamma_urb[,10]  
  cvar_iwls_glm_gamma_list_urb[[k]] = cvar_iwls_glm_gamma_urb[,10]  
  res_reg2_glm_gamma_var_urb <- res_dof_corr_urb[k]*var(res_reg1_2_iwls_glm_gamma_urb[,10]-cvar_iwls_glm_gamma_urb[,10]) 
  
  cond_var_at_n_iwls_glm_gamma_urb[k] = res_dof_corr_urb[k]*exp((cc0_iwls_glm_gamma_urb[k,10] + cc1_iwls_glm_gamma_urb[k,10] * max(bhc_ts_TIA_vec_station[[k]]))) 
  
  Q_99_at_n_iwls_glm_gamma_urb[k] = exp(cond_med_at_n_urb[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_iwls_glm_gamma_urb[k])) 
  rsquared_iwls_glm_gamma_urb[k] = 1-sum((res_reg1_2_iwls_glm_gamma_urb[,max(i)] - cvar_iwls_glm_gamma_urb[,max(i)])^2)/sum((res_reg1_2_iwls_glm_gamma_urb[,max(i)] - mean(res_reg1_2_iwls_glm_gamma_urb[,max(i)]))^2) 
  adj_rsquared_iwls_glm_gamma_urb[k] = 1 -adj_r2_corr_urb[k]*sum((res_reg1_2_iwls_glm_gamma_urb[,max(i)] - cvar_iwls_glm_gamma_urb[,max(i)])^2)/sum((res_reg1_2_iwls_glm_gamma_urb[,max(i)] - mean(res_reg1_2_iwls_glm_gamma_urb[,max(i)]))^2) 
  rmse_iwls_glm_gamma_urb[k] = sqrt(1/length(bhc_ts_TIA_vec_station[[k]])*sum((res_reg1_2_iwls_glm_gamma_urb[,max(i)] - cvar_iwls_glm_gamma_urb[,max(i)])^2))
  rrmse_iwls_glm_gamma_urb[k] = sqrt(1/length(bhc_ts_TIA_vec_station[[k]])*sum(((res_reg1_2_iwls_glm_gamma_urb[,max(i)] - cvar_iwls_glm_gamma_urb[,max(i)])/cvar_iwls_glm_gamma_urb[,max(i)])^2))
  
  se_cc1_iwls_glm_gamma_urb <- sqrt(sum((res_reg1_2_iwls_glm_gamma_urb[,max(i)] - cvar_iwls_glm_gamma_urb[,max(i)])^2)/((length(bhc_ts_TIA_vec_station[[k]]) - 2)*sum((bhc_ts_TIA_vec_station[[k]]-mean(bhc_ts_TIA_vec_station[[k]]))^2))) 
  t_cc1_iwls_glm_gamma_urb <- cc1_iwls_glm_gamma_urb[k]/se_cc1_iwls_glm_gamma_urb  
  p_cc1_iwls_glm_gamma_urb[k] = 2*(pt(-abs(t_cc1_iwls_glm_gamma_urb), df=length(bhc_ts_TIA_vec_station[[k]]) - 1)) 
  
  # Fit model using GLM with gamma regression with IWLS and quadratic model (7b)
  wt_iwls_glm_gamma_7b_urb  <- array(NA,dim=c(length(bhc_ts_TIA_vec_station[[k]]),10))
  res_reg1_iwls_glm_gamma_7b_urb  <- array(NA,dim=c(length(bhc_ts_TIA_vec_station[[k]]),10))
  res_reg1_2_iwls_glm_gamma_7b_urb  <- array(NA,dim=c(length(bhc_ts_TIA_vec_station[[k]]),10))
  cvar_iwls_glm_gamma_7b_urb <- array(NA,dim=c(length(bhc_ts_TIA_vec_station[[k]]),10))
  
  for (i in 1:10){
    if(i == 1){
      wt_iwls_glm_gamma_7b_urb[,i] = rep(1,length(bhc_ts_TIA_vec_station_2[[k]])) #Make a list
    } else {
      wt_iwls_glm_gamma_7b_urb[,i] = 1/(cvar_iwls_glm_gamma_7b_urb[,i-1]) #Make a list
    }
    flood_reg1_iwls_7b_urb.lm <- lm(log(peak_flow_urb + 0.01) ~ bhc_ts_TIA_vec_station_2[[k]], weights=wt_iwls_glm_gamma_7b_urb[,i]) # Make a list
    b0_reg1_iwls_glm_gamma_7b_urb[k,i] = as.numeric((flood_reg1_iwls_7b_urb.lm)$coefficients[1])
    b1_reg1_iwls_glm_gamma_7b_urb[k,i] = as.numeric((flood_reg1_iwls_7b_urb.lm)$coefficients[2])
    res_reg1_iwls_glm_gamma_7b_urb[,i] = as.numeric(residuals(flood_reg1_iwls_7b_urb.lm)) # Make a list
    res_reg1_2_iwls_glm_gamma_7b_urb[,i] = as.numeric((residuals(flood_reg1_iwls_7b_urb.lm))^2)# Make a list
    res_reg1_2_iwls_glm_gamma_7b_urb.lm <- glm2(res_reg1_2_iwls_glm_gamma_7b_urb[,i] ~ bhc_ts_TIA_vec_station_2[[k]],family=Gamma(link=log)) # Make a list
    cc0_iwls_glm_gamma_7b_urb[k,i] = as.numeric((res_reg1_2_iwls_glm_gamma_7b_urb.lm$coefficients[1]))
    cc1_iwls_glm_gamma_7b_urb[k,i] = as.numeric((res_reg1_2_iwls_glm_gamma_7b_urb.lm$coefficients[2]))
    cvar_iwls_glm_gamma_7b_urb[,i] = as.numeric(fitted(res_reg1_2_iwls_glm_gamma_7b_urb.lm)) # Make this a list
    for (l in 1:length(peak_flow_urb))
      if (cvar_iwls_glm_gamma_7b_urb[l,i] < 0.001){
        cvar_iwls_glm_gamma_7b_urb[l,i] <- 0.001
      } 
    
  } 
  
  wt_iwls_glm_gamma_7b_list_urb[[k]] = wt_iwls_glm_gamma_7b_urb[,10]  
  res_reg1_iwls_glm_gamma_7b_list_urb[[k]] = res_reg1_iwls_glm_gamma_7b_urb[,10]  
  res_reg1_2_iwls_glm_gamma_7b_list_urb[[k]] = res_reg1_2_iwls_glm_gamma_7b_urb[,10]  
  cvar_iwls_glm_gamma_7b_list_urb[[k]] = cvar_iwls_glm_gamma_7b_urb[,10]  
  res_reg2_glm_gamma_7b_var_urb <- res_dof_corr_urb[k]*var(res_reg1_2_iwls_glm_gamma_7b_urb[,10]-cvar_iwls_glm_gamma_7b_urb[,10]) 
  
  cond_var_at_n_iwls_glm_gamma_7b_urb[k] = res_dof_corr_urb[k]*exp((cc0_iwls_glm_gamma_7b_urb[k,10] + cc1_iwls_glm_gamma_7b_urb[k,10] * max(bhc_ts_TIA_vec_station_2[[k]]))) 
  
  Q_99_at_n_iwls_glm_gamma_7b_urb[k] = exp(cond_med_at_n_urb[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_iwls_glm_gamma_7b_urb[k])) 
  rsquared_iwls_glm_gamma_7b_urb[k] = 1-sum((res_reg1_2_iwls_glm_gamma_7b_urb[,max(i)] - cvar_iwls_glm_gamma_7b_urb[,max(i)])^2)/sum((res_reg1_2_iwls_glm_gamma_7b_urb[,max(i)] - mean(res_reg1_2_iwls_glm_gamma_7b_urb[,max(i)]))^2) 
  adj_rsquared_iwls_glm_gamma_7b_urb[k] = 1 -adj_r2_corr_urb[k]*sum((res_reg1_2_iwls_glm_gamma_7b_urb[,max(i)] - cvar_iwls_glm_gamma_7b_urb[,max(i)])^2)/sum((res_reg1_2_iwls_glm_gamma_7b_urb[,max(i)] - mean(res_reg1_2_iwls_glm_gamma_7b_urb[,max(i)]))^2) 
  rmse_iwls_glm_gamma_7b_urb[k] = sqrt(1/length(bhc_ts_TIA_vec_station_2[[k]])*sum((res_reg1_2_iwls_glm_gamma_7b_urb[,max(i)] - cvar_iwls_glm_gamma_7b_urb[,max(i)])^2))
  rrmse_iwls_glm_gamma_7b_urb[k] = sqrt(1/length(bhc_ts_TIA_vec_station_2[[k]])*sum(((res_reg1_2_iwls_glm_gamma_7b_urb[,max(i)] - cvar_iwls_glm_gamma_7b_urb[,max(i)])/cvar_iwls_glm_gamma_7b_urb[,max(i)])^2))
  
  se_cc1_iwls_glm_gamma_7b_urb <- sqrt(sum((res_reg1_2_iwls_glm_gamma_7b_urb[,max(i)] - cvar_iwls_glm_gamma_7b_urb[,max(i)])^2)/((length(bhc_ts_TIA_vec_station_2[[k]]) - 2)*sum((bhc_ts_TIA_vec_station_2[[k]]-mean(bhc_ts_TIA_vec_station_2[[k]]))^2))) 
  t_cc1_iwls_glm_gamma_7b_urb <- cc1_iwls_glm_gamma_7b_urb[k]/se_cc1_iwls_glm_gamma_7b_urb  
  p_cc1_iwls_glm_gamma_7b_urb[k] = 2*(pt(-abs(t_cc1_iwls_glm_gamma_7b_urb), df=length(bhc_ts_TIA_vec_station_2[[k]]) - 1)) 
  
  # Fit model using GLM with gamma regression with IWLS and logarithmic model
  wt_iwls_glm_gamma_7d_urb <- array(NA,dim=c(length(log(bhc_ts_TIA_vec_station[[k]])),10))
  res_reg1_iwls_glm_gamma_7d_urb <- array(NA,dim=c(length(log(bhc_ts_TIA_vec_station[[k]])),10))
  res_reg1_2_iwls_glm_gamma_7d_urb <- array(NA,dim=c(length(log(bhc_ts_TIA_vec_station[[k]])),10))
  cvar_iwls_glm_gamma_7d_urb <- array(NA,dim=c(length(log(bhc_ts_TIA_vec_station[[k]])),10))
  
  for (i in 1:10){
    if(i == 1){
      wt_iwls_glm_gamma_7d_urb[,i] = rep(1,length(log(bhc_ts_TIA_vec_station[[k]]))) #Make a list
    } else {
      wt_iwls_glm_gamma_7d_urb[,i] = 1/(cvar_iwls_glm_gamma_7d_urb[,i-1]) #Make a list
    }
    flood_reg1_iwls_7d_urb.lm <- lm(log(peak_flow_urb + 0.01) ~ bhc_ts_TIA_vec_station[[k]], weights=wt_iwls_glm_gamma_7d_urb[,i]) # Make a list
    b0_reg1_iwls_glm_gamma_7d_urb[k,i] = as.numeric((flood_reg1_iwls_7d_urb.lm)$coefficients[1])
    b1_reg1_iwls_glm_gamma_7d_urb[k,i] = as.numeric((flood_reg1_iwls_7d_urb.lm)$coefficients[2])
    res_reg1_iwls_glm_gamma_7d_urb[,i] = as.numeric(residuals(flood_reg1_iwls_7d_urb.lm)) # Make a list
    res_reg1_2_iwls_glm_gamma_7d_urb[,i] = as.numeric((residuals(flood_reg1_iwls_7d_urb.lm))^2)# Make a list
    res_reg1_2_iwls_glm_gamma_7d_urb.lm <- glm2(res_reg1_2_iwls_glm_gamma_7d_urb[,i] ~ log(bhc_ts_TIA_vec_station[[k]]),family=Gamma(link=log)) # Make a list
    cc0_iwls_glm_gamma_7d_urb[k,i] = as.numeric((res_reg1_2_iwls_glm_gamma_7d_urb.lm$coefficients[1]))
    cc1_iwls_glm_gamma_7d_urb[k,i] = as.numeric((res_reg1_2_iwls_glm_gamma_7d_urb.lm$coefficients[2]))
    cvar_iwls_glm_gamma_7d_urb[,i] = as.numeric(fitted(res_reg1_2_iwls_glm_gamma_7d_urb.lm)) # Make this a list
    for (l in 1:length(peak_flow_urb))
      if (cvar_iwls_glm_gamma_7d_urb[l,i] < 0.001){
        cvar_iwls_glm_gamma_7d_urb[l,i] <- 0.001
      } 
    
    
  } 
  
  wt_iwls_glm_gamma_7d_list_urb[[k]] = wt_iwls_glm_gamma_7d_urb[,10]  
  res_reg1_iwls_glm_gamma_7d_list_urb[[k]] = res_reg1_iwls_glm_gamma_7d_urb[,10]  
  res_reg1_2_iwls_glm_gamma_7d_list_urb[[k]] = res_reg1_2_iwls_glm_gamma_7d_urb[,10]  
  cvar_iwls_glm_gamma_7d_list_urb[[k]] = cvar_iwls_glm_gamma_7d_urb[,10]  
  res_reg2_glm_gamma_7d_var_urb <- res_dof_corr_urb[k]*var(res_reg1_2_iwls_glm_gamma_7d_urb[,10]-cvar_iwls_glm_gamma_7d_urb[,10]) 
  
  cond_var_at_n_iwls_glm_gamma_7d_urb[k] = res_dof_corr_urb[k]*exp((cc0_iwls_glm_gamma_7d_urb[k,10] + cc1_iwls_glm_gamma_7d_urb[k,10] * max(log(bhc_ts_TIA_vec_station[[k]])))) 
  
  Q_99_at_n_iwls_glm_gamma_7d_urb[k] = exp(cond_med_at_n_urb[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_iwls_glm_gamma_7d_urb[k])) 
  rsquared_iwls_glm_gamma_7d_urb[k] = 1-sum((res_reg1_2_iwls_glm_gamma_7d_urb[,max(i)] - cvar_iwls_glm_gamma_7d_urb[,max(i)])^2)/sum((res_reg1_2_iwls_glm_gamma_7d_urb[,max(i)] - mean(res_reg1_2_iwls_glm_gamma_7d_urb[,max(i)]))^2) 
  adj_rsquared_iwls_glm_gamma_7d_urb[k] = 1 -adj_r2_corr_urb[k]*sum((res_reg1_2_iwls_glm_gamma_7d_urb[,max(i)] - cvar_iwls_glm_gamma_7d_urb[,max(i)])^2)/sum((res_reg1_2_iwls_glm_gamma_7d_urb[,max(i)] - mean(res_reg1_2_iwls_glm_gamma_7d_urb[,max(i)]))^2) 
  rmse_iwls_glm_gamma_7d_urb[k] = sqrt(1/length(bhc_ts_TIA_vec_station[[k]])*sum((res_reg1_2_iwls_glm_gamma_7d_urb[,max(i)] - cvar_iwls_glm_gamma_7d_urb[,max(i)])^2)) 
  rrmse_iwls_glm_gamma_7d_urb[k] = sqrt(1/length(bhc_ts_TIA_vec_station[[k]])*sum(((res_reg1_2_iwls_glm_gamma_7d_urb[,max(i)] - cvar_iwls_glm_gamma_7d_urb[,max(i)])/cvar_iwls_glm_gamma_7d_urb[,max(i)])^2))
  
  se_cc1_iwls_glm_gamma_7d_urb <- sqrt(sum((res_reg1_2_iwls_glm_gamma_7d_urb[,max(i)] - cvar_iwls_glm_gamma_7d_urb[,max(i)])^2)/((length(log(bhc_ts_TIA_vec_station[[k]])) - 2)*sum((log(bhc_ts_TIA_vec_station[[k]])-mean(log(bhc_ts_TIA_vec_station[[k]])))^2))) 
  t_cc1_iwls_glm_gamma_7d_urb <- cc1_iwls_glm_gamma_7d_urb[k]/se_cc1_iwls_glm_gamma_7d_urb  
  p_cc1_iwls_glm_gamma_7d_urb[k] = 2*(pt(-abs(t_cc1_iwls_glm_gamma_7d_urb), df=length(log(bhc_ts_TIA_vec_station[[k]])) - 1)) 
  
  '
  # Fit model using GLM with gamma regression with IWLS and "standard deviation" model
  wt_iwls_glm_gamma_7f_time <- array(NA,dim=c(length(wy_order_time[[k]]),10))
  res_reg1_iwls_glm_gamma_7f_time <- array(NA,dim=c(length(wy_order_time[[k]]),10))
  res_reg1_abs_iwls_glm_gamma_7f_time <- array(NA,dim=c(length(wy_order_time[[k]]),10))
  cvar_iwls_glm_gamma_7f_time <- array(NA,dim=c(length(wy_order_time[[k]]),10))
  
  for (i in 1:10){
    if(i == 1){
      wt_iwls_glm_gamma_7f_time[,i] = rep(1,length(wy_order_time[[k]])) #Make a list
    } else {
      wt_iwls_glm_gamma_7f_time[,i] = 1/(cvar_iwls_glm_gamma_7f_time[,i-1]) #Make a list
    }
    flood_reg1_iwls_7f_time.lm <- lm(log(peak_flow_time + 0.01) ~ wy_order_time[[k]], weights=wt_iwls_glm_gamma_7f_time[,i]) # Make a list
    b0_reg1_iwls_glm_gamma_7f_time[k,i] = as.numeric((flood_reg1_iwls_7f_time.lm)$coefficients[1])
    b1_reg1_iwls_glm_gamma_7f_time[k,i] = as.numeric((flood_reg1_iwls_7f_time.lm)$coefficients[2])
    res_reg1_iwls_glm_gamma_7f_time[,i] = as.numeric(residuals(flood_reg1_iwls_7f_time.lm)) # Make a list
    res_reg1_abs_iwls_glm_gamma_7f_time[,i] = abs(as.numeric(residuals(flood_reg1_iwls_7f_time.lm))) # Make a list
    res_reg1_abs_iwls_glm_gamma_7f_time.lm <- glm2(res_reg1_abs_iwls_glm_gamma_7f_time[,i] ~ wy_order_time[[k]],family=Gamma(link=log)) # Make a list
    cc0_iwls_glm_gamma_7f_time[k,i] = as.numeric((res_reg1_abs_iwls_glm_gamma_7f_time.lm$coefficients[1])) 
    cc1_iwls_glm_gamma_7f_time[k,i] = as.numeric((res_reg1_abs_iwls_glm_gamma_7f_time.lm$coefficients[2])) 
    cvar_iwls_glm_gamma_7f_time[,i] = as.numeric(fitted(res_reg1_abs_iwls_glm_gamma_7f_time.lm)^2) # Make this a list
    for (l in 1:length(peak_flow_time))
      if (cvar_iwls_glm_gamma_7f_time[l,i] < 0.001){
        cvar_iwls_glm_gamma_7f_time[l,i] <- 0.001
      } 
    
    
  } 
  
  wt_iwls_glm_gamma_7f_list_time[[k]] = wt_iwls_glm_gamma_7f_time[,10]  
  res_reg1_iwls_glm_gamma_7f_list_time[[k]] = res_reg1_iwls_glm_gamma_7f_time[,10]  
  res_reg1_abs_iwls_glm_gamma_7f_list_time[[k]] = res_reg1_abs_iwls_glm_gamma_7f_time[,10] 
  cvar_iwls_glm_gamma_7f_list_time[[k]] = cvar_iwls_glm_gamma_7f_time[,10]  
  res_reg2_glm_gamma_7f_var_time <- res_dof_corr_time[k]*var(res_reg1_abs_iwls_glm_gamma_7f_time[,10] - sqrt(cvar_iwls_glm_gamma_7f_time[,10]))
  
  cond_var_at_n_iwls_glm_gamma_7f_time[k] = res_dof_corr_time[k]*exp(cc0_iwls_glm_gamma_7f_time[k,10] + cc1_iwls_glm_gamma_7f_time[k,10] * max(wy_order_time[[k]]))^2
  
  Q_99_at_n_iwls_glm_gamma_7f_time[k] = exp(cond_med_at_n_time[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_iwls_glm_gamma_7f_time[k])) 
  rsquared_iwls_glm_gamma_7f_time[k] = 1-sum((res_reg1_abs_iwls_glm_gamma_7f_time[,max(i)] - sqrt(cvar_iwls_glm_gamma_7f_time[,10])))/sum((res_reg1_abs_iwls_glm_gamma_7f_time[,max(i)] - mean(res_reg1_abs_iwls_glm_gamma_7f_time[,max(i)]))^2) 
  adj_rsquared_iwls_glm_gamma_7f_time[k] = 1 -adj_r2_corr_time[k]*sum((res_reg1_abs_iwls_glm_gamma_7f_time[,max(i)] - sqrt(cvar_iwls_glm_gamma_7f_time[,max(i)]))^2)/sum((res_reg1_abs_iwls_glm_gamma_7f_time[,max(i)] - mean(res_reg1_abs_iwls_glm_gamma_7f_time[,max(i)]))^2) 
  rmse_iwls_glm_gamma_7f_time[k] = sqrt(1/length(wy_order_time[[k]])*sum((res_reg1_abs_iwls_glm_gamma_7f_time[,max(i)] - sqrt(cvar_iwls_glm_gamma_7f_time[,max(i)]))^2)) 
  rrmse_iwls_glm_gamma_7f_time[k] = sqrt(1/length(wy_order_time[[k]])*sum(((res_reg1_abs_iwls_glm_gamma_7f_time[,max(i)] - sqrt(cvar_iwls_glm_gamma_7f_time[,max(i)]))/sqrt(cvar_iwls_glm_gamma_7f_time[,max(i)]))^2))
  
  se_cc1_iwls_glm_gamma_7f_time <- sqrt(sum((res_reg1_abs_iwls_glm_gamma_7f_time[,max(i)] - sqrt(cvar_iwls_glm_gamma_7f_time[,max(i)]))^2)/((length(wy_order_time[[k]]) - 2)*sum((wy_order_time[[k]]-mean(wy_order_time[[k]]))^2))) 
  t_cc1_iwls_glm_gamma_7f_time <- cc1_iwls_glm_gamma_7f_time[k]/se_cc1_iwls_glm_gamma_7f_time  
  p_cc1_iwls_glm_gamma_7f_time[k] = 2*(pt(-abs(t_cc1_iwls_glm_gamma_7f_time), df=length(wy_order_time[[k]]) - 1)) 
  '
  
  # MAIN LOOP ENDS
}

# Compare quantiles
# FIX: CODE AS MATRIX

# Compare different quantile estimates
Q99_compar <- c(Q_99_at_n_5a_urb,
                Q_99_at_n_5b_urb,
                Q_99_at_n_5c_urb,
                Q_99_at_n_5d_urb,
                Q_99_at_n_5e_urb,
                Q_99_at_n_6a_urb,
                Q_99_at_n_6b_urb,
                Q_99_at_n_6c_urb,
                Q_99_at_n_6d_urb,
                Q_99_at_n_glm_gamma_urb,
                Q_99_at_n_iwls_glm_gamma_urb)

# COMPARE R^2(adj) of different models for all stations
par(mfrow=c(1,3))

boxplot(adj_rsquared_5a_urb,
        adj_rsquared_iwls_glm_gamma_urb,
        names=c("OLS","GLM-IWLS"),
        ylab="Adjusted R^2",
        ylim=c(-0.1,0.7),
        main="Linear")

boxplot(adj_rsquared_5b_urb,
        adj_rsquared_iwls_glm_gamma_7b_urb,
        names=c("OLS","GLM-IWLS"),
        ylim=c(-0.1,0.7),
        main="Quadratic")

boxplot(adj_rsquared_5d_urb,
        adj_rsquared_iwls_glm_gamma_7d_urb,
        names=c("OLS","GLM-IWLS"),
        ylim=c(-0.1,0.7),
        main="Logarithmic")

par(mfrow=c(1,1))

# COMPARING LINEAR (5A), QUADRATIC (5B) AND LOG-LINEAR (5D) RESIDUAL FUNCTIONS 

# Compute range of estimates for each site for model 5
Q_99_at_n_5_urb <- matrix(data=c(Q_99_at_n_5a_urb,
                                  Q_99_at_n_5b_urb,
                                  Q_99_at_n_5d_urb),
                           nrow=3,ncol=length(Q_99_at_n_5a_urb),byrow=TRUE)
apply(Q_99_at_n_5_urb,2,min)/apply(Q_99_at_n_5_urb,2,max)

# Compare R^2(adj) for each site for model 5
adj_rsquared_5_urb <- matrix(c(adj_rsquared_5a_urb,adj_rsquared_5b_urb,adj_rsquared_5d_urb),nrow=3,ncol=length(Q_99_at_n_5a_urb),byrow=TRUE)
apply(adj_rsquared_5_urb,2,max)-apply(adj_rsquared_5_urb,2,min)

# Which model of model 5 is best? 
order(adj_rsquared_5_urb[,2],decreasing=T)[1]
adj_rsquared_5_best_urb <- apply(adj_rsquared_5_urb,2,function(x)order(x,decreasing=T)[1])

# COmpare OLS and GLM-IWLS R2adj for best models

par(mfrow=c(1,3))

boxplot(adj_rsquared_5a_urb[adj_rsquared_5_best_urb==1],
        adj_rsquared_iwls_glm_gamma_urb[adj_rsquared_5_best_urb==1],
        names=c("OLS","GLM-IWLS"),
        ylab="Adjusted R^2",
        ylim=c(-0.1,0.6),
        main="Linear")

boxplot(adj_rsquared_5b_urb[adj_rsquared_5_best_urb==2],
        adj_rsquared_iwls_glm_gamma_7b_urb[adj_rsquared_5_best_urb==2],
        names=c("OLS","GLM-IWLS"),
        ylim=c(-0.1,0.6),
        main="Quadratic")

boxplot(adj_rsquared_5d_urb[adj_rsquared_5_best_urb==3],
        adj_rsquared_iwls_glm_gamma_7d_urb[adj_rsquared_5_best_urb==3],
        names=c("OLS","GLM-IWLS"),
        ylim=c(-0.1,0.6),
        main="Logarithmic")

par(mfrow=c(1,1))

# Compare linear model quantiles with alternatives
par(mfrow=c(1,2))
plot(Q_99_at_n_5a_urb,Q_99_at_n_5b_urb,log="xy",xlab="100-year flood (Linear)",ylab="100-year flood (Quadratic)",cex.axis=0.85)
lines(seq(0,250000,100),seq(0,250000,100))
plot(Q_99_at_n_5a_urb,Q_99_at_n_5d_urb,log="xy",xlab="100-year flood (Linear)",ylab="100-year flood (Logarithmic)",cex.axis=0.85)
lines(seq(0,250000,100),seq(0,250000,100))
par(mfrow=c(1,1))

# Mean differences
mean((Q_99_at_n_5b_urb - Q_99_at_n_5a_urb)/Q_99_at_n_5a_urb)
mean((Q_99_at_n_5d_urb - Q_99_at_n_5a_urb)/Q_99_at_n_5a_urb)


# Compare OLS-IWLS vs. OLS in terms of bias (%) for each model structure
boxplot((Q_99_at_n_6a_urb-Q_99_at_n_5a_urb)/Q_99_at_n_5a_urb,
        (Q_99_at_n_6b_urb-Q_99_at_n_5b_urb)/Q_99_at_n_5b_urb,
        (Q_99_at_n_6d_urb-Q_99_at_n_5d_urb)/Q_99_at_n_5d_urb,
        names=c("Linear","Quadratic","Logarithmic"),ylab="Bias(%)")
abline(h=0,col="red",lty=2)

# Compare with stationary estimates and trends in the mean only----------------------

# Make comparison for linear model with histogram
par(mfrow=c(1,2),oma = c(0, 0, 3.5, 0))
hist(Q_99_at_n_5a_urb/Q_99_at_n_no_trend_urb,breaks=seq(0,3,0.2),main="Difference with \n stationary estimate",ylim=c(0,20),xlab="Ratio",cex.main=0.85,cex.axis=0.8)
hist(Q_99_at_n_5a_urb/Q_99_at_n_med_only_urb,breaks=seq(0,3,0.2),main="Difference with \n conditional mean estimate",ylim=c(0,20),xlab="Ratio",cex.main=0.85,cex.axis=0.8)
mtext("Current 100-Year Flood Estimates", outer = TRUE, cex = 1.3)
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

# Make comparison for linear model with 1:1 plot
par(mfrow=c(1,2),mar = c(5,6,3,1.5)+0.1)
plot(Q_99_at_n_no_trend_urb,Q_99_at_n_5a_urb,xlim=c(1000,50000),log="xy",xlab="Stationary Q100",ylab="Nonstationary Q100 \n (Trend in mean and Cv)")  
lines(seq(0,400000,100),seq(0,400000,100)) 
plot(Q_99_at_n_med_only_urb,Q_99_at_n_5a_urb,xlim=c(1000,50000),log="xy",xlab="Q100 for trend in mean only",ylab="Nonstationary Q100 \n (Trend in mean and Cv)") 
lines(seq(0,400000,100),seq(0,400000,100))  
mtext("Increasing trend in mean, \n Decreasing trend in Cv",outer = TRUE, cex = 1.3)
par(mfrow=c(1,1),mar = c(5,4,4,2)+0.1)

# Make comparison for best estimates

# Create vector with best quantile estimates
Q_99_at_n_5_best_urb <- c(Q_99_at_n_5a_urb[adj_rsquared_5_best_urb==1],
                           Q_99_at_n_5b_urb[adj_rsquared_5_best_urb==2],
                           Q_99_at_n_5d_urb[adj_rsquared_5_best_urb==3])

# QUICK AND DIRTY SOLUTION TO ALIGN WITH ORDERING OF Q_99_at_n_5_best 
Q_99_at_n_no_trend_best_urb <- c(Q_99_at_n_no_trend_urb[adj_rsquared_5_best_urb==1],
                                  Q_99_at_n_no_trend_urb[adj_rsquared_5_best_urb==2],
                                  Q_99_at_n_no_trend_urb[adj_rsquared_5_best_urb==3])

Q_99_at_n_med_only_best_urb <- c(Q_99_at_n_med_only_urb[adj_rsquared_5_best_urb==1], 
                                  Q_99_at_n_med_only_urb[adj_rsquared_5_best_urb==2],
                                  Q_99_at_n_med_only_urb[adj_rsquared_5_best_urb==3])


# Make comparison for all models
par(mfrow=c(1,2),oma = c(0, 0, 3.5, 0))
hist(Q_99_at_n_5_best_urb/Q_99_at_n_no_trend_best_urb,breaks=seq(1,3,0.2),main="Ratio with \n stationary estimate",ylim=c(0,20),xlab="Ratio",col="gray",cex.main=0.85,cex.axis=0.8)
hist(Q_99_at_n_5_best_urb/Q_99_at_n_med_only_best_urb,breaks=seq(1,3,0.2),main="Ratio with \n conditional mean estimate",ylim=c(0,20),xlab="Ratio",col="gray",cex.main=0.85,cex.axis=0.8)
mtext("Effect of modeling \n increasing trends in variance (N = 36)", outer = TRUE, cex = 1.3)
par(mfrow=c(1,1),oma = c(0, 0, 0, 0))

# Compare OLS and IWLS-GLM

par(mfrow=c(1,3))
#plot(Q_99_at_n_5a,Q_99_at_n_iwls_glm_gamma,xlab="OLS Linear",ylab="GLM-IWLS Linear")
plot(Q_99_at_n_5a_urb,Q_99_at_n_iwls_glm_gamma_urb,
     xlab="OLS",ylab="GLM-IWLS",log="xy",xlim=c(100,1000000),ylim=c(100,1000000))
lines(seq(100,1000000,100),seq(100,1000000,100))
title("Linear model") 

plot(Q_99_at_n_5b_urb,Q_99_at_n_iwls_glm_gamma_7b_urb,
     xlab="OLS",ylab="GLM-IWLS",log="xy",xlim=c(100,1000000),ylim=c(100,1000000))
lines(seq(100,1000000,100),seq(100,1000000,100))
title("Quadratic model") 

plot(Q_99_at_n_5d_urb,Q_99_at_n_iwls_glm_gamma_7d_urb,
     xlab="OLS",ylab="GLM-IWLS",log="xy",xlim=c(100,1000000),ylim=c(100,1000000))
lines(seq(100,1000000,100),seq(100,1000000,100))
title("Logarithmic model") 

par(mfrow=c(1,1))

# Median bias
median((Q_99_at_n_5a_urb - Q_99_at_n_iwls_glm_gamma_urb)/Q_99_at_n_iwls_glm_gamma_urb) 
# 0.07
median((Q_99_at_n_5b_urb - Q_99_at_n_iwls_glm_gamma_7b_urb)/Q_99_at_n_iwls_glm_gamma_7b_urb)
# 0.10
median((Q_99_at_n_5d_urb - Q_99_at_n_iwls_glm_gamma_7d_urb)/Q_99_at_n_iwls_glm_gamma_7d_urb)
# 0.06

# Boxplots of bias: OLS vs. GLM-IWLS for all stations
op <- par(mar=c(5, 6, 4, 2) + 0.1)
boxplot((Q_99_at_n_5a_urb - Q_99_at_n_iwls_glm_gamma_urb)/Q_99_at_n_iwls_glm_gamma_urb*100,
        (Q_99_at_n_5b_urb - Q_99_at_n_iwls_glm_gamma_7b_urb)/Q_99_at_n_iwls_glm_gamma_7b_urb*100,
        (Q_99_at_n_5d_urb - Q_99_at_n_iwls_glm_gamma_7d_urb)/Q_99_at_n_iwls_glm_gamma_7d_urb*100,
        names=c("Linear","Quadratic","Logarithmic"),
        ylab="Percent difference \n (OLS - GLM_IWLS)/GLM_IWLS",cex.lab=0.9,cex.axis=0.8)
abline(h=0,lty=2,col="gray")

IQR((Q_99_at_n_5a_urb - Q_99_at_n_iwls_glm_gamma_urb)/Q_99_at_n_iwls_glm_gamma_urb*100)

# Boxplots of bias: OLS vs. GLM-IWLS for stations where models perform best
op <- par(mar=c(5, 6, 4, 2) + 0.1)
boxplot((Q_99_at_n_5a_urb[adj_rsquared_5_best_urb==1] - Q_99_at_n_iwls_glm_gamma_urb[adj_rsquared_5_best_urb==1])/Q_99_at_n_iwls_glm_gamma_urb[adj_rsquared_5_best_urb==1]*100,
        (Q_99_at_n_5b_urb[adj_rsquared_5_best_urb==2] - Q_99_at_n_iwls_glm_gamma_7b_urb[adj_rsquared_5_best_urb==2])/Q_99_at_n_iwls_glm_gamma_7b_urb[adj_rsquared_5_best_urb==2]*100,
        (Q_99_at_n_5d_urb[adj_rsquared_5_best_urb==3] - Q_99_at_n_iwls_glm_gamma_7d_urb[adj_rsquared_5_best_urb==3])/Q_99_at_n_iwls_glm_gamma_7d_urb[adj_rsquared_5_best_urb==3]*100,
        names=c("Linear \n (N = 14)","Quadratic \n (N = 15)","Logarithmic \n (N = 22)"),
        ylab="Percent difference \n (OLS - GLM_IWLS)/GLM_IWLS",cex.lab=0.9,cex.axis=0.8)
abline(h=0,lty=2,col="gray")

median((Q_99_at_n_5a_urb[adj_rsquared_5_best_urb==1] - Q_99_at_n_iwls_glm_gamma_urb[adj_rsquared_5_best_urb==1])/Q_99_at_n_iwls_glm_gamma_urb[adj_rsquared_5_best_urb==1]*100)
median((Q_99_at_n_5b_urb[adj_rsquared_5_best_urb==2] - Q_99_at_n_iwls_glm_gamma_7b_urb[adj_rsquared_5_best_urb==2])/Q_99_at_n_iwls_glm_gamma_7b_urb[adj_rsquared_5_best_urb==2]*100)
median((Q_99_at_n_5d_urb[adj_rsquared_5_best_urb==3] - Q_99_at_n_iwls_glm_gamma_7d_urb[adj_rsquared_5_best_urb==3])/Q_99_at_n_iwls_glm_gamma_7d_urb[adj_rsquared_5_best_urb==3]*100)


# Boxplots of bias: OLS-IWLS vs. GLM-IWLS
boxplot((Q_99_at_n_6a_urb - Q_99_at_n_iwls_glm_gamma_urb)/Q_99_at_n_iwls_glm_gamma_urb*100,
        (Q_99_at_n_6b_urb - Q_99_at_n_iwls_glm_gamma_7b_urb)/Q_99_at_n_iwls_glm_gamma_7b_urb*100,
        (Q_99_at_n_6d_urb - Q_99_at_n_iwls_glm_gamma_7d_urb)/Q_99_at_n_iwls_glm_gamma_7d_urb*100,
        names=c("Linear","Quadratic","Logarithmic"),
        ylab="Percent bias \n (OLS - GLM-IWLS)/GLM-IWLS")
abline(h=0,lty=2,col="gray")

# FIX: ADD PLOTS FOR MODEL 5A

plot(Q_99_at_n_5b_urb,Q_99_at_n_iwls_glm_gamma_7b_urb)
plot(Q_99_at_n_5b_urb,Q_99_at_n_iwls_glm_gamma_7b_urb,xlim=c(0,20000),ylim=c(0,20000),
     xlab="OLS Estimate",
     ylab="IWLS-GLM Estimate") 
lines(seq(0,20000,1000),seq(0,20000,1000))
title("Quadratic model: OLS vs. IWLS-GLM")

plot(Q_99_at_n_5d_urb,Q_99_at_n_iwls_glm_gamma_7d_urb)
plot(Q_99_at_n_5d_urb,Q_99_at_n_iwls_glm_gamma_7d_urb,xlim=c(0,20000),ylim=c(0,20000))
lines(seq(0,20000,1000),seq(0,20000,1000))
title("Logarithmic model: OLS vs. IWLS-GLM")

# Compare RMSE for OLS and IWLS-GLM
par(mfrow=c(1,3))
boxplot(rrmse_5a_urb,rrmse_iwls_glm_gamma_urb,
        names=c("OLS","GLM-IWLS"),main="Linear",ylab = "Relative RMSE")
boxplot(rrmse_5b_urb,rrmse_iwls_glm_gamma_7b_urb,
        names=c("OLS","GLM-IWLS"),main="Quadratic",ylab = "Relative RMSE")
boxplot(rrmse_5d_urb,rrmse_iwls_glm_gamma_7d_urb,
        names=c("OLS","GLM-IWLS"),main="Logarithmic",ylab = "Relative RMSE")
par(mfrow=c(1,1)) 

# Compare RMSE for OLS and IWLS-OLS
par(mfrow=c(1,3))
boxplot(rrmse_5a_urb,rrmse_6a_urb,names=c("OLS","IWLS"),main="Linear")
boxplot(rrmse_5b_urb,rrmse_6b_urb,names=c("OLS","IWLS"),main="Quadratic")
boxplot(rrmse_5d_urb,rrmse_6d_urb,names=c("OLS","IWLS"),main="Logarithmic")
par(mfrow=c(1,1))

# Compare quantiles for OLS and IWLS-OLS
par(mfrow=c(1,3))
boxplot((Q_99_at_n_6a_urb - Q_99_at_n_5a_urb)/Q_99_at_n_5a_urb*100,
        ylim=c(-20,5),
        ylab="Percent Difference (OLS-IWLS vs. OLS)",
        main="Linear")
boxplot((Q_99_at_n_6b_urb - Q_99_at_n_5b_urb)/Q_99_at_n_5b_urb*100,
        ylim=c(-20,5),
        main="Quadratic")
boxplot((Q_99_at_n_6d_urb - Q_99_at_n_5d_urb)/Q_99_at_n_5d_urb*100,
        ylim=c(-20,5),
        main="Logarithmic")
par(mfrow=c(1,1))


# Explore stations without LN2 distribution
#station_group_urb_no_ln2 <- station_group_urb[which(station_group_urb$adq_ppcc_ln2==0),]

# Quantile estimate comparison
par(mfrow=c(1,2))
plot(Q_99_at_n_med_only_urb,Q_99_at_n_5a_no_trans_adj_urb,log="xy",xlim=c(100,100000),ylim=c(100,100000),
     xlab="Trend in mean only",ylab="Trend in mean and Cv",main="No transf bias \n adjustment",cex.main=0.9)
lines(seq(100,100000,100),seq(100,100000,100))

plot(Q_99_at_n_med_only_urb,Q_99_at_n_5a_urb,log="xy",xlim=c(100,100000),ylim=c(100,100000),
     xlab="Trend in mean only",ylab="Trend in mean and Cv",main="With transf bias \n adjustment",cex.main=0.9)
lines(seq(100,100000,100),seq(100,100000,100))
par(mfrow=c(1,1))

# Compare histograms to evaluate impact of transformation bias 
par(mfrow=c(1,2))
par(mar=c(5.1,5.1,4.1,2.1),mgp=c(4,1,0))
hist(Q_99_at_n_5a_no_trans_adj_urb/Q_99_at_n_med_only_urb,breaks=seq(0.8,1.8,0.1),xlab="100Y flood magnification \n (Trend in mean and Cv/ \n Trend in mean only)",ylim=c(0,10),
     cex.lab=0.8,
     main="No transf bias \n adjustment",
     cex.main=0.9,
     cex.axis=0.85)
hist(Q_99_at_n_5a_urb/Q_99_at_n_med_only_urb,breaks=seq(0.8,1.8,0.1),xlab="100Y flood magnification \n (Trend in mean and Cv/ \n Trend in mean only)",ylim=c(0,10),
     cex.lab=0.8,
     main="With transf bias \n adjustment",
     cex.main=0.9,
     cex.axis=0.85)
par(mfrow=c(1,1))
par(mar=c(5.1,4.1,4.1,2.1),mgp=c(3,1,0))

# Note: Removed magnification factors (see earlier versions)


# FIX: Compare R2 of +mean,-Cv and + mean,+Cv
boxplot(rsquared_5a_vminus_urb,rsquared_5a_urb,
        names=c(expression(paste("+",mu,", -",Cv)), 
                expression(paste("+",mu,", +",Cv))),
        ylab=expression(R^2))

