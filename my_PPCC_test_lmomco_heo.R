# Test residual normality -------

ppcc.test_jory <- function(y,res) {

  # Convert to L-moments in list format
  lmoms_res <- lmoms(res)
  
  # Estimate parameters from residuals
  par_nor <- parnor(lmoms_res)
  par_nor_valid <- are.parnor.valid(par_nor,nowarn=FALSE)
  
  # Compute Weibull probabiities for PPCC test
  y_ranked <- sort(y)
  p_y <- rank(y_ranked)/(length(y_ranked)+1)
  
  # Compute PPCC value
  F_nor_est <- quanor(p_y,par_nor)
  F_nor_QQ_PPCC <- cor(F_nor_est,y_ranked)
  plot(F_nor_est,y_ranked)
}

# Convert to L-moments in list format
lmom_residuals_fit_ols <- lmoms(residuals_fit_ols)

# Estimate parameters from residuals
par_nor <- parnor(lmom_residuals_fit_ols)
par_nor_valid <- are.parnor.valid(par_nor,nowarn=FALSE)

# Compute Weibull probabilities for PPCC test
qpeak_99_ranked <- sort(qpeak_99)
p_q99_peak <- rank(qpeak_99_ranked)/(length(qpeak_99_ranked)+1)

# Compute PPCC value
F_nor_est <- quanor(p_q99_peak,par_nor)
F_nor_QQ_PPCC <- cor(F_nor_est,qpeak_99_ranked)
plot(F_nor_est,qpeak_99_ranked)

# Find critical PPCC value
ppcc_sig_level = 0.10

# valid for sig level 0.005 - 0.10
rhs_heo_ppcc_1_10 <- 2.54*ppcc_sig_level^0.146*length(qpeak_99)^(0.152-0.00993*log(ppcc_sig_level))
ppcc_crit_1_10 <- 1-1/exp(rhs_heo_ppcc_1_10)

# valid for sig level 0.90 - 0.995
rhs_heo_ppcc_90_99 <- 1/(-5.22+12.3*ppcc_sig_level-6.81*ppcc_sig_level^2)*
length(qpeak_99)^(-2.42+5.63*ppcc_sig_level-3.11*ppcc_sig_level^2)
ppcc_crit_90_99 <- 1-1/exp(rhs_heo_ppcc_90_99)  
