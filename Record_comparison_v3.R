# Must have NWIS_Qpeak.R open while running this sheet

# Identify promising records
# Use 2/3 residuals test
# Use gamma regression

# EXTRACT RECORDS ---------------------------------------------------------------------

# Create vectors and lists to store output

station_group <- output_mplus_sig05_vminus_sig05

# Conditional mean model
wy_order_mean <- vector(length=nrow(station_group))
wy_order_max <- vector(length=nrow(station_group))
b0_reg1 <- vector(length=nrow(station_group))
b1_reg1 <- vector(length=nrow(station_group))
rho_reg1 <- vector(length=nrow(station_group))
pval_b0_reg1 <- vector(length=nrow(station_group))
pval_b1_reg1 <- vector(length=nrow(station_group))
rsquared_reg1 <- vector(length=nrow(station_group))
res_dof_corr <- vector(length=nrow(station_group))
res_reg1 <- list()
res_var <- vector(length=nrow(station_group))
cond_med_at_n <- vector(length=nrow(station_group))
cond_mean_at_n <- vector(length=nrow(station_group)) # Accounts for transformation bias
Q_99_at_n_med_only <- vector(length=nrow(station_group))
Q_99_at_n_no_trend <- vector(length=nrow(station_group))

# Conditional mean model with studentized residuals


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
cc0_4a <- vector(length=nrow(station_group))
cc1_4a <- vector(length=nrow(station_group))      
rsq_est_4a <- vector(length=nrow(station_group))

# Model 5A (Method of moments, linear trend in variance)
cc0_5a <- vector(length=nrow(station_group))
cc1_5a <- vector(length=nrow(station_group))
rsquared_5a <- vector(length=nrow(station_group))
cond_var_at_n_5a <- vector(length=nrow(station_group))
Q_99_at_n_5a <- vector(length=nrow(station_group))
p_cc0_5a <- vector(length=nrow(station_group))
p_cc1_5a <- vector(length=nrow(station_group)) 

pval_ppcc_res_reg2_5a <- vector(length=nrow(station_group))
pval_dw_res_reg2_5a <- vector(length=nrow(station_group))
p_dd1_5a <- vector(length=nrow(station_group))

# Model 5B (Method of moments, quadratic trend in variance)
cc0_5b <- vector(length=nrow(station_group))
cc1_5b <- vector(length=nrow(station_group))
rsquared_5b <- vector(length=nrow(station_group))
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
cond_var_at_n_5e <- vector(length=nrow(station_group))
Q_99_at_n_5e <- vector(length=nrow(station_group))
p_cc0_5e <- vector(length=nrow(station_group))
p_cc1_5e <- vector(length=nrow(station_group)) 

pval_ppcc_res_reg2_5e <- vector(length=nrow(station_group))
pval_dw_res_reg2_5e <- vector(length=nrow(station_group))
p_dd1_5e <- vector(length=nrow(station_group))

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
cond_var_at_n_6c <- vector(length=nrow(station_group))
Q_99_at_n_6c <- vector(length=nrow(station_group))
p_cc1_6c <- vector(length=nrow(station_group))

# Model 6d (IWLS, logistic trend in variance)
b0_reg1_iwls_6d <- matrix(nrow=nrow(station_group),ncol=10)
b1_reg1_iwls_6d <- matrix(nrow=nrow(station_group),ncol=10)
cc0_6d <- matrix(nrow=nrow(station_group),ncol=10)
cc1_6d <- matrix(nrow=nrow(station_group),ncol=10)

wt_6d_list <- list()
res_reg1_iwls_6d_list <- list()
res_reg1_iwls_2_3_6d_list <- list()
cond_var_6d_list <- list()

rsquared_6d <- vector(length=nrow(station_group))
cond_var_at_n_6d <- vector(length=nrow(station_group))
Q_99_at_n_6d <- vector(length=nrow(station_group))
p_cc1_6d <- vector(length=nrow(station_group))

# GLM - Model A 
b0_glm_cond_var_a <- vector(length=nrow(station_group))
b1_glm_cond_var_a <- vector(length=nrow(station_group)) 
pseudo_rsquared_glm_a <- vector(length=nrow(station_group)) 
cond_var_at_n_glm_a <- vector(length=nrow(station_group)) 
Q_99_at_n_glm_a <- vector(length=nrow(station_group))

# GLM - Model GAMMA
b0_glm_cond_var_gamma <- vector(length=nrow(station_group))
b1_glm_cond_var_gamma <- vector(length=nrow(station_group)) 
pseudo_rsquared_glm_gamma <- vector(length=nrow(station_group)) 
cond_var_at_n_glm_gamma <- vector(length=nrow(station_group)) 
Q_99_at_n_glm_gamma <- vector(length=nrow(station_group))

for (k in 1:nrow(station_group)){ 
wy <- wy_station[[station_group[k,1]]]
wy_order <- wy - min(wy) + 1
wy_order_mean[k] = mean(wy_order)
wy_order_max[k] = max(wy_order)
wyear <- wy_order # Add in water year stuff
peak_flow <- peaks[[station_group[k,1]]]
flood_reg1.lm = lm(log(peak_flow+0.01) ~ wy_order) # Add constant to avoid problems with ephemeral streams
b0_reg1[k] = summary(flood_reg1.lm)$coefficients[1,1]
b1_reg1[k] = summary(flood_reg1.lm)$coefficients[2,1]
rho_reg1[k] = cor(wy_order,log(peak_flow+0.01))
pval_b0_reg1[k] = summary(flood_reg1.lm)$coefficients[1,4]
pval_b1_reg1[k] = summary(flood_reg1.lm)$coefficients[2,4]
rsquared_reg1[k] = as.numeric(summary(flood_reg1.lm)$r.squared)
res_reg1[[k]] = residuals(flood_reg1.lm)
res_dof_corr[k] = length(wy_order)/(length(wy_order)-2)
res_var[k] = res_dof_corr[k] * var(res_reg1[[k]])
cond_med_at_n[k] = b0_reg1[k] + b1_reg1[k] * max(wy_order) 
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
#plot(wy_order,(res_reg1[[k]]^2)^(1/3),xlab="Water year in record (t)",ylab= "Residual Variance ^ (2/3)")
#title(c(site_names[station_group[k,1]],
#     as.character(site_id[station_group[k,1]]),
#     "DA",
#     as.numeric(as.character(site_info[station_group[k,1],6]))),cex.main=0.8)
#text(min(wy_order)+5,max((res_reg1[[k]]^2)^(1/3))-0.08,bquote(~R^2 ==.(round(rsquared_5a[k],3))),cex=0.7)

# View transformed residuals assuming wy_order~ln(res^2)
#log_res_reg1_2 <- log(res_reg1[[k]]^2)
#plot(wy_order,log_res_reg1_2,xlab="Water year in record (t)",ylab= "Residual Variance",pch=4)
#title(c(site_info[station_group[k,1]],
#        as.character(site_info[station_group[k,1]]),
#        "DA",
#        as.numeric(as.character(site_info[station_group[k,1],6]))),cex.main=0.8)

# Fit conditional variance model A using two-stage least squares with an Anscombe transformation
res_reg1_2_3 <- (res_reg1[[k]]^2)^(1/3)
wy_order_1_3 <- wy_order^(1/3)
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
wy_order_2_3 <- wy_order^(2/3)
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
flood_reg2_1a.lm <- lm(res_reg1_2 ~ wy_order)
b0_reg2_1a[k] = summary(flood_reg2_1a.lm)$coefficients[1,1]
b1_reg2_1a[k] = summary(flood_reg2_1a.lm)$coefficients[2,1]
pval_b0_reg2_1a[k] = summary(flood_reg2_1a.lm)$coefficients[1,4]
pval_b1_reg2_1a[k] = summary(flood_reg2_1a.lm)$coefficients[2,4]
adj_rsquared_reg2_1a[k] = as.numeric(summary(flood_reg2_1a.lm)$adj.r.squared)
res_reg2_1a[[k]] = residuals(flood_reg2_1a.lm)

cond_var_at_n_2sls_1a[k] = b0_reg2_1a[k] + b1_reg2_1a[k] * max(wy_order)
Q_99_at_n_2sls_1a[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_2sls_1a[k]))

# Fit conditional variance Model B with two-stage least squares without an Anscombe transformation
res_reg1_2 <- (res_reg1[[k]]^2)
wy_order_2 <- wy_order^2
flood_reg2_1b.lm <- lm(res_reg1_2 ~ wy_order_2)
b0_reg2_1b[k] = summary(flood_reg2_1b.lm)$coefficients[1,1]
b1_reg2_1b[k] = summary(flood_reg2_1b.lm)$coefficients[2,1]
pval_b0_reg2_1b[k] = summary(flood_reg2_1b.lm)$coefficients[1,4]
pval_b1_reg2_1b[k] = summary(flood_reg2_1b.lm)$coefficients[2,4]
adj_rsquared_reg2_1b[k] = as.numeric(summary(flood_reg2_1b.lm)$adj.r.squared)
res_reg2_1b[[k]] = residuals(flood_reg2_1b.lm)

cond_var_at_n_2sls_1b[k] = b0_reg2_1b[k] + b1_reg2_1b[k] * max(wy_order_2)
Q_99_at_n_2sls_1b[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_2sls_1b[k]))

# Fit Model A using 2-parameter derived equation
var_res_reg1 <- length(wy_order)/(length(wy_order)-2)*var(as.numeric(res_reg1[[k]]))
mu_t <- (length(wy_order)+1)/2
cond_var_reg2_2a <- var_res_reg1/mu_t*wy_order
cond_var_at_n_2a[k] = max(cond_var_reg2_2a)
Q_99_at_n_2a[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_2a[k]))

# Fit Model B using 2-parameter derived equation
var_res_reg1 <- length(wy_order)/(length(wy_order)-2)*var(as.numeric(res_reg1[[k]]))
mu_t_2 <- ((length(wy_order)+1)/2)^2
sigma_t_2 <- (length(wy_order)-1)^2/12
cond_var_reg2_2b <- var_res_reg1/(mu_t_2+sigma_t_2)*wy_order^2
cond_var_at_n_2b[k] = max(cond_var_reg2_2b)
Q_99_at_n_2b[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_2b[k]))
# NOT DONE USING ANSCOMBE TRANSFORMATION

# Fit Model A using 3-parameter derived equation
# set initial value of cc1
cc1 = 1
# estimate residual variance 
int <- (length(wy_order)-1)/length(wy_order)*var(as.numeric(res_reg1[[k]]))
# find optimal cc1 value
cc1_min_RSS <- function(cc1) {
  sum((res_reg1_2 - (int + cc1*(wy_order-mean(wy_order))))^2)
}

result <- optimize(cc1_min_RSS,interval=c(-1,1))

cc1_3a[k] <- result$minimum 
RSS <- result$objective
totSS <- sum((res_reg1_2 - mean(res_reg1_2))^2)
rsquared_3a[k] <- 1 - RSS/totSS

# Compute conditional variance and quantile
cond_var_3a <- int + cc1_3a[k]*(wy_order-mean(wy_order))
cond_var_at_n_3a[k] = max(cond_var_3a)
Q_99_at_n_3a[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_3a[k]))

# Fit Model C using 3-parameter derived equation
# set initial value of cc1
cc1 = 1
# estimate residual variance 
int <- (length(wy_order)-1)/length(wy_order)*var(as.numeric(res_reg1[[k]]))
# find optimal cc1 value
cc1_min_RSS <- function(cc1) {
  sum((res_reg1_2 - (int*exp(cc1*(wy_order-mean(wy_order)))))^2)
}

result <- optimize(cc1_min_RSS,interval=c(-1,1))

cc1_3c[k] <- result$minimum 
RSS <- result$objective
totSS <- sum((res_reg1_2 - mean(res_reg1_2))^2)
rsquared_3c[k] <- 1 - RSS/totSS

# Compute conditional variance and quantile
cond_var_3c <- int + cc1_3c[k]*(wy_order-mean(wy_order))
cond_var_at_n_3c[k] = max(cond_var_3c)
Q_99_at_n_3c[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_3c[k]))


# Fit Model A using 4-parameter derived equation
cc <- c(0.5,0.01) 
min.RSS <- function(cc) {
  sum((res_reg1_2_3 - (b1_reg1[k]^2*sigma_t_2*rho_reg1[k])/(cc[1]+cc[2]*mean(wy_order))*(cc[1]+cc[2]*wy_order))^2)
}
optim(cc,min.RSS)
tot.RSS <- sum((res_reg1_2_3 - mean(res_reg1_2_3))^2)

 # Set initial values of c0 and c1
  cc_param_est <- function(cc){
    cc0 <- cc[1]
    cc1 <- cc[2]
    cond_var_res_reg2_4a_init <- (b1_reg1[k]^2*sigma_t_2*rho_reg1[k])/(cc[1]+cc[2]*mean(wy_order))*(cc[1]+cc[2]*wy_order) 
    cond_var_4a_2_3_init <- cond_var_res_reg2_4a_init^(1/3)
    rsq_est <- 1 - sum((res_reg1_2_3-cond_var_4a_2_3_init)^2)/
      (sum((res_reg1_2_3-mean(res_reg1_2_3))^2)) 
  }

cc_param_optim <- optim(cc,cc_param_est,control=list(fnscale=-1)) #fnscale -1 for maximization
cc0_4a[k] = cc_param_optim$par[1]
cc1_4a[k] = cc_param_optim$par[2] 
cond_var_res_reg2_4a_init <- (b1_reg1[k]^2*sigma_t_2*rho_reg1[k])/(cc0_4a[k]+cc1_4a[2]*mean(wy_order))*wy_order 
cond_var_4a_2_3_init <- cond_var_res_reg2_4a_init^(1/3)

rsq_est_4a[k] = 1 - sum((res_reg1_2_3-cond_var_4a_2_3_init)^2)/
  (sum((res_reg1_2_3-mean(res_reg1_2_3))^2)) 

# Fit Model 5A 
rho_t_res_reg2 <- cor(wy_order,res_reg1_2_3)
sd_res_reg1_2_3 <- sd(res_reg1_2_3)
cc0_5a[k] = mean(res_reg1_2_3) - (rho_t_res_reg2 * sd_res_reg1_2_3  * mean(wy_order)) / sd(wy_order)
cc1_5a[k] = rho_t_res_reg2 * sd_res_reg1_2_3 / sd(wy_order)
res_reg2_5a_fit <- (cc0_5a[k] + cc1_5a[k] * wy_order) 
res_reg2_5a_var <- res_dof_corr[k]*var(res_reg1_2_3-res_reg2_5a_fit)
#cond_var_at_n_5a[k] = res_dof_corr[k]*(cc0_5a[k] + cc1_5a[k] * max(wy_order))^3 
cond_var_at_n_5a[k] = res_dof_corr[k]*((cc0_5a[k] + cc1_5a[k] * wy_order_mean[k])^3 + 3*res_reg2_5a_var*(cc0_5a[k] + cc1_5a[k] * wy_order_mean[k]) + mean((res_reg1_2_3-res_reg2_5a_fit)^3) + 
                                         (cov(wy_order,res_reg1_2)/var(wy_order) * (wy_order_max[k] - wy_order_mean[k])))
#plot(wy_order,res_reg1_2_3,col="blue")
#lines(wy_order,res_reg2_5a)
 
#Q_99_at_n_5a[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*(cond_var_at_n_5a[k]^(3/2) + 0.8604*(var(res_reg1_2_3 - res_reg2_5a_fit)^(3/4))))
Q_99_at_n_5a[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5a[k])) 

rsquared_5a[k] = 1-sum((res_reg1_2_3 - res_reg2_5a_fit)^2)/sum((res_reg1_2_3 - mean(res_reg1_2_3))^2)
se_cc0_5a <- sqrt(sum((res_reg1_2_3 - res_reg2_5a_fit)^2)/(length(wy_order) - 2)*(1/length(wy_order)+mean(wy_order)^2/sum((wy_order-mean(wy_order))^2))) 
se_cc1_5a <- sqrt(sum((res_reg1_2_3 - res_reg2_5a_fit)^2)/((length(wy_order) - 2)*sum((wy_order-mean(wy_order))^2)))
t_cc0_5a <- cc0_5a[k]/se_cc0_5a 
t_cc1_5a <- cc1_5a[k]/se_cc1_5a
p_cc0_5a[k] = 2*(pt(-abs(t_cc0_5a), df=length(wy_order) - 1))
p_cc1_5a[k] = 2*(pt(-abs(t_cc1_5a), df=length(wy_order) - 1)) 

# Test residual adequacy 
res_reg2_5a <- as.numeric(res_reg1_2_3 - res_reg2_5a_fit)
ppcc_res_reg2 <- ppcc.test(res_reg2_5a)
pval_ppcc_res_reg2_5a[k] = as.numeric(ppcc_res_reg2[2])
dw_res_reg2 <- dwtest(res_reg2_5a ~ wy_order)
pval_dw_res_reg2_5a[k] = as.numeric(dw_res_reg2[2]) 

# Test for heteroscedasticity of residuals of second regression
res_reg2_2_3_5a <- ((res_reg1_2_3 - res_reg2_5a_fit)^2)^(1/3)
res_reg2_2_3_5a.lm <- lm(res_reg2_2_3_5a ~ wy_order)
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
rho_t_res_reg2 <- cor(wy_order_2,res_reg1_2_3) 
sd_res_reg1_2_3 <- sd(res_reg1_2_3)
cc0_5b[k] = mean(res_reg1_2_3) - (rho_t_res_reg2 * sd_res_reg1_2_3  * mean(wy_order_2) / sd(wy_order_2))
cc1_5b[k] = rho_t_res_reg2 * sd_res_reg1_2_3 / sd(wy_order_2)
res_reg2_5b_fit <- cc0_5b[k] + cc1_5b[k] * wy_order_2
res_reg2_5b_var <- res_dof_corr[k]*var(res_reg1_2_3-res_reg2_5b_fit)
#cond_var_at_n_5b[k] = (cc0_5b[k] + cc1_5b[k] * max(wy_order^2))^3

res_reg2_5b_trans_mean <- res_dof_corr[k]*(cc0_5b[k] + cc1_5b[k] * mean(wy_order_2))^3 + 3*res_reg2_5b_var*(cc0_5b[k] + cc1_5b[k] * mean(wy_order_2)) + mean((res_reg1_2_3-res_reg2_5b_fit)^3)

qq <- var(res_reg1_2_3-res_reg2_5b_fit) #0.0315
rr <- 0.5*(res_reg2_5b_trans_mean - mean((res_reg1_2_3-res_reg2_5b_fit)^3)) #0.5*(0.08086177 - 0.00227) = 0.0382209
ss <- (rr + sqrt(qq^3+rr^2))^(1/3)
tt <- sign(rr - sqrt(qq^3+rr^2))*abs((rr - sqrt(qq^3+rr^2)))^(1/3)
ss_tt <- (ss+tt) 

res_reg2_5b_trans_mean_t <- sqrt(((ss + tt) - cc0_5b[k])/cc1_5b[k]) 
  
res_reg2_5b_trans_cond <- cov(wy_order_2,res_reg1_2)/var(wy_order_2) * (max(wy_order_2) - res_reg2_5b_trans_mean_t^2)
#cond_var_at_n_5b[k] = res_dof_corr[k]*(res_reg2_5b_trans_mean + res_reg2_5b_trans_cond)

cond_var_at_n_5b[k] = res_dof_corr[k]*((cc0_5b[k] + cc1_5b[k] * mean(wy_order_2))^3 + 3*res_reg2_5b_var*(cc0_5b[k] + cc1_5b[k] * mean(wy_order_2)) + res_dof_corr[k]*mean((res_reg1_2_3-res_reg2_5b_fit)^3) +
  cov(wy_order_2,res_reg1_2)/var(wy_order_2) * (max(wy_order_2) - mean(wy_order_2)))



#plot(wy_order,res_reg1_2_3,col="blue")
#lines(wy_order,res_reg2_5b)
Q_99_at_n_5b[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5b[k])) 
rsquared_5b[k] = 1-sum((res_reg1_2_3 - res_reg2_5b_fit)^2)/sum((res_reg1_2_3 - mean(res_reg1_2_3))^2) 
se_cc0_5b <- sqrt(sum((res_reg1_2_3 - res_reg2_5b_fit)^2)/(length(wy_order_2) - 2)*(1/length(wy_order_2)+mean(wy_order_2)^2/sum((wy_order_2-mean(wy_order_2))^2)))
se_cc1_5b <- sqrt(sum((res_reg1_2_3 - res_reg2_5b_fit)^2)/((length(wy_order_2) - 2)*sum((wy_order_2-mean(wy_order_2))^2)))
t_cc0_5b <- cc0_5b[k]/se_cc0_5b 
t_cc1_5b <- cc1_5b[k]/se_cc1_5b 
p_cc0_5b[k] = 2*(pt(-abs(t_cc0_5b), df=length(wy_order^2) - 1))
p_cc1_5b[k] = 2*(pt(-abs(t_cc1_5b), df=length(wy_order^2) - 1)) 

# Test residual adequacy 
res_reg2_5b <- as.numeric(res_reg1_2_3 - res_reg2_5b_fit)
ppcc_res_reg2 <- ppcc.test(res_reg2_5b)
pval_ppcc_res_reg2_5b[k] = as.numeric(ppcc_res_reg2[2])
dw_res_reg2 <- dwtest(res_reg2_5b ~ wy_order)
pval_dw_res_reg2_5b[k] = as.numeric(dw_res_reg2[2]) 

# Test for heteroscedasticity of residuals of second regression
res_reg2_2_3_5b <- ((res_reg1_2_3 - res_reg2_5b_fit)^2)^(1/3)
res_reg2_2_3_5b.lm <- lm(res_reg2_2_3_5b ~ wy_order)
p_dd1_5b[k] = summary(res_reg2_2_3_5b.lm)$coefficients[2,4]

# Fit Model 5C (exponential)
exp_wy_order_nrmlz <- exp(wy_order/max(wy_order))
rho_t_res_reg2 <- cor(exp_wy_order_nrmlz,res_reg1_2_3)
sd_res_reg1_2_3 <- sd(res_reg1_2_3)
cc0_5c[k] = mean(res_reg1_2_3) - (rho_t_res_reg2 * sd_res_reg1_2_3  * mean(exp_wy_order_nrmlz) / sd(exp_wy_order_nrmlz))
cc1_5c[k] = rho_t_res_reg2 * sd_res_reg1_2_3 / sd(exp_wy_order_nrmlz)
res_reg2_5c_fit <- cc0_5c[k] + cc1_5c[k] * exp_wy_order_nrmlz
res_reg2_5c_var <- res_dof_corr[k] * var(res_reg1_2_3-res_reg2_5c_fit) 
#cond_var_at_n_5c[k] = cc0_5c[k] + cc1_5c[k] * max(wy_order)^3 # ADD BIAS CORRECTION 
cond_var_at_n_5c[k] = res_dof_corr[k]*((cc0_5c[k] + cc1_5c[k] * mean(exp_wy_order_nrmlz))^3 + 3*res_reg2_5c_var*(cc0_5c[k] + cc1_5c[k] * mean(exp_wy_order_nrmlz)) + mean((res_reg1_2_3-res_reg2_5c_fit)^3) +
                                         cov(exp_wy_order_nrmlz,res_reg1_2)/var(exp_wy_order_nrmlz) * (max(exp_wy_order_nrmlz) - mean(exp_wy_order_nrmlz)))


Q_99_at_n_5c[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5c[k])) 
rsquared_5c[k] = 1-sum((res_reg1_2_3 - res_reg2_5c_fit)^2)/sum((res_reg1_2_3 - mean(res_reg1_2_3))^2) 
se_cc0_5c <- sqrt(sum((res_reg1_2_3 - res_reg2_5c_fit)^2)/(length(wy_order) - 2)*(1/length(wy_order)+mean(exp_wy_order_nrmlz)^2/sum((exp_wy_order_nrmlz-mean(exp_wy_order_nrmlz))^2)))
se_cc1_5c <- sqrt(sum((res_reg1_2_3 - res_reg2_5c_fit)^2)/((length(wy_order) - 2)*sum((exp_wy_order_nrmlz-mean(exp_wy_order_nrmlz))^2)))
t_cc0_5c <- cc0_5c[k]/se_cc0_5c
t_cc1_5c <- cc1_5c[k]/se_cc1_5c
p_cc0_5c[k] = 2*(pt(-abs(t_cc0_5c), df=length(wy_order) - 1))
p_cc1_5c[k] = 2*(pt(-abs(t_cc1_5c), df=length(wy_order) - 1)) 

# Test residual adequacy 
res_reg2_5c <- as.numeric(res_reg1_2_3 - res_reg2_5c_fit)
ppcc_res_reg2 <- ppcc.test(res_reg2_5c)
pval_ppcc_res_reg2_5c[k] = as.numeric(ppcc_res_reg2[2])
dw_res_reg2 <- dwtest(res_reg2_5c ~ wy_order)
pval_dw_res_reg2_5c[k] = as.numeric(dw_res_reg2[2]) 

# Test for heteroscedasticity of residuals of second regression
res_reg2_2_3_5c <- ((res_reg1_2_3 - res_reg2_5c_fit)^2)^(1/3)
res_reg2_2_3_5c.lm <- lm(res_reg2_2_3_5c ~ wy_order)
p_dd1_5c[k] = summary(res_reg2_2_3_5c.lm)$coefficients[2,4]

# Fit Model 5d (logarithmic)
rho_t_res_reg2 <- cor(log(wy_order),res_reg1_2_3)
sd_res_reg1_2_3 <- sd(res_reg1_2_3)
cc0_5d[k] = mean(res_reg1_2_3) - (rho_t_res_reg2 * sd_res_reg1_2_3  * mean(log(wy_order)) / sd(log(wy_order)))
cc1_5d[k] = rho_t_res_reg2 * sd_res_reg1_2_3 / sd(log(wy_order))
res_reg2_5d_fit <- cc0_5d[k] + cc1_5d[k] * log(wy_order)
res_reg2_5d_var <- length(wy_order)/(length(wy_order)-2)*var(res_reg1_2_3-res_reg2_5d_fit)
cond_var_at_n_5d[k] = res_dof_corr[k]*(cc0_5d[k] + cc1_5d[k] * mean(log(wy_order)))^3 + res_reg2_5d_var*(3*cc0_5d[k] + 3*cc1_5d[k] * mean(log(wy_order))) + mean((res_reg1_2_3-res_reg2_5d_fit)^3) +
  cov(log(wy_order),res_reg1_2)/var(log(wy_order)) * (max(log(wy_order)) - mean(log(wy_order)))
#plot(wy_order,res_reg1_2_3,col="blue")
#lines(wy_order,res_reg2_5d)
Q_99_at_n_5d[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5d[k]))
rsquared_5d[k] = 1-sum((res_reg1_2_3 - res_reg2_5d_fit)^2)/sum((res_reg1_2_3 - mean(res_reg1_2_3))^2)
se_cc0_5d <- sqrt(sum((res_reg1_2_3 - res_reg2_5d_fit)^2)/(length(wy_order) - 2)*(1/length(wy_order)+mean(log(wy_order))^2/sum((log(wy_order)-mean(log(wy_order)))^2))) 
se_cc1_5d <- sqrt(sum((res_reg1_2_3 - res_reg2_5d_fit)^2)/((length(wy_order) - 2)*sum((log(wy_order)-mean(log(wy_order)))^2)))
t_cc0_5d <- cc0_5d[k]/se_cc0_5d
t_cc1_5d <- cc1_5d[k]/se_cc1_5d
p_cc0_5d[k] = 2*(pt(-abs(t_cc0_5d), df=length(wy_order) - 1))
p_cc1_5d[k] = 2*(pt(-abs(t_cc1_5d), df=length(wy_order) - 1)) 

# Test residual adequacy 
res_reg2_5d <- as.numeric(res_reg1_2_3 - res_reg2_5d_fit)
ppcc_res_reg2 <- ppcc.test(res_reg2_5d)
pval_ppcc_res_reg2_5d[k] = as.numeric(ppcc_res_reg2[2])
dw_res_reg2 <- dwtest(res_reg2_5d ~ wy_order)
pval_dw_res_reg2_5d[k] = as.numeric(dw_res_reg2[2]) 

# Test for heteroscedasticity of residuals of second regression
res_reg2_2_3_5d <- ((res_reg1_2_3 - res_reg2_5d_fit)^2)^(1/3)
res_reg2_2_3_5d.lm <- lm(res_reg2_2_3_5d ~ wy_order)
p_dd1_5d[k] = summary(res_reg2_2_3_5d.lm)$coefficients[2,4]

# Model 5E (Method of moments, log-transformed squared residuals)
log_res_reg1_2 <- log(res_reg1[[k]]^2)
rho_t_res_reg2 <- cor(wy_order,log_res_reg1_2)
sd_log_res_reg1_2 <- sd(log_res_reg1_2)
cc0_5e[k] = mean(log_res_reg1_2) - (rho_t_res_reg2 * sd_log_res_reg1_2  * mean(wy_order) / sd(wy_order))
cc1_5e[k] = rho_t_res_reg2 * sd_log_res_reg1_2 / sd(wy_order)
res_reg2_5e_fit <- cc0_5e[k] + cc1_5e[k] * wy_order
res_reg2_5e_var <- var(log_res_reg1_2-res_reg2_5e_fit)
cond_var_at_n_5e[k] = exp(res_dof_corr[k]*((cc0_5e[k] + cc1_5e[k] * max(wy_order) + 0.5*res_reg2_5e_var)
                                           + cov(wy_order,log_res_reg1_2)/var(wy_order)*(max(wy_order)-mean(wy_order))))  
#plot(wy_order,log_res_reg1_2,col="blue")
#lines(wy_order,log_res_reg2_5e) 
Q_99_at_n_5e[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5e[k])) 
rsquared_5e[k] = 1-sum((log_res_reg1_2 - res_reg2_5e_fit)^2)/sum((log_res_reg1_2 - mean(log_res_reg1_2))^2) 
se_cc0_5e <- sqrt(sum((log_res_reg1_2 - res_reg2_5e_fit)^2)/(length(wy_order) - 2)*(1/length(wy_order)+mean(wy_order)^2/sum((wy_order-mean(wy_order))^2))) 
se_cc1_5e <- sqrt(sum((log_res_reg1_2 - res_reg2_5e_fit)^2)/((length(wy_order) - 2)*sum((wy_order-mean(wy_order))^2))) 
t_cc0_5e <- cc0_5e[k]/se_cc0_5e 
t_cc1_5e <- cc1_5e[k]/se_cc1_5e
p_cc0_5e[k] = 2*(pt(-abs(t_cc0_5e), df=length(wy_order) - 1)) 
p_cc1_5e[k] = 2*(pt(-abs(t_cc1_5e), df=length(wy_order) - 1)) 

# Test residual adequacy 
res_reg2_5e <- as.numeric(log_res_reg1_2 - res_reg2_5e_fit) 
ppcc_res_reg2 <- ppcc.test(res_reg2_5e) 
pval_ppcc_res_reg2_5e[k] = as.numeric(ppcc_res_reg2[2]) 
dw_res_reg2 <- dwtest(res_reg2_5e ~ wy_order) 
pval_dw_res_reg2_5e[k] = as.numeric(dw_res_reg2[2]) 

# Test for heteroscedasticity of residuals of second regression
log_res_reg2_5e <- ((log_res_reg1_2 - res_reg2_5e_fit)^2)^(1/3)
log_res_reg2_5e.lm <- lm(log_res_reg2_5e ~ wy_order)
p_dd1_5e[k] = summary(log_res_reg2_5e.lm)$coefficients[2,4]

# MODEL 6: ITERATIVE RE(WEIGHTED) LEAST SQUARES --------------------------------------------------

# Fit Model 6A 

# Establish output arrays for each record (length-dependent)
wt_6a <- array(NA,dim=c(length(wy_order),10))
res_reg1_iwls_6a <- array(NA,dim=c(length(wy_order),10))
res_reg1_iwls_2_3_6a <- array(NA,dim=c(length(wy_order),10))
cond_var_6a <- array(NA,dim=c(length(wy_order),10))

for (i in 1:10){
  if(i == 1){
    wt_6a[,i] = rep(1,length(wy_order)) #Make a list
  } else {
    wt_6a[,i] = 1/(cond_var_6a[,i-1]) #Make a list
  }
  flood_reg1_iwls.lm <- lm(log(peak_flow + 0.01) ~ wy_order, weights=wt_6a[,i]) # Make a list
  b0_reg1_iwls_6a[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[1])
  b1_reg1_iwls_6a[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[2])
  res_reg1_iwls_6a[,i] = as.numeric(residuals(flood_reg1_iwls.lm)) # Make a list
  res_reg1_iwls_2_3_6a[,i] = (as.numeric(residuals(flood_reg1_iwls.lm)^2))^(1/3)# Make a list
  res_reg1_iwls_2_3_6a.lm <- lm(res_reg1_iwls_2_3_6a[,i] ~ wy_order) # Make a list
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

cond_var_at_n_6a[k] = res_dof_corr[k]*(cc0_6a[k,10] + cc1_6a[k,10] * mean(wy_order))^3 + 3*res_reg2_6a_var*(cc0_6a[k,10] + cc1_6a[k,10] * mean(wy_order)) + mean((res_reg1_2_3-cond_var_6a[,10])^3) + 
  cov(wy_order,res_reg1_iwls_6a[,10]^2)/var(wy_order) * (max(wy_order) - mean(wy_order))

#plot(wy_order,res_reg1_2_3,col="blue")
#lines(wy_order,res_reg2_5a)
Q_99_at_n_6a[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_6a[k]))
rsquared_6a[k] = 1-sum((res_reg1_iwls_2_3_6a[,max(i)] - cond_var_6a[,max(i)])^2)/sum((res_reg1_iwls_2_3_6a[,max(i)] - mean(res_reg1_iwls_2_3_6a[,max(i)]))^2)
se_cc1_6a <- sqrt(sum((res_reg1_iwls_2_3_6a[,max(i)] - cond_var_6a[,max(i)])^2)/((length(wy_order) - 2)*sum((wy_order-mean(wy_order))^2)))
t_cc1_6a <- cc1_6a[k]/se_cc1_6a
p_cc1_6a[k] = 2*(pt(-abs(t_cc1_6a), df=length(wy_order) - 1)) 

# Fit model 6B (Quadratic with IWLS)

# Establish output arrays for each record (length-dependent)
wt_6b <- array(NA,dim=c(length(wy_order),10))
res_reg1_iwls_6b <- array(NA,dim=c(length(wy_order),10))
res_reg1_iwls_2_3_6b <- array(NA,dim=c(length(wy_order),10))
cond_var_6b <- array(NA,dim=c(length(wy_order),10))

for (i in 1:10){
  if(i == 1){
    wt_6b[,i] = rep(1,length(wy_order)) #Make a list
  } else {
    wt_6b[,i] = 1/(cond_var_6b[,i-1]) #Make a list
  }
  flood_reg1_iwls.lm <- lm(log(peak_flow + 0.01) ~ wy_order, weights=wt_6b[,i]) # Make a list
  b0_reg1_iwls_6b[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[1])
  b1_reg1_iwls_6b[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[2])
  res_reg1_iwls_6b[,i] = as.numeric(residuals(flood_reg1_iwls.lm)) # Make a list
  res_reg1_iwls_2_3_6b[,i] = (as.numeric(residuals(flood_reg1_iwls.lm)^2))^(1/3)# Make a list
  res_reg1_iwls_2_3_6b.lm <- lm(res_reg1_iwls_2_3_6b[,i] ~ wy_order_2) # Make a list
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

cond_var_at_n_6b[k] = (cc0_6b[k,10] + cc1_6b[k,10] * mean(wy_order_2))^3 + 3*res_reg2_6b_var*(cc0_6b[k,10] + cc1_6b[k,10] * mean(wy_order_2)) + mean(((res_reg1_iwls_2_3_6b[,10]-cond_var_6b[,10]))^3) +
  cov(wy_order_2,res_reg1_iwls_6b[,10]^2)/var(wy_order_2) * (max(wy_order_2) - mean(wy_order_2))

#plot(wy_order,res_reg1_2_3,col="blue")
#lines(wy_order,res_reg2_5a)
Q_99_at_n_6b[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_6b[k]))
rsquared_6b[k] = 1-sum((res_reg1_iwls_2_3_6b[,max(i)] - cond_var_6b[,max(i)])^2)/sum((res_reg1_iwls_2_3_6b[,max(i)] - mean(res_reg1_iwls_2_3_6b[,max(i)]))^2)
se_cc1_6b <- sqrt(sum((res_reg1_iwls_2_3_6b[,max(i)] - cond_var_6b[,max(i)])^2)/((length(wy_order) - 2)*sum((wy_order^2-mean(wy_order^2))^2)))
t_cc1_6b <- cc1_6b[k]/se_cc1_6b
p_cc1_6b[k] = 2*(pt(-abs(t_cc1_6b), df=length(wy_order) - 1)) 

# Fit model 6c (Exponential with IWLS)

# Establish output arrays for each record (length-dependent)
wt_6c <- array(NA,dim=c(length(wy_order),10))
res_reg1_iwls_6c <- array(NA,dim=c(length(wy_order),10))
res_reg1_iwls_2_3_6c <- array(NA,dim=c(length(wy_order),10))
cond_var_6c <- array(NA,dim=c(length(wy_order),10))

for (i in 1:10){
  if(i == 1){
    wt_6c[,i] = rep(1,length(wy_order)) #Make a list
  } else {
    wt_6c[,i] = 1/(cond_var_6c[,i-1]) #Make a list
  }
  flood_reg1_iwls.lm <- lm(log(peak_flow + 0.01) ~ wy_order, weights=wt_6c[,i]) # Make a list
  b0_reg1_iwls_6c[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[1])
  b1_reg1_iwls_6c[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[2])
  res_reg1_iwls_6c[,i] = as.numeric(residuals(flood_reg1_iwls.lm)) # Make a list
  res_reg1_iwls_2_3_6c[,i] = (as.numeric(residuals(flood_reg1_iwls.lm)^2))^(1/3)# Make a list
  res_reg1_iwls_2_3_6c.lm <- lm(res_reg1_iwls_2_3_6c[,i] ~ exp(wy_order)) # Make a list
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
cond_var_at_n_6c[k] = (exp(cc0_6c[k,10] + cc1_6c[k,10] * max(wy_order)))^3
#cond_var_at_n_6c[k] = (cc0_6c[k,10] + cc1_6c[k,10] * max(exp(wy_order)))^3 + 
#  var(res_reg1_iwls_2_3_6c[,10]-cond_var_6c[,10])*(3*cc0_6c[k,10] + 3*cc1_6c[k,10] * max(exp(wy_order))) + mean((res_reg1_iwls_2_3_6c[,10]-cond_var_6c[,10])^3)
#plot(wy_order,res_reg1_2_3,col="blue")
#lines(wy_order,res_reg2_5a)
Q_99_at_n_6c[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_6c[k]))
rsquared_6c[k] = 1-sum((res_reg1_iwls_2_3_6c[,max(i)] - cond_var_6c[,max(i)])^2)/sum((res_reg1_iwls_2_3_6c[,max(i)] - mean(res_reg1_iwls_2_3_6c[,max(i)]))^2)
se_cc1_6c <- sqrt(sum((res_reg1_iwls_2_3_6c[,max(i)] - cond_var_6c[,max(i)])^2)/((length(wy_order) - 2)*sum((exp(wy_order)-mean(exp(wy_order)))^2)))
t_cc1_6c <- cc1_6c[k]/se_cc1_6c
p_cc1_6c[k] = 2*(pt(-abs(t_cc1_6c), df=length(wy_order) - 1)) 


# Fit model 6d (Logarithmic with IWLS)

# Establish output arrays for each record (length-dependent)
wt_6d <- array(NA,dim=c(length(wy_order),10))
res_reg1_iwls_6d <- array(NA,dim=c(length(wy_order),10))
res_reg1_iwls_2_3_6d <- array(NA,dim=c(length(wy_order),10))
cond_var_6d <- array(NA,dim=c(length(wy_order),10))

for (i in 1:10){
  if(i == 1){
    wt_6d[,i] = rep(1,length(wy_order)) #Make a list
  } else {
    wt_6d[,i] = 1/(cond_var_6d[,i-1]) #Make a list
  }
  flood_reg1_iwls.lm <- lm(log(peak_flow + 0.01) ~ wy_order, weights=wt_6d[,i]) # Make a list
  b0_reg1_iwls_6d[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[1])
  b1_reg1_iwls_6d[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[2])
  res_reg1_iwls_6d[,i] = as.numeric(residuals(flood_reg1_iwls.lm)) # Make a list
  res_reg1_iwls_2_3_6d[,i] = (as.numeric(residuals(flood_reg1_iwls.lm)^2))^(1/3)# Make a list
  res_reg1_iwls_2_3_6d.lm <- lm(res_reg1_iwls_2_3_6d[,i] ~ log(wy_order)) # Make a list
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

cond_var_at_n_6d[k] = (cc0_6d[k,10] + cc1_6d[k,10] * mean(log(wy_order)))^3 + res_reg2_6d_var*(3*cc0_6d[k,10] + 3*cc1_6d[k,10] * mean(log(wy_order))) + mean((res_reg1_iwls_2_3_6d[,10]-cond_var_6d[,10])^3) +
  cov(log(wy_order),res_reg1_iwls_6d[,10]^2)/var(log(wy_order)) * (max(log(wy_order)) - mean(log(wy_order)))

#plot(wy_order,res_reg1_2_3,col="blue")
#lines(wy_order,res_reg2_5a)
Q_99_at_n_6d[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_6d[k]))
rsquared_6d[k] = 1-sum((res_reg1_iwls_2_3_6d[,max(i)] - cond_var_6d[,max(i)])^2)/sum((res_reg1_iwls_2_3_6d[,max(i)] - mean(res_reg1_iwls_2_3_6d[,max(i)]))^2)
se_cc1_6d <- sqrt(sum((res_reg1_iwls_2_3_6d[,max(i)] - cond_var_6d[,max(i)])^2)/((length(wy_order) - 2)*sum((log(wy_order)-mean(log(wy_order)))^2)))
t_cc1_6d <- cc1_6d[k]/se_cc1_6d 
p_cc1_6d[k] = 2*(pt(-abs(t_cc1_6d), df=length(wy_order) - 1))  

# Fit Model A using GLM 
glm.cond_var_a <- glm(formula = res_reg1_2_3 ~ wy_order_1_3, family=gaussian) 
b0_glm_cond_var_a[k] = as.numeric(glm.cond_var_a$coefficients[1]) 
b1_glm_cond_var_a[k] = as.numeric(glm.cond_var_a$coefficients[2]) 
res_dev <- as.numeric(glm.cond_var_a$deviance) 
null_dev <- as.numeric(glm.cond_var_a$null.deviance) 
pseudo_rsquared_glm_a[k] = 1 - res_dev/null_dev 
cond_var_at_n_glm_a[k] = (b0_glm_cond_var_a[k] + b1_glm_cond_var_a[k] * max(wy_order_1_3))^3
Q_99_at_n_glm_a[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_glm_a[k]))

# Fit model using GLM with gamma regression 
res_reg1_2 <- as.numeric((res_reg1[[k]])^2)
glm.cond_var_gamma <- glm(formula = res_reg1_2 ~ wy_order, family=Gamma(link=log),start=c(log(mean(res_reg1_2)),0))
b0_glm_cond_var_gamma[k] = as.numeric(glm.cond_var_gamma$coefficients[1])
b1_glm_cond_var_gamma[k] = as.numeric(glm.cond_var_gamma$coefficients[2])
res_dev <- as.numeric(glm.cond_var_gamma$deviance)
null_dev <- as.numeric(glm.cond_var_gamma$null.deviance) # intercept-only model
pseudo_rsquared_glm_gamma[k] = 1 - res_dev/null_dev
cond_var_at_n_glm_gamma[k] = exp(b0_glm_cond_var_gamma[k] + b1_glm_cond_var_gamma[k] * max(wy_order)) #Transformation bias?
Q_99_at_n_glm_gamma[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_glm_gamma[k]))

}



# Compare models by R^2

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
              Q_99_at_n_glm_gamma)

boxplot(rsquared_5a,
        rsquared_5b,
        #rsquared_5c,
        rsquared_5d,
        #rsquared_5e,
        pseudo_rsquared_glm_gamma)

boxplot((Q_99_at_n_6a-Q_99_at_n_5a)/Q_99_at_n_5a,
        (Q_99_at_n_6b-Q_99_at_n_5b)/Q_99_at_n_5b,
        (Q_99_at_n_6d-Q_99_at_n_5d)/Q_99_at_n_5d)