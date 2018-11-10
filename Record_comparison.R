# Must have NWIS_Qpeak.R open while running this sheet

# Identify promising records
# Use 2/3 residuals test
# Use gamma regression

# EXTRACT RECORDS ---------------------------------------------------------------------

# Create vectors and lists to store output

station_group <- output_mplus_sig05_vminus_sig05
# Conditional mean model
b0_reg1 <- vector(length=nrow(station_group))
b1_reg1 <- vector(length=nrow(station_group))
rho_reg1 <- vector(length=nrow(station_group))
pval_b0_reg1 <- vector(length=nrow(station_group))
pval_b1_reg1 <- vector(length=nrow(station_group))
rsquared_reg1 <- vector(length=nrow(station_group))
res_reg1 <- list()
cond_med_at_n <- vector(length=nrow(station_group))
cond_mean_at_n <- vector(length=nrow(station_group)) # Accounts for transformation bias
Q_99_at_n_med_only <- vector(length=nrow(station_group))
Q_99_at_n_no_trend <- vector(length=nrow(station_group))

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
p_dd1_5a <- vector(length=nrow(station_group))

# Model 5B (Method of moments, quadratic trend in variance)
cc0_5b <- vector(length=nrow(station_group))
cc1_5b <- vector(length=nrow(station_group))
rsquared_5b <- vector(length=nrow(station_group))
cond_var_at_n_5b <- vector(length=nrow(station_group))
Q_99_at_n_5b <- vector(length=nrow(station_group))
p_cc0_5b <- vector(length=nrow(station_group))
p_cc1_5b <- vector(length=nrow(station_group)) 

# Model 5C (Method of moments, log trend in variance)
cc0_5c<- vector(length=nrow(station_group))
cc1_5c <- vector(length=nrow(station_group))
rsquared_5c <- vector(length=nrow(station_group))
cond_var_at_n_5c <- vector(length=nrow(station_group))
Q_99_at_n_5c <- vector(length=nrow(station_group))
p_cc0_5c <- vector(length=nrow(station_group))
p_cc1_5c <- vector(length=nrow(station_group)) 

# Model 6A (IWLS, linear trend in variance)
b0_reg1_iwls <- matrix(nrow=nrow(station_group),ncol=10)
b1_reg1_iwls <- matrix(nrow=nrow(station_group),ncol=10)
cc0_6a <- matrix(nrow=nrow(station_group),ncol=10)
cc1_6a <- matrix(nrow=nrow(station_group),ncol=10)

wt_list <- list()
res_reg1_iwls_list <- list()
res_reg1_iwls_2_list <- list()
cond_var_6a_list <- list()

rsquared_6a <- vector(length=nrow(station_group))
cond_var_at_n_6a <- vector(length=nrow(station_group))
Q_99_at_n_6a <- vector(length=nrow(station_group))
p_cc1_6a <- vector(length=nrow(station_group)) 



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
cond_med_at_n[k] = b0_reg1[k] + b1_reg1[k] * max(wy_order) 
cond_mean_at_n[k] = cond_med_at_n[k] + 0.5*var(res_reg1[[k]])
Q_99_at_n_med_only[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1) * sd(res_reg1[[k]]))
Q_99_at_n_no_trend[k] = exp(mean(log(peak_flow + 0.01)) + qnorm(0.99,0,1) * sd(log(peak_flow+0.01)))

# View conditional mean model
plot(wy,log(peak_flow+0.01),xlab="Water year",ylab= "ln(Annual peak flow, cfs)")
title(c(site_names[station_group[k,1]],
      as.character(sites[station_group[k,1]]),
      "DA",
      as.numeric(as.character(site_info[station_group[k,1],6]))),cex.main=0.8)
text(min(wy)+5,max(log(peak_flow+0.01))-0.5,bquote(~R^2 ==. (round(rsquared_reg1[k],3))),cex=0.7)


# View transformed residuals assuming wy_order~res
#plot(wy_order,(res_reg1[[k]]^2)^(1/3),xlab="Water year in record (t)",ylab= "Residual Variance ^ (2/3)")
#title(c(site_names[station_group[k,1]],
#      as.character(sites[station_group[k,1]]),
#     "DA",
#     as.numeric(as.character(site_info[station_group[k,1],6]))),cex.main=0.8)
#text(min(wy_order)+5,max((res_reg1[[k]]^2)^(1/3)-0.08,bquote(~R^2 ==. (round(rsquared_5a[k],3))),cex=0.7)

# View transformed residuals assuming wy_order~ln(res^2)
#log_res_reg1_2 <- log(res_reg1[[k]]^2)
#plot(wy_order,log_res_reg1_2,xlab="Water year in record (t)",ylab= "Residual Variance",pch=4)
#title(c(site_names[station_group[k,1]],
#        as.character(sites[station_group[k,1]]),
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
#        as.character(sites[station_group[k,1]]),
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

# View transformed residuals
#plot(wy_order_2_3,(res_reg1[[k]]^2)^(2/3),xlab="Water year in record (t)",ylab= "Residual Variance",pch=2)
#title(c(site_names[station_group[k,1]],
#        as.character(sites[station_group[k,1]]),
#        "DA",
#        as.numeric(as.character(site_info[station_group[k,1],6]))),cex.main=0.8)

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
var(as.numeric(res_reg2_1a[[k]]))

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
var_res_reg1 <- var(as.numeric(res_reg1[[k]]))
mu_t <- (length(wy_order)+1)/2
cond_var_reg2_2a <- var_res_reg1/mu_t*wy_order
cond_var_at_n_2a[k] = max(cond_var_reg2_2a)
Q_99_at_n_2a[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_2a[k]))

# Fit Model B using 2-parameter derived equation
var_res_reg1 <- var(as.numeric(res_reg1[[k]]))
mu_t_2 <- ((length(wy_order)+1)/2)^2
sigma_t_2 <- (length(wy_order)-1)^2/12
cond_var_reg2_2b <- var_res_reg1/(mu_t_2+sigma_t_2)*wy_order^2
cond_var_at_n_2b[k] = max(cond_var_reg2_2b)
Q_99_at_n_2b[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_2b[k]))


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
res_reg2_5a_fit <- cc0_5a[k] + cc1_5a[k] * wy_order
cond_var_at_n_5a[k] = (cc0_5a[k] + cc1_5a[k] * max(wy_order))^3 + var(res_reg1_2_3-res_reg2_5a_fit)*(3*cc0_5a[k] + 3*cc1_5a[k] * max(wy_order)) + mean((res_reg1_2_3-res_reg2_5a_fit)^3)
#plot(wy_order,res_reg1_2_3,col="blue")
#lines(wy_order,res_reg2_5a)
#Q_99_at_n_5a[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*cond_var_at_n_5a[k]^(3/2))
#Q_99_at_n_5a[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*(cond_var_at_n_5a[k]^(3/2) + 0.8604*(var(res_reg1_2_3 - res_reg2_5a_fit)^(3/4))))
Q_99_at_n_5a[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5a[k]))

rsquared_5a[k] = 1-sum((res_reg1_2_3 - res_reg2_5a_fit)^2)/sum((res_reg1_2_3 - mean(res_reg1_2_3))^2)
se_cc0_5a <- sqrt(sum((res_reg1_2_3 - res_reg2_5a_fit)^2)/(length(wy_order) - 2)*(1/length(wy_order)+mean(wy_order)^2/sum((wy_order-mean(wy_order))^2)))
se_cc1_5a <- sqrt(sum((res_reg1_2_3 - res_reg2_5a_fit)^2)/((length(wy_order) - 2)*sum((wy_order-mean(wy_order))^2)))
t_cc0_5a <- cc0_5a[k]/se_cc0_5a
t_cc1_5a <- cc1_5a[k]/se_cc1_5a
p_cc0_5a[k] = 2*(pt(-abs(t_cc0_5a), df=length(wy_order) - 1))
p_cc1_5a[k] = 2*(pt(-abs(t_cc1_5a), df=length(wy_order) - 1)) 

# Test for heteroscedasticity of residuals of second regression
res_reg2_2_3_5a <- ((res_reg1_2_3 - res_reg2_5a_fit)^2)^(1/3)
res_reg2_2_3_5a.lm <- lm(res_reg2_2_3_5a ~ wy_order)
p_dd1_5a[k] = summary(res_reg2_2_3_5a.lm)$coefficients[2,4]


# View transformed residuals assuming wy_order~res
#plot(wy_order,(res_reg1[[k]]^2)^(1/3),xlab="Water year in record (t)",ylab= "Residual Variance ^ (2/3)")
# title(c(site_names[station_group[k,1]],
 #      as.character(sites[station_group[k,1]]),
 #      "DA",
 #      as.numeric(as.character(site_info[station_group[k,1],6]))),cex.main=0.8) 
#text(max(wy_order)-10,max((res_reg1[[k]]^2)^(1/3))-0.08,bquote(~R^2 ==. (round(rsquared_5a[k],3))),cex=0.7)
# max((res_reg1[[k]]^2)^(1/3))-0.08   


# Fit Model 5B 
rho_t_res_reg2 <- cor(wy_order^2,res_reg1_2_3) 
sd_res_reg1_2_3 <- sd(res_reg1_2_3)
cc0_5b[k] = mean(res_reg1_2_3) - (rho_t_res_reg2 * sd_res_reg1_2_3  * mean(wy_order^2) / sd(wy_order^2))
cc1_5b[k] = rho_t_res_reg2 * sd_res_reg1_2_3 / sd(wy_order^2)
res_reg2_5b_fit <- cc0_5b[k] + cc1_5b[k] * (wy_order^2)
cond_var_at_n_5b[k] = cc0_5b[k] + cc1_5b[k] * max(wy_order^2)
#plot(wy_order,res_reg1_2_3,col="blue")
#lines(wy_order,res_reg2_5c)
Q_99_at_n_5b[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5b[k]^3))
rsquared_5b[k] = 1-sum((res_reg1_2_3 - res_reg2_5b_fit)^2)/sum((res_reg1_2_3 - mean(res_reg1_2_3))^2)
se_cc0_5b <- sqrt(sum((res_reg1_2_3 - res_reg2_5b_fit)^2)/(length(wy_order_2) - 2)*(1/length(wy_order_2)+mean(wy_order_2)^2/sum((wy_order_2-mean(wy_order_2))^2)))
se_cc1_5b <- sqrt(sum((res_reg1_2_3 - res_reg2_5b_fit)^2)/((length(wy_order_2) - 2)*sum((wy_order_2-mean(wy_order_2))^2)))
t_cc0_5b <- cc0_5b[k]/se_cc0_5b 
t_cc1_5b <- cc1_5b[k]/se_cc1_5b 
p_cc0_5b[k] = 2*(pt(-abs(t_cc0_5b), df=length(wy_order^2) - 1))
p_cc1_5b[k] = 2*(pt(-abs(t_cc1_5b), df=length(wy_order^2) - 1)) 

# Fit Model 5C
rho_t_res_reg2 <- cor(log(wy_order),res_reg1_2_3)
sd_res_reg1_2_3 <- sd(res_reg1_2_3)
cc0_5c[k] = mean(res_reg1_2_3) - (rho_t_res_reg2 * sd_res_reg1_2_3  * mean(log(wy_order)) / sd(log(wy_order)))
cc1_5c[k] = rho_t_res_reg2 * sd_res_reg1_2_3 / sd(log(wy_order))
res_reg2_5c_fit <- cc0_5c[k] + cc1_5c[k] * log(wy_order)
cond_var_at_n_5c[k] = cc0_5c[k] + cc1_5c[k] * max(log(wy_order))
#plot(wy_order,res_reg1_2_3,col="blue")
#lines(wy_order,res_reg2_5c)
Q_99_at_n_5c[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_5c[k]^3))
rsquared_5c[k] = 1-sum((res_reg1_2_3 - res_reg2_5c_fit)^2)/sum((res_reg1_2_3 - mean(res_reg1_2_3))^2)
se_cc0_5c <- sqrt(sum((res_reg1_2_3 - res_reg2_5c_fit)^2)/(length(wy_order) - 2)*(1/length(wy_order)+mean(log(wy_order))^2/sum((log(wy_order)-mean(log(wy_order)))^2))) 
se_cc1_5c <- sqrt(sum((res_reg1_2_3 - res_reg2_5c_fit)^2)/((length(wy_order) - 2)*sum((log(wy_order)-mean(log(wy_order)))^2)))
t_cc0_5c <- cc0_5c[k]/se_cc0_5c
t_cc1_5c <- cc1_5c[k]/se_cc1_5c
p_cc0_5c[k] = 2*(pt(-abs(t_cc0_5c), df=length(wy_order) - 1))
p_cc1_5c[k] = 2*(pt(-abs(t_cc1_5c), df=length(wy_order) - 1)) 

# Fit Model 6A 

# Establish output arrays for each record (length-dependent)
wt <- array(NA,dim=c(length(wy_order),10))
res_reg1_iwls <- array(NA,dim=c(length(wy_order),10))
res_reg1_iwls_2 <- array(NA,dim=c(length(wy_order),10))
cond_var_6a <- array(NA,dim=c(length(wy_order),10))

for (i in 1:10){
  if(i == 1){
    wt[,i] = rep(1,length(wy_order)) #Make a list
  } else {
    wt[,i] = 1/(cond_var_6a[,i-1]) #Make a list
  }
  flood_reg1_iwls.lm <- lm(log(peak_flow + 0.01) ~ wy_order, weights=wt[,i]) # Make a list
  b0_reg1_iwls[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[1])
  b1_reg1_iwls[k,i] = as.numeric((flood_reg1_iwls.lm)$coefficients[2])
  res_reg1_iwls[,i] = as.numeric(residuals(flood_reg1_iwls.lm)) # Make a list
  res_reg1_iwls_2[,i] = (as.numeric(residuals(flood_reg1_iwls.lm)^2))^(1/3)# Make a list
  res_reg1_iwls_2.lm <- lm(res_reg1_iwls_2[,i] ~ wy_order) # Make a list
  cc0_6a[k,i] = as.numeric((res_reg1_iwls_2.lm$coefficients[1]))
  cc1_6a[k,i] = as.numeric((res_reg1_iwls_2.lm$coefficients[2]))
  cond_var_6a[,i] = as.numeric(fitted(res_reg1_iwls_2.lm)) # Make this a list
  for (l in 1:length(peak_flow))
  if (cond_var_6a[l,i] < 0.001){
    cond_var_6a[l,i] <- 0.001
  } 
 
} 

wt_list[[k]] = wt[,10]
res_reg1_iwls_list[[k]] = res_reg1_iwls[,10]
res_reg1_iwls_2_list[[k]] = res_reg1_iwls_2[,10]
cond_var_6a_list[[k]] = cond_var_6a[,10]

cond_var_at_n_6a[k] = cond_var_6a[length(peak_flow),max(i)]

cond_var_at_n_6a[k] = (cc0_6a[k,10] + cc1_6a[k,10] * max(wy_order))^3 + 
  var(res_reg1_iwls_2[,10]-cond_var_6a[,10])*(3*cc0_6a[k,10] + 3*cc1_6a[k,10] * max(wy_order)) + mean((res_reg1_iwls_2[,10]-cond_var_6a[,10])^3)
#plot(wy_order,res_reg1_2_3,col="blue")
#lines(wy_order,res_reg2_5a)
Q_99_at_n_6a[k] = exp(cond_med_at_n[k] + qnorm(0.99,0,1)*sqrt(cond_var_at_n_6a[k]))
rsquared_6a[k] = 1-sum((res_reg1_iwls[,max(i)] - cond_var_6a[,max(i)])^2)/sum((res_reg1_iwls[,max(i)] - mean(res_reg1_iwls[,max(i)]))^2)
se_cc1_6a <- sqrt(sum((res_reg1_iwls[,max(i)] - cond_var_6a[,max(i)])^2)/((length(wy_order) - 2)*sum((wy_order-mean(wy_order))^2)))
t_cc1_6a <- cc1_6a[k]/se_cc1_6a
p_cc1_6a[k] = 2*(pt(-abs(t_cc1_6a), df=length(wy_order) - 1)) 

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
res_reg1_2 <- as.numeric(res_reg1[[k]])^2
glm.cond_var_gamma <- glm(formula = res_reg1_2 ~ wy_order, family=Gamma(link=log))
b0_glm_cond_var_gamma[k] = as.numeric(glm.cond_var_gamma$coefficients[1])
b1_glm_cond_var_gamma[k] = as.numeric(glm.cond_var_gamma$coefficients[2])
res_dev <- as.numeric(glm.cond_var_gamma$deviance)
null_dev <- as.numeric(glm.cond_var_gamma$null.deviance)
pseudo_rsquared_glm_gamma[k] = 1 - res_dev/null_dev
cond_var_at_n_glm_gamma[k] = b0_glm_cond_var_gamma[k] + b1_glm_cond_var_gamma[k] * max(wy_order) #Transformation bias?
Q_99_at_n_glm_gamma[k] = exp(cond_med_at_n[k]+qnorm(0.99,0,1)*sqrt(cond_var_at_n_glm_gamma[k]))
}

# Model residuals with gamma regression and log link function
#glm.cond_var_a <- glm(formula = res_reg1_2_3 ~ wy_order_1_3, family=Gamma(link=log))
#b0_glm_cond_var_a <- as.numeric(glm.cond_var_a$coefficients[1])
#b1_glm_cond_var_a <- as.numeric(glm.cond_var_a$coefficients[2])
#Q_99_tt_50 <- b0_cond_med + b1_cond_med*nn + qnorm(0.99,0,1)*
#  sqrt(exp(b0_cond_var_glm+b1_cond_var_glm*nn))
# Fit model using GLM 2 

# Fit model using GLM 3 

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
              Q_99_at_n_glm_gamma)
