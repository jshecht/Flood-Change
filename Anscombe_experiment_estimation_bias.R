# Experiment with two-stage least squares fit with Anscombe residuals

# This experiment show

nn <- 76
tt <- runif(1000000,1,nn)
qq_base <- rnorm(1000000,10,1)
beta <- 0.20
qq <- qq_base + beta*tt
rho_qq_tt <- cor(tt,qq)

# Estimate the linear regression model
qq_tt.lm <- lm(qq~tt)
summary(qq_tt.lm)
fit_qq_tt <- fitted(qq_tt.lm)
b0_qq_tt <- summary(qq_tt.lm)$coefficients[1,1]
b1_qq_tt <- summary(qq_tt.lm)$coefficients[2,1] 
qq_tt_final <- b0_qq_tt + b1_qq_tt*nn


# Estimate 100-year flood with nonstationary quantile function
Q_99_mean_only <- qq_tt_final + qnorm(0.99,0,1)*sqrt(var(qq)-b1_qq_tt^2*var(tt))
print(Q_99_mean_only)

# Residual analysis 
res_qq_tt <- residuals(qq_tt.lm)
res_2_3_qq_tt <- (res_qq_tt^2)^(1/3)

# Try fit with tt_1_3
tt_1_3 <- tt^(1/3)
tt_1_3_res_2_3_qq_tt <- lm(res_2_3_qq_tt~tt_1_3)
fit_tt_1_3_res_2_3_qq_tt <- fitted(tt_1_3_res_2_3_qq_tt)
Q_99_tt_1_3 <- qq_tt_final + qnorm(0.99,0,1)*fit_tt_1_3_res_2_3_qq_tt^(3/2)
mean(Q_99_tt_1_3) 

# Try fit with tt_2_3
tt_2_3 <- (tt^2)^(1/3)
tt_2_3_res_2_3_qq_tt <- lm(res_2_3_qq_tt~tt_2_3)
fit_tt_2_3_res_2_3_qq_tt <- fitted(tt_2_3_res_2_3_qq_tt)
Q_99_tt_2_3 <- qq_tt_final + qnorm(0.99,0,1)*fit_tt_2_3_res_2_3_qq_tt^(3/2)
mean(Q_99_tt_2_3)
# 19.27346

# Try fit with absolute residuals
res_1_qq_tt <- (res_qq_tt^2)^(1/2)
tt_1_2 <- tt^(1/2)
tt_1_2_res_1_qq_tt <- lm(res_1_qq_tt~tt_1_2)
fit_tt_1_2_res_1_qq_tt <- fitted(tt_1_2_res_1_qq_tt)
Q_99_tt_1_2_res_1 <- qq_tt_final + qnorm(0.99,0,1)*fit_tt_1_2_res_1_qq_tt
mean(Q_99_tt_1_2_res_1) 
# 19.4578

# Try fit with square root residuals
res_1_2_qq_tt <- (res_qq_tt^2)^(1/4)
tt_1_4 <- tt^(1/4)
tt_1_4_res_1_2_qq_tt <- lm(res_1_2_qq_tt~tt_1_4)
fit_tt_1_4_res_1_2_qq_tt <- fitted(tt_1_4_res_1_2_qq_tt)
Q_99_tt_1_4_res_1_2 <- qq_tt_final + qnorm(0.99,0,1)*fit_tt_1_4_res_1_2_qq_tt^2
mean(Q_99_tt_1_4_res_1_2)
# 19.17384

# Try fit with square residuals
res_2_qq_tt <- res_qq_tt^2
tt_1_res_2_qq_tt <- lm(res_2_qq_tt~tt)
fit_tt_1_res_2_qq_tt <- fitted(tt_1_res_2_qq_tt)
Q_99_tt_1_res_2 <- qq_tt_final + qnorm(0.99,0,1)*fit_tt_1_res_2_qq_tt^(1/2)
mean(Q_99_tt_1_res_2) #Unbiased! 

# Compute fit for different powers k assuming tt^(k/2) ~ res_qq_tt^k
nn <- 50
tt <- runif(1000000,1,nn)
alphas <- 10
s_base <- 2
b_trend <- seq(0,0.05,0.005)
s_trend <- seq(0,0.02,0.002)
pwr <- seq(0.25,2,0.25)
Q_99_mean <-array(NA,dim=c(length(b_trend),length(s_trend),length(pwr)))
pbias <-array(NA,dim=c(length(b_trend),length(s_trend),length(pwr)))

for (b in 1:length(b_trend)){
  qq_trend = alphas + b_trend[b]*tt
  for (s in 1:length(s_trend)){
    ss <- s_base + s_trend[s]*sqrt(tt)
    qq <- qq_trend + rnorm(1000000,0,ss)
    qq_tt.lm <- lm(qq~tt)
    summary(qq_tt.lm)
    fit_qq_tt <- fitted(qq_tt.lm)
    b0_qq_tt <- summary(qq_tt.lm)$coefficients[1,1]
    b1_qq_tt <- summary(qq_tt.lm)$coefficients[2,1] #Assumes beta not affected by heterosc
    qq_tt_final <- b0_qq_tt + b1_qq_tt*nn
    res_qq_tt <- residuals(qq_tt.lm)
    for (k in 1:length(pwr)){
      res_pwr_qq_tt <- (res_qq_tt^2)^(pwr[k]/2)
      tt_pwr <- tt^(pwr[k]/2)
      tt_res_qq_tt_pwr <- lm(res_pwr_qq_tt ~ tt_pwr)
      fit_tt_res_qq_tt_pwr <- fitted(tt_res_qq_tt_pwr)
      Q_99_tt_res_pwr <- qq_tt_final + qnorm(0.99,0,1)*fit_tt_res_qq_tt_pwr^(1/pwr[k])
      Q_99_mean[b,s,k] = mean(Q_99_tt_res_pwr)
    }
  }
}

for (l in 1:length(pwr)){
  pbias[,,l] = (Q_99_mean[,,l]-Q_99_mean[,,4])/Q_99_mean[,,4]
}

# Convert to single column
dim(pbias) <- c(prod(dim(pbias)), 1) 

preds <- expand.grid(b_trend,s_trend,pwr)
pred1_btrend <- preds[,1]
pred2_strend <- preds[,2]
pred3_pwr <- preds[,3]
pwr_trnsf_reg_mod1 <- lm(pbias ~ pred1_btrend + pred2_strend + pred3_pwr)
summary(pwr_trnsf_reg_mod1)

---------------------------------------------------------------------------
  


#Boxplot comparison

Q_99_means <- c(mean(Q_99_tt_1_3),mean(Q_99_tt_2_3),mean(Q_99_tt_1_3_c),mean(Q_99_tt_2_3_c))
pbias_Q_99 <- (Q_99_means-Q_99_mean_only)/Q_99_mean_only
barplot(pbias_Q_99*100,names=c("2-stage \n OLS t^(1/3)","2-stage \n OLS t^(2/3)","Var Model 1","Var Model 2"),ylim=c(-4,6),ylab="Bias(%)",cex.names=0.75)
abline(h=0)



