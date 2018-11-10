# Anscombe residuals of perfectly normal distributions

# Generate normal distributions

# Normality

# Create vectors for storing results
skew_x_2_3 <- array(NA,dim=100)
pval_ppcc_test_x <- array(NA,dim=100)
pval_ppcc_test_x_2_3 <- array(NA,dim=100)
pval_s_w_test_x_2_3 <- array(NA,dim=100)
pval_ppcc_t_x_2_3 <- array(NA,dim=100)
pval_ppcc_ansc_gamma <- array(NA,dim=100)
pval_ppcc_ansc_gamma2 <- array(NA,dim=100)

for (i in 1:100){
  x <- rnorm(1000000,0,1) 
  ppcc_test_x = ppcc.test(x)
  pval_ppcc_test_x[i] <- as.numeric(ppcc_test_x[2])
  x_2_3 <- (x^2)^(1/3)
  #skew_x_2_3[i] = skewness(x_2_3)
  ppcc_test_x_2_3 = ppcc.test(x_2_3)
  pval_ppcc_test_x_2_3[i] = as.numeric(ppcc_test_x_2_3[2])
  s_w_test_x_2_3 = shapiro.test(x_2_3)
  pval_s_w_test_x_2_3[i] = as.numeric(s_w_test_x_2_3[2])
  
  t_x_2_3.lm <- lm(x_2_3~t)
  res_t_x_2_3 <- residuals(t_x_2_3.lm)
  ppcc_test_t_x_2_3 <- ppcc.test(res_t_x_2_3)
  pval_ppcc_t_x_2_3[i] = as.numeric(ppcc_test_t_x_2_3[2])
  
  x_2 <- x^2
  t_x_2.lm <- lm(x_2~t)
  res_t_x_2 <- residuals(t_x_2.lm)
  fit_t_x_2 <- fitted(t_x_2.lm)
  ansc_gamma <- 3*(x_2^(1/3)-fit_t_x_2^(1/3))/fit_t_x_2^(1/3)
  ppcc_test_ansc_gamma <- ppcc.test(ansc_gamma)
  pval_ppcc_ansc_gamma[i] = as.numeric(ppcc_test_ansc_gamma[2])
  
  t_x.lm <- lm(x~t)
  res_t_x <- residuals(t_x.lm)
  fit_t_x <- fitted(t_x.lm)
  ansc_gamma2 <- 3*((x^2)^(1/3)-(fit_t_x^2)^(1/3))/(fit_t_x^2)^(1/3)
  ppcc_test_ansc_gamma <- ppcc.test(ansc_gamma)
  pval_ppcc_ansc_gamma2[i] = as.numeric(ppcc_test_ansc_gamma[2])
 
}

boxplot(skew_x_2_3,main="Sample skewness coefficient")
boxplot(pval_ppcc_test_x,names=c("Sample from normal dist (N=76)"),ylab="PPCC value")
boxplot(pval_s_w_test_x_2_3)
boxplot(pval_ppcc_t_x_2_3)
boxplot(pval_ppcc_test_x_2_3,pval_ppcc_ansc_gamma2,names=c(expression(~epsilon^(2/3)),"SAS formula"),ylab="PPCC p-value")

boxplot(pval_ppcc_test_x_2_3,names=c("r^(2/3)"),main=expression("PPCC test for" ~ epsilon^(2/3)),ylab="Normal PPCC p-value")

hist(x_2_3)
hist(ansc_gamma2)





