# MONTE CARLO EXPERIMENTS

# Basic data
nn <- 50
tt <- runif(1000000,1,nn)
aa <- 10
bb <- 0.05
ss <- 0.01
qq <- aa + bb*tt + rnorm(1000000,0,1+ss*tt)
Q_99_tt_50 <- aa + bb*nn + qnorm(0.99,0,1)*(1+ss*nn) #15.98952

# Fit conditional model of the mean
lm.cond_mean <- lm(qq~tt)
summary(lm.cond_mean)
b0_cond_mean <- as.numeric(lm.cond_mean$coefficients[1])
b1_cond_mean <- as.numeric(lm.cond_mean$coefficients[2])
ee <- residuals(lm.cond_mean)
ee_2 <- ee^2

# Model residuals with gamma regression and log link function
glm.cond_var <- glm(formula = ee^2 ~ tt, family=Gamma(link=log))
summary(glm.cond_var)
b0_cond_var_glm <- as.numeric(glm.cond_var$coefficients[1])
b1_cond_var_glm <- as.numeric(glm.cond_var$coefficients[2])
Q_99_tt_50 <- b0_cond_mean + b1_cond_mean*nn + qnorm(0.99,0,1)*
  sqrt(exp(b0_cond_var_glm+b1_cond_var_glm*nn))
# Slight overestimate 16.02943

# Model absolute residuals with gamma regression and no link function
glm.cond_var <- glm(formula = abs(ee) ~ tt, family=Gamma)
b0_cond_var_glm <- as.numeric(glm.cond_var$coefficients[1])
b1_cond_var_glm <- as.numeric(glm.cond_var$coefficients[2])
Q_99_tt_50 <- b0_cond_mean + b1_cond_mean*nn + qnorm(0.99,0,1)*
  (b0_cond_var_glm+b1_cond_var_glm*nn)
# Underestimate 14.39404

# Model residuals with gamma regression 
glm.cond_var <- glm(formula = ee^2 ~ tt, family=Gamma)
b0_cond_var_glm <- as.numeric(glm.cond_var$coefficients[1])
b1_cond_var_glm <- as.numeric(glm.cond_var$coefficients[2])
Q_99_tt_50 <- b0_cond_mean + b1_cond_mean*nn + qnorm(0.99,0,1)*
  sqrt(b0_cond_var_glm+b1_cond_var_glm*nn)
# Underestimate 13.99359

# Model residuals with approximate Anscombe transformation (2/3 power)
ee_2_3 <- (ee^2)^(1/3)
tt_1_3 <- tt^(1/3)
lm_cond_var_ansc <- lm(ee_2_3 ~ tt_1_3) 
b0_cond_var_ansc <- as.numeric(lm_cond_var_ansc$coefficients[1])
b1_cond_var_ansc <- as.numeric(lm_cond_var_ansc$coefficients[2])
Q_99_tt_50 <- b0_cond_mean + b1_cond_mean*nn + qnorm(0.99,0,1)*
  (b0_cond_var_ansc+b1_cond_var_ansc*nn^(1/3))^(3/2)
# Underestimate 14.90895
# Fix with bias correction factor

ee_2_3 <- (ee^2)^(1/3)
tt_2_3 <- tt^(2/3)
lm_cond_var_ansc <- lm(ee_2_3 ~ tt_2_3) 
b0_cond_var_ansc <- as.numeric(lm_cond_var_ansc$coefficients[1])
b1_cond_var_ansc <- as.numeric(lm_cond_var_ansc$coefficients[2])
Q_99_tt_50 <- b0_cond_mean + b1_cond_mean*nn + qnorm(0.99,0,1)*
  (b0_cond_var_ansc+b1_cond_var_ansc*nn^(1/3))^(3/2)
# Underestimate 14.27874
# Fix with bias correction factor


# Model squared residuals 
lm_cond_var_sq <- lm(ee_2 ~ tt)
b0_cond_var_sq <- as.numeric(lm_cond_var_sq$coefficients[1])
b1_cond_var_sq <- as.numeric(lm_cond_var_sq$coefficients[2])
Q_99_tt_50 <- b0_cond_mean + b1_cond_mean*nn + qnorm(0.99,0,1)*
  sqrt(b0_cond_var_sq+b1_cond_var_sq*nn)
# Slight underestimate 15.95838


# Use glm.fit to do iteratively reweighted least squares


