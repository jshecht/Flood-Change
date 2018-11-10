# Anscombe probability plot 

w_ansc <- rnorm(1000000,0,1)
w_ansc_2_3 <- (w_ansc^2)^(1/3)
w_ansc_2_3_std <- (w_ansc_2_3-mean(w_ansc_2_3))/sd(w_ansc_2_3)

lmoms_w_ansc_2_3_std <- lmoms(w_ansc_2_3_std)

# Estimate parameters from residuals
par_nor <- parnor(lmoms_w_ansc_2_3_std)
par_nor_valid <- are.parnor.valid(par_nor,nowarn=FALSE)

# Compute Weibull probabiities for PPCC test
y_ranked <- sort(w_ansc_2_3_std)
p_y <- rank(y_ranked)/(length(y_ranked)+1)

# Compute PPCC value
F_nor_est <- quanor(p_y,par_nor)
F_nor_QQ_PPCC <- cor(F_nor_est,y_ranked)
plot(F_nor_est,y_ranked)

qqnorm_plot <- qqnorm(w_ansc_2_3_std,xlab="Normal Scores",ylab="Standardized Residuals",main="Anscombe Residual Normality",pch=1,col="gray2",cex=0.5)
qqline(w_ansc_2_3_std,col="black",lwd=0.75)


# Test significance using Heo et al. (2008)
q = 0.10
x = 1.29 + 0.283*log(q) + (0.887-0.751*q + 3.21*q^2)*log(100)
rq = 1 - 1/exp(x)

q = 0.95
x = (1/(-5.22+12.3*q-6.81*q^2))*1000000^(-2.42+5.63*q+-3.11*q^2)
rq = 1 - 1/exp(x)

# Compare with log transformation
w_log <- log(abs(rnorm(1000000,0,1)))
w_log_std <- (w_log - mean(w_log))/sd(w_log)
lmoms_w_log_std <- lmoms(w_log_std)

# Estimate parameters from residuals
par_nor <- parnor(lmoms_w_log_std)
par_nor_valid <- are.parnor.valid(par_nor,nowarn=FALSE)

# Compute Weibull probabiities for PPCC test
y_ranked <- sort(w_log_std)
p_y <- rank(y_ranked)/(length(y_ranked)+1)

# Compute PPCC value
F_nor_est <- quanor(p_y,par_nor)
F_nor_QQ_PPCC <- cor(F_nor_est,y_ranked)
plot(F_nor_est,y_ranked)

qqnorm_plot <- qqnorm(w_log_std,xlab="Normal Scores",ylab="Standardized Residuals",main="Anscombe residual normality",pch=1,col="gray2",cex=0.5)
qqline(w_ansc_2_3_std,col="black",lwd=0.75)

