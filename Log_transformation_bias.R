# Transformation bias check

x <- runif(1000000,1,10)
y <- rlnorm(1000000,5,1)
quantile(y,0.50)  #148.5744
quantile(y,0.99)  #1523.723

log_y <- log(y)
lm_test <- lm(log_y~x)
log_y_pred <- lm_test$coefficients[1] + lm_test$coefficients[2]*x
mean(exp(log_y_pred))
res <- log_y-log_y_pred

mean(exp(log_y_pred)) #148.4663
mean(exp(log_y_pred+0.5*var(res))) #244.8382

exp(mean(log_y_pred)+qnorm(0.99,0,1)*sd(res)) #1521.213

exp(mean(log_y_pred)+0.5*var(res)+qnorm(0.99,0,1)*sd(res)) #2508.656

