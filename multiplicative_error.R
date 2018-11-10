#New model

# Generate fake data
ss_ee <- 0.25
ee <- rnorm(10000,1,ss_ee*(1+0.01*xx)) #Possible problem with ee < 0
xx <- runif(10000,1,50)
yy <- (5 + 0.1*xx)*ee
plot(xx,yy)

# Fit a model to the fake data
lm_reg1 <- lm(yy ~ xx)
aa <- as.numeric(lm_reg1$coefficients[1])
bb <- as.numeric(lm_reg1$coefficients[2])
lines(xx,aa+bb*xx,col="red")
summary(lm_reg1)
vv <- residuals(lm_reg1)
#rse <- summary(lm_test)$sigma
ee_est <- (vv + (aa+bb*xx))/(aa+bb*xx) # ee_est have a constant scaling factor
plot(xx,vv)
plot(xx,ee_est)
yy_est <- (aa+bb*xx)*ee_est
plot(xx,yy_est)

ee_est_2_3 <- (ee_est^2)^(1/3)
plot(xx,ee_est_2_3)
lm_reg2 <- lm(ee_est_2_3 ~ xx)
cc0 <- as.numeric(lm_reg2$coefficients[1])
cc1 <- as.numeric(lm_reg2$coefficients[2]) # ~0 because there is no trend in the multiplicative error! 
lines(xx,cc0+cc1*xx,col="green")
ww <- residuals(lm_reg2)
uu_est <- (ww + cc0 + cc1*xx)/(cc0 + cc1*xx)
ee_est_2_3_est <- (cc0 + cc1*xx)*uu_est
plot(xx,ww) #Homosceadstic residuals

Z <- ((ee_est - 1)^2)^(1/3)
lm_reg3 <- lm(Z ~ xx)
ww2 <- residuals(lm_reg3)



#vv_est_2_3 <- (vv^2)^(1/3)
#plot(xx,vv_est_2_3)
#lm_reg2a <- lm(vv_est_2_3 ~ xx)
#dd0 <- as.numeric(lm_reg2a$coefficients[1])
#dd1 <- as.numeric(lm_reg2a$coefficients[2])
#lines(xx,dd0 + dd1*xx,col="blue")
#ww2 <- residuals(lm_reg2a)
#var_ww2 <- (dd0 + dd1*xx)^3 + 3*var(ww2)*(dd0+dd1*xx)


# Estimate 99th percentile
Q99 <- (aa + bb*50) + qnorm(0.99,0,1)*sqrt(max(var_ww2))
