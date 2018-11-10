c0 <- 0.5
c1 <- 0.05
vv <- rnorm(10000,0,0.5)
tt <- runif(10000,1,50)
ee <- vv*(c0+c1*tt)

vv_2_3 <- (vv^2)^(1/3)
ee_2_3 <- (ee^2)^(1/3)

uu <- rnorm(10000,1,0.5)

lm_test <- lm(ee_2_3~tt)
c0_ols <- lm_test$coefficients[1]
c1_ols <- lm_test$coefficients[2]

c1_2 <- (var(ee_2_3)-0.292*mean(ee_2_3))/(1.29*var(tt))

# Anscombe residuals
vv <- rnorm(1000000,0,1)
zz <- rnorm(1000000,0,0.5)

gg <- (vv^2)^(1/3)
hist(gg,xlab="e^(2/3)",main="Anscombe Residuals")

hh <- sign(aa)*((vv^2)^(1/3)*1 + (zz^2))^2^(3/4)
hist(hh,xlab="(e^(2/3)+u)^(3/2)",main="Histogram from Equation 4")
ppcc.test(hh)
hh_std <- (hh-mean(hh))/sd(hh)

# Compute PPCC plot
par_nor <- parnor(hh_std)
par_nor_valid <- are.parnor.valid(par_nor,nowarn=FALSE)

# Compute Weibull probabiities for PPCC test
y_ranked <- sort(hh_std)
p_y <- rank(y_ranked)/(length(y_ranked)+1)

# Compute PPCC value
F_nor_est <- quanor(p_y,par_nor)
F_nor_QQ_PPCC <- cor(F_nor_est,y_ranked)
plot(F_nor_est,y_ranked)

