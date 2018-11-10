# Test new residual derivations

tt <- runif(10000,1,50)
ee <- rnorm(10000,0,1)*(0.5+0.01*tt)
plot(tt,ee)
plot(tt,ee^2)
cc0 <- mean(ee^2) - cor(tt,ee^2)*sd(ee^2)*mean(tt)/sd(tt)
cc1 <- cor(tt,ee^2)*sd(ee^2)/sd(tt)

ee_2_3_fit <- 0.802*((cc0+cc1*tt)^2)^(1/3)
ww_2_3 <- (ee^2)^(1/3) - 0.802*(cc0+cc1*tt)^(1/3)

ee_2_est <- 0.516*(cc0+cc1*tt)+1.929*((cc0+cc1*tt)^2)^(1/3)*(ww_2_3)+2.406*(cc0+cc1*tt)^(1/3)*(ww_2_3)^2+(ww_2_3)^3

plot(ee_2_est - ee^2)