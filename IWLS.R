# Try IWLS

nn <- 50
tt <- seq(1,50,1)
aa <- 5
bb <- 0.02
ss <- 0.02

bb0_results <- array(NA,dim=c(10000,2))
bb1_results <- array(NA,dim=c(10000,2))
cc0_results <- array(NA,dim=c(10000,2))
cc1_results <- array(NA,dim=c(10000,2))

for (j in 1:10000){
  
  ee <- rnorm(nn,0,1-0.01*ss*tt)
  
  # Compute "true" values of y
  yy <- aa + bb*tt + ee
  
  cv_yy <- sd(yy)/mean(yy)
  
  bb0 <- vector(mode="numeric",length=10)
  bb1 <- vector(mode="numeric",length=10)
  ee_2 <- array(NA,dim=c(50,10))
  cc0 <- vector(mode="numeric",length=10)
  cc1 <- vector(mode="numeric",length=10)
  ee_2_tt_fit <- array(NA,dim=c(50,10))
  wt <- array(NA,dim=c(50,10))
  
  for (i in 1:10){
    if(i == 1){
      wt[,i] = rep(1,50)
    } else {
      wt[,i] = 1/ee_2_tt_fit[,i-1]
    }
    yy_tt.lm <- lm(yy~tt,weights=wt[,i])
    bb0[i] = as.numeric(yy_tt.lm$coefficients[1])
    bb1[i] = as.numeric(yy_tt.lm$coefficients[2])
    ee_2[,i] = residuals(yy_tt.lm)^2
    ee_2_tt.lm <- lm(ee_2[,i] + 1 ~ tt)
    cc0[i] = as.numeric(ee_2_tt.lm$coefficients[1])
    cc1[i] = as.numeric(ee_2_tt.lm$coefficients[2])
    ee_2_tt_fit[,i] = fitted(ee_2_tt.lm)
  }
  bb0_results[j,1] = bb0[1]
  bb1_results[j,1] = bb1[1]
  cc0_results[j,1] = cc0[1]
  cc1_results[j,1] = cc1[1]
  bb0_results[j,2] = bb0[10]
  bb1_results[j,2] = bb1[10]
  cc0_results[j,2] = cc0[10]
  cc1_results[j,2] = cc1[10]
}

