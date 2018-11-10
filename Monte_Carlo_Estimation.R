# New Monte Carlo experiments

# Number of runs
mm <- 100

# Record length
nn <- 100

# Select range of coefficients to test
bb0_true <- seq(2,4,2)
bb1_true <- seq(-0.01,0.01,0.01)
cc0_true <- seq(0.5,2,0.5)
cc1_true <- seq(-0.0025,0.0025,0.0025)

# Add codebreak for infeasible combinations

# Create storage arrays
bb0_est <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
bb1_est <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
cc0_est_ee_2 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
cc1_est_ee_2 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
cc0_est_ee_2_3 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
cc1_est_ee_2_3 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
cc0_est_ee_2_3_MoM <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
cc1_est_ee_2_3_MoM <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
cc0_est_glm_ee_2 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
cc1_est_glm_ee_2 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
var_max_tt_ee_2 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
var_max_tt_ee_2_3 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
var_max_tt_ee_2_3_MoM <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
var_max_tt_glm_ee_2 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))

# For IWLS runs

# Anscombe res
wt <- array(NA,dim=c(nn,10))
bb0_iwls_ee_2_3w <- vector(length=10)
bb1_iwls_ee_2_3w <- vector(length=10)
cc0_iwls_ee_2_3w<- vector(length=10)
cc1_iwls_ee_2_3w <- vector(length=10)
cvar_iwls_ee_2_3w <- array(NA,dim=c(nn,10))

bb0_iwls_ee_2_3 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
bb1_iwls_ee_2_3 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
cc0_iwls_ee_2_3 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
cc1_iwls_ee_2_3 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
cvar_iwls_ee_2_3 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
var_max_tt_iwls_ee_2_3 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))


# GLM
wt <- array(NA,dim=c(nn,10))
bb0_iwls_glm_ee_2w <- vector(length=10)
bb1_iwls_glm_ee_2w <- vector(length=10)
cc0_iwls_glm_ee_2w<- vector(length=10)
cc1_iwls_glm_ee_2w <- vector(length=10)
cvar_iwls_glm_ee_2w <- array(NA,dim=c(nn,10))

bb0_iwls_glm_ee_2 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
bb1_iwls_glm_ee_2 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
cc0_iwls_glm_ee_2 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
cc1_iwls_glm_ee_2 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
cvar_iwls_glm_ee_2 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
var_max_tt_iwls_glm_ee_2 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))

# For quantile estimates
yy_100_at_tt_max_TRUE <- array(NA,dim=c(length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
yy_100_at_tt_max_ee_2 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
yy_100_at_tt_max_ee_2_3 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
yy_100_at_tt_max_ee_2_3_MoM <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
yy_100_at_tt_max_iwls_ee_2_3 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
yy_100_at_tt_max_glm_ee_2 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))
yy_100_at_tt_max_iwls_glm_ee_2 <- array(NA,dim=c(length=mm,length(bb0_true),length(bb1_true),length(cc0_true),length(cc1_true)))


  # Set values
for (j in 1:length(bb0_true)){  
  for (k in 1:length(bb1_true)){
    for (l in 1:length(cc0_true)){
      for (m in 1:length(cc1_true)){
          for (i in 1:mm){
  
          #bb0 <- 5
          #bb1 <- 0.01
  
          tt_max <- 50
          tt <- runif(nn,1,tt_max)
          
          #cc0_true <- 1
          #cc1_true <- -0.01
  
          # Generate normally distributed errors
          error <- rnorm(nn,0,cc0_true[l] + cc1_true[m]*tt)
  
          # Computes yy
          yy <- bb0_true[j] + bb1_true[k]*tt + error
  
          # Compute TRUE VALUE of 100-year flood
          yy_100_at_tt_max_TRUE[j,k,l,m] = bb0_true[j] + bb1_true[k]*tt_max + qnorm(0.99,0,cc0_true[l]+cc1_true[m]*tt_max)
  
          # Estimate 100-year flood conditional mean
          lm.yy_tt <- lm(yy ~ tt)
          bb0_est[i,j,k,l,m] = as.numeric(lm.yy_tt$coefficients[1])
          bb1_est[i,j,k,l,m] = as.numeric(lm.yy_tt$coefficients[2])
          
          # Estimate 100-year flood using squared residuals
          ee <- as.numeric(residuals(lm.yy_tt))
          ee_2 <- ee^2
          lm.ee2_tt <- lm(ee_2 ~ tt)
          cc0_est_ee_2[i,j,k,l,m] = lm.ee2_tt$coefficients[1]
          cc1_est_ee_2[i,j,k,l,m] = lm.ee2_tt$coefficients[2]
          var_max_tt_ee_2[i,j,k,l,m] = cc0_est_ee_2[i,j,k,l,m] + cc1_est_ee_2[i,j,k,l,m] * tt_max
          yy_100_at_tt_max_ee_2[i,j,k,l,m] = (bb0_est[i,j,k,l,m] + bb1_est[i,j,k,l,m]*tt_max) + qnorm(0.99,0,1)*sqrt(var_max_tt_ee_2[i,j,k,l,m])
          
          # Estimate 100-year flood using Anscombe residuals
          ee_2_3 <- ee_2^(1/3)
          lm.ee2_3_tt <- lm(ee_2_3 ~ tt)
          cc0_est_ee_2_3[i,j,k,l,m] = lm.ee2_3_tt$coefficients[1]
          cc1_est_ee_2_3[i,j,k,l,m] = lm.ee2_3_tt$coefficients[2]
          ww <- residuals(lm.ee2_3_tt)
          var_max_tt_ee_2_3[i,j,k,l,m] = (cc0_est_ee_2_3[i,j,k,l,m] + cc1_est_ee_2_3[i,j,k,l,m]*tt_max)^3 + 3*var(ww)*(cc0_est_ee_2_3[i,j,k,l,m] + cc1_est_ee_2_3[i,j,k,l,m]*tt_max) + mean(ww^3)
          yy_100_at_tt_max_ee_2_3[i,j,k,l,m] = (bb0_est[i,j,k,l,m] + bb1_est[i,j,k,l,m]*tt_max) + qnorm(0.99,0,1)*sqrt(var_max_tt_ee_2_3[i,j,k,l,m])
          
          # Estimate 100-year flood using Anscombe residuals and method of moments
          ee_2_3 <- ee_2^(1/3)
          cc0_est_ee_2_3_MoM[i,j,k,l,m] = mean(ee_2_3) - cor(tt,ee_2_3)*sd(ee_2_3)*mean(tt)/sd(tt)
          cc1_est_ee_2_3_MoM[i,j,k,l,m] = cor(tt,ee_2_3)*sd(ee_2_3)/sd(tt)
          ww <- ee_2_3 - (cc0_est_ee_2_3_MoM[i,j,k,l,m] + cc1_est_ee_2_3_MoM[i,j,k,l,m] * tt)
          var_max_tt_ee_2_3_MoM[i,j,k,l,m] = (cc0_est_ee_2_3_MoM[i,j,k,l,m] + cc1_est_ee_2_3_MoM[i,j,k,l,m]*tt_max)^3 + 3*var(ww)*(cc0_est_ee_2_3_MoM[i,j,k,l,m] + cc1_est_ee_2_3_MoM[i,j,k,l,m]*tt_max) + mean(ww^3)
          yy_100_at_tt_max_ee_2_3_MoM[i,j,k,l,m] = (bb0_est[i,j,k,l,m] + bb1_est[i,j,k,l,m]*tt_max) + qnorm(0.99,0,1)*sqrt(var_max_tt_ee_2_3_MoM[i,j,k,l,m])
          
          # Estimate 100-year flood using Anscombe residuals and IWLS
          for (w in 1:10){
            if(w == 1){
              wt[,w] = rep(1,nn) #Make a list
            } else {
              wt[,w] = 1/(cvar_iwls_ee_2_3w[,w-1]) 
            }
            reg1_iwls_ee_2_3.lm <- lm(yy ~ tt, weights=wt[,w]) 
            bb0_iwls_ee_2_3w[w] = as.numeric((reg1_iwls_ee_2_3.lm)$coefficients[1])
            bb1_iwls_ee_2_3w[w] = as.numeric((reg1_iwls_ee_2_3.lm)$coefficients[2])
            iwls_ee <- as.numeric(residuals(reg1_iwls_ee_2_3.lm)) 
            iwls_ee_2_3 <- (iwls_ee^2)^(1/3)
            reg2_iwls_ee_2_3.lm <- lm(iwls_ee_2_3 ~ tt) 
            cc0_iwls_ee_2_3w[w] = as.numeric((reg2_iwls_ee_2_3.lm$coefficients[1])) 
            cc1_iwls_ee_2_3w[w] = as.numeric((reg2_iwls_ee_2_3.lm$coefficients[2])) 
            cvar_iwls_ee_2_3w[,w] = as.numeric(fitted(reg2_iwls_ee_2_3.lm))
          }  
          
          bb0_iwls_ee_2_3[i,j,k,l,m] = bb0_iwls_ee_2_3w[10] 
          bb1_iwls_ee_2_3[i,j,k,l,m] = bb1_iwls_ee_2_3w[10]
          cc0_iwls_ee_2_3[i,j,k,l,m] = cc0_iwls_ee_2_3w[10]
          cc1_iwls_ee_2_3[i,j,k,l,m] = cc1_iwls_ee_2_3w[10]
          cvar_iwls_ee_2_3[i,j,k,l,m] = cvar_iwls_ee_2_3w[nn,10]
          
          var_max_tt_iwls_ee_2_3[i,j,k,l,m] = (cc0_iwls_ee_2_3[i,j,k,l,m] + cc1_iwls_ee_2_3[i,j,k,l,m]*tt_max)^3 + 3*var(iwls_ee_2_3)*(cc0_iwls_ee_2_3[i,j,k,l,m] + cc1_iwls_ee_2_3[i,j,k,l,m]*tt_max) + mean((cc0_iwls_ee_2_3[i,j,k,l,m] + cc1_iwls_ee_2_3[i,j,k,l,m]*tt_max)^3)
          yy_100_at_tt_max_iwls_ee_2_3[i,j,k,l,m] = (bb0_iwls_ee_2_3[i,j,k,l,m] + bb1_iwls_ee_2_3[i,j,k,l,m]*tt_max) + qnorm(0.99,0,1)*sqrt(var_max_tt_iwls_ee_2_3[i,j,k,l,m])
          
          
          # Estimate 100-year flood using a GLM to estimate the conditional variance
          glm_ee_2 <- glm2(ee_2 ~ tt,family=Gamma(link=log)) 
          cc0_est_glm_ee_2[i,j,k,l,m] = as.numeric(glm_ee_2$coefficients[1])
          cc1_est_glm_ee_2[i,j,k,l,m] = as.numeric(glm_ee_2$coefficients[2])
          var_max_tt_glm_ee_2[i,j,k,l,m] = exp(cc0_est_glm_ee_2[i,j,k,l,m]*cc1_est_glm_ee_2[i,j,k,l,m]*tt_max)
          yy_100_at_tt_max_glm_ee_2[i,j,k,l,m] = (bb0_est[i,j,k,l,m] + bb1_est[i,j,k,l,m]*tt_max) + qnorm(0.99,0,1)*sqrt(var_max_tt_glm_ee_2[i,j,k,l,m])
  
          
          # Estimate using GLM-IWLS
          for (w in 1:10){
            if(w == 1){
              wt[,w] = rep(1,nn) #Make a list
            } else {
              wt[,w] = 1/(cvar_iwls_glm_ee_2w[,w-1]) 
            }
            reg1_iwls_glm_ee_2.lm <- lm(yy ~ tt, weights=wt[,w]) 
            bb0_iwls_glm_ee_2w[w] = as.numeric((reg1_iwls_glm_ee_2.lm)$coefficients[1])
            bb1_iwls_glm_ee_2w[w] = as.numeric((reg1_iwls_glm_ee_2.lm)$coefficients[2])
            iwls_glm_ee <- as.numeric(residuals(reg1_iwls_glm_ee_2.lm)) 
            iwls_glm_ee_2 <- iwls_glm_ee^2
            reg2_iwls_glm_ee_2.lm <- glm2(iwls_glm_ee_2 ~ tt,family=Gamma(link=log)) 
            cc0_iwls_glm_ee_2w[w] = as.numeric((reg2_iwls_glm_ee_2.lm$coefficients[1])) 
            cc1_iwls_glm_ee_2w[w] = as.numeric((reg2_iwls_glm_ee_2.lm$coefficients[2])) 
            cvar_iwls_glm_ee_2w[,w] = as.numeric(fitted(reg2_iwls_glm_ee_2.lm))
          }  
           
          bb0_iwls_glm_ee_2[i,j,k,l,m] = bb0_iwls_glm_ee_2w[10] 
          bb1_iwls_glm_ee_2[i,j,k,l,m] = bb1_iwls_glm_ee_2w[10]
          cc0_iwls_glm_ee_2[i,j,k,l,m] = cc0_iwls_glm_ee_2w[10]
          cc1_iwls_glm_ee_2[i,j,k,l,m] = cc1_iwls_glm_ee_2w[10]
          cvar_iwls_glm_ee_2[i,j,k,l,m] = cvar_iwls_glm_ee_2w[nn,10] 
          
          var_max_tt_iwls_glm_ee_2[i,j,k,l,m] = exp(cc0_iwls_glm_ee_2[i,j,k,l,m] + cc1_iwls_glm_ee_2[i,j,k,l,m]*tt_max)
          yy_100_at_tt_max_iwls_glm_ee_2[i,j,k,l,m] = (bb0_iwls_glm_ee_2[i,j,k,l,m] + bb1_iwls_glm_ee_2[i,j,k,l,m]*tt_max) + qnorm(0.99,0,1)*sqrt(var_max_tt_iwls_glm_ee_2[i,j,k,l,m])
          

        }   
      }
    }
  }
}

# ADD DEGREES OF FREEDOM CORRECTION

# Compare quantile estimates
mean_yy_100_at_tt_max_ee_2 <- apply(yy_100_at_tt_max_ee_2,c(2,3,4,5),mean)
mean_yy_100_at_tt_max_ee_2_3 <- apply(yy_100_at_tt_max_ee_2_3,c(2,3,4,5),mean)
mean_yy_100_at_tt_max_ee_2_3_MoM <- apply(yy_100_at_tt_max_ee_2_3_MoM,c(2,3,4,5),mean)
mean_yy_100_at_tt_max_iwls_ee_2_3 <- apply(yy_100_at_tt_max_iwls_ee_2_3,c(2,3,4,5),mean)
mean_yy_100_at_tt_max_glm_ee_2 <- apply(yy_100_at_tt_max_glm_ee_2,c(2,3,4,5),mean)
mean_yy_100_at_tt_max_iwls_glm_ee_2 <- apply(yy_100_at_tt_max_iwls_glm_ee_2,c(2,3,4,5),mean)

# Compute relative bias of quantile estimates (%)
rbias_ee_2 <- (mean_yy_100_at_tt_max_ee_2 - yy_100_at_tt_max_TRUE)/yy_100_at_tt_max_TRUE
rbias_ee_2_3 <- (mean_yy_100_at_tt_max_ee_2_3 - yy_100_at_tt_max_TRUE)/yy_100_at_tt_max_TRUE
rbias_ee_2_3_MoM <- (mean_yy_100_at_tt_max_ee_2_3_MoM - yy_100_at_tt_max_TRUE)/yy_100_at_tt_max_TRUE
rbias_iwls_ee_2_3 <- (mean_yy_100_at_tt_max_iwls_ee_2_3 - yy_100_at_tt_max_TRUE)/yy_100_at_tt_max_TRUE
rbias_glm_ee_2 <- (mean_yy_100_at_tt_max_glm_ee_2 - yy_100_at_tt_max_TRUE)/yy_100_at_tt_max_TRUE
rbias_iwls_glm_ee_2 <- (mean_yy_100_at_tt_max_iwls_glm_ee_2 - yy_100_at_tt_max_TRUE)/yy_100_at_tt_max_TRUE
boxplot(rbias_ee_2,rbias_ee_2_3,rbias_ee_2_3_MoM,rbias_iwls_ee_2_3,rbias_glm_ee_2,rbias_iwls_glm_ee_2)

# Compute relative RMSE of quantile estimates (%)
sq_ee_2 <- apply(yy_100_at_tt_max_ee_2,1,function(yy_100_at_tt_max_ee_2) ((yy_100_at_tt_max_ee_2-yy_100_at_tt_max_TRUE)/yy_100_at_tt_max_TRUE)^2)
sq_ee_2_3 <- apply(yy_100_at_tt_max_ee_2_3,1,function(yy_100_at_tt_max_ee_2_3) ((yy_100_at_tt_max_ee_2_3-yy_100_at_tt_max_TRUE)/yy_100_at_tt_max_TRUE)^2)
sq_ee_2_3_MoM <- apply(yy_100_at_tt_max_ee_2_3_MoM,1,function(yy_100_at_tt_max_ee_2_3_MoM) ((yy_100_at_tt_max_ee_2_3_MoM-yy_100_at_tt_max_TRUE)/yy_100_at_tt_max_TRUE)^2)
sq_iwls_ee_2_3 <- apply(yy_100_at_tt_max_iwls_ee_2_3,1,function(yy_100_at_tt_max_iwls_ee_2_3) ((yy_100_at_tt_max_iwls_ee_2_3-yy_100_at_tt_max_TRUE)/yy_100_at_tt_max_TRUE)^2)
sq_glm_ee_2 <- apply(yy_100_at_tt_max_glm_ee_2,1,function(yy_100_at_tt_max_glm_ee_2) ((yy_100_at_tt_max_glm_ee_2-yy_100_at_tt_max_TRUE)/yy_100_at_tt_max_TRUE)^2)
sq_iwls_glm_ee_2 <- apply(yy_100_at_tt_max_iwls_glm_ee_2,1,function(yy_100_at_tt_max_iwls_glm_ee_2) ((yy_100_at_tt_max_iwls_glm_ee_2-yy_100_at_tt_max_TRUE)/yy_100_at_tt_max_TRUE)^2)

rrmse_ee_2 <- sqrt(1/mm*apply(sq_ee_2,1,sum))
rrmse_ee_2_3 <- sqrt(1/mm*apply(sq_ee_2_3,1,sum))
rrmse_ee_2_3_MoM <- sqrt(1/mm*apply(sq_ee_2_3_MoM,1,sum))
rrmse_iwls_ee_2_3 <- sqrt(1/mm*apply(sq_iwls_ee_2_3,1,sum))
rrmse_glm_ee_2 <- sqrt(1/mm*apply(sq_glm_ee_2,1,sum))
rrmse_iwls_glm_ee_2 <- sqrt(1/mm*apply(sq_iwls_glm_ee_2,1,sum))

boxplot(rrmse_ee_2,rrmse_ee_2_3,rrmse_ee_2_3_MoM,rrmse_iwls_ee_2_3,rrmse_glm_ee_2,rrmse_iwls_glm_ee_2,
        names=c("ee_2","ee_2_3","ee_2_3_MoM","iwls_ee_2_3","glm_ee_2","iwls_glm_ee_2"),ylab="Relative RMSE",cex.axis=0.8)


# Compute relative bias of regression coefficients
rbias_bb0 <- (bb0_est - bb0_true)/bb0_true
rbias_bb0_iwls_ee_2_3 <- (bb0_iwls_ee_2_3 - bb0_true)/bb0_true
rbias_bb0_iwls_glm_ee_2 <- (bb0_iwls_glm_ee_2 - bb0_true)/bb0_true
boxplot(rbias_bb0,rbias_bb0_iwls_ee_2_3,rbias_bb0_iwls_glm_ee_2)
#Unbiased but large variance, IWLS does NOT reduce variance

bias_bb1 <- (bb1_est - bb1_true)
bias_bb1_iwls_ee_2_3 <- (bb1_iwls_ee_2_3 - bb1_true)
bias_bb1_iwls_glm_ee_2 <- (bb1_iwls_glm_ee_2 - bb1_true)
boxplot(bias_bb1,bias_bb1_iwls_ee_2_3,bias_bb1_iwls_glm_ee_2)

# Add for conditional variance
bias_cc0_est_ee_2 <- (cc0_est_ee_2 - cc0_true)
bias_cc0_est_ee_2_3 <- (cc0_est_ee_2_3 - cc0_true)
bias_cc0_est_ee_2_3_MoM <- (cc0_est_ee_2_3_MoM - cc0_true)
bias_cc0_iwls_ee_2_3 <- (cc0_iwls_ee_2_3  - cc0_true)
bias_cc0_est_glm_ee_2 <- (cc0_est_glm_ee_2 - cc0_true)
bias_cc0_iwls_glm_ee_2 <- (cc0_iwls_glm_ee_2 - cc0_true)

boxplot(bias_cc0_est_ee_2,
        bias_cc0_est_ee_2_3,
        bias_cc0_est_ee_2_3_MoM,
        bias_cc0_iwls_ee_2_3,
        bias_cc0_est_glm_ee_2,
        bias_cc0_iwls_glm_ee_2)

bias_cc1_est_ee_2 <- (cc1_est_ee_2 - cc1_true)
bias_cc1_est_ee_2_3 <- (cc1_est_ee_2_3 - cc1_true)
bias_cc1_est_ee_2_3_MoM <- (cc1_est_ee_2_3_MoM - cc1_true)
bias_cc1_iwls_ee_2_3 <- (cc1_iwls_ee_2_3  - cc1_true)
bias_cc1_est_glm_ee_2 <- (cc1_est_glm_ee_2 - cc1_true)
bias_cc1_iwls_glm_ee_2 <- (cc1_iwls_glm_ee_2 - cc1_true)

boxplot(bias_cc1_est_ee_2,
        bias_cc1_est_ee_2_3,
        bias_cc1_est_ee_2_3_MoM,
        bias_cc1_iwls_ee_2_3,
        bias_cc1_est_glm_ee_2,
        bias_cc1_iwls_glm_ee_2)



# Compute relative RMSE of regeesion coefficients


# Compute coverage of known value

#ci_lo_res_2 <- 
#ci_hi_res_2 <-
#ci_lo_res_2_3 <-
#ci_hi_res_2_3 <- 

#if (ci_lo_res_2[i] <= yy_100_at_tt_max_TRUE & ci_hi_res_2[i] >= yy_100_at_tt_max_TRUE == TRUE){
#  covrg_res_2[i] == 1
#} else {
#  covrg_res_2[i] == 0
#}
#sum(covrg_res2)

# 

# Compare regression coefficients
#mean(bb0_est)
#mean(bb1_est)
#mean(cc0_est_ee_2)
#mean(cc1_est_ee_2)
#mean(cc0_est_ee_2_3)
#mean(cc1_est_ee_2_3)



# Compute RRMSEs


# Compare 95% region overlap


