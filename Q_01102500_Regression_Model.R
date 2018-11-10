# Aberjona analysis 

# DATA PROCESSING -----------------------------------------------------------------

# Set working directory
setwd("C:/Users/Jory/Box Sync/USACE/R_Analysis")

# Import annual peak flows as .csv
Q_01102500_imp <- read.csv("Qpeak_01102500_dl_20160507.csv")

# Compute two-column matrix with just water years (Oct 1-Sep 30) and flows
Q_01102500_a <- Q_01102500_imp[,2:3]

# Add data from water year 2011 (observation on 08-10-2011), which was not available
# in recent USGS download
data_2011 <- c(2011,418)

Q_01102500_b <- as.matrix(rbind(Q_01102500_a[1:71,],data_2011,Q_01102500_a[72:75,]))

# Compute natural logs of annual peak flows
Q_01102500_c <- cbind(Q_01102500_b[,1],log(Q_01102500_b[,2]))
t <- 1:nrow(Q_01102500_c)
ln_Q <- Q_01102500_c[,2]

ln_Q <- as.numeric(log(peaks[[6327]])[3:80])
t <- 1:length(ln_Q)

plot(t+1939,Q_01102500_b[,2],typ="l",col="blue",xlab="Water Year in Record (t)",ylab="Annual Peak Flow (cfs)")
title("Annual peak flows in Aberjona River \n at Winchester, MA (USGS 01102500)",cex=0.95)
text("Annual Maximum Series",cex =0.8)

# TIME-DEPENDENT MODEL OF THE MEAN --------------------------------------------

# Compute regression model
ln_Q_mean_reg <-lm(ln_Q~t)
summary(ln_Q_mean_reg)
ln_Q_mean_reg_coeff <- as.numeric(coefficients(ln_Q_mean_reg))
ln_Q_mean_reg_b0 <- ln_Q_mean_reg_coeff[1]
ln_Q_mean_reg_b1 <- ln_Q_mean_reg_coeff[2]
ln_Q_mean_reg_fitted <- fitted(ln_Q_mean_reg)
ln_Q_mean_reg_residuals <- residuals(ln_Q_mean_reg)
ln_Q_mean_reg_influence <- influence(ln_Q_mean_reg)
ln_Q_mean_reg_cor <- cor(t,ln_Q)

# Transfer to Mathcad
# a, b, rho, n, sigma(y)
mcad_mean_reg <- c(ln_Q_mean_reg_b0,
                   ln_Q_mean_reg_b1,
                   ln_Q_mean_reg_cor,
                   length(t),
                   sd(ln_Q))

# Plot regression model
plot(t,ln_Q,xlab="Water year in record (t)",ylab="ln(Annual Peak Flow, cfs)")
lines(t,ln_Q_mean_reg_fitted,col="blue",lwd=2)
title("Regression model of conditional mean")
# Add R^2 and p-value
text(8,7.2,expression(~R^2 == 0.202),cex=0.6)
text(9,6.9,expression(p < 0.001),cex=0.6)
text(57,4.8,"ln(Q)=5.469+0.0124*t",cex=0.8)


# Plot residuals
plot(t+1939,ln_Q_mean_reg_residuals,xlab="Year",ylab="Residuals")
#abline(h=0,lty=2,col="red")
abline(v=1955.5,lty=2,col="red")
title("Residuals of conditional mean model")

# Check normality of residuals with PPCC test
ppcc_Q_mean_reg_residuals <- ppcc.test(ln_Q_mean_reg_residuals)
pval_ppcc_Q_mean_reg_residuals <- as.numeric(ppcc_Q_mean_reg_residuals[2])
if (pval_ppcc_Q_mean_reg_residuals < 0.05) {
  Q_mean_res_norm = FALSE
  #stop("Residuals not normally distributed")
} else {
  Q_mean_res_norm = TRUE
}

# TRANsFORM RESIDUALS ---------------------------------------------------------------

# Apply Anscombe transformation to residuals
res_2 <- (ln_Q_mean_reg_residuals^2)
res_2_3 <- res_2^(1/3)
t_1_3 <- t^(1/3)
t_2_3 <- t^(2/3)

# Plot transformed residuals
plot(t,res_2_3,xlab="Water Year in Record (t)",ylab=expression(~epsilon^(2/3)))
title("Anscombe Transformed Residuals")
lines(t,fitted(lm(res_2_3~t)))
text(6,1.1,"p = 0.025",cex=0.6,col="blue")
text(7,0.8,"R^2 = 0.07",cex=0.6,col="blue")
summary(lm(res_2_3~t))


plot(t_1_3,res_2_3,xlab="t^(1/3)",ylab="r^(2/3)")


# Plot squared residuals
plot(t,res_2,xlab="Water Year in Record (t)",ylab=expression(~epsilon^2))
title("Squared Residuals")

PPCC.test(res_2_3)
hist(res_2)
hist(res_2_3)

# RUN CONDITIONAL VARIANCE MODELS ----------------------

# Initialize arrays for storing outputs
model_rsquare <- array(NA,dim=6)
model_pval <- array(NA,dim=6)

# Model 1: OLS fit with t^1/3
ols_2stage_t_1_3 <- lm(res_2_3 ~ t_1_3)
summary(ols_2stage_t_1_3)
ols_2stage_t_1_3_coeff <- as.matrix(coefficients(ols_2stage_t_1_3))
ols_2stage_t_1_3_b0 <- ols_2stage_t_1_3_coeff[1,]
ols_2stage_t_1_3_b1 <- ols_2stage_t_1_3_coeff[2,]
ols_2stage_t_1_3_fit <- fitted(ols_2stage_t_1_3)
plot(t_1_3,res_2_3,xlab=expression(~t^(1/3),ylab=expression(~epsilon^(2/3))))
lines(t_1_3,ols_2stage_t_1_3_fit)
# Compute Rsq
ss_tot <- sum((res_2_3-mean(res_2_3))^2)
ss_res <- sum((res_2_3-ols_2stage_t_1_3_fit)^2)
model_rsquare[1] = 1-ss_res/ss_tot
# Compute p-value
se_b <- sqrt(sum((res_2_3-ols_2stage_t_1_3_fit)^2)/(length(t)-2)/(sum((t_1_3-mean(t_1_3))^2)))
t.value <- ols_2stage_t_1_3_b1/se_b
model_pval[1] = dt(t.value, df=length(t) - 2)

# Model 2: oLS fit with t^2/3
ols_2stage_t_2_3 <- lm(res_2_3 ~ t_2_3)
summary(ols_2stage_t_2_3)
ols_2stage_t_2_3_coeff <- as.matrix(coefficients(ols_2stage_t_2_3))
ols_2stage_t_2_3_int <- ols_2stage_t_2_3_coeff[1]
ols_2stage_t_2_3_b <- ols_2stage_t_2_3_coeff[2]
ols_2stage_t_2_3_fit <- fitted(ols_2stage_t_2_3)
plot(t_2_3,res_2_3,xlab=expression(~t^(2/3),ylab=expression(~epsilon^(2/3))))
lines(t_2_3,ols_2stage_t_2_3_fit)
# Compute Rsq
ss_tot <- sum((res_2_3-mean(res_2_3))^2)
ss_res <- sum((res_2_3-ols_2stage_t_1_3_fit)^2)
model_rsquare[2] = 1-ss_res/ss_tot
# Compute p-value
se_b <- sqrt(sum((res_2_3-ols_2stage_t_2_3_fit)^2)/(length(t)-2)/(sum((t_2_3-mean(t_2_3))^2)))
t.value <- ols_2stage_t_1_3_b1/se_b
model_pval[2] = dt(t.value, df=length(t) - 2)

# Model 3: Derived fit with t^1/3
nu_2_t_1_3 <- (ln_Q_mean_reg_b1^2*(length(t)-1)^2)/(6*(length(t)+1))*(1-ln_Q_mean_reg_cor^2)/ln_Q_mean_reg_cor^2
res_var_t_1_3 <- nu_2_t_1_3^(1/3)*t_1_3
plot(t_1_3,res_2_3)
lines(t_1_3,res_var_t_1_3)
# Compute rsquare
ss_tot <- sum((res_2_3-mean(res_2_3))^2)
ss_res <- sum((res_2_3-res_var_t_1_3)^2)
model_rsquare[3]= 1-ss_res/ss_tot
# Compute p-value
se_b <- sqrt(sum((res_2_3-res_var_t_1_3)^2)/(length(t)-2)/(sum((t_1_3-mean(t_1_3))^2)))
#t.value <- res_var_t_1_3/se_b
model_pval[3] = dt(t.value, df=length(t) - 2)

# Model 4: Derived fit with t^2/3
nu_2_t_2_3 <- (ln_Q_mean_reg_b1^2*(length(t)-1)^2)/(4*(length(t)^2+length(t)+1))*(1-ln_Q_mean_reg_cor^2)/ln_Q_mean_reg_cor^2
res_var_t_2_3 <- nu_2_t_2_3^(1/3)*t_2_3
plot(t_2_3,res_2_3)
lines(t_2_3,res_var_t_2_3)
ss_tot <- sum((res_2_3-mean(res_2_3))^2)
ss_res <- sum((res_2_3-res_var_t_2_3)^2)
model_rsquare[4] = 1-ss_res/ss_tot
# Compute p-value
#se_b <- sqrt(sum((res_2_3-res_var_t_2_3)^2)/(length(t)-2)/(sum((t_2_3-mean(t_2_3))^2)))
#t.value <- res_var_t_2_3/se_b
#model_pval[4] = dt(t.value, df=length(t) - 2)


# Model 5: Derived fit with t^1/3, c0, c1

#Aberjona
#c0_1_3 = 312800
#c1_1_3 = 8125

#Skokie
#c0_1_3 = 812800
#c1_1_3 = 33180

# Whiteoak Bayou (constrain c1 < 0)
c0_1_3 = 1352000
c1_1_3 = -5000
t_1_3_c <- (c0_1_3+c1_1_3*t)^(1/3)

#For sqrt(c0+c1*t)
nu_2_t_1_3_c <- (ln_Q_mean_reg_b1^2*(length(t)-1)^2)/(6*(2*c0_1_3+c1_1_3*(length(t)+1)))*(1-ln_Q_mean_reg_cor^2)/ln_Q_mean_reg_cor^2

#For c0+c1*sqrt(t)
#nu_2_t_1_3_c <- (1/12)*(ln_Q_mean_reg_b1^2*(length(t)-1)^2)/(c0_1_3^2+2*c0_1_3*c1_1_3*(2*(length(t)^(3/2)-1))/(3*(length(t)-1))+c1_1_3^2*((length(t)+1)/2))*(1-ln_Q_mean_reg_cor^2)/ln_Q_mean_reg_cor^2 

res_var_t_1_3_c <- nu_2_t_1_3_c^(1/3)*t_1_3_c
# nu^2^(1/3)=nu^(2/3)

par(mar=c(5.1,5.1,4.1,2.1))
plot(t_1_3_c,res_2_3,main="Residual Variance Model 1",xlab=expression(~(c0+c1*t)^(1/3)),ylab=expression(~sigma[epsilon]^(2/3)))
lines(t_1_3_c,res_var_t_1_3_c)

# Compute rsquare
ss_tot <- sum((res_2_3-mean(res_2_3))^2)
ss_res <- sum((res_2_3-res_var_t_1_3_c)^2)
model_rsquare[5] = 1-ss_res/ss_tot
print(model_rsquare[5])
# Compute p-value
#se_b <- sqrt(sum((res_2_3-res_var_t_1_3_c)^2)/(length(t)-2)/(sum((t_1_3_c-mean(t_1_3_c))^2)))
#t.value <- res_var_t_1_3_c/se_b
#model_pval[5] = dt(t.value, df=length(t) - 2)

mae <- mean(abs(res_2_3-res_var_t_1_3_c))
mae2 <- mean(abs(res_2_3-mean(res_2_3)))


# Model 6: Derived fit with t^2/3,c0,c1

#Aberjona 
#c0_2_3 = 18860
#c1_2_3 = 393.426

#Skokie
#c0_2_3 = 18690
#c1_2_3 = 611
  
#Whiteoak Bayou
c0_2_3 = 458100
c1_2_3 = -2053

t_2_3_c <- (c0_2_3+c1_2_3*t)^(2/3)

nu_2_t_2_3_c <- (ln_Q_mean_reg_b1^2*(length(t)-1)^2)/(4*((3*c0_2_3^2)+(c1_2_3^2*(length(t)^2))+(3*c0_2_3*c1_2_3+c1_2_3^2)*(length(t)+1)))*(1-ln_Q_mean_reg_cor^2)/ln_Q_mean_reg_cor^2
res_var_t_2_3_c <- nu_2_t_2_3_c^(1/3)*(c0_2_3+c1_2_3*t)^(2/3)

# Plot
plot(t_2_3_c,res_2_3)
lines(t_2_3_c,res_var_t_2_3_c)

par(mar=c(5.1,5.1,4.1,2.1))
plot(t_2_3_c,res_2_3,main="Residual Variance Model 2",xlab=expression(~(c0+c1*t)^(2/3)),ylab=expression(~sigma[epsilon]^(2/3)))
lines(t_2_3_c,res_var_t_2_3_c)
(max(res_var_t_2_3_c)-min(res_var_t_2_3_c))/(max(t_2_3_c)-min(t_2_3_c))

# Compute R^2
ss_tot <- sum((res_2_3-mean(res_2_3))^2)
ss_res <- sum((res_2_3-res_var_t_2_3_c)^2)
model_rsquare[6] = 1-ss_res/ss_tot
print(model_rsquare[6])

mae <- mean(abs(res_2_3-res_var_t_2_3_c))
mae2 <- mean(abs(res_2_3-mean(res_2_3)))

# Compute p-value 
#se_b <- sqrt(sum((res_2_3-res_var_t_2_3_c)^2)/(length(t)-2)/(sum((t_2_3_c-mean(t_2_3_c))^2)))
#b_res_var_t_2_3_c <- sum(t_2_3_c*res_var_t_2_3_c)/sum(t_2_3_c^2)
#t.value <- b_res_var_t_2_3_c/se_b
#model_pval[6] = dt(t.value, df=length(t) - 2)

# SUB-PERIOD MODELS
t_post_hwy <- t[17:76]
ln_Q_post_hwy <- ln_Q[17:76]

# Linear model
post_hwy.lm <- lm(ln_Q_post_hwy~t_post_hwy)
summary(post_hwy.lm)
post_hwy_b0 <- summary(post_hwy.lm)$coefficients[1,1]
post_hwy_b1 <- summary(post_hwy.lm)$coefficients[2,1]
post_hwy_fitted <- fitted(post_hwy.lm)

# Check residuals
post_hwy_res <- residuals(post_hwy.lm)
plot(post_hwy_res)
post_hwy_res_2_3 <- (post_hwy_res^2)^(1/3)
plot(post_hwy_res_2_3)

# Estimate quantile
post_hwy_var_res <- sqrt(var(ln_Q_post_hwy)-post_hwy_b1^2*var(t_post_hwy))
post_hwy_Q99 <- max(post_hwy_fitted) + qnorm(0.99,0,1)*post_hwy_var_res
exp(post_hwy_Q99)
#2302 cfs


# ESTIMATE FLOOD QUANTILES -----------------------------------------------------------

# Create standard normal variates at 0.2% intervals
z_p <- qnorm(seq(1,99,1)/100,0,1)

# Create output vectors  
qtile_stnry <- vector(mode="numeric",length=length(z_p))
qtile_mean <- vector(mode="numeric",length=length(z_p))
qtile_t_1_3_c <- vector(mode="numeric",length=length(z_p))
qtile_t_2_3_c <- vector(mode="numeric",length=length(z_p))
qtile_post_hwy <- vector(mode="numeric",length=length(z_p))
qtile_ols_1_3 <- vector(mode="numeric",length=length(z_p))
qtile_ols_2_3 <- vector(mode="numeric",length=length(z_p))

for (z in 1:length(z_p)){
  qtile_stnry[z] = exp(mean(ln_Q)+z_p[z]*sd(ln_Q))
  qtile_mean[z] = exp(ln_Q_mean_reg_b0+ln_Q_mean_reg_b1*t[length(t)]+z_p[z]*sqrt(var(ln_Q)-ln_Q_mean_reg_b1^2*sd(t)^2))
  qtile_t_1_3_c[z] = exp(ln_Q_mean_reg_b0+ln_Q_mean_reg_b1*t[length(t)]+z_p[z]*res_var_t_1_3_c[length(t)]^(3/2))
  qtile_t_2_3_c[z] = exp(ln_Q_mean_reg_b0+ln_Q_mean_reg_b1*t[length(t)]+z_p[z]*res_var_t_2_3_c[length(t)]^(3/2))
  qtile_post_hwy[z]= exp(max(post_hwy_fitted) + z_p[z]*post_hwy_var_res)
  qtile_ols_1_3[z] = exp(ln_Q_mean_reg_b0+ln_Q_mean_reg_b1*t[length(t)]+z_p[z]*ols_2stage_t_1_3_fit[length(t)]^(3/2))
  qtile_ols_2_3[z] = exp(ln_Q_mean_reg_b0+ln_Q_mean_reg_b1*t[length(t)]+z_p[z]*ols_2stage_t_2_3_fit[length(t)]^(3/2))
}

qtile_100y <- c(qtile_stnry[99],
                qtile_mean[99],
                qtile_t_1_3_c[99],
                qtile_t_2_3_c[99],
                qtile_post_hwy[99],
                qtile_ols_1_3[99],
                qtile_ols_1_3[99])


plot(1/(1-pnorm(z_p,0,1)),qtile_t_1_3_c,typ="l",xlim=c(1,100),ylim=c(0,5000),log="x",xlab="Return Period (Yrs)",
     ylab="Peak Flow (cfs)",col="red",cex.axis=0.8)
title(main="Current flood frequency curves: \n Aberjona River at Winchester, MA", cex.main=1.0)
#mtext("Aberjona River at Winchester, MA (USGS 01102500)",cex=0.8)
lines(1/(1-pnorm(z_p,0,1)),qtile_t_2_3_c,col="blue")
lines(1/(1-pnorm(z_p,0,1)),qtile_stnry,col="black")
lines(1/(1-pnorm(z_p,0,1)),qtile_mean,col="green4")
lines(1/(1-pnorm(z_p,0,1)),qtile_post_hwy,col="purple",lwd=2.5)
#lines(1/(1-pnorm(z_p,0,1)),qtile_ols_1_3,col="red",lty=2)
#lines(1/(1-pnorm(z_p,0,1)),qtile_ols_2_3,col="blue",lty=2)
abline(h=1590,lty=2,col="gray")
leg.txt <- c("Mean + Cv (Model 2)", "Mean + Cv (Model 1)", "Mean Only","Stationary","Mean Only (1956-2015)")
legend(x="topleft",legend=leg.txt,lty=1,cex=0.75,col=c("blue","red","green4","black","purple"),bty="n")
text(2,1900,"Flood of record",cex=0.7)
text(40,100,"USGS 01102500",cex=0.6)

# Step change analysis 
ln_Q.ts <- ts(ln_Q)
vvalue <- cpt.var(ln_Q.ts,penalty="None",test.stat="Normal")

# Compute nonstationary Cv

# Results analysis
#barplot(qtile_100y,names=c("Stationary","Mean Only","Var Model 1","Var Model 2"),cex.axis= 1, cex.lab=1)
#barplot(model_rsquare,cex.axis= 1, cex.lab=1)
#barplot(model_pval)

# GRAPH WITH PREDICTION INTERVALS
ln_Q_mean_reg_se <- sqrt(var(ln_Q) - ln_Q_mean_reg_b1^2*var(t))
t_PI <- seq(1940,2030,1)

PI_low_ln_Q <- (ln_Q_mean_reg_b0 + ln_Q_mean_reg_b1*(t_PI-1939)) + ln_Q_mean_reg_se * qt(0.025,length(t_PI)-2)*sqrt((1 + 1/length(t_PI)) + (t_PI-mean(t_PI))^2/var(t_PI))
PI_high_ln_Q <- (ln_Q_mean_reg_b0 + ln_Q_mean_reg_b1*(t_PI-1939)) + ln_Q_mean_reg_se * qt(0.975,length(t_PI)-2)*sqrt((1 + 1/length(t_PI)) + (t_PI-mean(t_PI))^2/var(t_PI))

plot(t+1939,Q_01102500_b[,2],log="y",xlab="Year",ylab="Annual Peak Flow (cfs)",xlim=c(1940,2030),ylim=c(10,10000))
lines(t+1939,exp(ln_Q_mean_reg_fitted))
lines(t_PI,exp(PI_low_ln_Q),lty=2)
lines(t_PI,exp(PI_high_ln_Q),lty=2)
abline(v=2015,lty=2,col="red")
text(1960,8000,"Type I error prob < 0.001%",cex=0.6)
text(1958,5000,"Type II error prob = 48%",cex=0.6)
text(2008,40,"Current",cex=0.75)
text(2006,25,"conditions",cex=0.75)

# ESTIMATE MULTIVARIATE MODEL
read.csv("Aberjona_Analysis_Model_Prep.csv")

