#Heteroscedastic error

u <- runif(100,0,50)
w <- rnorm(100,0,1)
plot(u,0.01*w*u,xlab="Time or Covariate",ylab="Residuals")
abline(h=0,col="gray")

k_time <- seq(1,100,1)
k_rnorm <- rnorm(100,5,1)
k_beta <- 0.05

k_q_time <- k_rnorm + k_beta * k_time
plot(k_time,k_q_time,xlab="Time or Covariate",ylab="Annual Peak Flow",cex.axis=0.8)
lines(k_time,fitted(lm(k_q_time~k_time)),col="blue",lwd=2)

#xaxt="n",yaxt="n",ann=FALSE