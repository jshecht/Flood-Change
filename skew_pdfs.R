qp <- seq(0.0005,0.9995,0.0005)
d1 <- vector(mode = "numeric",length=length(qp))
d2 <- vector(mode = "numeric",length=length(qp))

for (v in 1:length(qp)){
  d1[v] <- qnorm(qp[v],3,1)
  d2[v] <- qnorm(qp[v],5,0.5)
}
  
skew_res <- skewness(c(d1,d2))

# Plot different pdfs
par(mfrow=c(3,1))

x <- seq(0,8,length = 100)
hx1 <- dnorm(x,3,1)
hx2 <- dnorm(x,5,1.5)
hx3 <- (hx1 + hx2)/2
par(mar=c(2,4.1,2,2))
plot(x, hx3, type="l", lty=2, ylim=c(0,1), yaxt='n', xlab="",ylab="Density")
text(0.5,0.75,"COMBINED")
plot(x, hx2, type="l", lty=2, ylim=c(0,1), yaxt='n', xlab="",ylab="Density")
text(0.5,0.75,"PERIOD 2")
plot(x, hx1, type="l", lty=2, ylim=c(0,1), yaxt='n', xlab="",ylab="Density")
text(0.5,0.75,"PERIOD 1")

par(mfrow=c(1,1))
#