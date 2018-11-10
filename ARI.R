# Magnification factors

# Givens
nn <- 50
tt <- seq(1,nn,1)

cv_real <- c(0.01,seq(0.25,2.5,0.25)) # Test Cv's from 0.25 to 0.50, also 0.01
rho_1 <- seq(-1,1,0.1) 
rho_2 <- seq(-1,1,0.1)

zp0_10 <- qnorm(0.99,0,1) 

# Create arrays to store outputs
var_yy <- array(NA,dim=c(length(cv_real)))

bb <- array(NA,dim=c(length(cv_real),length(rho_1)))
var_ee <- array(NA,dim=c(length(cv_real),length(rho_1)))
var_ee_2 <- array(NA,dim=c(length(cv_real),length(rho_1)))

cc0 <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))
cc1 <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))

M_10 <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))
M_10_mean <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))
M_10_var <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))
M_10_prod <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))

var_term_1 <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))
var_term_2 <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))

zpc_10 <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))
Tc_10 <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))

for (i in 1:length(cv_real)){
  for (j in 1:length(rho_1)){
    for (k in 1:length(rho_2)){
    
    # Derived conditional mean terms
    var_yy[i] = log(1+cv_real[i]^2)
    
    bb[i,j] = sign(rho_1[j])*sqrt(rho_1[j]^2 * var_yy[i]/var(tt))
    var_ee[i,j] = var_yy[i] * ((1 - rho_1[j]^2))
    
    # Since ee^2 arises from a gamma distribution
    var_ee_2[i,j] <- 2*var_ee[i,j]^4
    
    # How well does this perform when ee is not normal? 

    # Derived conditional variance terms
    cc1[i,j,k] = rho_2[k] * sqrt(var_ee_2[i,j]/var(tt))
    
    cc0[i,j,k] = var_ee[i,j] - (rho_2[k] * sqrt(var_ee_2[i,j]) * mean(tt))/sd(tt)
    
    delta_t <- 10 
    t_cur <- max(tt) 
    t_ref <- max(tt) - 10 
    t_bar <- mean(tt) 
    
    zpc_10[i,j,k] = (zp0_10 * sqrt(var_ee[i,j]) * sqrt(1 + cc1[i,j,k]*(t_ref - t_bar)) - bb[i,j]*(max(tt) - t_ref))/(sqrt(var_ee[i,j]) * sqrt(1 + cc1[i,j,k]*(t_cur - t_bar)))  
    
    #zpc_10[i,j,k] = ((zp0_10 * sqrt(var_ee[i,j]) + sign(cc1[i,j,k])*sqrt(abs(cc1[i,j,k]*(t_ref - t_bar)))) - bb[i,j]*(max(tt) - t_ref))/(sqrt(var_ee[i,j]) + sign(cc1[i,j,k])*sqrt(abs(cc1[i,j,k]*(t_cur - t_bar))))  
    
    Tc_10[i,j,k] = 1/(1-pnorm(zpc_10[i,j,k],0,1))
    
    }     
  } 
}

# PLOTS ----------------------------------------


# ADJUSTED RECURRENCE INTERVAL FOR DIFFERENT CVs FOR THE LINEAR MODEL

# Set margins for first plot
par(mar=c(5.1, 4.1, 4.1, 2.1),xpd=FALSE) #Creates more space, allows for legend to the right of plot

plot(cv_real,Tc_10[,12,12],typ="l",lty=1,xlim=c(0,2.5),ylim=c(0,250),xlab="Real-space Cv",ylab="Adj. Recurrence Interval (ARI)")
title("Linear model")
lines(cv_real,Tc_10[,14,14],lty=2) 
lines(cv_real,Tc_10[,10,10],lty=1,col="gray")
lines(cv_real,Tc_10[,8,8],lty=2,col="gray") 
abline(h=100,lty=3,lwd=0.5)
leg.txt <- c(expression(~ rho[1] * " = 0.1," ~ rho[2] * " = 0.1" ),
             expression(~ rho[1] * " = 0.3," ~ rho[2] * " = 0.3" ),
             expression(~ rho[1] * " = 0.1," ~ rho[2] * " = -0.1" ), 
             expression(~ rho[1] * " = 0.3," ~ rho[2] * " = -0.3" ))
legend(x="topleft",legend=leg.txt,lty=c(1,2,1,2),col=c("black","black","gray","gray"),cex=0.6,bty="n")
abline(h=100,lty=6)

# Reset plot margins
par(mar=c(5.1, 4.1, 4.1, 2.1))


