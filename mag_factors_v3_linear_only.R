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
var_ee_abs <- array(NA,dim=c(length(cv_real),length(rho_1)))
var_ee_2 <- array(NA,dim=c(length(cv_real),length(rho_1)))

cc0 <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))
cc1 <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))

M_10 <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))
M_10_mean <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))
M_10_var <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))
M_10_prod <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))
M_10_t2 <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))
var_term_1 <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))
var_term_2 <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))

zpc_10 <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))
Tc_10 <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))

zpc_10_t2 <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))
Tc_10_t2 <- array(NA,dim=c(length(cv_real),length(rho_1),length(rho_2)))

for (i in 1:length(cv_real)){
  for (j in 1:length(rho_1)){
    for (k in 1:length(rho_2)){
    
    # Derived conditional mean terms
    var_yy[i] = log(1+cv_real[i]^2)
    
    bb[i,j] = sign(rho_1[j])*sqrt(rho_1[j]^2 * var_yy[i]/var(tt))
    var_ee[i,j] = var_yy[i] * ((1 - rho_1[j]^2))
    
    # Since ee^2 arises from a gamma distribution
    var_ee_2[i,j] <- 2*var_ee[i,j]^4
    
    var_ee_abs[i,j] <- var_ee[i,j]*(1 - 2/pi)
    # How well does this perform when ee is not normal? 

    # Derived conditional variance terms
    cc1[i,j,k] = rho_2[k] * sqrt(var_ee_2[i,j]/var(tt))
    
    cc0[i,j,k] = var_ee[i,j] - (rho_2[k] * sqrt(var_ee_2[i,j]) * mean(tt))/sd(tt)
    
    cc1_t2[i,j,k] = rho_2[k] * sqrt(var_ee_abs[i,j])/sd(tt)
    
    cc0_t2[i,j,k] = var_ee[i,j]*sqrt(2/pi) - (rho_2[k] * sqrt(var_ee_abs[i,j]) * mean(tt))/sd(tt)
    
    #cc1_t2[i,j,k] = rho_2[k] * sqrt(var_ee_2[i,j]/var(tt^2))
    #cc0_t2[i,j,k] = var_ee[i,j] - (rho_2[k] * sqrt(var_ee_2[i,j]) * mean(tt^2))/sd(tt^2)
    
    # Make all negative values of cc0 infeasible since the variance cannot have a negative intercept
    if (cc0[i,j,k] < 0) {
      cc0[i,j,k] = sqrt(cc0[i,j,k]) # sets infeasible negative value
    } 
    
    if (cc0_t2[i,j,k] < 0){
      cc0_t2[i,j,k] = sqrt(cc0_t2[i,j,k]) # sets infeasible negative value
    }
    
    delta_t <- 10 
    t_cur <- max(tt) 
    t_ref <- max(tt) - 10 
    t_bar <- mean(tt) 
    
    
    v_sign1 <- sign(cc0[i,j,k] + cc1[i,j,k] *(t_ref+delta_t))
    v_sign2 <- sign(cc0[i,j,k] + cc1[i,j,k] * t_ref)
    
    v_sign1_t2 <- sign(cc0[i,j,k] + cc1[i,j,k] *(t_ref+delta_t)^2)
    v_sign2_t2 <- sign(cc0[i,j,k] + cc1[i,j,k] * t_ref^2)
    #var_term_1[i,j,k] = sqrt(cc0[i,j,k] + cc1[i,j,k] * (t_ref+delta_t))
    #var_term_2[i,j,k] = sqrt(cc0[i,j,k] + cc1[i,j,k] * t_ref)
    
    M_10[i,j,k] = exp(bb[i,j]*delta_t + qnorm(0.99,0,1)*(sqrt(cc0[i,j,k] + cc1[i,j,k] * (t_ref+delta_t))
                                                       - sqrt(cc0[i,j,k] + cc1[i,j,k] * t_ref)))
    M_10_mean[i,j,k] = exp(bb[i,j]*delta_t)
    M_10_var[i,j,k] = exp(qnorm(0.99,0,1)*(v_sign1 * sqrt(abs(cc0[i,j,k] + cc1[i,j,k]*(t_ref+delta_t)))
                                    - v_sign2 * sqrt(abs(cc0[i,j,k]+cc1[i,j,k]*t_ref))))
    M_10_prod[i,j,k] = M_10_mean[i,j,k] * M_10_var[i,j,k]
    
    M_10_t2[i,j,k] = exp(bb[i,j]*delta_t + qnorm(0.99,0,1)*(cc1_t2[i,j,k] * delta_t))
    
    
    
    #M_10_t2[i,j,k] = exp(bb[i,j]*delta_t + qnorm(0.99,0,1)*(sqrt(cc0_t2[i,j,k] + cc1_t2[i,j,k] * (t_ref+delta_t)^2)
    #                                                    - sqrt(cc0_t2[i,j,k] + cc1_t2[i,j,k] * t_ref^2)))
    
    
    zpc_10[i,j,k] = (zp0_10 * sqrt(var_ee[i,j]) * sqrt(1 + cc1[i,j,k]*(t_ref - t_bar)) - bb[i,j]*(max(tt) - t_ref))/(sqrt(var_ee[i,j]) * sqrt(1 + cc1[i,j,k]*(t_cur - t_bar)))  
    
    #zpc_10[i,j,k] = ((zp0_10 * sqrt(var_ee[i,j]) + sign(cc1[i,j,k])*sqrt(abs(cc1[i,j,k]*(t_ref - t_bar)))) - bb[i,j]*(max(tt) - t_ref))/(sqrt(var_ee[i,j]) + sign(cc1[i,j,k])*sqrt(abs(cc1[i,j,k]*(t_cur - t_bar))))  
    
    Tc_10[i,j,k] = 1/(1-pnorm(zpc_10[i,j,k],0,1))
    
    zpc_10_t2[i,j,k] = (zp0_10 * sqrt(var_ee[i,j]) * (1 + cc1_t2[i,j,k]*(t_ref - t_bar)) - bb[i,j]*(max(tt) - t_ref))/(sqrt(var_ee[i,j]) * (1 + cc1_t2[i,j,k]*(t_cur - t_bar)))
      
    Tc_10_t2[i,j,k] = 1/(1-pnorm(zpc_10_t2[i,j,k],0,1))
    
    }     
  } 
}

# PLOTS ----------------------------------------

# Make correlation space plots
par(mfrow=c(2,2))

# CORRELATION CONTOURS (Cv = 0.5)
contour(rho_1,rho_2,M_10[2,,],levels=c(0.5,0.6,0.75,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8),
        labels = NULL,
        xlim = range(rho_1, finite = TRUE),
        ylim = range(rho_2, finite = TRUE),
        zlim = range(M_10[2,,], finite = TRUE),
        xlab="Corr. trend in mean",ylab="Corr. trend in Cv",main="Linear Model (Cv = 0.5)",cex.main=0.95,
        labcex = 0.6, drawlabels = TRUE, method = "flattest", 
        col = gray.colors(5), lty = par("lty"), lwd = par("lwd"),
        add = FALSE)
#points(1,0.01,pch=8)
#text(0.8,0.05,"M",cex=0.75) 

# CORRELATION CONTOURS (Cv = 1.5)
contour(rho_1,rho_2,M_10[6,,],levels=c(0,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3),
        labels = NULL,
        xlim = range(rho_1, finite = TRUE),
        ylim = range(rho_2, finite = TRUE),
        zlim = range(M_10[6,,], finite = TRUE),
        xlab="Corr. trend in mean",ylab="Corr. Trend in Cv",main="Linear Model (Cv = 1.5)",cex.main=0.95,
        labcex = 0.6, drawlabels = TRUE, method = "flattest", 
        col = gray.colors(5), lty = par("lty"), lwd = par("lwd"),
        add = FALSE)
#points(1,0.01,pch=8)
#text(0.8,0.05,"M",cex=0.75) 

# CORRELATION CONTOURS t2 (Cv = 0.5)
contour(rho_1,rho_2,M_10_t2[3,,],levels=c(0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.8),
        labels = NULL,
        xlim = range(rho_1, finite = TRUE),
        ylim = range(rho_2, finite = TRUE),
        zlim = range(M_10_t2[3,,], finite = TRUE),
        xlab="Corr. trend in mean",ylab="Corr. trend in Cv",main="Quadratic Model (Cv = 0.5)",cex.main=0.95,
        labcex = 0.6, drawlabels = TRUE, method = "flattest", 
        col = gray.colors(5), lty = par("lty"), lwd = par("lwd"),
        add = FALSE)
#points(1,0.01,pch=8)
#text(0.8,0.05,"M",cex=0.75) 

# CORRELATION CONTOURS t2 (Cv = 1.5)
contour(rho_1,rho_2,M_10_t2[7,,],levels=c(0,0.1,0.25,0.5,0.75,1,1.25,1.5,1.75,2,2.25,2.5,2.75,3),
        labels = NULL,
        xlim = range(rho_1, finite = TRUE),
        ylim = range(rho_2, finite = TRUE),
        zlim = range(M_10_t2[7,,], finite = TRUE),
        xlab="Corr. trend in mean",ylab="Corr. trend in Cv",main="Quadratic Model (Cv = 1.5)",cex.main=0.95,
        labcex = 0.6, drawlabels = TRUE, method = "flattest", 
        col = gray.colors(5), lty = par("lty"), lwd = par("lwd"),
        add = FALSE)
#points(1,0.01,pch=8)
#text(0.8,0.05,"M",cex=0.75) 

mtext("Magnification Factors", side = 3, line = -1, outer = TRUE)

par(mfrow=c(1,2))

# AVERAGE RECURRENCE INTERVAL (Cv = 0.5)
contour(rho_1,rho_2,Tc_10[3,,],levels=c(1.5,2,5,10,25,50,100,200,500,1000,10000),
        labels = NULL,
        xlim = range(rho_1, finite = TRUE),
        ylim = range(rho_2, finite = TRUE),
        zlim = range(Tc_10[3,,], finite = TRUE),
        xlab="Corr. trend in Mean",ylab="Corr. trend in Cv",main="Linear Model (Cv = 0.5)",cex.main=0.95,
        labcex = 0.6, drawlabels = TRUE, method = "flattest", 
        col = gray.colors(5), lty = par("lty"), lwd = par("lwd"),
        add = FALSE)
#points(1,0.01,pch=8)
#text(0.8,0.05,"M",cex=0.75) 

# AVERAGE RECURRENCE INTERVAL (Cv = 1.5)
contour(rho_1,rho_2,Tc_10[7,,],levels=c(1.5,2,5,10,25,50,100,200,500,1000,10000),
        labels = NULL,
        xlim = range(rho_1, finite = TRUE),
        ylim = range(rho_2, finite = TRUE),
        zlim = range(Tc_10[7,,], finite = TRUE),
        xlab="Corr. trend in mean",ylab="Corr. trend in Cv",main="Linear Model (Cv = 1.5)",cex.main=0.95,
        labcex = 0.6, drawlabels = TRUE, method = "flattest", 
        col = gray.colors(5), lty = par("lty"), lwd = par("lwd"),
        add = FALSE)
#points(1,0.01,pch=8)
#text(0.8,0.05,"M",cex=0.75) 

mtext("Adjusted recurrence intervals for 100-year flood", side = 3, line = -1, outer = TRUE)

par(mfrow=c(1,1))

# ADJUSTED RECURRENCE INTERVAL FOR DIFFERENT CVs FOR THE LINEAR MODEL

# Set margins for first plot
par(mar=c(5.1, 4.1, 4.1, 2.1),xpd=FALSE) #Creates more space, allows for legend to the right of plot

plot(cv_real,Tc_10[,12,12],typ="l",lty=1,xlim=c(0,2.5),ylim=c(0,250),xlab="Real-space Cv",ylab="Adj. Recurrence Interval (ARI)")
title("Linear model")
lines(cv_real,Tc_10[,14,14],lty=2) 
lines(cv_real,Tc_10[,10,10],lty=1,col="gray")
lines(cv_real,Tc_10[,8,8],lty=2,col="gray") 
#abline(h=100,lty=3,lwd=0.5)
leg.txt <- c(expression(~ rho[1] * " = 0.1," ~ rho[2] * " = 0.1" ),
             expression(~ rho[1] * " = 0.3," ~ rho[2] * " = 0.3" ),
             expression(~ rho[1] * " = 0.1," ~ rho[2] * " = -0.1" ), 
             expression(~ rho[1] * " = 0.3," ~ rho[2] * " = -0.3" ))
legend(x="topleft",legend=leg.txt,lty=c(1,2,1,2),col=c("black","black","gray","gray"),cex=0.6,bty="n")
abline(h=100,lty=6)

# Reset plot margins
par(mar=c(5.1, 4.1, 4.1, 2.1))

# ADJUSTED RECURRENCE INTERVAL FOR DIFFERENT CVs FOR THE QUADRATIC MODEL

par(mar=c(5.1, 4.1, 4.1, 2.1),xpd=FALSE) #Creates more space, allows for legend to the right of plot

plot(cv_real,Tc_10_t2[,12,12],typ="l",lty=1,xlim=c(0,2.5),ylim=c(0,250),xlab="Real-space Cv",ylab="")
title("Quadratic model")
lines(cv_real,Tc_10_t2[,14,14],lty=2) 
lines(cv_real,Tc_10_t2[,12,10],lty=1,col="gray")
lines(cv_real,Tc_10_t2[,14,9],lty=2,col="gray") 
#abline(h=100,lty=3,lwd=0.5)
leg.txt <- c(expression(~ rho[1] * " = 0.1," ~ rho[2] * " = 0.1" ),
             expression(~ rho[1] * " = 0.3," ~ rho[2] * " = 0.3" ),
             expression(~ rho[1] * " = 0.1," ~ rho[2] * " = -0.1" ), 
             expression(~ rho[1] * " = 0.3," ~ rho[2] * " = -0.3" ))
legend(x="topleft",legend=leg.txt,lty=c(1,2,1,2),col=c("black","black","gray","gray"),cex=0.6,bty="n")

# Reset plot margins
par(mar=c(5.1, 4.1, 4.1, 2.1))

# Reset to 1 x 1 panel


