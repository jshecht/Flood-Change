# Magnification factors

# Givens
nn <- 50
tt <- seq(1,nn,1)

cv <- seq(0.1,1,0.1)
mu <- 5
rho_1 <- seq(0,1,0.1)
rho_2 <- seq(0,1,0.1)

zp0_10 <- qnorm(0.99,0,1)

# Create arrays to store outputs

bb <- array(NA,dim=c(length(cv),length(rho_1)))
var_ee <- array(NA,dim=c(length(cv),length(rho_1)))
var_ee_abs <- array(NA,dim=c(length(cv),length(rho_1)))
var_ee_2 <- array(NA,dim=c(length(cv),length(rho_1)))

cc0 <- array(NA,dim=c(length(cv),length(rho_1),length(rho_2)))
cc1 <- array(NA,dim=c(length(cv),length(rho_1),length(rho_2)))

cc0_t2 <- array(NA,dim=c(length(cv),length(rho_1),length(rho_2)))
cc1_t2 <- array(NA,dim=c(length(cv),length(rho_1),length(rho_2)))

M_10 <- array(NA,dim=c(length(cv),length(rho_1),length(rho_2)))
M_10_mean <- array(NA,dim=c(length(cv),length(rho_1),length(rho_2)))
M_10_var <- array(NA,dim=c(length(cv),length(rho_1),length(rho_2)))
M_10_prod <- array(NA,dim=c(length(cv),length(rho_1),length(rho_2)))
M_10_t2 <- array(NA,dim=c(length(cv),length(rho_1),length(rho_2)))
var_term_1 <- array(NA,dim=c(length(cv),length(rho_1),length(rho_2)))
var_term_2 <- array(NA,dim=c(length(cv),length(rho_1),length(rho_2)))

zpc_10 <- array(NA,dim=c(length(cv),length(rho_1),length(rho_2)))
Tc_10 <- array(NA,dim=c(length(cv),length(rho_1),length(rho_2)))

for (i in 1:length(cv)){
  for (j in 1:length(rho_1)){
    for (k in 1:length(rho_2)){
    
    # Derived conditional mean terms
    var_yy <- (mu * cv[i])^2
    bb[i,j] = sqrt(rho_1[j]^2 * var_yy/var(tt))
    var_ee[i,j] = var_yy * ((1 - rho_1[j]^2))
    
    # Since ee^2 arises from a gamma distribution
    var_ee_2[i,j] <- 2*var_ee[i,j]^4
    
    var_ee_abs[i,j] <- var_ee[i,j]*(1 - 2/pi)
    # How well does this perform when ee is not normal? 

    # Derived conditional variance terms
    cc1[i,j,k] = sqrt(rho_2[k]^2 * var_ee_2[i,j]/var(tt))
    
    cc0[i,j,k] = var_ee[i,j] - (rho_2[k] * sqrt(var_ee_2[i,j]) * mean(tt))/sd(tt)
    
    #cc1_t2[i,j,k] = rho_2[k] * sqrt(var_ee_abs[i,j])/sd(tt)
  
    #cc0_t2[i,j,k] = var_ee[i,j]*sqrt(2/pi) - (rho_2[k] * sqrt(var_ee_abs[i,j]) * mean(tt))/sd(tt)
    
    cc1_t2[i,j,k] = sqrt(rho_2[k]^2 * var_ee_2[i,j]/var(tt^2))
    cc0_t2[i,j,k] = var_ee[i,j] - (rho_2[k] * sqrt(var_ee_2[i,j]) * mean(tt^2))/sd(tt^2)
      
      
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
    
    #M_10_t2[i,j,k] = exp(bb[i,j]*delta_t + qnorm(0.99,0,1)*(cc1_t2[i,j,k] * delta_t))
    
    
    
    M_10_t2[i,j,k] = exp(bb[i,j]*delta_t + qnorm(0.99,0,1)*(v_sign1_t2 * sqrt(abs(cc0_t2[i,j,k] + cc1_t2[i,j,k] * (t_ref+delta_t)^2))
                                                         - v_sign2_t2 * sqrt(abs(cc0_t2[i,j,k] + cc1_t2[i,j,k] * t_ref^2)))) 
    
    
    #M_10_t2[i,j,k] = exp(bb[i,j]*delta_t + qnorm(0.99,0,1)*(cc1_t2[i,j,k]*delta_t))
    
    zpc_10[i,j,k] = zp0_10 * (sqrt(var_ee[i,j]) + sqrt(cc1[i,j,k]*(t_ref - t_bar)) - bb[i,j]*(max(tt) - t_ref))/(sqrt(var_ee[i,j]) + sqrt(cc1[i,j,k]*(t_cur - t_bar))) 
    
    Tc_10[i,j,k] = 1/(1-pnorm(zpc_10[i,j,k],0,1))
    
    }  
  }
}

# CORRELATION CONTOURS
contour(rho_1,rho_2,M_10[3,,],levels=c(1.25,1.5,2,2.5,3,4,5,10),
        labels = NULL,
        xlim = range(rho_1, finite = TRUE),
        ylim = range(rho_2, finite = TRUE),
        zlim = range(M_10[3,,], finite = TRUE),
        xlab="Correlation of Trend in Mean",ylab="Correlation of Trend in Cv",main="Effects of correlation coefficients \n on magnification factors \n Linear Model (mu = 5, Cv = 0.3)",cex.main=0.9,
        labcex = 0.6, drawlabels = TRUE, method = "flattest", 
        col = gray.colors(5), lty = par("lty"), lwd = par("lwd"),
        add = FALSE)
#points(1,0.01,pch=8)
#text(0.8,0.05,"M",cex=0.75) 

# CORRELATION CONTOURS t2
contour(rho_1,rho_2,M_10_t2[3,,],levels=c(1.25,1.5,2,3,4,5,10,20),
        labels = NULL,
        xlim = range(rho_1, finite = TRUE),
        ylim = range(rho_2, finite = TRUE),
        zlim = range(M_10_t2[3,,], finite = TRUE),
        xlab="Correlation of Trend in Mean",ylab="Correlation of Trend in Cv",main="Effects of correlation coefficients \n on magnification factors \n Quadratic Model (mu = 5, Cv = 0.3)",cex.main=0.9,
        labcex = 0.6, drawlabels = TRUE, method = "flattest", 
        col = gray.colors(5), lty = par("lty"), lwd = par("lwd"),
        add = FALSE)
#points(1,0.01,pch=8)
#text(0.8,0.05,"M",cex=0.75) 

# AVERAGE RECURRENCE INTERVAL
contour(rho_1,rho_2,Tc_10[3,,],levels=c(2,5,10,15,20,25,30,35,40,45,50,75,100),
        labels = NULL,
        xlim = range(rho_1, finite = TRUE),
        ylim = range(rho_2, finite = TRUE),
        zlim = range(Tc_10[3,,], finite = TRUE),
        xlab="Correlation of Trend in Mean",ylab="Correlation of Trend in Cv",main="Effects of correlation coefficients \n on average recurrence interval \n Linear Model (mu = 5, Cv = 0.3)",cex.main=0.9,
        labcex = 0.6, drawlabels = TRUE, method = "flattest", 
        col = gray.colors(5), lty = par("lty"), lwd = par("lwd"),
        add = FALSE)
#points(1,0.01,pch=8)
#text(0.8,0.05,"M",cex=0.75) 

# CORRELATION CONTOURS M_10
contour(M_10_mean[1,,1],M_10_var[1,1,],M_10_prod[1,,],levels=c(1.05,1.1,1.2,1.3,1.4,1.5,2,5,10),
        labels = NULL,
        xlim = range(M_10_mean[1,,1], finite = TRUE),
        ylim = range(M_10_var[1,1,], finite = TRUE),
        zlim = range(M_10_prod[1,,], finite = TRUE),
        xlab = "Magnification from Trend in Mean",ylab="Magnification from Trend in Cv",main="Decadal Magnification for Cv = 0.5",
        labcex = 0.6, drawlabels = TRUE, method = "flattest", 
        col = gray.colors(5), lty = par("lty"), lwd = par("lwd"),
        add = FALSE)


#points(1,0.01,pch=8)
#text(0.8,0.05,"M",cex=0.75) 

par(mfrow=c(2,2))
par(mfrow=c(1,1))
