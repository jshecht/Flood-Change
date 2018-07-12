# ESTIMATE FLOOD QUANTILES -----------------------------------------------------------

# Extract the time series for a particular station
k = 11 # Set k
yr_start <- 1951 # set start year (water year)

peak_flow <- peaks[[station_group[k,1]]]/35.3146662127
year <- seq(yr_start,yr_start+length(peak_flow),1)

# Create standard normal variates at 0.2% intervals
z_p <- qnorm(seq(1,99,1)/100,0,1)

# Create output vectors  
qtile_stnry <- vector(mode="numeric",length=length(z_p))
qtile_mean <- vector(mode="numeric",length=length(z_p))
qtile_linear <- vector(mode="numeric",length=length(z_p))
qtile_logt <- vector(mode="numeric",length=length(z_p))
qtile_quad <- vector(mode="numeric",length=length(z_p))

for (z in 1:length(z_p)){
 qtile_stnry[z] = exp(mean(log(peak_flow + 0.00001)) + z_p[z] * sd(log(peak_flow+0.01)))
 qtile_mean[z] = exp(cond_med_at_n[k] + z_p[z]* sd(res_reg1[[k]]))
 qtile_logt[z] = exp(cond_med_at_n[k] + z_p[z] *sqrt(cond_var_at_n_5d[k]))
 qtile_linear[z] = exp(cond_med_at_n[k] + z_p[z] *sqrt(cond_var_at_n_5a[k]))
 qtile_quad[z] = exp(cond_med_at_n[k] + z_p[z] *sqrt(cond_var_at_n_5b[k]))
}
 
 qtile_100y <- c(qtile_stnry[99],
               qtile_mean[99],
               qtile_linear[99],
               qtile_logt[99],
               qtile_quad[99])
 
# # STONY BROOK NEAR WEST SUFFIELD, CT 
# plot(1/(1-pnorm(z_p,0,1)),qtile_linear,typ="l",xlim=c(1,100),ylim=c(0.0001,1.2*max(qtile_100y)),log="x",xlab="Return Period (Yrs)",
#    ylab="Peak Flow (m3s1)",lty=3,cex.axis=0.8)
# title(main="Current flood frequency curves \n Stony Brook near West Suffield, CT", cex.main=1.0)
# #mtext("Stony Brook near West Suffield, CT",cex=0.8)
# lines(1/(1-pnorm(z_p,0,1)),qtile_stnry,lty=1)
# lines(1/(1-pnorm(z_p,0,1)),qtile_mean,lty=2)
# #lines(1/(1-pnorm(z_p,0,1)),qtile_logt,lty=4)
# #lines(1/(1-pnorm(z_p,0,1)),qtile_quad,lty=5)
# abline(h=max(peak_flow),lty=2,col="gray")
# leg.txt <- c("Mean + Cv (Linear Model)","Mean Only","Stationary")
# legend(x="topleft",legend=leg.txt,lty=c(3,2,1),cex=0.75,bty="n")
# text(5,2100,"Flood of record",cex=0.7)
# text(40,100,"USGS 04166000",cex=0.6) 
# 
# 
# # Create panel with trend in mean, trend in variance, quantile
# Start year

 
# Conditional mean
par(mfrow=c(1,3))
plot(seq(yr_start,(yr_start-1)+length(peak_flow),1),log(peak_flow+0.0001),
     xlab="Year in record",
     ylab=expression(ln~"(Annual" ~ peak ~ "flow," ~ m^{3}/s~")"))  
lines(seq(yr_start,(yr_start-1)+length(peak_flow),1),
      b0_reg1[k] + b1_reg1[k]*seq(1,length(peak_flow),1))
title("Conditional mean")

# Conditional Cv plot
plot(seq(yr_start,(yr_start-1)+length(peak_flow),1),(res_reg1[[k]]^2)^(1/3),xlab="Year in record",ylab="")
title("Conditional Cv",
      ylab=expression(paste("[Conditional mean residuals,"~epsilon,"]"^(2/3))),
      mgp=c(2.5,1,0))
lines(seq(yr_start,(yr_start-1)+length(peak_flow),1),
      cc0_5a[k] + cc1_5a[k]*(seq(1,length(peak_flow),1))) 
  
 
# Quantile plot
par(mgp=c(3,1,0)) 
plot(1/(1-pnorm(z_p,0,1)),qtile_mean,typ="l",xlim=c(1,100),ylim=c(0.0001,1.2*max(qtile_100y)),
     log="x",xlab="Return Period (Yrs)",
     ylab=expression(~Annual ~ peak ~ flow ~ (m^{3}/s)),lty=2)
title("Frequency comparison")  
lines(1/(1-pnorm(z_p,0,1)),qtile_stnry,lty=1) 
lines(1/(1-pnorm(z_p,0,1)),qtile_linear,lty=3,col="gray50") 
lines(1/(1-pnorm(z_p,0,1)),qtile_logt,lty=6,col="gray50") 
#lines(1/(1-pnorm(z_p,0,1)),qtile_quad,lty=5)
abline(h=max(peak_flow),lty=2,col="gray") 

# For increasing trend in the mean, increasing trend in the Cv
#leg.txt <- c("Mean + Cv","Mean + Cv (ln(t))","Mean Only","Stationary")
#legend(x="topleft",legend=leg.txt,lty=c(3,6,2,1),col=c("gray50","gray50","black","black"),cex=0.85,bty="n")
#text(3,60,"Flood of record",cex=0.8)
#text(30,5,"USGS 01184100",cex=0.8)

# For increasing trend in the mean, decreasing trend in the Cv
leg.txt <- c("Mean Only","Mean + Cv (ln(t))","Stationary","Mean + Cv") 
legend(x="topleft",legend=leg.txt,lty=c(2,6,1,3),cex=0.85,bty="n",col=c("black","gray50","black","gray50"))
text(3,41,"Flood of record",cex=0.8)
text(30,5,"USGS 04166000",cex=0.8)



par(mfrow=c(1,1))
 


# Sub-period analysis for Birmingham, Michigan
#par(mfrow=c(1,3))
#hist(log(peak_flow),xlab="ln(Annual peak flow, m3s1)",ylim=c(0,16),
#    breaks=seq(4,8,0.25),main="1951-2014",
#    cex.main=0.9,cex.axis=0.8,cex.lab=0.9,col="gray")
#hist(log(peak_flow[1:32]),xlab="ln(Annual peak flow, , m3s1)",ylim=c(0,16),
#    breaks=seq(4,8,0.25),main="1951-1982",
#    cex.main=0.9,cex.axis=0.8,cex.lab=0.9,col="gray")
##hist(log(peak_flow[33:64]),xlab="ln(Annual peak flow, , m3s1)",ylim=c(0,16),
#    breaks=seq(4,8,0.25),main="1983-2014",
#    cex.main=0.9,cex.axis=0.8,cex.lab=0.9,col="gray")
#par(mfrow=c(1,1))
 

