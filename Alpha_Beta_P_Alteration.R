aa <- seq(0.01,0.99,0.01)
bb <- seq(0.01,0.99,0.01)

P_NA_CNA <- array(NA,dim=c(length(aa),length(bb)))
P_NA_CA <- array(NA,dim=c(length(aa),length(bb)))
P_A_CNA <- array(NA,dim=c(length(aa),length(bb)))
P_A_CA <- array(NA,dim=c(length(aa),length(bb)))

# For a non-informative prior
k = 0.5

for (i in 1:length(aa)){
  for (j in 1:length(bb)){
    P_NA_CNA[i,j] = (1-k)*(1 - aa[i])/((1-k)*(1-aa[i])+k*bb[j])
    P_NA_CA[i,j] = (1-k)*aa[i]/((1-k)*aa[i]+k*(1-bb[j]))
    P_A_CNA[i,j] = k*bb[j]/((1-k)*(1-aa[i])+k*bb[j])
    P_A_CA[i,j] = k*(1-bb[j])/((1-k)*aa[i]+k*(1-bb[j]))
  }
}

par(mfrow=c(2,2))
contour(aa,bb,P_NA_CNA,xlab="Type I Error Prob. (alpha)",ylab="Type II Error Prob. (beta)"
        ,main="P(NA|CNA)")
contour(aa,bb,P_NA_CA,xlab="Type I Error Prob. (alpha)",ylab="Type II Error Prob. (beta)"
        ,main="P(NA|CA)")
contour(aa,bb,P_A_CNA,xlab="Type I Error Prob. (alpha)",ylab="Type II Error Prob. (beta)"
        ,main="P(A|CNA)")
contour(aa,bb,P_A_CA,xlab="Type I Error Prob. (alpha)",ylab="Type II Error Prob. (beta)"
        ,main="P(A|CA)")
par(mfrow=c(1,1))