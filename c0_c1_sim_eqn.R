# Tester

c0 <- 0.5
c1 <- 0.05
tt <- runif(10000,1,50)
vt <- rnorm(10000,0,1)

et <- vt*(c0+c1*tt)^(3/2)
vt_2_3 <- (vt^2)^(1/3)
et_2_3 <- vt_2_3*(c0+c1*tt)
mean(et_2_3) 
var(et_2_3) 

# Analytical value of mean(et_2_3)
0.802*var(vt)^(1/3)*(c0+c1*mean(tt)) 

# Analytical value of var(et_2_3)
var(vt_2_3)*var(c0+c1*tt) + mean(vt_2_3)^2*var(c0+c1*tt) + var(vt_2_3)*mean(c0+c1*tt)^2 #3.408

# Compute mean of LHS
0.802*var(et)^(1/3)
mean(et_2_3)

# Compute mean of RHS
0.802*var(vt)^(1/3)*(c0+c1*mean(tt))


# Compute variance of RHS
0.188*var(vt)^(2/3)*var(c0+c1*tt)+
  0.188*var(vt)^(2/3)*mean(c0+c1*tt)^2 +
  (0.802*var(vt)^(1/3))^2*var(c0+c1*tt)

# Compute variance of LHS
0.188*var(et)^(2/3)
