# CHECK WITH ANA'S PAPER

# Inputs
aa <- 0.114  #Type I error prob
nn <- 73    #Record length
rr <- 0.142  #Correlation coefficient

# Computations
dd <- 1/(sqrt(1/rr^2-1)) 
tt <- qt(1-aa,nn-2)
t2_error <- pt(tt - dd*sqrt(nn),nn-2)
t2_error

# As the correlation increases and the 