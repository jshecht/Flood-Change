# Compute power of a regression model

# Set working directory
setwd("C:/Users/jhecht01/Box Sync/USACE/R_Analysis")

# Import data
Q_table = read.csv("01102500_R_Import.csv")

t_year<-Q_table[1:74,3]
lnQ_amax<-Q_table[1:74,6]

ns_trend.lm <-lm(lnQ_amax~t_year)
cor(t_year,lnQ_amax)

pwr.r.test(r=0.45,n=74,sig.level=0.000028,alternative="greater")