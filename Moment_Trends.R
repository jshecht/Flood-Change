# Moment trends

bb_mean_0 = 600
bb_sd_0 = 150
bb_trend_mean = 10
bb_trend_sd = 0.5

tt = seq(1,100,1)
bb_mean_tt <- bb_mean_0 + bb_trend_mean*tt
bb_sd_tt <- bb_sd_0 + bb_trend_sd*tt

q100y_tt = bb_mean_tt + 2.326*bb_sd_tt 

raw2nd_tt = bb_sd_tt^2 + bb_mean_tt^2

cv_tt = bb_sd_tt/bb_mean_tt

plot(tt,raw2nd_tt)

# Regression implies

