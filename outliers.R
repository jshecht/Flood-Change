# Conditional mean model with studentized residuals and outlier diagnostics
#dfbetas(flood_reg1.lm)
#dffits(flood_reg1.lm)
#cooks.distance(flood_reg1.lm)
#hatvalues(flood_reg1.lm)

# Compute studentized residuals
#rstudent(flood_reg1.lm)

# Untransform fitted values
#rstudent * var() * sqrt(1-hatvalues)