# GUIDE TO SCRIPTS FOR NONSTATIONARY FLOOD FREQUENCY ANALYSIS 

# Packages for analysis 
install.packages(car,glm2,
                 lmomco,Lmoments,lmtest)


# car is essential for computing the Durbin-Watson statistic, also POSSIBLY other regression things
# e1071 enables skewness (surprised it is not a base package function)
# glm2 estimates glm models including the gamma GLM with a log link function that we use
  # this function is better than glm at negotiating model convergence issues
# lmomco and Lmoments are useful for distribution function estimate ???
# lmtest has some auxiliary tools for linear regression models ???

# Packages for mapping
install.packages(GISTools,mapproj,maps,maptools)

# DATA ANALYSIS SCRIPTS

# NWIS_QPeak 
  # Imports annual maximum series from USGS stations stored in a .csv file
  # Get actual uSGS --> R downloading tools
# Record comparison 
  # Main script for comparing records at different stations
  # Takes subset of stations from NWIS_QPeak and produces quantile estimates from
  # (i) a linear conditional mean model and (ii) different conditional residual variance models
  # Includes basic two-step regressio models with squared residuals, the original derivation method
  # with the derived expression for sigma(v), the current two-step OLS models, OLS models with IWLS,
  # GLM gamma models with IWLS
  # Also processes results of the analysis 
# Flood_quantile_chart 
  # For case study analysis 
  # Produces a 1 x 3 plot of the (i) conditional mean model, (ii) conditional 
  # variance model and (iii) a chart comparing the flood frequency estimated with different 
  # Also examines sub-period differences 
# Station_maps
  # Produces map of stations with significant increasing trends in the mean with concurrent trends
  # in the mean and variance 

# Older sheets
  # Q_01102500_Aberjona
    # Analysis of different derived conditional variance models for the Aberjona River, also estimation comparison

# MODEL EXPERIMENTS

# Monte Carlo estimation
  # Preliminary Monte Carlo experiments for different values of bb0, bb1, cc0 and cc1 with different record lengths

# IWLS 
  # Template for setting up an iteratively weighted least squares model 
  # glm.fit also estimates this way, but easier to make my own code 

# Leverage 
  # Explores effects that leverage may have on residual estimates

# mag_factors
  # Explores effects of different stationary Cv and trend model correlation coefficient on flood magnification factors
  # Source of contour plots for flood magnification factors and adjusted recurrence intervals

# Power_transf_bias
  # Evaluates performance of different methods for dealing with transformation bias associated with Anscombe 
  # residual model (e^(2/3))

# my_PPCC_test_lmomco_heo
  # My own PPCC test using the lmomco package and critical values from Heo et al. (2008)

# PPCC_test_norm_USGS_Lorenz
  # PPCC test from the USGS

# Ansc_ppcc_plot
  # Generates PPCC chart with 1,000,000 samples 



# GLM_test is an OLDER sheet used to learn about different GLM models with different link functions

