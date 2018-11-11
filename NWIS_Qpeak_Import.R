# Import NWIS data

install.packages(c("car","e1071","EnvStats","glm2","lmtest","sandwich","trend"))
library("car")
library("e1071")
library("EnvStats")
library("glm2")
library("lmtest")
library("sandwich")
library("trend")
library("WriteXLS")

# testing

setwd("C:/Users/joryh/ownCloud/documents/Tufts/WRR_Flood_Nonstationarity/Q_Data_Analysis/NSFFA")

# FIX: SOURCING OF PPCC TEST
source("PPCC_test_norm_USGS_Lorenz.R")

############################## IMPORT GAGES II STATIONS #########################################

# Import list of stations with information about urbanization attributes from GAGES II
gages2 <- read.csv("gages2_UrbanAttrib.csv",colClasses = c(rep("character",2),rep("numeric",8)))

# Check column names and classes
colnames(gages2)
colnames(gages2) <- c("STAID",colnames(gages2[2:10])) # Fix first column name

lapply(gages2,class)

# Create a subset of urban gages
urbGages <- gages2[which(gages2$STO_RATIO <= 0.1 & 
                         gages2$IMPNLCD06 > 10),]

# Export this subset for work in Excel and ArcGIS
write.csv(urbGages, file="urbGages.csv")

# Import GIS results
bhc_ts <- read.csv("bhc_time_series.csv")
colnames(bhc_ts) <- c("STAID_TEXT",colnames(bhc_ts[2:14]))

# Put back zeros in text station IDs removed through .csv format
for(i in 1:nrow(bhc_ts)){
  if (i <= 534 | i == 537){
    bhc_ts$STAID_TEXT[i] = paste("0",bhc_ts$STAID_TEXT[i],sep="")
  } 
}

# Remove stations with less than 5 km2 basin
bhc_ts <- bhc_ts[!bhc_ts$DRAIN_SQKM < 5, ]

# Convert TIAs to numeric values
bhc_ts$IC1940 <- as.numeric(bhc_ts$IC1940) 
bhc_ts$IC1950 <- as.numeric(bhc_ts$IC1950)
bhc_ts$IC1960 <- as.numeric(bhc_ts$IC1960)
bhc_ts$IC1970 <- as.numeric(bhc_ts$IC1970)
bhc_ts$IC1980 <- as.numeric(bhc_ts$IC1980)
bhc_ts$IC1990 <- as.numeric(bhc_ts$IC1990)
bhc_ts$IC2000 <- as.numeric(bhc_ts$IC2000)
bhc_ts$IC2010 <- as.numeric(bhc_ts$IC2010)

# Create matrix for storing interpolated data
bhc_ts_ann <- array(NA,dim=c(nrow(bhc_ts),76))

# Interpolate TIA estimates for 1941-2016
for (j in 1:nrow(bhc_ts)){
  for (k in 1:7){
    for (i in 1:10){
      bhc_ts_ann[j,(k-1)*10+i] = bhc_ts[j,k+5] + (bhc_ts[j,k+5+1] - bhc_ts[j,k+5])*i/10
    }
  }
  # Interpolate six years past 2010
  for (i in 1:6){
    bhc_ts_ann[j,k*10+i] = bhc_ts[j,k+5] + (bhc_ts[j,k+5+1] - bhc_ts[j,k+5])*(10+i)/10
  }
}

# Append site info and 1940 data
df_bhc_ts_ann <- cbind(bhc_ts[,1:6],bhc_ts_ann)
bhc_ts_ann_mat <- as.matrix(df_bhc_ts_ann[,6:82])

# Check work with a plot
plot(seq(1940,2016,1),bhc_ts_ann_mat[1,],typ="l",
     ylim=c(0,50),
     xlab="Year",
     ylab="TIA (%)")
for(i in 502:510){
  lines(seq(1940,2016,1),bhc_ts_ann_mat[i,])
}


############################### IMPORT ANNUAL FLOOD SERIES ######################################

# Import list of stations with at least 30 years of data (file too big to store in git repo)
Qpeak_30yrs <- read.csv("C:/Users/joryh/ownCloud/documents/Tufts/WRR_Flood_Nonstationarity/Q_data_analysis/NSFFA_backup/Qpeak_TimeSeries_v3.csv")
#Qpeak_30yrs <- read.csv("Export_30yr_Feb2018.csv",colClasses=c(rep("character",2)))
Qpeak_30yrs <- Qpeak_30yrs[1:505399,]
nrow(Qpeak_30yrs )

# Check column names and classes
colnames(Qpeak_30yrs)
colnames(Qpeak_30yrs) <- c("agency_cd","30yr_Min",colnames(Qpeak_30yrs[3:19]))
Qpeak_30yrs$site_no <- as.character(Qpeak_30yrs$site_no)

# Add zeros to stations with HUCs beginning with zero
Qpeak_30yrs$site_no <- c(paste("0",Qpeak_30yrs$site_no[1:389160],sep=""),
                         Qpeak_30yrs$site_no[389161:400646],
                         paste("0",Qpeak_30yrs$site_no[400647:400676],sep=""),
                         Qpeak_30yrs$site_no[400677:nrow(Qpeak_30yrs)])

lapply(Qpeak_30yrs,class)

############################### IDENTIFY URBAN STATIONS WITH >30 YRS DATA ##############################
# Find stations in both databases

bhc_ts_30yrQ <- bhc_ts[bhc_ts$STAID_TEXT %in% Qpeak_30yrs$site_no,]
Qpeak_urb <- Qpeak_30yrs[Qpeak_30yrs$site_no %in% bhc_ts$STAID_TEXT,]
Qpeak_urb <- Qpeak_urb[!is.na(Qpeak_urb$peak_va)==TRUE,]

# Create an empty list to store results. List format needed for records of different lengths
main_list <- list()
for (i in 1:nrow(bhc_ts_30yrQ)){
  main_list[[i]] = as.matrix(Qpeak_urb[Qpeak_urb$site_no == bhc_ts_30yrQ$STAID_TEXT[i],])
}

# Determine stations with fewer than 30 years of data
# FIX: Problem with input data (shouldn't be any stations with < 30 years)
yrs_Qpeak <- vector(mode="numeric",length=length(main_list))
main_list_2 <- list()
k = 0
for (i in 1:length(main_list)){
  yrs_Qpeak[i] = nrow(main_list[[i]])
  if (yrs_Qpeak[i] >= 30){
    k = k + 1
    main_list_2[[k]] = main_list[[i]]
  }
}

bhc_ts_urb_2 <- bhc_ts_30yrQ[yrs_Qpeak>=30,]

# Determine stations with intermittent records that are too short at start or end of record

# For each year, compute number of obs in sliding window (15-yr on each side).
wy_cont <- vector(mode="logical",length=length(main_list_2))
for(i in 1:length(main_list_2)){
  wy <- as.numeric(main_list_2[[i]][,9])
  wy_window <- vector(mode="numeric",length=length(wy))
  for(m in 1:length(wy)){
    wy_window[m] = length(wy[wy >= max(min(wy),wy[m]-15) & wy <= min(wy[m]+15,max(wy))])
  }
  if(min(wy_window) < 10){
    wy_cont[i] = TRUE
  } else {
    wy_cont[i] = FALSE 
  }
}

# Remove stations that have excessively intermittent records
main_list_3 <- list()
k = 0 
for(i in 1:length(main_list_2)){
  if(wy_cont[i] == FALSE){
    k = k + 1
    main_list_3[[k]] = main_list_2[[i]]
  }
}

# Update urban station list
bhc_ts_urb_3 <- bhc_ts_urb_2[wy_cont == FALSE,]
  
# 4. Determine stations with at least 1% TIA growth/decade during period of record

# Determine start and end year of each record
min_wy <- vector(mode="numeric",length(main_list_3))
max_wy <- vector(mode="numeric",length(main_list_3))
bhc_ts_rec_len <- array(NA,dim=c(length(main_list_3),77))
bhc_ts_Qpeak_rec <- list()
bhc_ts_TIA_rate <- vector(mode="numeric",length(main_list_3))

for(j in 1:length(main_list_3)){
  wy <- as.numeric(main_list_3[[j]][,9]) 
  min_wy[j] = max(min(wy),1940)
  max_wy[j] = max(wy)
  bhc_ts_rec_len <- bhc_ts_ann_mat[j,(min_wy[j]-1939):(max_wy[j]-1939)]
  bhc_ts_Qpeak_rec[[j]] = bhc_ts_rec_len
  bhc_ts_TIA_rate[j] = (max(bhc_ts_rec_len) - min(bhc_ts_rec_len))/length(bhc_ts_rec_len) * 10
  # FIX: Doesn't account for intermittent years
}

# Examine histogram of TIA growth rates
hist(bhc_ts_TIA_rate)

# Remove stations with less than 1% increase over gauging record
# Record which stations have 1% TIA or greater
main_list_4 <- list()
urb_abv_1 <- vector(mode="logical",length=length(main_list_3))
k = 0  # initializes index for new list
for(i in 1:length(main_list_3)){
  if(bhc_ts_TIA_rate[i] >= 1.00){
    k = k + 1
    urb_abv_1[i] = TRUE
    main_list_4[[k]] = main_list_3[[i]]
  } else {
    urb_abv_1[i] = FALSE
  }
}

# Create list of urban stations
bhc_ts_urb_4 <- bhc_ts_urb_3[urb_abv_1==TRUE,]

# Create urban TIA time series 
bhc_ts_TIA <- list()
k = 0
for(j in 1:nrow(bhc_ts_urb_3)){
  if(bhc_ts_TIA_rate[j] >= 1.00){
    k = k + 1
    bhc_ts_TIA[[k]] = bhc_ts_Qpeak_rec[[j]]
  }
}

# Find stations in both databases

# Check to make sure all stations are also in Qpeak_30yrs
bhc_ts_urb_final <- bhc_ts_urb_4[bhc_ts_urb_4$STAID_TEXT %in% Qpeak_30yrs$site_no,]
if(nrow(bhc_ts_urb_final)!=nrow(bhc_ts_urb_4)){
  stop("Streamflow import data error")
}

# Create annual time series of TIA
df_bhc_ts_urb_final_ann <- df_bhc_ts_ann[df_bhc_ts_ann$STAID_TEXT %in% bhc_ts_urb_final$STAID_TEXT,]

# Create list of flow records for the final list of urban stations
Qpeak_urb <- Qpeak_30yrs[Qpeak_30yrs$site_no %in% bhc_ts_urb_final$STAID_TEXT,]
Qpeak_urb <- Qpeak_urb[!is.na(Qpeak_urb$peak_va)==TRUE,]

#################################### COMPUTE FLOOD STATISTICS ######################################

# Create vectors to store output
rec_id <- seq(1,nrow(bhc_ts_urb_final),1)
peaks <- list()
peaks_urb <- list()
wy_station <- list()
bhc_ts_TIA_vec_station <- list()

beta_reg1_time <- vector(length=nrow(bhc_ts_urb_final))
beta_reg1_urb <- vector(length=nrow(bhc_ts_urb_final))

pval_beta_reg1_time <- vector(length=nrow(bhc_ts_urb_final))
pval_beta_reg1_urb <- vector(length=nrow(bhc_ts_urb_final))

pval_beta_reg1_hc3_time <- vector(length=nrow(bhc_ts_urb_final))
pval_beta_reg1_hc3_urb <- vector(length=nrow(bhc_ts_urb_final))

adq_ppcc_ln2 <- vector(length=nrow(bhc_ts_urb_final))

adq_ppcc_res_flood_reg1_time <- vector(length=nrow(bhc_ts_urb_final))
adq_dw_res_flood_reg1_time <- vector(length=nrow(bhc_ts_urb_final))
adq_res_flood_reg1_time <- vector(length=nrow(bhc_ts_urb_final))
beta_reg2_max_time <- vector(length=nrow(bhc_ts_urb_final))
pval_beta_reg2_max_time <- vector(length=nrow(bhc_ts_urb_final))
pval_beta_reg2_max_hc3_time <- vector(length=nrow(bhc_ts_urb_final))


adq_ppcc_res_flood_reg1_urb <- vector(length=nrow(bhc_ts_urb_final))
adq_dw_res_flood_reg1_urb <- vector(length=nrow(bhc_ts_urb_final))
adq_res_flood_reg1_urb <- vector(length=nrow(bhc_ts_urb_final))
beta_reg2_max_urb <- vector(length=nrow(bhc_ts_urb_final))
pval_beta_reg2_max_urb <- vector(length=nrow(bhc_ts_urb_final))
pval_beta_reg2_max_hc3_urb <- vector(length=nrow(bhc_ts_urb_final))


cc1_LIN_time <- vector(length=nrow(bhc_ts_urb_final))
cc1_QUA_time <- vector(length=nrow(bhc_ts_urb_final))
cc1_LOGT_time <- vector(length=nrow(bhc_ts_urb_final))
pval_cc1_LIN_time <- vector(length=nrow(bhc_ts_urb_final))
pval_cc1_QUA_time <- vector(length=nrow(bhc_ts_urb_final))
pval_cc1_LOGT_time <- vector(length=nrow(bhc_ts_urb_final))
pval_cc1_LIN_hc3_time <- vector(length=nrow(bhc_ts_urb_final))
pval_cc1_QUA_hc3_time <- vector(length=nrow(bhc_ts_urb_final))
pval_cc1_LOGT_hc3_time <- vector(length=nrow(bhc_ts_urb_final))

cc1_LIN_urb <- vector(length=nrow(bhc_ts_urb_final))
cc1_QUA_urb <- vector(length=nrow(bhc_ts_urb_final))
cc1_LOGT_urb <- vector(length=nrow(bhc_ts_urb_final))
pval_cc1_LIN_urb <- vector(length=nrow(bhc_ts_urb_final))
pval_cc1_QUA_urb <- vector(length=nrow(bhc_ts_urb_final))
pval_cc1_LOGT_urb <- vector(length=nrow(bhc_ts_urb_final))
pval_cc1_LIN_hc3_urb <- vector(length=nrow(bhc_ts_urb_final))
pval_cc1_QUA_hc3_urb <- vector(length=nrow(bhc_ts_urb_final))
pval_cc1_LOGT_hc3_urb <- vector(length=nrow(bhc_ts_urb_final))

cc1_GLM_time <- vector(length=nrow(bhc_ts_urb_final))
pval_cc1_GLM_time <- vector(length=nrow(bhc_ts_urb_final))
pval_cc1_GLM_hc3_time <- vector(length=nrow(bhc_ts_urb_final))

cc1_GLM_urb <- vector(length=nrow(bhc_ts_urb_final))
pval_cc1_GLM_urb <- vector(length=nrow(bhc_ts_urb_final))
pval_cc1_GLM_hc3_urb <- vector(length=nrow(bhc_ts_urb_final))

for (j in 1:nrow(bhc_ts_urb_final)){
  wy <- as.numeric(main_list_4[[j]][,9]) 
  wy_station[[j]] <- wy
  wy_order <- wy - min(wy) + 1
  wy_order_2 <- wy_order^2
  
  peak_va <- as.numeric(main_list_4[[j]][,11])
  peaks[[j]] <- peak_va
  
  # Create vector of TIA to match water years
  # FIX: bhc_ts_ann_mat
  wy_in_1940_2016 <- vector(mode="logical",length=ncol(bhc_ts_ann_mat))
  for(i in 1:77){
    if((1939 + i) %in% wy == TRUE){
      wy_in_1940_2016[i] = TRUE
    } else {
      wy_in_1940_2016[i] = FALSE
    }
  }
  
  # Determine water year order for urban equation
  wy_order_urb <- wy_in_1940_2016*seq(1,77,1)
  wy_order_urb <- (wy_order_urb[wy_order_urb!=0]-min(wy_order_urb))+1

  bhc_ts_TIA_vec <- as.numeric(df_bhc_ts_urb_final_ann [j,6:82][wy_in_1940_2016 == TRUE])
  
  bhc_ts_TIA_vec_station[[j]] <- bhc_ts_TIA_vec
  
  # Remove peak flows earlier than 1940 (in cfs)
  wy_pre1940 <- vector(mode="logical",length(wy))
  for(i in 1:length(wy)){
    if(wy[i] < 1940){
      wy_pre1940[i] = TRUE
    } else {
      wy_pre1940[i] = FALSE
    }
  }
  
  peak_va_urb <- peak_va[wy_pre1940==FALSE]
  peaks_urb[[j]] <- peak_va_urb
  
  # Build models
  flood_reg1_time.lm <- lm(log(peak_va+0.01) ~ wy_order) # Add constant to avoid problems with ephemeral streams
  flood_reg1_urb.lm <- lm(log(peak_va_urb+0.01) ~ bhc_ts_TIA_vec)
  
  # Extract coefficients
  int_reg1_time <- summary(flood_reg1_time.lm)$coefficients[1,1]
  beta_reg1_time[j] = summary(flood_reg1_time.lm)$coefficients[2,1]

  int_reg1_urb <- summary(flood_reg1_urb.lm)$coefficients[1,1]
  beta_reg1_urb[j] = summary(flood_reg1_urb.lm)$coefficients[2,1]
  
  # p-values divided by 2 since we're applying a one-sided test
  pval_int_reg1_time <- summary(flood_reg1_time.lm)$coefficients[1,4]/2
  pval_beta_reg1_time[j] = summary(flood_reg1_time.lm)$coefficients[2,4]/2
  pval_beta_reg1_hc3_time[j] = coeftest(flood_reg1_time.lm,vcov=vcovHC(flood_reg1_time.lm,"HC3"))[2,4]/2
  res_flood_reg1_time <- residuals(flood_reg1_time.lm)
  
  pval_int_reg1_urb <- summary(flood_reg1_urb.lm)$coefficients[1,4]/2
  pval_beta_reg1_urb[j] = summary(flood_reg1_urb.lm)$coefficients[2,4]/2
  pval_beta_reg1_hc3_urb[j] = coeftest(flood_reg1_urb.lm,vcov=vcovHC(flood_reg1_urb.lm,"HC3"))[2,4]/2
  res_flood_reg1_urb <- residuals(flood_reg1_urb.lm)

  # Check plausibility of LN2 assumption of AFS
  # FIX: Make sure this isn't used to filter out stations
  # Emphasize this in the text
  ppcc_ln2 <- ppcc.test(log(peak_va+0.01))
  if (is.na(as.numeric(ppcc_ln2[1]))==TRUE){
    stat_ppcc_ln2 = -9999
  } else {
    stat_ppcc_ln2 <- as.numeric(ppcc_ln2[1])
  }
  
  if (is.finite(as.numeric(ppcc_ln2[2]))==FALSE){
    pval_ppcc_ln2 = -9999
  } else {
    pval_ppcc_ln2 <- as.numeric(ppcc_ln2[2])
    if (pval_ppcc_ln2 >= 0.05){
      adq_ppcc_ln2[j] = 1
    } else {
    adq_ppcc_ln2[j] = 0
    } 
  }

  #CHECK NORMALITY OF RESIDUALS WITH PPCC TEST
  ppcc_res_flood_reg1_time <- ppcc.test(res_flood_reg1_time)
  if (is.na(as.numeric(ppcc_res_flood_reg1_time[1]))==TRUE){
    stat_ppcc_res_flood_reg1_time = -9999
  } else {
    stat_ppcc_res_flood_reg1_time <- as.numeric(ppcc_res_flood_reg1_time[1])
  }
  
  if (is.finite(as.numeric(ppcc_res_flood_reg1_time[2]))==FALSE){
    pval_ppcc_res_flood_reg1_time = -9999
  } else {pval_ppcc_res_flood_reg1_time <- as.numeric(ppcc_res_flood_reg1_time[2])
    if (pval_ppcc_res_flood_reg1_time >= 0.05){
    adq_ppcc_res_flood_reg1_time[j] = 1
    } else {
    adq_ppcc_res_flood_reg1_time[j] = 0
    } 
  }
  
  
  ppcc_res_flood_reg1_urb <- ppcc.test(res_flood_reg1_urb)
  if (is.na(as.numeric(ppcc_res_flood_reg1_urb[1]))==TRUE){
    stat_ppcc_res_flood_reg1_urb = -9999
  } else {
    stat_ppcc_res_flood_reg1_urb <- as.numeric(ppcc_res_flood_reg1_urb[1])
  }
  
  if (is.finite(as.numeric(ppcc_res_flood_reg1_urb[2]))==FALSE){
    pval_ppcc_res_flood_reg1_urb = -9999
  } else {pval_ppcc_res_flood_reg1_urb <- as.numeric(ppcc_res_flood_reg1_urb[2])
  if (pval_ppcc_res_flood_reg1_urb >= 0.05){
    adq_ppcc_res_flood_reg1_urb[j] = 1
  } else {
    adq_ppcc_res_flood_reg1_urb[j] = 0
  } 
  }
  
  # CHECK AUTOCORRELATION OF RESIDUALS WITH THE DURBIN-WATSON TEST
  
  if (max(res_flood_reg1_time)==0){
    stat_dw_res_flood_reg1_time = -9999
  } else { 
    dw_res_flood_reg1_time <- dwtest(res_flood_reg1_time~wy_order)
      if (is.finite(as.numeric(dw_res_flood_reg1_time[[1]]))==FALSE){
        stat_dw_res_flood_reg1_time = -9999
      } else {
      stat_dw_res_flood_reg1_time <- as.numeric(dw_res_flood_reg1_time[[1]])
      }
  }
  
  
  if (max(res_flood_reg1_time)==0){
    pval_dw_res_flood_reg1_time = -9999
  } else if (is.finite(as.numeric(dw_res_flood_reg1_time[[4]]))==FALSE){
    pval_dw_res_flood_reg1_time = -9999
  } else {pval_dw_res_flood_reg1_time <- dw_res_flood_reg1_time[[4]]
    if (pval_dw_res_flood_reg1_time >= 0.05){
      adq_dw_res_flood_reg1_time[j] = 1
    } else {
      adq_dw_res_flood_reg1_time[j] = 0
    }
  }
  
  
  if (max(res_flood_reg1_urb)==0){
    stat_dw_res_flood_reg1_urb = -9999
  } else { 
    dw_res_flood_reg1_urb <- dwtest(res_flood_reg1_urb~wy_order_urb)
    if (is.finite(as.numeric(dw_res_flood_reg1_urb[[1]]))==FALSE){
      stat_dw_res_flood_reg1_urb = -9999
    } else {
      stat_dw_res_flood_reg1_urb <- as.numeric(dw_res_flood_reg1_urb[[1]])
    }
  }
  
  
  if(max(res_flood_reg1_urb)==0){
    pval_dw_res_flood_reg1_urb = -9999
  } else if (is.finite(as.numeric(dw_res_flood_reg1_urb[[4]]))==FALSE){
    pval_dw_res_flood_reg1_urb = -9999
  } else {pval_dw_res_flood_reg1_urb <- dw_res_flood_reg1_urb[[4]]
  if (pval_dw_res_flood_reg1_urb >= 0.05){
    adq_dw_res_flood_reg1_urb[j] = 1
  } else {
    adq_dw_res_flood_reg1_urb[j] = 0
    }
  }

  # COMPUTE OVERALL ADEQUACY OF RESIDUALS
  adq_res_flood_reg1_time[j] = adq_ppcc_res_flood_reg1_time[j]*adq_dw_res_flood_reg1_time[j]
  
  adq_res_flood_reg1_urb[j] = adq_ppcc_res_flood_reg1_urb[j]*adq_dw_res_flood_reg1_urb[j]
  
  # Transform residuals to the two-thirds power
  res_flood_reg_time_2_3 <- (res_flood_reg1_time^2)^(1/3)
  res_flood_reg_urb_2_3 <- (res_flood_reg1_urb^2)^(1/3)
  
  # Linear model
  lm_reg2_LIN_time <- lm(res_flood_reg_time_2_3 ~ wy_order) 
  cc1_LIN_time[j] = coeftest(lm_reg2_LIN_time,vcovHC(lm_reg2_LIN_time,"HC3"))[2,1]
  pval_cc1_LIN_hc3_time[j] = coeftest(lm_reg2_LIN_time,vcovHC(lm_reg2_LIN_time,"HC3"))[2,4]
  
  lm_reg2_LIN_urb <- lm(res_flood_reg_urb_2_3 ~ wy_order_urb) 
  cc1_LIN_urb[j] = coeftest(lm_reg2_LIN_urb,vcovHC(lm_reg2_LIN_urb,"HC3"))[2,1]
  pval_cc1_LIN_hc3_urb[j] = coeftest(lm_reg2_LIN_urb,vcovHC(lm_reg2_LIN_urb,"HC3"))[2,4]
  
  
  #rho_t_res_reg2 <- cor(wy_order,res_flood_reg_2_3)
  #sd_res_flood_reg_2_3 <- sd(res_flood_reg_2_3)
  #cc0_LIN <- mean(res_flood_reg_2_3) - (rho_t_res_reg2 * sd_res_flood_reg_2_3  * mean(wy_order)) / sd(wy_order)
  #cc1_LIN[j] = rho_t_res_reg2 * sd_res_flood_reg_2_3 / sd(wy_order)
  #res_reg2_LIN_fit <- cc0_LIN + cc1_LIN[j] * wy_order
  #res_reg2_LIN_var <- length(wy_order)/(length(wy_order)-2)*var(res_flood_reg_2_3-res_reg2_LIN_fit)
  #se_cc0_LIN <- sqrt(sum((res_flood_reg_2_3 - res_reg2_LIN_fit)^2)/(length(wy_order) - 2)*(1/length(wy_order)+mean(wy_order)^2/sum((wy_order-mean(wy_order))^2)))
  #se_cc1_LIN <- sqrt(sum((res_flood_reg_2_3 - res_reg2_LIN_fit)^2)/((length(wy_order) - 2)*sum((wy_order-mean(wy_order))^2)))
  #t_cc0_LIN <- cc0_LIN/se_cc0_LIN
  #t_cc1_LIN <- cc1_LIN[j]/se_cc1_LIN
  #pval_cc0_LIN <- 2*(pt(-abs(t_cc0_LIN), df=length(wy_order) - 1))
  #pval_cc1_LIN[j] = 2*(pt(-abs(t_cc1_LIN), df=length(wy_order) - 1)) 
  
  # QUADRATIC MODEL
  lm_reg2_QUA_time <- lm(res_flood_reg_time_2_3 ~ wy_order^2) 
  cc1_QUA_time[j] = coeftest(lm_reg2_QUA_time,vcovHC(lm_reg2_QUA_time,"HC3"))[2,1]
  pval_cc1_QUA_hc3_time[j] = coeftest(lm_reg2_QUA_time,vcovHC(lm_reg2_QUA_time,"HC3"))[2,4]
  
  lm_reg2_QUA_urb <- lm(res_flood_reg_urb_2_3 ~ wy_order_urb^2) 
  cc1_QUA_urb[j] = coeftest(lm_reg2_QUA_urb,vcovHC(lm_reg2_QUA_urb,"HC3"))[2,1]
  pval_cc1_QUA_hc3_urb[j] = coeftest(lm_reg2_QUA_urb,vcovHC(lm_reg2_QUA_urb,"HC3"))[2,4]
  
  #rho_t_res_reg2 <- cor(wy_order^2,res_flood_reg_2_3) 
  #sd_res_flood_reg_2_3 <- sd(res_flood_reg_2_3)
  #cc0_QUA <- mean(res_flood_reg_2_3) - (rho_t_res_reg2 * sd_res_flood_reg_2_3  * mean(wy_order^2) / sd(wy_order^2))
  #cc1_QUA[j] = rho_t_res_reg2 * sd_res_flood_reg_2_3 / sd(wy_order_2)
  #res_reg2_QUA_fit <- cc0_QUA + cc1_QUA[j] * (wy_order_2)
  #res_reg2_QUA_var <- length(wy_order)/(length(wy_order)-2)*var(res_flood_reg_2_3-res_reg2_QUA_fit)
  #se_cc0_QUA <- sqrt(sum((res_flood_reg_2_3 - res_reg2_QUA_fit)^2)/(length(wy_order_2) - 2)*(1/length(wy_order_2)+mean(wy_order_2)^2/sum((wy_order_2-mean(wy_order_2))^2)))
  #se_cc1_QUA <- sqrt(sum((res_flood_reg_2_3 - res_reg2_QUA_fit)^2)/((length(wy_order_2) - 2)*sum((wy_order_2-mean(wy_order_2))^2)))
  #t_cc0_QUA <- cc0_QUA/se_cc0_QUA 
  #t_cc1_QUA <- cc1_QUA[j]/se_cc1_QUA 
  #pval_cc0_QUA <- 2*(pt(-abs(t_cc0_QUA), df=length(wy_order_2) - 1))
  #pval_cc1_QUA[j] = 2*(pt(-abs(t_cc1_QUA), df=length(wy_order_2) - 1)) 
  
  
  
  # LOGARITHMIC MODEL (LOG OF T, NOT LOG-LINEAR MODEL)
  lm_reg2_LOGT_time <- lm(res_flood_reg_time_2_3 ~ wy_order) 
  cc1_LOGT_time[j] = coeftest(lm_reg2_LOGT_time,vcovHC(lm_reg2_LOGT_time,"HC3"))[2,1]
  pval_cc1_LOGT_hc3_time[j] = coeftest(lm_reg2_LOGT_time,vcovHC(lm_reg2_LOGT_time,"HC3"))[2,4]
  
  lm_reg2_LOGT_urb <- lm(res_flood_reg_urb_2_3 ~ wy_order_urb) 
  cc1_LOGT_urb[j] = coeftest(lm_reg2_LOGT_urb,vcovHC(lm_reg2_LOGT_urb,"HC3"))[2,1]
  pval_cc1_LOGT_hc3_urb[j] = coeftest(lm_reg2_LOGT_urb,vcovHC(lm_reg2_LOGT_urb,"HC3"))[2,4]
  
  #rho_t_res_reg2 <- cor(log(wy_order),res_flood_reg_2_3)
  #sd_res_flood_reg_2_3 <- sd(res_flood_reg_2_3)
  #cc0_LOGT <- mean(res_flood_reg_2_3) - (rho_t_res_reg2 * sd_res_flood_reg_2_3  * mean(log(wy_order)) / sd(log(wy_order)))
  #cc1_LOGT[j] = rho_t_res_reg2 * sd_res_flood_reg_2_3 / sd(log(wy_order))
  #res_reg2_LOGT_fit <- cc0_LOGT + cc1_LOGT[j]* log(wy_order)
  #res_reg2_LOGT_var <- length(wy_order)/(length(wy_order)-2)*var(res_flood_reg_2_3-res_reg2_LOGT_fit)
  #se_cc0_LOGT <- sqrt(sum((res_flood_reg_2_3 - res_reg2_LOGT_fit)^2)/(length(wy_order) - 2)*(1/length(wy_order)+mean(log(wy_order))^2/sum((log(wy_order)-mean(log(wy_order)))^2))) 
  #se_cc1_LOGT <- sqrt(sum((res_flood_reg_2_3 - res_reg2_LOGT_fit)^2)/((length(wy_order) - 2)*sum((log(wy_order)-mean(log(wy_order)))^2)))
  #t_cc0_LOGT <- cc0_LOGT/se_cc0_LOGT
  #t_cc1_LOGT <- cc1_LOGT[j]/se_cc1_LOGT 
  #pval_cc0_LOGT <- 2*(pt(-abs(t_cc0_LOGT), df=length(wy_order) - 1)) 
  #pval_cc1_LOGT[j] = 2*(pt(-abs(t_cc1_LOGT), df=length(wy_order) - 1))  
  
  #flood_reg2.lm <- lm(res_flood_reg_2_3~wy_order)
  # Replace with models 5A,B,C,D
  #int_reg2 <- summary(flood_reg2.lm)$coefficients[1,1]
  #beta_reg2_max[j] = summary(flood_reg2.lm)$coefficients[2,1]
  #pval_int_reg2 <- summary(flood_reg2.lm)$coefficients[1,4]
  #pval_beta_reg2_max[j] = summary(flood_reg2.lm)$coefficients[2,4]
  
  # Identify the maximum slope among the three OLS models
  beta_reg2_max_time[j] = max(cc1_LIN_time[j],cc1_QUA_time[j],cc1_LOGT_time[j])
  beta_reg2_max_urb[j] = max(cc1_LIN_urb[j],cc1_QUA_urb[j],cc1_LOGT_urb[j])
  
  # Identify the minimum p-value for slope trends among the three OLS models
  pval_beta_reg2_max_time[j] = min(pval_cc1_LIN_time[j],pval_cc1_QUA_time[j],pval_cc1_LOGT_time[j])
  pval_beta_reg2_max_hc3_time[j] = min(pval_cc1_LIN_hc3_time[j],pval_cc1_QUA_hc3_time[j],pval_cc1_LOGT_hc3_time[j])
  
  pval_beta_reg2_max_urb[j] = min(pval_cc1_LIN_urb[j],pval_cc1_QUA_urb[j],pval_cc1_LOGT_urb[j])
  pval_beta_reg2_max_hc3_urb[j] = min(pval_cc1_LIN_hc3_urb[j],pval_cc1_QUA_hc3_urb[j],pval_cc1_LOGT_hc3_urb[j])
  
  # GLM model
  res_flood_reg1_time_2 <- (res_flood_reg1_time^2) #Compute squared residuals
  
  glm_results_time <- glm2(res_flood_reg1_time_2 ~ wy_order, family=Gamma(link=log),start=c(log(mean(res_flood_reg1_time_2)),0)) 
  
  cc1_GLM_time[j] = as.numeric(summary(glm_results_time)$coefficients[2,1])
  pval_cc1_GLM_time[j] = as.numeric(summary(glm_results_time)$coefficients[2,4])
  pval_cc1_GLM_hc3_time[j] = coeftest(glm_results_time,vcov=vcovHC(glm_results_time,"HC3"))[2,4]
  
  
  res_flood_reg1_urb_2 <- (res_flood_reg1_urb^2) #Compute squared residuals
  
  glm_results_urb <- glm2(res_flood_reg1_urb_2 ~ wy_order_urb, family=Gamma(link=log),start=c(log(mean(res_flood_reg1_urb_2)),0)) 
  
  cc1_GLM_urb[j] = as.numeric(summary(glm_results_urb)$coefficients[2,1])
  pval_cc1_GLM_urb[j] = as.numeric(summary(glm_results_urb)$coefficients[2,4])
  pval_cc1_GLM_hc3_urb[j] = coeftest(glm_results_urb,vcov=vcovHC(glm_results_urb,"HC3"))[2,4]

}


# Find table that can store diff data types
output_table_time <- data.frame(rec_id,
                         bhc_ts_urb_final$STAID_TEXT,
                         bhc_ts_urb_final$STANAME,
                         bhc_ts_urb_final$DRAIN_SQKM,
                         beta_reg1_time,
                       # pval_beta_reg1,
                         pval_beta_reg1_hc3_time,
                         adq_ppcc_ln2,
                         adq_ppcc_res_flood_reg1_time,
                         adq_dw_res_flood_reg1_time,
                         adq_res_flood_reg1_time,
                         cc1_LIN_time,
                         cc1_QUA_time,
                         cc1_LOGT_time,
                         beta_reg2_max_time,
                         #pval_cc1_LIN,
                         #pval_cc1_QUA,
                         #pval_cc1_LOGT,
                         pval_cc1_LIN_hc3_time,
                         pval_cc1_QUA_hc3_time,
                         pval_cc1_LOGT_hc3_time,
                         pval_beta_reg2_max_time,
                         pval_beta_reg2_max_hc3_time,
                         cc1_GLM_time,
                         pval_cc1_GLM_time,
                         pval_cc1_GLM_hc3_time)

output_table_urb <- data.frame(rec_id,
                              bhc_ts_urb_final$STAID_TEXT,
                              bhc_ts_urb_final$STANAME,
                              bhc_ts_urb_final$DRAIN_SQKM,
                              beta_reg1_urb,
                              # pval_beta_reg1,
                              pval_beta_reg1_hc3_urb,
                              adq_ppcc_ln2,
                              adq_ppcc_res_flood_reg1_urb,
                              adq_dw_res_flood_reg1_urb,
                              adq_res_flood_reg1_urb,
                              cc1_LIN_urb,
                              cc1_QUA_urb,
                              cc1_LOGT_urb,
                              beta_reg2_max_urb,
                              #pval_cc1_LIN,
                              #pval_cc1_QUA,
                              #pval_cc1_LOGT,
                              pval_cc1_LIN_hc3_urb,
                              pval_cc1_QUA_hc3_urb,
                              pval_cc1_LOGT_hc3_urb,
                              pval_beta_reg2_max_urb,
                              pval_beta_reg2_max_hc3_urb,
                              cc1_GLM_urb,
                              pval_cc1_GLM_urb,
                              pval_cc1_GLM_hc3_urb)

# Set p-value threshold for significance
pval_thresh <- 0.05

# Extract subset of stations with normally distributed and serially independent residuals (4751/8192)  
output_adq_time <- output_table_time[which(output_table_time$adq_res_flood_reg1_time==1),]   
output_adq_urb <- output_table_urb[which(output_table_urb$adq_res_flood_reg1_urb==1),]   

# Find increasing trends in the mean
output_mplus_time <- output_adq_time[which(output_adq_time$beta_reg1_time>0),] 
output_mplus_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb>0),] 

# Find statistically significant (p <0.05) increasing trends in the mean
output_mplus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time > 0
                               & output_adq_time$pval_beta_reg1_hc3_time < pval_thresh),]
output_mplus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb > 0
                                                 & output_adq_urb$pval_beta_reg1_hc3_urb < pval_thresh),]

# Find statistically significant (p <0.05) decreasing trends in the mean
output_mminus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time < 0
                                     & output_adq_time$pval_beta_reg1_hc3_time < pval_thresh),]
output_mminus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb < 0
                                                  & output_adq_urb$pval_beta_reg1_hc3_urb < pval_thresh),]

# Find statistically significant (p < 0.05) increasing trends in the variance
output_vplus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg2_max_time > 0
                                      & output_adq_time$pval_beta_reg2_max_hc3_time < pval_thresh),]
output_vplus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg2_max_urb > 0
                                      & output_adq_urb$pval_beta_reg2_max_hc3_urb < pval_thresh),]

# Find statistically significant (p < 0.05) decreasing trends in the variance
output_vminus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg2_max_time < 0
                                      & output_adq_time$pval_beta_reg2_max_hc3_time < pval_thresh),]

output_vminus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg2_max_urb < 0
                                      & output_adq_urb$pval_beta_reg2_max_hc3_urb < pval_thresh),]

# Find increasing trends in BOTH the mean and variance

output_mplus_sig05_vplus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time > 0
                                           & output_adq_time$pval_beta_reg1_hc3_time < pval_thresh
                                           & output_adq_time$beta_reg2_max_time > 0
                                           & output_adq_time$pval_beta_reg2_max_hc3_time < pval_thresh),]

output_mplus_sig05_vplus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb > 0
                                           & output_adq_urb$pval_beta_reg1_hc3_urb < pval_thresh
                                           & output_adq_urb$beta_reg2_max_urb > 0
                                           & output_adq_urb$pval_beta_reg2_max_hc3_urb < pval_thresh),]

rownames_df_time <- as.numeric(rownames(output_mplus_sig05_vplus_sig05_time))
for (k in 1:length(rownames_df_time)){
plot(wy_station[[rownames_df_time[k]]],log(peaks[[rownames_df_time[k]]]),xlab="Water Year",ylab="ln(Annual Peak Flow, cfs)",
     main=c(as.character(bhc_ts_urb_final$STANAME[rownames_df_time[k]]),as.character(bhc_ts_urb_final$STAID_TEXT[rownames_df_time[k]]),
            "DA",as.numeric(as.character(bhc_ts_urb_final$DRAIN_SQKM[rownames_df_time[k]]))),cex.main=0.8)
}

rownames_df_urb <- as.numeric(rownames(output_mplus_sig05_vplus_sig05_urb))
for (k in 1:length(rownames_df_urb)){
  plot(bhc_ts_TIA_vec_station[[rownames_df_urb[k]]],log(peaks_urb[[rownames_df_urb[k]]]),xlab="Total Impervious Area (%)",ylab="ln(Annual Peak Flow, cfs)",
       main=c(as.character(bhc_ts_urb_final$STANAME[rownames_df_urb[k]]),as.character(bhc_ts_urb_final$STAID_TEXT[rownames_df_urb[k]]),
              "DA",as.numeric(bhc_ts_urb_final$DRAIN_SQKM[rownames_df_urb[k]])),cex.main=0.8)
}


# Find increasing trends in the mean with decrease trends in the variance

output_mplus_sig05_vminus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time > 0
                                                 & output_adq_time$pval_beta_reg1_hc3_time < pval_thresh
                                                 & output_adq_time$beta_reg2_max_time < 0
                                                 & output_adq_time$pval_beta_reg2_max_hc3_time < pval_thresh),]

rownames_df_time <- as.numeric(rownames(output_mplus_sig05_vminus_sig05_time))
for (k in 1:length(rownames_df_time)){
plot(wy_station[[rownames_df_time[k]]],log(peaks[[rownames_df_time[k]]]),xlab="Water Year",ylab="ln(Annual Peak Flow, cfs)",
       main=c(as.character(bhc_ts_urb_final$STANAME[rownames_df_time[k]]),as.character(bhc_ts_urb_final$STAID_TEXT[rownames_df_time[k]]),
              "DA",as.numeric(as.character(bhc_ts_urb_final$DRAIN_SQKM[rownames_df_time[k]]))),cex.main=0.8)
}

output_mplus_sig05_vminus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb > 0
                                                              & output_adq_urb$pval_beta_reg1_hc3_urb < pval_thresh
                                                              & output_adq_urb$beta_reg2_max_urb < 0
                                                              & output_adq_urb$pval_beta_reg2_max_hc3_urb < pval_thresh),]

rownames_df_urb <- as.numeric(rownames(output_mplus_sig05_vminus_sig05_urb))
for (k in 1:length(rownames_df_urb)){
  plot(bhc_ts_TIA_vec_station[[rownames_df_urb[k]]],log(peaks_urb[[rownames_df_urb[k]]]),xlab="Total Impervious Area (%)",ylab="ln(Annual Peak Flow, cfs)",
       main=c(as.character(bhc_ts_urb_final$STANAME[rownames_df_urb[k]]),as.character(bhc_ts_urb_final$STAID_TEXT[rownames_df_urb[k]]),
              "DA",as.numeric(as.character(bhc_ts_urb_final$DRAIN_SQKM[rownames_df_urb[k]]))),cex.main=0.8)
}



# Negative trend in mean and positive trend in variance

output_mminus_sig05_vplus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time  < 0
                                                 & output_adq_time$pval_beta_reg1_hc3_time < pval_thresh
                                                 & output_adq_time$beta_reg2_max_time > 0
                                                 & output_adq_time$pval_beta_reg2_max_hc3_time < pval_thresh),]


output_mminus_sig05_vplus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb  < 0
                                                               & output_adq_urb$pval_beta_reg1_hc3_urb < pval_thresh
                                                               & output_adq_urb$beta_reg2_max_urb > 0
                                                               & output_adq_urb$pval_beta_reg2_max_hc3_urb < pval_thresh),]

# Negative trends in both the mean and variance
output_mminus_sig05_vminus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time < 0
                                                  & output_adq_time$pval_beta_reg1_hc3_time < pval_thresh
                                                  & output_adq_time$beta_reg2_max_time < 0
                                                  & output_adq_time$pval_beta_reg2_max_hc3_time < pval_thresh),]


output_mminus_sig05_vminus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb < 0
                                                               & output_adq_urb$pval_beta_reg1_hc3_urb < pval_thresh
                                                               & output_adq_urb$beta_reg2_max_urb < 0
                                                               & output_adq_urb$pval_beta_reg2_max_hc3_urb < pval_thresh),]



# No trend in mean, increase in variance
output_mnone_sig05_vplus_sig05_time <- output_adq_time [which(output_adq_time$pval_beta_reg1_hc3_time  > pval_thresh
                                                   & output_adq_time$beta_reg2_max_time  > 0
                                                   & output_adq_time$pval_beta_reg2_max_hc3_time  < pval_thresh),]

output_mnone_sig05_vplus_sig05_urb <- output_adq_urb [which(output_adq_urb$pval_beta_reg1_hc3_urb  > pval_thresh
                                                              & output_adq_urb$beta_reg2_max_urb  > 0
                                                              & output_adq_urb$pval_beta_reg2_max_hc3_urb  < pval_thresh),]




# COMPUTE NUMBER OF TRENDS OF EACH TYPE

# Linear - time
output_lin_mminus_sig05_vminus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time<0
                                                         & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh
                                                         & output_adq_time$cc1_LIN_time<0
                                                         & output_adq_time$pval_cc1_LIN_hc3_time<pval_thresh),]

output_lin_mminus_sig05_vplus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time<0
                                                         & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh
                                                         & output_adq_time$cc1_LIN_time>0
                                                         & output_adq_time$pval_cc1_LIN_hc3_time<pval_thresh),]

output_lin_mplus_sig05_vminus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time>0
                                                         & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh
                                                         & output_adq_time$cc1_LIN_time<0
                                                         & output_adq_time$pval_cc1_LIN_hc3_time<pval_thresh),]

output_lin_mplus_sig05_vplus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time>0
                                                        & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh
                                                        & output_adq_time$cc1_LIN_time>0
                                                        & output_adq_time$pval_cc1_LIN_hc3_time<pval_thresh),]

output_lin_mminus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time<0
                                            & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh),]

output_lin_mplus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1>0
                                                       & output_adq_time$pval_beta_reg1_hc3<pval_thresh),]

output_lin_vminus_sig05_time <- output_adq_time[which(output_adq_time$cc1_LIN_time<0
                                          & output_adq_time$pval_cc1_LIN_hc3_time<pval_thresh),]

output_lin_vplus_sig05_time <- output_adq_time[which(output_adq_time$cc1_LIN_time>0
                                            & output_adq_time$pval_cc1_LIN_hc3_time<pval_thresh),]


# Linear - urban
output_lin_mminus_sig05_vminus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb<0
                                                                   & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh
                                                                   & output_adq_urb$cc1_LIN_urb<0
                                                                   & output_adq_urb$pval_cc1_LIN_hc3_urb<pval_thresh),]

output_lin_mminus_sig05_vplus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb<0
                                                             & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh
                                                             & output_adq_urb$cc1_LIN_urb>0
                                                             & output_adq_urb$pval_cc1_LIN_hc3_urb<pval_thresh),]

output_lin_mplus_sig05_vminus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb>0
                                                                  & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh
                                                                  & output_adq_urb$cc1_LIN_urb<0
                                                                  & output_adq_urb$pval_cc1_LIN_hc3_urb<pval_thresh),]

output_lin_mplus_sig05_vplus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb>0
                                                                 & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh
                                                                 & output_adq_urb$cc1_LIN_urb>0
                                                                 & output_adq_urb$pval_cc1_LIN_hc3_urb<pval_thresh),]

output_lin_mminus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb<0
                                                      & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh),]

output_lin_mplus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1>0
                                                     & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh),]

output_lin_vminus_sig05_urb <- output_adq_urb[which(output_adq_urb$cc1_LIN_urb<0
                                                      & output_adq_urb$pval_cc1_LIN_hc3_urb<pval_thresh),]

output_lin_vplus_sig05_urb <- output_adq_urb[which(output_adq_urb$cc1_LIN_urb>0
                                                     & output_adq_urb$pval_cc1_LIN_hc3_urb<pval_thresh),]
# Quadratic - Time

output_qua_mminus_sig05_vminus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time<0
                                                         & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh
                                                         & output_adq_time$cc1_QUA_time<0
                                                         & output_adq_time$pval_cc1_QUA_hc3_time<pval_thresh),]

output_qua_mminus_sig05_vplus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time<0
                                                        & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh
                                                        & output_adq_time$cc1_QUA_time>0
                                                        & output_adq_time$pval_cc1_QUA_hc3_time<pval_thresh),]

output_qua_mplus_sig05_vminus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time>0
                                                        & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh
                                                        & output_adq_time$cc1_QUA_time<0
                                                        & output_adq_time$pval_cc1_QUA_hc3_time<pval_thresh),]

output_qua_mplus_sig05_vplus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time>0
                                                       & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh
                                                       & output_adq_time$cc1_QUA_time>0
                                                       & output_adq_time$pval_cc1_QUA_hc3_time<pval_thresh),]

output_qua_mminus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time<0
                                            & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh),]

output_qua_mplus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time>0
                                           & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh),]

output_qua_vminus_sig05_time <- output_adq_time[which(output_adq_time$cc1_QUA_time<0
                                            & output_adq_time$pval_cc1_QUA_hc3_time<pval_thresh),]

output_qua_vplus_sig05_time <- output_adq_time[which(output_adq_time$cc1_QUA_time>0
                                           & output_adq_time$pval_cc1_QUA_hc3_time<pval_thresh),]

# Quadratic - Urban

output_qua_mminus_sig05_vminus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb<0
                                                                   & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh
                                                                   & output_adq_urb$cc1_QUA_urb<0
                                                                   & output_adq_urb$pval_cc1_QUA_hc3_urb<pval_thresh),]

output_qua_mminus_sig05_vplus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb<0
                                                                  & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh
                                                                  & output_adq_urb$cc1_QUA_urb>0
                                                                  & output_adq_urb$pval_cc1_QUA_hc3_urb<pval_thresh),]

output_qua_mplus_sig05_vminus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb>0
                                                                  & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh
                                                                  & output_adq_urb$cc1_QUA_urb<0
                                                                  & output_adq_urb$pval_cc1_QUA_hc3_urb<pval_thresh),]

output_qua_mplus_sig05_vplus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb>0
                                                                 & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh
                                                                 & output_adq_urb$cc1_QUA_urb>0
                                                                 & output_adq_urb$pval_cc1_QUA_hc3_urb<pval_thresh),]

output_qua_mminus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb<0
                                                      & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh),]

output_qua_mplus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb>0
                                                     & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh),]

output_qua_vminus_sig05_urb <- output_adq_urb[which(output_adq_urb$cc1_QUA_urb<0
                                                      & output_adq_urb$pval_cc1_QUA_hc3_urb<pval_thresh),]

output_qua_vplus_sig05_urb <- output_adq_urb[which(output_adq_urb$cc1_QUA_urb>0
                                                     & output_adq_urb$pval_cc1_QUA_hc3_urb<pval_thresh),]

# Log(t) - Time

output_logt_mminus_sig05_vminus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time<0
                                                         & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh
                                                         & output_adq_time$cc1_LOGT_time<0
                                                         & output_adq_time$pval_cc1_LOGT_hc3_time<pval_thresh),]

output_logt_mminus_sig05_vplus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time<0
                                                        & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh
                                                        & output_adq_time$cc1_LOGT_time>0
                                                        & output_adq_time$pval_cc1_LOGT_hc3_time<pval_thresh),]

output_logt_mplus_sig05_vminus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time>0
                                                        & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh
                                                        & output_adq_time$cc1_LOGT_time<0
                                                        & output_adq_time$pval_cc1_LOGT_hc3_time<pval_thresh),]

output_logt_mplus_sig05_vplus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time>0
                                                       & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh
                                                       & output_adq_time$cc1_LOGT_time>0
                                                       & output_adq_time$pval_cc1_LOGT_hc3_time<pval_thresh),]

output_logt_mminus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time<0
                                             & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh),]

output_logt_mplus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time>0
                                           & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh),]

output_logt_vminus_sig05_time <- output_adq_time[which(output_adq_time$cc1_LOGT_time<0
                                            & output_adq_time$pval_cc1_LOGT_hc3_time<pval_thresh),]

output_logt_vplus_sig05_time <- output_adq_time[which(output_adq_time$cc1_LOGT_time>0
                                           & output_adq_time$pval_cc1_LOGT_hc3_time<pval_thresh),]


# Log(t) - Urban

output_logt_mminus_sig05_vminus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb<0
                                                                    & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh
                                                                    & output_adq_urb$cc1_LOGT_urb<0
                                                                    & output_adq_urb$pval_cc1_LOGT_hc3_urb<pval_thresh),]

output_logt_mminus_sig05_vplus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb<0
                                                                   & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh
                                                                   & output_adq_urb$cc1_LOGT_urb>0
                                                                   & output_adq_urb$pval_cc1_LOGT_hc3_urb<pval_thresh),]

output_logt_mplus_sig05_vminus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb>0
                                                                   & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh
                                                                   & output_adq_urb$cc1_LOGT_urb<0
                                                                   & output_adq_urb$pval_cc1_LOGT_hc3_urb<pval_thresh),]

output_logt_mplus_sig05_vplus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb>0
                                                                  & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh
                                                                  & output_adq_urb$cc1_LOGT_urb>0
                                                                  & output_adq_urb$pval_cc1_LOGT_hc3_urb<pval_thresh),]

output_logt_mminus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb<0
                                                       & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh),]

output_logt_mplus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb>0
                                                      & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh),]

output_logt_vminus_sig05_urb <- output_adq_urb[which(output_adq_urb$cc1_LOGT_urb<0
                                                       & output_adq_urb$pval_cc1_LOGT_hc3_urb<pval_thresh),]

output_logt_vplus_sig05_urb <- output_adq_urb[which(output_adq_urb$cc1_LOGT_urb>0
                                                      & output_adq_urb$pval_cc1_LOGT_hc3_urb<pval_thresh),]

# GLM_hc3 - urb

output_GLM_mminus_sig05_vminus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time<0
                                                          & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh
                                                          & output_adq_time$cc1_GLM_time<0
                                                          & output_adq_time$pval_cc1_GLM_hc3_time<pval_thresh),]

output_GLM_mminus_sig05_vplus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time<0
                                                         & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh
                                                         & output_adq_time$cc1_GLM_time>0
                                                         & output_adq_time$pval_cc1_GLM_hc3_time<pval_thresh),]

output_GLM_mplus_sig05_vminus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time>0
                                                         & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh
                                                         & output_adq_time$cc1_GLM_time<0
                                                         & output_adq_time$pval_cc1_GLM_hc3_time<pval_thresh),]

output_GLM_mplus_sig05_vplus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time>0
                                                        & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh
                                                        & output_adq_time$cc1_GLM_time>0
                                                        & output_adq_time$pval_cc1_GLM_hc3_time<pval_thresh),]

output_GLM_mminus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time<0
                                            & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh),]

output_GLM_mplus_sig05_time <- output_adq_time[which(output_adq_time$beta_reg1_time>0
                                            & output_adq_time$pval_beta_reg1_hc3_time<pval_thresh),]

output_GLM_vminus_sig05_time <- output_adq_time[which(output_adq_time$cc1_GLM_time<0
                                             & output_adq_time$pval_cc1_GLM_hc3_time<pval_thresh),]

output_GLM_vplus_sig05_time <- output_adq_time[which(output_adq_time$cc1_GLM_time>0
                                            & output_adq_time$pval_cc1_GLM_hc3_time<pval_thresh),]

# GLM_hc3 - urb

output_GLM_mminus_sig05_vminus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb<0
                                                                   & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh
                                                                   & output_adq_urb$cc1_GLM_urb<0
                                                                   & output_adq_urb$pval_cc1_GLM_hc3_urb<pval_thresh),]

output_GLM_mminus_sig05_vplus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb<0
                                                                  & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh
                                                                  & output_adq_urb$cc1_GLM_urb>0
                                                                  & output_adq_urb$pval_cc1_GLM_hc3_urb<pval_thresh),]

output_GLM_mplus_sig05_vminus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb>0
                                                                  & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh
                                                                  & output_adq_urb$cc1_GLM_urb<0
                                                                  & output_adq_urb$pval_cc1_GLM_hc3_urb<pval_thresh),]

output_GLM_mplus_sig05_vplus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb>0
                                                                 & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh
                                                                 & output_adq_urb$cc1_GLM_urb>0
                                                                 & output_adq_urb$pval_cc1_GLM_hc3_urb<pval_thresh),]

output_GLM_mminus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb<0
                                                      & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh),]

output_GLM_mplus_sig05_urb <- output_adq_urb[which(output_adq_urb$beta_reg1_urb>0
                                                     & output_adq_urb$pval_beta_reg1_hc3_urb<pval_thresh),]

output_GLM_vminus_sig05_urb <- output_adq_urb[which(output_adq_urb$cc1_GLM_urb<0
                                                      & output_adq_urb$pval_cc1_GLM_hc3_urb<pval_thresh),]

output_GLM_vplus_sig05_urb <- output_adq_urb[which(output_adq_urb$cc1_GLM_urb>0
                                                     & output_adq_urb$pval_cc1_GLM_hc3_urb<pval_thresh),]

# Count rows and make tables

nrows_lin_time <- c(nrow(output_lin_mminus_sig05_vminus_sig05_time),
               nrow(output_lin_mminus_sig05_vplus_sig05_time),
               nrow(output_lin_mplus_sig05_vminus_sig05_time),
               nrow(output_lin_mplus_sig05_vplus_sig05_time),
               nrow(output_lin_mminus_sig05_time),
               nrow(output_lin_mplus_sig05_time),
               nrow(output_lin_vminus_sig05_time),
               nrow(output_lin_vplus_sig05_time))  


nrows_lin_urb <- c(nrow(output_lin_mminus_sig05_vminus_sig05_urb),
                    nrow(output_lin_mminus_sig05_vplus_sig05_urb),
                    nrow(output_lin_mplus_sig05_vminus_sig05_urb),
                    nrow(output_lin_mplus_sig05_vplus_sig05_urb),
                    nrow(output_lin_mminus_sig05_urb),
                    nrow(output_lin_mplus_sig05_urb),
                    nrow(output_lin_vminus_sig05_urb),
                    nrow(output_lin_vplus_sig05_urb))  


trend_tbl_lin_time <- t(matrix(c(nrow(output_lin_mplus_sig05_vminus_sig05_time),
                   nrow(output_lin_mplus_sig05_time)-(nrow(output_lin_mplus_sig05_vminus_sig05_time)+nrow(output_lin_mplus_sig05_vplus_sig05_time)),
                   nrow(output_lin_mplus_sig05_vplus_sig05_time),
                   nrow(output_lin_vminus_sig05_time)-(nrow(output_lin_mminus_sig05_vminus_sig05_time)+nrow(output_lin_mplus_sig05_vminus_sig05_time)),
                   nrow(output_adq_time)-sum(nrows_lin_time),
                   nrow(output_lin_vplus_sig05_time)-(nrow(output_lin_mminus_sig05_vplus_sig05_time)+nrow(output_lin_mplus_sig05_vplus_sig05_time)),
                   nrow(output_lin_mminus_sig05_vminus_sig05_time),  
                   nrow(output_lin_mminus_sig05_time)-(nrow(output_lin_mminus_sig05_vminus_sig05_time)+nrow(output_lin_mminus_sig05_vplus_sig05_time)),
                   nrow(output_lin_mminus_sig05_vplus_sig05_time)),3,3))

trend_tbl_lin_urb <- t(matrix(c(nrow(output_lin_mplus_sig05_vminus_sig05_urb),
                                 nrow(output_lin_mplus_sig05_urb)-(nrow(output_lin_mplus_sig05_vminus_sig05_urb)+nrow(output_lin_mplus_sig05_vplus_sig05_urb)),
                                 nrow(output_lin_mplus_sig05_vplus_sig05_urb),
                                 nrow(output_lin_vminus_sig05_urb)-(nrow(output_lin_mminus_sig05_vminus_sig05_urb)+nrow(output_lin_mplus_sig05_vminus_sig05_urb)),
                                 nrow(output_adq_urb)-sum(nrows_lin_urb),
                                 nrow(output_lin_vplus_sig05_urb)-(nrow(output_lin_mminus_sig05_vplus_sig05_urb)+nrow(output_lin_mplus_sig05_vplus_sig05_urb)),
                                 nrow(output_lin_mminus_sig05_vminus_sig05_urb),  
                                 nrow(output_lin_mminus_sig05_urb)-(nrow(output_lin_mminus_sig05_vminus_sig05_urb)+nrow(output_lin_mminus_sig05_vplus_sig05_urb)),
                                 nrow(output_lin_mminus_sig05_vplus_sig05_urb)),3,3))


nrows_qua_time <- c(nrow(output_qua_mminus_sig05_vminus_sig05_time),
               nrow(output_qua_mminus_sig05_vplus_sig05_time),
               nrow(output_qua_mplus_sig05_vminus_sig05_time),
               nrow(output_qua_mplus_sig05_vplus_sig05_time),
               nrow(output_qua_mminus_sig05_time),
               nrow(output_qua_mplus_sig05_time),
               nrow(output_qua_vminus_sig05_time),
               nrow(output_qua_vplus_sig05_time))  

nrows_qua_urb <- c(nrow(output_qua_mminus_sig05_vminus_sig05_urb),
                    nrow(output_qua_mminus_sig05_vplus_sig05_urb),
                    nrow(output_qua_mplus_sig05_vminus_sig05_urb),
                    nrow(output_qua_mplus_sig05_vplus_sig05_urb),
                    nrow(output_qua_mminus_sig05_urb),
                    nrow(output_qua_mplus_sig05_urb),
                    nrow(output_qua_vminus_sig05_urb),
                    nrow(output_qua_vplus_sig05_urb))  

trend_tbl_qua_time <- t(matrix(c(nrow(output_qua_mplus_sig05_vminus_sig05_time),
                          nrow(output_qua_mplus_sig05_time)-(nrow(output_qua_mplus_sig05_vminus_sig05_time)+nrow(output_qua_mplus_sig05_vplus_sig05_time)),
                          nrow(output_qua_mplus_sig05_vplus_sig05_time),
                          nrow(output_qua_vminus_sig05_time)-(nrow(output_qua_mminus_sig05_vminus_sig05_time)+nrow(output_qua_mplus_sig05_vminus_sig05_time)),
                          nrow(output_adq_time)-sum(nrows_qua_time),
                          nrow(output_qua_vplus_sig05_time)-(nrow(output_qua_mminus_sig05_vplus_sig05_time)+nrow(output_qua_mplus_sig05_vplus_sig05_time)),
                          nrow(output_qua_mminus_sig05_vminus_sig05_time),  
                          nrow(output_qua_mminus_sig05_time)-(nrow(output_qua_mminus_sig05_vminus_sig05_time)+nrow(output_qua_mminus_sig05_vplus_sig05_time)),
                          nrow(output_qua_mminus_sig05_vplus_sig05_time)),3,3))

trend_tbl_qua_urb <- t(matrix(c(nrow(output_qua_mplus_sig05_vminus_sig05_urb),
                                 nrow(output_qua_mplus_sig05_urb)-(nrow(output_qua_mplus_sig05_vminus_sig05_urb)+nrow(output_qua_mplus_sig05_vplus_sig05_urb)),
                                 nrow(output_qua_mplus_sig05_vplus_sig05_urb),
                                 nrow(output_qua_vminus_sig05_urb)-(nrow(output_qua_mminus_sig05_vminus_sig05_urb)+nrow(output_qua_mplus_sig05_vminus_sig05_urb)),
                                 nrow(output_adq_urb)-sum(nrows_qua_urb),
                                 nrow(output_qua_vplus_sig05_urb)-(nrow(output_qua_mminus_sig05_vplus_sig05_urb)+nrow(output_qua_mplus_sig05_vplus_sig05_urb)),
                                 nrow(output_qua_mminus_sig05_vminus_sig05_urb),  
                                 nrow(output_qua_mminus_sig05_urb)-(nrow(output_qua_mminus_sig05_vminus_sig05_urb)+nrow(output_qua_mminus_sig05_vplus_sig05_urb)),
                                 nrow(output_qua_mminus_sig05_vplus_sig05_urb)),3,3))

nrows_logt_time <- c(nrow(output_logt_mminus_sig05_vminus_sig05_time),
               nrow(output_logt_mminus_sig05_vplus_sig05_time),
               nrow(output_logt_mplus_sig05_vminus_sig05_time),
               nrow(output_logt_mplus_sig05_vplus_sig05_time),
               nrow(output_logt_mminus_sig05_time),
               nrow(output_logt_mplus_sig05_time),
               nrow(output_logt_vminus_sig05_time),
               nrow(output_logt_vplus_sig05_time))

nrows_logt_urb <- c(nrow(output_logt_mminus_sig05_vminus_sig05_urb),
                     nrow(output_logt_mminus_sig05_vplus_sig05_urb),
                     nrow(output_logt_mplus_sig05_vminus_sig05_urb),
                     nrow(output_logt_mplus_sig05_vplus_sig05_urb),
                     nrow(output_logt_mminus_sig05_urb),
                     nrow(output_logt_mplus_sig05_urb),
                     nrow(output_logt_vminus_sig05_urb),
                     nrow(output_logt_vplus_sig05_urb))

trend_tbl_logt_time <- t(matrix(c(nrow(output_logt_mplus_sig05_vminus_sig05_time),
                          nrow(output_logt_mplus_sig05_time)-(nrow(output_logt_mplus_sig05_vminus_sig05_time)+nrow(output_logt_mplus_sig05_vplus_sig05_time)),
                          nrow(output_logt_mplus_sig05_vplus_sig05_time),
                          nrow(output_logt_vminus_sig05_time)-(nrow(output_logt_mminus_sig05_vminus_sig05_time)+nrow(output_logt_mplus_sig05_vminus_sig05_time)),
                          nrow(output_adq_time)-sum(nrows_logt_time),
                          nrow(output_logt_vplus_sig05_time)-(nrow(output_logt_mminus_sig05_vplus_sig05_time)+nrow(output_logt_mplus_sig05_vplus_sig05_time)),
                          nrow(output_logt_mminus_sig05_vminus_sig05_time),  
                          nrow(output_logt_mminus_sig05_time)-(nrow(output_logt_mminus_sig05_vminus_sig05_time)+nrow(output_logt_mminus_sig05_vplus_sig05_time)),
                          nrow(output_logt_mminus_sig05_vplus_sig05_time)),3,3))

trend_tbl_logt_urb <- t(matrix(c(nrow(output_logt_mplus_sig05_vminus_sig05_urb),
                                  nrow(output_logt_mplus_sig05_urb)-(nrow(output_logt_mplus_sig05_vminus_sig05_urb)+nrow(output_logt_mplus_sig05_vplus_sig05_urb)),
                                  nrow(output_logt_mplus_sig05_vplus_sig05_urb),
                                  nrow(output_logt_vminus_sig05_urb)-(nrow(output_logt_mminus_sig05_vminus_sig05_urb)+nrow(output_logt_mplus_sig05_vminus_sig05_urb)),
                                  nrow(output_adq_urb)-sum(nrows_logt_urb),
                                  nrow(output_logt_vplus_sig05_urb)-(nrow(output_logt_mminus_sig05_vplus_sig05_urb)+nrow(output_logt_mplus_sig05_vplus_sig05_urb)),
                                  nrow(output_logt_mminus_sig05_vminus_sig05_urb),  
                                  nrow(output_logt_mminus_sig05_urb)-(nrow(output_logt_mminus_sig05_vminus_sig05_urb)+nrow(output_logt_mminus_sig05_vplus_sig05_urb)),
                                  nrow(output_logt_mminus_sig05_vplus_sig05_urb)),3,3))

nrows_GLM_time <- c(nrow(output_GLM_mminus_sig05_vminus_sig05_time),
               nrow(output_GLM_mminus_sig05_vplus_sig05_time),
               nrow(output_GLM_mplus_sig05_vminus_sig05_time),
               nrow(output_GLM_mplus_sig05_vplus_sig05_time),
               nrow(output_GLM_mminus_sig05_time),
               nrow(output_GLM_mplus_sig05_time),
               nrow(output_GLM_vminus_sig05_time),
               nrow(output_GLM_vplus_sig05_time)) 

nrows_GLM_urb <- c(nrow(output_GLM_mminus_sig05_vminus_sig05_urb),
                    nrow(output_GLM_mminus_sig05_vplus_sig05_urb),
                    nrow(output_GLM_mplus_sig05_vminus_sig05_urb),
                    nrow(output_GLM_mplus_sig05_vplus_sig05_urb),
                    nrow(output_GLM_mminus_sig05_urb),
                    nrow(output_GLM_mplus_sig05_urb),
                    nrow(output_GLM_vminus_sig05_urb),
                    nrow(output_GLM_vplus_sig05_urb))  

trend_tbl_GLM_time <- t(matrix(c(nrow(output_GLM_mplus_sig05_vminus_sig05_time),
                          nrow(output_GLM_mplus_sig05_time)-(nrow(output_GLM_mplus_sig05_vminus_sig05_time)+nrow(output_GLM_mplus_sig05_vplus_sig05_time)),
                          nrow(output_GLM_mplus_sig05_vplus_sig05_time),
                          nrow(output_GLM_vminus_sig05_time)-(nrow(output_GLM_mminus_sig05_vminus_sig05_time)+nrow(output_GLM_mplus_sig05_vminus_sig05_time)),
                          nrow(output_adq_time)-sum(nrows_GLM_time),
                          nrow(output_GLM_vplus_sig05_time)-(nrow(output_GLM_mminus_sig05_vplus_sig05_time)+nrow(output_GLM_mplus_sig05_vplus_sig05_time)),
                          nrow(output_GLM_mminus_sig05_vminus_sig05_time),  
                          nrow(output_GLM_mminus_sig05_time)-(nrow(output_GLM_mminus_sig05_vminus_sig05_time)+nrow(output_GLM_mminus_sig05_vplus_sig05_time)),
                          nrow(output_GLM_mminus_sig05_vplus_sig05_time)),3,3))

trend_tbl_GLM_urb <- t(matrix(c(nrow(output_GLM_mplus_sig05_vminus_sig05_urb),
                                 nrow(output_GLM_mplus_sig05_urb)-(nrow(output_GLM_mplus_sig05_vminus_sig05_urb)+nrow(output_GLM_mplus_sig05_vplus_sig05_urb)),
                                 nrow(output_GLM_mplus_sig05_vplus_sig05_urb),
                                 nrow(output_GLM_vminus_sig05_urb)-(nrow(output_GLM_mminus_sig05_vminus_sig05_urb)+nrow(output_GLM_mplus_sig05_vminus_sig05_urb)),
                                 nrow(output_adq_urb)-sum(nrows_GLM_urb),
                                 nrow(output_GLM_vplus_sig05_urb)-(nrow(output_GLM_mminus_sig05_vplus_sig05_urb)+nrow(output_GLM_mplus_sig05_vplus_sig05_urb)),
                                 nrow(output_GLM_mminus_sig05_vminus_sig05_urb),  
                                 nrow(output_GLM_mminus_sig05_urb)-(nrow(output_GLM_mminus_sig05_vminus_sig05_urb)+nrow(output_GLM_mminus_sig05_vplus_sig05_urb)),
                                 nrow(output_GLM_mminus_sig05_vplus_sig05_urb)),3,3))

# COMPUTE OVERLAP BETWEEN TRENDS
length(unique(c(output_lin_mminus_sig05_vminus_sig05_time[,1],
       output_qua_mminus_sig05_vminus_sig05_time[,1],
       output_logt_mminus_sig05_vminus_sig05_time[,1],
       output_GLM_mminus_sig05_vminus_sig05_time[,1])))  

length(unique(c(output_lin_mminus_sig05_vminus_sig05_urb[,1],
                output_qua_mminus_sig05_vminus_sig05_urb[,1],
                output_logt_mminus_sig05_vminus_sig05_urb[,1],
                output_GLM_mminus_sig05_vminus_sig05_urb[,1])))  

length(unique(c(output_lin_mminus_sig05_vplus_sig05_time[,1],
         output_qua_mminus_sig05_vplus_sig05_time[,1],
         output_logt_mminus_sig05_vplus_sig05_time[,1],
         output_GLM_mminus_sig05_vplus_sig05_time[,1]))) 

length(unique(c(output_lin_mminus_sig05_vplus_sig05_urb[,1],
                output_qua_mminus_sig05_vplus_sig05_urb[,1],
                output_logt_mminus_sig05_vplus_sig05_urb[,1],
                output_GLM_mminus_sig05_vplus_sig05_urb[,1]))) 

length(unique(c(output_lin_mplus_sig05_vminus_sig05_time[,1],
         output_qua_mplus_sig05_vminus_sig05_time[,1],
         output_logt_mplus_sig05_vminus_sig05_time[,1],
         output_GLM_mplus_sig05_vminus_sig05_time[,1]))) 

length(unique(c(output_lin_mplus_sig05_vminus_sig05_urb[,1],
                output_qua_mplus_sig05_vminus_sig05_urb[,1],
                output_logt_mplus_sig05_vminus_sig05_urb[,1],
                output_GLM_mplus_sig05_vminus_sig05_urb[,1]))) 

length(unique(c(output_lin_mplus_sig05_vplus_sig05_time[,1],
         output_qua_mplus_sig05_vplus_sig05_time[,1],
         output_logt_mplus_sig05_vplus_sig05_time[,1],
         output_GLM_mplus_sig05_vplus_sig05_time[,1])))  

length(unique(c(output_lin_mplus_sig05_vplus_sig05_urb[,1],
                output_qua_mplus_sig05_vplus_sig05_urb[,1],
                output_logt_mplus_sig05_vplus_sig05_urb[,1],
                output_GLM_mplus_sig05_vplus_sig05_urb[,1])))