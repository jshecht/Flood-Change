# EXPORT SITES WITH TRENDS -------------------------------------------------------------

# Extract matches 

# Increasing trend in mean, decreasing trend in Cv
# FIX: Site info
index_test <- match(output_lin_mplus_sig05_vminus_sig05_urb$site_id,site_info$site_no)
site_info_mplus_vminus <- site_info[na.omit(index_test),]

# Increasing trend in mean, increasing trend in Cv
index_test_2 <- match(output_lin_mplus_sig05_vplus_sig05_clean$site_id,site_info$site_no)
site_info_mplus_vplus <- site_info[na.omit(index_test_2),]

# Combine subsets of interest
output_mplus_sig05_vtrend_sig05 <- rbind(output_lin_mplus_sig05_vminus_sig05_clean,output_lin_mplus_sig05_vplus_sig05_clean)

# Identify records 
index_test_3 <- match(output_mplus_sig05_vtrend_sig05$site_id,site_info$site_no)
row_index <- seq(1,nrow(site_info),1) # Creates row index
row_index[index_test_3] <- NA # Assign NA to rows with matches in index_test_3
site_info_no_mplus_vtrend <- site_info[na.omit(row_index),]

# MAPS -----------------------------------------------------
# Basic map
map("state",col="gray50",xlim=c(-125,-65),ylim=c(25,50),mar=rep(0.05,4))
points(site_info_no_mplus_vtrend[,5],site_info_no_mplus_vtrend[,4],pch=46,cex=0.1,col="gray80")
points(site_info_mplus_vplus[,5],site_info_mplus_vplus[,4],pch=2,cex=0.5)
points(site_info_mplus_vminus[,5],site_info_mplus_vminus[,4],pch=6,cex=0.5)
leg.txt <- c(expression(paste("+",mu,", +",Cv," ")),
             expression(paste("+",mu,", -",Cv," ")),
             "Other")    
legend(-120,29.5, leg.txt, cex=0.5, horiz = FALSE,pch=c(2,6,46),bty="n",col=c("black","black","gray80"))
abline(h=c(25,30,35,40,45),lwd=0.5,col="gray90")
abline(v=c(-70,-80,-90,-100,-110,-120),lwd=0.5,col="gray90")
box(which="plot",lty="solid")
deg_long <- c(120,110,100,90,80,70)
axis(1, at = -deg_long, labels = parse(text=paste(deg_long,"*degree ~ W")),cex.axis=0.6)
axis(2, at = seq(25,50,5),labels = parse(text=paste(seq(25,50,5),"*degree ~ N")),cex.axis=0.6)

#grid(nx = NULL, ny = nx, col = "lightgray", lty = "dotted",
 #    lwd = par("lwd"), equilogs = TRUE)
#north.arrow(xb=-70, yb=34, len=0.5, lab="N")
par(mar=c(5.1,4.1,4.1,2.1)) # Set margins back to default


# OTHER ATTEMPTS

#Albers Equal Area projection of CONUS
map("state",proj="albers",par=c(29.5,45.5))
points(site_info_mplus_vminus[,5],site_info_mplus_vminus[,4])



map("state",projection="albers",par=c(29.5,45.5),col="gray")
#plot(coordinates(site_info_mplus_vminus),add=TRUE)

site_info_mplus_vminus_coords <- spTransform(site_info_mplus_vminus, "+proj=longlat", "+proj=albers +lat0=29.5 +lat1=45.5")
plot(site_info_mplus_vminus_coords)

plot(site_info_mplus_vminus[,5],site_info_mplus_vminus[,4],pch=1,cex=0.8)
points(site_info_mplus_vplus[,5],site_info_mplus_vplus[,4],pch=2,cex=0.8)
points(site_info_no_mplus_vtrend[,5],site_info_no_mplus_vtrend[,4],pch=46,cex=0.3,col="gray")