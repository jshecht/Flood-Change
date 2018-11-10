# List of stations to remove

# Increasing trends in mean and Cv----------------------------------------
output_lin_mplus_sig05_vplus_sig05_clean <- output_lin_mplus_sig05_vplus_sig05[-c(7,20,23),]

# Oquaga Creek at Deposit, NY (7)
# Kaskaskia Ditch and Bondville, IL (20)
# Calawah River near Forks, WA (23)

# Questionable ones to keep

# Wissahickon Creek at Fort Washington, PA (11)
# Little Patuxent River at Savage, MD (13)
# Scott Creek above Sylva, NC (16)
# Lamoille Creek near Lamoille, NV (22)

# Increasing trends in mean and decreasing trends in Cv-------------------
output_lin_mplus_sig05_vminus_sig05_clean <- output_lin_mplus_sig05_vminus_sig05[-c(7,18,22),]

# Brooker Creek near Lake Fern, FL (7)
# Platte River at Royalton, MN (18)
# Missouri River near Great Falls, MT (22)


# Questionable ones to keep

# Sprague Creek near Sprague, MB, Canada (17)
# South Platte River at Fort Lupton, CO (23)

# Export data
write.table(output_lin_mplus_sig05_vplus_sig05_clean,
            "c:/Users/Jory/Box Sync/USACE/R_Analysis/output_lin_mplus_sig05_vplus_sig05_clean.txt",
            sep="\t")
write.table(output_lin_mplus_sig05_vminus_sig05_clean,
            "c:/Users/Jory/Box Sync/USACE/R_Analysis/output_lin_mplus_sig05_vminus_sig05_clean.txt",
            sep="\t")

