library(raster)   
library(dplyr)
library(corrplot)
library(ggplot2)

setwd("D:/NCF/SDM_try-GH")
raster_files <- list.files("final-rasters/", pattern = ".tif$", full.names = TRUE)
env_stack <- stack(raster_files)
                   
nk <- layerStats(env_stack, stat = 'pearson', na.rm = TRUE)
cor_matrix <- nk$'pearson correlation coefficient'
plot.new()
# Plot
# Set output file path
outfile <- "correlation_plot.png"

# Open PNG graphics device (600 dpi, 6x6 inches)
png(filename = outfile, width = 6, height = 9, units = "in", res = 600)

# Generate the plot
corrplot(cor_matrix, method = "color", type = "lower",
         tl.col = "black", tl.cex = 0.3)

# Close the device
dev.off()


# Check if any device is open and close it safely
while (dev.cur() != 1) dev.off()

# Open new device (600 dpi PNG)
png("correlation_plot.png", width = 6, height = 9, units = "in", res = 600)

# Try the plot
corrplot(cor_matrix, method = "number", type = "upper", tl.cex = 0.3, number.cex = 0.3)

# Close device
dev.off()

write.csv(cor_matrix, "correlation_matrix.csv", row.names = TRUE)
#### filtering correlated variables####

## use th correlation matrix

# Read the CSV file (replace with your path)
cor_df <- read.csv("cor_matrenamed.csv", row.names = 1)

# Convert to a numeric matrix
cor_matrix <- as.matrix(cor_df)

# Optional: Check structure
str(cor_matrix)
# Use the same correlation matrix
dist_mat <- as.dist(1 - abs(cor_matrix))  # correlation distance
clust <- hclust(dist_mat, method = "average")
plot(clust)

library(corrplot)
corrplot(cor_matrix, method = "circle", type = "upper", order = "hclust", 
         tl.cex = 0.4, number.cex = 0.3, addrect = 4)

#### now do with variables selected by aparajta ####
#before that let us check correlation with asc format


#load cropped and aligned rasters

raster_files <- list.files("final-rasters/", pattern = ".tif$", full.names = TRUE)
env_stack <- stack(raster_files)
env_selected <- stack(env_stack$SRTM_Elevation_AP_BUFF_resampled_filled,
                      env_stack$GHS_POP_BUFF_resampled_filled,
                      env_stack$TreeCover_AP_resampled_filled,
                      env_stack$NVDI_MAY_AP_2015_resampled_filled,
                      env_stack$AP_CHELSA_bio1_1981.2010_V.2.1,
                      env_stack$AP_CHELSA_bio2_1981.2010_V.2.1,
                      env_stack$AP_CHELSA_bio4_1981.2010_V.2.1,
                      env_stack$AP_CHELSA_bio7_1981.2010_V.2.1,
                      env_stack$AP_CHELSA_bio12_1981.2010_V.2.1,
                      env_stack$AP_CHELSA_bio15_1981.2010_V.2.1,
                      env_stack$AP_CHELSA_bio17_1981.2010_V.2.1)


nk <- layerStats(env_selected, stat = 'pearson', na.rm = TRUE)
cor_matrix <- nk$'pearson correlation coefficient'
corrplot(cor_matrix, method = "circle", type = "upper", order = "hclust", 
         tl.cex = 0.5, number.cex = 0.5, addrect = 4)
####now my turn ####
env_selected <- stack(env_stack$GHS_POP_BUFF_resampled_filled,
                      env_stack$AP_CHELSA_bio17_1981.2010_V.2.1,
                      env_stack$AP_CHELSA_bio13_1981.2010_V.2.1,
                      env_stack$AP_CHELSA_bio4_1981.2010_V.2.1,
                      env_stack$AP_CHELSA_bio2_1981.2010_V.2.1,
                      env_stack$forest_height_BUFF_resampled_filled)

nk <- layerStats(env_selected, stat = 'pearson', na.rm = TRUE)
cor_matrix <- nk$'pearson correlation coefficient'
corrplot(cor_matrix, method = "circle", type = "upper", order = "hclust", 
         tl.cex = 0.5, number.cex = 0.5, addrect = 4)


