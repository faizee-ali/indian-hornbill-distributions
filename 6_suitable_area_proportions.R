#### finding the area ####
library(SDMtune)
library(terra)
library(sf)
library(dplyr)
library(rJava)
library(readr)


### load models

bubi <- readRDS("Buceros bicornis_final_model.rds")

ths_bubi <- thresholds(bubi, 
                  type = "cloglog")

ths_bubi$species <- rep("bubi")

acni <- readRDS("Aceros nipalensis_final_model.rds")

ths_acni <- thresholds(acni, 
                       type = "cloglog")
ths_acni$species <- rep("acni")

anau <- readRDS("Anorrhinus austeni_final_model.rds")

ths_anau <- thresholds(anau, 
                       type = "cloglog")
ths_anau$species <- rep("anau")

anal <- readRDS("Anthracoceros albirostris_final_model.rds")

ths_anal <- thresholds(anal, 
                       type = "cloglog")
ths_anal$species <- rep("anal")

rhun <- readRDS("Rhyticeros undulatus_final_model.rds")

ths_rhun <- thresholds(rhun, 
                       type = "cloglog")

ths_rhun$species <- rep("rhun")

ths <- rbind(ths_acni, ths_anal, ths_anau, ths_bubi, ths_rhun)

write_tsv(ths, "threshold_area_props.txt")



### custom threshold ####

library(terra)

# Load or use your raster
# prediction_raster <- rast("path/to/prediction.tif")

# Define threshold
threshold <- 0.5

prediction_raster <- rast("aceros-nipalensis/Aceros nipalensisgraduated_suitability.tif")
# Create binary raster (1 = suitable, 0 = not suitable)
binary_map <- prediction_raster > threshold

# Calculate cell areas (in m² or km² depending on CRS)
area_raster <- cellSize(prediction_raster)

# Total area (sum of all non-NA cells)
total_area <- global(area_raster, fun = "sum", na.rm = TRUE)[1]

# Area above threshold
suitable_area <- global(area_raster * binary_map, fun = "sum", na.rm = TRUE)[1]

# Proportion
prop_suitable <- suitable_area / total_area

# Print result
cat("Proportion of area above threshold:", round(prop_suitable$sum, 4)*100, "\n")

