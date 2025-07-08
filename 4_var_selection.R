#### variable selection ####
#install.packages("SDMtune")
library(SDMtune)
library(terra)
library(sf)
library(dplyr)
library(rJava)
library(readr)

# Environmental layers
env <- rast(list.files("final-rasters/", pattern = ".tif$", full.names = TRUE))

##Remove absoluterly useless variables
env <- env[[!names(env) %in% c("GHS_BUILT_S_BUFF_resampled_filled", 
                               "SRTM_Elevation_AP_BUFF_resampled_filled", 
                               "NVDI_NOV_AP_2015_resampled_filled",
                               "ET0_AP_V3_resampled_filled", 
                               "TreeCover_AP_resampled_filled",
                               "AP_CHELSA_bio1_1981-2010_V.2.1",
                               "AP_CHELSA_bio12_1981-2010_V.2.1")]]
# eBird presence and bias-corrected background
pres <- read_tsv("data/combined_curated_all_sp_data.txt")
pres <- pres %>% filter(latitude != 'NA' & longitude != 'NA')

pres <- pres %>% filter(scientific.name == "Anorrhinus austeni")  

bg_bias <- st_read("data/biased_background_200m.gpkg", layer = "biased_background_200m")

# Ensure it's in WGS84
bg_bias <- st_transform(bg_bias, crs = crs(env))

# Extract coordinates
bg_bias_df <- cbind(st_coordinates(bg_bias), st_drop_geometry(bg_bias))
colnames(bg_bias_df)[1:2] <- c("longitude", "latitude")
#correct names (remove objectionable caracters from predictor names)
names(env) <- make.names(names(env), unique = TRUE)
#PREPARE PRESENCE AND BACKGROUND LOCATIONS

p_coords <- select(pres, c("scientific.name", "longitude", "latitude"))
bg_coords <- select(bg_bias_df, c("checklist_id", "longitude", "latitude"))
# Create SWD object
data <- prepareSWD(species = "Buceros bicornis",
                   p = p_coords[, c("longitude", "latitude")],
                   a = bg_coords[, c("longitude", "latitude")],
                   env = env, categorical = "LULC_raster_aligned_resampled_filled")

#split into test and train datasets
datasets <- trainValTest(data,
                         test = 0.2)
train <- datasets[[1]]
test <- datasets[[2]]
# Basic model with all variables
model_full <- train(method = "Maxent", data = train, removeduplicates=TRUE, addsamplestobackground=TRUE, categorical = "LULC_raster_aligned_resampled_filled")


bg <- prepareSWD(species = "Anorrhinus austeni",
                 a = bg_coords[, c("longitude", "latitude")],
                 env = env)
# Filter with correlation threshold and importance metric

vs <- varSel(model_full,
             metric = "auc",
             bg4cor = bg,
             test = test,
             cor_th = 0.7,
             permut = 50,
             use_pc = TRUE)

#### with smaller dataset ####


save(model_full, file = "output/var_sel_model_trial.Rdata")


#### trial plot ####

#train a new model with only selected variables
#The variables AP_CHELSA_bio11_1981.2010_V.2.1, AP_CHELSA_bio16_1981.2010_V.2.1, AP_CHELSA_bio17_1981.2010_V.2.1, AP_CHELSA_bio18_1981.2010_V.2.1, AP_CHELSA_bio19_1981.2010_V.2.1, AP_CHELSA_bio2_1981.2010_V.2.1, AP_CHELSA_bio4_1981.2010_V.2.1, AP_CHELSA_bio5_1981.2010_V.2.1, AP_CHELSA_bio8_1981.2010_V.2.1, AP_CHELSA_bio9_1981.2010_V.2.1, and CHELSA_bio6_1981.2010_V.2.1 have been removed

env <- rast(list.files("final-rasters/", pattern = ".tif$", full.names = TRUE))

##Remove absoluterly useless variables
env <- env[[!names(env) %in% c("GHS_BUILT_S_BUFF_resampled_filled", 
                               "SRTM_Elevation_AP_BUFF_resampled_filled", 
                               "NVDI_NOV_AP_2015_resampled_filled",
                               "ET0_AP_V3_resampled_filled", 
                               "TreeCover_AP_resampled_filled",
                               "AP_CHELSA_bio1_1981-2010_V.2.1",
                               "AP_CHELSA_bio12_1981-2010_V.2.1",
                               "AP_CHELSA_bio11_1981.2010_V.2.1", 
                               "AP_CHELSA_bio16_1981.2010_V.2.1", 
                               "AP_CHELSA_bio17_1981.2010_V.2.1", 
                               "AP_CHELSA_bio18_1981.2010_V.2.1", 
                               "AP_CHELSA_bio19_1981.2010_V.2.1", 
                               "AP_CHELSA_bio2_1981.2010_V.2.1", 
                               "AP_CHELSA_bio4_1981.2010_V.2.1", 
                               "AP_CHELSA_bio5_1981.2010_V.2.1", 
                               "AP_CHELSA_bio8_1981.2010_V.2.1", 
                               "AP_CHELSA_bio9_1981.2010_V.2.1", 
                               "CHELSA_bio6_1981.2010_V.2.1")]]
# eBird presence and bias-corrected background
pres <- read_tsv("data/combined_curated_all_sp_data.txt")
pres <- pres %>% filter(latitude != 'NA' & longitude != 'NA')

pres <- pres %>% filter(scientific.name == "Buceros bicornis")  

bg_bias <- st_read("data/biased_background_200m.gpkg", layer = "biased_background_200m", categorical = "")

# Ensure it's in WGS84
bg_bias <- st_transform(bg_bias, crs = crs(env))

# Extract coordinates
bg_bias_df <- cbind(st_coordinates(bg_bias), st_drop_geometry(bg_bias))
colnames(bg_bias_df)[1:2] <- c("longitude", "latitude")
#correct names (remove objectionable caracters from predictor names)
names(env) <- make.names(names(env), unique = TRUE)
#PREPARE PRESENCE AND BACKGROUND LOCATIONS

p_coords <- select(pres, c("scientific.name", "longitude", "latitude"))
bg_coords <- select(bg_bias_df, c("checklist_id", "longitude", "latitude"))
# Create SWD object
data <- prepareSWD(species = "Buceros bicornis",
                   p = p_coords[, c("longitude", "latitude")],
                   a = bg_coords[, c("longitude", "latitude")],
                   env = env,
                   categorical = "LULC_raster_aligned_resampled_filled")

#split into test and train datasets
datasets <- trainValTest(data,
                         test = 0.2)
train <- datasets[[1]]
test <- datasets[[2]]
# Basic model with all variables
model_full <- train(method = "Maxent", data = train, removeduplicates=TRUE, addsamplestobackground=TRUE, categorical = "LULC_raster_aligned_resampled_filled", fc = "lqph")


# Get all variable names from the model
vs<- varImp(model_full)
auc(model_full)
auc(model_full, test = test)
# Loop through and plot each response
for (v in vars) {
    cat("Plotting:", v, "\n")
    plotResponse(model_full, var = v)
}
