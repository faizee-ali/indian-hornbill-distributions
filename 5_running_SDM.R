#### train eval and test a model for Buceros bicornis ####
#with spatial folds cross validation
library(SDMtune)
library(terra)
library(sf)
library(dplyr)
library(rJava)
library(readr)
library(ENMeval)
library(ggplot2)
library(plotROC)
library(rasterVis)

env <- rast(list.files("env/", pattern = ".tif$", full.names = TRUE))
#correct names (remove objectionable caracters from predictor names)
names(env) <- make.names(names(env), unique = TRUE)

bg_bias <- st_read("data/biased_background_200m.gpkg", layer = "biased_background_200m")

#Ensure it's in WGS84
bg_bias <- st_transform(bg_bias, crs = crs(env))

# Extract coordinates
bg_bias_df <- cbind(st_coordinates(bg_bias), st_drop_geometry(bg_bias))
colnames(bg_bias_df)[1:2] <- c("longitude", "latitude")
pres <- read_tsv("data/combined_curated_all_sp_data.txt")
pres <- pres %>% filter(latitude != 'NA' & longitude != 'NA')
species <- c("Buceros bicornis", "Rhyticeros undulatus", "Anthracoceros albirostris",
             "Anorrhinus austeni", "Aceros nipalensis")
i <- 1
sp = species[i]
pres <- pres %>% filter(scientific.name == sp)  


# Create SWD object
data <- prepareSWD(species = sp,
                   p = pres[, c("longitude", "latitude")],
                   a = bg_bias_df[, c("longitude", "latitude")],
                   env = env, categorical = "LULC_raster_aligned_resampled_filled")

#### we will do spatial cross validation - no we will do k fold cross validation
## because our data is clustered and we do not want to predict the model to other spaces
# rather we want You're trying to understand how well the model fits your study area, not extrapolate beyond it.
#	Your predictions will be made within the same area and sampling effort as your training data.
# 	Because test data drawn nearby reflects local ecological relationships, not necessarily overfitting.

##### split into train val test, then split train into k folds!

data_split <- trainValTest(data, test = 0.2, val = 0.2, seed = 69, only_presence = FALSE)

train <- data_split[[1]]
test <- data_split[[3]]
val <- data_split[[2]]
random_folds <- randomFolds(data = train, 
                            seed = 69, 
                            only_presence = FALSE,
                            k = 3)

model_cv <- train(method = "Maxent", data = train, removeduplicates=TRUE, addsamplestobackground=TRUE, categorical = "LULC_raster_aligned_resampled_filled", fc = "lqph", folds = random_folds)
print(1)
#check AUC
print(cat("Training AUC for:", sp, auc(model_cv)))
#Training AUC:  0.9658108
print(cat("Testing AUC: ", auc(model_cv, test = TRUE)))


## tune across combinations of hyper parameters
h <- list(reg = seq(0.2, 1.8, 0.4), 
          fc = c("l", "q", "p", "h", "t", "lq", "lqp", "lqph", "lqpht"),
          iter = seq(200, 2000, 400))

exp_cv <- gridSearch(model_cv, 
                     hypers = h, 
                     metric = "auc", 
                     test = val)
print(2)
index <- which.max(exp_cv@results$test_AUC)
hyp<-as.data.frame(rbind(exp_cv@results[index, 1], exp_cv@results[index, 2], exp_cv@results[index, 3]))
print(index)
write_tsv(hyp, file = paste0(sp, "_tunedhypers.txt"))
# New train dataset containing only the selected variables
new_train <- exp_cv@models[[index]]@data 

# Merge only presence data
merged_data <- mergeSWD(new_train,
                        val,
                        only_presence = TRUE) 


final_model <- train("Maxent", 
                     data = merged_data, 
                     fc = exp_cv@results[index, 1], 
                     reg = exp_cv@results[index, 2],
                     iter = exp_cv@results[index, 3])
print(auc(final_model, 
          test = test))
saveRDS(final_model, file = paste0(sp,"_final_model.rds"))
readRDS("Buceros bicornis_final_model.rds")

# Plot ROC curve
png(paste0(sp,"ROC_plot.png"))
plotROC(final_model, 
        test = test)
dev.off()
print(3)

### plot the prediction
pred <- SDMtune::predict(final_model,
                         data = env,
                         type = "cloglog")

map <- predict(final_model,
               data = env,
               type = "cloglog")
print(4)
p <- plotPred(map,
              lt = "Habitat\nsuitability",
              colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))
ggsave(paste0(sp,"_pred.jpeg"), p, width = 9, height = 6, units = "in", dpi = 900)

imp <- varImp(final_model, permut = 10)

pv <- ggplot(imp, aes(x = reorder(Variable, Permutation_importance), y = Permutation_importance)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "Percentage Contribution of Variables",
         x = "Variable",
         y = "Contribution (%)") +
    theme_minimal()
print(5)
ggsave(paste0(sp,"variable_contribution_plot.jpeg"), pv, width = 8, height = 6, dpi = 300)
### plot suitability ###

# Reclassify suitability raster into categories
rcl <- matrix(c(
    -Inf, 0.1, 0,  # Exactly 0 (unsuitable)
    0.1,   0.4,   1,  # Low
    0.4, 0.7,   2,  # Medium
    0.7, Inf,   3   # High
), ncol = 3, byrow = TRUE)

classified_map <- classify(map, rcl)

# Convert raster to data.frame for ggplot
df <- as.data.frame(classified_map, xy = TRUE)
names(df)[3] <- "class"

# Assign class labels
df$class <- factor(df$class,
                   levels = c(0, 1, 2, 3),
                   labels = c("Unsuitable", "Low", "Medium", "High"))
print(6)
# Plot using ggplot2
p <- ggplot(df) +
    geom_raster(aes(x = x, y = y, fill = class)) +
    scale_fill_manual(values = c("gray80", "yellow", "orange", "darkgreen")) +
    coord_equal() +
    labs(title = "Graduated Suitability Map", fill = "Suitability") +
    theme_minimal()
print(7)
# Save plot
ggsave(paste0(sp,"graduated_suitability.jpeg"), p, width = 9, height = 6, units = "in", dpi = 600)

writeRaster(classified_map, 
            filename = paste0(sp,"graduated_suitability.tif"), 
            overwrite = TRUE, 
            datatype = "INT1U")  # 1-byte integer, suitable for class maps

i <- i+1

#### species 2 ####
sp = species[i]
print(sp)
sp = species[i]
pres <- read_tsv("data/combined_curated_all_sp_data.txt")
pres <- pres %>% filter(latitude != 'NA' & longitude != 'NA')

pres <- pres %>% filter(scientific.name == sp)  


# Create SWD object
data <- prepareSWD(species = sp,
                   p = pres[, c("longitude", "latitude")],
                   a = bg_bias_df[, c("longitude", "latitude")],
                   env = env, categorical = "LULC_raster_aligned_resampled_filled")

#### we will do spatial cross validation - no we will do k fold cross validation
## because our data is clustered and we do not want to predict the model to other spaces
# rather we want You're trying to understand how well the model fits your study area, not extrapolate beyond it.
#	Your predictions will be made within the same area and sampling effort as your training data.
# 	Because test data drawn nearby reflects local ecological relationships, not necessarily overfitting.

##### split into train val test, then split train into k folds!

data_split <- trainValTest(data, test = 0.2, val = 0.2, seed = 69, only_presence = FALSE)

train <- data_split[[1]]
test <- data_split[[3]]
val <- data_split[[2]]
random_folds <- randomFolds(data = train, 
                            seed = 69, 
                            only_presence = FALSE,
                            k = 3)

model_cv <- train(method = "Maxent", data = train, removeduplicates=TRUE, addsamplestobackground=TRUE, categorical = "LULC_raster_aligned_resampled_filled", fc = "lqph", folds = random_folds)
print(1)
#check AUC
print(cat("Training AUC for:", sp, auc(model_cv)))
#Training AUC:  0.9658108
print(cat("Testing AUC: ", auc(model_cv, test = TRUE)))


## tune across combinations of hyper parameters
h <- list(reg = seq(0.2, 1.8, 0.4), 
          fc = c("l", "q", "p", "h", "t", "lq", "lqp", "lqph", "lqpht"),
          iter = seq(200, 2000, 400))

exp_cv <- gridSearch(model_cv, 
                     hypers = h, 
                     metric = "auc", 
                     test = val)
print(2)
index <- which.max(exp_cv@results$test_AUC)
hyp<-as.data.frame(rbind(exp_cv@results[index, 1], exp_cv@results[index, 2], exp_cv@results[index, 3]))
print(index)
write_tsv(hyp, file = paste0(sp, "_tunedhypers.txt"))
# New train dataset containing only the selected variables
new_train <- exp_cv@models[[index]]@data 

# Merge only presence data
merged_data <- mergeSWD(new_train,
                        val,
                        only_presence = TRUE) 


final_model <- train("Maxent", 
                     data = merged_data, 
                     fc = exp_cv@results[index, 1], 
                     reg = exp_cv@results[index, 2],
                     iter = exp_cv@results[index, 3])
print(auc(final_model, 
          test = test))
saveRDS(final_model, file = paste0(sp,"_final_model.rds"))


# Plot ROC curve
png(paste0(sp,"ROC_plot.png"))
plotROC(final_model, 
        test = test)
dev.off()
print(3)

### plot the prediction
pred <- SDMtune::predict(final_model,
                         data = env,
                         type = "cloglog")

map <- predict(final_model,
               data = env,
               type = "cloglog")
print(4)
p <- plotPred(map,
              lt = "Habitat\nsuitability",
              colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))
ggsave(paste0(sp,"_pred.jpeg"), p, width = 9, height = 6, units = "in", dpi = 900)

imp <- varImp(final_model, permut = 10)

pv <- ggplot(imp, aes(x = reorder(Variable, Permutation_importance), y = Permutation_importance)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "Percentage Contribution of Variables",
         x = "Variable",
         y = "Contribution (%)") +
    theme_minimal()
print(5)
ggsave(paste0(sp,"variable_contribution_plot.jpeg"), pv, width = 8, height = 6, dpi = 300)
### plot suitability ###

# Reclassify suitability raster into categories
rcl <- matrix(c(
    -Inf, 0.1, 0,  # Exactly 0 (unsuitable)
    0.1,   0.4,   1,  # Low
    0.4, 0.7,   2,  # Medium
    0.7, Inf,   3   # High
), ncol = 3, byrow = TRUE)

classified_map <- classify(map, rcl)

# Convert raster to data.frame for ggplot
df <- as.data.frame(classified_map, xy = TRUE)
names(df)[3] <- "class"

# Assign class labels
df$class <- factor(df$class,
                   levels = c(0, 1, 2, 3),
                   labels = c("Unsuitable", "Low", "Medium", "High"))
print(6)
# Plot using ggplot2
p <- ggplot(df) +
    geom_raster(aes(x = x, y = y, fill = class)) +
    scale_fill_manual(values = c("gray80", "yellow", "orange", "darkgreen")) +
    coord_equal() +
    labs(title = "Graduated Suitability Map", fill = "Suitability") +
    theme_minimal()
print(7)
# Save plot
ggsave(paste0(sp,"graduated_suitability.jpeg"), p, width = 9, height = 6, units = "in", dpi = 600)

writeRaster(classified_map, 
            filename = paste0(sp,"graduated_suitability.tif"), 
            overwrite = TRUE, 
            datatype = "INT1U")  # 1-byte integer, suitable for class maps
i <- i+1
#### species 3 ####

sp = species[i]
print(sp)
sp = species[i]
pres <- read_tsv("data/combined_curated_all_sp_data.txt")
pres <- pres %>% filter(latitude != 'NA' & longitude != 'NA')

pres <- pres %>% filter(scientific.name == sp)  


# Create SWD object
data <- prepareSWD(species = sp,
                   p = pres[, c("longitude", "latitude")],
                   a = bg_bias_df[, c("longitude", "latitude")],
                   env = env, categorical = "LULC_raster_aligned_resampled_filled")

#### we will do spatial cross validation - no we will do k fold cross validation
## because our data is clustered and we do not want to predict the model to other spaces
# rather we want You're trying to understand how well the model fits your study area, not extrapolate beyond it.
#	Your predictions will be made within the same area and sampling effort as your training data.
# 	Because test data drawn nearby reflects local ecological relationships, not necessarily overfitting.

##### split into train val test, then split train into k folds!

data_split <- trainValTest(data, test = 0.2, val = 0.2, seed = 69, only_presence = FALSE)

train <- data_split[[1]]
test <- data_split[[3]]
val <- data_split[[2]]
random_folds <- randomFolds(data = train, 
                            seed = 69, 
                            only_presence = FALSE,
                            k = 3)

model_cv <- train(method = "Maxent", data = train, removeduplicates=TRUE, addsamplestobackground=TRUE, categorical = "LULC_raster_aligned_resampled_filled", fc = "lqph", folds = random_folds)
print(1)
#check AUC
print(cat("Training AUC for:", sp, auc(model_cv)))
#Training AUC:  0.9658108
print(cat("Testing AUC: ", auc(model_cv, test = TRUE)))


## tune across combinations of hyper parameters
h <- list(reg = seq(0.2, 1.8, 0.4), 
          fc = c("l", "q", "p", "h", "t", "lq", "lqp", "lqph", "lqpht"),
          iter = seq(200, 2000, 400))

exp_cv <- gridSearch(model_cv, 
                     hypers = h, 
                     metric = "auc", 
                     test = val)
print(2)
index <- which.max(exp_cv@results$test_AUC)
hyp<-as.data.frame(rbind(exp_cv@results[index, 1], exp_cv@results[index, 2], exp_cv@results[index, 3]))
print(index)
write_tsv(hyp, file = paste0(sp, "_tunedhypers.txt"))
# New train dataset containing only the selected variables
new_train <- exp_cv@models[[index]]@data 

# Merge only presence data
merged_data <- mergeSWD(new_train,
                        val,
                        only_presence = TRUE) 


final_model <- train("Maxent", 
                     data = merged_data, 
                     fc = exp_cv@results[index, 1], 
                     reg = exp_cv@results[index, 2],
                     iter = exp_cv@results[index, 3])
print(auc(final_model, 
          test = test))
saveRDS(final_model, file = paste0(sp,"_final_model.rds"))


# Plot ROC curve
png(paste0(sp,"ROC_plot.png"))
plotROC(final_model, 
        test = test)
dev.off()
print(3)

### plot the prediction
pred <- SDMtune::predict(final_model,
                         data = env,
                         type = "cloglog")

map <- predict(final_model,
               data = env,
               type = "cloglog")
print(4)
p <- plotPred(map,
              lt = "Habitat\nsuitability",
              colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))
ggsave(paste0(sp,"_pred.jpeg"), p, width = 9, height = 6, units = "in", dpi = 900)

imp <- varImp(final_model, permut = 10)

pv <- ggplot(imp, aes(x = reorder(Variable, Permutation_importance), y = Permutation_importance)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "Percentage Contribution of Variables",
         x = "Variable",
         y = "Contribution (%)") +
    theme_minimal()
print(5)
ggsave(paste0(sp,"variable_contribution_plot.jpeg"), pv, width = 8, height = 6, dpi = 300)
### plot suitability ###

# Reclassify suitability raster into categories
rcl <- matrix(c(
    -Inf, 0.1, 0,  # Exactly 0 (unsuitable)
    0.1,   0.4,   1,  # Low
    0.4, 0.7,   2,  # Medium
    0.7, Inf,   3   # High
), ncol = 3, byrow = TRUE)

classified_map <- classify(map, rcl)

# Convert raster to data.frame for ggplot
df <- as.data.frame(classified_map, xy = TRUE)
names(df)[3] <- "class"

# Assign class labels
df$class <- factor(df$class,
                   levels = c(0, 1, 2, 3),
                   labels = c("Unsuitable", "Low", "Medium", "High"))
print(6)
# Plot using ggplot2
p <- ggplot(df) +
    geom_raster(aes(x = x, y = y, fill = class)) +
    scale_fill_manual(values = c("gray80", "yellow", "orange", "darkgreen")) +
    coord_equal() +
    labs(title = "Graduated Suitability Map", fill = "Suitability") +
    theme_minimal()
print(7)
# Save plot
ggsave(paste0(sp,"graduated_suitability.jpeg"), p, width = 9, height = 6, units = "in", dpi = 600)

writeRaster(classified_map, 
            filename = paste0(sp,"graduated_suitability.tif"), 
            overwrite = TRUE, 
            datatype = "INT1U")  # 1-byte integer, suitable for class maps

i <- i+1
#### species 4 ####
sp = species[i]
sp = species[i]
pres <- read_tsv("data/combined_curated_all_sp_data.txt")
pres <- pres %>% filter(latitude != 'NA' & longitude != 'NA')

pres <- pres %>% filter(scientific.name == sp)  


# Create SWD object
data <- prepareSWD(species = sp,
                   p = pres[, c("longitude", "latitude")],
                   a = bg_bias_df[, c("longitude", "latitude")],
                   env = env, categorical = "LULC_raster_aligned_resampled_filled")

#### we will do spatial cross validation - no we will do k fold cross validation
## because our data is clustered and we do not want to predict the model to other spaces
# rather we want You're trying to understand how well the model fits your study area, not extrapolate beyond it.
#	Your predictions will be made within the same area and sampling effort as your training data.
# 	Because test data drawn nearby reflects local ecological relationships, not necessarily overfitting.

##### split into train val test, then split train into k folds!

data_split <- trainValTest(data, test = 0.2, val = 0.2, seed = 69, only_presence = FALSE)

train <- data_split[[1]]
test <- data_split[[3]]
val <- data_split[[2]]
random_folds <- randomFolds(data = train, 
                            seed = 69, 
                            only_presence = FALSE,
                            k = 3)

model_cv <- train(method = "Maxent", data = train, removeduplicates=TRUE, addsamplestobackground=TRUE, categorical = "LULC_raster_aligned_resampled_filled", fc = "lqph", folds = random_folds)
print(1)
#check AUC
print(cat("Training AUC for:", sp, auc(model_cv)))
#Training AUC:  0.9658108
print(cat("Testing AUC: ", auc(model_cv, test = TRUE)))


## tune across combinations of hyper parameters
h <- list(reg = seq(0.2, 1.8, 0.4), 
          fc = c("l", "q", "p", "h", "t", "lq", "lqp", "lqph", "lqpht"),
          iter = seq(200, 2000, 400))

exp_cv <- gridSearch(model_cv, 
                     hypers = h, 
                     metric = "auc", 
                     test = val)
print(2)
index <- which.max(exp_cv@results$test_AUC)
hyp<-as.data.frame(rbind(exp_cv@results[index, 1], exp_cv@results[index, 2], exp_cv@results[index, 3]))
print(index)
write_tsv(hyp, file = paste0(sp, "_tunedhypers.txt"))
# New train dataset containing only the selected variables
new_train <- exp_cv@models[[index]]@data 

# Merge only presence data
merged_data <- mergeSWD(new_train,
                        val,
                        only_presence = TRUE) 


final_model <- train("Maxent", 
                     data = merged_data, 
                     fc = exp_cv@results[index, 1], 
                     reg = exp_cv@results[index, 2],
                     iter = exp_cv@results[index, 3])
print(auc(final_model, 
          test = test))
saveRDS(final_model, file = paste0(sp,"_final_model.rds"))


# Plot ROC curve
png(paste0(sp,"ROC_plot.png"))
plotROC(final_model, 
        test = test)
dev.off()
print(3)

### plot the prediction
pred <- SDMtune::predict(final_model,
                         data = env,
                         type = "cloglog")

map <- predict(final_model,
               data = env,
               type = "cloglog")
print(4)
p <- plotPred(map,
              lt = "Habitat\nsuitability",
              colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))
ggsave(paste0(sp,"_pred.jpeg"), p, width = 9, height = 6, units = "in", dpi = 900)

imp <- varImp(final_model, permut = 10)

pv <- ggplot(imp, aes(x = reorder(Variable, Permutation_importance), y = Permutation_importance)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "Percentage Contribution of Variables",
         x = "Variable",
         y = "Contribution (%)") +
    theme_minimal()
print(5)
ggsave(paste0(sp,"variable_contribution_plot.jpeg"), pv, width = 8, height = 6, dpi = 300)
### plot suitability ###

# Reclassify suitability raster into categories
rcl <- matrix(c(
    -Inf, 0.1, 0,  # Exactly 0 (unsuitable)
    0.1,   0.4,   1,  # Low
    0.4, 0.7,   2,  # Medium
    0.7, Inf,   3   # High
), ncol = 3, byrow = TRUE)

classified_map <- classify(map, rcl)

# Convert raster to data.frame for ggplot
df <- as.data.frame(classified_map, xy = TRUE)
names(df)[3] <- "class"

# Assign class labels
df$class <- factor(df$class,
                   levels = c(0, 1, 2, 3),
                   labels = c("Unsuitable", "Low", "Medium", "High"))
print(6)
# Plot using ggplot2
p <- ggplot(df) +
    geom_raster(aes(x = x, y = y, fill = class)) +
    scale_fill_manual(values = c("gray80", "yellow", "orange", "darkgreen")) +
    coord_equal() +
    labs(title = "Graduated Suitability Map", fill = "Suitability") +
    theme_minimal()
print(7)
# Save plot
ggsave(paste0(sp,"graduated_suitability.jpeg"), p, width = 9, height = 6, units = "in", dpi = 600)

writeRaster(classified_map, 
            filename = paste0(sp,"graduated_suitability.tif"), 
            overwrite = TRUE, 
            datatype = "INT1U")  # 1-byte integer, suitable for class maps

i <- i+1

#### species 5 ####
sp = species[i]
pres <- read_tsv("data/combined_curated_all_sp_data.txt")
pres <- pres %>% filter(latitude != 'NA' & longitude != 'NA')

pres <- pres %>% filter(scientific.name == sp)  


# Create SWD object
data <- prepareSWD(species = sp,
                   p = pres[, c("longitude", "latitude")],
                   a = bg_bias_df[, c("longitude", "latitude")],
                   env = env, categorical = "LULC_raster_aligned_resampled_filled")

#### we will do spatial cross validation - no we will do k fold cross validation
## because our data is clustered and we do not want to predict the model to other spaces
# rather we want You're trying to understand how well the model fits your study area, not extrapolate beyond it.
#	Your predictions will be made within the same area and sampling effort as your training data.
# 	Because test data drawn nearby reflects local ecological relationships, not necessarily overfitting.

##### split into train val test, then split train into k folds!

data_split <- trainValTest(data, test = 0.2, val = 0.2, seed = 69, only_presence = FALSE)

train <- data_split[[1]]
test <- data_split[[3]]
val <- data_split[[2]]
random_folds <- randomFolds(data = train, 
                            seed = 69, 
                            only_presence = FALSE,
                            k = 3)

model_cv <- train(method = "Maxent", data = train, removeduplicates=TRUE, addsamplestobackground=TRUE, categorical = "LULC_raster_aligned_resampled_filled", fc = "lqph", folds = random_folds)
print(1)
#check AUC
print(cat("Training AUC for:", sp, auc(model_cv)))
#Training AUC:  0.9658108
print(cat("Testing AUC: ", auc(model_cv, test = TRUE)))


## tune across combinations of hyper parameters
h <- list(reg = seq(0.2, 1.8, 0.4), 
          fc = c("l", "q", "p", "h", "t", "lq", "lqp", "lqph", "lqpht"),
          iter = seq(200, 2000, 400))

exp_cv <- gridSearch(model_cv, 
                     hypers = h, 
                     metric = "auc", 
                     test = val)
print(2)
index <- which.max(exp_cv@results$test_AUC)
hyp<-as.data.frame(rbind(exp_cv@results[index, 1], exp_cv@results[index, 2], exp_cv@results[index, 3]))
print(index)
write_tsv(hyp, file = paste0(sp, "_tunedhypers.txt"))
# New train dataset containing only the selected variables
new_train <- exp_cv@models[[index]]@data 

# Merge only presence data
merged_data <- mergeSWD(new_train,
                        val,
                        only_presence = TRUE) 


final_model <- train("Maxent", 
                     data = merged_data, 
                     fc = exp_cv@results[index, 1], 
                     reg = exp_cv@results[index, 2],
                     iter = exp_cv@results[index, 3])
print(auc(final_model, 
          test = test))
saveRDS(final_model, file = paste0(sp,"_final_model.rds"))


# Plot ROC curve
png(paste0(sp,"ROC_plot.png"))
plotROC(final_model, 
        test = test)
dev.off()
print(3)

### plot the prediction
pred <- SDMtune::predict(final_model,
                         data = env,
                         type = "cloglog")

map <- predict(final_model,
               data = env,
               type = "cloglog")
print(4)
p <- plotPred(map,
              lt = "Habitat\nsuitability",
              colorramp = c("#2c7bb6", "#abd9e9", "#ffffbf", "#fdae61", "#d7191c"))
ggsave(paste0(sp,"_pred.jpeg"), p, width = 9, height = 6, units = "in", dpi = 900)

imp <- varImp(final_model, permut = 10)

pv <- ggplot(imp, aes(x = reorder(Variable, Permutation_importance), y = Permutation_importance)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(title = "Percentage Contribution of Variables",
         x = "Variable",
         y = "Contribution (%)") +
    theme_minimal()
print(5)
ggsave(paste0(sp,"variable_contribution_plot.jpeg"), pv, width = 8, height = 6, dpi = 300)
### plot suitability ###

# Reclassify suitability raster into categories
rcl <- matrix(c(
    -Inf, 0.1, 0,  # Exactly 0 (unsuitable)
    0.1,   0.4,   1,  # Low
    0.4, 0.7,   2,  # Medium
    0.7, Inf,   3   # High
), ncol = 3, byrow = TRUE)

classified_map <- classify(map, rcl)

# Convert raster to data.frame for ggplot
df <- as.data.frame(classified_map, xy = TRUE)
names(df)[3] <- "class"

# Assign class labels
df$class <- factor(df$class,
                   levels = c(0, 1, 2, 3),
                   labels = c("Unsuitable", "Low", "Medium", "High"))
print(6)
# Plot using ggplot2
p <- ggplot(df) +
    geom_raster(aes(x = x, y = y, fill = class)) +
    scale_fill_manual(values = c("gray80", "yellow", "orange", "darkgreen")) +
    coord_equal() +
    labs(title = "Graduated Suitability Map", fill = "Suitability") +
    theme_minimal()
print(7)
# Save plot
ggsave(paste0(sp,"graduated_suitability.jpeg"), p, width = 9, height = 6, units = "in", dpi = 600)

writeRaster(classified_map, 
            filename = paste0(sp,"graduated_suitability.tif"), 
            overwrite = TRUE, 
            datatype = "INT1U")  # 1-byte integer, suitable for class maps
