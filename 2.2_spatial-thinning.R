library(sf)
library(dplyr)
library(readr)

# Set grid size
gridsize <- 200  # in meters

# 1. Read Arunachal boundary and transform to UTM Zone 44N (EPSG:32644)
ap <- st_read("data/BUFFERED_StateBoundary_Arunachal_Pradesh.kml") %>%
    st_transform(32644)

# 2. Create 500 m grid across entire Arunachal
grid <- st_make_grid(ap, cellsize = gridsize, what = "polygons")
grid_sf <- st_sf(grid_id = 1:length(grid), geometry = grid)

# 3. Read eBird sampling data from TSV and convert to spatial points
ebird_sampling <- read_tsv("data/ebd_sampling_ap.txt")  # replace with actual path
ebird_sampling_sf <- st_as_sf(ebird_sampling,
                              coords = c("longitude", "latitude"),
                              crs = 4326)

# 4. Transform sampling points to the same CRS as grid (UTM Zone 44N)
ebird_sampling_proj <- st_transform(ebird_sampling_sf, crs = 32644)

# 5. Spatial join: assign each checklist to a grid cell
ebird_joined <- st_join(ebird_sampling_proj, grid_sf, join = st_intersects)

# 6. Create unique site IDs using coordinate pairs
ebird_joined$site_id <- paste0(round(st_coordinates(ebird_joined)[,1], 2), "_",
                               round(st_coordinates(ebird_joined)[,2], 2))
# 7. Count visits per site within each grid cell
most_visited <- ebird_joined %>%
    count(grid_id, site_id) %>%
    group_by(grid_id) %>%
    slice_max(n, with_ties = FALSE)
# 8 Filter original data to retain only those sites
background_bias <- ebird_joined %>%
    filter(site_id %in% most_visited$site_id & grid_id %in% most_visited$grid_id) %>%
    group_by(grid_id) %>%
    slice(1) %>%
    ungroup()

# 9. Save to file (optional)
st_write(background_bias, "biased_background_200m.gpkg")
background_bias_coords <- background_bias %>%
    mutate(
        longitude = st_coordinates(.)[, 1],
        latitude = st_coordinates(.)[, 2]
    ) %>%
    st_drop_geometry()

# Save as CSV
write_tsv(background_bias_coords, "output/biased_background_200m.txt")

save(ebird_sampling_sf, file = "ebird_sam.Rdata")
rm(ebird_sampling_sf)
####need to calibrate thinnign according to the average distance to another site 
#This needs- a mush ram heavier processor

library(sf)
library(geosphere)  # for accurate distances on Earth
library(dplyr)

load("ebird_sam.Rdata")
# Assume 'ebird_sampling_sf' is your unthinned sf object with WGS84 CRS
coords <- st_coordinates(ebird_sampling_sf)

# Compute pairwise great-circle distances (in meters)
dist_matrix <- distm(coords, fun = distHaversine)  # or distVincentyEllipsoid

# Set diagonal to NA so self-distances don't affect mean
diag(dist_matrix) <- NA

# Calculate average distance to all other points for each site
avg_distances <- rowMeans(dist_matrix, na.rm = TRUE)

# Add to the sf object
ebird_sampling_sf$avg_distance_m <- avg_distances



### try####
library(FNN)

coords <- st_coordinates(ebird_sampling_sf)
nn <- get.knn(coords, k = 2)  # 1st is self, 2nd is nearest neighbor
ebird_sampling_sf$nearest_neighbor_dist_m <- nn$nn.dist[, 2]

summary(ebird_sampling_sf$nearest_neighbor_dist_m)

mean_dist <- mean(ebird_sampling_sf$nearest_neighbor_dist_m)
median_dist <- median(ebird_sampling_sf$nearest_neighbor_dist_m)

cat("Mean NN distance:", round(mean_dist, 2), "m\n")
cat("Median NN distance:", round(median_dist, 2), "m\n")


### effort distance claculations ####

distances <- ebird_sampling$effort_distance_km
distances <- distances[!is.na(distances)]
avg_dist <- mean(distances)
med_dist <- median(distances)
rm(ebird_sampling, ebird_sampling_proj, ebird_sampling_sf)

# library
library(ggplot2)

# dataset:
data=data.frame(value=distances)
# basic histogram
p <- ggplot(data, aes(x=value)) + 
    geom_histogram(binwidth = 0.1) + xlim(0, 100)
p
