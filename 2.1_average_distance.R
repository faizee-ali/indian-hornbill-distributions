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
