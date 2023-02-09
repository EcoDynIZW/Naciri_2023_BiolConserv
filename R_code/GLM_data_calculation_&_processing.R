#==============================================================================#
#                                                                              #
#                  Linear model of spatial roadkill predictors                 #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(rgdal)
library(broom)
library(raster)
library(sp)
library(sf)
library(geosphere)
library(rgeos)

Sys.setenv(LANG = "en")

source(file = "R_code/functions_calculate_GLM_variables.R")

# A. Calculate response and dependent variables ================================

#==============================================================================#
# STEP 1: Import cropped road data --------------------------------------------

roads <- readOGR(dsn  =  "data/spatial",
                 layer = "road_network_for_GLM")
roads@data[["osm_id"]] <- seq(from = 1, to = length(roads@data[["osm_id"]]))

# Convert the shapefile into a data frame (because the road segmentation function
# takes data frames as arguments)
roads.df <- CreateDataFrameRoads(roads) %>%
  mutate(ID = as.factor(ID))


#==============================================================================#
# STEP 2: Segment the roads ----------------------------------------------------

roads_segmented_df <- SegmentationDf(roads.df,     # ~20 seconds
                                     segment.length = 2500)

write_csv(roads_segmented_df, "data/GLM/roads_segmented_df.csv")


#==============================================================================#
# STEP 3: Find the midpoint and measure the length of each segment ---------------------------------------

roads_segmented_df <- read_csv("data/GLM/roads_segmented_df.csv")

start <- Sys.time()
midpoint_and_length <- MidpointAndLengthAllSegments(roads_segmented_df) # ~30sec
end <- Sys.time() ; end-start

write_csv(midpoint_and_length, "data/GLM/midpoint_and_length.csv")

# Create table with all the info
table_glm_3.0 <- roads_segmented_df %>%
  distinct(ID_road_seg,
           road_importance) %>%
  mutate(nbr_carcasses = 0) # This column will be filled in step 6

# Add the midpoint and length of each segment
table_glm_3 <- table_glm_3.0 %>%
  left_join(x = .,
            y = midpoint_and_length,
            by = "ID_road_seg")

write_csv(table_glm_3, "data/GLM/table_glm_3.csv")


#==============================================================================#
# STEP 4: Delete segments that are too short -----------------------------------------

table_glm_3.2 <- read_csv("data/GLM/table_glm_3.csv")

table_glm_4 <- table_glm_3.2 %>%
  filter(length > 2000)

write_csv(table_glm_4, "data/GLM/table_glm_4.csv")


# Remove the coordinates of the points of the road segments that are shorter than
# 2000 meters from the road dataframe
roads_segmented_cropped_no_short_df <- roads_segmented_df %>%
  filter(ID_road_seg %in% table_glm_4$ID_road_seg)
write_csv(roads_segmented_cropped_no_short_df,
          "data/GLM/roads_segmented_cropped_no_short_df.csv")


#==============================================================================#
# STEP 5: Count the number of carcasses on each segment ------------------------

table_glm_4 <- read_csv("data/GLM/table_glm_4.csv")
roads_segmented_cropped_no_short_df <- read_csv("data/GLM/roads_segmented_cropped_no_short_df.csv")


# Import the file with the coordinates of the carcasses
hy_carcasses_with_GPS <- read_delim("data/01_hy_carcasses.csv",
                                             ";", escape_double = FALSE, trim_ws = TRUE) %>%
  filter(location_certainty_score >= 0.75)

start <- Sys.time()
closest_segment <- GetClosestRoadSegment(carcasses.df = hy_carcasses_with_GPS, # ~1min
                                         roads.df = roads_segmented_cropped_no_short_df)
end <- Sys.time() ; end-start

table_glm_5 <- AssignCarcassesToSegments(table_glm_4,   # 5 sec
                                         closest_segment)


sum(table_glm_5$nbr_carcasses)
write_csv(table_glm_5, "data/GLM/table_glm_5.csv")


#==============================================================================#
# STEP 6: Calculate distance to the closest water source -----------------------

### Load the data
rivers <- readOGR("data/spatial",
                      layer = "rivers_GLM")
waterbodies <- readOGR("data/spatial",
                           layer = "waterbodies_GLM")

table_glm_5 <- read_csv("data/GLM/table_glm_5.csv")

start <- Sys.time()
closest_water <- ShortestDistanceToRiverAndWaterbody(table_glm_5,     # ~ 12min
                                                     waterbodies,
                                                     rivers) %>%
  mutate_at(c("segment_ID", "ID_closest_water", "distance", "long_water",
              "lat_water"), unlist)
end <- Sys.time() ; end - start


write_csv(closest_water, "data/GLM/closest_water.csv")

### Add the columns to the table
table_glm_6 <- table_glm_5 %>%
  left_join(x = table_glm_5,
            y = closest_water[ , c("segment_ID", "distance")],
            by = "segment_ID") %>%
  mutate(distance_water_km = as.numeric(distance)/1000) %>%
  dplyr::select(-distance)

write_csv(table_glm_6, "data/GLM/table_glm_6.csv")



#==============================================================================#
# STEP 7: Calculate distance to the closest amenity ----------------------------

### Load the data
amenity_polygons <- readOGR("data/spatial",
                                layer = "amenity_polygons_GLM")
amenity_points <- readOGR("data/spatial",
                              layer = "amenity_points_GLM")

table_glm_6 <- read_csv("data/GLM/table_glm_6.csv")

start <- Sys.time()
closest_amenity <- ShortestDistanceAmenity(table_glm_6,     # ~ 6min
                                           amenity_polygons,
                                           amenity_points)
end <- Sys.time() ; end - start

write_csv(closest_amenity, "data/GLM/closest_amenity.csv")
# closest_amenity <- read_csv("data/GLM/closest_amenity.csv")

### Add the columns to the table
table_glm_7 <- table_glm_6 %>%
  left_join(x = table_glm_6,
            y = closest_amenity[ , c(1, 3)],
            by = "segment_ID") %>%
  mutate(distance_amenity_km = as.numeric(distance)/1000) %>%
  dplyr::select(-distance)

write_csv(table_glm_7, "data/GLM/table_glm_7.csv")


#==============================================================================#
# STEP 8: Extract land cover percentages --------------------------------------

table_glm_7 <- read_csv("data/GLM/table_glm_7.csv")
land_cover <- raster("data/spatial/serengeti_land_cover_Reed_2009.tif")



# ~ 1. Crop the land cover data to reduce the computational burden -------------

crop_square <- extent(665000, 750000, # Determined using the extent of the road network, with a margin
                      9660000, 9785000)
land_cover_cropped <- raster::crop(land_cover, crop_square)



# ~ 2. Run the function --------------------------------------------------------

roads_segmented_cropped_no_short_df <- read_csv("data/GLM/roads_segmented_cropped_no_short_df.csv")


start <- Sys.time()
land_cover_seg <- ExtractLandCover(roads_segmented_cropped_no_short_df, # ~ 4min
                                   land_cover_cropped,
                                   500) # 500m-radius buffer
percentage_land_cover <- PercentageLandCover(land_cover_seg)
end <- Sys.time() ; end - start

# Convert columns to numeric
percentage_land_cover <- percentage_land_cover %>%
  mutate_at(c("unclassified", "bare_ground","water","clouds", "sparse_grassland",
              "open_grassland", "dense_grassland", "closed_grassland",
              "sparse_shrubbed_grassland", "open_shrubbed_grassland",
              "dense_shrubbed_grassland", "closed_shrubbed_grassland",
              "open_treed_grassland", "dense_treed_grassland","closed_treed_grassland",
              "dense_shrubland","open_grassed_shrubland", "dense_grassed_shrubland",
              "closed_grassed_shrubland","open_treed_shrubland","dense_treed_shrubland",
              "closed_treed_shrubland", "dense_forest", "open_grassed_woodland",
              "dense_grassed_woodland", "closed_grassed_woodland", "dense_shrubbed_forest",
              "closed_shrubbed_forest"), as.numeric)

# Merge some land cover types
percentage_land_cover_processed <- percentage_land_cover %>%
  mutate(sparse_open_grassland = sparse_grassland + open_grassland + sparse_shrubbed_grassland +
           open_shrubbed_grassland + open_treed_grassland,

         dense_closed_grassland = dense_grassland + closed_grassland + dense_shrubbed_grassland +
           closed_shrubbed_grassland + dense_treed_grassland + closed_treed_grassland,

         sparse_open_shrubland = open_grassed_shrubland + open_treed_shrubland,

         dense_closed_shrubland = dense_shrubland + dense_grassed_shrubland + closed_grassed_shrubland
         + dense_treed_shrubland + closed_treed_shrubland,

         woodland = closed_grassed_woodland + open_grassed_woodland + dense_grassed_woodland,

         forest = dense_shrubbed_forest + closed_shrubbed_forest + dense_forest) %>%

  dplyr::select(segment_ID, sparse_open_grassland, dense_closed_grassland, sparse_open_shrubland, dense_closed_shrubland,
                woodland, forest, water, bare_ground) %>%
  # Combine all the grassland covers, and all the shrubland covers
  mutate(grassland = sparse_open_grassland + dense_closed_grassland,
         shrubland = sparse_open_shrubland + dense_closed_shrubland)

# Save
save(land_cover_seg,
     file = "data/GLM/extracted_cell_cover_values.RData")
write_csv(percentage_land_cover,
          file =  "data/GLM/percentage_land_cover.csv")
write_csv(percentage_land_cover_processed,
          file =  "data/GLM/percentage_land_cover.processed.csv")


table_glm_8 <- table_glm_7 %>%
  left_join(x = table_glm_7,
            y = percentage_land_cover_processed,
            by = "segment_ID")

write_csv(table_glm_8, file =  "data/GLM/table_glm_8.csv")

#==============================================================================#
# STEP 9: Remove ducplicated row -----------------------------------------------

# It is unclear why, but one road segment is duplicated. Perhaps in the map, there were
# two lines for this specific road portion. I need to remove it

table_glm_8 <- read_csv("data/GLM/table_glm_8.csv")

# Find the index of the duplicated row
which(duplicated(table_glm_8[, -1], incomparables = FALSE)) # Here I removed the first column
# (segment ID) because the first column is different for all rows, including the
# duplicated row. i.e. even the duplicated segment has two different IDs)

# Check that the value of e.g. longitude appear twice
which(table_glm_8$long_midpoint == table_glm_8$long_midpoint[which(duplicated(table_glm_8[, -1],
                                                                              incomparables = FALSE))])

# Remove the second time that the row appears (i.e. line 411)
table_glm_9 <- table_glm_8[-which(duplicated(table_glm_8[, -1], incomparables = FALSE)), ]

write_csv(table_glm_9, file =  "data/GLM/table_glm_9.csv")

#==============================================================================#
# STEP 10: Scale the variables -------------------------------------------------
table_glm_9 <- read_csv("data/GLM/table_glm_9.csv")


colnames(table_glm_9)
table_glm_10 <- table_glm_9 %>%
  mutate_at(c("length", "distance_water_km", "distance_amenity_km", "sparse_open_grassland",
              "dense_closed_grassland", "sparse_open_shrubland", "dense_closed_shrubland",
              "woodland", "forest", "grassland", "shrubland", "water", "bare_ground"), ~ (scale(.) %>% as.vector))

write_csv(table_glm_10, file =  "data/GLM/table_glm_10.csv")


# Save the mean and sd in a csv
table_glm_9 <- table_glm_9 %>%
  mutate_at(c("length", "distance_water_km", "distance_amenity_km", "sparse_open_grassland",
              "dense_closed_grassland", "sparse_open_shrubland", "dense_closed_shrubland",
              "woodland", "forest", "grassland", "shrubland", "water", "bare_ground"), as.numeric)

predictor_mean_sd <- cbind(colnames(table_glm_9[sapply(table_glm_9, is.numeric)]),
                           as.data.frame(sapply(table_glm_9[sapply(table_glm_9, is.numeric)],
                                                mean,
                                                na.rm = T)),
                           as.data.frame(sapply(table_glm_9[sapply(table_glm_9, is.numeric)],
                                                sd,
                                                na.rm = T)))
rownames(predictor_mean_sd) <- NULL
colnames(predictor_mean_sd) <- c("variable", "mean", "sd")

write_csv(predictor_mean_sd,
          file =  "data/GLM/predictor_mean_sd.csv")


#__________________________________________________________________________-----
# B. Visualize the data and covariates ===================================================

table_glm_10 <- read_csv("data/GLM/table_glm_10.csv")

# ~ 1. Response variable -------------------------------------------------------

ggplot(table_glm_10, aes(x = nbr_carcasses)) +
  geom_histogram()


# ~ 2. Covariates --------------------------------------------------------------

library(patchwork)

table_glm_9 <- read_csv("data/GLM/table_glm_9.csv")
table_glm_10 <- read_csv("data/GLM/table_glm_10.csv")

roads_segmented_df_covariates <- read_csv("data/GLM/roads_segmented_df.csv") %>%
  left_join(x = .,
            y = table_glm_9,
            by = c("ID_road_seg" = "segment_ID"))

# ~~~ a. Road type -------------------------------------------------------------

VS_amenity <- ggplot(table_glm_9, aes(x = road_importance, y = distance_amenity_km)) +
  geom_violin() +
  theme_bw()

VS_water <- ggplot(table_glm_9, aes(x = road_importance, y = distance_water_km)) +
  geom_violin() +
  theme_bw()

VS_woodland <- ggplot(table_glm_9, aes(x = road_importance, y = woodland)) +
  geom_violin() +
  theme_bw()

VS_amenity + VS_water + VS_woodland


# ~~~ b. Distance to amenities -------------------------------------------------

ggplot(table_glm_9, aes(x = distance_amenity_km)) +
  geom_histogram()

distance_amenity <- ggplot(roads_segmented_df_covariates,
       aes(x = long, y = lat, group = ID_road_seg, color = distance_amenity_km)) +
  geom_line(linewidth = 1)  +
  theme_bw() +
  viridis::scale_color_viridis() +
  labs(color = "DistanceAmenity (km)")

# ~~~ c. Distance to water -----------------------------------------------------

ggplot(table_glm_9, aes(x = distance_water_km)) +
  geom_histogram()

distance_water <- ggplot(roads_segmented_df_covariates,
       aes(x = long, y = lat, group = ID_road_seg, color = distance_water_km)) +
  geom_line(linewidth = 1)  +
  theme_bw() +
  viridis::scale_color_viridis() +
  labs(color = "DistanceWater (km)")

# ~~~ d. Woodland cover --------------------------------------------------------

ggplot(table_glm_9, aes(x = woodland)) +
  geom_histogram()

woodland <- ggplot(roads_segmented_df_covariates,
       aes(x = long, y = lat, group = ID_road_seg, color = woodland)) +
  geom_line(linewidth = 1)  +
  theme_bw() +
  viridis::scale_color_viridis() +
  labs(color = "Woodland cover")


library(patchwork)
distance_amenity + distance_water + woodland &
  theme(legend.position = "bottom")

ggplot(table_glm_9, aes(x = road_importance, y = grassland)) +
  geom_boxplot()

wilcox.test(formula = grassland ~ road_importance, data = table_glm_9, alternative = "two.sided", exact = TRUE)

ggplot(table_glm_9, aes(x = road_importance, y = woodland)) +
  geom_boxplot()

wilcox.test(formula = woodland ~ road_importance, data = table_glm_9, alternative = "two.sided", exact = TRUE)

ggplot(table_glm_9, aes(x = lat_midpoint, y = woodland)) +
  geom_point()








