#==============================================================================#
#                                                                              #
#                 Test GLM with negative binomial distribution                 #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(MASS)
library(DHARMa)
library(sf)
library(spdep)
Sys.setenv(LANG = "en")

# A. High-certainty carcasses --------------------------------------------------
table_glm <- read_csv(paste0("data/GLM/table_glm_10.csv")) %>%
  rename(RoadType = road_importance,
         DistanceAmenity = distance_amenity_km,
         DistanceWater = distance_water_km,
         Woodland = woodland) %>%
  mutate(RoadType = factor(RoadType, levels = c("minor_road", "major_road")))

# ~ 1. Fit the models -----------------------------------------------------------

glm_nb <- glm.nb(formula = nbr_carcasses ~ RoadType + DistanceAmenity
                 + DistanceWater + Woodland,
                 data = table_glm,
                 link = log)

summary(glm_nb)

# ~ 2. Check the assumptions of the model using DHARMA -------------------------

simulationOutput <- DHARMa::simulateResiduals(fittedModel = glm_nb, n = 500)
plot(simulationOutput)

testDispersion(glm_nb)
testZeroInflation(glm_nb)
testQuantiles(glm_nb)
testUniformity(glm_nb)
testResiduals(glm_nb)

plotResiduals(glm_nb, as.factor(table_glm$RoadType))
plotResiduals(glm_nb, table_glm$DistanceAmenity)
plotResiduals(glm_nb, table_glm$DistanceWater)
plotResiduals(glm_nb, table_glm$Woodland)


# Spatial autocorrelation ++++++++++++++++++++++++++++++++++++++++++++++++++++++
testSpatialAutocorrelation(simulationOutput = simulationOutput,
                           x = table_glm$long_midpoint, y = table_glm$lat_midpoint)
# Using DHARMa, I don't find spatial autocorrelation

# Test for spatial autocorrelation by hand
table_glm_resids <- add_column(table_glm, resids = residuals(glm_nb))

table_glm_resids_sf <- st_as_sf(table_glm_resids, coords = c("long_midpoint", "lat_midpoint"),
                                crs = 4326)
table_glm_resids_sf_proj <- st_transform(table_glm_resids_sf, crs = 32636)

moran_resids_df <- table_glm_resids_sf_proj %>%
  mutate(x_coord = st_coordinates(table_glm_resids_sf_proj)[,1],
         y_coord = st_coordinates(table_glm_resids_sf_proj)[,2]) %>%
  st_drop_geometry() %>%
  dplyr::select(x_coord,y_coord,resids)

mycoord <- matrix(cbind(moran_resids_df$x_coord,moran_resids_df$y_coord), ncol=2)
mycoord.nb <- knn2nb(knearneigh(mycoord, k = 4))
mycoord.lw <- nb2listw(mycoord.nb, style="W")
plot(mycoord.nb,mycoord) # Represents each point linked to the nearest 4 neighbors

moran.test(moran_resids_df$resids, mycoord.lw)
# Autocorrelation here. This may have to do with Aimara's way of using the 4 closest
# neighbors, etc.


# ~ 3. Deviance explained ------------------------------------------------------

car::Anova(glm_nb, type = "III")


