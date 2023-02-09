#==============================================================================#
#                                                                              #
#                   Generalized linear model of roadkills                      #
#                                                                              #
#==============================================================================#

library(tidyverse)
library(rgdal)
library(raster)
library(sp)
library(sf)
library(tmap)
library(geosphere)
library(rgeos)
library(effects)
library(glmmTMB)

Sys.setenv(LANG = "en")

# 1. Format road data -------------------------------------------------------------

CreateDataFrameRoads <- function(Shp) { # Converts the road shp into a df
  roads.df <- c()
  for (i in 1:length(Shp@lines)) {
    temp.i <- c()
    for (j in 1:length(Shp@lines[[i]]@Lines)) {

      temp.j <- data.frame(ID = Shp@data[["osm_id"]][i],
                           class = Shp@data[["fclass"]][i],
                           long = Shp@lines[[i]]@Lines[[j]]@coords[,1],
                           lat = Shp@lines[[i]]@Lines[[j]]@coords[,2],
                           type = ifelse(Shp@data[["fclass"]][i] %in% c("primary", "secondary"), "major_road", "minor_road"))
      temp.i <- rbind(temp.i, temp.j)
    }
    roads.df <- rbind(roads.df,temp.i)
  }
  return(roads.df)
}


# 2. Segment roads ----------------------------------------------------------------

CreateSegment <- function(coords, from, to) {
  distance <- 0
  coordsOut <- c()
  biggerThanFrom <- FALSE
  for (i in 1:(nrow(coords) - 1)) {
    d <- sqrt((coords[i, 1] - coords[i + 1, 1])^2 + (coords[i, 2] - coords[i + 1, 2])^2)   #distance between two successive points
    distance <- distance + d
    if (!biggerThanFrom && (distance > from)) {
      w <- 1 - (distance - from)/d
      x <- coords[i, 1] + w * (coords[i + 1, 1] - coords[i, 1])
      y <- coords[i, 2] + w * (coords[i + 1, 2] - coords[i, 2])
      coordsOut <- rbind(coordsOut, c(x, y))
      biggerThanFrom <- TRUE
    }
    if (biggerThanFrom) {
      if (distance > to) {
        w <- 1 - (distance - to)/d
        x <- coords[i, 1] + w * (coords[i + 1, 1] - coords[i, 1])
        y <- coords[i, 2] + w * (coords[i + 1, 2] - coords[i, 2])
        coordsOut <- rbind(coordsOut, c(x, y))
        break
      }
      coordsOut <- rbind(coordsOut, c(coords[i + 1, 1], coords[i + 1, 2]))
    }
  }
  return(coordsOut)
}


CreateSegments <- function(coords, length = 0, n.parts = 0) {
  stopifnot((length > 0 || n.parts > 0))
  # calculate total length of the road
  total_length <- 0
  for (i in 1:(nrow(coords) - 1)) {
    d <- sqrt((coords[i, 1] - coords[i + 1, 1])^2 +
                (coords[i, 2] - coords[i + 1, 2])^2) # sqrt( (xB-xA)^2 + (yB-yA)^2 )
    total_length <- total_length + d
  }

  # Create a vector with the position of segment ends/beginings along a straight line
  # and starting from 0 (the begining of the road)/
  # E.g. for a road that is 6 km long, and with segments 2.5km long, the following lines
  # would return (0, 2.5, 5, 6), with the last segment shorter than the others since
  # 6 is not a multiple of 2.5.
  if (length > 0) {
    stationing <- c(seq(from = 0, to = total_length, by = length),
                    total_length)
  } else {
    stationing <- c(seq(from = 0, to = total_length, length.out = n.parts),
                    total_length)
  }

  # Create the segments and store them in a list
  newlines <- list()
  for (i in 1:(length(stationing) - 1)) {
    newlines[[i]] <- CreateSegment(coords, stationing[i], stationing[i + 1])
  }
  return(newlines)
}








SegmentationDf <- function(roads.df, segment.length) {
  # Convert the real desired segment length into the weird length that the
  # functions understands
  segment.length.for.function <- -1.382e-06 + 8.989e-06 * segment.length

  # Create a list with the IDs
  list.ID <- c(levels(roads.df$ID))

  coords.all.roads.segmented.df <- c() # This dataframe will contain the info about
  # all the points of all the segments

  for (i in 1:length(list.ID)) {
    roads.df.i <- roads.df %>%          # Keep only the points that belong to road i
      filter(ID == list.ID[i])
    coords.road.i <- roads.df.i %>%     # Coords of the points that belong to road i
      dplyr::select(long,lat)
    coords.road.i.segmented <- CreateSegments(coords.road.i,   #segmentation de la route qui correspond au subset, renvoie une liste (=une route) de listes(=les segments) de coordonnées (
                                              length = segment.length.for.function)
    ### Merge the last segment with the punultimate segment if the last segment is shorter than half the desired length.
    # L.last.seg <- 0
    # if (length(coords.road.i.segmented) > 1) {
    #   for (k in 1:(length(coords.road.i.segmented[[length(coords.road.i.segmented)]])/2-1)) { # Calculate the length of the last segment (last member of the list of lists coords.road.i.segmented)
    #     d <- distGeo(coords.road.i.segmented[[length(coords.road.i.segmented)]][k,],
    #                        coords.road.i.segmented[[length(coords.road.i.segmented)]][k+1,])
    #     L.last.seg <- L.last.seg + d
    #   }
    #   if (L.last.seg < (real.segment.length/2)) {  # Merge the last two segments (i.e. the last segments
    #
    #     coords.road.i.segmented[[length(coords.road.i.segmented)-1]] <- rbind(coords.road.i.segmented[[length(coords.road.i.segmented)-1]],
    #                                                                             coords.road.i.segmented[[length(coords.road.i.segmented)]])
    #     coords.road.i.segmented <- coords.road.i.segmented[1:(length(coords.road.i.segmented)-1)]
    #   }
    # }

    # Create a data frame with all the desired info for the segments of road i
    for (j in 1:length(coords.road.i.segmented)) {
      #fait un dataframe avec une colonne ID de la route + une colonne ID du segment de la route et les coordonnées GPS
      npoint <- seq(1:(length(coords.road.i.segmented[[j]])/2))
      coords.road.i.segmented.df <- data.frame(ID_point = npoint,
                                               ID_seg = j,
                                               ID_road = i,
                                               ID_road_seg = paste0(i, "_", j),
                                               ID_road_seg_point = paste0(i,"_",j, "_", npoint),
                                               long = coords.road.i.segmented[[j]][,1],
                                               lat = coords.road.i.segmented[[j]][,2],
                                               osm_ID = list.ID[i],
                                               class = roads.df.i$class[j],
                                               road_importance = (ifelse(roads.df.i$class[j] %in% c("primary", "secondary"),
                                                                         "major_road", "minor_road")),
                                               color = ifelse((j %% 2) == 0, "black","red"))
      # Add the info about segments of road i to the dataframe that contains the info
      # about segments of all the roads.
      coords.all.roads.segmented.df <- rbind(coords.all.roads.segmented.df,
                                             coords.road.i.segmented.df)
    }
  }
  return(coords.all.roads.segmented.df)
}


# 3. Find the midpoint and measure the length of each segment ---------------------

# Function to calculate the midpoint and length of a single segment
MidpointAndLengthSegment <- function(coords) {
  L <- 0
  for (i in 1:(nrow(coords) - 1)) {
    d <- distGeo(coords[i,], coords[i+1,])  # distHarvesine calculates the shortest ditstance between
    # two points assuming a spherical earth.
    L <- L + d           # Here we sum the distances between
                         # each two successive points of a segment
                         # to get the total segment length
  }
  # Below, we sum up the distance between two successive points until we reach
  # half of the total length of the segment. As moqt often, the midpoint of the
  # segment is located between two points, we take only part of the length between
  # these two points.
  j <- 1
  sum.intervals <- 0
  interval <- distGeo(coords[j,], coords[j+1,])
  midpoint <- c()
  while ((sum.intervals + interval) < L/2) {
    sum.intervals <- sum.intervals + interval
    j <- j + 1
    interval <- distGeo(coords[j,], coords[j+1,])

  }
  proportion <- ((L/2) - sum.intervals)/(interval) # How much of the distance between the point just before
                                                   # the midpoint and the point just after separates the point
                                                   # just before the midpoint from the midpoint ?
  x <- as.numeric(coords[j, 1] + proportion * (coords[j + 1, 1] - coords[j, 1]))
  y <- as.numeric(coords[j, 2] + proportion * (coords[j + 1, 2] - coords[j, 2]))
  midpoint <- c(x, y)
  return(c(midpoint, L))
}

# MidpointAndLengthSegment <- function(coords) {
#   L <- 0
#   for (i in 1:(nrow(coords) - 1)) {
#     d <- distGeo(coords[i,], coords[i+1,])  # distHarvesine calculates the shortest ditstance between
#                                                   # two points assuming a spherical earth.
#     L <- L + d           # Here we sum the distances between
#                          # each two successive points of a segment
#                          # to get the total segment length
#   }
#   sum.intervals <- 0
#   j <- 0
#   continue <- TRUE
#   midpoint <- c()
#   while (continue) {
#     j = j + 1
#     interval <- distGeo(coords[j,], coords[j+1,])
#     if ((sum.intervals + interval) < L/2) {
#       sum.intervals <- sum.intervals + interval
#     } else {
#       proportion = ((L/2) - sum.intervals)/(interval)
#       x <- coords[j, 1] + proportion * (coords[j + 1, 1] - coords[j, 1])
#       y <- coords[j, 2] + proportion * (coords[j + 1, 2] - coords[j, 2])
#       midpoint <- c(x, y)
#       continue <- FALSE
#     }
#   }
#   return(c(midpoint, L))
# }



# Function to calculate the midpoint and length of multiple segments
MidpointAndLengthAllSegments <- function(roads.segmented.df) {
  new.columns <- c()
  list.seg.ID <- roads.segmented.df %>%
    distinct(ID_road_seg)
  for (i in 1:nrow(list.seg.ID)) {            # For each segment
    coords.seg.i <- roads.segmented.df %>%    # Get coordinates of segment i
      filter(ID_road_seg == list.seg.ID$ID_road_seg[i]) %>%
      dplyr::select(long, lat)
    new.row <- MidpointAndLengthSegment(coords.seg.i) # calculate the midpoint and length of segment i
    new.columns <- rbind(new.columns, new.row)
  }
  new.columns <- as.data.frame(new.columns)
  colnames(new.columns) <- c("long_midpoint", "lat_midpoint", "length")
  midpoint_length_df <- data.frame(ID_road_seg = list.seg.ID$ID_road_seg,
                                   long_midpoint = new.columns[, 1],
                                   lat_midpoint = new.columns[, 2],
                                   length = new.columns[, 3])
  return(midpoint_length_df)
}



# # 4. Delete roads that are too far from the center of the study area --------------

# This functions calculates the distance between the midpoint of each segment, and
# the center of the 65km-radius circle that delimitates our study area. It removes
# the segments that have a midpoint that is further than 65km from the center of the
# circle
DeleteSegmentFarFromCenter <- function(table.glm, coord.center, lim.dist) {
  distances.from.center <- c()
  for (i in 1:nrow(table.glm)) {
    D <- distGeo(coord.center,
                       table.glm[i, c("long_midpoint", "lat_midpoint")])
    distances.from.center <- rbind(distances.from.center, D)
  }
  table.glm <- as.data.frame(cbind(table.glm, distances.from.center))
  table.glm <- table.glm %>%
    filter(distances.from.center < lim.dist)
  return(table.glm)
}



# 5. Count the number of carcasses on each segment ------------------------

### Locate the closest segment to each carcass
GetClosestRoadSegment <- function(carcasses.df, roads.df) {
  dist.min <- c()
  for (i in 1:nrow(carcasses.df)) {
    roads.df.i <- roads.df %>%
      filter(abs(long - as.numeric(carcasses.df$long[i])) < 0.05,
             abs(lat - as.numeric(carcasses.df$lat[i])) < 0.05)
    dist.carcass.i <- c()
    d.min <- 1000000000
    ID <- carcasses.df[i, 1]
    for (j in 1:nrow(roads.df.i)) {
      d <- distGeo(carcasses.df[i , c('long', 'lat')],
                         roads.df.i[j, c('long', 'lat')])
      if (d < d.min) { # if the distance between carcass j and road i is smaller
                       # than the minimum  distance until now for carcass i
        ID_closest_seg <- roads.df.i$ID_road_seg[j]
        d.min <- d
        long <- roads.df.i$long[j]
        lat <- roads.df.i$lat[j]
        ID_road_seg_point <- roads.df.i$ID_road_seg_point[j]
      }
    }
    dist.min <- rbind(dist.min, c(ID, ID_closest_seg, d.min, long, lat, ID_road_seg_point))
  }
  dist.min <- as.data.frame(dist.min)
  colnames(dist.min) <- c("ID", "ID_closest_seg", "distance", "long", "lat", "ID_full")
  dist.min$long <- as.numeric(dist.min$long)
  dist.min$lat <- as.numeric(dist.min$lat)
  return(as.data.frame(dist.min))
}


### Assign each carcass to the closest segment
AssignCarcassesToSegments <- function(table.glm, closest.segment) {
  table.glm <- as.data.frame(table.glm)
  for (i in 1:nrow(closest.segment)) {
    for (j in 1:nrow(table.glm)) {
      if (as.factor(table.glm[j, 1]) %in% as.factor(closest.segment[i, 2])) {
        table.glm$nbr_carcasses[j] <- table.glm$nbr_carcasses[j] + 1
      }
    }
  }
  colnames(table.glm) <- c("segment_ID", "road_importance", "nbr_carcasses",
                           "long_midpoint", "lat_midpoint", "length") # "dist_from_center") # I removed this since I don't use the cirlce anymore
  return(table.glm)
}


# 7. Find the closest waterbody or river ------------------------------------------

ShortestDistanceToRiverAndWaterbody <- function(table.glm, waterbodies.shp, rivers.shp) {
  ### Create a dataframe with the coordinates, osm Id, names, class of waterbodies
  waterbodies.df <- tidy(waterbodies.shp) %>%
    dplyr::select(ID = id,
           long, lat,
           order)
  ### Create a dataframe with the coordinates, osm ID, names, class of rivers
  rivers.df <- tidy(rivers.shp) %>%
    dplyr::select(ID = id,
           long, lat,
           order)
  ### Fuse the waterbody df and the river df
  water.df <- rbind(waterbodies.df, rivers.df)
  ### Calculate the pair wise distances between each segment and all the points and keep the minimum for each road segment
  closest.water <- c()
  for (i in 1:nrow(table.glm)) {
    water.df.i <- water.df %>%                                # This is just to reduce the computational burden:
      filter(abs(long - table.glm$long_midpoint[i]) < 0.122,  # for each segment I only consider water points that
             abs(lat - table.glm$lat_midpoint[i]) < 0.122)    # are less than 0.122 degrees (latitude and longitude) away
    # from the midpoint. 0.164 degrees is a conservative value.
    dist.seg.i <- c()
    ID_segment <- table.glm$segment_ID[i]
    d.min <- 1000000000                    # Initialize the minimum distance to a water
    # source with a  purposefully inflated value
    for (j in 1:nrow(water.df.i)) {
      d <- distGeo(table.glm[i, c("long_midpoint", "lat_midpoint")],
                         water.df.i[j, c("long", "lat")])
      if (d < d.min) {
        d.min <- d
        j.min <- j
      }
    }
    # closest.water <- rbind(closest.water, c(ID_segment,
    #                                         water.df.i$ID[j.min],
    #                                         d.min,
    #                                         water.df.i$long[j.min], water.df.i$lat[j.min]))
    closest.water <- rbind(closest.water, c(ID_segment,
                                            water.df.i[j.min, "ID"],
                                            d.min,
                                            water.df.i[j.min, c("long", "lat")]))
    if (i %% 10 == 0) {
      print(i)
    }
  }
  closest.water <- as.data.frame(closest.water)
  colnames(closest.water) <- c("segment_ID", "ID_closest_water", "distance",
                               "long_water", "lat_water")
  return(closest.water)
}



# 8. Find closest amenity ----------------------------------------------------

ShortestDistanceAmenity <- function(table.glm, amenity.polygons.shp, amenity.points.shp) {
  ### Create a dataframe with the osm ID, name, type and coordinates of polygon amenity
  amenity.polygons.df.temp <- tidy(amenity.polygons.shp) %>%
    dplyr::select(ID = id,
                  long, lat)
  amenity.polygons.df <- c()
  for (k in 1:length(unique(amenity.polygons.df.temp$ID))) {
    amenity.polygons.k <- amenity.polygons.df.temp %>%
      filter(ID == unique(amenity.polygons.df.temp$ID)[k])
    amenity.polygons.df <- rbind(amenity.polygons.df,
                                   amenity.polygons.k[1,])
  }

  ### Create a dataframe with the osm ID, name, type and coordinates of points amenity
  amenity.points.df <- c()
  for (i in 1:nrow(amenity.points.shp@coords)) {
    row.i <- data.frame(ID = i + nrow(amenity.polygons.df),
                        long = amenity.points.shp@coords[i,1],
                        lat = amenity.points.shp@coords[i,2])
    amenity.points.df <- rbind(amenity.points.df, row.i)
  }

  ### Fuse the polygon building df and the point building df
  amenity.df <- rbind(amenity.polygons.df, amenity.points.df)

  ### Calculate the pair wise distances between all the points, and keep the minimum
  closest.amenity <- c()
  for (i in 1:nrow(table.glm)) {
    amenity.df.i <- amenity.df %>%             # This is just to reduce the computational burden:
      filter(abs(long - table.glm$long_midpoint[i]) < 0.20,  # for each segment I only consider amenity that
             abs(lat - table.glm$lat_midpoint[i]) < 0.20)    # are less than 0.20 degrees (latitude and longitude) away
    # from the midpoint. 0.28 degrees is a conservative value.

    dist.seg.i <- c()
    ID_segment <- table.glm$segment_ID[i]
    d.min <- 1000000000                                      # Initialize the minimum distance to a building
    # with a  purposefully inflated value
    for (j in 1:nrow(amenity.df.i)) {
      d <- distGeo(table.glm[i, c("long_midpoint", "lat_midpoint")],
                         amenity.df.i[j, c("long", "lat")])
      if (d < d.min) {
        d.min <- d
        j.min <- j
      }
    }
    closest.amenity <- rbind(closest.amenity,
                                    c(ID_segment,
                                      amenity.df.i$ID[j.min],
                                      d.min,
                                      amenity.df.i$long[j.min],
                                      amenity.df.i$lat[j.min]))
    if (i %% 10 == 0) {
      print(i)
    }
  }
  colnames(closest.amenity) <- c("segment_ID", "ID_closest_amenity",
                                        "distance", "long_amenity",
                                        "lat_amenity") #, "type", "name")
  closest.amenity <- as.data.frame(closest.amenity)
  return(closest.amenity)
}



# 10. Extract land cover percentages -----------------------------------------------

ExtractLandCover <- function(roads.segmented.df, land.cover, buffer.size) {
  # Convert the segmented road dataframe to shp
  coordinates(roads.segmented.df) =~ long+lat # transform into spatial points

  x <- lapply(split(roads.segmented.df, roads.segmented.df$ID_road_seg), # classifies points according to the segment to which they belong.
              function(x) Lines(list(Line(coordinates(x))),
                                x$ID_road_seg[1L]))
  lines <- SpatialLines(x) # transform points into SpatialLInes
  proj4string(lines) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
  data.for.shp <- data.frame(roads.segmented.df)  %>% # Dataframe containing the info that will be incorporated into the shp
    dplyr::select(ID = ID_road_seg) %>%
    distinct()
  rownames(data.for.shp) <- data.for.shp$ID
  roads.segmented.shp <- SpatialLinesDataFrame(lines,
                                               data.for.shp)
  roads.segmented.shp.planar <- spTransform(roads.segmented.shp, # Project to planar coodinates
                                            CRS("+proj=utm +zone=36 +south +a=6378249.145 +rf=293.46500000303 +towgs84=-160,-6,-302,0,0,0,0 +units=m +no_defs"))

  # Generate buffers around/along each road segment
  segment.buffers.planar <- gBuffer(roads.segmented.shp.planar, # generate a buffer around each road segment
                                    byid = TRUE, id = NULL,
                                    width = buffer.size,
                                    quadsegs = 5,
                                    capStyle = "ROUND",
                                    joinStyle = "ROUND",
                                    mitreLimit = 1.0)

  # Extract the land cover in this buffer
  land.cover.seg <- raster::extract(x = land.cover,             # Returns a list of n list (n = nbr of segments). Each of these sublists
                                    y = segment.buffers.planar, # contains the landcover value of each cell present in the buffer of a
                                    weights = TRUE,             # given segment, along with coefficients that corresponds to the proportion
                                    normalizeWeights=  TRUE,    # of each cell that was present within the buffer polygon
                                    sp = FALSE)
  # List of levels in the land cover data:
  land.cover.types <- c(0, 1, 2, 3, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 20,
                        23, 24, 25, 27, 28, 29, 32, 35, 36, 39, 40, 41)
  # For each segment, and for each land cover type, count the number of cells
  # with the land cover type.
  final <- c()
  for (i in 1:length(land.cover.seg)) {          # For each segment
    table <- c()
    for (j in 1:(length(land.cover.types))) {    # For each land cover type
      s <- 0
      for (k in 1:(length(land.cover.seg[[i]])/2)) {   # For each cell of segment i
        if (land.cover.seg[[i]][k, 1] == land.cover.types[j]) {
          s <- s + land.cover.seg[[i]][k, 2]
        }
      }
      table <- rbind(table, c(segment.buffers.planar@polygons[[i]]@ID, # Segment ID
                              land.cover.types[j], s))
    }
    final <- rbind(final, table) # table with: the segment number (i),
                                 # each landcover value, and the
                                 # proportion of each land cover type.
    if (i %% 10 == 0) {
      print(i)
    }
  }
  colnames(final) <- c("segment_ID", "land_cover", "weight")
  return(as.data.frame(final))
}



# Argument: the df returned by the ExtractLandCover function
# Output: df with a row for each segment, and a column for the proportion of each land cover along a road segment
# (colummns corresponding to land covers that are completely absent from the buffers are deleted).
PercentageLandCover <- function(df) {
  segment_IDs <- levels(as.factor(df$segment_ID))
  table.land.covers <- c()
  for (i in 1:length(segment_IDs)) {
    temp.i <- df %>%
      filter(df$segment_ID == segment_IDs[i]) %>%
      dplyr::select(weight) %>%
      mutate(weight = 100 * as.numeric(weight))
    temp.i.transposed <- c(segment_IDs[i], as.numeric(t(temp.i)))
    table.land.covers <- rbind(table.land.covers, temp.i.transposed)
  }
  table.land.covers <- as.data.frame(table.land.covers)
  colnames(table.land.covers) <- c("segment_ID",
                                   "unclassified", "bare_ground","water","clouds", "sparse_grassland",
                                   "open_grassland", "dense_grassland", "closed_grassland",
                                   "sparse_shrubbed_grassland", "open_shrubbed_grassland",
                                   "dense_shrubbed_grassland", "closed_shrubbed_grassland",
                                   "open_treed_grassland", "dense_treed_grassland","closed_treed_grassland",
                                   "dense_shrubland","open_grassed_shrubland",
                                   "dense_grassed_shrubland","closed_grassed_shrubland",
                                   "open_treed_shrubland","dense_treed_shrubland", "closed_treed_shrubland",
                                   "dense_forest", "open_grassed_woodland", "dense_grassed_woodland",
                                   "closed_grassed_woodland", "dense_shrubbed_forest",
                                   "closed_shrubbed_forest")

  return(table.land.covers)
}
