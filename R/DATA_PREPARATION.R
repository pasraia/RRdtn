#' @importFrom raster projection crs
#' @importFrom sp proj4string SpatialPoints
#' @importFrom raster extract

DATA_PREPARATION<-function (species_input_data,
                            obs_col = "OBS",
                            input_mask,
                            time_col = NULL){
  if (!projection(input_mask) == proj4string(species_input_data) |
      is.na(projection(input_mask)) | is.na(proj4string(species_input_data)))
    stop("some missing or not-matching projections in input data or background area")
  species_input_data <- as.data.frame(species_input_data)
  coords_cols<- grep("coords", colnames(species_input_data))
  ones <- subset(species_input_data, species_input_data[, obs_col] ==
                   1)
  back <- subset(species_input_data, species_input_data[, obs_col] ==
                   0)
  input_ones <- data.frame(ones[, -coords_cols], geoID = extract(input_mask,
                                                                 ones[, coords_cols], cellnumbers = TRUE)[, "cells"])
  input_back <- data.frame(back[, -coords_cols], geoID = extract(input_mask,
                                                                 back[, coords_cols], cellnumbers = TRUE)[, "cells"])
  if (!all(input_ones$geoID %in% input_back$geoID))
    stop("some points fall outside of the background area")
  ones_coords <- SpatialPoints(ones[, coords_cols])
  raster::crs(ones_coords) <- projection(input_mask)
  OUTPUT <- list(call = "input_data", input_ones = input_ones,
                 input_back = input_back, obs_col = obs_col, geoID_col = "geoID",
                 time_col = time_col, ones_coords = ones_coords,
                 study_area = input_mask)
  return(OUTPUT)
}
