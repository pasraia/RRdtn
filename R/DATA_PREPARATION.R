#' @importFrom sf st_crs st_coordinates st_drop_geometry
#' @importFrom grDevices dev.new
#' @importFrom graphics points
#' @importFrom terra extract crs cellFromXY

DATA_PREPARATION<-function(species_input_data,
                           obs_col = "OBS",
                           input_mask,
                           time_col = NULL){
  if (!crs(input_mask, proj = TRUE) == st_crs(species_input_data)$proj4string |
      is.na(crs(input_mask, proj = TRUE)) | is.na(st_crs(species_input_data)$proj4string))
    stop("some missing or not-matching projections in input data or background area")
  occ.desaggregation.RASTER <- function(df, colxy, rast, plot = T) {
    df_ini <- df
    cacca <- cellFromXY(rast, df[, colxy])
    if (any(is.na(cacca))) {
      stop("no NA admitted in species occurrences!")
    }
    cacca1 <- split(cacca, cacca)
    l <- sapply(cacca1, length)
    if (max(l) > 1) {
      l1 <- l[l > 1]
      df_ok <- df[which(cacca %in% as.numeric(names(l[l ==
                                                        1]))), ]
      found <- lapply(1:length(l1), function(j) {
        w <- which(cacca %in% as.numeric(names(l1[j])))
        s <- sample(w, 1)
        df[s, ]
      })
      df_final <- rbind(df_ok, do.call(rbind, found))
      if (plot == T) {
        dev.new()
        plot(df_ini[, colxy], main = "distribution of occurences",
             sub = paste("# initial (black):", nrow(df_ini),
                         " | # kept (red): ", nrow(df_final)), pch = 19,
             col = "black", cex = 0.5)
        points(df_final[, colxy], pch = 19, col = "red",
               cex = 0.2)
      }
      return(df_final)
    }
    if (max(l) == 1)
      return(df)
  }
  species_input_data <- cbind(st_coordinates(species_input_data),
                              as.data.frame(st_drop_geometry(species_input_data)))
  species_input_data0 <- subset(species_input_data, species_input_data[,
                                                                       obs_col] == 0)
  species_input_data <- subset(species_input_data, species_input_data[,
                                                                      obs_col] == 1)
  if (is.null(time_col)) {
    species_input_data2 <- occ.desaggregation.RASTER(species_input_data,
                                                     1:2, rast = input_mask, plot = F)
    if (nrow(species_input_data2) < nrow(species_input_data))
      warning(paste(nrow(species_input_data) - nrow(species_input_data2),
                    "duplicated points into raster cells were found and removed!"))
    species_input_data <- species_input_data2
  }else {
    ss <- split(species_input_data, species_input_data[,
                                                       time_col])
    species_input_data2 <- lapply(ss, function(x) {
      occ.desaggregation.RASTER(x, 1:2, rast = input_mask,
                                plot = F)
    })
    species_input_data2 <- do.call(rbind, species_input_data2)
    if (nrow(species_input_data2) < nrow(species_input_data))
      warning(paste(nrow(species_input_data) - nrow(species_input_data2),
                    "duplicated points into raster cells were found and removed!"))
    species_input_data <- species_input_data2
  }
  species_input_data <- rbind(species_input_data, species_input_data0)
  colnames(species_input_data)[1:2] <- c("coordsX", "coordsY")
  coords_cols <- grep("coords", colnames(species_input_data))
  ones <- subset(species_input_data, species_input_data[, obs_col] ==
                   1)
  back <- subset(species_input_data, species_input_data[, obs_col] ==
                   0)
  input_ones <- cbind.data.frame(ones[, -coords_cols], geoID = extract(input_mask,
                                                                       ones[, coords_cols], cells = TRUE)[, "cell"])
  input_back <- cbind.data.frame(back[, -coords_cols], geoID = extract(input_mask,
                                                                       back[, coords_cols], cells = TRUE)[, "cell"])
  if (!all(input_ones$geoID %in% input_back$geoID))
    stop("some points fall outside of the background area")
  ones_coords <- ones[, coords_cols]
  rownames(ones_coords) <- NULL
  OUTPUT <- list(call = "input_data", input_ones = input_ones,
                 input_back = input_back, obs_col = obs_col, geoID_col = "geoID",
                 time_col = time_col, ones_coords = ones_coords, study_area = input_mask)
  return(OUTPUT)
}
