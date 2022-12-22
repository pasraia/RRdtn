#' @importFrom methods slot
#' @importFrom sp gridded gridparameters coordinates
#' @importFrom adehabitatMA .adehabitatMAEnv join
#' @importFrom ade4 acm.disjonctif

mahasuhab.custom<-function (x, pts)
{
  # if (class(x) == "SpatialPixelsDataFrame" & extends(class(pts),
  #                                                    "SpatialPoints")) {
  if (inherits(x,"SpatialPixelsDataFrame") & extends(class(pts),
                                                     "SpatialPoints")) {
    if (!inherits(x, "SpatialPixelsDataFrame"))
      stop("should be an object of class SpatialPixelsDataFrame")
    sp::gridded(x) <- TRUE
    gr <- gridparameters(x)
    if (nrow(gr) > 2)
      stop("x should be defined in two dimensions")
    if ((gr[1, 2] - gr[2, 2]) > get(".adeoptions",
                                    envir = .adehabitatMAEnv)$epsilon)
      stop("the cellsize should be the same in x and y directions")
    if (!inherits(pts, "SpatialPoints"))
      stop("should inherit from class \"SpatialPoints\"")
    hihi <- join(pts, x)
    if (is.null(dim(hihi)))
      hihi <- cbind(hihi)
    hihi <- hihi[!is.na(hihi[, 1]), , drop = FALSE]
    used <- list()
    for (i in 1:ncol(hihi)) {
      if (is.factor(hihi[, i]))
        used[[i]] <- acm.disjonctif(data.frame(hihi[,
                                                    i]))[, -1] else used[[i]] <- hihi[, i]
    }
    used[[i + 1]] <- rep(1, nrow(hihi))
    hihi <- as.data.frame(used)
    hihi <- hihi[!is.na(hihi[, 1]), ]
    mu <- apply(hihi, 2, function(x) mean(x, na.rm = TRUE))
    varcov <- t(as.matrix(hihi)) %*% as.matrix(hihi)/nrow(hihi)
    kasc <- slot(x, "data")
    ava <- list()
    for (i in 1:ncol(kasc)) {
      if (is.factor(kasc[, i]))
        ava[[i]] <- acm.disjonctif(data.frame(kasc[,
                                                   i]))[, -1] else ava[[i]] <- kasc[, i]
    }
    ava[[i + 1]] <- rep(1, nrow(kasc))
    df <- as.data.frame(ava)
    map <- mahalanobis(as.matrix(df), mu, varcov)
    map <- 1 - pchisq(map, ncol(hihi) - 1)
    map <- data.frame(MD = map)
    sp::coordinates(map) <- sp::coordinates(x)
    sp::gridded(map) <- TRUE
  } else {
    hihi <- pts
    if (is.null(dim(hihi)))
      hihi <- cbind(hihi)
    hihi <- hihi[!is.na(hihi[, 1]), , drop = FALSE]
    used <- list()
    for (i in 1:ncol(hihi)) {
      if (is.factor(hihi[, i]))
        used[[i]] <- acm.disjonctif(data.frame(hihi[,
                                                    i]))[, -1] else used[[i]] <- hihi[, i]
    }
    used[[i + 1]] <- rep(1, nrow(hihi))
    hihi <- as.data.frame(used)
    hihi <- hihi[!is.na(hihi[, 1]), ]
    mu <- apply(hihi, 2, function(xx) mean(xx, na.rm = TRUE))
    varcov <- t(as.matrix(hihi)) %*% as.matrix(hihi)/nrow(hihi)
    # if (class(x) == "SpatialPixelsDataFrame")
    if (inherits(x,"SpatialPixelsDataFrame"))
      kasc <- slot(x, "data") else kasc <- x
    ava <- list()
    for (i in 1:ncol(kasc)) {
      if (is.factor(kasc[, i]))
        ava[[i]] <- acm.disjonctif(data.frame(kasc[,
                                                   i]))[, -1] else ava[[i]] <- kasc[, i]
    }
    ava[[i + 1]] <- rep(1, nrow(kasc))
    df <- as.data.frame(ava)
    map <- mahalanobis(as.matrix(df), mu, varcov)
    map <- 1 - pchisq(map, ncol(hihi))
    map <- data.frame(MD = map)
  }
  return(map)
}
