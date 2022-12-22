#'@title Project the ENFA and ENphylo models into new geographical space and
#'  time interval
#'@description The function projects species marginality and specialization
#'  factors in other areas and other time scales. The function is able to
#'  convert marginality and specialization factors in habitat suitability values
#'  by using the Mahalanobis distances method.
#'@usage ENphylo_prediction(object, newdata, calc_method=c("raster", "terra"),
#'  convert.to.suitability=FALSE)
#'@param object a \code{list} returned by the \code{\link{ENphylo_modeling}}
#'  function where all the species are modelled with ENFA or ENphylo models.
#'@param newdata a set of explanatory variables onto which ENFA or ENphylo
#'  models will be projected. It could be a \code{data.frame} or a
#'  \code{rasterStack} object.
#'@param calc_method character. The R package and related functions which
#'  \code{ENphylo_prediction} uses to obtain model predictions.
#'@param convert.to.suitability logical. If \code{TRUE},
#'  \code{ENphylo_prediction} projects ENFA or ENphylo model predictions in
#'  other areas and other time scales.
#'@author Alessandro Mondanaro, Mirko Di Febbraro, Silvia Castiglione,
#'  Carmela Serio, Marina Melchionna, Pasquale Raia
#'@details If \code{convert.to.suitability} is set as \code{TRUE},
#'  \code{ENphylo_prediction} uses the function
#'  \code{\link[adehabitatHS]{mahasuhab}} from the \pkg{adehabitatHS} R
#'  package (\cite{Calenge, 2006}) to compute the habitat suitability map of the
#'  species over a given area. The conversion of Mahalanobis distances into
#'  probabilities follows the chi-squared distribution. Specifically, we set the
#'  degree of freedom equal to \emph{n} rather than \emph{n-1} following
#'  \cite{Etherington (2019)}. To convert habitat suitability values into binary
#'  presence/absence values, \code{ENphylo_prediction} relies on
#'  \code{\link{optimal.thresholds}} function in the package
#'  \pkg{PresenceAbsence} (\cite{Freeman & Moisen, 2008}).
#'@return A list of length equal to the number of modelled species. For each
#'  element of the list, the output objects differ depending on whether ENFA or
#'  ENphylo is implemented and if a raster or data.frame object is provided by
#'  the user as new environments where to project the model. Furthermore, if the
#'  \code{convert.to.suitability} is set to \code{TRUE}, each element of the
#'  list includes the habitat suitability values in addition to model
#'  predictions for marginality and specialization. In any case, each element of
#'  the output list contains two elements:
#'@return ENFA/raster
#'@return \enumerate{ \item\strong{$call}: a character specifying the algorithm
#'  used for the model predictions. \item\strong{$enfa_prediction}: a list of
#'  RasterLayer objects including marginality and specialization model
#'  predictions plus the habitat suitability and binary presence/absence maps
#'  obtained by using three different predefined thresholds: MaxSensSpec (i.e.
#'  maximize TSS), SensSpec (i.e. equalize sensitivity and specificity) and 10th
#'  percentile of predicted probability. Habitat suitability estimates and
#'  binary values are calculated only if \code{convert.to.suitability = TRUE}.}
#'@return ENFA/data.frame
#'@return \enumerate{ \item\strong{$call}: a character specifying the algorithm
#'  used for the model predictions.\item\strong{$enfa_prediction}: a data.frame
#'  containing marginality and specialization model predictions and the habitat
#'  suitability values (only if \code{convert.to.suitability = TRUE}).}
#'@return ENphylo/raster
#'@return \enumerate{ \item\strong{$call}: a character specifying the algorithm
#'  used for the model predictions. \item\strong{$imputed_prediction}: a list of
#'  raster objects including marginality and specialization model predictions.
#'  The layer number of each raster is equal to the integer number chosen for
#'  testing alternative phylogenies (i.e. \code{nsim} argument). If
#'  \code{convert.to.suitability = TRUE}, in addition to marginality and
#'  specialization, the function provides raster objects containing habitat
#'  suitability and binary presence/absence maps (see above) for each tested
#'  phylogeny.}
#'@return ENphylo/data.frame
#'@return \enumerate{ \item\strong{$call}: a character specifying the algorithm
#'  used for the model predictions. \item\strong{$imputed_prediction}: a list of
#'  data.frame objects of length equal to the integer number chosen for testing
#'  alternative phylogenies (i.e. \code{nsim} argument) containing marginality
#'  and specialization model predictions and the habitat suitability values
#'  (only if \code{convert.to.suitability argument = TRUE}).}
#'@importFrom raster stack nlayers
#'@importFrom terra app rast
#'@importFrom biomod2 bm_BinaryTransformation
#'@importFrom stats quantile
#'@export
#'@references Calenge, C. (2006) The package adehabitat for the R software: a
#'  tool for the analysis of space and habitat use by animals. \emph{Ecological
#'  Modelling}, 197, 516-519.
#'@references Etherington, T. R. (2019). Mahalanobis distances and ecological
#'  niche modelling: correcting a chi-squared probability error. \emph{PeerJ},
#'  7, e6678.
#'@references Freeman, E. A. & Moisen, G. (2008). PresenceAbsence: An R Package
#'  for Presence-Absence Model Analysis. \emph{Journal of Statistical Software},
#'  23(11):1-31.
#'@examples
#' \dontrun{
#'setwd("YOUR_DIRECTORY")
#'getwd()->main.dir
#'url<-"https://www.dropbox.com/s/pdcf1psidyzjkjk/ENphylo%20code%26data.zip?dl=1"
#'download.file(url,file.path(main.dir,"ENphylo code&data.zip"),mode="wb")
#'unzip("ENphylo code&data.zip")
#'load("example_data.RData")
#'MASK_FULL<-raster::raster("variable_bio1.tif")
#'external_data<-raster::stack(list.files("external_data",full.names=TRUE))
#'
#' ## NOTE: Given the size of the data, running the function is time-comsuming
#' ### CASE 1
#' ENmod<-ENphylo_modeling(input_data=DATA_FULL,
#'                        tree=tree,
#'                        input_mask=MASK_FULL,
#'                        obs_col="OBS",
#'                        time_col="TIME_factor",
#'                        eval_metric_for_imputation = "AUC",
#'                        output_options = "best")
#' gc()
#'
#' # Predicting for the first species in ENmod output (available 1-31)
#' ENmod_pred<-ENphylo_prediction(object=ENmod[1],
#'                               newdata=external_data,
#'                               convert.to.suitability=TRUE,
#'                               calc_method="terra")
#'
#' ### CASE 2
#' ## Simulating a rare species by subsampling its occurrences to force run ENphylo algorithm
#' k=10 # index of the species to subsample as in DATA_FULL (k available 1-31)
#' npoints=10 # Number of randomly selected presence data points for the k species in DATA_FULL
#'
#' # Subsampling procedure
#' Species_data_start<-DATA_FULL[[k]]
#' Species_data1_start<-subset(Species_data_start, OBS==1)
#' Species_data0_start<-subset(Species_data_start, OBS==0)
#' sample(nrow(Species_data1_start),npoints)->zz
#' Species_data1<-Species_data1_start[zz,]
#' Species_data0<-subset(Species_data0_start,TIME_factor%in%Species_data1$TIME_factor)
#' Species_data<-rbind(Species_data1,Species_data0)
#' Species_data->DATA_FULL[[k]]
#'
#' # Running ENphylo
#' ENmod_subsam<-ENphylo_modeling(input_data=DATA_FULL,
#'                               tree=tree,
#'                               input_mask=MASK_FULL,
#'                               obs_col="OBS",
#'                               time_col="TIME_factor",
#'                               min_occ_enfa=npoints+1,
#'                               eval_metric_for_imputation="AUC",
#'                               eval_threshold=0.7,
#'                               output_options = "best")
#' gc()
#'
#' # Predicting for the k species in ENmod_subsam output
#' ENmod_subsam_pred<-ENphylo_prediction(object=ENmod_subsam[k],
#'                                      newdata=external_data,
#'                                      convert.to.suitability=TRUE,
#'                                      calc_method="terra")
#'
#'
#' ### CASE 3
#' ## In case you do not trust an ENFA-derived SDM (regardless of the evaluation score)
#' ## you might want to run ENphylo anyway, retrieving the ENFA models and specifing
#' ## species k you what to look at ENphylo prediction for
#'
#' # Extract a species modelled with ENFA under ENphylo_modeling and 'force' ENphylo calibration on it
#' sample(which(sapply(ENmod,function(x) x$call)=="calibrated_enfa"),1)->k
#'
#'
#' ENmod3<-ENphylo_modeling(external_enfa_models=ENmod,
#'                              spec_for_imputation=names(ENmod)[k],
#'                              tree=tree,
#'                              output_options="full",
#'                              eval_metric_for_imputation="AUC")
#'
#' gc()
#'
#' # Transform the raster object to a data.frame object to predict on data.frame
#' external_dataframe<-raster::rasterToPoints(external_data)[,-c(1,2)]
#'
#' # Predicting for the k species in ENmod_enfa1 output
#' ENmod3_pred<-ENphylo_prediction(object=ENmod3[k],
#'                                     newdata=external_dataframe,
#'                                     convert.to.suitability=TRUE,
#'                                     calc_method="raster")
#'
#'
#' ### CASE 4
#' # Run ENphylo under alternative model evaluation metric and strategy
#' ENmod4<-ENphylo_modeling(external_enfa_models=ENmod,
#'                              spec_for_imputation=names(ENmod)[k],
#'                              tree=tree,
#'                              output_options="best",
#'                              eval_metric_for_imputation="TSS",
#'                              eval_threshold = 0.4)
#' gc()
#'
#' # Predicting for the k species in ENmod_enfa2 output
#' ENmod4_pred<-ENphylo_prediction(object=ENmod4[k],
#'                                     newdata=external_dataframe,
#'                                     convert.to.suitability=TRUE,
#'                                     calc_method="terra")
#'}

ENphylo_prediction<-function (object,
                              newdata,
                              calc_method = c("raster", "terra"),
                              convert.to.suitability = FALSE){

  if(!(extends(class(newdata),"Raster")|is.matrix(newdata)|is.data.frame(newdata)))
    stop("Please, provide newdata as a Raster* object or a data.frame")

  x <- newdata
  if (is.null(names(object))) {
    if (object[[1]]$call == "enfa") {
      U <- object[[1]]$co
      f1 <- function(y) y %*% U
      if (extends(class(x), "Raster")) {
        if(!(all(names(x)%in%rownames(U))&all(rownames(U)%in%names(x))))
          stop("Variable names in newdata must match with those used to calibrate models") else{
          x<-x[[match(rownames(U),names(x))]]
        }
        if (calc_method == "raster")
          ras <- raster::calc(x, fun = f1)
        if (calc_method == "terra")
          ras <- stack(app(rast(x), fun = f1))
        names(ras) <- names(object[[1]]$sf)
      }else {
        if(!(all(colnames(x)%in%rownames(U))&all(rownames(U)%in%colnames(x))))
          stop("Variable names in newdata must match with those used to calibrate models") else{
            x<-x[,match(rownames(U),colnames(x))]
          }

        ras <- t(apply(x, 1, f1))
        colnames(ras) <- names(object[[1]]$sf)
        ras <- as.data.frame(ras)
      }
    }
    if (object[[1]]$call == "calibrated_enfa") {
      U <- object[[1]]$full_model$co
      f1 <- function(y) y %*% U
      if (extends(class(x), "Raster")) {
        if(!(all(names(x)%in%rownames(U))&all(rownames(U)%in%names(x))))
          stop("Variable names in newdata must match with those used to calibrate models") else{
            x<-x[[match(rownames(U),names(x))]]
          }
        if (calc_method == "raster")
          ras <- raster::calc(x, fun = f1)
        if (calc_method == "terra")
          ras <- stack(app(rast(x), fun = f1))
        names(ras) <- names(object[[1]]$full_model$sf)
      }else {
        if(!(all(colnames(x)%in%rownames(U))&all(rownames(U)%in%colnames(x))))
          stop("Variable names in newdata must match with those used to calibrate models") else{
            x<-x[,match(rownames(U),colnames(x))]
          }
        ras <- t(apply(x, 1, f1))
        colnames(ras) <- names(object[[1]]$full_model$sf)
        ras <- as.data.frame(ras)
      }
    }
    ras <- list(call = "enfa_prediction", enfa_prediction = ras)
  }else {
    cat(paste("\n", "PREDICTING ENFA/IMPUTATION",
              "\n"))
    ras <- pblapply(object, function(ob) {
      if (ob$call == "calibrated_enfa") {
        U <- ob$calibrated_model$full_model$co
        f1 <- function(y) y %*% U
        if (extends(class(x), "Raster")) {
          if(!(all(names(x)%in%rownames(U))&all(rownames(U)%in%names(x))))
            stop("Variable names in newdata must match with those used to calibrate models") else{
              x<-x[[match(rownames(U),names(x))]]
            }
          if (calc_method == "raster")
            ras <- raster::calc(x, fun = f1)
          if (calc_method == "terra")
            ras <- stack(app(rast(x), fun = f1))
          names(ras) <- names(ob$calibrated_model$full_model$sf)
          ras <- ras[[1:ob$calibrated_model$full_model$significant_axes]]
          ras
        }else {
          if(!(all(colnames(x)%in%rownames(U))&all(rownames(U)%in%colnames(x))))
            stop("Variable names in newdata must match with those used to calibrate models") else{
              x<-x[,match(rownames(U),colnames(x))]
            }
          ras <- t(apply(x, 1, f1))
          colnames(ras) <- names(ob$calibrated_model$full_model$sf)
          ras <- as.data.frame(ras[, 1:ob$calibrated_model$full_model$significant_axes])
          ras
        }
        ras <- list(call = "enfa_prediction", enfa_prediction = ras)
      }
      if (ob$call == "calibrated_imputed") {
        ras <- lapply(ob$calibrated_model$co, function(co) {
          U <- co
          f1 <- function(y) y %*% U
          if (extends(class(x), "Raster")) {
            if(!(all(names(x)%in%rownames(U))&all(rownames(U)%in%names(x))))
              stop("Variable names in newdata must match with those used to calibrate models") else{
                x<-x[[match(rownames(U),names(x))]]
              }
            if (calc_method == "raster")
              ras <- raster::calc(x, fun = f1)
            if (calc_method == "terra")
              ras <- stack(app(rast(x), fun = f1))
            names(ras) <- colnames(U)
            ras
          }else {
            if(!(all(colnames(x)%in%rownames(U))&all(rownames(U)%in%colnames(x))))
              stop("Variable names in newdata must match with those used to calibrate models") else{
                x<-x[,match(rownames(U),colnames(x))]
              }
            ras <- t(apply(x, 1, f1))
            colnames(ras) <- colnames(U)
            ras <- as.data.frame(ras)
          }
        })
        if (extends(class(x), "Raster")) {
          ras <- lapply(1:nlayers(ras[[1]]), function(i) stack(lapply(ras,
                                                                      "[[", i)))
          names(ras)[1] <- "Marg"
          names(ras)[-1] <- paste("Spec", 1:(length(ras) -
                                               1), sep = "")
          ras <- lapply(ras, function(o) {
            if (any(grepl("evaluation", names(ob$calibrated_model)))) {
              names(o) <- rownames(ob$calibrated_model$evaluation)
            }else {
              names(o) <- paste("swap", 1:nlayers(o),
                                sep = "_")
            }
            o
          })
        }else {
          if (any(grepl("evaluation", names(ob$calibrated_model)))) {
            ras <- mapply(function(x, y) {
              x$swap_rep <- y
              x
            }, x = ras, y = rownames(ob$calibrated_model$evaluation),
            SIMPLIFY = FALSE)
          }else {
            ras <- mapply(function(x, y) {
              x$swap_rep <- y
              x
            }, x = ras, y = paste("swap", 1:length(ras),
                                  sep = "_"), SIMPLIFY = FALSE)
          }
        }
        if (ob$calibrated_model$output_options[1] ==
            "best" & extends(class(x), "Raster")) {
          ras2 <- ras
          ras <- stack(ras)
          names(ras) <- paste(names(ras2), strsplit(names(ras[[1]]),
                                                    "[.]")[[1]][1], sep = "_")
          ras
        }

        if (ob$calibrated_model$output_options[1] ==
            "weighted.mean" & ob$calibrated_model$output_options[2]=="OMR") warning("The weighted mean of model predictions using OMR provides misleading results")

        if (ob$calibrated_model$output_options[1] ==
            "weighted.mean") {
          if (extends(class(x), "Raster")) {
            ras <- list(Reduce("+", Map("*",
                                        lapply(1:nlayers(ras[[1]]), function(i) stack(lapply(ras,
                                                                                             "[[", i))), ob$calibrated_model$evaluation[,
                                                                                                                                        ob$calibrated_model$output_options[2]]))/nlayers(ras[[1]]))
            rr <- lapply(1:nlayers(ras[[1]]), function(xx) {
              dd <- ras[[1]][[xx]]
              dd
            })
            rr <- stack(rr)
            names(rr) <- paste(names(ras[[1]]), "swap_weighted",
                               sep = "_")
            ras <- rr
            ras
          }else {
            ras <- list(Reduce("+", Map("*",
                                        lapply(ras, function(xx) xx[, !grepl("swap_rep",
                                                                             colnames(xx))]), ob$calibrated_model$evaluation[,
                                                                                                                             ob$calibrated_model$output_options[2]]))/length(ras))
            ras[[1]]$swap_rep <- "swap_weighted"
            ras
          }
        }
        ras <- list(call = "imputed_prediction",
                    imputed_prediction = ras)
      }
      return(ras)
    })
  }
  if (convert.to.suitability) {
    cat(paste("\n", "CONVERTING PREDICTED VALUES TO SUITABILITY",
              "\n"))
    if (extends(class(x), "Raster")) {
      reference <- rasterToPoints(x)[, -c(1:2)]
    }else reference <- x
    ras_suitability <- pblapply(names(object), function(sp) {
      cat(paste("\n", sp, "\n"))
      mydata <- object[[sp]]$formatted_data
      obs_col <- mydata$obs_col
      time_col <- mydata$time_col
      geoID_col <- mydata$geoID_col
      reference <- as.data.frame(reference)
      reference$type <- factor("reference", levels = c("obs_background",
                                                       "reference"))
      reference[, obs_col] <- 0
      if (nrow(mydata$input_back) > 10000) {
        independent_data_for_pred <- rbind(mydata$input_ones,
                                           mydata$input_back[sample(nrow(mydata$input_back),
                                                                    10000), ])
      }else independent_data_for_pred <- rbind(mydata$input_ones,
                                               mydata$input_back)
      obs_background <- independent_data_for_pred
      obs_background <- obs_background[, !grepl(geoID_col,
                                                colnames(obs_background))]
      obs_background$type <- factor("obs_background",
                                    levels = c("obs_background", "reference"))
      obs_background <- obs_background[, match(colnames(reference),
                                               colnames(obs_background))]
      ras_for_proj <- rbind(reference, obs_background)
      ras_obs <- ras_for_proj[, obs_col, drop = FALSE]
      ras_type <- ras_for_proj[, "type", drop = FALSE]
      ras_for_proj <- ras_for_proj[, !grepl(paste(c("type",
                                                    obs_col), collapse = "|"), colnames(ras_for_proj))]
      if (object[[sp]]$call == "calibrated_enfa") {
        U <- object[[sp]]$calibrated_model$full_model$co
        f1 <- function(y) y %*% U
        ma1 <- t(apply(as.matrix(ras_for_proj), 1, f1))
        colnames(ma1) <- colnames(U)
        ma1 <- as.data.frame(scale(ma1[, 1:object[[sp]]$calibrated_model$full_model$significant_axes]))
        ma1 <- list(ma1)
      }
      if (object[[sp]]$call == "calibrated_imputed") {
        ma1 <- lapply(1:length(object[[sp]]$calibrated_model$co),
                      function(kk) {
                        U <- object[[sp]]$calibrated_model$co[[kk]]
                        f1 <- function(y) y %*% U
                        ma1 <- t(apply(as.matrix(ras_for_proj), 1,
                                       f1))
                        colnames(ma1) <- colnames(U)
                        ma1 <- as.data.frame(scale(ma1))
                      })
        if (object[[sp]]$calibrated_model$output_options[1] ==
            "weighted.mean") {
          ma1 <- list(Reduce("+", Map("*",
                                      ma1, object[[sp]]$calibrated_model$evaluation[,
                                                                                    object[[sp]]$calibrated_model$output_options[2]]))/length(ma1))
        }else ma1 <- ma1
      }
      suitability_final <- lapply(1:length(ma1), function(kk) {
        repeat {
          hsm1 <- suppressWarnings(try(mahasuhab.custom(x = ma1[[kk]],
                                                        pts = ma1[[kk]][which(ras_obs == 1), ]),
                                       silent = TRUE))
          # if (class(hsm1) != "try-error")
          #   break else {
          if(inherits(hsm1,"try-error")){
              ma1[[kk]] <- ma1[[kk]][, -ncol(ma1[[kk]]),
                                     drop = FALSE]
              warning(paste("Dimensionality was reduced to",
                            ncol(ma1[[kk]]), "axes"))
            }else break
        }
        ma1[[kk]]$Suitability <- hsm1$MD
        reference <- ma1[[kk]][which(ras_type == "reference"),
        ]
        pred_for_th <- ma1[[kk]][which(ras_type == "obs_background"),
        ]
        data_for_th <- data.frame(ID = 1:nrow(pred_for_th),
                                  observed = obs_background$OBS, pred = pred_for_th$Suitability)
        suppressWarnings(TH <- optimal.thresholds(data_for_th)[2:3,
                                                               2])
        TH <- c(TH, quantile(data_for_th[which(data_for_th$observed ==
                                                 1), ]$pred, 0.1))
        if (extends(class(x), "Raster")) {
          hsm1 <- x[[1]]
          hsm1[!is.na(hsm1)] <- reference$Suitability
          names(hsm1) <- "Suitability"
          if (object[[sp]]$call == "calibrated_imputed") {
            if (any(grepl("evaluation", names(object[[sp]]$calibrated_model)))) {
              if (object[[sp]]$calibrated_model$output_options[1] ==
                  "best") {
                names(hsm1) <- paste("Suitability",
                                     rownames(object[[sp]]$calibrated_model$evaluation)[kk],
                                     sep = "_")
              }
              if (object[[sp]]$calibrated_model$output_options[1] ==
                  "weighted.mean") {
                names(hsm1) <- paste(names(hsm1), "swap_weighted",
                                     sep = "_")
              }
              if (object[[sp]]$calibrated_model$output_options[1] ==
                  "full") {
                names(hsm1) <- paste("Suitability",
                                     "swap", kk, sep = "_")
              }
            }else names(hsm1) <- paste("Suitability",
                                       "swap", kk, sep = "_")
          }
        }else {
          hsm1 <- reference[, "Suitability", drop = FALSE]
        }
        return(list(TH, hsm1))
      })
      if (object[[sp]]$call == "calibrated_enfa") {
        if (extends(class(x), "Raster")) {
          thh <- lapply(suitability_final, function(xx) {
            xx <- stack(sapply(xx[[1]], function(yy) bm_BinaryTransformation(xx[[2]],
                                                                          yy)))
            names(xx) <- c("SensSpec", "MaxSensSpec",
                           "TenPerc")
            xx
          })
          names(thh[[1]]) <- paste("Binary", names(thh[[1]]),
                                   sep = "_")
          ras[[sp]]$enfa_prediction <- stack(ras[[sp]]$enfa_prediction,
                                             suitability_final[[1]][[2]], thh[[1]])
        }else {
          ras[[sp]]$enfa_prediction <- mapply(function(a,
                                                       b) {
            cbind(a, b[[2]])
          }, a = list(ras[[sp]]$enfa_prediction), b = suitability_final,
          SIMPLIFY = FALSE)
        }
        ras_final <- list(call = "enfa_prediction",
                          enfa_prediction = ras[[sp]]$enfa_prediction)
      }
      if (object[[sp]]$call == "calibrated_imputed") {
        if (extends(class(x), "Raster")) {
          thh <- lapply(suitability_final, function(xx) {
            xx <- stack(sapply(xx[[1]], function(yy) bm_BinaryTransformation(xx[[2]],
                                                                          yy)))
            names(xx) <- c("SensSpec", "MaxSensSpec",
                           "TenPerc")
            xx
          })
          thh <- lapply(1:3, function(i) stack(lapply(thh,
                                                      "[[", i)))
          if (any(grepl("evaluation", names(object[[sp]]$calibrated_model)))) {
            if (object[[sp]]$calibrated_model$output_options[1] ==
                "best") {
              thh <- lapply(thh, function(jj) {
                names(jj) <- paste("Binary", strsplit(names(jj[[1]]),
                                                      "[.]")[[1]][1], rownames(object[[sp]]$calibrated_model$evaluation),
                                   sep = "_")
                jj
              })
            }
            if (object[[sp]]$calibrated_model$output_options[1] ==
                "weighted.mean") {
              thh <- lapply(thh, function(jj) {
                names(jj) <- paste("Binary", strsplit(names(jj[[1]]),
                                                      "[.]")[[1]][1], "swap_weighted",
                                   sep = "_")
                jj
              })
            }
            if (object[[sp]]$calibrated_model$output_options[1] ==
                "full") {
              thh <- lapply(thh, function(jj) {
                names(jj) <- paste(strsplit(names(jj[[1]]),
                                            "[.]")[[1]][1], "swap",
                                   1:nlayers(jj), sep = "_")
                jj
              })
            }
          }else {
            thh <- lapply(thh, function(jj) {
              names(jj) <- paste(strsplit(names(jj[[1]]),
                                          "[.]")[[1]][1], "swap", 1:nlayers(jj),
                                 sep = "_")
              jj
            })
          }
          if (object[[sp]]$calibrated_model$output_options[1] ==
              "full") {
            ras[[sp]]$imputed_prediction <- c(ras[[sp]]$imputed_prediction,
                                              Suitability = stack(lapply(suitability_final,
                                                                         "[[", 2)), Binary_SensSpec = thh[[1]],
                                              Binary_MaxSensSpec = thh[[2]], Binary_TenPerc = thh[[3]])
          }else {
            ras[[sp]]$imputed_prediction <- stack(ras[[sp]]$imputed_prediction,
                                                  stack(lapply(suitability_final, "[[",
                                                               2)), thh[[1]], thh[[2]], thh[[3]])
          }
        }else {
          ras[[sp]]$imputed_prediction <- mapply(function(a,
                                                          b) {
            cc <- cbind(a, b[[2]])
            cc
          }, a = ras[[sp]]$imputed_prediction, b = suitability_final,
          SIMPLIFY = FALSE)
        }
        ras_final <- list(call = "imputed_prediction",
                          imputed_prediction = ras[[sp]]$imputed_prediction)
      }
      return(ras_final)
    })
    names(ras_suitability) <- names(ras)
    ras <- ras_suitability
  }
  return(ras)
}



