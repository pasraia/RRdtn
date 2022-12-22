ecospat.ESM.EnsembleModeling_CUSTOM<-function (ESM.modeling.output, weighting.score, threshold = NULL, 
          models) 
{
  require(PresenceAbsence)
  require(gtools)
  if (!weighting.score %in% c("AUC", "TSS", "Boyce", 
                              "Kappa", "SomersD")) {
    stop("weighting score not supported! Choose one of the following: AUC, TSS, Boyce, Kappa or SomersD")
  }
  iniwd <- getwd()
  on.exit(setwd(iniwd))
  setwd(ESM.modeling.output$wd)
  data <- ESM.modeling.output$data
  models. <- ESM.modeling.output$models.
  NbRunEval <- ESM.modeling.output$NbRunEval
  models <- ESM.modeling.output$models
  calib.lines <- as.data.frame(ESM.modeling.output$calib.lines)
  ESM_Projection <- ESM.modeling.output$ESM_Projection
  new.env <- ESM.modeling.output$new.env
  for (i in 1:length(models.)) load(models.[i])
  models. <- NULL
  for (n in 1:length(ESM.modeling.output$modeling.id)) {
    xx<-grep(ESM.modeling.output$modeling.id[n], 
             grep(paste(".", ESM.modeling.output$modeling.id[n], 
                        ".models.out", sep = ""), ls(), value = TRUE, 
                  fixed = TRUE), value = TRUE)
    yy<-as.numeric(gsub("ESM.BIOMOD.", "", unlist(strsplit(xx, paste(".", ESM.modeling.output$modeling.id[n], 
                                                                     ".models.out", sep = "")))))
    
    models. <- c(models., xx[order(yy)])
  }
  mymodel <- list()
  for (i in 1:length(models.)) mymodel[[i]] <- get(models.[i])
  weights <- unlist(lapply(mymodel, function(x) {
    y <- x
    if (weighting.score == "Boyce") {
      z <- get_predictions(y, as.data.frame = TRUE)
      z <- z[, grep(paste("RUN", NbRunEval + 1, sep = ""), 
                    colnames(z), invert = TRUE)]
      x <- get_evaluations(y)[, "Testing.data", , 
                              , ]
      if (length(models) > 1) {
        x[, ] <- NA
        x <- x[, colnames(x) != "Full" & colnames(x) != 
                 paste("RUN", NbRunEval + 1, sep = "")]
      }
      else {
        x[] <- NA
        x <- x[names(x) != "Full" & names(x) != 
                 paste("RUN", NbRunEval + 1, sep = "")]
      }
      for (n in 1:(length(models))) {
        model <- models[n]
        if (!TRUE %in% c(grepl("Full", y@models.failed), 
                         grepl(paste("RUN", NbRunEval + 1, sep = ""), 
                               y@models.failed))) {
          for (i in 1:NbRunEval) {
            if (sum(is.na(z[, grep(paste("RUN", 
                                         i, "_", model, sep = ""), colnames(z))])) != 
                nrow(z)) {
              if (length(models) > 1) {
                if (grepl("_", names(calib.lines))) {
                  x[rownames(x) == model, colnames(x) == 
                      paste("RUN", i, sep = "")] <- ecospat.boyce(z[!calib.lines[, 
                                                                                 paste("_RUN", i, sep = "")], 
                                                                    grep(paste("RUN", i, "_", 
                                                                               model, sep = ""), colnames(z))], 
                                                                  z[!calib.lines[, paste("_RUN", 
                                                                                         i, sep = "")] & data@data.species == 
                                                                      1, grep(paste("RUN", i, "_", 
                                                                                    model, sep = ""), colnames(z))], 
                                                                  PEplot = F)$cor
                }
                else {
                  x[rownames(x) == model, colnames(x) == 
                      paste("RUN", i, sep = "")] <- ecospat.boyce(z[!calib.lines[, 
                                                                                 paste("RUN", i, sep = "")], 
                                                                    grep(paste("RUN", i, "_", 
                                                                               model, sep = ""), colnames(z))], 
                                                                  z[!calib.lines[, paste("RUN", 
                                                                                         i, sep = "")] & data@data.species == 
                                                                      1, grep(paste("RUN", i, "_", 
                                                                                    model, sep = ""), colnames(z))], 
                                                                  PEplot = F)$cor
                }
              }
              else {
                if (grepl("_", names(calib.lines))) {
                  x[names(x) == paste("RUN", i, 
                                      sep = "")] <- ecospat.boyce(z[!calib.lines[, 
                                                                                 paste("_RUN", i, sep = "")], 
                                                                    grep(paste("RUN", i, "_", 
                                                                               model, sep = ""), colnames(z))], 
                                                                  z[!calib.lines[, paste("_RUN", 
                                                                                         i, sep = "")] & data@data.species == 
                                                                      1, grep(paste("RUN", i, "_", 
                                                                                    model, sep = ""), colnames(z))], 
                                                                  PEplot = F)$cor
                }
                else {
                  x[names(x) == paste("RUN", i, 
                                      sep = "")] <- ecospat.boyce(z[!calib.lines[, 
                                                                                 paste("RUN", i, sep = "")], 
                                                                    grep(paste("RUN", i, "_", 
                                                                               model, sep = ""), colnames(z))], 
                                                                  z[!calib.lines[, paste("RUN", 
                                                                                         i, sep = "")] & data@data.species == 
                                                                      1, grep(paste("RUN", i, "_", 
                                                                                    model, sep = ""), colnames(z))], 
                                                                  PEplot = F)$cor
                }
              }
            }
          }
        }
      }
      if (length(models) > 1) {
        x <- round(apply(x, 1, mean, na.rm = TRUE), 4)
      }
      else {
        x <- round(mean(x, na.rm = TRUE), 4)
      }
    }
    else {
      x <- get_evaluations(y)[, "Testing.data", , 
                              , ]
      if (length(models) > 1) {
        for (row in models) {
          if (ncol(calib.lines) == NbRunEval + 1) {
            if (is.na(x[row, NbRunEval + 1])) {
              x[row, ] <- NA
            }
          }
        }
        x <- x[, colnames(x) != "Full" & colnames(x) != 
                 paste("RUN", NbRunEval + 1, sep = "")]
        if (NbRunEval > 1) {
          x <- round(apply(x, 1, mean, na.rm = TRUE), 
                     4)
        }
      }
      else {
        if (ncol(calib.lines) == NbRunEval + 1) {
          if (is.na(x[NbRunEval + 1])) {
            x <- NA
          }
        }
        x <- x[names(x) != "Full" & names(x) != 
                 paste("RUN", NbRunEval + 1, sep = "")]
        x <- round(mean(x, na.rm = TRUE), 4)
      }
      if (weighting.score == "SomersD") {
        x <- x * 2 - 1
      }
    }
    names(x) <- paste(names(x), y@sp.name, sep = ".")
    return(x)
  }), recursive = TRUE)
  failed.mod <- "none"
  failed.mod <- lapply(mymodel, function(x) {
    return(if (x@models.failed[1] == "none") {
      paste(x@sp.name, "none", sep = ": ")
    } else {
      x@models.failed
    })
  })
  failed <- "none"
  if (sum(is.na(weights)) > 0) {
    warning(cat(paste("\n\n\n################################\n", 
                      paste("The following bivariate model(s) failed and are not included for building the ESMs\n", 
                            sep = " "))), print(names(which(is.na(weights)))))
    failed <- names(which(is.na(weights)))
  }
  if (is.null(threshold)) {
    if (weighting.score == "AUC") {
      threshold <- 0.5
    }
    else {
      threshold <- 0
    }
  }
  if (sum(weights <= threshold & !is.na(weights)) > 0) {
    warning(cat(paste("\n\n\n################################\n", 
                      paste("The following bivariate model(s) is (are) not included for building the ESMs\n because evaluation score is smaller or equal to the given threshold of", 
                            weighting.score, "=", threshold, ":\n", 
                            sep = " "))), print(names(weights[weights <= 
                                                                threshold & !is.na(weights)])))
  }
  weights[(weights <= threshold) | is.na(weights)] <- 0
  if (length(models) == 1) {
    mymodel <- mymodel[which(weights > 0)]
    weights <- weights[which(weights > 0)]
  }
  test.pred <- lapply(mymodel, function(x) {
    x <- get_predictions(x, as.data.frame = TRUE)
    return(x)
  })
  test.ESM <- NULL
  biva.st2 <- do.call(cbind, test.pred)
  for (i in 1:length(models)) {
    for (run in 1:ncol(calib.lines)) {
      if (length(models) > 1) {
        test.ESM1 <- apply(biva.st2[, grep(paste("RUN", 
                                                 run, "_", models[i], sep = ""), 
                                           colnames(biva.st2))], 1, function(x) weighted.mean(x, 
                                                                                              weights[grep(models[i], names(weights))], na.rm = TRUE))
      }
      else {
        test.ESM1 <- apply(biva.st2[, grep(paste("RUN", 
                                                 run, "_", models[i], sep = ""), 
                                           colnames(biva.st2))], 1, function(x) weighted.mean(x, 
                                                                                              weights, na.rm = TRUE))
      }
      test.ESM <- cbind(test.ESM, test.ESM1)
      colnames(test.ESM)[ncol(test.ESM)] <- paste("RUN", 
                                                  run, "_", models[i], sep = "")
    }
  }
  test.ESM <- as.data.frame(test.ESM)
  DATA <- cbind(1:length(data@data.species), resp.var = data@data.species, 
                test.ESM/1000)
  if (any(DATA$resp.var==1)&any(is.na(DATA$resp.var))) {
    DATA$resp.var[is.na(DATA$resp.var)] <- 0
  }
  EVAL <- NULL
  for (i in 1:NbRunEval) {
    DATA1 <- cbind(DATA[, 1:2], DATA[, grep(paste("RUN", 
                                                  i, "_", sep = ""), colnames(DATA))])
    colnames(DATA1)[3] <- colnames(DATA)[2 + i]
    if (length(models) > 1) {
      if (nrow(DATA1) %in% apply(DATA1[, 3:ncol(DATA1)], 
                                 2, function(x) {
                                   sum(is.na(x))
                                 })) {
        failed <- which(nrow(DATA1) == apply(DATA1, 2, 
                                             function(x) {
                                               sum(is.na(x))
                                             }))
        DATA1[, nrow(DATA1) == apply(DATA1, 2, function(x) {
          sum(is.na(x))
        })] <- 0
      }
    }
    else {
      if (nrow(DATA1) == sum(is.na(DATA1[, 3]))) {
        failed <- paste("RUN", i, sep = "")
        DATA1[, 3] <- 0
      }
    }
    EVAL1 <- presence.absence.accuracy(DATA1[!calib.lines[, 
                                                          i], ], threshold = as.vector(optimal.thresholds(DATA1[!calib.lines[, 
                                                                                                                             i], ], opt.methods = "MaxSens+Spec"), mode = "numeric")[-1])
    EVAL1 <- EVAL1[c(1, 2, 4:7, 9:12)]
    EVAL1$TSS <- EVAL1$sensitivity + EVAL1$specificity - 
      1
    if (length(models) > 1) {
      if (nrow(DATA1) %in% apply(DATA1[, 3:ncol(DATA1)], 
                                 2, function(x) {
                                   sum(is.na(x))
                                 })) {
        EVAL1[EVAL1$model == names(failed), 2:11] <- NA
      }
    }
    else {
      if (nrow(DATA1) == sum(is.na(DATA1[, 3]))) {
        EVAL1[EVAL1$model == names(failed), 2:11] <- NA
      }
    }
    EVAL1$SomersD <- EVAL1$AUC * 2 - 1
    EVAL1$Boyce <- EVAL1$MPA <- NA
    for (n in 1:nrow(EVAL1)) {
      EVAL1$MPA[EVAL1$model == EVAL1$model[n]] <- ecospat.mpa(DATA1[!calib.lines[, 
                                                                                 i] & DATA1[, 2] == 1, EVAL1$model[n]])
      EVAL1$Boyce[EVAL1$model == EVAL1$model[n]] <- ecospat.boyce(DATA1[!calib.lines[, 
                                                                                     i], EVAL1$model[n]], DATA1[!calib.lines[, i] & 
                                                                                                                  DATA1[, 2] == 1, EVAL1$model[n]], PEplot = F)$cor
    }
    EVAL1$technique <- unlist(strsplit(EVAL1$model, split = "_"))[seq(2, 
                                                                      nrow(EVAL1) * 2, 2)]
    EVAL1$RUN <- paste("RUN", i, sep = "")
    EVAL <- rbind(EVAL, EVAL1)
  }
  if (length(models) > 1) {
    weights.double <- aggregate(EVAL[, weighting.score], 
                                by = list(EVAL$technique), FUN = mean, na.rm = TRUE)
    weights.double <- weights.double[order(weights.double[, 
                                                          1]), ]
    if (!is.null(threshold)) {
      weights.double[weights.double <= threshold] <- 0
      if (sum(weights.double == 0) > 0) {
        warning(cat(paste("\n\n\n################################\n", 
                          paste("The following ESM model(s) is (are) not included for building the 'double ensemble'\n because evaluation score is smaller or equal to the given threshold:\n", 
                                list(weights.double[weights.double[, 2] == 
                                                      0, 1]), sep = " ", "\n################################\n"))))
      }
    }
    for (run in 1:ncol(calib.lines)) {
      test.ESM$EF <- apply(test.ESM[, grep(paste("RUN", 
                                                 run, "_", sep = ""), colnames(test.ESM))], 
                           1, function(x) weighted.mean(x, weights.double[, 
                                                                          2], na.rm = TRUE))
      colnames(test.ESM)[ncol(test.ESM)] <- paste("RUN", 
                                                  run, "_EF", sep = "")
    }
    DATA <- cbind(1:length(data@data.species), resp.var = data@data.species, 
                  test.ESM/1000)
    if (any(DATA$resp.var==1)&any(is.na(DATA$resp.var))) {
      DATA$resp.var[is.na(DATA$resp.var)] <- 0
    }
    EVAL <- NULL
    for (i in 1:NbRunEval) {
      DATA1 <- cbind(DATA[, 1:2], DATA[, grep(paste("RUN", 
                                                    i, "_", sep = ""), colnames(DATA))])
      if (nrow(DATA1) %in% apply(DATA1[, 3:ncol(DATA1)], 
                                 2, function(x) {
                                   sum(is.na(x))
                                 })) {
        failed <- which(nrow(DATA1) == apply(DATA1, 2, 
                                             function(x) {
                                               sum(is.na(x))
                                             }))
        DATA1[, nrow(DATA1) == apply(DATA1, 2, function(x) {
          sum(is.na(x))
        })] <- 0
      }
      EVAL1 <- presence.absence.accuracy(DATA1[!calib.lines[, 
                                                            i], ], threshold = as.vector(optimal.thresholds(DATA1[!calib.lines[, 
                                                                                                                               i], ], opt.methods = "MaxSens+Spec"), mode = "numeric")[-1])
      EVAL1 <- EVAL1[c(1, 2, 4:7, 9:12)]
      EVAL1$TSS <- EVAL1$sensitivity + EVAL1$specificity - 
        1
      if (nrow(DATA1) %in% apply(DATA1[, 3:ncol(DATA1)], 
                                 2, function(x) {
                                   sum(is.na(x))
                                 })) {
        EVAL1[EVAL1$model == names(failed), 2:11] <- NA
      }
      EVAL1$SomersD <- EVAL1$AUC * 2 - 1
      EVAL1$Boyce <- EVAL1$MPA <- NA
      for (n in 1:nrow(EVAL1)) {
        EVAL1$MPA[EVAL1$model == EVAL1$model[n]] <- ecospat.mpa(DATA1[!calib.lines[, 
                                                                                   i] & DATA1[, 2] == 1, EVAL1$model[n]])
        EVAL1$Boyce[EVAL1$model == EVAL1$model[n]] <- ecospat.boyce(DATA1[!calib.lines[, 
                                                                                       i], EVAL1$model[n]], DATA1[!calib.lines[, i] & 
                                                                                                                    DATA1[, 2] == 1, EVAL1$model[n]], PEplot = F)$cor
      }
      EVAL1$technique <- unlist(strsplit(EVAL1$model, split = "_"))[seq(2, 
                                                                        nrow(EVAL1) * 2, 2)]
      EVAL1$RUN <- paste("RUN", i, sep = "")
      EVAL <- rbind(EVAL, EVAL1)
    }
  }
  colnames(DATA) <- gsub(paste("RUN", NbRunEval + 1, 
                               sep = ""), "Full", colnames(DATA))
  colnames(DATA)[3:ncol(DATA)] <- paste(colnames(DATA)[3:ncol(DATA)], 
                                        "ESM", sep = "_")
  DATA[, -c(1, 2)] <- DATA[, -c(1, 2)] * 1000
  EVAL[, 2:14] <- round(EVAL[, 2:14], 3)
  EVAL[, 2] <- EVAL[, 2] * 1000
  if (length(models) > 1) {
    output <- list(species = data@sp.name, ESM.fit = round(DATA[, 
                                                                -1]), ESM.evaluations = EVAL, weights = weights, 
                   weights.EF = weights.double, failed = failed.mod)
  }
  else {
    output <- list(species = data@sp.name, ESM.fit = round(DATA[, 
                                                                -1]), ESM.evaluations = EVAL, weights = weights, 
                   failed = failed.mod)
  }
  save(output, file = paste("ESM_EnsembleModeling..", 
                            weighting.score, threshold, ESM.modeling.output$modeling.id, 
                            "out", sep = "."))
  setwd(iniwd)
  return(output)
}