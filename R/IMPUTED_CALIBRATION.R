#' @importFrom stats median aggregate cov mahalanobis median pchisq
#' @importFrom RRphylo swapONE RRphylo
#' @importFrom ape Ntip
#' @importFrom pbapply pblapply
#' @importFrom gtools mixedorder
#' @importFrom Rphylopars phylopars

IMPUTED_CALIBRATION<-function (ENFA_output,
                               tree,
                               nsim = 10,
                               si = 0.2,
                               si2 = 0.2,
                               clust,
                               evaluate = FALSE,
                               boot_test_perc = 20,
                               boot_reps = 20,
                               output_options,
                               eval_threshold,
                               eval.metric)
{

  all_models_NN <- ENFA_output[!sapply(ENFA_output, function(k) is.null(k$calibrated_model))]

  if (!is.null(clust))
    parallel = TRUE
  else parallel = FALSE
  all_M <- all_models_NN[which(sapply(all_models_NN, function(x) x$call) ==
                                 "calibrated_enfa")]
  axx <- sapply(all_M, function(y) {
    y$calibrated_model$full_model$significant_axes
  })
  matri_test <- lapply(all_M, function(x) {
    ax <- x$calibrated_model$full_model$co[, 1:median(axx)]
    ax
  })
  temp <- lapply(1:nrow(matri_test[[1]]), function(x) {
    vv <- lapply(matri_test, function(k) {
      vv <- as.data.frame(t(k[x, ]))
      rownames(vv) <- rownames(k)[x]
      vv
    })
    do.call(rbind, vv)
  })
  names(temp) <- rownames(matri_test[[1]])
  tree1 <- tree
  if (parallel) {
    cl <- makeCluster(clust, type = "SOCK")
    registerDoParallel(cl)
    output <- foreach(i = 1:nsim, .packages = c("RRphylo",
                                                "Rphylopars", "ape")) %dopar% {
                                                  tree <- swapONE(tree1, si, si2)[[1]]
                                                  tree$edge.length[which(tree$edge.length < max(diag(vcv(tree)))/10000)] <- tree$edge.length[which(tree$edge.length <
                                                                                                                                                     max(diag(vcv(tree)))/10000)] + max(diag(vcv(tree)))/10000
                                                  RR_sh <- lapply(1:length(temp), function(x) {
                                                    m <- matrix(NA, ncol = ncol(temp[[x]]), nrow = sum(!tree$tip.label %in%
                                                                                                         rownames(temp[[x]])))
                                                    colnames(m) <- colnames(temp[[x]])
                                                    rownames(m) <- tree$tip.label[!tree$tip.label %in%
                                                                                    rownames(temp[[x]])]
                                                    temp[[x]] <- rbind(temp[[x]], m)
                                                    xx <- cbind(species = rownames(temp[[x]]), temp[[x]])
                                                    xx <- xx[match(tree$tip.label, xx$species), ]
                                                    rownames(xx) <- NULL
                                                    ress <- phylopars(xx, tree)
                                                    y <- ress$anc_recon[1:Ntip(tree), 1:(ncol(xx) -
                                                                                           1)]
                                                    RR_sh <- RRphylo(tree, y)
                                                    RR_sh
                                                  })
                                                  res_phylo <- lapply(RR_sh, function(x) x$predicted.phenotype)
                                                  names(res_phylo) <- rownames(matri_test[[1]])
                                                  dat_tot <- lapply(1:nrow(res_phylo[[1]]), function(x) {
                                                    bb <- lapply(res_phylo, function(y) {
                                                      y[x, ]
                                                    })
                                                    do.call(rbind, bb)
                                                  })
                                                  names(dat_tot) <- rownames(res_phylo[[1]])
                                                  data_final <- dat_tot[!names(dat_tot) %in% names(matri_test)]

                                                }
    stopCluster(cl)
    closeAllConnections()
    gc()
    output <- unlist(output, recursive = FALSE)
    output <- split(output, names(output))
    output <- lapply(output, function(f) {
      names(f) <- NULL
      c(list(call = "calibrated_imputed"), list(co = f))
    })
  }
  if (!parallel) {
    output <- replicate(nsim, {
      tree <- swapONE(tree1, si, si2)[[1]]
      RR_sh <- lapply(1:length(temp), function(x) {
        m <- matrix(NA, ncol = ncol(temp[[x]]), nrow = sum(!tree$tip.label %in%
                                                             rownames(temp[[x]])))
        colnames(m) <- colnames(temp[[x]])
        rownames(m) <- tree$tip.label[!tree$tip.label %in%
                                        rownames(temp[[x]])]
        temp[[x]] <- rbind(temp[[x]], m)
        xx <- cbind(species = rownames(temp[[x]]), temp[[x]])
        xx <- xx[match(tree$tip.label, xx$species), ]
        rownames(xx) <- NULL
        ress <- phylopars(xx, tree)
        y <- ress$anc_recon[1:Ntip(tree), 1:(ncol(xx) -
                                               1)]
        RR_sh <- RRphylo(tree, y)
        RR_sh
      })
      res_phylo <- lapply(RR_sh, function(x) x$predicted.phenotype)
      names(res_phylo) <- rownames(matri_test[[1]])
      dat_tot <- lapply(1:nrow(res_phylo[[1]]), function(x) {
        bb <- lapply(res_phylo, function(y) {
          y[x, ]
        })
        do.call(rbind, bb)
      })
      names(dat_tot) <- rownames(res_phylo[[1]])
      data_final <- dat_tot[!names(dat_tot) %in% names(matri_test)]
    }, simplify = FALSE)

    gc()

    output <- unlist(output, recursive = FALSE)
    output <- split(output, names(output))
    output <- lapply(output, function(f) {
      names(f) <- NULL
      c(list(call = "calibrated_imputed"), list(co = f))
    })
  }
  if (evaluate) {
    cat(paste("\n", "EVALUATING IMPUTATION", "\n"))
    tt <- ENFA_output[match(names(output), names(ENFA_output))]
    tt <- mapply(function(xx, yy) {
      yy$calibrated_model$co <- xx$co
      yy
    }, xx = output, yy = tt, SIMPLIFY = FALSE)
    EVAL_ALL <- lapply(tt, function(xx) {
      aa <- xx$formatted_data$input_ones
      if (nrow(aa) >= 10) {
        ss <- replicate(boot_reps, {
          sample(1:nrow(aa), nrow(aa) * boot_test_perc/100)
        }, simplify = FALSE)
        b_train <- lapply(ss, function(bb) aa[-bb, ])
        b_test <- lapply(ss, function(bb) aa[bb, , drop = FALSE])
        if (parallel) {
          cl <- makeCluster(clust, type = "SOCK")
          registerDoParallel(cl)
          s = NULL
          EVAL <- foreach(s = 1:nsim, .packages = c("dismo",
                                                    "ecospat", "PresenceAbsence"), .export = c("mahasuhab.custom",
                                                                                               "aa", "boot_reps")) %dopar% {
                                                                                                 U <- xx$calibrated_model$co[[s]]
                                                                                                 f1 <- function(y) y %*% U
                                                                                                 EVAL <- mapply(function(train, test) {
                                                                                                   if (is.null(xx$formatted_data$time_col)) {
                                                                                                     train$globalID <- train[, xx$formatted_data$geoID_col]
                                                                                                   }
                                                                                                   else {
                                                                                                     train$globalID <- paste(train[, xx$formatted_data$geoID_col],
                                                                                                                             train[, xx$formatted_data$time_col],
                                                                                                                             sep = "_")
                                                                                                   }
                                                                                                   if (is.null(xx$formatted_data$time_col)) {
                                                                                                     test$globalID <- test[, xx$formatted_data$geoID_col]
                                                                                                   }
                                                                                                   else {
                                                                                                     test$globalID <- paste(test[, xx$formatted_data$geoID_col],
                                                                                                                            test[, xx$formatted_data$time_col],
                                                                                                                            sep = "_")
                                                                                                   }
                                                                                                   ras_back <- t(apply(xx$formatted_data$input_back[,
                                                                                                                                                    rownames(U)], 1, f1))
                                                                                                   gg <- ras_back <- as.data.frame(ras_back)
                                                                                                   ras_back[, xx$formatted_data$geoID_col] <- xx$formatted_data$input_back[,
                                                                                                                                                                           xx$formatted_data$geoID_col]
                                                                                                   if (!is.null(xx$formatted_data$time_col)) {
                                                                                                     ras_back[, xx$formatted_data$time_col] <- xx$formatted_data$input_back[,
                                                                                                                                                                            xx$formatted_data$time_col]
                                                                                                     ras_back$globalID <- paste(ras_back[,
                                                                                                                                         xx$formatted_data$geoID_col], ras_back[,
                                                                                                                                                                                xx$formatted_data$time_col], sep = "_")
                                                                                                   }
                                                                                                   else {
                                                                                                     ras_back$globalID <- ras_back[, xx$formatted_data$geoID_col]
                                                                                                   }
                                                                                                   PRED_sites_train <- ras_back[which(ras_back$globalID %in%
                                                                                                                                        train$globalID), ]
                                                                                                   PRED_sites_test <- ras_back[which(ras_back$globalID %in%
                                                                                                                                       test$globalID), ]
                                                                                                   a <- ras_back
                                                                                                   if (nrow(a) > 10000) {
                                                                                                     k <- sample(1:nrow(a), 10000)
                                                                                                     a <- a[k, ]
                                                                                                   }
                                                                                                   else k <- 1:nrow(a)
                                                                                                   a_scaled <- dudi.pca(gg, scannf = FALSE)$tab
                                                                                                   a_scaled$globalID <- ras_back$globalID
                                                                                                   p_scaled <- a_scaled[which(a_scaled$globalID %in%
                                                                                                                                PRED_sites_train$globalID), ]
                                                                                                   maha_prob <- mahasuhab.custom(a_scaled[,
                                                                                                                                          1:ncol(gg)], p_scaled[, 1:ncol(gg)])
                                                                                                   maha_prob$globalID <- ras_back$globalID
                                                                                                   fit <- maha_prob$MD[k]
                                                                                                   obs <- maha_prob[which(maha_prob$globalID %in%
                                                                                                                            PRED_sites_test$globalID), ]$MD
                                                                                                   DDATA <- data.frame(ID = 1:(length(fit) +
                                                                                                                                 length(obs)), OBS = c(rep(1, length(obs)),
                                                                                                                                                       rep(0, length(fit))), PRED = c(obs, fit))
                                                                                                   th <- optimal.thresholds(DDATA)[3, 2]
                                                                                                   ev0 <- presence.absence.accuracy(DDATA,
                                                                                                                                    th)
                                                                                                   TPR <- ev0$sensitivity
                                                                                                   TNR <- ev0$specificity
                                                                                                   FPR <- 1 - TPR
                                                                                                   FNR <- 1 - TNR
                                                                                                   TSS <- TPR + TNR - 1
                                                                                                   AUC <- ev0$AUC
                                                                                                   SORENSEN <- 2 * TPR/(FNR + (2 * TPR) +
                                                                                                                          FPR)
                                                                                                   CBI <- ecospat.boyce(DDATA[which(DDATA$OBS ==
                                                                                                                                      0), ]$PRED, DDATA[which(DDATA$OBS ==
                                                                                                                                                                1), ]$PRED, PEplot = F)$cor
                                                                                                   maha_prob_train <- maha_prob[which(maha_prob$globalID %in%
                                                                                                                                        PRED_sites_train$globalID), ]$MD
                                                                                                   maha_prob_test <- obs
                                                                                                   n_vals <- length(maha_prob_test)
                                                                                                   suit_at_th <- quantile(maha_prob_train,
                                                                                                                          0.1)
                                                                                                   omr <- 1 - length(which(maha_prob_test >=
                                                                                                                             suit_at_th))/n_vals
                                                                                                   return(c(AUC = AUC, TSS = TSS, CBI = CBI,
                                                                                                            SORENSEN = SORENSEN, OMR = omr))
                                                                                                 }, train = b_train, test = b_test, SIMPLIFY = FALSE)
                                                                                               }
          stopCluster(cl)
          closeAllConnections()
          gc()
          EVAL_ALL <- lapply(1:length(EVAL), function(x) {
            dd <- as.data.frame(do.call(rbind, EVAL[[x]]))
            dd$fold <- paste("fold", 1:nrow(dd), sep = "_")
            dd$swap <- paste("swap", x, sep = "_")
            dd
          })
          EVAL_ALL <- do.call(rbind, EVAL_ALL)
          return(EVAL_ALL)
        }
        else {
          EVAL_ALL <- pblapply(xx$calibrated_model$co,
                               function(co) {
                                 U <- co
                                 f1 <- function(y) y %*% U
                                 EVAL <- mapply(function(train, test) {
                                   if (is.null(xx$formatted_data$time_col)) {
                                     train$globalID <- train[, xx$formatted_data$geoID_col]
                                   }
                                   else {
                                     train$globalID <- paste(train[, xx$formatted_data$geoID_col],
                                                             train[, xx$formatted_data$time_col],
                                                             sep = "_")
                                   }
                                   if (is.null(xx$formatted_data$time_col)) {
                                     test$globalID <- test[, xx$formatted_data$geoID_col]
                                   }
                                   else {
                                     test$globalID <- paste(test[, xx$formatted_data$geoID_col],
                                                            test[, xx$formatted_data$time_col],
                                                            sep = "_")
                                   }
                                   ras_back <- t(apply(xx$formatted_data$input_back[,
                                                                                    rownames(U)], 1, f1))
                                   gg <- ras_back <- as.data.frame(ras_back)
                                   ras_back[, xx$formatted_data$geoID_col] <- xx$formatted_data$input_back[,
                                                                                                           xx$formatted_data$geoID_col]
                                   if (!is.null(xx$formatted_data$time_col)) {
                                     ras_back[, xx$formatted_data$time_col] <- xx$formatted_data$input_back[,
                                                                                                            xx$formatted_data$time_col]
                                     ras_back$globalID <- paste(ras_back[,
                                                                         xx$formatted_data$geoID_col], ras_back[,
                                                                                                                xx$formatted_data$time_col], sep = "_")
                                   }
                                   else {
                                     ras_back$globalID <- ras_back[, xx$formatted_data$geoID_col]
                                   }
                                   PRED_sites_train <- ras_back[which(ras_back$globalID %in%
                                                                        train$globalID), ]
                                   PRED_sites_test <- ras_back[which(ras_back$globalID %in%
                                                                       test$globalID), ]
                                   a <- ras_back
                                   if (nrow(a) > 10000) {
                                     k <- sample(1:nrow(a), 10000)
                                     a <- a[k, ]
                                   }
                                   else k <- 1:nrow(a)
                                   a_scaled <- dudi.pca(gg, scannf = FALSE)$tab
                                   a_scaled$globalID <- ras_back$globalID
                                   p_scaled <- a_scaled[which(a_scaled$globalID %in%
                                                                PRED_sites_train$globalID), ]
                                   maha_prob <- mahasuhab.custom(a_scaled[,
                                                                          1:ncol(gg)], p_scaled[, 1:ncol(gg)])
                                   maha_prob$globalID <- ras_back$globalID
                                   fit <- maha_prob$MD[k]
                                   obs <- maha_prob[which(maha_prob$globalID %in%
                                                            PRED_sites_test$globalID), ]$MD
                                   DDATA <- data.frame(ID = 1:(length(fit) +
                                                                 length(obs)), OBS = c(rep(1, length(obs)),
                                                                                       rep(0, length(fit))), PRED = c(obs,
                                                                                                                      fit))
                                   th <- optimal.thresholds(DDATA)[3, 2]
                                   ev0 <- presence.absence.accuracy(DDATA,
                                                                    th)
                                   TPR <- ev0$sensitivity
                                   TNR <- ev0$specificity
                                   FPR <- 1 - TPR
                                   FNR <- 1 - TNR
                                   TSS <- TPR + TNR - 1
                                   AUC <- ev0$AUC
                                   SORENSEN <- 2 * TPR/(FNR + (2 * TPR) +
                                                          FPR)
                                   CBI <- ecospat.boyce(DDATA[which(DDATA$OBS ==
                                                                      0), ]$PRED, DDATA[which(DDATA$OBS ==
                                                                                                1), ]$PRED, PEplot = F)$cor
                                   maha_prob_train <- maha_prob[which(maha_prob$globalID %in%
                                                                        PRED_sites_train$globalID), ]$MD
                                   maha_prob_test <- obs
                                   n_vals <- length(maha_prob_test)
                                   suit_at_th <- quantile(maha_prob_train,
                                                          0.1)
                                   omr <- 1 - length(which(maha_prob_test >=
                                                             suit_at_th))/n_vals
                                   return(c(AUC = AUC, TSS = TSS, CBI = CBI,
                                            SORENSEN = SORENSEN, OMR = omr))
                                 }, train = b_train, test = b_test, SIMPLIFY = FALSE)
                                 EVAL <- as.data.frame(do.call(rbind, EVAL))
                                 EVAL$fold <- paste("fold", 1:nrow(EVAL),
                                                    sep = "_")
                                 EVAL
                               })
          EVAL_ALL <- as.data.frame(do.call(rbind, EVAL_ALL))
          swap <- paste("swap", 1:length(xx$calibrated_model$co),
                        sep = "_")
          EVAL_ALL$swap <- unlist(lapply(swap, rep, boot_reps))
          EVAL_ALL
        }
      }
      if (nrow(aa) < 10)
        EVAL_ALL <- "NULL"
      return(EVAL_ALL)
    })
    EVAL_ALL <- suppressWarnings(lapply(EVAL_ALL, function(x) {
      if (length(x) > 1) {
        aggregate(cbind(AUC, TSS, CBI, SORENSEN, OMR) ~
                    swap, data = x, FUN = "mean")
      }
      else NULL
    }))
    EVAL_ALL <- lapply(EVAL_ALL, function(zz) {
      if (!is.null(zz)) {
        zz[mixedorder(zz[, 1]), ]
      }
      else NULL
    })
    output_opt <- rep(output_options, length(EVAL_ALL))
    if (output_options == "full") {
      EVAL_ALL <- EVAL_ALL
    }
    if (output_options == "best") {
      EVAL <- lapply(EVAL_ALL, function(x) {
        if (!is.null(x)) {
          if (eval.metric != "OMR") {
            xx <- x[which(x[, eval.metric] >= eval_threshold),
            ]
            xx[which.max(xx[, eval.metric]), ]
          }
          else {
            xx <- x[which(x[, eval.metric] <= eval_threshold),
            ]
            xx[which.min(xx[, eval.metric]), ]
          }
        }
        else NULL
      })
      EVAL_ALL <- lapply(1:length(EVAL), function(p) {
        if (is.null(EVAL[[p]]) || nrow(EVAL[[p]]) ==
            0) {
          warning(paste(names(EVAL_ALL[p]), "model validation never exceeds the threshold value. ENphylo will return the full model",
                        sep = " "))
          output_opt[p] <<- "full"
          EVAL_ALL[[p]] <- EVAL_ALL[[p]]
        }
        else EVAL_ALL[[p]] <- EVAL[[p]]
      })
    }
    if (output_options == "weighted.mean") {
      EVAL <- lapply(EVAL_ALL, function(x) {
        if (!is.null(x)) {
          if (eval.metric != "OMR") {
            xx <- x[which(x[, eval.metric] >= eval_threshold),
            ]
          }
          else {
            xx <- x[which(x[, eval.metric] <= eval_threshold),
            ]
          }
          xx
        }
        else NULL
      })
      EVAL_ALL <- lapply(1:length(EVAL), function(p) {
        if (is.null(EVAL[[p]]) || nrow(EVAL[[p]]) ==
            0) {
          warning(paste(names(EVAL_ALL[p]), "model validation never exceeds the threshold value. ENphylo will return the full model",
                        sep = " "))
          output_opt[p] <<- "full"
          EVAL_ALL[[p]] <- EVAL_ALL[[p]]
        }
        else EVAL_ALL[[p]] <- EVAL[[p]]
      })
    }
    output <- mapply(function(x, y) {
      if (!is.null(y))
        x$co <- x$co[as.numeric(gsub("swap_", "", y$swap))]
      return(x)
    }, x = output, y = EVAL_ALL, SIMPLIFY = F)
  }
  else {
    EVAL_ALL <- output_opt <- replicate(length(output), NULL,
                                        simplify = FALSE)
  }
  output <- mapply(function(jj, kk, ll, nam) {
    if (!is.null(kk)) {
      jj$evaluation <- kk
      rownames(jj$evaluation) <- jj$evaluation$swap
      jj$evaluation <- jj$evaluation[, -grep("swap", colnames(jj$evaluation))]
    }
    if (evaluate & !is.null(kk)) {
      jj$output_options <- c(ll, eval.metric)
    }
    else {
      jj$output_options = "full"
      jj <- c(jj, list(evaluation = NULL))
    }
    dir.create(file.path("ENphylo_imputed_models", nam),
               recursive = TRUE)
    model_outputs<- list(jj)
    names(model_outputs) <- nam
    save(model_outputs, file = paste0("ENphylo_imputed_models/", nam,
                                      "/output_model.RData"))
    jj
  }, jj = output, kk = EVAL_ALL, ll = output_opt, nam = names(output),
  SIMPLIFY = FALSE)
  return(output)
}
