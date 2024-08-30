#' @importFrom foreach foreach %dopar%
#' @importFrom parallel makeCluster stopCluster
#' @importFrom ade4 dudi.pca
#' @importFrom dismo mahal
#' @importFrom PresenceAbsence optimal.thresholds presence.absence.accuracy
#' @importFrom ecospat ecospat.boyce
#' @importFrom doParallel registerDoParallel

ENFA_CALIBRATION<-function (formatted_data,
                            boot_test_perc = 20,
                            boot_reps = 20,
                            sig_axes_selection = c("brStick", "user.defined"),
                            user_defined_sig_axes = NULL,
                            clust){
  if (!is.null(clust))
    parallel = TRUE else parallel=FALSE
  CENFA_brStick <- function(eigs) {
    if (max(Im(eigs)) > 1e-05)
      stop("broken-stick method does not work for complex eigenvalues")
    eigs <- Re(eigs)
    p <- length(eigs)
    a <- NULL
    r <- NULL
    for (j in 1:p) {
      a[j] <- 1/p * sum(1/(j:p))
      r[j] <- eigs[j]/(sum(eigs))
    }
    length(which(r > a))
  }
  if (sig_axes_selection == "user.defined" & is.null(user_defined_sig_axes))
    stop("Please, specify a number of ENFA axes to retain")
  FULL_MODEL <- enfa.custom(x = formatted_data$input_back[,!grepl(formatted_data$obs_col, colnames(formatted_data$input_back))],
                            s.dat = formatted_data$input_ones[, grepl(paste(c(formatted_data$obs_col,
                                                                              formatted_data$time_col, formatted_data$geoID_col),
                                                                            collapse = "|"), colnames(formatted_data$input_ones))],
                            obs_col = formatted_data$obs_col, geoID_col = formatted_data$geoID_col,
                            time_col = formatted_data$time_col)
  if (sig_axes_selection == "brStick")
    FULL_MODEL$significant_axes <- ifelse(CENFA_brStick(FULL_MODEL$sf) <
                                            2, 2, CENFA_brStick(FULL_MODEL$sf))
  if (sig_axes_selection == "user.defined")
    FULL_MODEL$significant_axes <- user_defined_sig_axes
  FULL_MODEL$ones_coords <- formatted_data$ones_coords
  ss <- replicate(boot_reps, {
    sample(1:nrow(formatted_data$input_ones), nrow(formatted_data$input_ones) *
             boot_test_perc/100)
  }, simplify = FALSE)
  if (boot_reps > 0) {
    if (parallel) {
      cl <- makeCluster(clust, type = "SOCK")
      registerDoParallel(cl)
      EVALUATIONS <- foreach(s = ss, .packages = c("dismo",
                                                   "ecospat", "PresenceAbsence")) %dopar%
        {
          s.dat.train <- formatted_data$input_ones[-s,
                                                   , drop = FALSE]
          if (is.null(formatted_data$time_col)) {
            s.dat.train$globalID <- s.dat.train[, formatted_data$geoID_col]
          }else {
            s.dat.train$globalID <- paste(s.dat.train[,
                                                      formatted_data$geoID_col], s.dat.train[,
                                                                                             formatted_data$time_col], sep = "_")
          }
          s.dat.test <- formatted_data$input_ones[s,
                                                  , drop = FALSE]
          if (is.null(formatted_data$time_col)) {
            s.dat.test$globalID <- s.dat.test[, formatted_data$geoID_col]
          }else {
            s.dat.test$globalID <- paste(s.dat.test[,
                                                    formatted_data$geoID_col], s.dat.test[,
                                                                                          formatted_data$time_col], sep = "_")
          }
          mod1 <- enfa.custom(x = formatted_data$input_back[,
                                                            !grepl(formatted_data$obs_col, colnames(formatted_data$input_back))],
                              s.dat = s.dat.train[, grepl(paste(c(formatted_data$obs_col,
                                                                  formatted_data$time_col, formatted_data$geoID_col),
                                                                collapse = "|"), colnames(s.dat.train))],
                              obs_col = formatted_data$obs_col, geoID_col = formatted_data$geoID_col,
                              time_col = formatted_data$time_col)
          if (sig_axes_selection == "brStick")
            sig <- ifelse(CENFA_brStick(mod1$sf) < 2,
                          2, CENFA_brStick(mod1$sf))
          if (sig_axes_selection == "user.defined")
            sig <- user_defined_sig_axes
          if (is.null(formatted_data$time_col)) {
            PRED <- ENphylo_prediction(object = list(mod1),
                                       newdata = formatted_data$input_back[, !grepl(paste(c(formatted_data$geoID_col,
                                                                                            formatted_data$obs_col), collapse = "|"),
                                                                                    colnames(formatted_data$input_back))])$enfa_prediction
          }else {
            PRED <- ENphylo_prediction(object = list(mod1),
                                       newdata = formatted_data$input_back[, !grepl(paste(c(formatted_data$geoID_col,
                                                                                            formatted_data$obs_col, formatted_data$time_col),
                                                                                          collapse = "|"), colnames(formatted_data$input_back))])$enfa_prediction
          }
          PRED[, formatted_data$geoID_col] <- formatted_data$input_back[,
                                                                        formatted_data$geoID_col]
          if (!is.null(formatted_data$time_col))
            PRED[, formatted_data$time_col] <- formatted_data$input_back[,
                                                                         formatted_data$time_col]
          if (is.null(formatted_data$time_col)) {
            PRED$globalID <- PRED[, formatted_data$geoID_col]
          }else {
            PRED$globalID <- paste(PRED[, formatted_data$geoID_col],
                                   PRED[, formatted_data$time_col],
                                   sep = "_")
          }
          PRED_sites_train <- PRED[which(PRED$globalID %in%
                                           s.dat.train$globalID), ]
          PRED_sites_test <- PRED[which(PRED$globalID %in%
                                          s.dat.test$globalID), ]
          a <- PRED
          if (nrow(a) > 10000) {
            k <- sample(1:nrow(a), 10000)
            a <- a[k, ]
          }else k <- 1:nrow(a)
          a_scaled <- dudi.pca(PRED[, 1:sig], scannf = FALSE)$tab
          a_scaled$globalID <- PRED$globalID
          p_scaled <- a_scaled[which(a_scaled$globalID %in%
                                       PRED_sites_train$globalID), ]
          maha_prob <- mahasuhab.custom(a_scaled[, 1:sig,
                                                 drop = FALSE], p_scaled[, 1:sig, drop = FALSE])
          maha_prob$globalID <- PRED$globalID
          fit <- maha_prob$MD[k]
          obs <- maha_prob[which(maha_prob$globalID %in%
                                   PRED_sites_test$globalID), ]$MD
          DDATA <- data.frame(ID = 1:(length(fit) + length(obs)),
                              OBS = c(rep(1, length(obs)), rep(0, length(fit))),
                              PRED = c(obs, fit))
          th <- optimal.thresholds(DDATA)[3, 2]

          ev0 <- presence.absence.accuracy(DDATA, th)
          TPR <- ev0$sensitivity
          TNR <- ev0$specificity
          FPR <- 1 - TPR
          FNR <- 1 - TNR
          TSS <- TPR + TNR - 1
          AUC <- ev0$AUC
          SORENSEN <- 2 * TPR/(FNR + (2 * TPR) + FPR)
          CBI <- ecospat.boyce(DDATA[which(DDATA$OBS ==
                                             0), ]$PRED, DDATA[which(DDATA$OBS == 1),
                                             ]$PRED, PEplot = FALSE)$cor

          maha_prob_train<-maha_prob[which(maha_prob$globalID %in%
                                             PRED_sites_train$globalID), ]$MD
          maha_prob_test<-obs
          n_vals <- length(maha_prob_test)

          suit_at_th <- quantile(maha_prob_train,0.1)
          omr <- 1 - length(which(maha_prob_test>=suit_at_th))/n_vals

          c(AUC = AUC, TSS = TSS, CBI = CBI, SORENSEN = SORENSEN,OMR=omr)
        }
      stopCluster(cl)
    }else {
      EVALUATIONS <- list()
      for (s in ss) {
        s.dat.train <- formatted_data$input_ones[-s,
                                                 , drop = FALSE]
        if (is.null(formatted_data$time_col)) {
          s.dat.train$globalID <- s.dat.train[, formatted_data$geoID_col]
        }else {
          s.dat.train$globalID <- paste(s.dat.train[,
                                                    formatted_data$geoID_col], s.dat.train[,
                                                                                           formatted_data$time_col], sep = "_")
        }
        s.dat.test <- formatted_data$input_ones[s, ,
                                                drop = FALSE]
        if (is.null(formatted_data$time_col)) {
          s.dat.test$globalID <- s.dat.test[, formatted_data$geoID_col]
        }else {
          s.dat.test$globalID <- paste(s.dat.test[, formatted_data$geoID_col],
                                       s.dat.test[, formatted_data$time_col],
                                       sep = "_")
        }
        mod1 <- enfa.custom(x = formatted_data$input_back[,
                                                          !grepl(formatted_data$obs_col, colnames(formatted_data$input_back))],
                            s.dat = s.dat.train[, grepl(paste(c(formatted_data$obs_col,
                                                                formatted_data$time_col, formatted_data$geoID_col),
                                                              collapse = "|"), colnames(s.dat.train))],
                            obs_col = formatted_data$obs_col, geoID_col = formatted_data$geoID_col,
                            time_col = formatted_data$time_col)
        if (sig_axes_selection == "brStick")
          sig <- ifelse(CENFA_brStick(mod1$sf) < 2, 2,
                        CENFA_brStick(mod1$sf))
        if (sig_axes_selection == "user.defined")
          sig <- user_defined_sig_axes
        if (is.null(formatted_data$time_col)) {
          PRED <- ENphylo_prediction(object = list(mod1),
                                     newdata = formatted_data$input_back[, !grepl(paste(c(formatted_data$geoID_col,
                                                                                          formatted_data$obs_col), collapse = "|"),
                                                                                  colnames(formatted_data$input_back))])$enfa_prediction
        }else {
          PRED <- ENphylo_prediction(object = list(mod1),
                                     newdata = formatted_data$input_back[, !grepl(paste(c(formatted_data$geoID_col,
                                                                                          formatted_data$obs_col, formatted_data$time_col),
                                                                                        collapse = "|"), colnames(formatted_data$input_back))])$enfa_prediction
        }
        PRED[, formatted_data$geoID_col] <- formatted_data$input_back[,
                                                                      formatted_data$geoID_col]
        if (!is.null(formatted_data$time_col)) {
          PRED[, formatted_data$time_col] <- formatted_data$input_back[,
                                                                       formatted_data$time_col]
          PRED$globalID <- paste(PRED[, formatted_data$geoID_col],
                                 PRED[, formatted_data$time_col], sep = "_")
        }else{
          PRED$globalID <- PRED[, formatted_data$geoID_col]
        }


        PRED_sites_train <- PRED[which(PRED$globalID %in%
                                         s.dat.train$globalID), ]
        PRED_sites_test <- PRED[which(PRED$globalID %in%
                                        s.dat.test$globalID), ]
        a <- PRED
        if (nrow(a) > 10000) {
          k <- sample(1:nrow(a), 10000)
          a <- a[k, ]
        }else k <- 1:nrow(a)
        a_scaled <- dudi.pca(PRED[, 1:sig], scannf = FALSE)$tab
        a_scaled$globalID <- PRED$globalID
        p_scaled <- a_scaled[which(a_scaled$globalID %in%
                                     PRED_sites_train$globalID), ]
        maha_prob <- mahasuhab.custom(a_scaled[, 1:sig,
                                               drop = FALSE], p_scaled[, 1:sig, drop = FALSE])
        maha_prob$globalID <- PRED$globalID
        fit <- maha_prob$MD[k]
        obs <- maha_prob[which(maha_prob$globalID %in%
                                 PRED_sites_test$globalID), ]$MD
        DDATA <- data.frame(ID = 1:(length(fit) + length(obs)),
                            OBS = c(rep(1, length(obs)), rep(0, length(fit))),
                            PRED = c(obs, fit))
        th <- optimal.thresholds(DDATA)[3, 2]
        ev0 <- presence.absence.accuracy(DDATA, th)
        TPR <- ev0$sensitivity
        TNR <- ev0$specificity
        FPR <- 1 - TPR
        FNR <- 1 - TNR
        TSS <- TPR + TNR - 1
        AUC <- ev0$AUC
        SORENSEN <- 2 * TPR/(FNR + (2 * TPR) + FPR)
        CBI <- ecospat.boyce(DDATA[which(DDATA$OBS ==
                                           0), ]$PRED, DDATA[which(DDATA$OBS == 1), ]$PRED,
                             PEplot = F)$cor

        maha_prob_train<-maha_prob[which(maha_prob$globalID %in%
                                           PRED_sites_train$globalID), ]$MD
        maha_prob_test<-obs
        n_vals <- length(maha_prob_test)

        suit_at_th <- quantile(maha_prob_train,0.1)
        omr <- 1 - length(which(maha_prob_test>=suit_at_th ))/n_vals

        EVALUATIONS[[which(sapply(ss, function(x) all(s ==
                                                        x)))]] <- c(AUC = AUC, TSS = TSS, CBI = CBI,
                                                                    SORENSEN = SORENSEN,OMR=omr)
      }
    }
    closeAllConnections()
    gc()
    EVALUATIONS <- t(EVALUATIONS)
    EVALUATIONS <- do.call(rbind, EVALUATIONS)
    rownames(EVALUATIONS) <- paste("replicate", 1:length(ss),
                                   sep = "_")
  } else EVALUATIONS <- NULL
  output <- list(call = "calibrated_enfa", full_model = FULL_MODEL,
                 evaluation = EVALUATIONS)
}
