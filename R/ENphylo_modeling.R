#'@title Calculating species marginality and specialization via ENFA and
#'  phylogenetic imputation
#'@description The function computes vectors of marginality and specialization
#'  according to \cite{Rinnan & Lawler (2019)} via Environmental Niche Factor
#'  Analysis (ENFA) and phylogenetic imputation (\cite{Garland & Ives, 2000}).
#'  It takes a list of \code{Spatial} or \code{Simple Features} (\pkg{sp} or
#'  \pkg{sf} objects) and a phylogenenetic tree to train ENFA and/or ENphylo
#'  models. This function allows to calibrate and evaluate both model techniques
#'  for a given species while accounting for phylogenetic uncertainty.
#'  Calibrations are made on a random subset of the data under the bootstrap
#'  cross-validation scheme. The predictive power of the different models is
#'  estimated using four different evaluation metrics.
#'@usage ENphylo_modeling(input_data, tree, input_mask, obs_col, time_col=NULL,
#'  min_occ_enfa=50, boot_test_perc=20, boot_reps=10, nsim=10, si=0.2, si2=0.2,
#'  eval_metric_for_imputation=c("AUC","TSS","CBI","SORENSEN","OMR"),
#'  output_options=c("full","weighted.mean","best"), eval_threshold=0.7,
#'  clust=0.5, external_enfa_models=NULL, spec_for_imputation=NULL)
#'@param input_data a list of \code{SpatialPointsDataFrame} or
#'  \code{sf::data.frame} objects containing species occurrence data in binary
#'  format (ones for presence, zero for background points) along with the
#'  explanatory variables to be used in ENFA or ENphylo. Each element of the
#'  list must be named (i.e. the species name). If \code{external_enfa_models}
#'  are provided, \code{input_data} should not be supplied.
#'@param tree an object of class \code{phylo} containing a phylogenetic tree
#'  written in Newick format with all the species considered in the
#'  \code{INPUT_DATA} argument. The tree needs not to be ultrametric or fully
#'  dichotomous.
#'@param input_mask the geographical mask defining the spatial domain
#'  encompassing the background area enclosing all the species in the tree.
#'@param obs_col character. Name of the column containing the vector of species
#'  occurrence data in binary format.
#'@param time_col character. Name of column containing the time intervals
#'  associated to each species presence and background point (optional).
#'@param min_occ_enfa numeric. The minimum number of occurrence data required
#'  for a species to be modelled with ENFA.
#'@param boot_test_perc numeric. Percentage of data (comprised between 0 and
#'  100) used to calibrate ENFA and/or ENphylo models within a bootstrap
#'  cross-validation scheme (the \code{100-boot_test_perc} remaining percentage
#'  will be used for evaluating model performances).
#'@param boot_reps numeric. Number of evaluation runs performed within the
#'  bootstrap cross-validation scheme to evaluate ENFA and/or ENphylo models. If
#'  set to 0, models evaluation is skipped.
#'@param nsim numeric. Number of alternative phylogenies with altered topology
#'  and branch lengths (as compared to the reference tree) to be tested by means
#'  of function \code{\link[RRphylo]{swapONE}}
#'@param si numeric. As in \code{\link[RRphylo]{swapONE}}, the proportion of
#'  tips whose topologic arrangement will be swapped to provide alternative
#'  trees.
#'@param si2 numeric. As in \code{\link[RRphylo]{swapONE}}, the proportion of
#'  nodes whose age will be changed to provide alternative trees.
#'@param eval_metric_for_imputation character. Evaluation metric used to select
#'  the most accurate ENphylo models. The viable options are: \code{"AUC"},
#'  \code{"TSS"}, \code{"CBI"}, \code{"SORENSEN"}, or \code{"OMR"}.
#'@param output_options character. The strategy adopted to return ENphylo models
#'  results (see details). The viable options are: \code{"full"},
#'  \code{"weighted.mean"}, and \code{"best"}.
#'@param eval_threshold numeric. The minimum evaluation score below which
#'  ENphylo models derived from the swapped trees are excluded from the function
#'  output (only used if \code{output_options = "weighted.mean"} or
#'  \code{"best"}).
#'@param clust numeric. The proportion of cores used to train ENFA and ENphylo
#'  models. If \code{NULL}, parallel computing is disabled. It is set at 0.5 by
#'  default.
#'@param external_enfa_models list. If not \code{NULL}, an object returned by
#'  the \code{ENphylo_modeling} function where all the species are modelled with
#'  ENFA.
#'@param spec_for_imputation character. The names of the species among those
#'  provided in the \code{external_enfa_models} argument, whose ENFA models will
#'  be replaced with a new ENphylo model (used only if
#'  \code{external_enfa_models} argument is not \code{NULL}).
#'@author Alessandro Mondanaro, Mirko Di Febbraro, Silvia Castiglione, Carmela
#'  Serio, Marina Melchionna, Pasquale Raia
#'@details \code{ENphylo_modeling} automatically arranges the input data in a
#'  convenient format to run ENFA or ENphylo. The internal call of the function
#'  is \code{"calibrated_enfa"} for ENFA and \code{"calibrated_imputed"} for
#'  ENphylo, respectively. The function requires a list as the input data. Each
#'  element of the list must be provided with a name (e.g. species names). Model
#'  evaluation can be switched off by setting \code{boot_rep} =0. In this case,
#'  the internal evaluation element will return \code{NULL}. The function cannot
#'  work with \code{nsim} < 1 since one of the strongest points of
#'  \code{ENphylo_modeling} is to test alternative phylogenies for providing the
#'  most accurate possible reconstruction of the species environmental
#'  preferences.
#'
#'  ENphylo automatically switches from ENFA to ENphylo algorithm for any
#'  species having ENFA model accuracy below the minimum evaluation score set by
#'  the user (i.e. \code{eval_threshold}). In this case, both models are
#'  performed but only the one performing best according to validation metrics
#'  is retained. Overall, the function is able to phylogenetically impute up to
#'  30 percent of the species on the tree. If the sum of species with a low
#'  number of occurrences (i.e. a number below \code{min_occ_enfa}) or poorly
#'  modeled with ENFA (i.e. evaluation score below \code{eval_threshold})
#'  exceeds 30\%, ENphylo is set to use lower \code{min_occ_enfa} or lower
#'  \code{eval_threshold} progressively adding poorly modeled species until the
#'  remaining set amounts to less than 30\% of the total. In this case, the
#'  function reports a warning where the new \code{min_occ_enfa} and/or
#'  \code{eval_threshold} value is expressed. If \code{ENphylo_modeling} runs
#'  the ENphylo algorithm, the outputs depends on the strategy adopted by the
#'  user through the \code{output_options} argument. If
#'  \code{output_options="full"}, all CO matrices and evaluation metrics for all
#'  the swapped trees are reported. With \code{output_options="weighted.mean"},
#'  the output takes a subset of CO matrices and evaluation metrics only for
#'  those tree swapping iterations, achieving a predictive accuracy above the
#'  user-defined threshold, as indicated in the
#'  \code{eval_metric_for_imputation} and \code{eval_threshold} arguments.
#'  Finally, if \code{output_options="best"}, a single CO matrix and evaluation
#'  scores vector corresponding to the most accurate swapped tree is reported.
#'  If any of the tree swapping iterations under either \code{"best"} or
#'  \code{"weighted.mean"}  provides above-thethreshold accuracy, the function
#'  automatically switches to \code{"full"} strategy. In any case, the function
#'  returns a list element where model options and metric used to validate model
#'  are expressed. Eventually, the function creates two new folders,
#'  "ENphylo_enfa_models" and "ENphylo_imputed_models", in the current working
#'  directory. In each of these folders, a number of new named folders equal to
#'  the number of modelled species are created. Model outputs will be saved
#'  therein as "model_output.RData".
#'@return A list of length equal to the number of modelled species. For each
#'  element of this list, the output contains three elements: \enumerate{ \item
#'  \strong{$call} a character specifying the algorithm used to model the
#'  species (i.e. ENFA or ENphylo). \item \strong{$formatted data} a list of
#'  input data formatted to run either ENFA or ENphylo algorithms. Specifically,
#'  the list reports \code{$input_ones} the presence data points,
#'  \code{$input_back} the background points, \code{$one_COORDS} the coordinates
#'  of presence data only, \code{$study_area} the background area, and the name
#'  of the columns associated to the arguments \code{OBS_col} and
#'  \code{TIME_factor_col}. \item \strong{$calibrated_model} a list of three
#'  elements. The output objects are different depending on whether ENFA or
#'  ENphylo is used to model the species:}
#'@return  ENFA
#'@return \itemize{\item$call: a character specifying the algorithm used.
#'  \item$full_ model: a list containing marginality and specialization factors,
#'  the 'co' matrix, the number of significant axes, and all the other objects
#'  generated by applying ENFA on the entire occurrence dataset (see
#'  \cite{Rinnan et al. 2019} for additional details). \item$evaluation: a
#'  matrix containing the evaluation scores of the ENFA model assessed by all
#'  possible evaluation metrics (i.e. Area Under the Curve (AUC), True Skill
#'  Statistic (TSS), Boyce Index, Sorensen Index, and Omission Rate (OMR)) for
#'  each model evaluations run.}
#'@return ENphylo
#'@return \itemize{\item$call: a character specifying the algorithm used.
#'  \item$co: a list of the 'co' matrices of length equal to the number of
#'  alternative phylogenies tested (i.e. \code{nsim} argument).
#'  \item$evaluation: a data.frame object containing the evaluation scores of
#'  ENphylo model assessed by the evaluation metrics for each alternative
#'  phylogeny. The output of this object depends on the strategy adopted by the
#'  user through the \code{output_options} argument (see
#'  details).\item$output_options:  a character vector where the model options
#'  to run ENphylo and the evaluation metric used to assess the accuracy of
#'  ENphylo model are expressed.\item$evaluation: a matrix containing the
#'  evaluation scores of the ENFA model assessed by all possible evaluation
#'  metrics (i.e. Area Under the Curve (AUC), True Skill Statistic (TSS), Boyce
#'  Index, Sorensen Index, and Omission Rate (OMR)) for each model evaluations
#'  run.}
#'@importFrom raster rasterToPoints crop extent extend raster res
#'@importFrom ape drop.tip vcv
#'@importFrom methods extends as
#'@importFrom parallel detectCores
#'@importFrom sf as_Spatial
#'@importClassesFrom sp SpatialPixelsDataFrame
#'@export
#'@references Rinnan, D. S., &  Lawler, J. (2019). Climate-niche factor
#'  analysis: a spatial approach to quantifying species vulnerability to climate
#'  change. \emph{Ecography}, 42(9), 1494–1503. doi/full/10.1111/ecog.03937
#'@references Garland, T., & Ives, A. R. (2000). Using the past to predict the
#'  present: Confidence intervals for regression equations in phylogenetic
#'  comparative methods. \emph{American Naturalist}, 155(3),346–364.
#'  doi.org/10.1086/303327
#'@examples
#' \dontrun{
#'setwd("YOUR_DIRECTORY")
#'getwd()->main.dir
#'url<-"https://www.dropbox.com/s/wm7qcuwm0nm35w2/ENphylo%20code%26data%202.zip?dl=1"
#'download.file(url,file.path(main.dir,"ENphylo code&data 2.zip"),mode="wb")
#'unzip("ENphylo code&data 2.zip")
#'load("ENphylo code&data/example_data.RData")
#'MASK_FULL<-raster::raster("ENphylo code&data/variable_bio1.tif")
#'external_data<-raster::stack(list.files("ENphylo code&data/external_data",full.names=TRUE))
#'
#'## NOTE: Given the size of the data, running the function is time-comsuming
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
#' ### CASE 4
#' # Run ENphylo under alternative model evaluation metric and strategy
#' ENmod4<-ENphylo_modeling(external_enfa_models=ENmod,
#'                              spec_for_imputation=names(ENmod)[k],
#'                              tree=tree,
#'                              output_options="best",
#'                              eval_metric_for_imputation="TSS",
#'                              eval_threshold = 0.4)
#' gc()
#'}


ENphylo_modeling<-function (input_data = NULL,
                            tree,
                            input_mask,
                            obs_col, time_col = NULL,
                            min_occ_enfa = 50,
                            boot_test_perc = 20,
                            boot_reps = 10,
                            nsim = 10,
                            si = 0.2,
                            si2 = 0.2,
                            eval_metric_for_imputation = c("AUC","TSS", "CBI", "SORENSEN", "OMR"),
                            output_options = c("full", "weighted.mean", "best"),
                            eval_threshold = 0.7,
                            clust = 0.5,
                            external_enfa_models = NULL,
                            spec_for_imputation = NULL)
{
  if (is.null(external_enfa_models) && any(is.null(names(input_data))))
    stop("all the elements in input_data list must be provided with a species name")
  if (is.null(external_enfa_models) && any(!names(input_data) %in%
                                           tree$tip.label))
    stop("all species in input_data must be on the phylogenetic tree")
  if (boot_reps > 0)
    evaluate_imputed = TRUE else evaluate_imputed = FALSE
    if (evaluate_imputed & length(eval_metric_for_imputation) >
        1)
      stop("Please, set just one evaluation metric for imputation")
    if (evaluate_imputed & length(output_options) > 1)
      stop("Please, set just one output option")
    if (any(tree$edge.length == 0))
      tree$edge.length[which(tree$edge.length < max(diag(vcv(tree)))/1000)] <- tree$edge.length[which(tree$edge.length <
                                                                                                        max(diag(vcv(tree)))/1000)] + max(diag(vcv(tree))) *
        0.001
    if (!is.null(clust))
      clust <- detectCores() * clust
    if (is.null(external_enfa_models)) {
      if (sum(sapply(input_data, function(x) sum(as.data.frame(x)[,
                                                                  obs_col])) < min_occ_enfa) >= length(input_data) *
          0.3) {
        xx <- sort(sapply(input_data, function(x) sum(as.data.frame(x)[,
                                                                       obs_col])))[ceiling(length(input_data) * 0.1)]
        min_occ_enfa <- xx + 1
        warning(paste("The number of species with occurrences below min_occ_enfa exceeds the 30% of total species.\n                      min_occ_enfa value was reduced to",
                      xx, sep = " "))
      }
      all_models <- mapply(function(data, nam) {
        cat(paste("\n", "Modelling", nam, "\n"))
        if (is.null(time_col)) {
          prova <- subset(data, as.data.frame(data)[, obs_col] ==
                            0)
        } else {
          prova <- subset(data, as.data.frame(data)[, obs_col] ==
                            0 & as.data.frame(data)[, time_col] == unique(as.data.frame(data)[,
                                                                                              time_col])[1])
        }
        qq <- input_mask
        vv <- prova[, obs_col]
        if (!extends(class(vv), "Spatial"))
          vv <- as_Spatial(vv)
        qq <- as(qq, "SpatialPixelsDataFrame")
        qq <- qq[vv, ]
        qq <- raster(qq)
        qq_RTP <- apply(rasterToPoints(qq)[, 1:2], 2, range)
        qq_RTP <- extent(qq_RTP[, 1], qq_RTP[, 2])
        qq_RTP <- extend(qq_RTP, res(qq)[1])
        maskk <- crop(qq, qq_RTP)
        if (!extends(class(data), "Spatial"))
          data <- as_Spatial(data)
        mydata <- DATA_PREPARATION(species_input_data = data,
                                   obs_col = obs_col, input_mask = maskk, time_col = time_col)
        if (length(mydata$ones_coords) >= min_occ_enfa) {
          mymodel <- ENFA_CALIBRATION(formatted_data = mydata,
                                      boot_test_perc = boot_test_perc, boot_reps = boot_reps,
                                      sig_axes_selection = "brStick", clust = clust)
        } else mymodel = NULL

        if(!is.null(mymodel)){
          model_outputs <-list(call = "calibrated_enfa", formatted_data = mydata,
                               calibrated_model = mymodel)
          list(model_outputs)->model_outputs
          names(model_outputs)<-nam
          dir.create(file.path("ENphylo_enfa_models", nam),
                     recursive = TRUE)
          save(model_outputs, file = paste0("ENphylo_enfa_models/",
                                            nam, "/output_model.RData"))
        }

        return(list(call = "calibrated_enfa", formatted_data = mydata,
                    calibrated_model = mymodel))

      }, data = input_data, nam = names(input_data), SIMPLIFY = FALSE)
      names(all_models) <- names(input_data)
      enf_mod <- all_models
      if (evaluate_imputed) {
        check <- sapply(all_models, function(k) {
          if (!is.null(k$calibrated_model)) {
            mean(k$calibrated_model$evaluation[, eval_metric_for_imputation]) <
              eval_threshold
          } else TRUE
        })
        vv <- eval_threshold
        if (sum(check) > (length(input_data) * 30/100)) {
          repeat {
            vv <- vv * 0.95
            check <- sapply(all_models, function(k) {
              if (!is.null(k$calibrated_model)) {
                mean(k$calibrated_model$evaluation[, eval_metric_for_imputation]) <
                  vv
              } else TRUE
            })
            if (sum(check) <= (length(input_data) * 30/100))
              break
          }
          all_models <- lapply(all_models, function(k) {
            if (!is.null(k$calibrated_model)) {
              if (mean(k$calibrated_model$evaluation[,
                                                     eval_metric_for_imputation]) < vv) {
                k <- list(call = k$call, formatted_data = k$formatted_data,
                          calibrated_model = NULL)
              }
            }
            k
          })
          warning(paste("eval_threshold value for ENFA models evaluation was reduced to",
                        vv, "because more than 30% of species have to be modelled with phylogenetic imputation",
                        sep = " "))
        } else {
          all_models <- lapply(all_models, function(k) {
            if (!is.null(k$calibrated_model)) {
              if (mean(k$calibrated_model$evaluation[,
                                                     eval_metric_for_imputation]) < vv) {
                k <- list(call = k$call, formatted_data = k$formatted_data,
                          calibrated_model = NULL)
              }
            }
            k
          })
        }
      }
      if (all(sapply(all_models, function(k) !is.null(k$calibrated_model)) ==
              TRUE)) {
        warning("All species are modelled with Enfa")
      } else {
        if (any(tree$tip.label %in% names(all_models) ==
                FALSE)) {
          sp <- tree$tip.label[!tree$tip.label %in% names(all_models)]
          tree <- drop.tip(tree, sp)
          warning(paste(as.data.frame(sp), " not present in input_data. They will be dropped from the tree",
                        sep = ""))
        }
        if (nsim == 0)
          stop("Please, set nsim as an integer >0")
        myimputed <- IMPUTED_CALIBRATION(ENFA_output = all_models,
                                         tree = tree, nsim = nsim, si = si, si2 = si2,
                                         clust = clust, evaluate = evaluate_imputed, boot_test_perc = boot_test_perc,
                                         boot_reps = boot_reps, output_options = output_options,
                                         eval.metric = eval_metric_for_imputation, eval_threshold = eval_threshold)
        sel <- lapply(names(myimputed), function(x) {
          if (is.null(enf_mod[x][[1]]$calibrated_model)) {
            a = 0
          } else {
            a <- mean(enf_mod[x][[1]]$calibrated_model$evaluation[,
                                                                  eval_metric_for_imputation])
          }
          b <- mean(myimputed[x][[1]]$evaluation[, eval_metric_for_imputation])
          b > a
        })
        myimputed <- myimputed[which(sel == TRUE)]
        myimputed <- mapply(function(x, y) {
          list(call = "calibrated_imputed", formatted_data = y$formatted_data,
               calibrated_model = x)
        }, x = myimputed, y = all_models[names(myimputed)],
        SIMPLIFY = FALSE)
        w <- match(names(myimputed), names(enf_mod))
        enf_mod[w] <- myimputed
        all_models <- enf_mod
      }
    } else {
      all_models <- external_enfa_models
      if (any(tree$tip.label %in% names(all_models) == FALSE)) {
        sp <- tree$tip.label[!tree$tip.label %in% names(all_models)]
        tree <- drop.tip(tree, sp)
      }
      if (is.null(spec_for_imputation))
        stop("Please, specify which species to impute")
      if (nsim == 0)
        stop("Please, set nsim as an integer >0")
      myimputed <- IMPUTED_CALIBRATION(ENFA_output = all_models,
                                       tree = tree, nsim = nsim, si = si, si2 = si2, clust = clust,
                                       spec_for_imputation = spec_for_imputation, evaluate = evaluate_imputed,
                                       boot_test_perc = boot_test_perc, boot_reps = boot_reps,
                                       output_options = output_options, eval.metric = eval_metric_for_imputation,
                                       eval_threshold = eval_threshold)
      myimputed <- mapply(function(x, y) {
        list(call = "calibrated_imputed", formatted_data = y$formatted_data,
             calibrated_model = x)
      }, x = myimputed, y = all_models[names(myimputed)], SIMPLIFY = FALSE)
      w <- match(names(myimputed), names(all_models))
      all_models[w] <- myimputed
    }
    return(all_models)
}

