ecospat.ESM.Modeling_CUSTOM<-function(data, 
                                      NbRunEval = NULL, 
                                      DataSplit, 
                                      DataSplitTable = NULL, 
                                      Yweights=NULL, 
                                      VarImport=0,
                                      weighting.score, 
                                      rescal.all.models=T,
                                      models, 
                                      modeling.id = as.character(format(Sys.time(), "%s")), 
                                      models.options, 
                                      which.biva = NULL, 
                                      parallel, 
                                      ncores=NULL, 
                                      cleanup = F) 
{
  require(doParallel)
  require(foreach)
  require(gtools)
  if (!weighting.score %in% c("AUC", "TSS", "Boyce", "Kappa", 
                              "SomersD")) {
    stop("weighting score not supported! Choose one of the following: AUC, TSS, Boyce, Kappa or SomersD")
  }
  if (data@has.data.eval) {
    stop("Evaluation with independant data is not supported yet!")
  }
  if ("PA" %in% slotNames(data)) {
    if (ncol(data@PA) > 1) {
      stop("It is not possible to use more than one Pseudo Absences dataset")
    }
  }
  #models.eval.meth <- weighting.score
  if (weighting.score == "AUC" | weighting.score == "SomersD") {
    models.eval.meth <- "ROC"
  }
  if (weighting.score == "Kappa") {
    models.eval.meth <- "KAPPA"
  }
  if (weighting.score == "TSS") {
    models.eval.meth <- "TSS"
  }
  if (weighting.score == "Boyce") {
    models.eval.meth <- "ROC"
  }
  models <- sort(models)
  iniwd <- getwd()
  dir.create(paste("./ESM.BIOMOD.output", data@sp.name, sep = "_"))
  newwd <- paste(getwd(), "/ESM.BIOMOD.output_", data@sp.name, 
                 sep = "")
  setwd(newwd)
  combinations <- combn(colnames(data@data.env.var), 2)
  if (is.null(which.biva)) {
    which.biva <- 1:ncol(combinations)
  }
  mydata <- data
  if (is.null(DataSplitTable)) {
    mod.prep.dat <- .Models.prepare.data(mydata, NbRunEval, 
                                         DataSplit, Yweights = NULL, Prevalence = 0.5, do.full.models = T)
    if (length(dim(mod.prep.dat[[1]]$calibLines)) == 3) {
      calib.lines <- mod.prep.dat[[1]]$calibLines[, , 1]
    }
    if (length(dim(mod.prep.dat[[1]]$calibLines)) == 2) {
      calib.lines <- mod.prep.dat[[1]]$calibLines
    }
    rm(mod.prep.dat)
  }
  else {
    calib.lines <- DataSplitTable
  }
  if (is.null(NbRunEval)) {
    if (ncol(calib.lines > 1)) {
      if (sum(!calib.lines[, ncol(calib.lines)]) == 0) {
        NbRunEval <- ncol(calib.lines) - 1
      }
      else {
        NbRunEval <- ncol(calib.lines)
      }
    }
    else {
      if (sum(!calib.lines) == 0) {
        NbRunEval <- ncol(calib.lines) - 1
      }
      else {
        NbRunEval <- ncol(calib.lines)
      }
    }
  }
  mymodels <- list()
  if (parallel == F) {
    for (k in which.biva) {
      mydata@data.env.var <- data@data.env.var[, colnames(data@data.env.var) %in% 
                                                 combinations[, k]]
      mydata@sp.name <- paste("ESM.BIOMOD", k, sep = ".")
      mymodels[[k]] <- BIOMOD_Modeling(data = mydata, models = models, 
                                       models.options = models.options, models.eval.meth = models.eval.meth, 
                                       DataSplitTable = calib.lines, Yweights=Yweights, 
                                       rescal.all.models = rescal.all.models, do.full.models = T, VarImport = VarImport, 
                                       modeling.id = modeling.id)
      if (cleanup != F) {
        removeTmpFiles(h = cleanup)
      }
    }
  }
  if (parallel == T) {
    cl <- makeCluster(ncores, outfile="log.txt")
    registerDoParallel(cl)
    
    mymodels <- foreach(k = which.biva, .packages = c("biomod2", "raster")) %dopar% 
    {
      setwd(newwd)
      mydata@data.env.var <- data@data.env.var[, colnames(data@data.env.var) %in% 
                                                 combinations[, k]]
      mydata@sp.name <- paste("ESM.BIOMOD", k, sep = ".")
      if (cleanup != F) {
        removeTmpFiles(h = cleanup)
      }
      BIOMOD_Modeling(data = mydata, models = models, 
                      models.options = models.options, models.eval.meth = models.eval.meth, 
                      DataSplitTable = calib.lines, Yweights=Yweights, 
                      rescal.all.models = rescal.all.models, do.full.models = T, 
                      VarImport = VarImport, modeling.id = modeling.id)
    }
    stopCluster(cl)
  }
  output <- list(modeling.id = modeling.id, 
                 models. = grep(modeling.id, mixedsort(list.files(getwd(), "models.out", 
                                                                  recursive = T, 
                                                                  full.names = T)), 
                                value = T), 
                 models = models, 
                 calib.lines = calib.lines, 
                 NbRunEval = NbRunEval, 
                 data = data, 
                 wd = getwd(), 
                 which.biva = which.biva, 
                 mymodels = mymodels)
  save(output, file = paste("ESM_Modeling..models", modeling.id, 
                            "out", sep = "."))
  setwd(iniwd)
  return(output)
}