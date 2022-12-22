ecospat.ESM.Projection_CUSTOM<-function (ESM.modeling.output, 
                                         new.env, 
                                         proj.name, 
                                         parallel = F, 
                                         ncores=NULL){
  require(gtools)
  require(doParallel)
  iniwd <- getwd()
  setwd(ESM.modeling.output$wd)
  models <- ESM.modeling.output$models
  models. <- ESM.modeling.output$models.[mixedorder(as.numeric(gsub("/", "", sapply(strsplit(ESM.modeling.output$models., "ESM.BIOMOD."), "[[", 3))))]
  mymodels <- ESM.modeling.output$mymodels
  # FAILED=sapply(lapply(mymodels, function(x)x@models.failed), function(x)all(grepl(paste("RUN", ESM.modeling.output$NbRunEval+1, sep=""), x)))
  FAILED=sapply(lapply(mymodels, function(x)x@models.failed), function(x){
    length(grep(paste("RUN", ESM.modeling.output$NbRunEval+1, sep=""), x))==length(models)
  })
  mymodels<-mymodels[!FAILED]
  combinations <- combn(colnames(ESM.modeling.output$data@data.env.var),
                        2)[,!FAILED]

  which.biva <- ESM.modeling.output$which.biva
  NbRunEval <- ESM.modeling.output$NbRunEval
  modeling.id <- ESM.modeling.output$modeling.id
  name.env <- deparse(substitute(new.env))
  if (parallel == F) {
    for (k in 1:length(mymodels)) {
      mymodel <- mymodels[[k]]
      if (is.data.frame(new.env)) {
        BIOMOD_Projection(modeling.output = mymodel, 
                          new.env = new.env[, colnames(new.env) %in% 
                                              combinations[, k]], 
                          proj.name = paste(proj.name, "ESM.BIOMOD", k, mymodel@modeling.id, sep = "."), 
                          selected.models = c(grep("Full", mymodel@models.computed, 
                                                   value = T), grep(paste("RUN", NbRunEval + 
                                                                            1, sep = ""), mymodel@models.computed, value = T)), 
                          do.stack = F, build.clamping.mask = F)
      }
      if (class(new.env) == "RasterStack") {
        BIOMOD_Projection(modeling.output = mymodel, 
                          new.env = new.env[[which(names(new.env)%in%combinations[, k])]], 
                          proj.name = paste(proj.name, 
                                            "ESM.BIOMOD", 
                                            k, 
                                            mymodel@modeling.id, 
                                            sep = "."), 
                          selected.models = c(grep("Full", mymodel@models.computed, value = T), 
                                              grep(paste("RUN", NbRunEval + 1, sep = ""), 
                                                   mymodel@models.computed, value = T)), 
                          do.stack = T, 
                          build.clamping.mask = F)
      }
    }
  }
  if (parallel == T) {
    
    cl <- makeCluster(ncores)
    registerDoParallel(cl)
    
    foreach(k = 1:length(mymodels), .packages = c("biomod2", "raster")) %dopar% {
      mymodel <- mymodels[[k]]
      if (is.data.frame(new.env)) {
        BIOMOD_Projection(modeling.output = mymodel, 
                          new.env = new.env[, colnames(new.env) %in% 
                                              combinations[, k]], proj.name = paste(proj.name, 
                                                                                    strsplit(mymodel@models.computed[1], "_")[[1]][1], 
                                                                                    mymodel@modeling.id, sep = "."), 
                          selected.models = c(grep("Full", mymodel@models.computed, 
                                                   value = T), grep(paste("RUN", NbRunEval + 
                                                                            1, sep = ""), mymodel@models.computed, 
                                                                    value = T)), do.stack = F, build.clamping.mask = F)
      }
      if (class(new.env) == "RasterStack") {
        BIOMOD_Projection(modeling.output = mymodel, 
                          new.env = new.env[[which(names(new.env) %in% 
                                                     combinations[, k])]], proj.name = paste(proj.name, 
                                                                                             strsplit(mymodel@models.computed[1], "_")[[1]][1], 
                                                                                             mymodel@modeling.id, sep = "."), 
                          selected.models = c(grep("Full", mymodel@models.computed, 
                                                   value = T), grep(paste("RUN", NbRunEval + 
                                                                            1, sep = ""), mymodel@models.computed, 
                                                                    value = T)), do.stack = T, build.clamping.mask = F)
      }
    }
    stopCluster(cl)
  }
  removeTmpFiles(h = T)
  pred.biva <- grep(modeling.id, mixedsort(list.files(getwd(), 
                                                      paste("proj_", proj.name, sep = ""), recursive = T, full.names = T)), 
                    value = T)
  output <- list(modeling.id = modeling.id, models. = grep(modeling.id, 
                                                           mixedsort(list.files(getwd(), "models.out", recursive = T, 
                                                                                full.names = T)), value = T), models = models, pred.biva = grep(modeling.id, 
                                                                                                                                                mixedsort(list.files(getwd(), paste("proj_", proj.name, 
                                                                                                                                                                                    sep = ""), recursive = T, full.names = T)), value = T), 
                 NbRunEval = NbRunEval, name.env = name.env, new.env.raster = class(new.env) == 
                   "RasterStack", wd = getwd(), which.biva = which.biva)
  save(output, file = paste("ESM_Projections", modeling.id, 
                            "out", sep = "."))
  setwd(iniwd)
  return(output)
}