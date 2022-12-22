ESM.prediction.output=myproj
ESM.EnsembleModeling.output=myEM
models <- ESM.prediction.output$models
weights <- ESM.EnsembleModeling.output$weights
weights.EF <- ESM.EnsembleModeling.output$weights.EF
NbRunEval <- ESM.prediction.output$NbRunEval
pred.biva <- ESM.prediction.output$pred.biva
new.env.raster <- ESM.prediction.output$new.env.raster
failed.mod <- grep(paste("RUN", NbRunEval + 1, sep = ""), 
                   unlist(ESM.EnsembleModeling.output$failed), value = T)
if (!exists("new.env")) {
  stop("new.env object required!")
}
if (new.env.raster) 
  #pred.biva <- grep("gri", pred.biva, value = T)
  pred.biva <- pred.biva[-grep(".grd", pred.biva, value = F)]
if (!new.env.raster) 
  pred.biva <- grep("RData", pred.biva, value = T)
biva.proj <- list()
for (i in 1:length(pred.biva)){
  if (!new.env.raster)   {
    biva.proj[[i]] <- as.data.frame(get(load(pred.biva[i])))
    colnames(biva.proj[[i]]) <- gsub("AllData", "BIOMOD", colnames(biva.proj[[i]]))
    colnames(biva.proj[[i]]) <- gsub(paste("RUN", NbRunEval + 1, sep = ""), "ESM", colnames(biva.proj[[i]]))
    colnames(biva.proj[[i]]) <- paste(colnames(biva.proj[[i]]), i, sep = ".")
  }
  if (new.env.raster) {biva.proj[[i]] <- stack(pred.biva[i])}
}

if (!new.env.raster) {
  pred.ESM <- list()
  biva.st <- do.call(cbind, biva.proj)
  for (i in 1:length(models)) {
    if (length(models) > 1) {
      wm <- weights[grep(models[i], names(weights))][names(weights[grep(models[i], 
                                                                        names(weights))]) %in% colnames(biva.st[, grep(models[i], 
                                                                                                                       colnames(biva.st))])]
      pred.ESM[[i]] <- apply(biva.st[, grep(models[i], 
                                            colnames(biva.st))], 1, function(x) weighted.mean(x, 
                                                                                              wm, na.rm = T))
    }
    else {
      pred.ESM[[i]] <- apply(biva.st[, grep(models[i], 
                                            colnames(biva.st))], 1, function(x) weighted.mean(x, 
                                                                                              wm, na.rm = T))
    }
  }
  rm(biva.st)
  pred.ESM <- as.data.frame(do.call(cbind, pred.ESM))
  colnames(pred.ESM) <- models
}
if (new.env.raster) {
  pred.ESM <- stack(biva.proj)
  for (i in 1:length(models)){
    if(nlayers(pred.ESM[[grep(models[i], names(pred.ESM))]])!=0)assign(models[i], pred.ESM[[grep(models[i], names(pred.ESM))]])
    if(nlayers(pred.ESM[[grep(models[i], names(pred.ESM))]])==0){
      empty<-biva.proj[[1]][[1]]
      empty[]<-NA
      names(empty)<-"layer"
      assign(models[i], empty)
      next
    }
    
    n <- grep(models[i], failed.mod, value = T)
    if (length(models) > 1){
      weights.mod <- weights[grep(models[i], names(weights))]
      
      weights.mod <- weights.mod[gsub(paste(models[i], ".", sep=""), 
                       "", 
                       names(weights.mod))%in%sapply(strsplit(names(get(models[i])), 
                                                              "_"), "[[", 1)]
      
      
      if (length(n) > 0){
        weights.mod <- weights.mod[!names(weights.mod)%in%paste(models[i], 
                                                                ".", 
                                                                unlist(strsplit(n, "_"))[((1:length(n)) * 4) - 3], 
                                                                sep = "")]
        assign(models[i], round(weighted.mean(get(models[i]), weights.mod, na.rm = T)))
      } else {
        assign(models[i], round(weighted.mean(get(models[i]), weights.mod, na.rm = T)))
      }
    }else{
      if (length(n) > 0) {
        #weights.mod <- weights[grep(paste(unlist(strsplit(n, "_"))[1], "_", sep = ""), paste(names(weights), "_", sep = ""), invert = T)]
        weights.mod <- weights[-match(paste(".", 
                                            sapply(strsplit(n, "_"), "[[", c(1)), 
                                            "_", 
                                            sep = ""), 
                                      paste(names(weights), 
                                            "_", 
                                            sep = ""))]
        if(all(is.na(weights.mod)))weights.mod <- weights
        
        X<-get(models[i])
        Y<-names(X)
        Y<-sapply(strsplit(Y, "_"), "[[", c(1))
        Y<-paste(".", Y,  sep="")
        Y<-sapply(names(weights.mod), function(x)grep(paste(x, "$", sep=""), Y))
        
        assign(models[i], round(weighted.mean(X[[Y]], 
                                              weights.mod, na.rm = T)))
      }
      else {
        X<-get(models[i])
        Y<-names(X)
        Y<-sapply(strsplit(Y, "_"), "[[", c(1))
        Y<-paste(".", Y,  sep="")
        Y<-sapply(names(weights), function(x)grep(paste(x, "$", sep=""), Y))
        assign(models[i], round(weighted.mean(X[[Y]], weights, na.rm = T)))
      }
    }
  }
  pred.ESM <- stack(mget(models))[[order(models)]]
  do.call("rm", as.list(models))
}
if (length(models) > 1) {
  if (new.env.raster) {
    ESM.EF <- round(weighted.mean(pred.ESM[[names(pred.ESM)]], 
                                  weights.EF[order(weights.EF[, 1]), 2], na.rm = T))
    pred.ESM <- stack(pred.ESM, ESM.EF)
    names(pred.ESM) <- c(names(pred.ESM)[1:(nlayers(pred.ESM) - 
                                              1)], "EF")
    rm(ESM.EF)
  }
  if (!new.env.raster) {
    pred.ESM$EF <- apply(pred.ESM, 1, function(x) weighted.mean(x, 
                                                                weights.EF[order(weights.EF[, 1]), 2], na.rm = T))
  }
  if (!new.env.raster) {
    pred.ESM <- round(pred.ESM)
  }
}
myEF<-pred.ESM
