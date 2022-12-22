### Installing/Updating packages
required<-c("ape", "raster", "methods","CENFA","pbapply","Rphylopars","RRphylo","dismo","gtools",
            "adehabitatMA", "ecospat", "foreach", "doParallel", "PresenceAbsence", 
            "parallel", "ade4", "sp", "biomod2", "terra", "sf", "adehabitatHS",
            "caret","ENMeval","rJava")
if(any(!required%in%installed.packages()[,1]))
  install.packages(required[which(!required%in%installed.packages()[,1])])

latest_version<-c("ecospat","dismo","biomod2","ENMeval","rJava")
latest_version[which(latest_version%in%rownames(old.packages()))]->to_update
if(length(to_update)>0) install.packages(to_update)

sapply(required,require,character.only=TRUE)

### Loading Data and RRdtn package
setwd("YOUR_DIRECTORY")

getwd()->main.dir
url<-"https://www.dropbox.com/sh/fwaie31o9zhg0x7/AAC1AA7nhqKzLgQrc_7Nok8ea/ENphylo%20code%26data.zip?dl=1"

download.file(url,file.path(main.dir,"ENphylo code&data.zip"),mode="wb")
unzip("ENphylo code&data.zip")
load("example_data.RData")
source("ecospat.ESM.EnsembleModeling_CUSTOM.R")
source("ecospat.ESM.Projection_CUSTOM.R")
source("ecospat.ESM.Modeling_CUSTOM.R")
MASK_FULL<-raster::raster("variable_bio1.tif")
install.packages("RRdtn_0.2.0.tar.gz", repos = NULL, type = "source", INSTALL_opts = "--install-tests")
require(RRdtn)


  
###DO NOT RUN, WE HAVE ALREADY UPLOADED "enfa_all_models" OBJECT INTO example_data.RData 
###EMBEDDED IN ENphylo code&data.zip TO REDUCE THE COMPUTATIONAL TIME.

# ENphylo_modeling(input_data=DATA_FULL,
#                  tree=tree,
#                  input_mask=MASK_FULL,
#                  obs_col="OBS",
#                  time_col="TIME_factor",
#                  min_occ_enfa = 10, ## this argument is not arbitrarily set because no species embedded in Mondanaro et al.2021 dataset has a number of occurrences <=10
#                  eval_metric_for_imputation ="AUC",
#                  output_options ="best",
#                  eval_threshold=0.5)->enfa_all_species


### running analyses ###
npoints=10
dir.create(paste("test_", npoints, "_points", sep=""))

eval_function<-function(test_p, test_a, mm, train_p, model_name){
  if(model_name!="ens.biv"){
    ev0<- dismo::evaluate(test_p, test_a, mm) 
    
    TPR<-ev0@TPR[which.max(ev0@TPR+ev0@TNR)]
    TNR<-ev0@TNR[which.max(ev0@TPR+ev0@TNR)]
    FPR<-ev0@FPR[which.max(ev0@TPR+ev0@TNR)]
    FNR<-ev0@FNR[which.max(ev0@TPR+ev0@TNR)]
    
    TSS<-TPR+TNR-1
    AUC<-ev0@auc
    
    SORENSEN<-2*TPR/(FNR + (2*TPR) + FPR)
  
    
    pred_testing<-rbind(test_p, test_a)
    pca<-scale(pred_testing)
    pred_train_sites_scaled<-scale(train_p, 
                                   center = attr(pca, "scaled:center"), 
                                   scale = attr(pca, "scaled:scale"))
    repeat{
      maha_prob <- suppressWarnings(try(RRdtn:::mahasuhab.custom(as.data.frame(pca), 
                                                         as.data.frame(pred_train_sites_scaled)), 
                                        silent=T))
      if(class(maha_prob)!="try-error")break else {
        pca<-pca[,-ncol(pca),drop=F]
        pred_train_sites_scaled<-pred_train_sites_scaled[,-ncol(pred_train_sites_scaled),drop=F]
        warning(paste("Dimensionality was reduced to", ncol(pca), "axes"))
      }
    }
    
    DDATA<-data.frame(ID=1:nrow(maha_prob),
                      OBS=c(rep(1, nrow(test_p)), rep(0, nrow(test_a))),
                      PRED=maha_prob$MD)
    
    th<-optimal.thresholds(DDATA)[3,2]
    ev0<-presence.absence.accuracy(DDATA, th)
    TPR<-ev0$sensitivity
    TNR<-ev0$specificity
    FPR<-1-TPR
    FNR<-1-TNR
    
    TSS_suit<-TPR+TNR-1
    AUC_suit<-ev0$AUC
    SORENSEN_suit<-2*TPR/(FNR + (2*TPR) + FPR)
    
    CBI<-ecospat.boyce(subset(DDATA, OBS==0)$PRED,
                       subset(DDATA, OBS==1)$PRED,
                       PEplot = F)$cor
    
    outcome<-c(AUC=AUC, AUC_suit=AUC_suit, TSS=TSS, TSS_suit=TSS_suit, SORENSEN=SORENSEN, 
      SORENSEN_suit=SORENSEN_suit, CBI=CBI)
    
  } else {
      DDATA<-data.frame(ID=1:nrow(rbind(test_p,test_a)), 
                        OBS=c(rep(1, nrow(test_p)), rep(0, nrow(test_a))),
                        PRED=rbind(test_p, test_a))
      colnames(DDATA)[3]<-"PRED"
      
      th<-optimal.thresholds(DDATA)[3,2]
      ev0<-presence.absence.accuracy(DDATA, th)
      
      TPR<-ev0$sensitivity
      TNR<-ev0$specificity
      FPR<-1-TPR
      FNR<-1-TNR
      
      TSS<-TPR+TNR-1
      AUC<-ev0$AUC
      SORENSEN<-2*TPR/(FNR + (2*TPR) + FPR)
      
      CBI<-ecospat.boyce(subset(DDATA, OBS==0)$PRED,
                         subset(DDATA, OBS==1)$PRED,
                         PEplot = F)$cor
      
      outcome<-c(AUC=AUC, TSS=TSS, SORENSEN=SORENSEN, CBI=CBI)
    }
  
  return(outcome)
}

jj<-25 ### this sets Rangifer tarandus (the reindeer) only to be modelled. 
for(jj in 25){  ### run line 133 instead of lines 131 - 132 to model all species
###  for(jj in 1:length(DATA_FULL)){  
  dir.create(paste("test_", npoints, "_points/", names(enfa_all_species)[jj], sep=""),
             recursive=T)
  
  setwd(paste("test_", npoints, "_points/", names(enfa_all_species)[jj], sep=""))  
  
  Species_data_start<-DATA_FULL[[jj]]
  Species_data1_start<-subset(Species_data_start, OBS==1)
  Species_data0_start<-subset(Species_data_start, OBS==0)
  print(jj)
 
  list()->extreme_test   
  
  zz0<-replicate(1, sample(nrow(Species_data1_start),npoints)) 
  ### please notice that in our analyses the 
  ### number of reps is set to 20, so that you have to 
  ### type zz0<-replicate(20,...to reproduce them
  
  save(zz0, file="sampled_points.RData")

  
  for(s in 1){ ### run line 155 instead of line 154 to perform all reps
    ### loop for(s in c(1:20)){ in order to perform all 20 reps as in our analyses
    zz<-zz0[,s]
    Species_data1<-Species_data1_start[zz,] 
    Species_data0<-subset(Species_data0_start,TIME_factor%in%Species_data1$TIME_factor) 
    
    Species_data1_for_test<-Species_data1_start[-zz,] 
    Species_data0_for_test<-subset(Species_data0_start,TIME_factor%in%Species_data1_for_test$TIME_factor) 
    
    Species_data<-list(rbind(Species_data1,
                             Species_data0)) 
    names(Species_data)<-names(enfa_all_species)[[jj]]
  
    
    
    enfa_subsampled<-ENphylo_modeling(input_data=Species_data,
                                      tree=tree,
                                      input_mask=MASK_FULL,
                                      obs_col="OBS",
                                      time_col="TIME_factor",
                                      boot_test_perc=20,
                                      boot_reps=10,
                                      min_occ_enfa=npoints-1,
                                      eval_metric_for_imputation ="AUC",
                                      output_options ="best",
                                      eval_threshold=0.5)
    
    
    test_enfa<-enfa_all_species
    test_enfa[[jj]]<-enfa_subsampled[[1]]
    gc()
    
    save(enfa_subsampled, file=paste("enfa_subsampled_rep", s, ".RData", sep=""))
   
    
    enfa_all_species[jj][[1]]$formatted_data<-enfa_subsampled[[1]]$formatted_data
    
  
    test_imputed<-ENphylo_modeling(external_enfa_models=enfa_all_species,
                                   spec_for_imputation=names(enfa_all_species[jj]),
                                   tree=tree,
                                   nsim=10,  #set to 20 in our analyses
                                   si=0.5,
                                   si2=0.5,
                                   eval_metric_for_imputation="AUC",
                                   output_options="full",
                                   eval_threshold=0.7)
  
    
    closeAllConnections()
    gc()
    
    
    test_imputed[jj]->test_imputed
    
    save(test_imputed, file=paste("test_imputed_rep", s, ".RData", sep=""))
    
    
    wm<-which(test_imputed[[1]]$calibrated_model$evaluation$AUC>=0.7)
    
    test_imputed_th<-test_imputed
    test_imputed_th[[1]]$calibrated_model$co[wm]->test_imputed_th[[1]]$calibrated_model$co
    test_imputed_th[[1]]$calibrated_model$evaluation[wm,,drop=F]->test_imputed_th[[1]]$calibrated_model$evaluation
    
    
    save(test_imputed_th, file=paste("test_imputed_th_rep", s, ".RData", sep=""))
    
    
    wm<-which.max(test_imputed[[1]]$calibrated_model$evaluation$AUC)
    
    test_imputed_best<-test_imputed
    test_imputed_best[[1]]$calibrated_model$co[wm]->test_imputed_best[[1]]$calibrated_model$co
    test_imputed_best[[1]]$calibrated_model$evaluation[wm,,drop=F]->test_imputed_best[[1]]$calibrated_model$evaluation
    
    save(test_imputed_best, file=paste("test_imputed_best_rep", s, ".RData", sep=""))
    
    
    
    ## create training dataset  ################################################
    if(nrow(Species_data0)>10000){
      sample(nrow(Species_data0),10000)->k
      Species_data0[k,]->Species_data0
    }
    training_dataset<-rbind(Species_data1,
                            Species_data0)
    training_dataset<-as.data.frame(training_dataset)[,c("OBS","bio4","bio8","bio10","bio13","bio14","bio18")] ## variables selected by VIF
    
    
    ## predict imputed and enfa models on training dataset  ####################
    ENphylo_prediction(object=test_imputed,
                       newdata=training_dataset[,-1])[[1]]$imputed_prediction->pred_training_imputed
    lapply(pred_training_imputed,function(x){
      subset(x,training_dataset$OBS==1)->xx
      xx[,grep("swap",colnames(xx),invert=T)]->xx
    })->pred_train_sites_imputed
    
    
    if(nrow(test_imputed_th[[1]]$calibrated_model$evaluation)>0){
      ENphylo_prediction(object=test_imputed_th,
                         newdata=training_dataset[,-1])[[1]]$imputed_prediction->pred_training_imputed_th
      lapply(pred_training_imputed_th,function(x){
        subset(x,training_dataset$OBS==1)->xx
        xx[,grep("swap",colnames(xx),invert=T)]->xx
      })->pred_train_sites_imputed_th
      
      pred_train_sites_imputed_th<-list(Reduce("+", 
                                               Map("*", pred_train_sites_imputed_th, 
                                                   test_imputed_th[[1]]$calibrated_model$evaluation$AUC))/
                                          length(pred_train_sites_imputed_th))
    }
    
    

    ENphylo_prediction(object=test_imputed_best,
                       newdata=training_dataset[,-1])[[1]]$imputed_prediction->pred_training_imputed_best
    lapply(pred_training_imputed_best,function(x){
      subset(x,training_dataset$OBS==1)->xx
      xx[,grep("swap",colnames(xx),invert=T)]->xx
    })->pred_train_sites_imputed_best
    
    
    ENphylo_prediction(object=test_enfa[jj], 
                       newdata=training_dataset[,-1])[[1]]$enfa_prediction->pred_training_enfa
    pred_train_sites_enfa<-subset(pred_training_enfa, training_dataset$OBS==1)
    pred_train_sites_enfa[,1:test_enfa[[jj]]$calibrated_model$full_model$significant_axes]->pred_train_sites_enfa
    gc()
    

        
    ## create testing dataset  #################################################
    if (nrow(Species_data0_for_test)>10000){
      sample(nrow(Species_data0_for_test),10000)->k
      Species_data0_for_test[k,]->Species_data0_for_test
    }
    testing_dataset<-rbind(Species_data1_for_test,
                           Species_data0_for_test)
    testing_dataset<-as.data.frame(testing_dataset)[,c("OBS","bio4","bio8","bio10","bio13","bio14","bio18")]
    
    
    
    ## predict imputed and enfa models on testing dataset  #####################
    ENphylo_prediction(object=test_imputed, 
                       newdata=testing_dataset[,-1])[[1]]$imputed_prediction->pred_testing_imputed
    
    lapply(pred_testing_imputed,function(x){
      subset(x,testing_dataset$OBS==1)->xx
      xx[,grep("swap",colnames(xx),invert=T)]->xx
    })->pred_test_sites_imputed
    
    lapply(pred_testing_imputed,function(x){
      subset(x,testing_dataset$OBS==0)->xx
      xx[,grep("swap",colnames(xx),invert=T)]->xx
    })->pred_test_background_imputed
    
    
    if(nrow(test_imputed_th[[1]]$calibrated_model$evaluation)>0){
      ENphylo_prediction(object=test_imputed_th, 
                         newdata=testing_dataset[,-1])[[1]]$imputed_prediction->pred_testing_imputed_th
      
      lapply(pred_testing_imputed_th,function(x){
        subset(x,testing_dataset$OBS==1)->xx
        xx[,grep("swap",colnames(xx),invert=T)]->xx
      })->pred_test_sites_imputed_th
      pred_test_sites_imputed_th<-list(Reduce("+", 
                                              Map("*", pred_test_sites_imputed_th, 
                                                  test_imputed_th[[1]]$calibrated_model$evaluation$AUC))/
                                         length(pred_test_sites_imputed_th))
      
      lapply(pred_testing_imputed_th,function(x){
        subset(x,testing_dataset$OBS==0)->xx
        xx[,grep("swap",colnames(xx),invert=T)]->xx
      })->pred_test_background_imputed_th
      pred_test_background_imputed_th<-list(Reduce("+", 
                                                   Map("*", pred_test_background_imputed_th, 
                                                       test_imputed_th[[1]]$calibrated_model$evaluation$AUC))/
                                              length(pred_test_background_imputed_th))
    }
    
    
    
    
    ENphylo_prediction(object=test_imputed_best, 
                       newdata=testing_dataset[,-1])[[1]]$imputed_prediction->pred_testing_imputed_best
    
    lapply(pred_testing_imputed_best,function(x){
      subset(x,testing_dataset$OBS==1)->xx
      xx[,grep("swap",colnames(xx),invert=T)]->xx
    })->pred_test_sites_imputed_best
    
    lapply(pred_testing_imputed_best,function(x){
      subset(x,testing_dataset$OBS==0)->xx
      xx[,grep("swap",colnames(xx),invert=T)]->xx
    })->pred_test_background_imputed_best
    
    
    
    ENphylo_prediction(object=test_enfa[jj], 
                      newdata=testing_dataset[,-1])[[1]]$enfa_prediction->pred_testing_enfa
    pred_testing_enfa[,1:test_enfa[[jj]]$calibrated_model$full_model$significant_axes]->pred_testing_enfa
    
    pred_test_sites_enfa<-subset(pred_testing_enfa, testing_dataset$OBS==1)
    pred_test_background_enfa<-subset(pred_testing_enfa, testing_dataset$OBS==0)
    gc()

    
    
    ##  calculate Mahalanobis distance for training ones only ##################
    lapply(1:length(pred_train_sites_imputed),function(x){
      repeat{
        mm_imputed <- suppressWarnings(try(mahal(pred_train_sites_imputed[[x]]), 
                                           silent=T))
        if(class(mm_imputed)!="try-error")break else {
          pred_train_sites_imputed[[x]]<-pred_train_sites_imputed[[x]][,-ncol(pred_train_sites_imputed[[x]]),drop=F]
          
        }
      }
      mm_imputed
    })->mm_imputed
    
    
    min(sapply(mm_imputed,function(x)ncol(x@presence)))->min_ax
    
    lapply(pred_train_sites_imputed,function(x)x[,1:min_ax])->pred_train_sites_imputed
    lapply(pred_test_sites_imputed,function(x)x[,1:min_ax])->pred_test_sites_imputed
    lapply(pred_test_background_imputed,function(x)x[,1:min_ax])->pred_test_background_imputed
    
    
    if(nrow(test_imputed_th[[1]]$calibrated_model$evaluation)>0){
      lapply(1:length(pred_train_sites_imputed_th),function(x){
        repeat{
          mm_imputed_th <- suppressWarnings(try(mahal(pred_train_sites_imputed_th[[x]]), 
                                                silent=T))
          if(class(mm_imputed_th)!="try-error")break else {
            pred_train_sites_imputed_th[[x]]<-pred_train_sites_imputed_th[[x]][,-ncol(pred_train_sites_imputed_th[[x]]),drop=F]
            
          }
        }
        mm_imputed_th
      })->mm_imputed_th
      
      
      min(sapply(mm_imputed_th,function(x)ncol(x@presence)))->min_ax
      
      lapply(pred_train_sites_imputed_th,function(x)x[,1:min_ax])->pred_train_sites_imputed_th
      lapply(pred_test_sites_imputed_th,function(x)x[,1:min_ax])->pred_test_sites_imputed_th
      lapply(pred_test_background_imputed_th,function(x)x[,1:min_ax])->pred_test_background_imputed_th
    }
    
    
    
    
    lapply(1:length(pred_train_sites_imputed_best),function(x){
      repeat{
        mm_imputed_best <- suppressWarnings(try(mahal(pred_train_sites_imputed_best[[x]]), 
                                           silent=T))
        if(class(mm_imputed_best)!="try-error")break else {
          pred_train_sites_imputed_best[[x]]<-pred_train_sites_imputed_best[[x]][,-ncol(pred_train_sites_imputed_best[[x]]),drop=F]
          
        }
      }
      mm_imputed_best
    })->mm_imputed_best
    
    
    min(sapply(mm_imputed_best,function(x)ncol(x@presence)))->min_ax
    
    lapply(pred_train_sites_imputed_best,function(x)x[,1:min_ax])->pred_train_sites_imputed_best
    lapply(pred_test_sites_imputed_best,function(x)x[,1:min_ax])->pred_test_sites_imputed_best
    lapply(pred_test_background_imputed_best,function(x)x[,1:min_ax])->pred_test_background_imputed_best
    
    
    mm_enfa<-mahal(pred_train_sites_enfa)
    
    
    pblapply(1:length(pred_train_sites_imputed),function(j){
      params_eval<-list(imputed=list(test_p=pred_test_sites_imputed[[j]], 
                                     test_a=pred_test_background_imputed[[j]], 
                                     mm=mm_imputed[[j]],
                                     model_name="imputed",
                                     train_p=pred_train_sites_imputed[[j]]))
      
      
      lapply(params_eval, function(x){
        eval_function(test_p=x$test_p,
                      test_a=x$test_a,
                      mm=x$mm,
                      model_name=x$model_name,
                      train_p=x$train_p)
      })
      
    })->eval_imputed
    do.call(rbind,unlist(eval_imputed, recursive=FALSE))->eval_imputed
    
    
    if(nrow(test_imputed_th[[1]]$calibrated_model$evaluation)>0){
      pblapply(1:length(pred_train_sites_imputed_th),function(j){
        params_eval<-list(imputed_th=list(test_p=pred_test_sites_imputed_th[[j]], 
                                          test_a=pred_test_background_imputed_th[[j]], 
                                          mm=mm_imputed_th[[j]],
                                          model_name="imputed",
                                          train_p=pred_train_sites_imputed_th[[j]]))
        
        
        lapply(params_eval, function(x){
          eval_function(test_p=x$test_p,
                        test_a=x$test_a,
                        mm=x$mm,
                        model_name=x$model_name,
                        train_p=x$train_p)
        })
        
      })->eval_imputed_th
      do.call(rbind,unlist(eval_imputed_th, recursive=FALSE))->eval_imputed_th
    } else eval_imputed_th<-c(AUC=0.5,
                              AUC_suit=0.5,
                              TSS=0,
                              TSS_suit=0,
                              SORENSEN=0,
                              SORENSEN_suit=0,
                              CBI=0)
    

    
    pblapply(1:length(pred_train_sites_imputed_best),function(j){
      params_eval<-list(imputed_best=list(test_p=pred_test_sites_imputed_best[[j]], 
                                     test_a=pred_test_background_imputed_best[[j]], 
                                     mm=mm_imputed_best[[j]],
                                     model_name="imputed",
                                     train_p=pred_train_sites_imputed_best[[j]]))
      
      
      lapply(params_eval, function(x){
        eval_function(test_p=x$test_p,
                      test_a=x$test_a,
                      mm=x$mm,
                      model_name=x$model_name,
                      train_p=x$train_p)
      })
      
    })->eval_imputed_best
    do.call(rbind,unlist(eval_imputed_best, recursive=FALSE))->eval_imputed_best
    
    
    ## ESM CALIBRATION  ########################################################
    
    # Some recommendations before running the next lines:
    # 1) please ensure R and Java have matching architectures (i.e. 32 or 64 bit).
    # 2) please ensure "maxent.jar" file is in the your LibPath dismo/java folder.

    
    pp_maxent_full<-rbind(Species_data1,Species_data0) 
    covars<-grep("bio", colnames(pp_maxent_full), value=T)
    covars<-combn(covars, 2, simplify = F)

    cl<-makeCluster(detectCores()*0.5)
    clusterExport(cl, c("Species_data1",
                        "Species_data0"),
                  envir = environment())

    clusterEvalQ(cl, {
      library(sf)
      library(ecospat)
      library(ENMeval)
      library(biomod2)
      library(caret)
    })

    ESM_settings<-pblapply(covars, function(var){
      pp_maxent<-rbind(Species_data1[,c("OBS", "TIME_factor", var)],
                       Species_data0[,c("OBS", "TIME_factor", var)]) ### training dataset

      occs<-subset(cbind(st_coordinates(pp_maxent), pp_maxent), OBS==1)
      bg<-subset(cbind(st_coordinates(pp_maxent), pp_maxent), OBS==0)

      as.data.frame(pp_maxent)->pp_maxent


      enmeval_results<-ENMevaluate(occs=as.data.frame(occs)[,c("X", "Y", var)],
                                   bg=as.data.frame(bg)[,c("X", "Y", var)],
                                   tune.args=list(fc = c("L", "LQ", "H", "LQH", "LQHP", "LQHPT"),
                                                  rm = seq(0.5, 4, 0.5)),
                                   partitions="none",
                                   algorithm = "maxent.jar",
                                   other.settings = list(abs.auc.diff = FALSE,
                                                         pred.type = "logistic",
                                                         validation.bg = "full"),
                                   parallel=F,
                                   quiet=TRUE)

      best_model_aic<-subset(enmeval_results@results, !is.na(AICc))
      MAXENT<-enmeval_results@tune.settings[as.numeric(rownames(best_model_aic)[which.min(best_model_aic$AICc)]),]


      mydata<-BIOMOD_FormatingData(resp.var=ifelse(pp_maxent$OBS==0, NA, 1),
                                   expl.var=pp_maxent[,var],
                                   resp.name = paste(paste(var,collapse="."),sep="."),
                                   resp.xy = as.data.frame(rbind(occs, bg))[,c("X","Y")])


      trControl<-trainControl(method="boot",
                              p=0.80,
                              number=10,
                              summaryFunction = twoClassSummary,
                              classProbs = T,
                              returnData = F)

      models.options = BIOMOD_ModelingOptions()

      fm <- list()
      GLM.results <- NULL
      i <- 0
      resp <- ifelse(is.na(mydata@data.species), "Absence", "Presence")
      for (type in c("simple", "quadratic", "polynomial")) {
        for (IA in c(0,1)) {
          i <- i + 1
          try(tune.GLM <- caret::train(makeFormula("resp",
                                                   mydata@data.env.var,
                                                   type = type, interaction.level = IA),
                                       data = cbind(mydata@data.env.var, resp = resp),
                                       method = "glmStepAIC", trControl = trControl,
                                       weights = NULL))
          try(GLM.results <- rbind(GLM.results, cbind(tune.GLM$results,
                                                      il = IA, type = type)))
          try(fm[[i]] <- formula(tune.GLM$finalModel))
        }
      }
      glm.best <- which.max(GLM.results$ROC)
      models.options@GLM$interaction.level <- GLM.results[glm.best,
                                                          "il"]
      models.options@GLM$type <- as.character(GLM.results[glm.best,
                                                          "type"])



      models.options@MAXENT.Phillips$betamultiplier<-as.numeric(MAXENT$rm)

      if(MAXENT$fc=="L"){
        models.options@MAXENT.Phillips$linear<-T
        models.options@MAXENT.Phillips$quadratic<-F
        models.options@MAXENT.Phillips$product<-F
        models.options@MAXENT.Phillips$threshold<-F
        models.options@MAXENT.Phillips$hinge<-F
      }
      if(MAXENT$fc=="LQ"){
        models.options@MAXENT.Phillips$linear<-T
        models.options@MAXENT.Phillips$quadratic<-T
        models.options@MAXENT.Phillips$product<-F
        models.options@MAXENT.Phillips$threshold<-F
        models.options@MAXENT.Phillips$hinge<-F
      }
      if(MAXENT$fc=="H"){
        models.options@MAXENT.Phillips$linear<-F
        models.options@MAXENT.Phillips$quadratic<-F
        models.options@MAXENT.Phillips$product<-F
        models.options@MAXENT.Phillips$threshold<-F
        models.options@MAXENT.Phillips$hinge<-T
      }
      if(MAXENT$fc=="LQH"){
        models.options@MAXENT.Phillips$linear<-T
        models.options@MAXENT.Phillips$quadratic<-T
        models.options@MAXENT.Phillips$product<-F
        models.options@MAXENT.Phillips$threshold<-F
        models.options@MAXENT.Phillips$hinge<-T
      }
      if(MAXENT$fc=="LQHP"){
        models.options@MAXENT.Phillips$linear<-T
        models.options@MAXENT.Phillips$quadratic<-T
        models.options@MAXENT.Phillips$product<-T
        models.options@MAXENT.Phillips$threshold<-F
        models.options@MAXENT.Phillips$hinge<-T
      }
      if(MAXENT$fc=="LQHPT"){
        models.options@MAXENT.Phillips$linear<-T
        models.options@MAXENT.Phillips$quadratic<-T
        models.options@MAXENT.Phillips$product<-T
        models.options@MAXENT.Phillips$threshold<-T
        models.options@MAXENT.Phillips$hinge<-T
      }

      return(models.options)
    }, cl=cl)

    stopCluster(cl)
    closeAllConnections()
    gc()


    ESM_settings<-lapply(ESM_settings, function(x){
      list(MAXENT=data.frame(linear=x@MAXENT.Phillips$linear,
                             quadratic=x@MAXENT.Phillips$quadratic,
                             product=x@MAXENT.Phillips$product,
                             threshold=x@MAXENT.Phillips$threshold,
                             hinge=x@MAXENT.Phillips$hinge,
                             betamultiplier=x@MAXENT.Phillips$betamultiplier),
           GLM=data.frame(type=x@GLM$type,
                          int=x@GLM$interaction.level))

    })

    MAXENT_settings<-do.call(rbind, lapply(ESM_settings, "[[", 1))
    GLM_settings<-do.call(rbind, lapply(ESM_settings, "[[", 2))

    GLM_settings<-apply(GLM_settings, 2, function(x){
      ll<-as.data.frame(table(x))
      ll[which.max(ll$Freq),1]
    })
    MAXENT_settings<-apply(MAXENT_settings, 2, function(x){
      ll<-as.data.frame(table(x))
      ll[which.max(ll$Freq),1]
    })
    MAXENT_settings<-as.data.frame(t(as.data.frame(MAXENT_settings)))
    MAXENT_settings[1,1:5]<-MAXENT_settings[1,1:5]==1

    myBiomodOption<-Print_Default_ModelingOptions()
    myBiomodOption@GLM$test<-"none"
    myBiomodOption@GLM$type<-as.character(GLM_settings[1])
    myBiomodOption@GLM$interaction.level<-as.numeric(as.character(GLM_settings[2]))

    myBiomodOption@MAXENT.Phillips$path_to_maxent.jar=paste(installed.packages()["dismo",2],"dismo","java",sep="/") # please ensure this path contains "maxent.jar"
    myBiomodOption@MAXENT.Phillips$linear=as.logical(MAXENT_settings$linear)
    myBiomodOption@MAXENT.Phillips$quadratic=as.logical(MAXENT_settings$quadratic)
    myBiomodOption@MAXENT.Phillips$product=as.logical(MAXENT_settings$product)
    myBiomodOption@MAXENT.Phillips$threshold=as.logical(MAXENT_settings$threshold)
    myBiomodOption@MAXENT.Phillips$hinge=as.logical(MAXENT_settings$hinge)
    myBiomodOption@MAXENT.Phillips$betamultiplier=as.numeric(MAXENT_settings$betamultiplier)


    occs<-subset(cbind(st_coordinates(pp_maxent_full), pp_maxent_full), OBS==1)
    bg<-subset(cbind(st_coordinates(pp_maxent_full), pp_maxent_full), OBS==0)


    mydata<-BIOMOD_FormatingData(resp.var=ifelse(pp_maxent_full$OBS==0, NA, 1),
                                 expl.var=as.data.frame(pp_maxent_full)[,grep("bio",
                                                                              colnames(pp_maxent_full))],
                                 resp.name = paste("rep", s, sep="."),
                                 resp.xy = as.data.frame(rbind(occs, bg))[,c("X","Y")])


    mymodel<-ecospat.ESM.Modeling_CUSTOM(data=mydata,
                                         DataSplit = 80,
                                         NbRunEval = 10,
                                         DataSplitTable=NULL,
                                         models=c("GLM","RF","MAXENT.Phillips"),
                                         weighting.score = "AUC",
                                         models.options = myBiomodOption,
                                         cleanup = T,
                                         parallel = T,
                                         ncores = detectCores()*0.5)


    gc()

    myEM<-ecospat.ESM.EnsembleModeling_CUSTOM(ESM.modeling.output=mymodel,
                                              weighting.score = "AUC",
                                              threshold = 0.7,
                                              models=c("GLM", "RF", "MAXENT.Phillips"))

    if(any(myEM$weights.EF[,2]!=0)){
      
      ## preparation testing dataset for Ensemble bivariate
      myproj<-ecospat.ESM.Projection_CUSTOM(mymodel,
                                            proj.name=paste("rep", npoints, sep="."),
                                            new.env = testing_dataset,
                                            parallel=T,
                                            ncores=10)

      gc()

      source(paste(main.dir,"ecospat.ESM.EnsembleProjection_CUSTOM.R",sep="/"))
      
      r<-myEF$EF/1000

      pred_test_background_ENS<-data.frame(EF=subset(r, testing_dataset$OBS==0))
      pred_test_sites_ENS<-data.frame(EF=subset(r, testing_dataset$OBS==1))

      gc()

    #   #### EVALUATION ############################################################
      params_eval<-list(enfa=list(test_p=pred_test_sites_enfa,
                                  test_a=pred_test_background_enfa,
                                  mm=mm_enfa,
                                  model_name="enfa",
                                  train_p=pred_train_sites_enfa),
                        ens.biv=list(test_p=pred_test_sites_ENS,
                                     test_a=pred_test_background_ENS,
                                     train_p=NULL,
                                     mm=NULL,
                                     model_name="ens.biv"))


      lapply(params_eval, function(x){
        eval_function(test_p=x$test_p,
                      test_a=x$test_a,
                      mm=x$mm,
                      model_name=x$model_name,
                      train_p=x$train_p)
      })->eval_others

    } else {
      params_eval<-list(enfa=list(test_p=pred_test_sites_enfa,
                                  test_a=pred_test_background_enfa,
                                  mm=mm_enfa,
                                  model_name="enfa",
                                  train_p=pred_train_sites_enfa))


      lapply(params_eval, function(x){
        eval_function(test_p=x$test_p,
                      test_a=x$test_a,
                      mm=x$mm,
                      model_name=x$model_name,
                      train_p=x$train_p)
      })->eval_others

      eval_others$ens.biv<-c(AUC=0.5,
                             TSS=0,
                             SORENSEN=0,
                             CBI=0)
    }

    

    extreme_test[[s]]<-list(enfa=eval_others$enfa,
                            ens.biv=eval_others$ens.biv,
                            imputed.full=eval_imputed,
                            imputed.th=eval_imputed_th,
                            imputed.best=eval_imputed_best)
    
    
   }
  

  
  save(extreme_test, file="external_evaluation.RData")
  
  
  setwd(main.dir)
  
} 

