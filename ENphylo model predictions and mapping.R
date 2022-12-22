### Installing/Updating packages
required<-c("ape", "raster", "methods","CENFA","pbapply","Rphylopars","RRphylo","dismo","gtools",
            "adehabitatMA", "ecospat", "foreach", "doParallel", "PresenceAbsence", 
            "parallel", "ade4", "sp", "biomod2", "terra", "sf", "adehabitatHS")
if(any(!required%in%installed.packages()[,1]))
  install.packages(required[which(!required%in%installed.packages()[,1])])

latest_version<-c("ecospat","dismo","biomod2","terra")
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
MASK_FULL<-raster::raster("variable_bio1.tif")
external_data<-raster::stack(list.files("external_data",full.names=TRUE))
install.packages("RRdtn_0.2.0.tar.gz", repos = NULL, type = "source", INSTALL_opts = "--install-tests")
require(RRdtn)

### Running the analyses
#the species is indicated by the numeric value, here 'k=25' in the following line (available 1 - 31)
k=25

npoints=10 ### for subsampling the k species

Species_data_start<-DATA_FULL[[k]]
Species_data1_start<-subset(Species_data_start, OBS==1)
Species_data0_start<-subset(Species_data_start, OBS==0)
sample(nrow(Species_data1_start),npoints)->zz
Species_data1<-Species_data1_start[zz,]
Species_data0<-subset(Species_data0_start,TIME_factor%in%Species_data1$TIME_factor)
Species_data<-rbind(Species_data1,Species_data0)
Species_data->DATA_FULL[[k]]

##### ENFA ######
enfa_test<-ENphylo_modeling(input_data=DATA_FULL[k],
                            input_mask=MASK_FULL,
                            tree=tree,
                            obs_col="OBS",
                            time_col="TIME_factor",
                            min_occ_enfa=npoints-1,
                            eval_metric_for_imputation="AUC",
                            output_options = "best",
                            eval_threshold=0.5)

gc()


enfa_all_species[[k]]$formatted_data<-enfa_test[[1]]$formatted_data


##### ENphylo full model ######
enphylo_test1<-ENphylo_modeling(external_enfa_models=enfa_all_species,
                                spec_for_imputation=names(enfa_all_species)[k],
                                tree=tree,
                                nsim=50,
                                si=0.5,
                                si2=0.5,
                                eval_metric_for_imputation="AUC",
                                output_options="full",
                                eval_threshold=0.70)

gc()
closeAllConnections()

### check average model performance under "full" strategy ###
apply(enphylo_test1[[k]]$calibrated_model$evaluation,2,mean)  

##### Mapping ENphylo predictions; committee average of the full models #####
enphylo_pred1<-ENphylo_prediction(object=enphylo_test1[k],
                                  newdata=external_data,
                                  convert.to.suitability=TRUE,
                                  calc_method="terra")

grep("MaxSensSpec",names(enphylo_pred1[[1]]$imputed_prediction))->l
terra::rast(sum(enphylo_pred1[[1]]$imputed_prediction[[l]]))->rr
terra::vect(enphylo_test1[[k]]$formatted_data$ones_coords)->pp
terra::plot(rr)
plot(pp,add=TRUE,pch=3,col="red")



##### ENphylo best selection model #####
enphylo_test2<-ENphylo_modeling(external_enfa_models=enfa_all_species,
                                spec_for_imputation=names(enfa_all_species)[k],
                                tree=tree,
                                nsim=50,
                                si=0.5,
                                si2=0.5,
                                eval_metric_for_imputation="AUC",
                                output_options="best",
                                eval_threshold=0.70)

gc()
closeAllConnections()

### check average model performance under "best" strategy ###
apply(enphylo_test2[[k]]$calibrated_model$evaluation,2,mean)  


##### Mapping ENphylo prediction; best model #####
enphylo_pred2<-ENphylo_prediction(object=enphylo_test2[k],
                                  newdata=external_data,
                                  convert.to.suitability=TRUE,
                                  calc_method="terra")


grep("Suitability",names(enphylo_pred2[[1]]$imputed_prediction))->l
terra::rast(enphylo_pred2[[1]]$imputed_prediction[[l]])->rr
terra::vect(enphylo_test2[[k]]$formatted_data$ones_coords)->pp
terra::plot(rr)
plot(pp,add=TRUE,pch=3,col="red")



