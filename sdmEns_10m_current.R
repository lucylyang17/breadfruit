# load library
library(biomod2)
library(raster)
library(dplyr)
library(ggplot2)
library(usdm)


artocarpus <- read.csv("/projects/b1045/yang/gbif_points_10m.csv")



# the name
myRespName <- 'Artocarpus altilis'
# the XY coordinates of the presence
myRespXY <- artocarpus[,c(1,2)]
# and the presence data
myResp <-  as.numeric(artocarpus[,3])


# load the environmental raster layers from local computer (currently, all 19 bioclim variables)


r = stack( "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_1.asc",
           "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_2.asc",
           "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_3.asc",
           "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_4.asc",
           "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_5.asc",
           "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_6.asc",
           "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_7.asc",
           "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_8.asc",
           "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_9.asc",
           "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_10.asc",
           "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_11.asc",
           "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_12.asc",
           "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_13.asc",
           "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_14.asc",
           "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_15.asc",
           "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_16.asc",
           "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_17.asc",
           "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_18.asc",
           "/projects/b1045/yang/wc21/10m/wc2.1_10m_bio_19.asc")
names(r) <- c('bio1', 'bio2', 'bio3', 'bio4', 'bio5', 'bio6', 'bio7', 'bio8', 'bio9', 'bio10', 'bio11', 'bio12', 'bio13', 'bio14', 'bio15', 'bio16', 'bio17', 'bio18', 'bio19')


  
  myExpl <- r[[c(2,7,11,12,17:19)]]
  e <- extent(-180, 180, -45, 45)
  myExpl <- crop(myExpl, e)
  myExpl <- stack(myExpl)
  

  
    
#need to create pseudoabsence points - random
#test with biomod
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep = 3,
                                     PA.nb.absences = 200,
                                     PA.strategy = 'random')

get_PAtab <- function(bfd){
  dplyr::bind_cols(
    x = bfd@coord[, 1],
    y = bfd@coord[, 2],
    status = bfd@data.species,
    bfd@PA
  )
}

myBiomodData


#extract pseudoabsences, set to 0, add to DataSpecies dataframe with all ground-points
pres.xy <- get_PAtab(myBiomodData)
m <- as.data.frame(pres.xy, xy=TRUE)
m[is.na(m)] <- 0

write.csv(m, "PAloc.csv", row.names = TRUE)

env <- extract(myExpl, m[,c(1,2)])
DataSpecies <- cbind(m[,1:3], env)
#DataSpecies

# Defining Models Options using default options.
myBiomodOption <- BIOMOD_ModelingOptions()


# Computing the models
#early testing only look at GLM and GAM
myBiomodModelOut <- BIOMOD_Modeling(
  myBiomodData,
  models = c('MAXENT.Phillips.2','GLM','GAM','MARS','RF','GBM'),
  models.options = myBiomodOption,
  NbRunEval=10,
  DataSplit=70,
  Yweights=NULL,
  VarImport=3,
  models.eval.meth = c('TSS','ROC','KAPPA'),
  SaveObj = TRUE,
  rescal.all.models = TRUE)




myBiomodModelOut
get_evaluations(myBiomodModelOut)
modelscores <- get_evaluations(myBiomodModelOut, as.data.frame=TRUE)
write.csv(modelscores, paste("/projects/b1045/yang/PA",x1,".csv",sep=""), row.names = TRUE)
#dimnames(modelscores)

get_variables_importance(myBiomodModelOut)

##BOX PLOTS
modelscores
##TSS
scoresTSS <- modelscores[ which(modelscores$Eval.metric=='TSS'),]
scoresTSS$model <- substr(scoresTSS$Model.name, 0, 3)
##ROC
scoresROC <- modelscores[ which(modelscores$Eval.metric=='ROC'),]
scoresROC$model <- substr(scoresROC$Model.name, 0, 3)
##KAPPA
scoresKAPPA <- modelscores[ which(modelscores$Eval.metric=='KAPPA'),]
scoresKAPPA$model <- substr(scoresKAPPA$Model.name, 0, 3)

#PLOTS
pdf(paste("boxplot",x1,".pdf",sep=""))
p1 <- ggplot(scoresTSS, aes(model,Testing.data)) + geom_boxplot() + ggtitle("TSS")
p2 <- ggplot(scoresROC, aes(model,Testing.data)) + geom_boxplot() + ggtitle("ROC")
p3 <- ggplot(scoresKAPPA, aes(model,Testing.data)) + geom_boxplot() + ggtitle("KAPPA")
print(p1)
print(p2)
print(p3)
dev.off()


#Seeing variation of parameters
pdf(paste("modelscores",x1,".pdf",sep=""))
models_scores_graph(myBiomodModelOut, by = "models", metrics = c("ROC","TSS"), xlim = c(0,1), ylim=c(0,1))
models_scores_graph(myBiomodModelOut, by = "cv_run", metrics = c("ROC","TSS"), xlim = c(0,1), ylim=c(0,1))
models_scores_graph(myBiomodModelOut, by = "data_set", metrics = c("ROC","TSS"), xlim = c(0,1), ylim=c(0,1))
dev.off()


#Will this work with multiple???
#may need to spit out into csv file
#var_import <- 
#get_variables_importance(myBiomodModelOut)


##RESPONSE CURVES
#loadmodels
model_glm <- BIOMOD_LoadModels(myBiomodModelOut, models = 'GLM')
model_gam <- BIOMOD_LoadModels(myBiomodModelOut, models = 'GAM')
model_maxent <- BIOMOD_LoadModels(myBiomodModelOut, models = 'MAXENT.Phillips.2')
model_mars <- BIOMOD_LoadModels(myBiomodModelOut, models = 'MARS')
model_rf <- BIOMOD_LoadModels(myBiomodModelOut, models = 'RF')
model_gbm <- BIOMOD_LoadModels(myBiomodModelOut, models = 'GBM')

pdf(paste("responsecurves",x1,".pdf",sep=""))
glm_eval_strip <- biomod2::response.plot2(
  models = model_glm,
  Data = get_formal_data(myBiomodModelOut,'expl.var'), 
  show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
#  legend = FALSE,
#  display_title = FALSE,
  data_species = get_formal_data(myBiomodModelOut,'resp.var'))

gam_eval_strip <- biomod2::response.plot2(
  models = model_gam,
  Data = get_formal_data(myBiomodModelOut,'expl.var'), 
  show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  #  legend = FALSE,
  #  display_title = FALSE,
  data_species = get_formal_data(myBiomodModelOut,'resp.var'))

maxent_eval_strip <- biomod2::response.plot2(
  models = model_maxent,
  Data = get_formal_data(myBiomodModelOut,'expl.var'),
  show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  #  legend = FALSE,
  #  display_title = FALSE,
  data_species = get_formal_data(myBiomodModelOut,'resp.var'))

mars_eval_strip <- biomod2::response.plot2(
  models = model_mars,
  Data = get_formal_data(myBiomodModelOut,'expl.var'),
  show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  #  legend = FALSE,
  #  display_title = FALSE,
  data_species = get_formal_data(myBiomodModelOut,'resp.var'))

rf_eval_strip <- biomod2::response.plot2(
  models = model_rf,
  Data = get_formal_data(myBiomodModelOut,'expl.var'),
  show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  #  legend = FALSE,
  #  display_title = FALSE,
  data_species = get_formal_data(myBiomodModelOut,'resp.var'))

gbm_eval_strip <- biomod2::response.plot2(
  models = model_gbm,
  Data = get_formal_data(myBiomodModelOut,'expl.var'),
  show.variables = get_formal_data(myBiomodModelOut,'expl.var.names'),
  do.bivariate = FALSE,
  fixed.var.metric = 'median',
  #  legend = FALSE,
  #  display_title = FALSE,
  data_species = get_formal_data(myBiomodModelOut,'resp.var'))



dev.off()
###PROJECTION
#project the potential distribution of the species over space and time
# projection over the globe under current conditions
myBiomomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl,
  proj.name = 'current_wmax',
  selected.models = 'all',
  binary.meth = c('TSS','ROC','KAPPA'),
  compress = 'xz',
  clamping.mask = F,
  output.format = '.grd')

myBiomomodProj
current <- get_predictions(myBiomomodProj)
#current
#names(current)


GLMstack <- raster::subset(current, grep('Full_GLM', names(current), value = TRUE))
GLMavg <-  calc(GLMstack, fun = mean, na.rm = T)

GAMstack <- raster::subset(current, grep('Full_GAM', names(current), value = TRUE))
GAMavg <- calc(GAMstack, fun = mean, na.rm = T)

MAXENTstack <- raster::subset(current, grep('Full_MAXENT', names(current), value = TRUE))
MAXENTavg <- calc(MAXENTstack, fun = mean, na.rm = T)

MARSstack <- raster::subset(current, grep('Full_MARS', names(current), value = TRUE))
MARSavg <- calc(MARSstack, fun = mean, na.rm = T)

RFstack <- raster::subset(current, grep('Full_RF', names(current), value = TRUE))
RFavg <- calc(RFstack, fun = mean, na.rm = T)

GBMstack <- raster::subset(current, grep('Full_GBM', names(current), value = TRUE))
GBMavg <- calc(GBMstack, fun = mean, na.rm = T)

#GLMstack <- subset(current,c(21,43,65))
#GAMstack <- subset(current,c(22,44,66))
#GLMavg <-  calc(GLMstack, fun = mean, na.rm = T)
#GAMavg <- calc(GAMstack, fun = mean, na.rm = T)

pdf(paste("MAP",x1,".pdf",sep=""))
plot(myBiomomodProj, str.grep = 'Full_GLM')
plot(GLMavg)

plot(myBiomomodProj, str.grep = 'Full_GAM')
plot(GAMavg)

plot(myBiomomodProj, str.grep = 'Full_MAXENT')
plot(MAXENTavg)


plot(myBiomomodProj, str.grep = 'Full_MARS')
plot(MARSavg)

plot(myBiomomodProj, str.grep = 'Full_RF')
plot(RFavg)

plot(myBiomomodProj, str.grep = 'Full_GBM')
plot(GBMavg)

dev.off()

GLMscores <- extract(GLMavg, m[,c(1,2)])
GLMscores <- as.data.frame(GLMscores)
GLMscores <- cbind(m,GLMscores)

GAMscores <- extract(GAMavg, m[,c(1,2)])
GAMscores <- as.data.frame(GAMscores)
GAMscores <- cbind(m,GAMscores)

MAXENTscores <- extract(MAXENTavg, m[,c(1,2)])
MAXENTscores <- as.data.frame(MAXENTscores)
MAXENTscores <- cbind(m,MAXENTscores)

MARSscores <- extract(MARSavg, m[,c(1,2)])
MARSscores <- as.data.frame(MARSscores)
MARSscores <- cbind(m,MARSscores)

RFscores <- extract(RFavg, m[,c(1,2)])
RFscores <- as.data.frame(RFscores)
RFscores <- cbind(m,RFscores)

GBMscores <- extract(GBMavg, m[,c(1,2)])
GBMscores <- as.data.frame(GBMscores)
GBMscores <- cbind(m,GBMscores)


GLMscoressorted <- GLMscores
GLMscoressorted <- GLMscoressorted[order(GLMscoressorted$GLMscores),]
GLMscoressorted$observation <- 1:nrow(GLMscores)
GLMscoressorted$cumper <- GLMscoressorted$observation/nrow(GLMscoressorted)


GAMscoressorted <- GAMscores
GAMscoressorted <- GAMscoressorted[order(GAMscoressorted$GAMscores),]
GAMscoressorted$observation <- 1:nrow(GAMscores)
GAMscoressorted$cumper <- GAMscoressorted$observation/nrow(GAMscoressorted)

MAXENTscoressorted <- MAXENTscores
MAXENTscoressorted <- MAXENTscoressorted[order(MAXENTscoressorted$MAXENTscores),]
MAXENTscoressorted$observation <- 1:nrow(MAXENTscores)
MAXENTscoressorted$cumper <- MAXENTscoressorted$observation/nrow(MAXENTscoressorted)

MARSscoressorted <- MARSscores
MARSscoressorted <- MARSscoressorted[order(MARSscoressorted$MARSscores),]
MARSscoressorted$observation <- 1:nrow(MARSscores)
MARSscoressorted$cumper <- MARSscoressorted$observation/nrow(MARSscoressorted)

RFscoressorted <- RFscores
RFscoressorted <- RFscoressorted[order(RFscoressorted$RFscores),]
RFscoressorted$observation <- 1:nrow(RFscores)
RFscoressorted$cumper <- RFscoressorted$observation/nrow(RFscoressorted)

GBMscoressorted <- GBMscores
GBMscoressorted <- GBMscoressorted[order(GBMscoressorted$GBMscores),]
GBMscoressorted$observation <- 1:nrow(GBMscores)
GBMscoressorted$cumper <- GBMscoressorted$observation/nrow(GBMscoressorted)

pdf(paste("cum",x1,".pdf",sep=""))
plot(GLMscoressorted$GLMscores,GLMscoressorted$cumper,main="GLM",xlab="Model Value",ylab="Cumulative Percent")
plot(GAMscoressorted$GAMscores,GAMscoressorted$cumper, main="GAM",xlab="Model Value",ylab="Cumulative Percent")
plot(MAXENTscoressorted$MAXENTscores,MAXENTscoressorted$cumper, main="MAXENT",xlab="Model Value",ylab="Cumulative Percent")

plot(MARSscoressorted$MARSscores,MARSscoressorted$cumper, main="MARS",xlab="Model Value",ylab="Cumulative Percent")

plot(RFscoressorted$RFscores,RFscoressorted$cumper, main="RF",xlab="Model Value",ylab="Cumulative Percent")

plot(GBMscoressorted$GBMscores,GBMscoressorted$cumper, main="GBM",xlab="Model Value",ylab="Cumulative Percent")

dev.off()

##ANOVA
print("GLM")
res.aov <- aov(GLMscores ~ status, data = GLMscores)
aovGLM <- summary(res.aov)
print(aovGLM)
print("GAM")
res.aovGAM <- aov(GAMscores ~ status, data = GAMscores)
aovGAM <- summary(res.aovGAM)
print(aovGAM)

print("MAXENT")
res.aovMAXENT <- aov(MAXENTscores ~ status, data = MAXENTscores)
aovMAXENT <- summary(res.aovMAXENT)
print(aovMAXENT)

print("MARS")
res.aovMARS <- aov(MARSscores ~ status, data = MARSscores)
aovMARS <- summary(res.aovMARS)
print(aovMARS)


print("RF")
res.aovRF <- aov(RFscores ~ status, data = RFscores)
aovRF <- summary(res.aovRF)
print(aovRF)

print("GBM")
res.aovGBM <- aov(GBMscores ~ status, data = GBMscores)
aovGBM <- summary(res.aovGBM)
print(aovGBM)

#EM
myBiomodEM <- BIOMOD_EnsembleModeling(
  modeling.output = myBiomodModelOut,
  chosen.models = 'all',
  eval.metric = c('TSS'),
  eval.metric.quality.threshold = c(0.7),
  prob.mean = T,
  prob.cv = T,
  prob.ci = T,
  prob.ci.alpha = 0.05,
  prob.median = T,
  committee.averaging = T,
  prob.mean.weight = T,
  prob.mean.weight.decay = 'proportional' )

myBiomodEM
getEMeval(myBiomodEM)


# make some plots sub-selected by str.grep argument
pdf("individual.pdf")

plot(myBiomomodProj, str.grep = 'GLM')
plot(myBiomomodProj, str.grep = 'GAM')
plot(myBiomomodProj, str.grep = 'GBM')
plot(myBiomomodProj, str.grep = 'MARS')
plot(myBiomomodProj, str.grep = 'RF')
plot(myBiomomodProj, str.grep = 'MAXENT.Phillips.2')



dev.off()


myBiomodEMF <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomomodProj)






myExplBCC45 = stack( "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp245_2061_2080_2.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp245_2061_2080_7.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp245_2061_2080_11.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp245_2061_2080_12.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp245_2061_2080_17.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp245_2061_2080_18.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp245_2061_2080_19.asc")

names(myExplBCC45) <- c('bio2', 'bio7', 'bio11', 'bio12', 'bio17', 'bio18', 'bio19')


myBiomodProjBCC45 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExplBCC45,
  proj.name = '2070_SDM_4_5_BCC_wmax',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

myExplBCC85 = stack( "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp585_2061_2080_2.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp585_2061_2080_7.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp585_2061_2080_11.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp585_2061_2080_12.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp585_2061_2080_17.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp585_2061_2080_18.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp585_2061_2080_19.asc")

names(myExplBCC85) <- c('bio2', 'bio7', 'bio11', 'bio12', 'bio17', 'bio18', 'bio19')


myBiomodProjBCC85 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExplBCC85,
  proj.name = '2070_SDM_8_5_BCC_wmax',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')




#FUTURE FORECASTING
myBiomodEF_BCC45 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjBCC45
)

myBiomodEF_BCC45


pdf("2070_4_5_BCC_wmax.pdf")
plot(myBiomodEF_BCC45)
dev.off()

myBiomodEF_BCC85 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjBCC85
)

myBiomodEF_BCC85


pdf("2070_8_5_BCC_wmax.pdf")
plot(myBiomodEF_BCC85)
dev.off()


#####CanESM
myExpl45 = stack( "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CanESM5_ssp245_2061_2080_2.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CanESM5_ssp245_2061_2080_7.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CanESM5_ssp245_2061_2080_11.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CanESM5_ssp245_2061_2080_12.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CanESM5_ssp245_2061_2080_17.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CanESM5_ssp245_2061_2080_18.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CanESM5_ssp245_2061_2080_19.asc")

names(myExpl45) <- c('bio2', 'bio7', 'bio11', 'bio12', 'bio17', 'bio18', 'bio19')


myBiomodProjCan45 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl45,
  proj.name = '2070_SDM_4_5_Can_wmax',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

myExplCan85 = stack( "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CanESM5_ssp585_2061_2080_2.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CanESM5_ssp585_2061_2080_7.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CanESM5_ssp585_2061_2080_11.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CanESM5_ssp585_2061_2080_12.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CanESM5_ssp585_2061_2080_17.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CanESM5_ssp585_2061_2080_18.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CanESM5_ssp585_2061_2080_19.asc")

names(myExplCan85) <- c('bio2', 'bio7', 'bio11', 'bio12', 'bio17', 'bio18', 'bio19')


myBiomodProjCan85 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExplCan85,
  proj.name = '2070_SDM_8_5_Can_wmax',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

#FUTURE FORECASTING
myBiomodEF_Can45 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjCan45
)

myBiomodEF_Can45


pdf("2070_4_5_Can.pdf")
plot(myBiomodEF_Can45)
dev.off()

myBiomodEF_Can85 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjCan85
)

myBiomodEF_Can85


pdf("2070_8_5_Can.pdf")
plot(myBiomodEF_Can85)
dev.off()


#####CM
myExpl45 = stack( "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_CM6_1_ssp245_2061_2080_2.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_CM6_1_ssp245_2061_2080_7.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_CM6_1_ssp245_2061_2080_11.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_CM6_1_ssp245_2061_2080_12.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_CM6_1_ssp245_2061_2080_17.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_CM6_1_ssp245_2061_2080_18.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_CM6_1_ssp245_2061_2080_19.asc")

names(myExpl45) <- c('bio2', 'bio7', 'bio11', 'bio12', 'bio17', 'bio18', 'bio19')


myBiomodProjCM45 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl45,
  proj.name = '2070_SDM_4_5_CM_wmax',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

myExplCM85 = stack( "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_CM6_1_ssp585_2061_2080_2.asc",
                    "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_CM6_1_ssp585_2061_2080_7.asc",
                    "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_CM6_1_ssp585_2061_2080_11.asc",
                    "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_CM6_1_ssp585_2061_2080_12.asc",
                    "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_CM6_1_ssp585_2061_2080_17.asc",
                    "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_CM6_1_ssp585_2061_2080_18.asc",
                    "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_CM6_1_ssp585_2061_2080_19.asc")

names(myExplCM85) <- c('bio2', 'bio7', 'bio11', 'bio12', 'bio17', 'bio18', 'bio19')


myBiomodProjCM85 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExplCM85,
  proj.name = '2070_SDM_8_5_CM_wmax',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

#FUTURE FORECASTING
myBiomodEF_CM45 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjCM45
)

myBiomodEF_CM45



pdf("2070_4_5_CM.pdf")
plot(myBiomodEF_CM45)
dev.off()

myBiomodEF_CM85 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjCM85
)

myBiomodEF_CM85



pdf("2070_8_5_CM.pdf")
plot(myBiomodEF_CM85)
dev.off()

###ESM
myExpl45 = stack( "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_ESM2_1_ssp245_2061_2080_2.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_ESM2_1_ssp245_2061_2080_7.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_ESM2_1_ssp245_2061_2080_11.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_ESM2_1_ssp245_2061_2080_12.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_ESM2_1_ssp245_2061_2080_17.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_ESM2_1_ssp245_2061_2080_18.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_ESM2_1_ssp245_2061_2080_19.asc")

names(myExpl45) <- c('bio2', 'bio7', 'bio11', 'bio12', 'bio17', 'bio18', 'bio19')


myBiomodProjESM45 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl45,
  proj.name = '2070_SDM_4_5_ESM_wmax',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

myExplESM85 = stack( "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_ESM2_1_ssp585_2061_2080_2.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_ESM2_1_ssp585_2061_2080_7.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_ESM2_1_ssp585_2061_2080_11.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_ESM2_1_ssp585_2061_2080_12.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_ESM2_1_ssp585_2061_2080_17.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_ESM2_1_ssp585_2061_2080_18.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_CNRM_ESM2_1_ssp585_2061_2080_19.asc")

names(myExplESM85) <- c('bio2', 'bio7', 'bio11', 'bio12', 'bio17', 'bio18', 'bio19')


myBiomodProjESM85 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExplESM85,
  proj.name = '2070_SDM_8_5_ESM_wmax',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

#FUTURE FORECASTING
myBiomodEF_ESM45 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjESM45
)

myBiomodEF_ESM45



pdf("2070_4_5_ESM.pdf")
plot(myBiomodEF_ESM45)
dev.off()

myBiomodEF_ESM85 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjESM85
)

myBiomodEF_ESM85

#writeRaster(myBiomodEF, filename="future_SDM.grd", overwrite=TRUE, format="raster")


pdf("2070_8_5_ESM.pdf")
plot(myBiomodEF_ESM85)
dev.off()

###IPSL
myExpl45 = stack( "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_IPSL_CM6A_LR_ssp245_2061_2080_2.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_IPSL_CM6A_LR_ssp245_2061_2080_7.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_IPSL_CM6A_LR_ssp245_2061_2080_11.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_IPSL_CM6A_LR_ssp245_2061_2080_12.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_IPSL_CM6A_LR_ssp245_2061_2080_17.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_IPSL_CM6A_LR_ssp245_2061_2080_18.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_IPSL_CM6A_LR_ssp245_2061_2080_19.asc")

names(myExpl45) <- c('bio2', 'bio7', 'bio11', 'bio12', 'bio17', 'bio18', 'bio19')


myBiomodProjIPSL45 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl45,
  proj.name = '2070_SDM_4_5_IPSL_wmax',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

myExplIPSL85 = stack( "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_IPSL_CM6A_LR_ssp585_2061_2080_2.asc",
                      "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_IPSL_CM6A_LR_ssp585_2061_2080_7.asc",
                      "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_IPSL_CM6A_LR_ssp585_2061_2080_11.asc",
                      "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_IPSL_CM6A_LR_ssp585_2061_2080_12.asc",
                      "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_IPSL_CM6A_LR_ssp585_2061_2080_17.asc",
                      "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_IPSL_CM6A_LR_ssp585_2061_2080_18.asc",
                      "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_IPSL_CM6A_LR_ssp585_2061_2080_19.asc")

names(myExplIPSL85) <- c('bio2', 'bio7', 'bio11', 'bio12', 'bio17', 'bio18', 'bio19')


myBiomodProjIPSL85 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExplIPSL85,
  proj.name = '2070_SDM_8_5_IPSL_wmax',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

#FUTURE FORECASTING
myBiomodEF_IPSL45 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjIPSL45
)

myBiomodEF_IPSL45



pdf("2070_4_5_IPSL.pdf")
plot(myBiomodEF_IPSL45)
dev.off()

myBiomodEF_IPSL85 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjIPSL85
)

myBiomodEF_IPSL85



pdf("2070_8_5_IPSL.pdf")
plot(myBiomodEF_IPSL85)
dev.off()

####ES
myExpl45 = stack( "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC_ES2L_ssp245_2061_2080_2.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC_ES2L_ssp245_2061_2080_7.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC_ES2L_ssp245_2061_2080_11.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC_ES2L_ssp245_2061_2080_12.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC_ES2L_ssp245_2061_2080_17.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC_ES2L_ssp245_2061_2080_18.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC_ES2L_ssp245_2061_2080_19.asc")

names(myExpl45) <- c('bio2', 'bio7', 'bio11', 'bio12', 'bio17', 'bio18', 'bio19')


myBiomodProjES45 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl45,
  proj.name = '2070_SDM_4_5_ES_wmax',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

myExplES85 = stack( "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC_ES2L_ssp585_2061_2080_2.asc",
                    "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC_ES2L_ssp585_2061_2080_7.asc",
                    "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC_ES2L_ssp585_2061_2080_11.asc",
                    "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC_ES2L_ssp585_2061_2080_12.asc",
                    "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC_ES2L_ssp585_2061_2080_17.asc",
                    "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC_ES2L_ssp585_2061_2080_18.asc",
                    "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC_ES2L_ssp585_2061_2080_19.asc")

names(myExplES85) <- c('bio2', 'bio7', 'bio11', 'bio12', 'bio17', 'bio18', 'bio19')


myBiomodProjES85 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExplES85,
  proj.name = '2070_SDM_8_5_ES_wmax',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

#FUTURE FORECASTING
myBiomodEF_ES45 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjES45
)

myBiomodEF_ES45



pdf("2070_4_5_ES.pdf")
plot(myBiomodEF_ES45)
dev.off()

myBiomodEF_ES85 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjES85
)

myBiomodEF_ES85



pdf("2070_8_5_ES.pdf")
plot(myBiomodEF_ES85)
dev.off()


##MIROC6 (MIR)
myExpl45 = stack( "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC6_ssp245_2061_2080_2.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC6_ssp245_2061_2080_7.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC6_ssp245_2061_2080_11.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC6_ssp245_2061_2080_12.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC6_ssp245_2061_2080_17.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC6_ssp245_2061_2080_18.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC6_ssp245_2061_2080_19.asc")

names(myExpl45) <- c('bio2', 'bio7', 'bio11', 'bio12', 'bio17', 'bio18', 'bio19')


myBiomodProjMIR45 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl45,
  proj.name = '2070_SDM_4_5_MIR_wmax',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

myExplMIR85 = stack( "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC6_ssp585_2061_2080_2.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC6_ssp585_2061_2080_7.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC6_ssp585_2061_2080_11.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC6_ssp585_2061_2080_12.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC6_ssp585_2061_2080_17.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC6_ssp585_2061_2080_18.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MIROC6_ssp585_2061_2080_19.asc")

names(myExplMIR85) <- c('bio2', 'bio7', 'bio11', 'bio12', 'bio17', 'bio18', 'bio19')


myBiomodProjMIR85 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExplMIR85,
  proj.name = '2070_SDM_8_5_MIR_wmax',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

#FUTURE FORECASTING
myBiomodEF_MIR45 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjMIR45
)

myBiomodEF_MIR45



pdf("2070_4_5_MIR.pdf")
plot(myBiomodEF_MIR45)
dev.off()

myBiomodEF_MIR85 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjMIR85
)

myBiomodEF_MIR85



pdf("2070_8_5_MIR.pdf")
plot(myBiomodEF_MIR85)
dev.off()

###MRI
myExpl45 = stack( "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MRI_ESM2_0_ssp245_2061_2080_2.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MRI_ESM2_0_ssp245_2061_2080_7.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MRI_ESM2_0_ssp245_2061_2080_11.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MRI_ESM2_0_ssp245_2061_2080_12.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MRI_ESM2_0_ssp245_2061_2080_17.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MRI_ESM2_0_ssp245_2061_2080_18.asc",
                  "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MRI_ESM2_0_ssp245_2061_2080_19.asc")

names(myExpl45) <- c('bio2', 'bio7', 'bio11', 'bio12', 'bio17', 'bio18', 'bio19')


myBiomodProjMRI45 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl45,
  proj.name = '2070_SDM_4_5_MRI_wmax',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

myExplMRI85 = stack( "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MRI_ESM2_0_ssp585_2061_2080_2.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MRI_ESM2_0_ssp585_2061_2080_7.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MRI_ESM2_0_ssp585_2061_2080_11.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MRI_ESM2_0_ssp585_2061_2080_12.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MRI_ESM2_0_ssp585_2061_2080_17.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MRI_ESM2_0_ssp585_2061_2080_18.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_MRI_ESM2_0_ssp585_2061_2080_19.asc")

names(myExplMRI85) <- c('bio2', 'bio7', 'bio11', 'bio12', 'bio17', 'bio18', 'bio19')


myBiomodProjMRI85 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExplMRI85,
  proj.name = '2070_SDM_8_5_MRI_wmax',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')

#FUTURE FORECASTING
myBiomodEF_MRI45 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjMRI45
)

myBiomodEF_MRI45



pdf("2070_4_5_MRI.pdf")
plot(myBiomodEF_MRI45)
dev.off()

myBiomodEF_MRI85 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjMRI85
)

myBiomodEF_MRI85



pdf("2070_8_5_MRI.pdf")
plot(myBiomodEF_MRI85)
dev.off()





sink()
