# load library
library(biomod2)
library(raster)
library(dplyr)
library(ggplot2)
library(usdm)

#Ensemble Model
(myBiomodEM <-load("/projects/b1045/yang/Artocarpus.altilis/Artocarpus.altilis.1611612758ensemble.models.out"))
myBiomodEM <- get(myBiomodEM)


(myBiomodModelOut <-load("/projects/b1045/yang/Artocarpus.altilis/Artocarpus.altilis.1611612758.models.out"))
myBiomodModelOut <- get(myBiomodModelOut)



myExplBCC45 = stack( "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp245_2061_2080_2.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp245_2061_2080_7.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp245_2061_2080_11.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp245_2061_2080_12.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp245_2061_2080_17.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp245_2061_2080_18.asc",
                     "/projects/b1045/yang/wc21/10m/wc2.1_10m_bioc_BCC_CSM2_MR_ssp245_2061_2080_19.asc")

names(myExplBCC45) <- c('bio2', 'bio7', 'bio11', 'bio12', 'bio17', 'bio18', 'bio19')

e <- extent(-180, 180, -45, 45)
myExplBCC45 <- crop(myExplBCC45, e)
myExplBCC45 <- stack(myExplBCC45)



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

myExplBCC85 <- crop(myExplBCC85, e)
myExplBCC85 <- stack(myExplBCC85)

myBiomodProjBCC85 <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExplBCC85,
  proj.name = '2070_SDM_8_5_BCC_wmax',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = 'xz',
  clamping.mask = T,
  output.format = '.grd')



print(4)

#FUTURE FORECASTING
myBiomodEF_BCC45 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjBCC45,
  proj.name = 'BCC45'
)



myBiomodEF_BCC45



pdf("2070_4_5_BCC.pdf")
plot(myBiomodEF_BCC45)
dev.off()

myBiomodEF_BCC85 <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjBCC85
)

myBiomodEF_BCC85

#writeRaster(myBiomodEF, filename="future_SDM.grd", overwrite=TRUE, format="raster")
 
 
 pdf("2070_8_5_BCC.pdf")
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

myExpl45 <- crop(myExpl45, e)
myExpl45 <- stack(myExpl45)


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


myExplCan85 <- crop(myExplCan85, e)
myExplCan85 <- stack(myExplCan85)

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

myExpl45 <- crop(myExpl45, e)
myExpl45 <- stack(myExpl45)

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

myExplCM85 <- crop(myExplCM85, e)
myExplCM85 <- stack(myExplCM85)


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

myExpl45 <- crop(myExpl45, e)
myExpl45 <- stack(myExpl45)

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

myExplESM85 <- crop(myExplESM85, e)
myExplESM85 <- stack(myExplESM85)



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

myExpl45 <- crop(myExpl45, e)
myExpl45 <- stack(myExpl45)


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

myExplIPSL85 <- crop(myExplIPSL85, e)
myExplIPSL85 <- stack(myExplIPSL85)



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

#writeRaster(myBiomodEF, filename="future_SDM.grd", overwrite=TRUE, format="raster")


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

myExpl45 <- crop(myExpl45, e)
myExpl45 <- stack(myExpl45)


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

myExplES85 <- crop(myExplES85, e)
myExplES85 <- stack(myExplES85)



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

myExpl45 <- crop(myExpl45, e)
myExpl45 <- stack(myExpl45)

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

myExplMIR85 <- crop(myExplMIR85, e)
myExplMIR85 <- stack(myExplMIR85)



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

myExpl45 <- crop(myExpl45, e)
myExpl45 <- stack(myExpl45)

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

myExplMRI85 <- crop(myExplMRI85, e)
myExplMRI85 <- stack(myExplMRI85)

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
