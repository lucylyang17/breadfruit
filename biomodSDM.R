library("biomod2")

#Artocarpus altilis presence data
artocarpus <- read.csv("/home/lly2413/breadfruit/dat_full.csv")

# the name
myRespName <- 'Artocarpus altilis'
# the XY coordinates of the presence
myRespXY <- artocarpus[,c(1,2)]
# and the presence data
myResp <-  as.numeric(artocarpus[,3])


myExpl = stack( system.file( "external/bioclim/current/bio3.grd", 
                             package="biomod2"),
                system.file( "external/bioclim/current/bio4.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio7.grd", 
                             package="biomod2"),  
                system.file( "external/bioclim/current/bio11.grd", 
                             package="biomod2"), 
                system.file( "external/bioclim/current/bio12.grd", 
                             package="biomod2"))

myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep = 3,
                                     PA.nb.absences = 200,
                                     PA.strategy = 'random')


myBiomodData
plot(myBiomodData)

