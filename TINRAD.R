###################################################################################################################################################################################
# Unstructured-Mesh Terrain Analysis and Incident Solar Radiation for Continuos Hydrologic Modeling in Mountain Watersheds
#Please cite:
# Moreno, H.A., Ogden, F.L, Alvarez, L.V., 201x. Unstructured-Mesh Terrain Analysis and Incident Solar Radiation for Continuous Hydrologic Modeling in Mountain Watershed
# Computers and Geoscience, Issue, Vol, pages...
###################################################################################################################################################################################
################# LOAD LIBRARIES ##################################################################################################################################################
rm(list=ls())  # remove all objects before starting
library(plotrix)  # color legend if makeplots=1
library(insol)  # tools to compute declination and hour angle
options(digits=16)
################ INPUT MODULE #####################################################################################################################################################
WD="/media/D/Work/Software/Scripts/OU/SOLAR_RADIATION_COMPLETE/Data/" # Working directory
TINF<<-"TIN_Centers.txt"  # txt file with the points file x,y,z coordinates of the TIN center points from the TINGEOFIX script
Vert<<-"Vertices2.txt" ## txt file with the second list of triangle vertices from the TINGEOFIX script
clipTIN<<-"clipTIN.rda" ## Delaunay file produced with TINGEOFIX script
shading<<-"Shading.rda" ## Shading file produced with TINGEOFIX script
ordersky<<-"ordersky.rda" # Sky view fractions produced with TINGEOFIX script
nuv<<-"nu.txt"  # Normal unit vector from TINGEOFIX script
slices=16  # how many slices should use? (e.g. 4,8,16,32,64,128)..Must coincide with slices used in TINGEOFIX
makeplots=0 # 1 if want to see some plots; O if not
# Dates and sun parameters
iyear=2000    # initial year GMT 
imonth=3     # initial month in GMT
iday=21       # Initial day on GMT
ihour=0       # Initial hour in GMT
imin=0        # Initial minute in GMT
isec=0        # Initial second in GMT
fyear=2000    # final year GMT 
fmonth=3      # final month in GMT
fday=22       # final day on GMT
fhour=0       # final hour in GMT
fmin=0        # final minute in GMT
fsec=0        # final second in GMT
sconstant=1361  # solar constant in w/m2
## coordinate transformation
RE=6378137 		# in meters; Earth's radius at the equator GRS80 ellipsoid
Lambda00<<--109 	# reference meridian in degrees. Sign indicates if it's west or east of Greenwich
FalseEast<<-20000000 	# False East of desired coordinate system
FalseNorth<<-10000000	# False North of  desired coordinate system
#### Atmospheric parameters at GMT times coincident with the sequence of SUN parameters above ############################################################################################################
Ta=c(277,272,270,267,268,269,260,259,259,258,258,255,256,258,257,259,263,267,271,273,276,277,278,278,277)  # Hourly air temperatures in kELVIN.. this can be a vector with values at each triangle center of TINCENTERS
RH=c(24,32,36,50,52,48,49,49,50,50,50,48,54,52,52,50,51,44,39,41,33,25,23,23,24) # Hourly basin-average air relative humidities in %.. this can be a vector with values at each triangle center of TINCENTERS
albedo=0.65   # Abbedo map values. It can be one basin-average value for all the simulation or one for each element and time step 
################INPUT MODULE ENDS ########################################################################################################################################################################

#### MODULE TO COMPUTE PRECIPITABLE WATER CONTENT in cm from the Relative Humidity(%) and Air Temperature(K)
w=0.00493*(RH/Ta)*exp(26.23-(5416/Ta))  # Precipitable water in cm; when RH and Ta are matrices per hour
############### MODULE METADATA AND HEADERS #######################################################################################################################################
cat("Time started processing ",date(),fill=TRUE)
TINET=read.table(paste(WD,TINF,sep=""))
centers=TINET[,1:2]
elecenters=TINET[,3]
TINCENTERS=cbind(centers,elecenters)
clipTINO<<-load(file=paste(WD,clipTIN,sep="")) 
nu<<-read.table(paste(WD,nuv,sep=""))
load(file=paste(WD,shading,sep="")) 
load(file=paste(WD,ordersky,sep="")) 
#### UNIT SOLAR VECTOR COMPUTATION ###########################################################################################################################
## lATITUDE AND LONGITUDE computation from the centers vector
clonlat=matrix(NA,dim(centers)[1],dim(centers)[2])
colnames(clonlat)=c("Lon","Lat") 
Lambda0=Lambda00*pi/180
clonlat[,2]=(180/pi)*((centers[,2]-FalseNorth)/RE)
clonlat[,1]=(((centers[,1]-FalseEast)/(RE*cos(clonlat[,2]*pi/180)))+Lambda0)*(180/pi)
timezone=clonlat[,1]/15 # correction factor by longitude position to compute Julian date
cat("Maximum local aparent time difference between extreme longitude elements of the basin is ", 3600*(max(timezone)-min(timezone)), " seconds", fill=TRUE)

# Creates sequence of hours for GMT times
secuencia=as.character(seq(ISOdate(iyear,imonth,iday,ihour, imin, isec), ISOdate(fyear,fmonth,fday, fhour, fmin, fsec),"hour"))
maxdi=timezone[which(abs(timezone)==max(abs(timezone)))]
cat("Maximum time difference between GMT and local time to generate new sequence is ",maxdi," hours", fill=TRUE)
sunarr=matrix(NA,length(secuencia),3)  ## matrix containing the basin mean sun position vector i,j,k for each day of secuencia only for sun animation purposes
diradarr=matrix(NA, dim(clonlat)[1],length(secuencia))  # matrix containing the short wave radiatin for every element in the rows and every time in the columns
direm=matrix(NA, dim(clonlat)[1],length(secuencia))  # matrix containing the shadowed elements per every time step
Kdif=matrix(NA, dim(clonlat)[1],length(secuencia))  # matrix containing the diffuse radiation per every time step
Kbs=matrix(NA, dim(clonlat)[1],length(secuencia))  # matrix containing the backscattered radiation per every time step
Kal=matrix(NA,dim(clonlat)[1],length(secuencia))  # matrix containing the backscattered radiation from surrounding landscapes
Twa=matrix(NA, dim(clonlat)[1],length(secuencia)) # 
Tda=matrix(NA, dim(clonlat)[1],length(secuencia))
Tws=matrix(NA, dim(clonlat)[1],length(secuencia))
Trs=matrix(NA, dim(clonlat)[1],length(secuencia))
Tds=matrix(NA, dim(clonlat)[1],length(secuencia))
finalrad = matrix(NA,nrow(diradarr),ncol(diradarr))
marr=matrix(NA, dim(clonlat)[1],length(secuencia)) # optical air mass values per time step averaged over the basin

for (ti in 1:length(secuencia)){  # for loop for every time step we want to compute the sunvector
  cat("GMT time ", secuencia[ti],fill=TRUE)
  pa= substring(secuencia[ti] , 1:19, 1:19)
  yeGMT=as.numeric(paste(pa[1],pa[2],pa[3],pa[4],sep=""))
  moGMT=as.numeric(paste(pa[6],pa[7],sep=""))
  daGMT=as.numeric(paste(pa[9],pa[10],sep=""))
  hoGMT=as.numeric(paste(pa[12],pa[13],sep=""))
  miGMT=as.numeric(paste(pa[15],pa[16],sep=""))
  seGMT=as.numeric(paste(pa[18],pa[19],sep=""))
  JudayGMT=JDymd(yeGMT, moGMT, daGMT, hoGMT, miGMT, seGMT)
  localJuday=JudayGMT+(timezone/24) # Vector of julian days for each element of the TIN
  delta=declination(localJuday)*pi/180 # vector of declination of sun in radians.. south negative... north positive  
  ha=hourangle(localJuday,clonlat[,1],timezone) # return hour angle in radians and Juday is local time.. ha should be positive before noon (E) and negative afternoon (W), although this command provides the contrary
  sundis=sunr(localJuday) 
  ## Here are two ways to compute sunvector  SVX == sunvec[,1], SVY= sunvec[,2], SVZ= sunvec[,3] these should conduct to the same computations
  SVX=-sin(ha)*cos(delta)
  SVY=(sin(clonlat[,2]*pi/180)*cos(ha)*cos(delta))-(cos(clonlat[,2]*pi/180)*sin(delta))
  SVZ=(cos(clonlat[,2]*pi/180)*cos(ha)*cos(delta))+(sin(clonlat[,2]*pi/180)*sin(delta))
  
  sunvec=sunvector(localJuday,clonlat[,2],clonlat[,1],timezone)
  sunvec[,2]=-sunvec[,2] # change of sign for the south-north component to make it poisitive looking north
  sunarr[ti,]=colMeans(sunvec)  # sunvec[1,]  # takes the mean of sunvector
  zenithangle=(atan((((sunvec[,1]^2)+(sunvec[,2]^2))^0.5)/(sunvec[,3])))    # average spatial zenith angle in rads. Negative values are below horizon
  zenithangle[which(zenithangle<=0)]=NA
  zenithangle_degree=zenithangle*180/pi
  cosza=cos(zenithangle)  # cosine of zenith angle
  m=1/(cosza+(0.50572*(96.07995-(zenithangle*180/pi))^-1.6364))*exp(-TINCENTERS[,3]/7000) # Kasten and Young (1989)
  
  marr[,ti]=m
  ### Atmospheric Transmissivities
  # Transmissivity of water vapor
  Twa[,ti]=1-(0.077*(m*w[ti])^0.3)   # Transmissivity of water vapor 
  Tda[,ti]=0.965^m               # Transmissivity of dust and aerosols 
  Tws[,ti]=1-(0.0225*m*w[ti])        # Scattered by water vapor 
  Trs[,ti]=0.972-(0.08262*m)+(0.00933*(m^2))-(0.00095*(m^3))+(0.0000437*(m^4)) # Raileigh Scattering
  Tds[,ti]=Tda[,ti]                   # Scattering by dust and aerosols 
  
  ### calculation of cosine influence on every triangle
  nuv=as.matrix(nu[,1:3])
  costheta=rep(NA,dim(sunvec)[1])
  for (dp in 1:dim(sunvec)[1]){
  costheta[dp]=sunvec[dp,]%*%nuv[dp,]  # computes the dot product
  }
  
  dirad=sconstant*sundis*costheta*Twa[,ti]*Tda[,ti]*(Tws[,ti]*Trs[,ti]*Tds[,ti])   # solar constant * sundistance * costheta
  dirad[which(dirad<=0)]=0          # pixel is shading itself as cos(theta)<0, so angle between two vectors is pi/2 or greater
 
  belowhorizon=which(sunvec[,3]<0) # identify all triangles for which sun is already below the horizon (k<0)
  dirad[belowhorizon]=0
  diradarr[,ti]=dirad
  
  azsun=(((sunvec[,1]/abs(sunvec[,1]))*(pi/2))-atan(sunvec[,2]/sunvec[,1]))*180/pi  # gives angles between -180 and 180
  azsun[which(azsun<0)]=azsun[which(azsun<0)]+360   # gives angles between 0 and 360
  slicesun=trunc(slices*azsun/360)+1 #determines what slices is the sun shining to
  REM=rep(1,dim(sunvec)[1])
  
  for (cui in 1:dim(sunvec)[1]){
    if (is.na(zenithangle_degree[cui])==FALSE){
    if (Shading[cui,slicesun[cui]]<=zenithangle_degree[cui]){
    REM[cui]=1  # shaded
    } else {
    REM[cui]=0  # not shaded 
    }
  }
  }
   
  direm[,ti]=REM
  direm[belowhorizon,ti]=0
 
  Kdif[,ti]=0.05*sconstant*sundis*Twa[,ti]*Tda[,ti]*(Tws[,ti]*Trs[,ti]*Tds[,ti])*ordersky[,4]   
  Kdif[belowhorizon,ti]=0
  Kbs[,ti]=albedo*(diradarr[,ti]+Kdif[,ti])*(0.05*Twa[,ti]*Tda[,ti]*(Tws[,ti]*Trs[,ti]*Tds[,ti]))*ordersky[,4]  
  Kbs[belowhorizon,ti]=0
  Kal[,ti]=albedo*(1-ordersky[,4])*mean(diradarr[,ti])
  finalrad[,ti]=(diradarr[,ti]*(1-direm[,ti])) + Kdif[,ti]  + Kbs[,ti]  + Kal[,ti]
  }
if (makeplots==1){
########################################################################################################
# Snapshots of final_rad ##
for (ri in 1:length(secuencia)){
  cat("radiation plotting ",ri," of ", length(secuencia), fill=TRUE)
  pdf(file=paste(WD,"Final_Rad_GMT",secuencia[ri],".pdf",sep=""))
  testcol<-grey(seq(0,1, len=128))   # grey pallete
  #testcol = ygobb(256)  # topo.colors, terrain.colors, cm.colors, heat.colors, rainbow, ygobb, matlab.like2, primary.colors, diverge_hsv, cyan2yellow, green2red, blue2red, magenta2green 
  col.labels<-c(round(min(finalrad), digits=2),round((((max(finalrad)-min(finalrad)))/4)+min(finalrad),digits=2),
                round((((max(finalrad)-min(finalrad)))/2)+min(finalrad),digits=2),
                round((((max(finalrad)-min(finalrad)))*0.75)+min(finalrad), digits=2),round(max(finalrad), digits=2))
  plot(clipTIN, lwd=0.05,col=grey((finalrad[,ri]-min(finalrad))/(max(finalrad)-min(finalrad))))
  #plot(clipTIN, lwd=0.05,col=testcol[c(round(255*(diradarr[,ri]-min(diradarr[,ri]))/(max(diradarr[,ri])-min(diradarr[,ri]))))+1])
  color.legend(max(centers[,1]+1000),min(centers[,2]),max(centers[,1])+2000,max(centers[,2]), col.labels, testcol, gradient="y", cex=0.6, align="rb")
  dev.off()
}
}



