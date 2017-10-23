###################################################################################################################################################################################
# Code to produce a delaunay triangulation from a set of triangle vertices (vertices.txt)
# Please cite:
# Moreno, H.A., Ogden, F.L, Alvarez, L.V., 201x. Unstructured-Mesh Terrain Analysis and Incident Solar Radiation for Continuous Hydrologic Modeling in Mountain Watershed
# Computers and Geoscience, Issue, Vol, pages...
################# LOAD LIBRARIES ##################################################################################################################################################
rm(list=ls())  # remove all objects before starting
library(rgeos)
library(sp)
library(maptools)
library(spatstat)
library(fields)
library(colorRamps)
library(plotrix)
library(moments)
options(digits=16)

################ INPUT MODULE #####################################################################################################################################################
WD         = "./"           # set working directory
divisoria  = "catchments"   # Watershed divide shapefile to clip triangulation
vertices   = "vertices.txt" # textfile with the TIN triangle vertices X,Y,Z (elevation).. Coordinates in projected system in meters
maxres     = 60             # in length units, largest expected triangle edge within the TIN to find nearest three nodes. Rule of thumb: If TIN was derived from a 30m DEM then use 30m to begin
tolerancia = 0              # tolerance value for distance when simplifying the tesseltaion achieved with Delaunay.. begin with tolerance of 0, if processing is heavy then increase it
slices     = 16             # if calculateSKY=1, how many slices should use? (e.g. 4,8,16,32,64,128)..The more the more accurate the calculations
RE         = 6378137 		    # in meters; Earth's radius at the equator GRS80 ellipsoid
makeplots  = 0              # 0 for no making plots; 1 for making plots..May save time not making the plots.

###############Script processsing begins #######################################################################################################################################
cat("Time started processing ",date(),fill=TRUE)
TIN=read.table(paste(WD,vertices,sep=""))

############# Checking for duplicated nodes and removing them from the list of triangle nodes
DUPT     = duplicated(TIN)
dondedup = which(DUPT==TRUE)

if (length(dondedup)>0){
  cat("Eliminating ",length(dondedup)," points",fill=TRUE)
  TINET=TIN[-dondedup,]
} else {
  TINET=TIN
}

colnames(TINET)=c("x","y","z")

if (makeplots==1){
  pdf(file=paste(WD,"1_Basin_domain_Vertices.pdf",sep=""))
  plot(TINET[,1],TINET[,2],pch=19,cex=0.1,main="Basin domain", xlim=c(min(TINET[,1]),max(TINET[,1])),ylim=c(min(TINET[,2]),max(TINET[,2])))
  dev.off()
}

################################## Creating triangle tesselation ########################################################################
W <- ripras(TINET[,1:2], shape="convex") 
#W <- owin(c(0, 30), c(0, 30)) 
## create point pattern 
X <- as.ppp(TINET[,1:2], W=W) 

if (makeplots==1){
  pdf(file=paste(WD,"2_Point_pattern.pdf",sep=""))
  plot(X, pch=20, cex=0.2)
  dev.off()
}

## generate tiles
if (makeplots==1){
  Y <- delaunay(X) 
  pdf(file=paste(WD,"3_Delaunay.pdf",sep=""))
  plot(Y,lwd=0.2)
  dev.off()
}

## Read shape file with watershed divide to clip final map using maptools
divi=readShapePoly(paste(WD,divisoria,sep=""))  # using maptools library

## convert Y (Delaunay) to spatial polygons
Z <- as(Y, "SpatialPolygons") 
if (makeplots==1){
  pdf(file=paste(WD,"4_Delaunay_SpatialPolygons.pdf",sep=""))
  plot(Z, col=gray(TINET[,3]/max(TINET[,3])), lwd=0.1) ## polygons are not plotted in the original order
  dev.off()
}

### Clipping and plotting ####################################################
if (tolerancia==0){  # Use this if clipping and plotting takes less than an hour
  clipTIN=gIntersection(Z,divi, byid=TRUE, drop_lower_td=TRUE)   # using library sp and geos
  if (makeplots==1){
    pdf(file=paste(WD,"5_Delaunay_Clipped.pdf",sep=""))
    plot(clipTIN) #, col=grey(TINET[,3]/max(TINET[,3])), lwd=0.1)
    dev.off()
  }
} else {  # Simplify TIN segments by simplifing geometry using a tolerance value tolerancia
  SimplyTIN=gSimplify(Z, tol=tolerancia) # Function simplifies the given geometry using the Douglas-Peuker algorithm
  # The purpose of the algorithm is, given a curve composed of line segments, to find a similar curve with fewer points. 
  # The algorithm defines 'dissimilar' based on the maximum distance between the original curve and the simplified curve. 
  # The simplified curve consists of a subset of the points that defined the original curve.
  if (makeplots==1){
    pdf(file=paste(WD,"6_Delaunay_SpatialPolygons_Simplified_tol=",tolerancia,".pdf",sep=""))
    plot(SimplyTIN) #, col=grey(TINET[,3]/max(TINET[,3])), lwd=0.1) ## polygons are not plotted in the original order
    dev.off()
  }
  clipTIN=gIntersection(SimplyTIN,divi, byid=TRUE, drop_lower_td=TRUE)

  if (makeplots==1){
    pdf(file=paste(WD,"7_Delaunay_Clipped_Simplified_tol=",tolerancia,".pdf",sep=""))
    plot(clipTIN) #, col=grey(TINET[,3]/max(TINET[,3])), lwd=0.1)
    dev.off()
  }
}

centers    = coordinates(clipTIN)
elecenters = matrix(NA,dim(centers)[1])
vertex     = matrix(NA,dim(centers)[1],9) 

#### Computing the element centers and assigning one elevation to each Tile based on the average of the three closest triangle nodes ################
for (mi in 1:10){
  maxresi = maxres*mi  
  cat("Maxresi is ",maxres*(mi), fill=TRUE)
  dista  = fields.rdist.near(centers,TINET[,1:2], delta=maxresi)  # the parameter mean.neighbor is set to 100 by default
  dista1 = as.data.frame(dista[1])
  dista2 = as.data.frame(dista[2])
  for (ei in 1:dim(centers)[1]){
    cat("Finding mean elevation of node ", ei, " of ", dim(centers)[1], fill=TRUE)
    dondecenter=which(dista1[,1]==ei)  
    if (length(dondecenter)<3){
      cat("Increasing maxresi to ",maxres*(mi+1), fill=TRUE)
      if (mi==20) cat("Please increase maxresi to a higher value 2X (e.g double)",fill=TRUE)
      break
    } else {
      Onei           = which(dista2[dondecenter,1] == sort(dista2[dondecenter,1])[1])  
      Twoi           = which(dista2[dondecenter,1] == sort(dista2[dondecenter,1])[2])  
      Threei         = which(dista2[dondecenter,1] == sort(dista2[dondecenter,1])[3]) 
      Indexnear      = c(dista1[dondecenter[Onei],2], dista1[dondecenter[Twoi],2], dista1[dondecenter[Threei],2]) 
      elecenters[ei] = mean(TINET[Indexnear,3]) 
      vertex[ei,1:3] = TINET[Indexnear,1]
      vertex[ei,4:6] = TINET[Indexnear,2]
      vertex[ei,7:9] = TINET[Indexnear,3]
    } 
  }
  if (ei==dim(centers)[1])
    break
}

if (makeplots==1){
  pdf(file=paste(WD,"8_Delaunay_Elevations_matlab_like_colors",".pdf",sep=""))
  testcol = matlab.like2(256)  # topo.colors, terrain.colors, cm.colors, heat.colors, rainbow, ygobb, matlab.like2, primary.colors, , diverge_hsv, cyan2yellow, green2red, blue2red, magenta2green 
  col.labels<-c(round(min(elecenters)),round((((max(elecenters)-min(elecenters)))/4)+min(elecenters)),round((((max(elecenters)-min(elecenters)))/2)+min(elecenters)),
              round((((max(elecenters)-min(elecenters)))*0.75)+min(elecenters)),round(max(elecenters)))
  #plot(clipTIN, lwd=0.05,col=grey((elecenters-min(elecenters))/(max(elecenters)-min(elecenters))))
  plot(clipTIN, lwd=0.05,col=testcol[c(round(255*(elecenters-min(elecenters))/(max(elecenters)-min(elecenters))))+1])
  color.legend(max(centers[,1]+1000),min(centers[,2]),max(centers[,1])+2000,max(centers[,2]), col.labels, testcol, gradient="y", cex=0.6, align="rb")
  dev.off()
}

## Exporting the TIN element centers and elevations
TINCENTERS=cbind(centers, elecenters)
write.table(TINCENTERS,paste(WD,"TIN_Centers.txt",sep=""), quote=F, row.names=F,col.names=F)

## expoirting clipTIN for future access and plotting
save(clipTIN, file=paste(WD,"clipTIN.rda",sep=""))

## Exporting organized Vertices in txt
write.table(vertex,paste(WD,"vertices2.txt",sep=""), quote=F, row.names=F,col.names=F)

########################################## MODULE TO COMPUTE THE UNIT NORMAL VECTOR ##############################################################################
nv      = matrix(NA,dim(centers)[1],4)   # normal vector whose columns determine the components i,j,k, magnitude
#vertex = read.table(paste(WD,Vert,sep=""))
X1      = vertex[,1]
X2      = vertex[,2]
X3      = vertex[,3]
Y1      = vertex[,4]
Y2      = vertex[,5]
Y3      = vertex[,6]
Z1      = vertex[,7]
Z2      = vertex[,8]
Z3      = vertex[,9]
nv[,1]  = Z3*(Y2-Y1)+Z2*(Y1-Y3)+Z1*(Y3-Y2)  # i component
nv[,2]  = Z3*(X1-X2)+Z2*(X3-X1)+Z1*(X2-X3)  # j component
nv[,3]  = Y3*(X2-X1)+Y2*(X1-X3)+Y1*(X3-X2)  # k component
nv[,4]  = ((nv[,1]^2)+(nv[,2]^2)+(nv[,3]^2))^0.5
nu      = nv
difo    = which(nv[,4] != 0)

if (length(difo)!=dim(nu)[1]) {
  cat("careful..there are null triangles in the region",fill=TRUE)
  Sys.sleep(600)
}

nu[difo,]=nv[difo,]/nv[difo,4] # Normal unitary vector

## multiply any negative k component row by -1 to obtain the upwards direction vector
negnu=which(nu[,3]<0)
nu[negnu,]=nu[negnu,]*(-1)
write.table(nu,paste(WD,"nu.txt",sep=""),col.names=FALSE, row.names=FALSE)

##### Plotting the components of normal unit vector
### i component largest values (east-looking) are white color and smallest values are dark
if (makeplots==1){
  pdf(file=paste(WD,"9_X-coordinate_of_Normal_Unit_Vector.pdf",sep=""))
  #testcol<-grey(seq(0,1, len=128))   # grey pallete
  testcol = matlab.like2(256)  # topo.colors, terrain.colors, cm.colors, heat.colors, rainbow, ygobb, matlab.like2, primary.colors, diverge_hsv, cyan2yellow, green2red, blue2red, magenta2green 
  col.labels<-c(round(min(nu[,1]), digits=2),round((((max(nu[,1])-min(nu[,1])))/4)+min(nu[,1]),digits=2),round((((max(nu[,1])-min(nu[,1])))/2)+min(nu[,1]),digits=2),
              round((((max(nu[,1])-min(nu[,1])))*0.75)+min(nu[,1]), digits=2),round(max(nu[,1]), digits=2))
  #plot(clipTIN, lwd=0.05,col=grey((nu[,1]-min(nu[,1]))/(max(nu[,1])-min(nu[,1]))))
  plot(clipTIN, lwd=0.05,col=testcol[c(round(255*(nu[,1]-min(nu[,1]))/(max(nu[,1])-min(nu[,1]))))+1])
  color.legend(max(centers[,1]+1000),min(centers[,2]),max(centers[,1])+2000,max(centers[,2]), col.labels, testcol, gradient="y", cex=0.6, align="rb")
  dev.off()

  pdf(file=paste(WD,"10_X-coordinate_Normal_Unit_Vector_density.pdf",sep=""))
  hist(nu[,1], probability=T, ylim=c(0,5))
  lines(density(nu[,1]), col="red", ylim=c(0,5))
  dev.off()
}

cat("Moments of component x (i) east-west component", fill=TRUE)
all.moments(nu[,1], order.max=4)
summary(nu[,1])

### j component
if (makeplots==1){
  pdf(file=paste(WD,"11_Y-coordinate_of_Normal_Unit_Vector.pdf",sep=""))
  #testcol<-grey(seq(0,1, len=128))   # grey pallete
  testcol = matlab.like2(256)  # topo.colors, terrain.colors, cm.colors, heat.colors, rainbow, ygobb, matlab.like2, primary.colors, , diverge_hsv, cyan2yellow, green2red, blue2red, magenta2green 
  col.labels<-c(round(min(nu[,2]), digits=2),round((((max(nu[,2])-min(nu[,2])))/4)+min(nu[,2]),digits=2),round((((max(nu[,2])-min(nu[,2])))/2)+min(nu[,2]),digits=2),
              round((((max(nu[,2])-min(nu[,2])))*0.75)+min(nu[,2]), digits=2),round(max(nu[,2]), digits=2))
  #plot(clipTIN, lwd=0.05,col=grey((nu[,2]-min(nu[,2]))/(max(nu[,2])-min(nu[,2]))))
  plot(clipTIN, lwd=0.05,col=testcol[c(round(255*(nu[,2]-min(nu[,2]))/(max(nu[,2])-min(nu[,2]))))+1])
  color.legend(max(centers[,1]+1000),min(centers[,2]),max(centers[,1])+2000,max(centers[,2]), col.labels, testcol, gradient="y", cex=0.6, align="rb")
  dev.off()

  pdf(file=paste(WD,"12_Y-coordinate_of_Normal_Unit_Vector_density.pdf",sep=""))
  hist(nu[,2], probability=T, ylim=c(0,6))
  lines(density(nu[,2]), col="red", ylim=c(0,6))
  dev.off()
}
cat("Moments of component y (j) east-west component", fill=TRUE)
all.moments(nu[,2], order.max=4)
summary(nu[,2])

### k component
### k component largest values (upward-looking) are white color and smallest (vertical-looking) values are dark
if (makeplots==1){
  pdf(file=paste(WD,"13_Z-coordinate_of_Normal_Unit_Vector.pdf",sep=""))
  #testcol<-grey(seq(0,1, len=128))   # grey pallete
  testcol = matlab.like2(256)  # topo.colors, terrain.colors, cm.colors, heat.colors, rainbow, ygobb, matlab.like2, primary.colors, , diverge_hsv, cyan2yellow, green2red, blue2red, magenta2green 
  col.labels<-c(round(min(nu[,3]), digits=2),round((((max(nu[,3])-min(nu[,3])))/4)+min(nu[,3]),digits=2),round((((max(nu[,3])-min(nu[,3])))/2)+min(nu[,3]),digits=3),
              round((((max(nu[,3])-min(nu[,3])))*0.75)+min(nu[,3]), digits=2),round(max(nu[,3]), digits=2))
  #plot(clipTIN, lwd=0.05,col=grey((nu[,3]-min(nu[,3]))/(max(nu[,3])-min(nu[,3]))))
  plot(clipTIN, lwd=0.05,col=testcol[c(round(255*(nu[,3]-min(nu[,3]))/(max(nu[,3])-min(nu[,3]))))+1])
  color.legend(max(centers[,1]+1000),min(centers[,2]),max(centers[,1])+2000,max(centers[,2]), col.labels, testcol, gradient="y", cex=0.6, align="rb")
  dev.off()

  pdf(file=paste(WD,"14_Z-coordinate_of_Normal_Unit_Vector_density.pdf",sep=""))
  hist(nu[,3], probability=T, ylim=c(0,40))
  lines(density(nu[,3]), col="red", ylim=c(0,40))
  dev.off()
}
cat("Moments of component z (k upward component", fill=TRUE)
all.moments(nu[,3], order.max=4)
summary(nu[,3])

############################################################ Slope computation
## True slope
trueslope=180*acos(nu[,3])/pi
if (makeplots==1){
  pdf(file=paste(WD,"15_True_Slopes_matlablike.pdf",sep=""))
  #testcol<-grey(seq(0,1, len=128))   # grey pallete
  testcol = matlab.like(256)  # topo.colors, terrain.colors, cm.colors, heat.colors, rainbow, ygobb, matlab.like2, primary.colors, diverge_hsv, cyan2yellow, green2red, blue2red, magenta2green 
  col.labels<-c(round(min(trueslope)),round((((max(trueslope)-min(trueslope)))/4)+min(trueslope)),round((((max(trueslope)-min(trueslope)))/2)+min(trueslope)),
              round((((max(trueslope)-min(trueslope)))*0.75)+min(trueslope)),round(max(trueslope)))
  #plot(clipTIN, lwd=0.05,col=grey((trueslope-min(trueslope))/(max(trueslope)-min(trueslope))))
  plot(clipTIN, lwd=0.05,col=testcol[c(round(255*(trueslope-min(trueslope))/(max(trueslope)-min(trueslope))))+1])
  color.legend(max(centers[,1]+1000),min(centers[,2]),max(centers[,1])+2000,max(centers[,2]), col.labels, testcol, gradient="y", cex=0.6, align="rb")
  dev.off()
}
cat("Moments of component true slope", fill=TRUE)
all.moments(trueslope, order.max=4)
summary(trueslope)

############################################################ Aspect computation
novec=which(nu[,1]==0 & nu[,2]==0 & nu[,3]==0)        # Vector does not exist
cat("There is not normal vector possible at ",length(novec), " points",fill=TRUE)

flatt         = which(nu[,1]==0 & nu[,2]==0 & nu[,3] != 0)                    # trueAZ=NA
purenorth     = which(nu[,1] == 0 & nu[,2]>0)                 # trueAZ=0
purenorthwall = which(nu[,1] == 0 & nu[,2]>0 & nu[,3]==0) # trueAZ=0
pureeast      = which(nu[,2] == 0 & nu[,1]>0)                  # trueAZ=90
pureeastwall  = which(nu[,2] == 0 & nu[,1]>0 & nu[,3]==0)  # trueAZ=90
puresouth     = which(nu[,1] == 0 & nu[,2]<0)                 # trueAZ=180
puresouthwall = which(nu[,1] == 0 & nu[,2]<0 & nu[,3]==0) # trueAZ=180
purewest      = which(nu[,2] == 0 & nu[,1]<0)                  # trueAZ=-90
purewestwall  = which(nu[,2] == 0 & nu[,1]<0 & nu[,3]==0)  # trueAZ=-90

####### trueaspect
trueAZ                = (((nu[,1]/abs(nu[,1]))*(pi/2))-atan(nu[,2]/nu[,1]))*180/pi
trueAZ[purewestwall]  = -90
trueAZ[purewest]      = -90
trueAZ[puresouthwall] = 180
trueAZ[puresouth]     = 180
trueAZ[pureeastwall]  = 90
trueAZ[pureeast]      = 90
trueAZ[purenorthwall] = 0
trueAZ[purenorth]     = 0
trueAZ[flatt]         = NA
trueAZ[novec]         = NA

### reclassify values to 8 directins
# Direction    Values range     Class
# N          [0-22.5)           0
# NE         [22.5-67.5)        45
# E          [67.5-112.5)       90
# SE         [112.5-157.5)      135
# S          [157.5-180]        180
# S          [-180 - -157.5)    180
# SW         [-157.5 - -112.5) -135
# W          [-112.5 - -67.5)  -90
# NW         [-67.5 - -22.5)   -45
# N          [-22.5 - 0]         0

trueAZclass=trueAZ
trueAZclass[which(trueAZ>=0 & trueAZ<22.5)]=1 #0 # red
trueAZclass[which(trueAZ>=22.5 & trueAZ<67.5)]=2 #45 # orange
trueAZclass[which(trueAZ>=67.5 & trueAZ<112.5)]=3 #90 # yellow
trueAZclass[which(trueAZ>=112.5 & trueAZ<157.5)]=4 #135 # green
trueAZclass[which(trueAZ>=157.5 & trueAZ<=180)]=5 #180 # cyan
trueAZclass[which(trueAZ>=-180 & trueAZ < -157.5)]=5 # 180 # cyan
trueAZclass[which(trueAZ>=-157.5 & trueAZ < -112.5)]=6 #-135 # cornflowerblue
trueAZclass[which(trueAZ>=-112.5 & trueAZ < -67.5)]=7 #-90 # blue
trueAZclass[which(trueAZ>=-67.5 & trueAZ < -22.5)]=8 #-45 # pink
trueAZclass[which(trueAZ>=-22.5 & trueAZ <= 0)]=1 # red

if (makeplots==1){
  palette = c("red","orange","yellow","green","cyan","cornflowerblue","blue","pink")
  palette1=c("grey75","white","grey75","grey50","grey25","black","grey25","grey50")
  testcol = rainbow(8)  # topo.colors, terrain.colors, cm.colors, heat.colors, rainbow

  pdf(file=paste(WD,"20_True_AZ_color.pdf",sep=""))
  plot(clipTIN, lwd=0.05,col=palette[c(trueAZclass)])
  #[c(trueAZclass-min(trueAZclass,na.rm=TRUE))/(max(trueAZclass,na.rm=TRUE)-min(trueAZclass,na.rm=TRUE))])
  #color.legend(max(centers[,1]+1000),min(centers[,2]),max(centers[,1])+2000,max(centers[,2]), col.labels, testcol, gradient="y", cex=0.6, align="rb")
  dev.off()

  pdf(file=paste(WD,"21_Legend_pie.pdf",sep=""))
  pie(rep(1,8),col=palette,labels=c("N","NE","E","SE","S","SW","W","NW"), clockwise=T, init.angle=112.5)
  #[c(trueAZclass-min(trueAZclass,na.rm=TRUE))/(max(trueAZclass,na.rm=TRUE)-min(trueAZclass,na.rm=TRUE))])
  #color.legend(max(centers[,1]+1000),min(centers[,2]),max(centers[,1])+2000,max(centers[,2]), col.labels, testcol, gradient="y", cex=0.6, align="rb")
  dev.off()

  cat("Histogram of aspects of component true aspect", fill=TRUE)
  pdf(file=paste(WD,"21_HistogramAspects.pdf",sep=""))
  hist(trueAZclass, breaks=c(0.5,1.5,2.5,3.5,4.5,5.5,6.5,7.5,8.5), axes=FALSE, col=palette)
  axis(1,c(1,2,3,4,5,6,7,8), c("N","NE","E","SE","S","SW","W","NW"))
  axis(2)
  dev.off()

  ### GREYSCALE ASPECTS
  trueAZclass[which(trueAZ>=  0     & trueAZ <   11.25)] = 1   #GREY75
  trueAZclass[which(trueAZ>= 11.25  & trueAZ <   22.5 )] = 2  #GREY81
  trueAZclass[which(trueAZ>= 22.5   & trueAZ <   33.75)] = 3 #GREY88
  trueAZclass[which(trueAZ>= 33.75  & trueAZ <   45   )] = 4   #GREY94
  trueAZclass[which(trueAZ>= 45     & trueAZ <=  56.25)] = 5   # WHITE
  trueAZclass[which(trueAZ>= 56.25  & trueAZ <   67.5 )] = 6 #GREY94
  trueAZclass[which(trueAZ>= 67.5   & trueAZ <   78.75)] = 7 #GREY88
  trueAZclass[which(trueAZ>= 78.75  & trueAZ <   90   )] = 8   #GREY81
  trueAZclass[which(trueAZ>= 90     & trueAZ <  101.25)] = 9  #GREY75
  trueAZclass[which(trueAZ>= 101.25 & trueAZ <  112.5 )] = 10 #GREY69
  trueAZclass[which(trueAZ>= 112.5  & trueAZ <  123.75)] = 11 #GREY62
  trueAZclass[which(trueAZ>= 123.75 & trueAZ <  135   )] = 12   #GREY56
  trueAZclass[which(trueAZ>= 135    & trueAZ <  146.25)] = 13   #GREY50
  trueAZclass[which(trueAZ>= 146.25 & trueAZ <  157.5 )] = 14 #GREY44
  trueAZclass[which(trueAZ>= 157.5  & trueAZ <  168.75)] = 15 #GREY38
  trueAZclass[which(trueAZ>= 168.75 & trueAZ <  180   )] = 16   #GREY31
  trueAZclass[which(trueAZ>=-180    & trueAZ < -168.75)] = 17 #GREY25
  trueAZclass[which(trueAZ>=-168.75 & trueAZ < -157.5 )] = 18 #GREY19
  trueAZclass[which(trueAZ>=-157.5  & trueAZ < -146.25)] = 19 # GREY12
  trueAZclass[which(trueAZ>=-146.25 & trueAZ < -135   )] = 20   #GREY6
  trueAZclass[which(trueAZ>=-135    & trueAZ < -123.75)] = 21 # BLACK
  trueAZclass[which(trueAZ>=-123.75 & trueAZ < -112.5 )] = 22 #GREY6
  trueAZclass[which(trueAZ>=-112.5  & trueAZ < -101.25)] = 23 #GREY12
  trueAZclass[which(trueAZ>=-101.25 & trueAZ <  -90   )] = 24    #GREY19 
  trueAZclass[which(trueAZ>= -90    & trueAZ <  -78.75)] = 25   #GREY25
  trueAZclass[which(trueAZ>= -78.75 & trueAZ <  -67.5 )] = 26 #GREY31
  trueAZclass[which(trueAZ>= -67.5  & trueAZ <  -56.25)] = 27 #GREY38
  trueAZclass[which(trueAZ>= -56.25 & trueAZ <  -45   )] = 28 #GREY44
  trueAZclass[which(trueAZ>= -45    & trueAZ <  -33.75)] = 29 #GREY50
  trueAZclass[which(trueAZ>= -33.75 & trueAZ <  -22.5 )] = 30 #GREY56
  trueAZclass[which(trueAZ>= -22.5  & trueAZ <  -11.25)] = 31 #GREY62
  trueAZclass[which(trueAZ>= -11.25 & trueAZ <=   0   )] = 32  #GREY69

  palette1=c("grey75","grey81","grey88","grey94","white","grey94","grey88","grey81","grey75","grey69","grey62","grey56","grey50","grey44","grey38","grey31","grey25","grey19","grey12","grey6",
             "grey0","grey6","grey12","grey19","grey25","grey31","grey38","grey44","grey50","grey56","grey62","grey69")
  testcol = rainbow(8)  # topo.colors, terrain.colors, cm.colors, heat.colors, rainbow

  pdf(file=paste(WD,"20_True_AZ_color_grey32.pdf",sep=""))
  plot(clipTIN, lwd=0.05,col=palette1[c(trueAZclass)])
  #[c(trueAZclass-min(trueAZclass,na.rm=TRUE))/(max(trueAZclass,na.rm=TRUE)-min(trueAZclass,na.rm=TRUE))])
  #color.legend(max(centers[,1]+1000),min(centers[,2]),max(centers[,1])+2000,max(centers[,2]), col.labels, testcol, gradient="y", cex=0.6, align="rb")
  dev.off()
}

## SKY VIEW FRACTION + REMOTE SHADING MODULE
Shading=matrix(90,dim(TINCENTERS)[1],slices)  # it initially assumed all skyview fractions are 90 degree for the calculation of remote shadfing

Sys.time()
cat("Calculating Sky View Fraction for",dim(TINCENTERS)[1]," points", fill=TRUE)

SVF      = rep(NA,dim(TINCENTERS)[1])
lslices  = 360/slices
anslices = seq(0,360,lslices)

## inclusion of vertices
vertvec1    = cbind(vertex[,1],vertex[,4],vertex[,7])
vertvec2    = cbind(vertex[,2],vertex[,5],vertex[,8])
vertvec3    = cbind(vertex[,3],vertex[,6],vertex[,9])
vertvecalli = rbind(vertvec1,vertvec2, vertvec3)
indicesdup  = which(duplicated(vertvecalli) == TRUE)
vertvecall  = vertvecalli[-indicesdup,]   # removing duplicate vertices because one triangle shares some vertices with other traingles

for (tn in 1:dim(TINCENTERS)[1]){
  cat("Element",tn, " of ",dim(TINCENTERS)[1],fill=TRUE)
  skys           = rep(NA,slices)
  ctnx           = TINCENTERS[tn,1]
  ctny           = TINCENTERS[tn,2]
  ctnz           = TINCENTERS[tn,3]
  newTINCENTERSx = vertvecall[,1]-ctnx
  newTINCENTERSy = vertvecall[,2]-ctny
  newTINCENTERSz = vertvecall[,3]-ctnz
  
  angle                 = (((newTINCENTERSx/abs(newTINCENTERSx))*(pi/2))-atan(newTINCENTERSy/newTINCENTERSx))*180/pi  # gives angles between -180 and 180
  angle[which(angle<0)] = angle[which(angle<0)]+360   # gives angles between 0 and 360
  dist                  = sqrt(((newTINCENTERSx)^2)+(newTINCENTERSy)^2)
  radius                = max(dist)
  dhdx                  = newTINCENTERSz/dist
  
  for (angi in 1:slices){
    #cat("angi=",angi,fill=TRUE)
    inside=which(newTINCENTERSz>=0 & angle>anslices[angi] & angle<=anslices[angi+1])
    coordinside=cbind(newTINCENTERSx[inside],newTINCENTERSy[inside],newTINCENTERSz[inside])
    
    if (length(inside)>0){
      #plot(c(0,newTINCENTERSx[inside]),c(0,newTINCENTERSy[inside]),pch=20)
      #points(c(0,0),col="orange",pch=15)
      dhdxin = dhdx[inside]
      for (siga in 1:length(dhdxin)){
        maxdhdx   = sort(dhdxin,decreasing=TRUE)[siga]
        dondemaxi = which(dhdxin == maxdhdx)
        blocking  = inside[dondemaxi]
        dc        = sqrt(((vertvecall[blocking,3])^2)-(ctnz^2)+(2*RE*(vertvecall[blocking,3]-ctnz)))
        
        if (dist[blocking]<dc){
          break
        }
        cat("blocking element is so far it can't be viewed from the distace", fill=TRUE) 
      }    
      #cat("tallest element is at", dist[blocking],"m distance, at an elevation of ",TINCENTERS[blocking,3]," and has elevation difference of ", newTINCENTERSz[blocking]," m, for Dh/Dx=",dhdx[blocking],fill=TRUE)
      skys[angi]       = atan(dist[blocking]/newTINCENTERSz[blocking])*180/pi
      Shading[tn,angi] = atan(dist[blocking]/newTINCENTERSz[blocking])*180/pi
    } else {
      cat("Blocking below horizon SKY>1",fill=TRUE)
      insidebelowhorizon      = which(newTINCENTERSz<0 & angle>anslices[angi] & angle<=anslices[angi+1])
      coordinsidebelowhorizon = cbind(newTINCENTERSx[insidebelowhorizon],newTINCENTERSy[insidebelowhorizon],newTINCENTERSz[insidebelowhorizon])
      
      if (length(insidebelowhorizon)>0){
        dzdxin     = dhdx[insidebelowhorizon]  
        mindzdx    = sort(dzdxin,decreasing=FALSE)[1]
        dondemini  = which(dzdxin == mindzdx)
        blocking   = insidebelowhorizon[dondemini][1]
        skys[angi] = 90+ atan(abs(newTINCENTERSz[blocking])/dist[blocking])*180/pi
      } 
      
      if (length(insidebelowhorizon)==0){
        skys[angi]=90  # flat horizon
      }
    }
  }  
  SVF[tn]=mean(skys)/90
}
  
ordersky=cbind(TINCENTERS,SVF)
save(ordersky,file=paste(WD,"ordersky.rda",sep=""))  # uncomment this if need to save ordersky
save(Shading,file=paste(WD,"Shading.rda",sep=""))  # uncomment this if need to save ordersky

if (makeplots==1){
  pdf(file=paste(WD,"22_Sky_view_fractions1.pdf",sep=""))
  #testcol<-grey(seq(0,1, len=128))   # grey pallete
  testcol = terrain.colors(256) # topo.colors, terrain.colors, cm.colors, heat.colors, rainbow, ygobb, matlab.like2, primary.colors, diverge_hsv, cyan2yellow, green2red, blue2red, magenta2green 
  col.labels<-c(round(min(ordersky[,4]),digits=2),round((((max(ordersky[,4])-min(ordersky[,4])))/4)+min(ordersky[,4]), digits=2),round((((max(ordersky[,4])-min(ordersky[,4])))/2)+min(ordersky[,4]), digits=2),
              round((((max(ordersky[,4])-min(ordersky[,4])))*0.75)+min(ordersky[,4]), digits=2),round(max(ordersky[,4]), digits=2))
  #plot(clipTIN, lwd=0.05,col=grey((trueslope-min(trueslope))/(max(trueslope)-min(trueslope))))
  plot(clipTIN, lwd=0.05,col=testcol[c(round(255*(ordersky[,4]-min(ordersky[,4]))/(max(ordersky[,4])-min(ordersky[,4]))))+1])
  color.legend(max(centers[,1]+1000),min(centers[,2]),max(centers[,1])+2000,max(centers[,2]), col.labels, testcol, gradient="y", cex=0.6, align="rb")
  dev.off()
}
