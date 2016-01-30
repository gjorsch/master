#libs
require(raster)
require(rgdal)
require(landsat8)
require(sp)
require(rgeos)


#set work directory
setwd("D:/Users/gjors_000/Documents/Masterarbeit/Landsat Corrections")

#Data input

#Landsat 8 Thermalband

#Load MTL File and DWD Data
mtlfile = "LC81920252015157LGN00_MTL.txt"
emissitivityfile= "emiss_ETRS_align_clip_GRASS.tif"
LS8fileband10 = "LS8b10_ETRS_align_clip_GRASS.tif"
LS8fileband11 = "Band11_ETRS_align_clip.tif"

mtl = read.delim(paste0("data/LS8/",mtlfile), sep = '=', stringsAsFactors = F)
dwd = read.delim("data/DWD/produkt_temp_Terminwerte_20140712_20160112_01048.txt", sep = ';', stringsAsFactors = F)
dwdp = read.delim("data/DWD/pressure/produkt_synop_Terminwerte_20140712_20160112_01048.txt", sep = ';', stringsAsFactors = F)
date = toString(mtl[grep("FILE_DATE",mtl$GROUP),]["L1_METADATA_FILE"])

#substring date to match DWD format
date = gsub("-","",date)
date = sub("T","",date)
date = substr(date,2,11)

#Load TIRS files
ls8band10 = readGDAL(paste("data/LS8/",LS8fileband10, sep=""),band=1)
ls8band11 = readGDAL(paste("data/LS8/",LS8fileband11, sep=""),band=1)

#Load emmissitivity file
ls8emiss = readGDAL(paste("data/LS8/",emissitivityfile, sep=""),band=1)

#band specific multiplicative rescaling factor
Ml10 = as.numeric(mtl[grep("RADIANCE_MULT_BAND_10",mtl$GROUP),]["L1_METADATA_FILE"])
Ml11 = as.numeric(mtl[grep("RADIANCE_MULT_BAND_11",mtl$GROUP),]["L1_METADATA_FILE"])

#band specific additive rescaling factor 
Al10 = as.numeric(mtl[grep("RADIANCE_ADD_BAND_10",mtl$GROUP),]["L1_METADATA_FILE"])
Al11 = as.numeric(mtl[grep("RADIANCE_ADD_BAND_11",mtl$GROUP),]["L1_METADATA_FILE"])

#band specific thermal conversion constants
K110 = as.numeric(mtl[grep("K1_CONSTANT_BAND_10",mtl$GROUP),]["L1_METADATA_FILE"])
K111 = as.numeric(mtl[grep("K1_CONSTANT_BAND_11",mtl$GROUP),]["L1_METADATA_FILE"])
K210 = as.numeric(mtl[grep("K2_CONSTANT_BAND_10",mtl$GROUP),]["L1_METADATA_FILE"])
K211 = as.numeric(mtl[grep("K2_CONSTANT_BAND_11",mtl$GROUP),]["L1_METADATA_FILE"])

#Radiance transfer coeficients
C1 = 14387.7
C2 = 1.19104*10^8

#effective Wavelength for Landsat8 TIRS
LS8band10lambda = 10.896
LS8band11lambda = 12.006

#near LST temperature (°C) and relative humidity
Tnear = as.numeric(dwd[grep(date,dwd$MESS_DATUM),]["LUFTTEMPERATUR"])
Rhumid = as.numeric(dwd[grep(date,dwd$MESS_DATUM),]["REL_FEUCHTE"])
Rhumid = Rhumid/100

#air pressure
p = as.numeric(dwdp[grep(date,dwdp$MESS_DATUM),]["LUFTDRUCK_STATIONSHOEHE"])

#atmospheric water vapour content
w = Rhumid*((1.0007+(3.46*10^(-6)*p))*(6.1121*exp((17.502*Tnear)/(240.97+Tnear))))
w = 0.098 * w

#omega coefficients band 10 (model of atmosphere)
omega1_10 <- function(x){
  y = (0.0109*x^3)+(0.0079*x^2)+(0.0991*x)+1.0090
  return(y)
}
omega2_10 <- function(x){
  y = (-0.0620*w^3)+(-0.4671*w^2)+(-1.2105*w)+0.1176
  return(y)
}
omega3_10 <- function(x){
  y = (-0.0533*w^3)+(0.4013*w^2)+(0.8585*w)-0.0451
  return(y)
}

#omega coefficients band 11 (model of atmosphere)
omega1_11 <- function(x){
  y = (0.0405*x^3)+(-0.0809*x^2)+(0.2919*x)+0.9620
  return(y)
}
omega2_11 <- function(x){
  y = (-0.2960*w^3)+(0.3611*w^2)+(-1.0257*w)+0.4644
  return(y)
}
omega3_11 <- function(x){
  y = (-0.0443*w^3)+(0.2509*w^2)+(1.4573*w)-0.0854
  return(y)
}

#correction functions (part of main function)
gamma <- function(x,lambda,T){
  y = ((C1*x)/(T*T))*((((lambda^4)/C2)*x)+lambda^(-1))
  y = 1/y
  return(y)
}
delta <- function(x,lambda,T){
  y = (-gamma(x,lambda,T)*x)+T
  return(y)
}
#correction main function
corr <- function(x,lambda,T){
  y = (gamma(x,lambda,T)*((1/e)*(((omega1*x)+omega2)+omega3)))+delta(x,lambda,T)
  return(y)
}

#LST for Band 10
LS8LST <- function(x,lambda,T,e,w){
  y = (gamma(x,lambda,T)*((1/e)*((omega1_10(w)*x)+omega2_10(w))+omega3_10(w)))+delta(x,lambda,T)
  return(y)
}

#Conversion to At Satellite Brightness Temperature
BTls8band10=tempconv(ls8band10,Ml10,Al10,K110,K210)
BTls8band11=tempconv(ls8band11,Ml11,Al11,K111,K211)

#Jimenez et al, 2003
#finally calculate band 10 LST spatial frame
LST = LS8LST(BTls8band10@data,LS8band10lambda,Tnear,ls8emiss@data,w)
LST_plot = BTls8band10
LST_plot@data = LST

#buggy due to little LST differences ~0.01°C
#check units!?
#writeGDAL(LST_plot, "outLST_LS8b10.tif")
plot(raster(LST_plot))

#Alternate calc
#Weng et al, 2003

BTls8band10alt = BTls8band10
BTls8band10alt@data =BTls8band10alt@data - 273.15
BTls8band11alt = BTls8band11
BTls8band11alt@data =BTls8band11alt@data - 273.15

LSTband10alt  = BTls8band10alt
LSTband10alt@data = BTls8band10alt@data/1+ls8band10@data*(BTls8band10alt@data/14380)*log(ls8emiss@data)
writeGDAL(LSTband10alt, "out_alt_LST_LS8b10.tif")

LSTband11alt  = BTls8band11alt
LSTband11alt@data = BTls8band11alt@data/1+ls8band11@data*(BTls8band11alt@data/14380)*log(ls8emiss@data)
writeGDAL(LSTband11alt, "out_alt_LST_LS8b11.tif")

LST_alt_av = LSTband10alt
LST_alt_av@data = (LSTband10alt@data+LSTband11alt@data)/2
writeGDAL(LST_alt, "out_alt_LST_LS8av.tif")