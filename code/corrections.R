#required packages
require(raster)
require(rgdal)
require(landsat8)

#set work directory
setwd("D:/Users/gjors_000/Documents/Masterarbeit/Landsat Corrections")

#Data input

#Landsat 8 Thermalband

#Load MTL File
mtl = read.delim("data/LS8/LC81920252015237LGN00/LC81920252015237LGN00_MTL.txt", sep = '=', stringsAsFactors = F)

#Load TIRS files
ls8band10 = toString(mtl[grep("FILE_NAME_BAND_10",mtl$GROUP),]["L1_METADATA_FILE"])
ls8band11 = toString(mtl[grep("FILE_NAME_BAND_11",mtl$GROUP),]["L1_METADATA_FILE"])

ls8band10 = substr(ls8band10, 2, nchar(ls8band10))
ls8band11 = substr(ls8band11, 2, nchar(ls8band11))

ls8band10 = readGDAL(paste("data/LS8/LC81920252015237LGN00/",ls8band10, sep=""))
ls8band11 = readGDAL(paste("data/LS8/LC81920252015237LGN00/",ls8band11, sep=""))

#Load temperature conversion variables

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


#Conversion to At Satellite Brightness Temperature
BTls8band10=tempconv(ls8band10,Ml10,Al10,K110,K210)
BTls8band11=tempconv(ls8band11,Ml11,Al11,K111,K211)
