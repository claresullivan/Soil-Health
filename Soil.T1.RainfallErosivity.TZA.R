####################################################
# Formatting Rainfall Erosivity for use as RUSLE R Factor
# Date - August 12, 2015
####################################################
# resamples from 25 miles resolution to 250 m resolution
library(sp)
library(raster)
library(GISTools)
library(rgdal)

# data sources
dir_data1 <- "C:\\Data\\Vital Signs\\Soil Health Thread\\Inputs\\ErosivityAfrica"
dir_data2 <- "C:\\Data\\Vital Signs\\Soil Health Thread\\Inputs\\Boundaries"

# data output
dir_out1 <- "C:\\Data\\Vital Signs\\Soil Health Thread\\Outputs"

# load data 
setwd(dir_data1)
# erosivity dataset is from Vrieling 2014
ero <- raster("R_EI30rescaled_TMPA1998-2012_Africa_geotiff.tif")
setwd(dir_data2)
sagcot <- readOGR(dir_data2, "SAGCOT")
setwd(dir_out1)
sagcot_ph <- raster("sagcot_ph.tif")

# Reproject all files to Africa Albers Equal Area Conic (ESRI:102022)
aea <- '+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0'
# reproject SAGCOT boundary shapefile )
sagcot <- spTransform(sagcot, crs(aea))
# reproject erosivity: 
ero <- projectRaster(ero, crs = aea)
# reproject reference raster pH raster file for resampling
sagcot_aea <- projectRaster(sagcot_ph, crs = aea)

# Crop using extent, rasterize polygon and finally, create poly-raster
e <- extent(sagcot)
ero.sagcot <- crop(ero, e, snap ="out")
rsagcot <- rasterize(sagcot, ero.sagcot)
ero <- mask(x = ero.sagcot, mask = rsagcot)

# Set NA to "-999" values
ero[-999] <- NA

# use the pH layer to create a reference for resampling:
sagcot_aea
#dimensions  : 2177, 4083, 8888691  (nrow, ncol, ncell)
#resolution  : 235, 266  (x, y)
#extent      : 554167.7, 1513673, -1254644, -675562  (xmin, xmax, ymin, ymax)
#coord. ref. : +proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 
#+ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0 
sagcot_ref <- raster(nrow=2177, ncol=4083, crs=aea, res=c(235, 266), extent(c(554167.7, 1513673, -1254644, -675562)))
sagcot_ref <- setValues(sagcot_ref, NA)

# Note: same extent, same projection, but different resolution
ero 
sagcot_aea

# Set up higher resolution for erosivity (apply a grid with MORE pixels within the same extent) and resample:
#ero2 <- raster(nrow=2177, ncol=4083)
ero2 <- resample(ero, sagcot_ref, method='bilinear')
ero2[ero2<0] <- NA

#write outputs
setwd(dir_out1)
writeRaster(ero2, filename="erosivity_sagcot.tif", format='GTiff', overwrite=TRUE)

####################################################################################

