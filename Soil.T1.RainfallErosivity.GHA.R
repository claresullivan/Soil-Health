####################################################
# Formatting Rainfall Erosivity for use as RUSLE R Factor
# Date - August 12, 2015
####################################################
# 25 miles resolution to 250 m resolution
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
ero <- raster("R_EI30rescaled_TMPA1998-2012_Africa_geotiff.tif")
setwd(dir_data2)
gha <- readOGR(dir_data2, "GHA_adm1")
setwd(dir_out1)
gha_ph <- raster("gha_ph.tif")

# Reproject all files to Africa Albers Equal Area Conic (ESRI:102022)
aea <- '+proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 +ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0'
# reproject GHANA boundary shapefile )
gha <- spTransform(gha, crs(aea))
# reproject erosivity: 
ero <- projectRaster(ero, crs = aea)
# reproject reference raster pH raster file for resampling
gha_aea <- projectRaster(gha_ph, crs = aea)

# Crop using extent, rasterize polygon and finally, create poly-raster
e <- extent(gha)
ero.gha <- crop(ero, e, snap ="out")
rgha <- rasterize(gha, ero.gha)
ero <- mask(x = ero.gha, mask = rgha)

# Set NA to "-999" values
ero[-999] <- NA

# use the pH layer to create a reference for resampling:
gha_aea
#dimensions  : 2966, 2102, 6234532  (nrow, ncol, ncell)
#resolution  : 240, 260  (x, y)
#extent      : -2972559, -2468079, 538097.3, 1309257  (xmin, xmax, ymin, ymax)
#coord. ref. : +proj=aea +lat_1=20 +lat_2=-23 +lat_0=0 +lon_0=25 +x_0=0 +y_0=0 
#+ellps=WGS84 +datum=WGS84 +units=m +no_defs +towgs84=0,0,0 
gha_ref <- raster(nrow=2966, ncol=2102, crs=aea, res=c(240, 260), extent(c(-2972559, -2468079, 538097.3, 1309257)))
gha_ref <- setValues(gha_ref, NA)

# Note: same extent, same projection, but different resolution
ero 
gha_aea

# Set up higher resolution for erosivity (apply a grid with MORE pixels within the same extent) and resample:
#ero2 <- raster(nrow=2966, ncol=2102)
ero2 <- resample(ero, gha_ref, method='bilinear')
ero2[ero2<0] <- NA

#write outputs
setwd(dir_out1)
writeRaster(ero, filename="erosivity_gha_25miles.tif", format='GTiff', overwrite=TRUE)
writeRaster(ero2, filename="erosivity_gha.tif", format='GTiff', overwrite=TRUE)

####################################################################################

