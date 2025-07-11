#>> We will not include plotting functions due to their size
#>> Similar sections to those already included will be omitted

options(digits = 8, max.print = 100) # console format

# install.packages("terra")
library(terra)
# install.packages("performance")
library(performance)
# install.packages("gstat")
library(gstat)

## ==== Data processing ====

### ---- Covariables ----

#> Import data and create rasters
Altitude.frast <- rast("./features/GIMP/gimpdem_90m_v01.1.tif")
crs(Altitude.frast) <- crs("EPSG:3413")

gravity <- read.table("./features/gravity_anomaly/grid_bouguer.txt")
gravity_lonlat.rast <- rast(gravity, type = "xyz", crs = "EPSG:4326")
Gravity_Anomaly.frast <- project(gravity_lonlat.rast, "EPSG:3413")

Mass_Balance.frast <- rast("./features/surface_mass_balance/smb_rec-mean.nc")
crs(Mass_Balance.frast) <- crs("EPSG:3413")

bedmachine.rast <- rast("./measurements/BedMachine/BedMachineGreenland-v5.nc")
bedmachine_mask.rast <- bedmachine.rast$mask
bedmachine_mask.rast[bedmachine.rast$mask == 2] <- 1
Distance_to_Coast.frast <- distance(bedmachine_mask.rast, target = 1)/1e3
crs(Distance_to_Coast.frast) <- crs("EPSG:3413")

icevelx.rast <- rast("./features/ice_velocity/greenland_vel_mosaic250_vx.tif")
icevely.rast <- rast("./features/ice_velocity/greenland_vel_mosaic250_vy.tif")
Ice_Velocity.frast <- lapp(c(icevelx.rast, icevely.rast), 
                                  fun = function(x, y){ return(sqrt(x^2+y^2)) })
crs(Ice_Velocity.frast) <- crs("EPSG:3413")

#> Organize into a list for efficiency
raw_covariables.rastlist <- mget(ls(pattern = "\\.frast$"))
names(raw_covariables.rastlist) <- 
  gsub(".frast", "", names(raw_covariables.rastlist))

#> Transform to the extension and resolution of the one with smallest resolution
name_biggest <- names(which.max(sapply(raw_covariables.rastlist, 
                                       function(rast){ prod(res(rast)) })))
resampled_covariables.rastlist <- lapply(raw_covariables.rastlist,function(rast) 
  resample(rast, raw_covariables.rastlist[[name_biggest]], method = "bilinear"))

#> Create grid from which to calculate DistCoast before adjusting all to it
template.rast <- trim(!is.na(app(rast(resampled_covariables.rastlist), 
                                 fun = sum)))
covariables.rastlist <- lapply(resampled_covariables.rastlist, function(rast) 
  trim(mask(crop(rast, template.rast), template.rast, maskvalues = FALSE)) )

names(covariables.rastlist) <- 
  c("Altitude", "GravAnom", "MassBal", "DistCoast", "IceVel")
covariables.rast <- rast(covariables.rastlist)

### ---- Variable of interest ----

#> Import .csv files to a list and name them by expedition
temp <- list.files(path = "./measurements/CReSIS/csv_files", 
                   pattern = "\\.csv$", full.names = TRUE)
CReSIS.dflist <- lapply(temp, read.csv)
names(CReSIS.dflist) <- gsub("_Greenland|.csv", "", substring(temp, 40))

#> Reduce by contributed percentage
reduced_CReSIS.dflist <- 
  CReSIS.dflist[which(prop.table(sapply(CReSIS.dflist, nrow)) > .024)]

#> For each expedition, transform the relevant columns to SpatVectors, omitting
#> NAs and using the correct crs before projecting to polar stereographic.
CReSIS.vectlist <- lapply(reduced_CReSIS.dflist, function(dataframe){ 
  project(vect(na.omit(dataframe[c("LON", "LAT", "THICK")]), 
               geom = c("LON", "LAT"), crs = "EPSG:4326"), "EPSG:3413") })

#> Remove outliers
CReSIS_outliers <-performance::check_outliers(vect(CReSIS.vectlist)$Thick)
CReSIS_outliers # detected 7 outliers
split_index <- split(!CReSIS_outliers, 
                     unlist(sapply(1:length(CReSIS.vectlist), function(k) 
                       rep(k, length(CReSIS.vectlist[[k]])))))
CReSIS.vectlist <- lapply(
  setNames(1:length(CReSIS.vectlist), names(CReSIS.vectlist)), function(k) 
    CReSIS.vectlist[[k]][split_index[[k]]])


#> Remove ice thickness <= 0
CReSIS_mod.vect <- vect(CReSIS.vectlist)[vect(CReSIS.vectlist)$THICK > 0,]

#> Adjust to the grid (Altitude == template.rast)
CReSIS_mean.vect <- {
  rast <- rasterize(CReSIS_mod.vect, covariables.rast$Altitude, 
                    field = "THICK", fun = mean)
  vect <- as.points(mask(rast, covariables.rast$Altitude, maskvalues = NA))
  names(vect) <- "THICK"; vect
}

#> Create observation matrix
CReSIS_data.df <- 
  as.data.frame(extract(covariables.rast, CReSIS_mean.vect,
                        ID = FALSE, na.rm = TRUE, bind = TRUE, xy = TRUE))
names(CReSIS_data.df)[1] <- "Thickness"

## ==== Model training ====

#> Choice of data source
data.df <- CReSIS_data.df

#> Split into training and testing/evaluation data
proportion <- c(train = .7, test = .3)
n <- 5e4; n/nrow(data.df)

{set.seed(321)
datasample <- data.df[sample(1:nrow(data.df), min(n,nrow(data.df)), replace=F),]
assignation <- sample(cut(seq(nrow(datasample)),
                          nrow(datasample)*cumsum(c(0, proportion)),
                          labels = names(proportion) ))
split_datasample <- split(datasample, assignation)
TrainingData.df <- split_datasample$train
TestingData.df <- split_datasample$test}

#> Hipothesis testing with performance
performance::check_model(lm(
  formula="Thickness ~ Altitude + GravAnom + MassBal + DistCoast + log(IceVel)",
  data = TrainingData.df), size_title = 24, base_size = 20)

##### Ordinary Kriging 
formula = "Thickness ~ 1"; model = "OK"
psill <- NA; range <- 3e4; nugget <- 1e4; vmodel <- "Lin"; anis <- c(0, .8)

##### Universal Kriging
formula = "Thickness ~ Altitude + GravAnom + MassBal + DistCoast + log(IceVel)"
model = "UK"
psill <- NA; range <- 1e5; nugget <- 1e4; vmodel <- "Exc"; anis <- c(0, 1)

cutoff <- max(TrainingData.df$DistCoast)*1e3/2
nbins <- 20

#> Experimental semivariogram
variogram <- variogram(as.formula(formula), ~x+y, data = TrainingData.df,
                       cutoff = cutoff, width = cutoff/nbins)

variogram_map <- variogram(gstat(id = "Map", formula = as.formula(formula),
                                 locations = ~x+y, data = TrainingData.df), 
                           map = T, cutoff = cutoff, width = cutoff/nbins)

#> Adjust theoretical model
theoretical_variogram <- fit.variogram(variogram, debug.level = 2, fit.method=2,
                                       vgm(psill = psill, vmodel, range = range, 
                                           nugget = nugget, anis = anis))
anisotropy_variogram <- variogram(as.formula(formula), ~x+y, 
                                  data = TrainingData.df, cutoff = cutoff, 
                                  width = cutoff/nbins,
                                  alpha = rev(seq(0, 180, by = 30))[-1])

assign(paste0(model, ".thvar"), theoretical_variogram)


## ==== Prediction ====

#> Kriging model construction and prediction
nmin <- 1; nmax <- 100; maxdist <- 1e5
OK.model <- gstat(id = "OK", formula = as.formula("Thickness ~ 1"), 
                  data = TrainingData.df, locations = ~x+y, model = OK.thvar, 
                  nmin = nmin, nmax = nmax, maxdist = maxdist, set =list(gls=1))
OK.pred <- interpolate(covariables.rast, OK.model, cores = 4, cpkgs = "gstat", 
                       debug.level = 2, na.rm = TRUE, index = 3:4)

nmin <- 1; nmax <- 100; maxdist <- 2e5
UK.model <- gstat(id = "UK", formula = as.formula(
  "Thickness ~ Altitude + GravAnom + MassBal + DistCoast + log(IceVel)"), 
  data = TrainingData.df, locations = ~x+y, model = UK.thvar, nmin = nmin, 
  nmax = nmax, maxdist = maxdist, set = list(gls = 1))
UK.pred <- interpolate(covariables.rast, UK.model, cores = 4, cpkgs = "gstat", 
                       debug.level = 2, na.rm = TRUE, index = 3:4)

#> Other models construction and prediction
MLR.model <- lm(formula = as.formula(
  "Thickness ~ Altitude + GravAnom + MassBal + DistCoast + log(IceVel)"), 
  data = TrainingData.df)
MLR.pred <- predict(covariables.rast, MLR.model, na.rm = TRUE)

IDW.model <- gstat(id = "IDW", formula = as.formula(
  "Thickness ~ Altitude + GravAnom + MassBal + DistCoast + log(IceVel)"),
  data = TrainingData.df, locations = ~x+y, nmin = nmin, nmax = nmax, 
  maxdist = maxdist, set = list(gls = 1, idp = 2))
IDW.pred <- interpolate(covariables.rast, IDW.model, debug.level = 2, 
                        na.rm = TRUE, index = 3:4)

BM.pred <- trim(mask(
  resample(bedmachine.rast$thickness, covariables.rast$Altitude, 
           method = "bilinear"), covariables.rast$Altitude, maskvalues = NA))

### ---- Comparación ----
#> Metric functions
RMSE <- function(prediction, measurements = TestingData.df$Thickness){
  return(sqrt(mean((prediction - measurements)^2)))
}
R2 <- function(prediction, measurements = TestingData.df$Thickness){
  RSS <- sum((prediction - measurements)^2)
  TSS <- sum((measurements - mean(measurements))^2)
  
  return(1 - RSS/TSS)
}

#> Volume and area function
SLEandArea <- function(thickness.rast){
  require(terra)
  
  thick_pos.rast <- mask(thickness.rast, thickness.rast > 0, 
                         maskvalue = FALSE)
  
  area.km2 <- prod(res(thickness.rast)/1e3)
  extent.percent <- (area.km2 * global(!is.na(thick_pos.rast), sum)) / 
    (area.km2 * global(!is.na(thickness.rast), sum)) * 100
  
  volume.km3 <- as.numeric(global(area.km2 * thick_pos.rast/1e3, 
                                  sum, na.rm =T))
  print("Volume (km3):"); print(volume.km3)
  Gt <- volume.km3 * 0.9167 # 0.9167 density of ice
  print("Gt:"); print(Gt)
  
  SLE.m <- Gt * 1/361.8 # 361.8 Gt ice = 1 mm sea level elevation
  
  return(list("SLE" = SLE.m/1e3, "Extent" = extent.percent))
}

#> Evaluate all models
TestingData.vect <- vect(TestingData.df[c("Thickness", "x", "y")], 
                         geom = c("x","y"), crs= crs(covariables.rast$Altitude))
models <- c("UK_197X.pred","MLR.pred","IDW.pred","OK.pred","UK.pred","BM.pred")
Test_results.df <- cbind(sapply(mget(models), function(rast) 
  extract(rast[[1]], TestingData.vect, ID = FALSE)[,1] ))

#> Apply metrics
Metrics <- cbind(apply(Test_results.df, 2, RMSE), apply(Test_results.df, 2, R2))
colnames(Metrics) <- c("RMSE", "R2")
rownames(Metrics) <- gsub(".pred", "", models)
Metrics

#> Calculate volume and area
Calcs <- sapply(models, function(model) unlist(SLEandArea(get(model)[[1]])))
rownames(Calcs) <- c("SLE", "Extent")
colnames(Calcs) <- gsub(".pred", "", models)
t(Calcs)



### ----  Anexo I ----

# my.plot_terra_maps(list("Máscara BedMachine" = bedmachine.rast$mask), resolution = c(800, 1200), filename = "./LaTeX/Images/bedmachine-mask")


## ==== Examples =====
library(sp)
data(meuse)
coordinates(meuse) = ~x+y
v <- variogram(log(zinc) ~ 1, meuse)
v.fit <- fit.variogram(v, vgm(1, "Sph", 700, 1))
v.fit
png("./LaTeX/Images/meuse-variogram-example.png", width = 1200, height = 600)
plot(v, v.fit, ylab = list("Semivarianza", cex = 1.2), xlab = list("Distancia (m)", cex = 1.2), par.settings = list(fontsize = list(text = 25, points = 15)))
library(lattice)
trellis.focus("panel",1,1)
panel.arrows(0, sum(v.fit[,2]), v.fit[2,3], sum(v.fit[,2]), lwd = 3, length = .2, angle = 20)
panel.arrows(v.fit[2,3], sum(v.fit[,2]), 0, sum(v.fit[,2]), lwd = 3, length = .2, angle = 20)
ltext(x=v.fit[2,3]/2, y=1.04*sum(v.fit[,2]), "Alcance", cex=1)
panel.arrows(v.fit[2,3], sum(v.fit[,2]), v.fit[2,3], 0, lwd = 3, length = .2, angle = 20)
panel.arrows(v.fit[2,3], 0, v.fit[2,3], sum(v.fit[,2]), lwd = 3, length = .2, angle = 20)
ltext(x=v.fit[2,3]+10, y=sum(v.fit[,2])/2, "Umbral", cex=1, adj=0)
llines(x=c(0,30), y=v.fit[1,2], col = "black", lwd=1, lty=3)
panel.arrows(30, 0, 30, v.fit[1,2], lwd = 3, length = .15, angle = 20)
panel.arrows(30, v.fit[1,2], 30, 0, lwd = 3, length = .15, angle = 20)
ltext(x=40, y=v.fit[1,2]/1.8, "Pepita", cex=1, adj=0)
trellis.unfocus()
dev.off()

# show.vgms()
png("./LaTeX/Images/theoretical-variograms-example.png", width = 1200, height = 500)
show.vgms(models = c("Nug", "Lin", "Sph", "Exp"), max = 1,
          xlab = list("Distancia (m)"), ylab = list("Semivarianza"),
          par.settings = list(fontsize = list(text = 30, points = 20)))
dev.off()


# Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jre1.8.0_431/")
# library(glmulti)






######### cementery

# CReSIS_mean.vectlist <- lapply(CReSIS.vectlist, function(vect) {
#   rast <- c(rasterize(vect, template.rast, field = "THICK", fun = mean),
#             rasterize(vect, template.rast, fun = "count"),
#             rasterize(vect, template.rast, field = "y", fun = function(x) if(all(x%in%x)){x}else{NA}));
#   vect <- as.points(mask(rast, template.rast, maskvalues = FALSE)); names(vect) <- c("Thickness", "count"); vect
#   })

#> MERGED SEPARATELY
# my.plot_terra_maps_histograms_boxplots(list("1971-1979" = Seventies_mean.vect, CReSIS = CReSIS_mean.vect), "Thickness", 
#                                      background.rast = covariables.rast$Altitude, palette = "Temps", consistent_colors = FALSE, 
#                                      global_x = "Thickness del hielo", byrow = TRUE, resolution = c(1600, 1200), 
#                                      filename = "merged_Seventies_CReSIS_maps_histograms_boxplots")

# CReSIS_data.dflist <- lapply(CReSIS_mean.vectlist, function(vect) 
#   as.data.frame(na.omit(extract(covariables.rast, vect, ID = FALSE, na.rm = TRUE, bind = TRUE, xy = TRUE), field = "")) )
# head(CReSIS_data.dflist$`2005_TO`)

# xyplot(gamma ~ dist | id, data = v, layout = c(2, 1),
#        main = list("Semivariograma no estacionario"),
#        xlab = list("Distancia (m)"), ylab = list("Semivarianza"),
#        par.settings = list(fontsize = list(text = 20, points = 15)))

#> SCATTERPLOTS
# my.plot_scatterplots(TrainingData.df$Thickness, "Thickness", TrainingData.df[2:6], palette = "Plasma",
#                      resolution = c(1600, 1200), filename = "./LaTeX/Images/sample-covariable-scatterplots")
# dividir en dos modelos cortando distancia a la costa a la mitad?

# my.plot_crossplot_matrix(TrainingData.df[1:6], palette = "Plasma", resolution = c(1600,1200), filename = "./LaTeX/Images/sample-scatterplot-matrix")
