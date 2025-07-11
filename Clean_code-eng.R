source("my_plotting_functions.R") # funciones caseras de gráficos
options(digits = 8, max.print = 100) # formato de consola

# install.packages("terra")
library(terra)
# install.packages("performance")
library(performance)
# install.packages("gstat")
library(gstat)

# writeRaster(covariables.rast, "./data_export/covariables.rast.tif", overwrite = T)
covariables.rast <- rast("./data_export/covariables.rast.tif")
# write.table(CReSIS_data.df, "./data_export/CReSIS_data.df.txt")
CReSIS_data.df <- read.table("./data_export/CReSIS_data.df.txt")

## ==== Procesado de datos ====

### ---- Covariables ----

#> Importamos los datos y creamos los rasters
Altitude.frast <- rast("./features/GIMP/gimpdem_90m_v01.1.tif")
crs(Altitude.frast) <- crs("EPSG:3413")

gravity <- read.table("./features/gravity_anomaly/grid_bouguer_-74_-11_59_84.txt")
gravity_lonlat.rast <- rast(gravity, type = "xyz", crs = "EPSG:4326")
Gravity_Anomaly.frast <- project(gravity_lonlat.rast, "EPSG:3413")

Mass_Balance.frast <- rast("./features/surface_mass_balance/smb_rec.1958-1990.BN_RACMO2.3p2_ERA5_3h_FGRN055.1km.YY-mean.nc")
crs(Mass_Balance.frast) <- crs("EPSG:3413")

bedmachine.rast <- rast("./measurements/BedMachine/BedMachineGreenland-v5.nc")
bedmachine_mask.rast <- bedmachine.rast$mask
bedmachine_mask.rast[bedmachine.rast$mask == 2] <- 1
Distance_to_Coast.frast <- distance(bedmachine_mask.rast, target = 1)/1e3
crs(Distance_to_Coast.frast) <- crs("EPSG:3413")

icevelx.rast <- rast("./features/ice_velocity/greenland_vel_mosaic250_vx_v1.tif")
icevely.rast <- rast("./features/ice_velocity/greenland_vel_mosaic250_vy_v1.tif")
Ice_Velocity.frast <- lapp(c(icevelx.rast, icevely.rast), fun = function(x, y){ return(sqrt(x^2+y^2)) })
crs(Ice_Velocity.frast) <- crs("EPSG:3413")

#> Los organizamos en una lista por eficiencia
raw_covariables.rastlist <- mget(ls(pattern = "\\.frast$")[c(1,3,5,2,4)]) # [c(1,3,5,2,4)] for reordering english
names(raw_covariables.rastlist) <- gsub(".frast", "", names(raw_covariables.rastlist))

#> RAW, histograms for outliers?
my.plot_terra_maps(raw_covariables.rastlist, palette = "Temps", numbering = c(1,2,2,3,4,5), byrow = TRUE,
                   resolution = c(1200, 1200), filename = "./LaTeX/Images/all-raw-covariable-maps")

#> Transformamos los importados a la extensión y resolución del de menos resolución
name_biggest <- names(which.max(sapply(raw_covariables.rastlist, function(rast){ prod(res(rast)) })))
resampled_covariables.rastlist <- lapply(raw_covariables.rastlist, function(rast) 
  resample(rast, raw_covariables.rastlist[[name_biggest]], method = "bilinear") )

#> Construimos la plantilla y distancia a la costa a partir de este antes de juntarlos
template.rast <- trim(!is.na(app(rast(resampled_covariables.rastlist), fun = sum)))
covariables.rastlist <- lapply(resampled_covariables.rastlist, function(rast) 
  trim(mask(crop(rast, template.rast), template.rast, maskvalues = FALSE)) )

names(covariables.rastlist) <- c("Altitude", "GravAnom", "MassBal", "DistCoast", "IceVel")
covariables.rast <- rast(covariables.rastlist)

#> PROCESSED, WITH HISTOGRAMS
my.plot_terra_maps_histograms_boxplots(c(setNames(as.list(covariables.rast), names(covariables.rast)), "log(IceVel)" = log(covariables.rast$IceVel)),
                                       palette = "Temps", ncol = 4, byrow = TRUE, resolution = c(1200, 1200), filename = "./LaTeX/Images/all-covariable-maps-histograms-boxplots")

### ---- Variable de interés ----

#> Importamos los archivos .csv a una lista y los nombramos de acuerdo a la expedición de la que provienen
temp <- list.files(path = "./measurements/CReSIS/csv_files", pattern = "\\.csv$", full.names = TRUE)
CReSIS.dflist <- lapply(temp, read.csv)
names(CReSIS.dflist) <- gsub("_Greenland|.csv","", substring(temp, 40))

#> Reducir por porcentaje contribuída
sum(prop.table(sapply(CReSIS.dflist, nrow)) > .024)
png("./LaTeX/Images/percentage-of-contributed-observations.png", width = 1200, height = 600, pointsize = 18)
par(mar = c(6.1,4.1,4.1,0.1))
colors = rep("gray", length(CReSIS.dflist))
colors[which(prop.table(sapply(CReSIS.dflist, nrow)) < .024)] <- "indianred2"
barplot(prop.table(sapply(CReSIS.dflist, nrow)), las = 2, col = colors,
        main = "Percentage of contributed observations (CReSIS)")
dev.off()

reduced_CReSIS.dflist <- CReSIS.dflist[which(prop.table(sapply(CReSIS.dflist, nrow)) > .024)]
sum(prop.table(sapply(CReSIS.dflist, nrow))[which(prop.table(sapply(CReSIS.dflist, nrow)) > .024)])

#> Buen ejemplo
write.table(head(reduced_CReSIS.dflist$`2014_P3`), "./LaTeX/Overleaf/CReSIS_expedition.txt")

#> Para cada expedición, transformamos las columnas que nos interesan a vector espacial usando el paquete "terra" (project, vect),
#> omitiendo NAs y teniendo cuidado de construírlos en el sistema de coordenadas (crs) correcto antes de proyectarlos al que usaremos.
CReSIS.vectlist <- lapply(reduced_CReSIS.dflist, function(dataframe){ 
  project(vect(na.omit(dataframe[c("LON", "LAT", "THICK")]), geom = c("LON", "LAT"), crs = "EPSG:4326"), "EPSG:3413") })

#> RAW, with histograms for outliers (over altitude raster) & density representation
my.plot_terra_maps_histograms_boxplots(setNames(lapply(CReSIS.vectlist, function(vect) crop(vect, covariables.rast$Altitude)), names(CReSIS.vectlist)),
                                       "THICK", background.rast = covariables.rast$Altitude, palette = "Temps", ncol = 4, nrow = 5,
                                       consistent_colors = TRUE, global_x = "Ice Thickness", byrow = TRUE, pointsize = 2, orientation = "portrait",
                                       resolution = c(1200, 1600), filename = "./LaTeX/Images/all-raw-CReSIS-maps-histograms-boxplots")

#> Quitamos espesor <= 0
CReSIS_mod.vect <- vect(CReSIS.vectlist)[vect(CReSIS.vectlist)$THICK > 0,]
nrow(CReSIS_mod.vect)

#> Los transformamos a la extensión y resolución acordada mediante la rejilla (Altitude == template.rast)
CReSIS_mean.vect <- {
  rast <- rasterize(CReSIS_mod.vect, covariables.rast$Altitude, field = "THICK", fun = mean)
  vect <- as.points(mask(rast, covariables.rast$Altitude, maskvalues = NA)); names(vect) <- "THICK"; vect
}
nrow(CReSIS_mean.vect)

#> COMPARISON between raw, filtered and processed with histograms
my.plot_terra_maps_histograms_boxplots(list("Original CReSIS" = crop(vect(CReSIS.vectlist), covariables.rast$Altitude),
                                            "Filtered CReSIS" = crop(CReSIS_mod.vect, covariables.rast$Altitude),
                                            "Processed CReSIS" = CReSIS_mean.vect),
                                       "THICK", background.rast = covariables.rast$Altitude, palette = "Temps", consistent_colors = TRUE,
                                       global_x = "Ice Thickness", nrow = 2, orientation = "landscape",
                                       resolution = c(1200, 800), filename = "./LaTeX/Images/comparison-modified-CReSIS-maps-histograms-boxplots")

#> DENSITY in each cell (Altitude ~= template.rast)
my.plot_terra_maps(list("Number of points (CReSIS)" = as.points(mask(rasterize(CReSIS_mod.vect, covariables.rast$Altitude, fun = "count"),
                                                                     covariables.rast$Altitude, maskvalues = NA))),
                   "count", background.rast = covariables.rast$Altitude, palette = "Temps",
                   resolution = c(900, 1200), filename = "./LaTeX/Images/CReSIS-count-map")

#> Creamos las matrices con todos los datos necesarios, divididos aún en expediciones???
CReSIS_data.df <- as.data.frame(extract(covariables.rast, CReSIS_mean.vect,
                                        ID = FALSE, na.rm = TRUE, bind = TRUE, xy = TRUE)) # cell = TRUE
names(CReSIS_data.df)[1] <- "Thickness"
# table(is.na(CReSIS_data.df))
nrow(CReSIS_data.df)

#> Buenos ejemplos
write.table(head(CReSIS_data.df), "./LaTeX/Overleaf/CReSIS_data.txt")

head(CReSIS_data.df)

## ==== Entrenamiento de los modelos ====

data.df <- CReSIS_data.df

#> Distribución en conjuntos de datos de entrenamiento y comprobación
proportion <- c(train = .7, test = .3)
n <- 5e4; n/nrow(data.df)

{set.seed(321)
  datasample <- data.df[sample(1:nrow(data.df), min(n, nrow(data.df)), replace = F),]
  assignation <- sample(cut(
    seq(nrow(datasample)),
    nrow(datasample)*cumsum(c(0, proportion)),
    labels = names(proportion) ))
  split_datasample <- split(datasample, assignation)
  print(addmargins(as.table(sapply(split_datasample, nrow)/nrow(datasample))))
  TrainingData.df <- split_datasample$train
  TestingData.df <- split_datasample$test}
nrow(TrainingData.df)

## ---- Creación de Modelo, Hipótesis ----
#> CROSSPLOT MATRIX
my.plot_crossplot_matrix(cbind(TrainingData.df[1:5], "log(IceVel)" = log(TrainingData.df$IceVel)), palette = "Plasma",
                         resolution = c(1200, 900), filename = "./LaTeX/Images/sample-modified-scatterplot-matrix")

#> PERFORMANCE, log(IceVel)
png("./LaTeX/Images/performance-multiple-linear-regression.png", width = 1200, height = 1200)
performance::check_model(lm(formula = "Thickness ~ Altitude + GravAnom + MassBal + DistCoast + log(IceVel)",
                            data = TrainingData.df), size_title = 24, base_size = 20)
dev.off()
png("./LaTeX/Images/presentation/performance-multiple-linear-regression.png", width = 1200, height = 800)
performance::check_model(lm(formula = "Thickness ~ Altitude + GravAnom + MassBal + DistCoast + log(IceVel)",
                            data = TrainingData.df), check = c("pp_check", "ncv", "homogeneity", "qq"), size_title = 24, base_size = 20)
dev.off()

### ---- Variograma ----

window <- max(TrainingData.df$DistCoast)*1e3/2
nbins <- 20

Vplot <- function(..., resolution = c(800, 600), filename = "Vplot"){
  png(paste0("./LaTeX/Images/", filename, ".png"), width = resolution[1], height = resolution[2])
  print(plot(..., xlab = list("Distance (m)"), ylab = list("Semivariance"),
             par.settings = list(fontsize = list(text = 30, points = 20))))
  dev.off()
}

##### Ordinary Kriging 
formula = "Thickness ~ 1"
cutoff <- window; psill <- NA; range <- 3e4; nugget <- 1e4; vmodel <- "Lin"; anis <- c(0, .8)

#> Dependencia espacial, estacionariedad e isotropía
variogram <- variogram(as.formula(formula), ~x+y, data = TrainingData.df,
                       cutoff = window, width = window/nbins, cressie = FALSE)
Vplot(variogram, main = list("Semivariogram"), filename = "OK-isotropic-semivariogram")

variogram_map <- variogram(gstat(id = "Map", formula = as.formula(formula),
                                 locations = ~x+y, data = TrainingData.df), 
                           map = T, cutoff = window, width = window/nbins, cressie = FALSE)
png("./LaTeX/Images/OK-map-semivariogram.png", width = 800, height = 600)
print(plot(variogram_map, par.settings = list(fontsize = list(text = 30, points = 20))))
dev.off()

#> Ajustar modelo teórico
experimental_variogram <- variogram(as.formula(formula), ~x+y, data = TrainingData.df,
                                    cutoff = cutoff, width = cutoff/nbins, cressie = FALSE)
theoretical_variogram <- fit.variogram(experimental_variogram, debug.level = 2, fit.method = 2,
                                       vgm(psill = psill, vmodel, range = range, nugget = nugget, anis = anis))
print(theoretical_variogram)
Vplot(experimental_variogram, theoretical_variogram, main = list("Adjusted semivariogram"),
      filename = "OK-experimental-theoretical-semivariogram")
anisotropy_variogram <- variogram(as.formula(formula), ~x+y, data = TrainingData.df,
                                  cutoff = cutoff, width = cutoff/nbins,
                                  alpha = rev(seq(0, 180, by = 30))[-1], cressie = FALSE)
Vplot(anisotropy_variogram, theoretical_variogram, filename = "OK-anisotropic-theoretical-semivariogram")

OK.thvar <- theoretical_variogram

##### Universal Kriging
formula = "Thickness ~ Altitude + GravAnom + MassBal + DistCoast + log(IceVel)"
cutoff <- window; psill <- NA; range <- 1e5; nugget <- 1e4; vmodel <- "Exc"; anis <- c(0, 1)

#> Dependencia espacial, estacionariedad e isotropía
variogram <- variogram(as.formula(formula), ~x+y, data = TrainingData.df,
                       cutoff = window, width = window/nbins, cressie = FALSE)
Vplot(variogram, main = list("Semivariogram"), filename = "UK-isotropic-semivariogram")

variogram_map <- variogram(gstat(id = "Map", formula = as.formula(formula),
                                 locations = ~x+y, data = TrainingData.df), 
                           map = T, cutoff = window, width = window/nbins, cressie = FALSE)
png("./LaTeX/Images/UK-map-semivariogram.png", width = 800, height = 600)
print(plot(variogram_map, par.settings = list(fontsize = list(text = 30, points = 20))))
dev.off()

#> Ajustar modelo teórico
experimental_variogram <- variogram(as.formula(formula), ~x+y, data = TrainingData.df,
                                    cutoff = cutoff, width = cutoff/nbins, cressie = FALSE)
theoretical_variogram <- fit.variogram(experimental_variogram, debug.level = 2, fit.method = 2,
                                       vgm(psill = psill, vmodel, range = range, nugget = nugget, anis = anis))
print(theoretical_variogram)
Vplot(experimental_variogram, theoretical_variogram, main = list("Adjusted semivariogram"),
      filename = "UK-experimental-theoretical-semivariogram")
anisotropy_variogram <- variogram(as.formula(formula), ~x+y, data = TrainingData.df,
                                  cutoff = cutoff, width = cutoff/nbins,
                                  alpha = rev(seq(0, 180, by = 30))[-1], cressie = FALSE)
Vplot(anisotropy_variogram, theoretical_variogram, filename = "UK-anisotropic-theoretical-semivariogram")

UK.thvar <- theoretical_variogram


## ==== Predicción ====

#> Construcción de modelos Kriging
nmin <- 1; nmax <- 100; maxdist <- 1e5
OK.model <- gstat(id = "OK", formula = as.formula("Thickness ~ 1"), data = TrainingData.df, locations = ~x+y, 
                  model = OK.thvar, nmin = nmin, nmax = nmax, maxdist = maxdist, set = list(gls = 1))

OK.pred <- interpolate(covariables.rast, OK.model, cores = 4, cpkgs = "gstat", debug.level = 2, na.rm = TRUE, index = 3:4)
my.plot_terra_maps(list("OK ice thickness prediction" = OK.pred[[1]], "Standard deviation of errors" = sqrt(abs(OK.pred[[2]]))),
                   consistent_colors = FALSE, resolution = c(1600, 1200), filename = "./LaTeX/Images/OK-prediction")
OK.pred_pos <- mask(OK.pred, OK.pred[[1]] > 0, maskvalues = FALSE)

nmin <- 1; nmax <- 100; maxdist <- 2e5
UK.model <- gstat(id = "UK", formula = as.formula("Thickness ~ Altitude + GravAnom + MassBal + DistCoast + log(IceVel)"), data = TrainingData.df, locations = ~x+y, 
                  model = UK.thvar, nmin = nmin, nmax = nmax, maxdist = maxdist, set = list(gls = 1))

UK.pred <- interpolate(covariables.rast, UK.model, cores = 4, cpkgs = "gstat", debug.level = 2, na.rm = TRUE, index = 3:4)
my.plot_terra_maps(list("UK ice thickness prediction" = UK.pred[[1]], "Standard deviation of errors" = sqrt(abs(UK.pred[[2]]))),
                   resolution = c(1600, 1200), filename = "./LaTeX/Images/UK-prediction")

#> Eliminamos negativos y con alta variación
sum(as.data.frame(UK.pred[[1]]) <= 0)
UK.pred_pos <- mask(UK.pred, UK.pred[[1]] > 0, maskvalues = FALSE)
maxvar <- max(as.data.frame(OK.pred[[2]])); sqrt(maxvar)
sum(as.data.frame(UK.pred[[2]])$UK.var > maxvar)
UK.pred_pos <- mask(UK.pred_pos, UK.pred[[2]] <= maxvar, maskvalues = FALSE)
my.plot_terra_maps(list("UK ice thickness prediction" = UK.pred_pos[[1]], "Standard deviation of errors" = sqrt(abs(UK.pred_pos[[2]]))),
                   resolution = c(1600, 1200), filename = "./LaTeX/Images/UK-prediction-positive")

### ---- Comparación con otros modelos ----
Methods.rast <- as.numeric(bedmachine.rast$source == 1)
Methods.rast[bedmachine.rast$source == 2] <- 2
Methods.rast[bedmachine.rast$source == 6] <- 3
# my.plot_terra_maps(list("BedMachine Methods" = Methods.rast), resolution = c(800, 1200), filename = "./LaTeX/Images/bedmachine-methods")

#> Introducción de modelos y predicción
MLR.model <- lm(formula = as.formula("Thickness ~ Altitude + GravAnom + MassBal + DistCoast + log(IceVel)"), 
                data = TrainingData.df)
IDW.model <- gstat(id = "IDW", formula = as.formula("Thickness ~ Altitude + GravAnom + MassBal + DistCoast + log(IceVel)"),
                   data = TrainingData.df, locations = ~x+y, nmin = nmin, nmax = nmax, maxdist = maxdist, set = list(gls = 1, idp = 2))

MLR.pred <- predict(covariables.rast, MLR.model, na.rm = TRUE)
MLR.pred_pos <- mask(MLR.pred, MLR.pred > 0, maskvalues = FALSE)
my.plot_terra_maps(list("MLR ice thickness prediction" = MLR.pred_pos), 
                   resolution = c(800, 1200), filename = "./LaTeX/Images/MLR-thickness-positive")
IDW.pred <- interpolate(covariables.rast, IDW.model, debug.level = 2, na.rm = TRUE, index = 3:4)
IDW.pred_pos <- mask(IDW.pred, IDW.pred[[1]] > 0, maskvalues = FALSE)
IDW.pred_pos <- mask(IDW.pred_pos, IDW.pred[[2]] <= maxvar, maskvalues = FALSE)
my.plot_terra_maps(list("IDW ice thickness prediction" = IDW.pred_pos[[1]], "Standard deviation of errors" = sqrt(abs(IDW.pred_pos[[2]]))),
                   resolution = c(1600, 1200), filename = "./LaTeX/Images/IDW-thickness-positive")

bedmachine.rast <- rast("./measurements/BedMachine/BedMachineGreenland-v5.nc")
BM.pred <- trim(mask(resample(bedmachine.rast$thickness, covariables.rast$Altitude, method = "bilinear"), 
                     covariables.rast$Altitude, maskvalues = NA))
BM.pred_pos <- mask(BM.pred, BM.pred > 0, maskvalue = FALSE)
my.plot_terra_maps(list("BedMachine ice thickness prediction" = BM.pred_pos),
                   resolution = c(800, 1200), filename = "./LaTeX/Images/bedmachine-thickness-positive")

# no contenidos en UK + desviacion tipica
exceeded <- BM.pred_pos > UK.pred_pos[[1]] + UK.pred_pos[[2]]
table(as.data.frame(exceeded))

### ---- Compare with using Seventies data ----

temp <- list.files(path = "./measurements/1971-1979/Coord_csvfiles", pattern = "\\.csv$", full.names = TRUE)
Seventies.dflist <- lapply(temp, read.csv)
names(Seventies.dflist) <- gsub("_final_QC.csv","", substring(temp, 41))

write.table(head(Seventies.dflist$`1974_fl10`), "./LaTeX/Overleaf/Seventies_expedition.txt")

Seventies.vectlist <- lapply(Seventies.dflist, function(dataframe){
  dataframe <- transform(dataframe, Thick = srf_elev_m.asl - bed_elev_m.asl);
  if (sum(!is.nan(dataframe$Thick)) != 0){
    project(vect(na.omit(dataframe[c("lon_degrees", "lat_degrees", "Thick")]), geom = c("lon_degrees", "lat_degrees"), crs = "EPSG:4326"), "EPSG:3413")
  }else{NULL}
})
Seventies.vectlist <- Seventies.vectlist[!sapply(Seventies.vectlist, is.null)]  # remove null 1974_fl5

Seventies_outliers <- performance::check_outliers(vect(Seventies.vectlist)$Thick)
Seventies_outliers # tiene 7 valores atípicos
split_index <- split(!Seventies_outliers, unlist(sapply(1:length(Seventies.vectlist), function(k) rep(k, length(Seventies.vectlist[[k]])))))
Seventies.vectlist <- lapply(setNames(1:length(Seventies.vectlist), names(Seventies.vectlist)), function(k) Seventies.vectlist[[k]][split_index[[k]]])

Seventies_mod.vect <- vect(Seventies.vectlist)[vect(Seventies.vectlist)$Thick > 0,]
nrow(Seventies_mod.vect)

Seventies_mean.vect <- {
  rast <- rasterize(Seventies_mod.vect, covariables.rast$Altitude, field = "Thick", fun = mean);
  vect <- as.points(mask(rast, covariables.rast$Altitude, maskvalues = NA)); names(vect) <- "Thick"; vect
}

my.plot_terra_maps_histograms_boxplots(list("Original 1970s" = vect(Seventies.vectlist), "Filtered 1970s" = Seventies_mean.vect), 
                                       "Thick", background.rast = covariables.rast$Altitude, palette = "Temps", 
                                       consistent_colors = FALSE, global_x = "Ice Thickness", byrow = FALSE,
                                       resolution = c(1200, 1200), filename = "./LaTeX/Images/Seventies-comparison-modified-map-histogram-boxplot")

Seventies_data.df <- as.data.frame(extract(covariables.rast, Seventies_mean.vect,
                                           ID = FALSE, na.rm = TRUE, bind = TRUE, xy = TRUE))
names(Seventies_data.df)[1] <- "Thickness"
# table(is.na(Seventies_data.df))
nrow(Seventies_data.df)

write.table(head(Seventies_data.df), "./LaTeX/Overleaf/Seventies_data.txt")

png("./LaTeX/Images/Seventies-performance-multiple-linear-regression.png", width = 1200, height = 1200)
performance::check_model(lm(formula = "Thickness ~ Altitude + GravAnom + MassBal + DistCoast + log(IceVel)",
                            data = Seventies_data.df), size_title = 24, base_size = 20)
dev.off()

### ---- Variograma ----

window <- max(Seventies_data.df$DistCoast)*1e3/2
nbins <- 20

Vplot <- function(..., resolution = c(800, 600), filename = "Vplot"){
  png(paste0("./LaTeX/Images/Seventies-", filename, ".png"), width = resolution[1], height = resolution[2])
  print(plot(..., xlab = list("Distance (m)"), ylab = list("Semivariance"),
             par.settings = list(fontsize = list(text = 30, points = 20))))
  dev.off()
}


##### Universal Kriging on 197X
formula = "Thickness ~ Altitude + GravAnom + MassBal + DistCoast + log(IceVel)"
cutoff <- window; psill <- NA; range <- 1e5; nugget <- 1e4; vmodel <- "Exc"; anis <- c(0, 1)

#> Dependencia espacial, estacionariedad e isotropía
# variogram <- variogram(as.formula(formula), ~x+y, data = Seventies_data.df,
#                        cutoff = window, width = window/nbins, cressie = FALSE)
# Vplot(variogram, main = list("Semivariograma"), filename = "UK-isotropic-semivariogram")

variogram_map <- variogram(gstat(id = "Map", formula = as.formula(formula),
                                 locations = ~x+y, data = Seventies_data.df), 
                           map = T, cutoff = window, width = window/nbins, cressie = FALSE)
png("./LaTeX/Images/Seventies-UK-map-semivariogram.png", width = 800, height = 600)
print(plot(variogram_map, par.settings = list(fontsize = list(text = 30, points = 20))))
dev.off()

#> Ajustar modelo teórico
experimental_variogram <- variogram(as.formula(formula), ~x+y, data = Seventies_data.df,
                                    cutoff = cutoff, width = cutoff/nbins, cressie = FALSE)
theoretical_variogram <- fit.variogram(experimental_variogram, debug.level = 2, fit.method = 2,
                                       vgm(psill = psill, vmodel, range = range, nugget = nugget, anis = anis))
print(theoretical_variogram)
# Vplot(experimental_variogram, theoretical_variogram, main = list("Semivariograma ajustado"),
#       filename = "UK-experimental-theoretical-semivariogram")
anisotropy_variogram <- variogram(as.formula(formula), ~x+y, data = Seventies_data.df,
                                  cutoff = cutoff, width = cutoff/nbins,
                                  alpha = rev(seq(0, 180, by = 30))[-1], cressie = FALSE)
Vplot(anisotropy_variogram, theoretical_variogram, filename = "UK-anisotropic-theoretical-semivariogram")

UK_197X.thvar <- theoretical_variogram

nmin <- 1; nmax <- 100; maxdist <- Inf
UK_197X.model <- gstat(id = "UK_197X", formula = as.formula("Thickness ~ Altitude + GravAnom + MassBal + DistCoast + log(IceVel)"), data = Seventies_data.df, locations = ~x+y, 
                       model = UK_197X.thvar, nmin = nmin, nmax = nmax, maxdist = maxdist, set = list(gls = 1))

UK_197X.pred <- interpolate(covariables.rast, UK_197X.model, cores = 4, cpkgs = "gstat", debug.level = 2, na.rm = TRUE, index = 3:4)
my.plot_terra_maps(list("UK-197X ice thickness prediction" = UK_197X.pred[[1]], "Standard deviation of errors" = sqrt(abs(UK_197X.pred[[2]]))),
                   resolution = c(1600, 1200), filename = "./LaTeX/Images/Seventies-UK-prediction")

#> Eliminamos negativos y con alta variación
sum(as.data.frame(UK_197X.pred[[1]]) <= 0)
UK_197X.pred_pos <- mask(UK_197X.pred, UK_197X.pred[[1]] > 0, maskvalues = FALSE)
maxvar <- max(as.data.frame(OK.pred[[2]])); sqrt(maxvar)
sum(as.data.frame(UK_197X.pred[[2]])$UK.var > maxvar)
UK_197X.pred_pos <- mask(UK_197X.pred_pos, UK_197X.pred[[2]] <= maxvar, maskvalues = FALSE)
my.plot_terra_maps(list("UK-197X ice thickness prediction" = UK_197X.pred_pos[[1]], "Standard deviation of errors" = sqrt(abs(UK_197X.pred_pos[[2]]))),
                   resolution = c(1600, 1200), filename = "./LaTeX/Images/Seventies-UK-prediction-positive")

#> Comparison functions
RMSE <- function(prediction, measurements = TestingData.df$Thickness){
  return(sqrt(mean((prediction - measurements)^2)))
}
R2 <- function(prediction, measurements = TestingData.df$Thickness){
  RSS <- sum((prediction - measurements)^2)
  TSS <- sum((measurements - mean(measurements))^2)
  
  return(1 - RSS/TSS)
}
SLEandArea <- function(thickness.rast){
  require(terra)
  
  thickness_pos.rast <- mask(thickness.rast, thickness.rast > 0, maskvalue = FALSE)
  
  area.km2 <- prod(res(thickness.rast)/1e3)
  extent.percent <- (area.km2 * global(!is.na(thickness_pos.rast), sum)) / 
    (area.km2 * global(!is.na(thickness.rast), sum)) * 100
  
  volume.km3 <- as.numeric(global(area.km2 * thickness_pos.rast/1e3, sum, na.rm = TRUE))
  print("Volume (km3):"); print(volume.km3)
  Gt <- volume.km3 * 0.9167 # 0.9167 density of ice
  print("Gt:"); print(Gt)
  
  SLE.m <- Gt * 1/361.8 # 361.8 Gt ice = 1 mm sea level elevation
  
  return(list("SLE" = SLE.m/1e3, "Extent" = extent.percent))
}

#> Evaluate test results for all models
TestingData.vect <- vect(TestingData.df[c("Thickness", "x", "y")], geom = c("x", "y"), crs = crs(covariables.rast$Altitude))
models <- c("UK_197X.pred", "MLR.pred", "IDW.pred", "OK.pred", "UK.pred", "BM.pred")
Test_results.df <- cbind(sapply(mget(models), function(rast) extract(rast[[1]], TestingData.vect, ID = FALSE)[,1] ))

Metrics <- cbind(apply(Test_results.df, 2, RMSE), apply(Test_results.df, 2, R2))
colnames(Metrics) <- c("RMSE", "R2"); rownames(Metrics) <- gsub(".pred", "", models)
Metrics

Calcs <- sapply(models, function(model) unlist(SLEandArea(get(model)[[1]])))
rownames(Calcs) <- c("SLE", "Extent"); colnames(Calcs) <- gsub(".pred", "", models)
t(Calcs)

# SLEandArea(bedmachine.rast$thickness)

#> Plot all models for visual comparison
my.plot_terra_maps(list(
  "Predicción UK-197X del espesor" = UK_197X.pred_pos[[1]],
  "Predicción MLR del espesor" = MLR.pred_pos[[1]],
  "Predicción IDW del espesor" = IDW.pred_pos[[1]]
  ,
  "Predicción OK del espesor" = OK.pred_pos[[1]],
  "Predicción UK del espesor" = UK.pred_pos[[1]],
  "Predicción BedMachine del espesor" = BM.pred_pos[[1]]
), byrow = TRUE, ncol = 3,
resolution = c(1600, 1600), filename = "./LaTeX/Images/all-model-predictions")

### ----  Anexo I ----

# source("my_plotting_functions.R")
# my.plot_terra_maps(list("BedMachine Mask" = bedmachine.rast$mask), resolution = c(800, 1200), filename = "./LaTeX/Images/bedmachine-mask")


## ==== Examples =====
library(sp)
data(meuse)
coordinates(meuse) = ~x+y
v <- variogram(log(zinc) ~ 1, meuse)
v.fit <- fit.variogram(v, vgm(1, "Sph", 700, 1))
v.fit
png("./LaTeX/Images/meuse-experimental-variogram-example.png", width = 1200, height = 600)
plot(v, ylab = list("Semivariance", cex = 1.2), xlab = list("Distance (m)", cex = 1.2), par.settings = list(fontsize = list(text = 25, points = 15)))
dev.off()
png("./LaTeX/Images/meuse-variogram-example.png", width = 1200, height = 600)
plot(v, v.fit, ylab = list("Semivariance", cex = 1.2), xlab = list("Distance (m)", cex = 1.2), par.settings = list(fontsize = list(text = 25, points = 15)))
library(lattice)
trellis.focus("panel",1,1)
panel.arrows(0, sum(v.fit[,2]), v.fit[2,3], sum(v.fit[,2]), lwd = 3, length = .2, angle = 20)
panel.arrows(v.fit[2,3], sum(v.fit[,2]), 0, sum(v.fit[,2]), lwd = 3, length = .2, angle = 20)
ltext(x=v.fit[2,3]/2, y=1.04*sum(v.fit[,2]), "Range", cex=1)
panel.arrows(v.fit[2,3], sum(v.fit[,2]), v.fit[2,3], 0, lwd = 3, length = .2, angle = 20)
panel.arrows(v.fit[2,3], 0, v.fit[2,3], sum(v.fit[,2]), lwd = 3, length = .2, angle = 20)
ltext(x=v.fit[2,3]+10, y=sum(v.fit[,2])/2, "Sill", cex=1, adj=0)
llines(x=c(0,30), y=v.fit[1,2], col = "black", lwd=1, lty=3)
panel.arrows(30, 0, 30, v.fit[1,2], lwd = 3, length = .15, angle = 20)
panel.arrows(30, v.fit[1,2], 30, 0, lwd = 3, length = .15, angle = 20)
ltext(x=40, y=v.fit[1,2]/1.8, "Nugget", cex=1, adj=0)
trellis.unfocus()
dev.off()

# show.vgms()
png("./LaTeX/Images/theoretical-variograms-example.png", width = 1200, height = 500)
show.vgms(models = c("Nug", "Lin", "Sph", "Exp"), max = 1,
          xlab = list("Distance (m)"), ylab = list("Semivariance"),
          par.settings = list(fontsize = list(text = 30, points = 20)))
dev.off()


# Sys.setenv(JAVA_HOME = "C:/Program Files/Java/jre1.8.0_431/")
# library(glmulti)


v <- vect(CReSIS.vectlist)
r <- covariables.rast$Altitude

plot(covariables.rast$Altitude)
z <- zoom(covariables.rast$Altitude)
window(r) <- z

png("./LaTeX/Images/presentation/non-adjusted-measurements.png", height = 800, width = 800, pointsize = 18)
plot(r, col = gray(seq(0, .8, by = .01)), legend = F)
plot(crop(v, z), "THICK", type = "continuous", col = c("white", hcl.colors(256, "Temps")), add = T)
dev.off()

v2 <- CReSIS_mean.vect
png("./LaTeX/Images/presentation/adjusted-measurements.png", height = 800, width = 800, pointsize = 18)
plot(r, col = gray(seq(0, .8, by = .01)), legend = F)
plot(crop(v2, z), "THICK", type = "continuous", col = hcl.colors(256, "Temps"), add = T)
dev.off()


######### cementery

BMkriging.rast <- trim(mask(resample(Methods.rast, covariables.rast$Altitude, method = "bilinear"), 
                            covariables.rast$Altitude, maskvalues = NA))

models <- c("UK_197X.pred", "MLR.pred", "IDW.pred", "OK.pred", "UK.pred", "BM.pred")
Test_results.df <- cbind(sapply(mget(models), function(rast) extract(mask(rast[[1]], BMkriging.rast == 3, maskvalues = FALSE), TestingData.vect, ID = FALSE)[,1] ))

Metrics <- cbind(apply(na.omit(Test_results.df), 2, function(x) RMSE(x, measurements = TestingData.df$Thickness[!is.na(Test_results.df[,1])])), 
                 apply(na.omit(Test_results.df), 2, function(x) R2(x, measurements = TestingData.df$Thickness[!is.na(Test_results.df[,1])])))
colnames(Metrics) <- c("RMSE", "R2"); rownames(Metrics) <- gsub(".pred", "", models)
Metrics

Calcs <- sapply(models, function(model) unlist(SLEandArea(get(model)[[1]])))
rownames(Calcs) <- c("SLE", "Extent"); colnames(Calcs) <- gsub(".pred", "", models)
t(Calcs)


my.plot_terra_maps_histograms_boxplots(list("Datos 1971-1979" = vect(Seventies_data.df, geom = c("x", "y"), crs = crs(covariables.rast$Altitude))), 
                                       "Thickness", background.rast = covariables.rast$Altitude, palette = "Temps", 
                                       consistent_colors = FALSE, global_x = "Thickness del hielo", byrow = TRUE,
                                       resolution = c(1200, 600), filename = "./LaTeX/Images/presentation/Seventies-modified-map-histogram-boxplot")

################################

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
