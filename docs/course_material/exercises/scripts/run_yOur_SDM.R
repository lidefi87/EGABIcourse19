library(ncdf4)
library(raster)
library(dismo)
library(gbm)

#### Abre los records de ocurrencia ####
#---------------------------------------
data_folder <- "docs/course_material/exercises/data"
scripts_folder <- "docs/course_material/exercises/scripts"
occ.sterechinus <- read.csv(file.path(data_folder, "occurrences_sterechinus.csv"), 
                            header=T, sep=";")
# head(occ.sterechinus)

#### Abre las capas con información ambiental y únelas en un solo raster ####
#-----------------------------------------------------------------------
depth <- raster(file.path(data_folder, "environmental_layers/depth.nc"))
sediments <- raster(file.path(data_folder, "environmental_layers/sediments.nc"))
seafloor_temp_2005_2012_max <- raster(file.path(data_folder, "environmental_layers/seafloor_temp_2005_2012_max.nc"))
POC_2005_2012_max <- raster(file.path(data_folder, "environmental_layers/POC_2005_2012_max.nc"))
seafloor_current_speed <- raster(file.path(data_folder, "environmental_layers/seafloor_current_speed.nc"))

predictors_stack <- brick(depth,sediments,seafloor_temp_2005_2012_max,POC_2005_2012_max,seafloor_current_speed)

## miremos a las propiedades de las variables ambientales
#.......................................................
#predictors_stack
#plot(predictors_stack)

# miremos a la distribucion de las ocurrencias
#.............................................
#plot(depth)
#points(occ.sterechinus[,c(2,1)], pch=20) # longitude first, latitude second

#--------------------------------------------------------------------------------------------------------------
# Abramos la capar KDE con el esfuerzo de muestreo. Esta es la capa que usaremos para crear puntos background 
# (con ponderacion) 
# El KDE (Estimación de Densidad Kernel) es una herramienta estadística que ayuda a medir la probabilidad de
# encontrar una ocurrencia en un pixel basada en los records del conjunto de datos bénticos del Océano del Sur
# (datos del Atlas Biogeográfico del Océano del Sur)
#--------------------------------------------------------------------------------------------------------------
KDE <- raster(file.path(data_folder, "KDE.asc"))
# en esta capa, los continentes están definedos por pixeles con valor NA
# esto permite a R reconocer las áreas donde los puntos background no deberían
# ser muestreados
#-----------------------------------------------------------------------------------------

#### Iniciando
#-------------
cv.boot <- 2 # número de réplicas (# de veces que el modelo es entrenado)
source(file.path(scripts_folder, "Function_gbm.R"))

#### Alineando las matrices vacías
#---------------------------------

stack.pred<-subset(predictors_stack,1);values(stack.pred)<-NA
#testvaluesp<-rep(NA,nrow(fichier_data)) ; testvaluesa <- testvaluesp

model_stats <- matrix(NA, 6, cv.boot*4,dimnames = list(c("AUC", "COR", "TSS", "maxSSS", "valid_test_data","prop_test"), NULL))

# Guardando la contribución
contTr <- matrix(NA, dim(predictors_stack)[3], cv.boot, dimnames = list(names(predictors_stack), NULL))
n.seed <- seq(1,60,1)[-c(4,5,17,28,32,41,45,46,48,53)] # controla y fija el muestreo al azar (así podemos comparar luego)

for (j in 1:cv.boot){

  #-----------------------------------------
  #### crea la matriz de ocurrencia-ambiente 
  #-----------------------------------------
  envi.presences <- unique(extract (predictors_stack,occ.sterechinus[,c(2,1)]))
  presence.data <- occ.sterechinus[-which(duplicated(extract(predictors_stack,occ.sterechinus[,c(2,1)]))),c(2,1)]
  colnames(presence.data)<- c("longitude","latitude")
  # la función 'unique' permite remover duplicados en el conjunto de datos (ocurrencias que ocurren en un mismo pixel)
  # la funcion 'duplicated' permite identificar las filas con los mismos datos
  #head(envi.presences)
  # los datos relacionados a presencias estarán asociadas con ID=1
  
  set.seed(n.seed[j])
  # muestreo de datos background: en el bucle, cambios con cada réplica son aplicados
  # 1000 puntos background son muestreados al azar en el ambiente de acuerdo a la ponderación con la capa KDE
  background_data <- xyFromCell(KDE, sample(which(!is.na(values(KDE))), 200, prob=values(KDE)[!is.na(values(KDE))]))

  colnames(background_data) <- colnames(presence.data)
  # extrayendo condiciones ambientales en los puntos background
  envi.background <- extract(predictors_stack,background_data)
  # los puntos background serán asociados con ID=0

  # Inicializar la matriz que contiene datos de presencia, background y de condiciones ambientales
  id<-0;sdmdata.unique<-0;  id<-c(rep(1,nrow(envi.presences)),rep(0,nrow(envi.background))) 
  MATRIX_OCC_ENVI<-data.frame(cbind(id,rbind(envi.presences,envi.background)))
  #head(MATRIX_OCC_ENVI)

  # Dividir los datos de ocurrencia-background en pliegues de datos para prueba y entrenamiento con segregación espacial
  dat1 <- rbind(cbind(background_data, Isp=rep(0,nrow(background_data))), 
                cbind(presence.data,Isp=rep(1,nrow(presence.data))))
  colnames(dat1)<- c("longitude","latitude","Isp")
  #tail(dat1)
  idP <- which(dat1$Isp == 1) #  id de datos de presencia para dividirlos
  MyFold <- rep(NA, nrow(dat1)) # variable vacía para guardar los grupos del conjunto de datos (1 al 4)
  
  source(file.path(scripts_folder, "clock4_crossValidation.R"))
  clock4F <- clock4(dat1[idP, c("longitude", "latitude")], dat1[-idP, c("longitude", "latitude")])
  
  # Extrayendo los pliegues
  MyFold[idP] <- clock4F$occ.grp
  MyFold[-idP] <- clock4F$bg.coords.grp
  plot(dat1[,c("longitude", "latitude")], pch = 20, col = c("red", "blue","black","purple")[as.factor(MyFold)])
  
  #------------------------
  #### Entrenando el modelo 
  #------------------------
  model.res<- gbm.step_v2 (data=MATRIX_OCC_ENVI, 
                         gbm.x = 2:ncol(MATRIX_OCC_ENVI),
                         gbm.y = 1,
                         family = "bernoulli",
                         n.folds=4,
                         fold.vector = MyFold, 
                         tree.complexity = 3,
                         learning.rate = 0.015,
                         bag.fraction =0.5)

  #--------------------------------------------
  #### Extrayendo datos y resultados del modelo
  #--------------------------------------------
  # Predicciones 
  p<-predict(predictors_stack,model.res,n.trees=model.res$gbm.call$best.trees,type="response", na.rm=F)
  stack.pred<-stack(stack.pred,p) # apilamos todos los mapas
  
  ########
  ## CV ## (= pliegues de validación cruzada o CV)
  ########
    j_cv <- ((j-1) * 4+1):(j*4) # para contar los pliegues
  
  model_stats["AUC", j_cv] <- model.res$cv.roc.matrix
  model_stats["COR", j_cv] <- model.res$cv.cor.matrix
  model_stats["TSS", j_cv] <- model.res$tss.cv
  model_stats["maxSSS", j_cv] <- model.res$cv.th.matrix

  ## clasificando los datos de prueba correctamente
  model_stats["valid_test_data", j_cv] <- model.res$cv.corr.class*100
  model_stats["prop_test", j_cv] <- model.res$cv.length*100
  
  # Obtener contribuciones
  RI <- summary(model.res, plotit = F) # extrayendo las contribuciones
  contTr[match(RI$var, rownames(contTr)), j] <- RI[,"rel.inf"]

}

#------------------------
#### Mapas de predicción 
#------------------------
mean_stack <- raster::calc(stack.pred, mean, na.rm=T); mean_stack <- mask(mean_stack, depth)
#sd_stack <- raster::calc(stack.pred,sd, na.rm=T); sd_stack <- mask(sd_stack, depth)

# graficamos los resultados 
#continent <- read.csv("data/worldmap.csv") # añadimos a los continentes
#plot(mean_stack) ; points(continent, type="l")

# este es un mapa aproximado, si quisieras un mejor map, puedes exportar a un documento ascii 
# y abrirlo en otro software como Qgis o similar
# writeRaster(mean_stack, "results/mean_raster.asc")
# writeRaster(sd_stack, "results/sd_raster.asc")

#---------------------------
#### Estadísticas del modelo 
#---------------------------
ecM <- apply(model_stats, 1, mean, na.rm=T)
ecSD <- apply(model_stats, 1, sd, na.rm=T)
ecTot <- paste(round(ecM, 3), round(ecSD, 3), sep = " ± ")
names(ecTot) <- names(ecM)

ResF <- data.frame(c(ecTot["AUC"],ecTot["COR"],ecTot["TSS"],ecTot["maxSSS"], ecTot["valid_test_data"], ecTot["prop_test"]))
rownames(ResF) <- c("AUC","COR","TSS", "maxSSS",  "Correctly classified test data", "Test data (% of total dataset)")
colnames(ResF) <- "Estadísticas promedio del modelo"

# matriz con resultados en crudo
#model_stats
# valores promedio
#ResF

# ahora podemos exportar los resultados
#write.csv(model_stats,"results/model_stats.csv"))

## Contribución de variables ambientales
# raw data 
#contTr

# calculamos el promedio
CtM <- apply(contTr, 1, mean)
CtSD <- apply(contTr, 1, sd)
CtTot <- paste(round(CtM, 3), round(CtSD, 3), sep = " ± ")
names(CtTot) <- names(CtM)

CtTot <- data.frame(CtTot) ; colnames(CtTot) <- "Contribución de las variables ambientales al modelo (%)"
#write.csv(CtTot, "results/avg_contribution.csv")


####---------------------------------------------------------------------------
####---------------------------------------------------------------------------
### CALCULANDO EXTRAPOLACIÓN
# Superficie de Similaridad Ambinetal Multivariada (Elith et al. 2010) 
# envi.presences <- unique(extract (predictors_stack,occ.sterechinus[,c(2,1)]))
# x <- dismo::mess(predictors_stack, na.omit(envi.presences))
# 
# y <- x; values(y)<- values(x)>0  # mira a Elith et al. (2010): 
#cuando el valor calculado de MESS es negativo significa que hubo una extrapolación (datos predichos fuera 
#del rango de entrenamiento)
# y <- reclassify(y,cbind(FALSE,0)) # area de extrapolación
# y <- reclassify(y,cbind(TRUE,1))  # no hubo extrapolación, dentro de los rangos de calibración
# 
# plot(y)

####---------------------------------------------------------------------------
####---------------------------------------------------------------------------

