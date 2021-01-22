setwd("D:/Paper Thesis UQ/Statistical Analysis/Crown-based/fd_cleared")
#load package
library(raster)
library(MASS)
library(randomForest)
library(glmnet)
library(tidyverse)
library(broom)
library(glmnet)
library(caret)
library(e1071)
library(rsample)
library(Boruta)
library(snow)
library(doParallel)
library(rgdal)
library(rsample)
rasterOptions(progress = "Text", tmpdir = paste0(getwd(), "/tmp"))

#load data
setwd("H:/Paper Thesis UQ/Hyperspectral/Hawk/Original")
eagle.agg <- stack("eagle_mosaic_5m_spectrum.tif")
hawk.agg <- stack("hawk_mosaic_5m_spectrum.tif")
setwd("H:/Paper Thesis UQ/Lidar data/Merged")
chm <- raster("chm_pitfree_lidar_2m_des2019.tif")
setwd("H:/Paper Thesis UQ/Lidar data/Raster_Metrics")
lidar <- stack("2m_lidar_height_metrics_clip.tif")
grid <- raster("grid_100_sel.tif")
setwd("H:/Paper Thesis UQ/Statistical Analysis/Raster-based/field data (raster)")
fd.5m <- raster("04_FD_ht_max_5m.tif")

#PCA predictor
setwd("H:/Paper Thesis UQ/Statistical Analysis/Raster-based/Predictors")
blue.pca <- stack("00_Blue_PCA_Eagle_Spectrum.tif")
green.pca <- stack("01_Green_PCA_Eagle_Spectrum.tif")
re.pca <- stack("02_RE_PCA_Eagle_Spectrum.tif")
nir.pca <- stack("03_NIR2_PCA_EagleHawk_Spectrum.tif")
swir.pca <- stack("04_SWIR_PCA_Hawk_Spectrum.tif")
lidar.pca <- stack("05_HeightM_PCA_Lidar_Spectrum.tif")

blue.res <- resample(blue.pca, lidar.pca, "bilinear")
green.res <- resample(green.pca, lidar.pca, "bilinear")
re.res <- resample(re.pca, lidar.pca, "bilinear")
nir.res <- resample(nir.pca, lidar.pca, "bilinear")
swir.res <- resample(swir.pca, lidar.pca, "bilinear")
pca.bands <- stack(blue.res, green.res, re.res, nir.res, swir.res, lidar.pca)

#quantile fd.40m <- aggregate(fd.10m, fact = 4, fun = function(i, na.rm) quantile(i, probs=0.5, na.rm=na.rm))
#data resampling
fd.50m <- aggregate(fd.5m, fact = 10, fun = mean)
cl <- makeCluster(2)
beginCluster(cl)
chm.res <- resample(chm, fd.50m, "bilinear")
lidar.res <- resample(lidar, fd.50m, "bilinear")
eagle.res <- resample(eagle.agg, fd.50m, "bilinear")
hawk.res <- resample(hawk.agg, fd.50m, "bilinear")
pca.res <- resample(pca.bands, fd.50m, "bilinear")
all.pred <- stack(eagle.res, hawk.res, lidar.res, chm.res, pca.res)
endCluster()

#creating sample training and test
predv <- values(all.pred)
fdv <- values(fd.50m)
comb <- cbind(fdv, predv)
combc <- comb[complete.cases(comb),]
combc <- as.data.frame(combc)
combc_split <- initial_split(combc, prop = .7)
train <- training(combc_split)
train.ori <- data.frame(cbind(train$fdv, train[,2:517]))
names(train.ori)[1] <- c("fdv")
train.pca <- data.frame(cbind(train$fdv, train[,518:535]))
names(train.pca)[1] <- c("fdv")
test  <- testing(combc_split)
test.ori <- data.frame(cbind(test$fdv, test[,2:517]))
names(test.ori)[1] <- c("fdv")
test.pca <- data.frame(cbind(test$fdv, test[,518:535]))
names(test.pca)[1] <- c("fdv")

#for boruta
predv2 <- values(all.pred[[1:516]])
fdv <- values(fd.50m)
comb2 <- cbind(fdv, predv2)
combc2 <- comb2[complete.cases(comb2),]
combc2 <- as.data.frame(combc2)
names(combc2)[1] <- c("fdv")

#variable selection
boruta.fs <- Boruta(fdv~., data = combc2, pValue = 0.01, maxRun = 200, doTrace = 2, ntree = 1000)
boruta.fs <- TentativeRoughFix(boruta.fs)
#sa_ctrl <- safsControl(functions = rfSA, method = "repeatedcv", repeats = 3, improve = 5)
#sa_obj <- safs(x=combc[,2:503], y=combc[,1], safsControl = sa_ctrl)

#Tuning parameter
#random search ori
rand_ctrl <- trainControl(method = "repeatedcv", repeats = 5, search = "random")
rand_search.xgb <- train(getConfirmedFormula(boruta.fs), data = train.ori, method = "xgbTree",
                         tuneLength = 50, metric = "Rsquared", trControl = rand_ctrl)
rand_search.xgb.DART <- train(getConfirmedFormula(boruta.fs), data = train.ori, method = "xgbDART", 
                          tuneLength = 50, metric = "Rsquared", trControl = rand_ctrl)
rand_search.xgb.Linear <- train(getConfirmedFormula(boruta.fs), data = train.ori, method = "xgbLinear", 
                              tuneLength = 50, metric = "Rsquared", trControl = rand_ctrl)
rand_search.earth <- train(getConfirmedFormula(boruta.fs), data = train.ori, method = "earth",
                           tuneLength = 50, metric = "Rsquared", trControl = rand_ctrl)
rand_search.rf <- train(getConfirmedFormula(boruta.fs), data = train.ori, method = "extraTrees",
                        tuneLength = 50, metric = "Rsquared", trControl = rand_ctrl)
rand_search.svm <- train(getConfirmedFormula(boruta.fs), data = train.ori, method = "svmRadial",
                         tuneLength = 50, metric = "Rsquared", preProc = c("center", "scale"), trControl = rand_ctrl)
rand_search.earth.lidar <- train(getConfirmedFormula(boruta.fs.lidar), data = train.ori.lidar, method = "earth",
                           tuneLength = 50, metric = "Rsquared", trControl = rand_ctrl)
rand_search.rf.lidar <- train(getConfirmedFormula(boruta.fs.lidar), data = train.ori.lidar, method = "extraTrees",
                                 tuneLength = 50, metric = "Rsquared", trControl = rand_ctrl)
rand_search.earth.hs <- train(getConfirmedFormula(boruta.fs.hs), data = train.ori.hs, method = "earth",
                             tuneLength = 50, metric = "Rsquared", trControl = rand_ctrl)


#random search pca
rand_ctrl <- trainControl(method = "repeatedcv", repeats = 5, search = "random")
rand_search.xgb.pca <- train(train.fdv~., data = train.pca, method = "xgbTree",
                         tuneLength = 50, metric = "Rsquared", trControl = rand_ctrl)
rand_search.earth.pca <- train(train.fdv~., data = train.pca, method = "earth",
                           tuneLength = 50, metric = "Rsquared", trControl = rand_ctrl)
rand_search.rf.pca <- train(train.fdv~., data = train.pca, method = "extraTrees",
                        tuneLength = 50, metric = "Rsquared", trControl = rand_ctrl)
rand_search.svm.pca <- train(train.fdv~., data = train.pca, method = "svmRadial",
                         tuneLength = 50, metric = "Rsquared", preProc = c("center", "scale"), trControl = rand_ctrl)

#bayesian optimization
#see tuning text files for each algorithms

#validation for test data
#final_search.xgb.train <- train(getConfirmedFormula(boruta.fs), data = train, method = "xgbTree", 
                                #tuneGrid = data.frame(nrounds = ba_search$Best_Par["nrounds"], max_depth = ba_search$Best_Par["max_depth"], eta = ba_search$Best_Par["eta"], gamma = ba_search$Best_Par["gamma"], colsample_bytree = ba_search$Best_Par["colsample_bytree"], min_child_weight = ba_search$Best_Par["min_child_weight"], subsample = ba_search$Best_Par["subsample"]), 
#metric = "Rsquared",  trControl = ctrl)
#final_search.xgb.test <- train(getConfirmedFormula(boruta.fs), data = test, method = "xgbTree", 
#tuneGrid = data.frame(nrounds = ba_search$Best_Par["nrounds"], max_depth = ba_search$Best_Par["max_depth"], eta = ba_search$Best_Par["eta"], gamma = ba_search$Best_Par["gamma"], colsample_bytree = ba_search$Best_Par["colsample_bytree"], min_child_weight = ba_search$Best_Par["min_child_weight"], subsample = ba_search$Best_Par["subsample"]), 
#metric = "Rsquared",  trControl = ctrl)
#final_search.earth.train <- train(getConfirmedFormula(boruta.fs), data = train, method = "earth", tuneGrid = data.frame(degree = ba_search$Best_Par["degree"], nprune = ba_search$Best_Par["nprune"]),
#metric = "Rsquared",  trControl = ctrl)
#final_search.earth.test <- train(getConfirmedFormula(boruta.fs), data = test, method = "earth", tuneGrid = data.frame(degree = ba_search$Best_Par["degree"], nprune = ba_search$Best_Par["nprune"]),
#metric = "Rsquared",  trControl = ctrl)
#final_search.rf.train <- train(getConfirmedFormula(boruta.fs), data = train, method = "extraTrees", 
#tuneGrid = data.frame(mtry = ba_search$Best_Par["mtry"], numRandomCuts = ba_search$Best_Par["numRandomCuts"]),
#metric = "Rsquared",  trControl = ctrl)
#final_search.rf.test <- train(getConfirmedFormula(boruta.fs), data = test, method = "extraTrees", 
#tuneGrid = data.frame(mtry = ba_search$Best_Par["mtry"], numRandomCuts = ba_search$Best_Par["numRandomCuts"]),
#                              metric = "Rsquared",  trControl = ctrl)
#final_search.svm.train <- train(getConfirmedFormula(boruta.fs), data = train, method = "svmRadial", tuneGrid = data.frame(sigma = ba_search$Best_Par["sigma"], C = ba_search$Best_Par["c"]),
#                               metric = "Rsquared", preProc = c("center", "scale"), trControl = ctrl)
#final_search.svm.test <- train(getConfirmedFormula(boruta.fs), data = test, method = "svmRadial", tuneGrid = data.frame(sigma = ba_search$Best_Par["sigma"], C = ba_search$Best_Par["c"]),
#                           metric = "Rsquared", preProc = c("center", "scale"), trControl = ctrl)
final_search.earth.train <- train(getConfirmedFormula(boruta.fs), data = train.ori, method = "earth", 
                                  tuneGrid = rand_search.earth$bestTune, 
                                  metric = "Rsquared",  trControl = rand_ctrl)
final_search.earth.train.lidar <- train(getConfirmedFormula(boruta.fs.lidar), data = train.ori.lidar, method = "earth", 
                                  tuneGrid = rand_search.earth.lidar$bestTune, 
                                  metric = "Rsquared",  trControl = rand_ctrl)
final_search.rf.train.lidar <- train(getConfirmedFormula(boruta.fs.lidar), data = train.ori.lidar, method = "extraTrees", 
                                        tuneGrid = rand_search.rf.lidar$bestTune, 
                                        metric = "Rsquared",  trControl = rand_ctrl)
final_search.earth.train.hs <- train(getConfirmedFormula(boruta.fs.hs), data = train.ori.hs, method = "earth", 
                                        tuneGrid = rand_search.earth.hs$bestTune, 
                                        metric = "Rsquared",  trControl = rand_ctrl)
final_search.rf.train <- train(getConfirmedFormula(boruta.fs), data = train.ori, method = "extraTrees", 
                               tuneGrid = rand_search.rf$bestTune, 
                               metric = "Rsquared",  trControl = rand_ctrl)
final_search.xgb.train <- train(getConfirmedFormula(boruta.fs), data = train.ori, method = "xgbTree", 
                                tuneGrid = rand_search.xgb$bestTune, 
                                metric = "Rsquared",  trControl = rand_ctrl)
final_search.svm.train <- train(getConfirmedFormula(boruta.fs), data = train.ori, method = "svmRadial", 
                                tuneGrid = rand_search.svm$bestTune, preProc = c("center", "scale"),
                                metric = "Rsquared",  trControl = rand_ctrl)

#pca test and validation
final_search.earth.train.pca <- train(fdv~., data = train.pca, method = "earth", 
                                      tuneGrid = rand_search.earth.pca$bestTune, 
                                      metric = "Rsquared",  trControl = rand_ctrl)
final_search.rf.train.pca <- train(fdv~., data = train.pca, method = "extraTrees", tuneGrid = rand_search.rf.pca$bestTune, metric = "Rsquared",  trControl = rand_ctrl)
final_search.xgb.train.pca <- train(fdv~., data = train.pca, method = "xgbTree", 
                                    tuneGrid = rand_search.xgb.pca$bestTune, 
                                    metric = "Rsquared",  trControl = rand_ctrl)
final_search.xgb.Linear.train.pca <- train(fdv~., data = train.pca, method = "xgbLinear", 
                                    tuneGrid = rand_search.xgb.Linear.pca$bestTune, 
                                    metric = "Rsquared",  trControl = rand_ctrl)
final_search.svm.train.pca <- train(fdv~., data = train.pca, method = "svmRadial", tuneGrid = rand_search.svm.pca$bestTune, preProc = c("center", "scale"), metric = "Rsquared",  trControl = rand_ctrl)




############validation gak jadi dipake
final_search.earth.test <- train(getConfirmedFormula(boruta.fs), data = test.ori, method = "earth", 
                                 tuneGrid = rand_search.earth$bestTune, 
                                 metric = "Rsquared",  trControl = rand_ctrl)
final_search.rf.test <- train(getConfirmedFormula(boruta.fs), data = test.ori, method = "extraTrees", 
                              tuneGrid = rand_search.rf$bestTune, 
                              metric = "Rsquared",  trControl = rand_ctrl)
final_search.xgb.test <- train(getConfirmedFormula(boruta.fs), data = test.ori, method = "xgbTree", 
                               tuneGrid = rand_search.xgb$bestTune, 
                               metric = "Rsquared",  trControl = rand_ctrl)
final_search.svm.test <- train(getConfirmedFormula(boruta.fs), data = test.ori, method = "svmRadial", 
                               tuneGrid = rand_search.svm$bestTune, preProc = c("center", "scale"),
                               metric = "Rsquared",  trControl = rand_ctrl)
final_search.earth.test.pca <- train(fdv~., data = test.pca, method = "earth", 
                                     tuneGrid = rand_search.earth.pca$bestTune, 
                                     metric = "Rsquared",  trControl = rand_ctrl)
final_search.rf.test.pca <- train(fdv~., data = test.pca, method = "extraTrees", 
                                  tuneGrid = rand_search.rf.pca$bestTune, 
                                  metric = "Rsquared",  trControl = rand_ctrl)
final_search.xgb.test.pca <- train(fdv~., data = test.pca, method = "xgbTree", 
                                   tuneGrid = rand_search.xgb.pca$bestTune, 
                                   metric = "Rsquared",  trControl = rand_ctrl)
final_search.svm.test.pca <- train(fdv~., data = test.pca, method = "svmRadial", 
                                   tuneGrid = rand_search.svm.pca$bestTune, preProc = c("center", "scale"),
                                   metric = "Rsquared",  trControl = rand_ctrl)


#function machine learning all (ga jadi dipake)
fun.ml <- function(boruta.fs){
  model.rf <- randomForest(getConfirmedFormula(boruta.fs), ntree = 1000, data = train)
  predicted.rf <- predict(model.rf, test)
  trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)
  model.svm <- caret::train(getConfirmedFormula(boruta.fs), data = train, method = "svmRadial", preProc = c("center", "scale"), metric = "RMSE", trControl=trctrl)
  predicted.svm <- predict(model.svm, test)
  model.mars <- earth(getConfirmedFormula(boruta.fs), data = train)
  predicted.mars <- predict(model.mars, test)
  model.xgb <- caret::train(getConfirmedFormula(boruta.fs), data = train, method = 'xgbTree')
  predicted.xgb <- predict(model.xgb, test)
  rsq.xgb <- cor(predicted.xgb, test$fdv)^2
  rsq.rf <- cor(predicted.rf, test$fdv)^2
  rsq.svm <- cor(predicted.svm, test$fdv)^2
  rsq.mars <- cor(predicted.mars, test$fdv)^2
  rsq.list <- c(rsq.rf, rsq.svm, rsq.mars, rsq.xgb)
  return(rsq.list)
}

result <- fun.ml(boruta.fs)

library(MLmetrics)
library(Metrics)
library(hydroGOF)

nrmse(predict.test.xgb.DART, test$fdv, norm = "maxmin")
bias(test$fdv, predict.test.svm)