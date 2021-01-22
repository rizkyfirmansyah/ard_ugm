library(foreign)
dbf.files <- list.files(getwd(), ".dbf")
df <- lapply(dbf.files, read.dbf)
df.final <- do.call("rbind", df)
train <- data.frame(df.final)
rasterOptions(progress = "Text", tmpdir = paste0(getwd(), "/tmp"), memfrac = .6)

library(caret)
rand_search.rf <- train(factor(id)~., data = df.final, method = "rf", tuneLength = 50, metric = "Accuracy", trControl = rand_ctrl)
final_search.rf <- train(factor(id)~., data = df.final, method = "rf", tuneLength = rand_search.rf$bestTune, metric = "Accuracy")

setwd("G:/UGM_ARD/ToA/T49MHU")
list.dir.sen2 <- list.dirs(getwd(), full.names = T, recursive = F)

for(i in 1:length(list.dir.sen2)){
  setwd(list.dir.sen2[i])
  rst <- list.files(getwd(), "calibrated_cor_Geo")
  print(paste0(basename(rst)))
  rstr <- stack(rst)
  names(rstr) <- c("calibrat_1", "calibrat_2", "calibrat_3", "calibrat_4", "calibrat_5", "calibrat_6" )
  cm <- raster::predict(rstr, final_search.rf, ext = extent(rstr), progress = "text")
  cm[cm == 3] <- NA
  writeRaster(cm, paste0("CM_RF_", basename(rst)), format = "GTiff", datatype = "INT1U")
}

for(i in 1:length(list.dir.sen2)){
  setwd(list.dir.sen2[i])
  rst2 <- list.files(getwd(), "^calibrated_cor_Geo")
  cm <- list.files(getwd(), "CM_RF")
  print(paste0(basename(rst2)))
  print(paste0(basename(cm)))
  cm_rst <- raster(cm)
  cm_rst[is.na(cm_rst)] <- 0
  cm_rst[cm_rst > 0] <- NA
  cm_rst[cm_rst == 0] <- 1
  rst <- stack(rst2)
  rst_mask <- rst * cm_rst
  writeRaster(rst_mask, paste0("masked_CM_", basename(rst2)), format = "GTiff", datatype = "FLT4S")
}

