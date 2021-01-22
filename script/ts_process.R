library(raster)
library(RStoolbox)
library(xts)
library(zoo)
library(rts)
library(signal)
library(zoo)
library(tsbox)
library(imputeTS)
library(timeSeries)
library(rHarmonics)
rasterOptions(progress = "Text", tmpdir = paste0(getwd(), "/tmp"), memfrac = .6)
setwd("G:/UGM_ARD/Landsat-8_ToA-BRDF/T49MHU_118061")
l8_masked <- list.files(getwd(), "geo_mask.tif", full.names = T, recursive = T)

setwd("G:/UGM_ARD/ToA/T49MHU")
sen_masked <- list.files(getwd(), "^calibrated_cor", recursive = T, full.names = T)
ref <- stack(list.files(getwd(), "310118"))

vec_l8 <- list()
vec_s2 <- list()

for(i in 1:length(l8_list)){
  name_l8 <- strsplit(l8_list[i], "m_")[[1]][2]
  name_l8 <- strsplit(name_l8, "_geo")[[1]][1]
  date_l8 <- as.Date(as.character(name_l8), format = "%d%m%y")
  vec_l8[i] <- as.character(date_l8)
}

for(i in 1:length(sen_masked)){
  name_sen <- strsplit(sen_masked[i], "T0")[[1]][1]
  name_sen <- strsplit(name_sen, "MHU_")[[1]][2]
  date_sen <- as.Date(as.character(name_sen), format = "%Y%m%d")
  vec_s2[i] <- as.character(date_sen)
}

recap <- list.files(getwd(), "^recap_pif", full.names = T, recursive = T)
pv_b1 <- list()
pv_b2 <- list()
pv_b3 <- list()
pv_b4 <- list()
pv_b5 <- list()
pv_b6 <- list()

for(i in 1:length(recap)){
  a <- read.csv(recap[i], header = T)
  pv_b1[i] <- a$rsq2a[1]
  pv_b2[i] <- a$rsq2a[2]
  pv_b3[i] <- a$rsq2a[3]
  pv_b4[i] <- a$rsq2a[4]
  pv_b5[i] <- a$rsq2a[5]
  pv_b6[i] <- a$rsq2a[6]
}

#filtering stack
sen2_b1 <- sen_masked[pv_b1 >= mean(unlist(pv_b1))]
sen2_b2 <- sen_masked[pv_b2 >= mean(unlist(pv_b2))]
sen2_b3 <- sen_masked[pv_b3 >= mean(unlist(pv_b3))]
sen2_b4 <- sen_masked[pv_b4 >= mean(unlist(pv_b4))]
sen2_b5 <- sen_masked[pv_b5 >= mean(unlist(pv_b5))]
sen2_b6 <- sen_masked[pv_b6 >= mean(unlist(pv_b6))]
date_b1 <- vec_s2[pv_b1 >= mean(unlist(pv_b1))]
date_b2 <- vec_s2[pv_b2 >= mean(unlist(pv_b2))]
date_b3 <- vec_s2[pv_b3 >= mean(unlist(pv_b3))]
date_b4 <- vec_s2[pv_b4 >= mean(unlist(pv_b4))]
date_b5 <- vec_s2[pv_b5 >= mean(unlist(pv_b5))]
date_b6 <- vec_s2[pv_b6 >= mean(unlist(pv_b6))]

#making raster time series stack
l8_b1 <- list()
l8_b2 <- list()
l8_b3 <- list()
l8_b4 <- list()
l8_b5 <- list()
l8_b6 <- list()

for(i in 1:length(l8_masked)){
  rst <- stack(l8_masked[i])
  rst <- calc(rst, fun.maskrm)
  l8_b1[i] <- rst[[2]]
  l8_b2[i] <- rst[[3]]
  l8_b3[i] <- rst[[4]]
  l8_b4[i] <- rst[[5]]
  l8_b5[i] <- rst[[6]]
  l8_b6[i] <- rst[[7]]
}

sen2_b1_rst <- list()
sen2_b2_rst <- list()
sen2_b3_rst <- list()
sen2_b4_rst <- list()
sen2_b5_rst <- list()
sen2_b6_rst <- list()

for(i in 1:length(sen2_b1)){
  rst <- stack(sen2_b1[i])
  rst <- rst[[1]]
  #rst <- resample(rst, ref, "bilinear")
  sen2_b1_rst[i] <- rst
}

for(i in 1:length(sen2_b2)){
  rst <- stack(sen2_b2[i])
  rst <- rst[[2]]
  #rst <- resample(rst, ref, "bilinear")
  sen2_b2_rst[i] <- rst
}

for(i in 1:length(sen2_b3)){
  rst <- stack(sen2_b3[i])
  rst <- rst[[3]]
  #rst <- resample(rst, ref, "bilinear")
  sen2_b3_rst[i] <- rst
}

for(i in 1:length(sen2_b4)){
  rst <- stack(sen2_b4[i])
  rst <- rst[[4]]
  #rst <- resample(rst, ref, "bilinear")
  sen2_b4_rst[i] <- rst
}

for(i in 1:length(sen2_b5)){
  rst <- stack(sen2_b5[i])
  rst <- rst[[5]]
  #rst <- resample(rst, ref, "bilinear")
  sen2_b5_rst[i] <- rst
}

for(i in 1:length(sen2_b6)){
  rst <- stack(sen2_b6[i])
  rst <- rst[[6]]
  #rst <- resample(rst, ref, "bilinear")
  sen2_b6_rst[i] <- rst
}

fun.median <- function(x) median(x, na.rm = T)
    
comb_b1 <- append(l8_b1, sen2_b1_rst)
comb_b2 <- append(l8_b2, sen2_b2_rst)
comb_b3 <- append(l8_b3, sen2_b3_rst)
comb_b4 <- append(l8_b4, sen2_b4_rst)
comb_b5 <- append(l8_b5, sen2_b5_rst)
comb_b6 <- append(l8_b6, sen2_b6_rst)
comb_b1_dates <- append(vec_l8, date_b1)
comb_b2_dates <- append(vec_l8, date_b2)
comb_b3_dates <- append(vec_l8, date_b3)
comb_b4_dates <- append(vec_l8, date_b4)
comb_b5_dates <- append(vec_l8, date_b5)
comb_b6_dates <- append(vec_l8, date_b6)
rts_b1 <- rts(stack(comb_b1), as.Date(unlist(comb_b1_dates)))
rts_b2 <- rts(stack(comb_b2), as.Date(unlist(comb_b2_dates)))
rts_b3 <- rts(stack(comb_b3), as.Date(unlist(comb_b3_dates)))
rts_b4 <- rts(stack(comb_b4), as.Date(unlist(comb_b4_dates)))
rts_b5 <- rts(stack(comb_b5), as.Date(unlist(comb_b5_dates)))
rts_b6 <- rts(stack(comb_b6), as.Date(unlist(comb_b6_dates)))
rts_b1_monmed <- apply.monthly(rts_b1, fun.median)
rts_b2_monmed <- apply.monthly(rts_b2, fun.median)
rts_b3_monmed <- apply.monthly(rts_b3, fun.median)
rts_b4_monmed <- apply.monthly(rts_b4, fun.median)
rts_b5_monmed <- apply.monthly(rts_b5, fun.median)
rts_b6_monmed <- apply.monthly(rts_b6, fun.median)
rts_b1_qmed <- apply.quarterly(rts_b1, fun.median)
rts_b2_qmed <- apply.quarterly(rts_b2, fun.median)
rts_b3_qmed <- apply.quarterly(rts_b3, fun.median)
rts_b4_qmed <- apply.quarterly(rts_b4, fun.median)
rts_b5_qmed <- apply.quarterly(rts_b5, fun.median)
rts_b6_qmed <- apply.quarterly(rts_b6, fun.median)
removeTmpFiles(h=0.01)

writeRaster(rts_b1@raster, paste0("T49MHU_monmed_", comb_b1_dates[1], "_", comb_b1_dates[length(comb_b1_dates)], "_b1_raw.tif"), format = "GTiff", datatype = "FLT4S")
writeRaster(rts_b2@raster, paste0("T49MHU_monmed_", comb_b2_dates[1], "_", comb_b2_dates[length(comb_b2_dates)], "_b2_raw.tif"), format = "GTiff", datatype = "FLT4S")
writeRaster(rts_b3@raster, paste0("T49MHU_monmed_", comb_b3_dates[1], "_", comb_b3_dates[length(comb_b3_dates)], "_b3_raw.tif"), format = "GTiff", datatype = "FLT4S")
writeRaster(rts_b4@raster, paste0("T49MHU_monmed_", comb_b4_dates[1], "_", comb_b4_dates[length(comb_b4_dates)], "_b4_raw.tif"), format = "GTiff", datatype = "FLT4S")
writeRaster(rts_b5@raster, paste0("T49MHU_monmed_", comb_b5_dates[1], "_", comb_b5_dates[length(comb_b5_dates)], "_b5_raw.tif"), format = "GTiff", datatype = "FLT4S")
writeRaster(rts_b6@raster, paste0("T49MHU_monmed_", comb_b6_dates[1], "_", comb_b6_dates[length(comb_b6_dates)], "_b6_raw.tif"), format = "GTiff", datatype = "FLT4S")

writeRaster(rts_b1_monmed@raster, paste0("T49MHU_monmed_", comb_b1_dates[1], "_", comb_b1_dates[length(comb_b1_dates)], "_b1.tif"), format = "GTiff", datatype = "FLT4S")
writeRaster(rts_b2_monmed@raster, paste0("T49MHU_monmed_", comb_b2_dates[1], "_", comb_b2_dates[length(comb_b2_dates)], "_b2.tif"), format = "GTiff", datatype = "FLT4S")
writeRaster(rts_b3_monmed@raster, paste0("T49MHU_monmed_", comb_b3_dates[1], "_", comb_b3_dates[length(comb_b3_dates)], "_b3.tif"), format = "GTiff", datatype = "FLT4S")
writeRaster(rts_b4_monmed@raster, paste0("T49MHU_monmed_", comb_b4_dates[1], "_", comb_b4_dates[length(comb_b4_dates)], "_b4.tif"), format = "GTiff", datatype = "FLT4S")
writeRaster(rts_b5_monmed@raster, paste0("T49MHU_monmed_", comb_b5_dates[1], "_", comb_b5_dates[length(comb_b5_dates)], "_b5.tif"), format = "GTiff", datatype = "FLT4S")
writeRaster(rts_b6_monmed@raster, paste0("T49MHU_monmed_", comb_b6_dates[1], "_", comb_b6_dates[length(comb_b6_dates)], "_b6.tif"), format = "GTiff", datatype = "FLT4S")

writeRaster(rts_b1_qmed@raster, paste0("T49MHU_qmed_", comb_b1_dates[1], "_", comb_b1_dates[length(comb_b1_dates)], "_b1.tif"), format = "GTiff", datatype = "FLT4S")
writeRaster(rts_b2_qmed@raster, paste0("T49MHU_qmed_", comb_b2_dates[1], "_", comb_b2_dates[length(comb_b2_dates)], "_b2.tif"), format = "GTiff", datatype = "FLT4S")
writeRaster(rts_b3_qmed@raster, paste0("T49MHU_qmed_", comb_b3_dates[1], "_", comb_b3_dates[length(comb_b3_dates)], "_b3.tif"), format = "GTiff", datatype = "FLT4S")
writeRaster(rts_b4_qmed@raster, paste0("T49MHU_qmed_", comb_b4_dates[1], "_", comb_b4_dates[length(comb_b4_dates)], "_b4.tif"), format = "GTiff", datatype = "FLT4S")
writeRaster(rts_b5_qmed@raster, paste0("T49MHU_qmed_", comb_b5_dates[1], "_", comb_b5_dates[length(comb_b5_dates)], "_b5.tif"), format = "GTiff", datatype = "FLT4S")
writeRaster(rts_b6_qmed@raster, paste0("T49MHU_qmed_", comb_b6_dates[1], "_", comb_b6_dates[length(comb_b6_dates)], "_b6.tif"), format = "GTiff", datatype = "FLT4S")

fun.svg <- function(x){
  if (all(is.na(x))){
    return(rep(NA, length(x)))
  } else if (sum(!is.na(x)) == 1){
    return(rep(NA, length(x)))
  } else {
    x1 <- zoo(c(x), c(as.Date(unlist(date.list))))
    x1 <- as.vector(as.vector(unlist(x1)))
    x2 <- na_interpolation(x1, option = "linear")
    x1.sg <- sgolayfilt(x2, p = 2, n = 7)
    return(x1.sg)
  }
}

rts_b1 <- rts(stack(comb_b1), as.Date(unlist(comb_b1_dates)))
rts_b2 <- rts(stack(comb_b2), as.Date(unlist(comb_b2_dates)))
rts_b3 <- rts(stack(comb_b3), as.Date(unlist(comb_b3_dates)))
rts_b4 <- rts(stack(comb_b4), as.Date(unlist(comb_b4_dates)))
rts_b5 <- rts(stack(comb_b5), as.Date(unlist(comb_b5_dates)))
rts_b6 <- rts(stack(comb_b6), as.Date(unlist(comb_b6_dates)))

#rts_b5_sm <- calc(rts_b5@raster, fun.svg)
rts_b5_sm <- calc(stack(comb_b5), fun.svg)
rts_b5_sm2 <- rts(rts_b5_sm, timeSeries::sort(as.Date(unlist(comb_b5_dates))))
rts_b5_monmed2 <- apply.monthly(rts_b5_sm2, fun.median)
rts_b5_qmed2 <- apply.quarterly(rts_b5_sm2, fun.median)

date.list <- as.vector(unlist(comb_b4_dates))
rts_b4_sm <- calc(stack(comb_b4), fun.svg)
date.list <- as.vector(unlist(comb_b5_dates))
rts_b5_sm <- calc(stack(comb_b5), fun.svg)
date.list <- as.vector(unlist(comb_b6_dates))
rts_b6_sm <- calc(stack(comb_b6), fun.svg)

rts_b4_sm2 <- rts(rts_b4_sm, timeSeries::sort(as.Date(unlist(comb_b4_dates))))
rts_b4_monmed2 <- apply.monthly(rts_b4_sm2, fun.median)
rts_b4_qmed2 <- apply.quarterly(rts_b4_sm2, fun.median)

rts_b5_sm2 <- rts(rts_b5_sm, timeSeries::sort(as.Date(unlist(comb_b5_dates))))
rts_b5_monmed2 <- apply.monthly(rts_b5_sm2, fun.median)
rts_b5_qmed2 <- apply.quarterly(rts_b5_sm2, fun.median)

rts_b6_sm2 <- rts(rts_b6_sm, timeSeries::sort(as.Date(unlist(comb_b6_dates))))
rts_b6_monmed2 <- apply.monthly(rts_b6_sm2, fun.median)
rts_b6_qmed2 <- apply.quarterly(rts_b6_sm2, fun.median)


rts_b1_monmed <- apply.monthly(rts_b1, fun.median)
rts_b2_monmed <- apply.monthly(rts_b2, fun.median)
rts_b3_monmed <- apply.monthly(rts_b3, fun.median)
rts_b4_monmed <- apply.monthly(rts_b4, fun.median)
rts_b5_monmed <- apply.monthly(rts_b5, fun.median)
rts_b6_monmed <- apply.monthly(rts_b6, fun.median)
rts_b1_qmed <- apply.quarterly(rts_b1, fun.median)
rts_b2_qmed <- apply.quarterly(rts_b2, fun.median)
rts_b3_qmed <- apply.quarterly(rts_b3, fun.median)
rts_b4_qmed <- apply.quarterly(rts_b4, fun.median)
rts_b5_qmed <- apply.quarterly(rts_b5, fun.median)
rts_b6_qmed <- apply.quarterly(rts_b6, fun.median)
removeTmpFiles(h=0.01)

#new2
fun.svg2 <- function(x){
  if (all(is.na(x))){
    return(rep(NA, length(x)))
  } else if (sum(!is.na(x)) == 1){
    return(rep(NA, length(x)))
  } else {
    x1 <- as.vector(as.vector(unlist(x)))
    x2 <- na_interpolation(x1, option = "spline")
    x1.sg <- sgolayfilt(x2, p = 2, n = 7)
    return(x1.sg)
  }
}

#new
rts_b4_monmed2 <- calc(rts_b4_monmed@raster, fun.svg)
rts_b5_monmed2 <- calc(rts_b5_monmed@raster, fun.svg)
rts_b6_monmed2 <- calc(rts_b6_monmed@raster, fun.svg)

i2 <- rep(1:(nlayers(rts_b6_monmed2)/2), each =2)
i3 <- rep(1:(nlayers(rts_b6_monmed2)/3), each =3)
i4 <- rep(1:(nlayers(rts_b6_monmed2)/4), each =4)
i6 <- rep(1:(nlayers(rts_b6_monmed2)/6), each =6)
i12 <- rep(1:(nlayers(rts_b6_monmed2)/12), each =12)

fun.median <- function(x){
  b <- rnorm(length(x))
  b <- median(b, na.rm = T)
  if (all(is.na(x))){
    return(rep(NA, length(b)))
  } else {
    c <- median(x, na.rm = T)
    return(c)
  }
}

rts_b4_qmed2 <- stackApply(rts_b4_monmed2, i3, median)
rts_b5_qmed2 <- stackApply(rts_b5_monmed2, i3, median)
rts_b6_qmed2 <- stackApply(rts_b6_monmed2, i3, median)

rts_b4_half <- stackApply(rts_b4_monmed2, i6, median)
rts_b5_half <- stackApply(rts_b5_monmed2, i6, median)
rts_b6_half <- stackApply(rts_b6_monmed2, i6, median)

fun.qaflag <- function(x){
  if (all(is.na(x))){
    return(NA)
  } else if (sum(is.na(x)) > (length(x)*0.9)){
    return(9)
  } else if ((sum(is.na(x)) > (length(x)*0.8)) & (sum(is.na(x)) <= (length(x)*0.9))) {
    return(8)
  } else if ((sum(is.na(x)) > (length(x)*0.7)) & (sum(is.na(x)) <= (length(x)*0.8))) {
    return(7)
  } else if ((sum(is.na(x)) > (length(x)*0.6)) & (sum(is.na(x)) <= (length(x)*0.7))) {
    return(6)
  } else if ((sum(is.na(x)) > (length(x)*0.5)) & (sum(is.na(x)) <= (length(x)*0.6))) {
    return(5)
  } else if ((sum(is.na(x)) > (length(x)*0.4)) & (sum(is.na(x)) <= (length(x)*0.5))) {
    return(4)
  } else if ((sum(is.na(x)) > (length(x)*0.3)) & (sum(is.na(x)) <= (length(x)*0.4))) {
    return(3)
  } else if ((sum(is.na(x)) > (length(x)*0.2)) & (sum(is.na(x)) <= (length(x)*0.3))) {
    return(2)
  } else {
    return(1)
  }
}


fun.harmonics <- function(x){
  if (all(is.na(x))){
    return(rep(NA, length(x)))
  } else if (sum(is.na(x)) > (length(x)*0.3)){
    return(rep(NA, length(x)))
  } else {
    x1 <- zoo(c(x), c(as.Date(unlist(date.list))))
    fitted_3rd_deg <- harmonics_fun(user_vals = as.vector(x1),
                                    user_dates = timeSeries::sort(as.Date(unlist(date.list))),
                                    harmonic_deg = 3)
    return(c(fitted_3rd_deg))
  }
}


library(doParallel)
library(snow)
library(parallel)

cl <- makeCluster(8)
beginCluster(cl)

date.list <- as.vector(unlist(comb_b4_dates))
#rts_b4_sm5 <- calc(stack(comb_b4), fun.harmonics)
rts_b4_sm5 <- clusterR(stack(comb_b4), calc, args = list(fun.harmonics), export = c("date.list"))
date.list <- as.vector(unlist(comb_b5_dates))
#rts_b5_sm5 <- calc(stack(comb_b5), fun.harmonics)
rts_b5_sm5 <- clusterR(stack(comb_b5), calc, args = list(fun.harmonics), export = c("date.list"))
date.list <- as.vector(unlist(comb_b6_dates))
#rts_b6_sm5 <- calc(stack(comb_b6), fun.harmonics)
rts_b6_sm5 <- clusterR(stack(comb_b6), calc, args = list(fun.harmonics), export = c("date.list"))

endCluster()
stopCluster(cl)

rts_b4_sm5a <- rts(rts_b4_sm5, timeSeries::sort(as.Date(unlist(comb_b4_dates))))
rts_b5_sm5a <- rts(rts_b5_sm5, timeSeries::sort(as.Date(unlist(comb_b5_dates))))
rts_b6_sm5a <- rts(rts_b6_sm5, timeSeries::sort(as.Date(unlist(comb_b6_dates))))


rts_b4_monmed5 <- apply.monthly(rts_b4_sm5a, median)
rts_b5_monmed5 <- apply.monthly(rts_b5_sm5a, median)
rts_b6_monmed5 <- apply.monthly(rts_b6_sm5a, median)

rts_b4_qmed5 <- apply.quarterly(rts_b4_sm5a, median)
rts_b5_qmed5 <- apply.quarterly(rts_b5_sm5a, median)
rts_b6_qmed5 <- apply.quarterly(rts_b6_sm5a, median)
