library(raster)
library(RStoolbox)
library(parallel)
library(doParallel)
library(foreach)
library(snow)
library(XML)
rasterOptions(progress = "Text", tmpdir = paste0(getwd(), "/tmp"), memfrac = .6)

rsq1a <- list()
rsq2a <- list()
pv1a <- list()

rsq1b <- list()
rsq2b <- list()
pv1b <- list()

#list.sentinel2 <- list.files(getwd(), "_25.tif", recursive = T)
dir.lst <- list.dirs(getwd(), full.names = T)

#band pass adjustment function
fun.bpa.s2a <- function(x){
 if (all(is.na(x))){
   return(rep(NA, length(x)))
 } else {
   x1 <- (0.9778 * (x[1])) - (0.004 * 10000)
   x2 <- (1.0053 * (x[2])) - (0.0009 * 10000)
   x3 <- (0.9765 * (x[3])) + (0.0009 * 10000)
   x4 <- (0.9983 * (x[4])) - (0.0001 * 10000)
   x5 <- (0.9987 * (x[5])) - (0.0011 * 10000)
   x6 <- (1.003 * (x[6])) - (0.0012 * 10000)
  return(c(x1, x2, x3, x4, x5, x6))
 }
}

fun.bpa.s2b <- function(x){
  if (all(is.na(x))){
    return(rep(NA, length(x)))
  } else {
    x1 <- (0.9778 * (x[1])) - (0.004 * 10000)
    x2 <- (1.0075 * (x[2])) - (0.0008 * 10000)
    x3 <- (0.9761 * (x[3])) + (0.001 * 10000)
    x4 <- (0.9966 * (x[4]))
    x5 <- (1 * (x[5])) - (0.0003 * 10000)
    x6 <- (0.9867 * (x[6])) + (0.0004 * 10000)
    return(c(x1, x2, x3, x4, x5, x6))
  }
}

#for(i in 1:length(list.sentinel2)){
#parallel computation
cl <- makeCluster(2)
registerDoParallel(cl)

#setting reference file, ganti filenya untuk setiap tile
reference.file <- "L8LTP118061m_200318_geo_mask.tif"
ref <- stack(reference.file)
ref <- stack(ref[[2:7]])

foreach(i = 1:length(dir.lst), packages = c("raster", "RStoolbox")) %dopar% {
  setwd(dir.lst[i])
  data <- xmlParse("MTD_MSIL1C.xml")
  xml_data <- xmlToList(data)
  sen.type <- xml_data$General_Info$Product_Info$Datatake$SPACECRAFT_NAME
  list.sentinel2 <- list.files(getwd(), "_25.tif", recursive = T)
  sen.uncal <- stack(list.sentinel2)
  name <- strsplit(basename(list.sentinel2[i]), ".tif")[[1]][1]
  print(paste0("start processing ", name))
  if (sen.type == "Sentinel-2A"){
   sen.bpa <- calc(sen.uncal, fun.bpa.s2a)
  } else {
   sen.bpa <- calc(sen.uncal, fun.bpa.s2b)
  } 
  sen.pif1 <- pifMatch(sen.bpa, ref, method = "cor", quantile = 0.95, returnModels = T)
  sen.pif2 <- pifMatch(sen.bpa, ref, method = "ed", quantile = 0.95, returnModels = T)
  for(j in 1:nlayers(sen.uncal)){
    sink(paste0("pif1_cor_", name, "_b", j, ".txt"))
    print(summary(sen.pif1$model[[j]]))
    sink()
    b <- summary(sen.pif1$model[[j]])
    rsq1a[j] <- b$r.squared
    rsq2a[j] <- b$adj.r.squared
    pv1a[j] <- b$coefficients[2,4]  
    print(paste0("adj rsq of ", j, " is ", b$adj.r.squared))
    sink(paste0("pif2_ed_", name, "_b", j, ".txt"))
    print(summary(sen.pif2$model[[j]]))
    sink()
    b <- summary(sen.pif2$model[[j]])
    rsq1b[j] <- b$r.squared
    rsq2b[j] <- b$adj.r.squared
    pv1b[j] <- b$coefficients[2,4]  
    print(paste0("adj rsq of ", j, " is ", b$adj.r.squared))
  }
  writeRaster(sen.pif1$img, paste0("calibrated_cor_", name, ".tif"), format= "GTiff", datatype = "FLT4S")
  writeRaster(sen.pif2$img, paste0("calibrated_ed_", name, ".tif"), format= "GTiff", datatype = "FLT4S")
  recap <- cbind(rsq1a, rsq2a, pv1a, rsq1b, rsq2b, pv1b)
  names(recap) <- c("rsq1", "adj.rsq1", "pvalue1", "rsq2", "adj.rsq2", "pvalue2")
  write.csv(recap, paste0("recap_pif_", name, ".csv"))
  writeRaster(sen.pif1$pifMap, paste0("pif1_area_cor_", name, ".tif"), format= "GTiff", datatype = "FLT4S")
  writeRaster(sen.pif2$pifMap, paste0("pif2_area_ed_", name, ".tif"), format= "GTiff", datatype = "FLT4S")
}
stopCluster(cl)
