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

list.sentinel2 <- list.files(getwd(), "_25.tif", recursive = T)

#for(i in 1:length(list.sentinel2)){
#parallel computation
cl <- makeCluster(2)
registerDoParallel(cl)

#setting reference file, ganti filenya untuk setiap tile
reference.file <- "L8LTP118061m_200318_geo_mask.tif"
ref <- stack(reference.file)
ref <- stack(ref[[2:7]])

foreach(i = 1:length(list.sentinel2), packages = c("raster", "RStoolbox")) %dopar% {
  sen.uncal <- stack(list.sentinel2[i])
  name <- strsplit(basename(list.sentinel2[i]), ".tif")[[1]][1]
  print(paste0("start processing ", name))
  sen.pif1 <- pifMatch(sen.uncal, ref, method = "cor", quantile = 0.95, returnModels = T)
  sen.pif2 <- pifMatch(sen.uncal, ref, method = "ed", quantile = 0.95, returnModels = T)
  for(j in 1:nlayers(sen.uncal){
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
