rasterOptions(progress="Text", tmpdir=paste0(getwd(),"/tmp"))
dir <- list.dirs(getwd(), full.names = T, recursive = F)
for(i in 1:length(dir)){
  setwd(dir[i])
  ers_files<-list.files(getwd(),".ers")
  print("jumlah file = ", length(ers_files))
for(j in 1:length(ers_files)){
  rst1 <- stack(ers_files[j])
  print(paste0("start",j,"from", length(ers_files)))
  names <- strsplit(ers_files[j], ".ers")[[1]][1]
  writeRaster(rst1, paste0(names, ".tif"), format = "GTiff", datatype = "FLT4S", overwrite = T)
}
}
coba <- list.files(getwd(), full.names = T, recursive = F)
tiff_files <- list.files(getwd(),".tif")
file.remove(tiff_files)
for (i in 1:length(dir)){
  setwd(dir[i])
  tiff_files <- list.files(getwd(), ".tif")
  file.remove(tiff_files)
}
for(i in 1:length(dir)){
  setwd(dir[i])
  dataFiles <- list.files(getwd(), "CM.tif$")}
