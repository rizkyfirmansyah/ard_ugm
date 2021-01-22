
fun.maskrm1 <- function(x){
  if (all(is.na(x))){
    return(rep(NA, length(x)))
  } else if (sum(x) == 0){
    return(rep(NA, length(x)))
  } else {
    return(x)
  }
}

fun.maskrm <- function(x){
  if (all(is.na(x))){
    return(rep(NA, length(x)))
  } else if (sum(x) == 0){
    return(rep(NA, length(x)))
  } else if (any(x == 0)){
    return(rep(NA, length(x)))
  } else {
    return(x)
  }
}

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
sen2_b1 <- sen_masked[pv_b5 >= mean(unlist(pv_b5))]
sen2_b2 <- sen_masked[pv_b5 >= mean(unlist(pv_b5))]
sen2_b3 <- sen_masked[pv_b5 >= mean(unlist(pv_b5))]
sen2_b4 <- sen_masked[pv_b5 >= mean(unlist(pv_b5))]
sen2_b5 <- sen_masked[pv_b5 >= mean(unlist(pv_b5))]
sen2_b6 <- sen_masked[pv_b5 >= mean(unlist(pv_b5))]
date_b1 <- vec_s2[pv_b5 >= mean(unlist(pv_b5))]
date_b2 <- vec_s2[pv_b5 >= mean(unlist(pv_b5))]
date_b3 <- vec_s2[pv_b5 >= mean(unlist(pv_b5))]
date_b4 <- vec_s2[pv_b5 >= mean(unlist(pv_b5))]
date_b5 <- vec_s2[pv_b5 >= mean(unlist(pv_b5))]
date_b6 <- vec_s2[pv_b5 >= mean(unlist(pv_b5))]

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
  rst <- calc(rst, fun.maskrm)
  rst <- rst[[1]]
  #rst <- resample(rst, ref, "bilinear")
  sen2_b1_rst[i] <- rst
}

for(i in 1:length(sen2_b2)){
  rst <- stack(sen2_b2[i])
  rst <- calc(rst, fun.maskrm)
  rst <- rst[[2]]
  #rst <- resample(rst, ref, "bilinear")
  sen2_b2_rst[i] <- rst
}

for(i in 1:length(sen2_b3)){
  rst <- stack(sen2_b3[i])
  rst <- calc(rst, fun.maskrm)
  rst <- rst[[3]]
  #rst <- resample(rst, ref, "bilinear")
  sen2_b3_rst[i] <- rst
}

for(i in 1:length(sen2_b4)){
  rst <- stack(sen2_b4[i])
  rst <- calc(rst, fun.maskrm)
  rst <- rst[[4]]
  #rst <- resample(rst, ref, "bilinear")
  sen2_b4_rst[i] <- rst
}

for(i in 1:length(sen2_b5)){
  rst <- stack(sen2_b5[i])
  rst <- calc(rst, fun.maskrm)
  rst <- rst[[5]]
  #rst <- resample(rst, ref, "bilinear")
  sen2_b5_rst[i] <- rst
}

for(i in 1:length(sen2_b6)){
  rst <- stack(sen2_b6[i])
  rst <- calc(rst, fun.maskrm)
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

date.list <- as.vector(unlist(comb_b1_dates))
rts_b1_sm5 <- clusterR(stack(comb_b1), calc, args = list(fun.harmonics), export = c("date.list"))
date.list <- as.vector(unlist(comb_b2_dates))
rts_b2_sm5 <- clusterR(stack(comb_b2), calc, args = list(fun.harmonics), export = c("date.list"))
date.list <- as.vector(unlist(comb_b3_dates))
rts_b3_sm5 <- clusterR(stack(comb_b3), calc, args = list(fun.harmonics), export = c("date.list"))
date.list <- as.vector(unlist(comb_b4_dates))
rts_b4_sm5 <- clusterR(stack(comb_b4), calc, args = list(fun.harmonics), export = c("date.list"))
date.list <- as.vector(unlist(comb_b5_dates))
rts_b5_sm5 <- clusterR(stack(comb_b5), calc, args = list(fun.harmonics), export = c("date.list"))
date.list <- as.vector(unlist(comb_b6_dates))
rts_b6_sm5 <- clusterR(stack(comb_b6), calc, args = list(fun.harmonics), export = c("date.list"))