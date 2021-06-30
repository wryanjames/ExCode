#'""" Hypervolumes for CELA trophic niche
#'    @author: Ryan James
#'    Date: 1/19/21"""

library(tidyverse)
library(hypervolume)

# format for hv analysis ----
# updated on 3/23/2020****
# function to convert MixSIAR output to a table to generate random data
# to use for hypervolumes
# file = name of .txt of summary statistics of MixSIAR
# site = 'site' or 'location' of the data set
# ind = if true will output with columns as source values
# nest = if nested mixing model, T will return the sources of nested 
mixTable = function(file,site,ind = F,nest = F){
  require(tidyverse)
  cn = c('ID', 'Mean', 'SD', '2.5%', '5%', '25%', '50%', '75%', '95%', '97.5%')
  x = read_table(file, skip = 7,  col_names = T)
  names(x) = cn
  x$source = NA
  x$name = NA
  x$code = NA
  
  if (nest == F){
    for (i in 1:nrow(x)){
      temp = strsplit(x$ID, split = '.', fixed = T)
      x$source[i] = temp[[i]][3]
      x$name[i] = temp[[i]][2]
      
      x$site = site
      x$ymax = x$`75%` + 1.5*(x$`75%` - x$`25%`)
      x$ymin = x$`25%` - 1.5*(x$`75%` - x$`25%`)
      
      df = data.frame(x$name, x$site, x$source, x$Mean, x$SD, x$`2.5%`, x$`97.5%`,
                      x$`50%`, x$`25%`, x$`75%`, x$ymax, x$ymin)
      colnames(df) = c('name', 'site', 'source', 'mean', 'sd', 'lowend', 'highend',
                       'mid', 'low', 'up', 'ymax', 'ymin')
    }
  }else{
    for (i in 1:nrow(x)){
      temp = strsplit(x$ID, split = '.', fixed = T)
      x$source[i] = temp[[i]][4]
      x$code[i] = temp[[i]][3]
      x$name[i] = temp[[i]][2]
      
      x$site = site
      x$ymax = x$`75%` + 1.5*(x$`75%` - x$`25%`)
      x$ymin = x$`25%` - 1.5*(x$`75%` - x$`25%`)
      
      df = tibble(x$name, x$site, x$source, x$code, x$Mean, x$SD, x$`2.5%`, x$`97.5%`,
                  x$`50%`, x$`25%`, x$`75%`, x$ymax, x$ymin)
      colnames(df) = c('name', 'site', 'source', 'code', 'mean', 'sd', 'lowend', 'highend',
                       'mid', 'low', 'up', 'ymax', 'ymin')
    }
  }
  
  for (i in 1:nrow(df)){
    if (df$ymax[i] > df$highend[i]){
      df$ymax[i] = df$highend[i]
    }
    if (df$ymin[i] < df$lowend[i]){
      df$ymin[i] = df$lowend[i]
    }
  }
  df = na.omit(df)
  df = subset(df, name != 'global') 
  
  if (ind == T){
    if (nest == T){
      #df = data.frame(x$name, x$site, x$code, x$source, x$Mean)
      #colnames(df) = c('name', 'site', 'code', 'source', 'mean')
      df = df %>% select(name, site, code, source, mean)
      df = pivot_wider(df, names_from = 'source', values_from = 'mean')
    }else{
      df = df %>% select(name, site, source, mean)
      df = pivot_wider(df, names_from = 'source', values_from = 'mean')
    }
  }
  
  return(df)
}

# mixing model results 
ac = mixTable('data/ac4prey_ss.txt', 'ACS', T, T)
mc = mixTable('data/mc4prey_ss.txt', 'MCS', T, T)

df = bind_rows(ac, mc)%>% 
  mutate(ind = parse_number(code)) %>% 
  select(Species = name, System = site, ind, FB, MBD, MBZ, MP)


# generate z scores
df = df %>%  
  mutate_if(is.numeric, scale) 


# hypervolumes ----

# function to bootstrap hypervolumes to generate confidence intervals
# data = dataframe or tibble of z-scored values to make hvs
# n = number of bootstrapped hvs to make
# frac = fraction of data points to use to make bootstrapped hvs
#    ** can also be input as a whole number for how many  
#       data points to use in resampling of bootstraps
# Returns vector of bootstrapped hvs
HVboot = function(data, n, frac, vol = NULL){
  require(tidyverse)
  df = tibble(Volume = rep(NA, len = num))
  num = ceiling(frac*nrow(data))
  if(frac>1){
    num = frac
  }
  for (i in 1:n){
    rand = sample_n(data, num, replace = TRUE)
    require(hypervolume)
    randh = hypervolume_gaussian(rand,
                                 samples.per.point = 5000,
                                 kde.bandwidth = estimate_bandwidth(rand), 
                                 sd.count = 3, 
                                 quantile.requested = 0.95, 
                                 quantile.requested.type = "probability", 
                                 chunk.size = 1000, 
                                 verbose = F)
    df$Volume[i] = as.numeric(get_volume(randh))
    if (i/num == .1){
      cat('10% done\n')
    } else if (i/num == 0.25){
      cat('25% done\n')
    }else if (i/num == 0.50){
      cat('50% done\n')
    }else if (i/num == 0.75){
      cat('75% done\n')
    }else if (i/num == 0.9){
      cat('90% done\n')
    }
  }
  if (is.null(vol) == F){
    cat('Volume (95% CI) = ', vol, ' (', 
        quantile(df$Volume, 0.025, na.rm = T), '-', 
        quantile(df$Volume, 0.975, na.rm = T), ')\n')
  }
  return(df)
}


# generate hvs ----

# Snook in ACS
acSN = df %>% filter(Species == 'Snook', System == 'ACS') %>% 
  select(FB, MBD, MBZ, MP)

# generate hypervolume
acSNhv = hypervolume_gaussian(acSN, name = 'Snook ACS',
                              #samples.per.point = ceiling((10^(3 + sqrt(ncol(acSN))))/nrow(acSN)),
                              samples.per.point = 5000,
                              kde.bandwidth = estimate_bandwidth(acSN), 
                              sd.count = 3, 
                              quantile.requested = 0.95, 
                              quantile.requested.type = "probability", 
                              chunk.size = 1000, 
                              verbose = F)
# get volume 
# value = 121.38
get_volume(acSNhv)

# calculate CI
# Volume (95% CI) =  121.4  ( 64.50002 - 206.9598 )       
acSNhvCI = HVboot(acSN, 100, 0.75, vol = 121.4)
acSNhvCI$Species = 'Snook'
acSNhvCI$System = 'ACS'


# Snook in MCS
mcSN = df %>% filter(Species == 'Snook', System == 'MCS') %>% 
  select(FB, MBD, MBZ, MP)

# generate hypervolume
mcSNhv = hypervolume_gaussian(mcSN, name = 'Snook MCS',
                              #samples.per.point = ceiling((10^(3 + sqrt(ncol(mcSN))))/nrow(mcSN)),
                              samples.per.point = 5000,
                              kde.bandwidth = estimate_bandwidth(mcSN), 
                              sd.count = 3, 
                              quantile.requested = 0.95, 
                              quantile.requested.type = "probability", 
                              chunk.size = 1000, 
                              verbose = F)
# get volume 
# value = 3.04
get_volume(mcSNhv)

# calculate CI
# Volume (95% CI) =  3  ( 0.809157 - 6.022489 )       
mcSNhvCI = HVboot(mcSN, 100, 0.75, vol = 3.0)
mcSNhvCI$Species = 'Snook'
mcSNhvCI$System = 'MCS'

# Tarpon in ACS
acTR = df %>% filter(Species == 'Tarpon', System == 'ACS') %>% 
  select(FB, MBD, MBZ, MP)

# generate hypervolume
acTRhv = hypervolume_gaussian(acTR, name = 'Tarpon AC',
                              samples.per.point = 5000,
                              kde.bandwidth = estimate_bandwidth(acTR), 
                              sd.count = 3, 
                              quantile.requested = 0.95, 
                              quantile.requested.type = "probability", 
                              chunk.size = 1000, 
                              verbose = F)
# get volume 
# value = 32.1
get_volume(acTRhv)

# calculate CI
# Volume (95% CI) =  32.1  ( 4.668775 - 72.18363 )     
acTRhvCI = HVboot(acTR, 100, 0.75, vol = 32.1)
acTRhvCI$Species = 'Tarpon'
acTRhvCI$System = 'ACS'

# Tarpon in MCS
mcTR = df %>% filter(Species == 'Tarpon', System == 'MCS') %>% 
  select(FB, MBD, MBZ, MP)

# generate hypervolume
mcTRhv = hypervolume_gaussian(mcTR, name = 'Tarpon MCS',
                              samples.per.point = 5000,
                              kde.bandwidth = estimate_bandwidth(mcTR), 
                              sd.count = 3, 
                              quantile.requested = 0.95, 
                              quantile.requested.type = "probability", 
                              chunk.size = 1000, 
                              verbose = F)
# get volume 
# value = 1.23
get_volume(mcTRhv)

# calculate CI
# Volume (95% CI) =  1.2  ( 0.0753357 - 3.86419 )     
mcTRhvCI = HVboot(mcTR, 100, 0.75, vol = 1.2)
mcTRhvCI$Species = 'Tarpon'
mcTRhvCI$System = 'MCS'

# bind all data together for plotting 
hvCI = bind_rows(acSNhvCI, mcSNhvCI, acTRhvCI, mcTRhvCI)
write_csv(hvCI, 'data/hv_CI.csv')

# overlap 
# data1 = z scored data for hv 1
# data2 = z scored data for hv 2
# n = number of iterations to run for bootstrapping
# frac = fraction of data points to use to make bootstrapped hvs
#    ** can also be input as a whole number for how many  
#       data points to use in resampling of bootstraps
# returns df with overlap and amount unique of each hv
bootOverlap = function(data1,data2,n,frac){
  require(tidyverse)
  df = tibble(sorenson = rep(NA, n), unique1 = NA, unique2 = NA)
  n1 = ceiling(frac*nrow(data1))
  if(frac>1){
    n1 = frac
  }
  n2 = ceiling(frac*nrow(data2))
  if(frac>1){
    n2 = frac
  }
  for (i in 1:num){
    rand1 = sample_n(data1, n1, replace = T)
    rand2 = sample_n(data2, n2, replace = T)
    require(hypervolume)
    randh1 = hypervolume_gaussian(rand1,
                                  samples.per.point = 5000,
                                  kde.bandwidth = estimate_bandwidth(rand1), 
                                  sd.count = 3, 
                                  quantile.requested = 0.95, 
                                  quantile.requested.type = "probability", 
                                  chunk.size = 1000, 
                                  verbose = F)
    randh2 = hypervolume_gaussian(rand2,
                                  samples.per.point = 5000,
                                  kde.bandwidth = estimate_bandwidth(rand2), 
                                  sd.count = 3, 
                                  quantile.requested = 0.95, 
                                  quantile.requested.type = "probability", 
                                  chunk.size = 1000, 
                                  verbose = F)
    set = hypervolume_set(randh1,randh2,check.memory = F, verbose = F)
    df$sorenson[i] = hypervolume_overlap_statistics(set)[2]
    df$unique1[i] = hypervolume_overlap_statistics(set)[3]
    df$unique2[i] = hypervolume_overlap_statistics(set)[4]
    
    if (i/num == .1){
      cat('10% done\n')
    } else if (i/num == 0.25){
      cat('25% done\n')
    }else if (i/num == 0.50){
      cat('50% done\n')
    }else if (i/num == 0.75){
      cat('75% done\n')
    }else if (i/num == 0.9){
      cat('90% done\n')
    }

  }
  return(df)
}


# species comparisons
# Alligator
hvAC= hypervolume_join(acSNhv, acTRhv)
AC = hypervolume_set(acSNhv, acTRhv,  check.memory = F)
hypervolume_overlap_statistics(AC)

ovAC = bootOverlap(acSN, acTR, 100, 0.75)
ovAC$comp = 'Species'
ovAC$type = 'ACS'

# McCormick
hvMC= hypervolume_join(mcSNhv, mcTRhv)
MC = hypervolume_set(mcSNhv, mcTRhv, check.memory = F, verbose = F)

hypervolume_overlap_statistics(MC)

ovMC = bootOverlap(mcSN, mcTR, 100, 0.75)
ovMC$comp = 'Species'
ovMC$type = 'MCS'

# bind data together
OV = bind_rows(ovAC, ovMC)
write_csv(OV, 'data/Ov_CI.csv')


# overlap metrics

# hypervolume_overlap_statistics(SN)
# jaccard      sorensen frac_unique_1 frac_unique_2 
# 0.02065068    0.04046571    0.97925928    0.17370518

# hypervolume_overlap_statistics(TR)
# jaccard      sorensen frac_unique_1 frac_unique_2 
# 0.02106570    0.04126218    0.97857912    0.44043478 

# hypervolume_overlap_statistics(AC)
# jaccard      sorensen frac_unique_1 frac_unique_2 
# 0.07033968    0.13143431    0.91689751    0.68586981 

# hypervolume_overlap_statistics(MC)
# jaccard      sorensen frac_unique_1 frac_unique_2 
# 0.07389729    0.13762451    0.90342465    0.76063307 

# plot hypervolumes----
# plot snook

tiff("figs/SnookHV_sys.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
plot(hvSN, pairplot = T, colors=c('darkolivegreen3','deepskyblue3'),
     names= c("Freshwater\n Benthic",
              "Marine\n Benthic\n Detritivore",
              "Marine\n Benthic\n Zoobenthivore",
              "Marine\n Pelagic",
              "Trophic\n Position"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=T, limits=c(-6,6), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)
dev.off()

# plot tarpon
tiff("figs/TarponHV_sys.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

plot(hvTR, pairplot = T, colors=c('darkolivegreen3','deepskyblue3'),
     names= c("Freshwater\n Benthic",
              "Marine\n Benthic\n Detritivore",
              "Marine\n Benthic\n Zoobenthivore",
              "Marine\n Pelagic",
              "Trophic\n Position"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=T, limits=c(-6,6), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)
dev.off()

# plot ACS
tiff("figs/ACSHV_sp.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
plot(hvAC, pairplot = T, colors=c('goldenrod','gainsboro'),
     names= c("Freshwater\n Benthic",
              "Marine\n Benthic\n Detritivore",
              "Marine\n Benthic\n Zoobenthivore",
              "Marine\n Pelagic",
              "Trophic\n Position"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=T, limits=c(-6,6), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)
dev.off()

# plot MCS
tiff("figs/MCSHV_sp.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')

plot(hvMC, pairplot = T, colors=c('goldenrod','gainsboro'),
     names= c("Freshwater\n Benthic",
              "Marine\n Benthic\n Detritivore",
              "Marine\n Benthic\n Zoobenthivore",
              "Marine\n Pelagic",
              "Trophic\n Position"),
     show.3d=FALSE,plot.3d.axes.id=NULL,
     show.axes=TRUE, show.frame=TRUE,
     show.random=T, show.density=TRUE,show.data=F,
     show.legend=T, limits=c(-6,6), 
     show.contour=F, contour.lwd= 2, 
     contour.type='alphahull', 
     contour.alphahull.alpha=0.25,
     contour.ball.radius.factor=1, 
     contour.kde.level=0.01,
     contour.raster.resolution=100,
     show.centroid=TRUE, cex.centroid=2,
     point.alpha.min=0.2, point.dark.factor=0.5,
     cex.random=0.5,cex.data=1,cex.axis=1.5,cex.names=2,cex.legend=2,
     num.points.max.data = 100000, num.points.max.random = 200000, reshuffle=TRUE,
     plot.function.additional=NULL,
     verbose=FALSE
)
dev.off()

