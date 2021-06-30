#'""" Example code to make hv with raw/individual mm data
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

# Function to make random points of n length----
# data from random sample with mean and sd but 
# can be generated between a high and low value if end_points = T
# *** Note chose either column names or column numbers for ID_rows and names
# either work but must be the same
# df = dataframe or tibble with each row containing 
#        unique entry for making random points 
# ID_rows = vector of column names or numbers with id information
# names = column name or number of name of measure variables
# mean = column name or column number of df with mean 
# sd = column name or column number of df with sd 
# n = number of points to randomly generate
# z_score = T or F, if T z-scores values
# end_points = T or F for if random points need to be generated between
#        a high and low end point (e.g. 5% and 95% interval)
#        low and high required if end_points = T
# low = column name or column number of df with lower bound to sample in
# high = column name or column number of df with upper bound to sample in
HVvalues = function(df, ID_rows, names, mean, sd, n, z_score = F,
                    end_points = F, low = NULL, high = NULL){
  require(tidyverse)
  require(truncnorm)
  
  # check to see if information is needed to restrict where points are
  if (end_points){
    if (is_empty(df[,low]) | is_empty(df[,high])){
      return(cat('Warning: low and/or high columns not specified \n
                  Specific and run again or end_points = F \n'))
    }
  }
  
  # check to see if there are more 
  if(T %in% duplicated(df[,c(ID_rows,names)])){
    return(cat('Warning: some of the rows contain duplicated information \n
                make sure data is correct \n'))
  }
  
  # rename variables to make code work
  if (is.numeric(mean)){
    names(df)[mean] = 'Mean'
  }else {
    df = df %>% rename(Mean = mean)
  }
  
  if (is.numeric(sd)){
    names(df)[sd] = 'SD'
  }else {
    df = df %>% rename(SD = sd)
  }
  
  if (end_points){
    if (is.numeric(low)){
      names(df)[low] = 'lower'
    }else {
      df = df %>% rename(lower = low)
    }
    
    if (is.numeric(high)){
      names(df)[high] = 'upper'
    }else {
      df = df %>% rename(upper = high)
    }
  }
  
  # make sure the names is not numeric 
  if (is.numeric(names)){
    names = names(df)[names]
  }
  
  # generate random points within bounds
  if(end_points){
    
    df_tot = df %>% slice(rep(1:n(), each=n))%>% 
      mutate(point = 
               truncnorm::rtruncnorm(1, a = lower, b = upper,
                                     mean = Mean, sd = SD),
             num = rep(1:n, times=nrow(df))) %>%
      select(-Mean, -SD, -lower, -upper)%>%
      pivot_wider(names_from = names, values_from = point)%>% 
      select(-num)
  }else {
    # generate random points outside of bounds
    df_tot = df %>% slice(rep(1:n(), each=n))%>%
      mutate(point = 
               truncnorm::rtruncnorm(1, mean = Mean, sd = SD),
             num = rep(1:n, times=nrow(df))) %>%
      select(-Mean, -SD)%>%
      pivot_wider(names_from = names, values_from = point)%>% 
      select(-num)
  }
  if (z_score){
    df_tot = df_tot %>% 
      mutate_if(is.numeric, scale)
  }
  
  return(df_tot)
  
}

# Oyster reference and restored----
# mixing model results 
ref = mixTable(file = 'C:/Users/wrjam/Dropbox/WorkDocs/R/Github/ExCode/Hypervolumes/data/OysterRes_ss.txt', 
              site = 'Reference',ind = F, nest = F)
res = mixTable(file = 'C:/Users/wrjam/Dropbox/WorkDocs/R/Github/ExCode/Hypervolumes/data/OysterRef_ss.txt', 
               site = 'Restored', ind = F, nest = F)

# bind data together
df = bind_rows(ref, res) %>% 
  select(name, site, source, mean, sd, lowend, highend)


# number or iterations
n = 50


# create tibble to store data
df_hvAll = tibble(rep = rep(seq(1,n,1)),
                  ref_vol = NA, res_vol= NA,
                  sorenson = NA, randOV = NA)



for (i in 1:n){
  # generate random from mean and sd that are z-scored
  oy = HVvalues(df = df, ID_rows = c('name','site'),names = c('source'), 
                mean = 'mean', sd = 'sd', n = 20,
                end_points = T, low = 'lowend', high = 'highend', 
                z_score = T)
  
  # Reference hypervolume
  # subset reference data
  ref_df = oy %>% filter(site == 'Reference')%>%
    select(-name, -site)
  
  # generate hypervolume
  ref_hv = hypervolume_gaussian(ref_df, name = 'Reference',
                                #samples.per.point = ceiling((10^(3 + sqrt(ncol(ref_df))))/nrow(ref_df)),
                                samples.per.point = 5000,
                                kde.bandwidth = estimate_bandwidth(ref_df), 
                                sd.count = 3, 
                                quantile.requested = 0.95, 
                                quantile.requested.type = "probability", 
                                chunk.size = 1000, 
                                verbose = F)
  
  # get reference volume
  df_hvAll$ref_vol[i] = get_volume(ref_hv)

  
  # subset restoration data
  res_df = oy %>% filter(site == 'Restored') %>% 
    select(-name, -site)
  
  # generate restored hypervolume
  res_hv = hypervolume_gaussian(res_df, name = 'Restored',
                                #samples.per.point = ceiling((10^(3 + sqrt(ncol(res_df))))/nrow(res_df)),
                                samples.per.point = 5000,
                                kde.bandwidth = estimate_bandwidth(res_df), 
                                sd.count = 3, 
                                quantile.requested = 0.95, 
                                quantile.requested.type = "probability", 
                                chunk.size = 1000, 
                                verbose = F)
  # get restoration volume
  df_hvAll$res_vol[i] = get_volume(res_hv)
  
  # join hypervolumes
  set = hypervolume_set(ref_hv, res_hv, check.memory = F, verbose = F)
  
  # calculate sorenson overlap
  df_hvAll$sorenson[i] = hypervolume_overlap_statistics(set)[2]
  
  # create random overlap from data
  # randomly assign data to make hvs
  hvs = bind_rows(ref_df, res_df)
  sub = sample(nrow(hvs),nrow(res_df))
  df1 = hvs[sub,]
  df2 = hvs[-sub,]
  
  # generate random hv 1 
  hv1 = hypervolume_gaussian(df1, name = 'hv1',
                             #samples.per.point = ceiling((10^(3 + sqrt(ncol(df1))))/nrow(df1)),
                             samples.per.point = 5000,
                             kde.bandwidth = estimate_bandwidth(df1), 
                             sd.count = 3, 
                             quantile.requested = 0.95, 
                             quantile.requested.type = "probability", 
                             chunk.size = 1000, 
                             verbose = F)
  
  hv2 = hypervolume_gaussian(df2, name = 'hv2',
                             #samples.per.point = ceiling((10^(3 + sqrt(ncol(df2))))/nrow(df2)),
                             samples.per.point = 5000,
                             kde.bandwidth = estimate_bandwidth(df2), 
                             sd.count = 3, 
                             quantile.requested = 0.95, 
                             quantile.requested.type = "probability", 
                             chunk.size = 1000, 
                             verbose = F)
  # join hypervolumes
  set_rand = hypervolume_set(hv1, hv2, check.memory = F, verbose = F)
  
  # calculate overlap
  df_hvAll$randOV[i] = hypervolume_overlap_statistics(set_rand)[2]
  
}

# calculate percent recovery
df_hvALL = read_csv('oyster.csv') %>% 
  mutate(rec = sorenson/randOV)


# Reference hypervolume
# subset reference data
ref_df = oy %>% filter(site == 'Reference')%>%
  select(-name, -site)

# generate hypervolume
ref_hv = hypervolume_gaussian(ref_df, name = 'Reference',
                              #samples.per.point = ceiling((10^(3 + sqrt(ncol(ref_df))))/nrow(ref_df)),
                              samples.per.point = 5000,
                              kde.bandwidth = estimate_bandwidth(ref_df), 
                              sd.count = 3, 
                              quantile.requested = 0.95, 
                              quantile.requested.type = "probability", 
                              chunk.size = 1000, 
                              verbose = F)

# Restored hypervolume
# subset restoration data
res_df = oy %>% filter(site == 'Restored') %>% 
  select(-name, -site)

#generate restored hypervolume
res_hv = hypervolume_gaussian(res_df, name = 'Restored',
                              #samples.per.point = ceiling((10^(3 + sqrt(ncol(res_df))))/nrow(res_df)),
                              samples.per.point = 5000,
                              kde.bandwidth = estimate_bandwidth(res_df), 
                              sd.count = 3, 
                              quantile.requested = 0.95, 
                              quantile.requested.type = "probability", 
                              chunk.size = 1000, 
                              verbose = F)

# Plot oyerter data
hvOY= hypervolume_join(ref_hv, res_hv)

tiff("C:/Users/wrjam/Dropbox/WorkDocs/R/Github/ExCode/Hypervolumes/figs/HV_oy.tiff", width = 9, height = 9, units = 'in', 
     res = 600, compression = 'lzw')
plot(hvOY, pairplot = T, colors=c("#EBCC2A","#2ca1db"),
     names= c(expression(italic("Gracilaria"),"SPOM", "SSOM")),
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