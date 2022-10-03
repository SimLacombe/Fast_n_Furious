################################################################################
#
#
# THIS FILE CONTAINS FUNCTIONS TO GET SEVERAL FEATURES FROM AN
# ACCELEROMETRY DATASET
# 
# ----------------------
# 
# get_features.R - Simon Lacombe - Accelerometry
#
# To be launched from accelero.R
#
################################################################################


#######get_features ------------------------------------------------------------

get_features <- function(dat, window, res_orig, res_analysis){

  # dat : accelerometry dataset :
  #     dat$time : timestamp of the observation
  #     dat$ms : time in ms since begining of the data series
  #     dat$ms_to_on : time in ms since the deployment of the colar
  #     dat$x/y/z : acceleration in the x/y/z direction
  #     
  # window : time window used to calculate static acceleration
  #
  # res_orig : resolution of the original dataset
  #
  # res_analysis : resultion used to further conduct the analysis 
  #                (resolution of the output dataset)
  
#1. Get static and Dynamic acceleration ------
  dat <- dat%>%
    mutate(Raw.x = x, Raw.y = y, Raw.z = z,
           Static.x = rollmean(x, k = floor(window/res_orig), align = "center", fill = NA),
           Static.y = rollmean(y, k = floor(window/res_orig), align = "center", fill = NA),
           Static.z = rollmean(z, k = floor(window/res_orig), align = "center", fill = NA),
           Dynamic.x = x - Static.x,
           Dynamic.y = y - Static.y,
           Dynamic.z = z - Static.z)%>%
    filter(!is.na(Static.x))

#2. Get ODBA, VeDBA ------
  dat <- dat %>%
    mutate(ODBA = abs(Dynamic.x) + abs(Dynamic.y) + abs(Dynamic.z),
           VeDBA = sqrt(Dynamic.x**2 + Dynamic.y**2 + Dynamic.z**2)) 
  
#3. Change resolution and get mean and sd of the previous features ------
  dat_seg <- dat %>% 
    mutate(cut = ms_to_on%/%(res_analysis*1000)) %>%
    group_by(cut)%>%
    summarise(time=time[1],
              nsample = n(),
              Raw.x.mn = mean(Raw.x),
              Raw.y.mn = mean(Raw.y),
              Raw.z.mn = mean(Raw.z),
              Raw.x.sd = sd(Raw.x),
              Raw.y.sd = sd(Raw.y),
              Raw.z.sd = sd(Raw.z),
              Static.x.mn = mean(Static.x),
              Static.y.mn = mean(Static.y),
              Static.z.mn = mean(Static.z),
              Static.x.sd = sd(Static.x),
              Static.y.sd = sd(Static.y),
              Static.z.sd = sd(Static.z),
              Dynamic.x.mn = mean(Dynamic.x),
              Dynamic.y.mn = mean(Dynamic.y),
              Dynamic.z.mn = mean(Dynamic.z),
              Dynamic.x.sd = sd(Dynamic.x),
              Dynamic.y.sd = sd(Dynamic.y),
              Dynamic.z.sd = sd(Dynamic.z),
              ODBA.mn = mean(ODBA),
              VeDBA.mn = mean(VeDBA),
              ODBA.sd = sd(ODBA),
              VeDBA.sd = sd(VeDBA))
  
#4. Get spectral features ------
  dat_seg <- dat %>% 
    mutate(cut = ms_to_on%/%(res_analysis*1000))%>%
    group_by(cut)%>%
    summarise(Dps.x = get_spectrum(Dynamic.x, res_analysis)[[1]],
              FDps.x = get_spectrum(Dynamic.x, res_analysis)[[2]],
              Dps.y = get_spectrum(Dynamic.y, res_analysis)[[1]],
              FDps.y = get_spectrum(Dynamic.y, res_analysis)[[2]],
              Dps.z = get_spectrum(Dynamic.z, res_analysis)[[1]],
              FDps.z = get_spectrum(Dynamic.z, res_analysis)[[2]])%>%
    merge(dat_seg, .,  by = "cut")
  
  return(list(dat, dat_seg))
}

#######get spectrum ------------------------------------------------------------

get_spectrum <- function(x, T) {
  
  # Get two spectral features from a 1-D acceleration sample
  # ---
  # x : 1-D acceleration sample
  # T : length (in s) of the acceleration sample. Fixed to res_analysis
  #     in get_features() 
  #
  
  n <- length(x)
  freq <- c(0:(length(x)/2-1))/T
  
  fftx  <- Mod(fft(x))/n
  fftx[2:n] <- 2*fftx[2:n]
  fftx <- fftx[1:floor(n/2)]
  
  Dps <- max(fftx)
  FDps <- freq[which(fftx == Dps)]
  return(list(Dps, FDps))
}