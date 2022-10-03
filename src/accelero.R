################################################################################
#
#
# ANALYSIS OF ACCELEROMETRY DATA
# 
# ----------------------
# 
# accelero.R - Simon Lacombe - Accelerometry
#
#
################################################################################

rm(list = ls())

library(foreach)
library(tidyverse)
library(lubridate)
library(stringi)
library(zoo)
library(depmixS4)

source("src/get_features.R")

window = 60                                                                      # time window used to get static acceleration
res_analysis = 60                                                                # resolution used for the analysis (HMM time step)
files <- 1:24

####### LOAD DATA --------------------------------------------------------------
boar_dat.path <- "datafile/acc_sanglier_33395"
boar_dat.filenames <- list.files(path = boar_dat.path, pattern = ".rds")
boar_dat.filenames <-boar_dat.filenames[files]

boar_dat <- foreach(file = boar_dat.filenames, .combine = rbind) %do%{
  readRDS(paste0(boar_dat.path,"/",file))
}

res_orig <- boar_dat$ms[nrow(boar_dat)]/(nrow(boar_dat)-1) / 1000               # original resolution (s)

####### PROCESS DATA -----------------------------------------------------------
dat_processed.l <- get_features(boar_dat,
                                window = window,
                                res_orig = res_orig,
                                res_analysis = res_analysis )

boar_dat <- dat_processed.l[[1]]
boar_dat_seg <- dat_processed.l[[2]]
    

####### FIT HMM FROM DEPMIXS4 --------------------------------------------------

features_ok <- c("Dps.x", "ODBA.mn", "Dynamic.x.mn")

boar_dat_seg[,features_ok]%>%
  cor%>%
  corrplot::corrplot(method = "color", type = "upper")

features.l <- paste(features_ok, "~ 1")

features.l <- lapply(features.l, function(x){eval(parse(text = x))})

model <- depmix(features.l,
                   data = boar_dat_seg,
                   family = rep(list(gaussian()), length(features.l)),
                   ntimes = nrow(boar_dat_seg), 
                   nstates = 2)
HMM_boar <- fit(model)

summary(HMM_boar)

posterior(HMM_boar)$state

cbind(boar_dat_seg, data_frame(est.state = posterior(HMM_boar)$state))%>%
  pivot_longer(cols = c("Dynamic.x.mn","est.state"),
               names_to = "VAR", values_to = "accl")%>%
  ggplot()+
  geom_line(aes(x=time,y=accl, color = VAR))+
  facet_grid(VAR~., scale = "free")

