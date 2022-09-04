# 2022 07 Alec Erb on Behalf of Andrew Chen: Data download for PBR


# ENVIRONMENT ====

rm(list = ls())
library(tidyverse)
library(data.table)
library(googledrive)
library(lubridate)


dir.create('../data/')

# March 2022 cz data
url <- "https://drive.google.com/drive/folders/1O18scg9iBTiBaDiQFhoGxdn4FdsbMqGo" 


## Download Baseline CZ Data ====

if (!file.exists('../data/PredictorPortsFull.csv')){
  
  SUBDIR = 'Full Sets OP'; FILENAME = 'PredictorPortsFull.csv'
  
  url %>% drive_ls() %>%
    filter(name == "Portfolios") %>% drive_ls() %>% 
    filter(name == SUBDIR) %>% drive_ls() %>% 
    filter(name == FILENAME) %>% 
    drive_download(path = paste0("../data/",FILENAME), overwrite = TRUE)
  
  # signal doc 
  url %>% drive_ls() %>% 
    filter(name == "SignalDoc.csv") %>% 
    drive_download(path = "../data/SignalDoc.csv", overwrite = TRUE)
}



## Download Liquidity Tests ====
# only used for one figure

if (!file.exists('../data/PredictorAltPorts_LiqScreen_VWforce.csv')){
  
  FILENAME = 'PredictorAltPorts_LiqScreen_VWforce.zip'
  
  url %>% drive_ls() %>%
    filter(name == "Portfolios") %>% drive_ls() %>% 
    filter(name == 'Full Sets Alt') %>% drive_ls() %>% 
    filter(name == FILENAME) %>% 
    drive_download(path = paste0("../data/deleteme.zip"), overwrite = TRUE)
  unzip('../data/deleteme.zip', exdir = '../data')
  
  FILENAME = 'PredictorAltPorts_LiqScreen_ME_gt_NYSE20pct.zip'
  
  url %>% drive_ls() %>%
    filter(name == "Portfolios") %>% drive_ls() %>% 
    filter(name == 'Full Sets Alt') %>% drive_ls() %>% 
    filter(name == FILENAME) %>% 
    drive_download(path = paste0("../data/deleteme.zip"), overwrite = TRUE)
  unzip('../data/deleteme.zip', exdir = '../data')  
  
  
  FILENAME = 'PredictorAltPorts_HoldPer_6.zip'
  
  url %>% drive_ls() %>%
    filter(name == "Portfolios") %>% drive_ls() %>% 
    filter(name == 'Full Sets Alt') %>% drive_ls() %>% 
    filter(name == FILENAME) %>% 
    drive_download(path = paste0("../data/deleteme.zip"), overwrite = TRUE)
  unzip('../data/deleteme.zip', exdir = '../data')
  
  FILENAME = 'PredictorAltPorts_HoldPer_12.zip'
  
  url %>% drive_ls() %>%
    filter(name == "Portfolios") %>% drive_ls() %>% 
    filter(name == 'Full Sets Alt') %>% drive_ls() %>% 
    filter(name == FILENAME) %>% 
    drive_download(path = paste0("../data/deleteme.zip"), overwrite = TRUE)
  unzip('../data/deleteme.zip', exdir = '../data')  
  
  file.remove('../data/deleteme.zip')
  
} # end if file exists



## import YZ ====
# only used for one figure
# (YZ data is not public, can't auto download)


if (!file.exists('../data/yz_sum.csv')){
  temp = read_sas('../data-YZ/Yan_Zheng_RFS_Data.sas7bdat')
  
  yzsum = temp %>%
    mutate(
      signalname = paste(transformation, fsvariable, sep = '.')
    ) %>%
    transmute(
      signalname, date = DATE, ret = 100*ddiff_ew
    ) %>% 
    filter(!is.na(ret)) %>% 
    group_by(signalname) %>% 
    summarize(
      rbar = mean(ret)
      , vol = sd(ret)
      , tstat = mean(ret)/sd(ret)*sqrt(dplyr::n()) %>% abs()
    ) %>%
    mutate(tstat = abs(tstat)) %>%
    write.csv(.,file = "../data/yz_sum.csv")
}