# 2022 07 Alec Erb on Behalf of Andrew Chen: Data download for PBR


# ENVIRONMENT ====

rm(list = ls())
library(tidyverse)
library(data.table)
library(googledrive)
library(lubridate)

dir.create('./data/')

SUBDIR = 'Full Sets OP'; FILENAME = 'PredictorPortsFull.csv'


## Download ====

if (!file.exists('../data/new-PredictorPortsFull.csv')){
  
  # March 2022 cz data
  url <- "https://drive.google.com/drive/folders/1O18scg9iBTiBaDiQFhoGxdn4FdsbMqGo" 
  
  url %>% drive_ls() %>%
    filter(name == "Portfolios") %>% drive_ls() %>% 
    filter(name == SUBDIR) %>% drive_ls() %>% 
    filter(name == FILENAME) %>% 
    drive_download(path = paste0("./data/",FILENAME), overwrite = TRUE)
  
  # signal doc 
  url %>% drive_ls() %>% 
    filter(name == "SignalDoc.csv") %>% 
    drive_download(path = "./data/SignalDoc.csv", overwrite = TRUE)
}


