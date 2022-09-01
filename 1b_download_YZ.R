## Environment
library(tidyverse)
library(data.table)
library(googledrive)
library(lubridate)

## import YZ ====
if (!file.exists('../data/yz_sum.csv')){
  temp = read_sas('../data/Yan_Zheng_RFS_Data.sas7bdat')
  
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