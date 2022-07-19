# 2021 12 Andrew Chen: replicate the main result of McLean-Pontiff 

# ENVIRONMENT ====

rm(list = ls())
library(tidyverse)
library(data.table)
library(googledrive)
library(readxl)
library(RColorBrewer)
library(lubridate)

### USER ENTRY

# root of March 2022 release
pathRelease = 'https://drive.google.com/drive/folders/1O18scg9iBTiBaDiQFhoGxdn4FdsbMqGo'

# login to gdrive
# this prompts a login
pathRelease %>% drive_ls()

# create temporary directory 
dir.create('temp/')                                                             ## ?Where is this temp directory located?

# DOWNLOAD DATA =====

## reads in temp SignalDoc.csv to signaldoc df and isolates to relevant date cols ====
target_dribble = pathRelease %>% drive_ls() %>% 
  filter(name=='SignalDoc.csv')

drive_download(target_dribble, path = 'temp/deleteme.csv', overwrite = T)

signaldoc = fread('temp/deleteme.csv') %>% 
  mutate(
    signalname = Acronym
    , pubdate = as.Date(paste0(Year, '-12-31'))
    , sampend = as.Date(paste0(SampleEndYear, '-12-31'))
    , sampstart = as.Date(paste0(SampleStartYear, '-01-01'))
  ) %>% 
  arrange(signalname) %>% 
  select(signalname, pubdate, sampend, sampstart)





## Download all long-short returns (OP) ====

# use this for original papers
SUBDIR = 'Full Sets OP'; FILENAME = 'PredictorPortsFull.csv'

if (!file.exists('../temp/deleteme.csv')){
  
  # March 2022 cz data
  url <- "https://drive.google.com/drive/folders/1O18scg9iBTiBaDiQFhoGxdn4FdsbMqGo" 
  
  url %>% drive_ls() %>%
    filter(name == "Portfolios") %>% drive_ls() %>% 
    filter(name == SUBDIR) %>% drive_ls() %>% 
    filter(name == FILENAME) %>% 
    drive_download(target_dribble[1,], path = 'temp/deleteme.csv', overwrite = T)



ret0 = fread('temp/deleteme.csv') %>%                                           ## reshape returns data from wide to long and removes NaNs
  # pivot_longer(-c(date),names_to = 'signalname', values_to = 'ret') %>%
  filter(!is.na(ret))


# MERGE AND FIND OUT OF SAMPLE RETURNS ====
ret1 = ret0 %>%                                                                 ## merging, on signalname, returns and date (sampend, sampstart, pubdate) data
  left_join(signaldoc) %>% 
  mutate(
    samptype = case_when(
      (date >= sampstart) & (date <= sampend) ~ 'in-samp'
      , (date > sampend) & (date <= pubdate) ~ 'out-of-samp'
      , (date > pubdate) ~ 'post-pub'
      , TRUE ~ NA_character_
    )
  ) %>% 
  filter(!is.na(samptype))

sumsignal = ret1 %>%                                                            ## collapses on sampletype then signalname to find mean return and t-stat of each strategy in all three sample types
  group_by(signalname, samptype) %>% 
  summarize(
    rbar = mean(ret)
    , tstat = rbar/sd(ret)*sqrt(n())
  ) 

# remove bad reproductions / bad predictors (if desired)
signalok = sumsignal %>% 
  filter(samptype=='in-samp') %>% 
  filter(tstat > -Inf) %>% 
  transmute(signalname)
sumsignal = sumsignal %>% inner_join(signalok)

# check out of sample decay simple way ====
sumsamp = sumsignal %>% 
  group_by(samptype) %>% 
  summarize(rbar = mean(rbar), nsignal = n())

baseline = sumsamp %>% 
  filter(samptype == 'in-samp') %>% 
  select(rbar) %>% 
  as.numeric()

sumsamp = sumsamp %>%                                                           ## finds % decay of returns in out-of-samp and post-pub using in-samp as baseline 
  mutate(
    decay = (baseline - rbar)/baseline
  )


# BOOTSTRAP MEAN DISTRIBUTIONS ====
# clustered by month

set.seed(6)

nboot = 500

# make wide dataset, use NA if not correct sample

bootfun = function(sampname, portnum){
  
  # remove unnecessary portfolios
  ret1 <- ret1[ret1$port == portnum, ]
  
  # make wide dataset, use NA if not correct sample
  wide_is = ret1 %>%
    filter(samptype == sampname) %>% 
    pivot_wider(
      c(date, ret, signalname), names_from = signalname, values_from = ret
    ) %>% 
    select(-date) %>% 
    as.matrix()
  
  # make array that only has enough signals in each month (10)
  tgood = rowSums(!is.na(wide_is), na.rm=T) > 10
  mat = wide_is[tgood, ]
  T = dim(mat)[1]
  
  # bootstrap pooled mean
  rboot = rep(NA_real_, nboot)
  for (i in 1:nboot){
    tempt = sample(1:T, replace = T)
    rboot[i] = mat[tempt,]  %>% as.vector %>% mean(na.rm=T)
  }
  
  return(rboot)
  
} # end bootfun


# bootstrap for each sample type
rboot1 = bootfun('in-samp', "01")
rboot2 = bootfun('out-of-samp', "01")

# compile and plot
bootdat = data.frame(
  pooled_mean_ret = rboot1, samptype = 'in-samp' 
) %>% 
  rbind(
    data.frame(
      pooled_mean_ret = rboot2, samptype = 'out-of-samp' 
    )
  )


bootdat  %>% 
  ggplot(aes(x=pooled_mean_ret, fill=samptype)) +
  geom_histogram(
    alpha = 0.6, position = 'identity', breaks = seq(0,1,0.025), aes(y=..density..)
  ) +
  ggtitle('bootstrapped distribution') +
  labs(x='pooled mean return (% monthly)') +
  geom_vline(xintercept = 0)


