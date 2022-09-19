# libraries ---------------------------------------------------------------


library(tidyverse)

library(data.table)
library(googledrive)
library(ggplot2)
library(lubridate)
library(zoo)
library(extrafont)
library(latex2exp)

dir.create('../results/')

# globals ---------------------------------------------------------------



MATBLUE = rgb(0,0.4470,0.7410)
MATRED = rgb(0.8500, 0.3250, 0.0980)
MATYELLOW = rgb(0.9290, 0.6940, 0.1250)

NICEBLUE = "#619CFF"
NICEGREEN = "#00BA38"
NICERED = "#F8766D"

chen_theme =   theme_minimal() +
  theme(
    text = element_text(family = "Palatino Linotype")
    , panel.border = element_rect(colour = "black", fill=NA, size=1)
    
    # Font sizes
    , axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.text = element_text(size = 18),
    
    # Tweaking legend
    legend.position = c(0.7, 0.8),
    legend.text.align = 0,
    legend.background = element_rect(fill = "white", color = "black"),
    legend.margin = margin(t = 5, r = 20, b = 5, l = 5), 
    legend.key.size = unit(1.5, "cm")
    , legend.title = element_blank()    
  ) 


# shared data cleaning ----------------------------------------------------

# read from disc
cz_all = fread("../data/PredictorPortsFull.csv")
signaldoc = fread('../data/SignalDoc.csv')
yzsum = fread('../data/yz_sum.csv')


signaldoc = fread('../data/SignalDoc.csv') %>% 
  rename(signalname = Acronym) %>% 
  mutate(
    pubdate = as.Date(paste0(Year, '-12-31'))
    , sampend = as.Date(paste0(SampleEndYear, '-12-31'))
    , sampstart = as.Date(paste0(SampleStartYear, '-01-01'))
  ) %>% 
  select(-c(Notes,`Detailed Definition`))

# czret (monthly returns)
czret = cz_all %>%                                         
  filter(!is.na(ret), port == 'LS') %>%                                                           
  left_join(signaldoc) %>% 
  mutate(
    samptype = case_when(
      (date >= sampstart) & (date <= sampend) ~ 'in-samp'
      , (date > sampend) & (date <= sampend %m+% months(36)) ~ 'out-of-samp'
      , (date > pubdate) ~ 'post-pub'
      , TRUE ~ 'NA_character_'
    )
  ) %>% 
  select(signalname, date, ret, samptype, sampstart, sampend)


# cz sum (sum stats by signal-samptype)
czsum = czret %>%
  group_by(signalname, samptype) %>%
  summarize(
    rbar = mean(ret)
    , vol = sd(ret)
    , tstat = mean(ret)/sd(ret)*sqrt(dplyr::n())
  )


# czretmat (return matrix for bootstraps)
# used in correlations and one of the bootstraps
# use full sample for most overlap
czretmat = czret %>% 
  select(signalname, date, ret) %>% 
  pivot_wider(names_from = signalname, values_from = ret) %>%
  select(-date) %>% 
  as.matrix()




# for liquidity screens
cz_alt = rbind(
  fread('../data/PredictorAltPorts_HoldPer_12.csv') %>% 
    mutate(group = 'holdper12')
  ,   fread('../data/PredictorAltPorts_HoldPer_6.csv') %>% 
    mutate(group = 'holdper06')
  , fread('../data/PredictorAltPorts_LiqScreen_ME_gt_NYSE20pct.csv') %>% 
    mutate(group = 'mescreen')  
  ,   fread('../data/PredictorAltPorts_LiqScreen_VWforce.csv') %>% 
    mutate(group = 'VWforce')
) %>% 
  filter(port == 'LS', !is.na(ret)) %>% 
  select(group, signalname, date, ret) %>% 
  left_join(signaldoc) %>% 
  mutate(
    samptype = case_when(
      (date >= sampstart) & (date <= sampend) ~ 'in-samp'
      , (date > sampend) & (date <= sampend %m+% months(36)) ~ 'out-of-samp'
      , (date > pubdate) ~ 'post-pub'
      , TRUE ~ 'NA_character_'
    )
  )


