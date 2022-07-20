# 2022 07: 
# Generates figure comparing Chen & Zimmerman's reproduced t-stats and original paper t-stats 


# SETUP ====

rm(list = ls())

## Globals ====
library(tidyverse)
library(data.table)
library(googledrive)
library(gridExtra)
library(nlme)

dir.create('../data/')


# use this for original papers
SUBDIR = 'Full Sets OP'; FILENAME = 'PredictorPortsFull.csv'

## Download ====

if (!file.exists('../data/new-PredictorPortsFull.csv')){

  # March 2022 cz data
  url <- "https://drive.google.com/drive/folders/1O18scg9iBTiBaDiQFhoGxdn4FdsbMqGo" 

  url %>% drive_ls() %>%
    filter(name == "Portfolios") %>% drive_ls() %>% 
    filter(name == SUBDIR) %>% drive_ls() %>% 
    filter(name == FILENAME) %>% 
    drive_download(path = paste0("../data/new-",FILENAME), overwrite = TRUE)
  
  # signal doc 
  url %>% drive_ls() %>% 
    filter(name == "SignalDoc.csv") %>% 
    drive_download(path = "../data/SignalDoc.csv", overwrite = TRUE)
}


## Import ====

# All signals
doc =  fread("../data/SignalDoc.csv") %>%
  rename(signalname = Acronym) %>%
  select(-c(`Detailed Definition`, `Notes`))


# Chen-Zimmerman dataset
cz_all = fread(paste0("../data/new-",FILENAME)) %>%
  mutate(vint = 2022) %>%
  rbind(
    fread(paste0("../data/new-",FILENAME)) %>%
      mutate(vint = 2021)
  ) %>%
  select(vint, signalname, port, date, ret) %>%
  left_join(
    doc %>% select(signalname, SampleStartYear, SampleEndYear)
    , by = 'signalname'
  ) %>%
  mutate(
    insamp = (year(date) >= SampleStartYear) &  (year(date) <= SampleEndYear)
  )


## Performance Measures ====

target_data = cz_all %>% 
  filter(vint == 2022 & insamp & port == 'LS') %>% 
  mutate(yearm = year(date)*100 + month(date)) %>% 
  select(signalname, yearm, ret)

# no adjustment
fit_all = target_data[
  , list(
      alpha = summary(lm(ret~1))$coefficients['(Intercept)' , 'Estimate']
      , tstat = summary(lm(ret~1))$coefficients['(Intercept)' , 't value']
  )
  , by=signalname
] %>% 
  mutate(
    model = 'raw'
  )


# Name Scatter cleaner ====


# select comparable t-stats
fit_OP = doc %>% 
  mutate(
    tstat_OP = abs(as.numeric(`T-Stat`))
  ) %>% 
  select(
    signalname, tstat_OP, `Predictability in OP`, `Signal Rep Quality`, `Test in OP`
  ) %>% 
  filter(
    `Signal Rep Quality` %in% c('1_good','2_fair')
    , grepl('port', `Test in OP`)
    , `Predictability in OP` != 'indirect'
  ) 

# merge 
fitcomp = fit_all %>% 
  filter(model == 'raw') %>% 
  rename(tstat_CZ = tstat) %>% 
  inner_join(
    fit_OP
    , by = 'signalname'
  ) %>% 
  filter(!is.na(tstat_OP)) %>%  # some port sorts have only point estimates 
  filter(tstat_CZ>0)  # for activism you can get a negative ff alpha

# plot
ablines = tibble(slope = 1, 
                 intercept = 0,
                 group = factor(x = c('45 degree line'),
                                levels = c('45 degree line')))

fitcomp %>% 
  ggplot(aes(y=tstat_CZ, x = tstat_OP)) +
  geom_line(aes(y=tstat_CZ, x= tstat_OP), colour=NA) +
  geom_smooth(method = lm, se = FALSE, color = "#666666", size=.6) +
  geom_point(size=3) +
  coord_trans(x='log10', y='log10', xlim = c(1.5, 17), ylim = c(1.0, 15)) +
  ggrepel::geom_text_repel(aes(label=signalname), max.overlaps = Inf, box.padding = 0.5) +
  scale_x_continuous(breaks=c(2, 5, 10, 15)) +
  scale_y_continuous(breaks=c(2, 5, 10, 15)) +
  theme_minimal(
    base_size = 20
  ) +
  theme(
    legend.position = c(.9, .1), 
    legend.title = element_blank(),

    axis.title.x = element_text(family = "sans", size = 12, face="bold"),
    axis.title.y = element_text(family = "sans", size = 12, face="bold"),
    axis.text.x = element_text(family = "sans", size = 12),
    axis.text.y = element_text(family = "sans", size = 12),
    plot.caption = element_text(family = "sans", hjust = 0),
  ) +
  geom_abline(
    data = ablines, aes(slope = slope, intercept = intercept, linetype = group)
    , color = 'black', size = 1
  ) +
  labs(y = 't-stat reproduction', 
       x = 't-stat original paper')  


