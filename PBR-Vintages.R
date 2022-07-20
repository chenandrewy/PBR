# 2022 05: 

# compare main scatter using no factors, ff3 factors, old ff3 factors
# compare new and old 

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
  
  
  # download FF (new)
  url = 'http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Research_Data_Factors_CSV.zip'
  download.file(url,'../data/deleteme.zip')
  unzip('../data/deleteme.zip', exdir = '../data')
  
  # dl from web.archive.org like this doesn't work
  # temp_url = 'https://web.archive.org/web/20011218003540/http://mba.tuck.dartmouth.edu/pages/faculty/ken.french/ftp/F-F_Benchmark_Factors_Monthly.zip'
  # download.file(temp_url, '../data/deleteme.zip')
  # unzip('../data/deleteme.zip', exdir = '../data')


} # if exists

## Import ====

# cz
doc =  fread("../data/SignalDoc.csv") %>%
  rename(signalname = Acronym) %>%
  select(-c(`Detailed Definition`, `Notes`))


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


  
## ff
ff_all = fread('../data/F-F_Research_Data_Factors.CSV')  %>% 
  mutate(vint = 2022) %>% 
  rbind(
    fread('data/FF-Vintage/F-F_Research_Data_Factors_2005.txt') %>% mutate(vint = 2005)
  ) %>% 
  rbind(
    fread('data/FF-Vintage/F-F_Research_Data_Factors_2012.txt') %>% mutate(vint = 2012)    
  )  
  
colnames(ff_all) = c('yearm', 'mktrf', 'smb', 'hml', 'rf', 'vint')

## Performance Measures ====

target_data = cz_all %>% 
  filter(vint == 2022 & insamp & port == 'LS') %>% 
  mutate(yearm = year(date)*100 + month(date)) %>% 
  select(signalname, yearm, ret)

# no adjustment
fit_raW = target_data[
  , list(
      alpha = summary(lm(ret~1))$coefficients['(Intercept)' , 'Estimate']
      , tstat = summary(lm(ret~1))$coefficients['(Intercept)' , 't value']
  )
  , by=signalname
] %>% 
  mutate(
    model = 'raw'
  )
 
## ff3 model 2022
fitme = target_data %>% 
  left_join(
    ff_all %>% filter(vint == 2022), by = 'yearm'
  )

fit_ff_2022 = fitme[
  , list(
    alpha = summary(lm(ret ~ mktrf+smb+hml))$coefficients['(Intercept)' , 'Estimate']
    , tstat = summary(lm(ret ~ mktrf+smb+hml))$coefficients['(Intercept)' , 't value']
  )
  , by=signalname
] %>% 
  mutate(
    model = 'ff3_2022'
  )

## ff3 model 2012
fitme = target_data %>% 
  left_join(
    ff_all %>% filter(vint == 2012), by = 'yearm'
  )
fit_ff_2012 = fitme[
  , list(
    alpha = summary(lm(ret ~ mktrf+smb+hml))$coefficients['(Intercept)' , 'Estimate']
    , tstat = summary(lm(ret ~ mktrf+smb+hml))$coefficients['(Intercept)' , 't value']
  )
  , by=signalname
] %>% 
  mutate(
    model = 'ff3_2012'
  )


## ff3 model 2005
fitme = target_data %>% 
  left_join(
    ff_all %>% filter(vint == 2005), by = 'yearm'
  )
fit_ff_2005 = fitme[
  , list(
    alpha = summary(lm(ret ~ mktrf+smb+hml))$coefficients['(Intercept)' , 'Estimate']
    , tstat = summary(lm(ret ~ mktrf+smb+hml))$coefficients['(Intercept)' , 't value']
  )
  , by=signalname
] %>% 
  mutate(
    model = 'ff3_2005'
  )

fit_all = fit_raW %>% rbind(fit_ff_2012) %>% rbind(fit_ff_2022) %>% rbind(fit_ff_2005)

# HML revisions rep ====

ff_all %>% 
  group_by(vint) %>% 
  summarize_at(.vars =vars(mktrf,smb,hml), sd)

temp = ff_all %>% 
  select(yearm,vint,hml) %>% 
  pivot_wider(
    names_from = vint, values_from = hml, names_prefix = 'vint'
  ) %>% 
  mutate(
    rev = vint2022- vint2012, time = floor(yearm/100) + (yearm/100 - floor(yearm/100))/0.12
  ) %>% 
  filter(
    !is.na(rev)
  )


ggplot(temp) +
  geom_line(aes(x=time,y=rev)) +
  theme_minimal(
    base_size = 20
  ) +
  ylab('HML revision (ppt)') +
  annotate(
    geom = 'text', x = 1980, y = 4.2
    , label = paste0(
      'SD(revision) = ', round(sd(temp$rev),3), '\n'
      , 'SD(HML 2022) = ', round(sd(temp$vint2022), 3)
    )
    , size = 6
  ) +
  xlab(NULL)



ggsave('../results/hml-rev.png', width = 8, height = 5)


# Name Scatter cleaner ====

ablines = tibble(slope = 1, 
                 intercept = 0,
                 group = factor(x = c('45 degree line'),
                                levels = c('45 degree line')))


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
    # plot.title = element_text(text_family = "sans", face="bold"),
    # plot.subtitle = element_text(family = "sans"),
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


