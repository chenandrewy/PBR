# 2022 07 Alec Erb on behalf of Andrew Chen: Generate Simple Figures for PBR

# SETUP ====

rm(list = ls())

## Globals ====
library(tidyverse)
library(data.table)
library(googledrive)
library(gridExtra)
library(nlme)

SUBDIR = 'Full Sets OP'; FILENAME = 'PredictorPortsFull.csv'

dir.create('./results/')

# Generate McLean-Pontiff bootstrapped mean returns figure ====
# read in PredictorPortsFull.csv
ret0 = fread('data/PredictorPortsFull.csv') %>%                                           
  filter(!is.na(ret), port == 'LS')

# read in SignalDoc.csv and isolates to relevant date cols
signaldoc = fread('data/SignalDoc.csv') %>% 
  mutate(
    signalname = Acronym
    , pubdate = as.Date(paste0(Year, '-12-31'))
    , sampend = as.Date(paste0(SampleEndYear, '-12-31'))
    , sampstart = as.Date(paste0(SampleStartYear, '-01-01'))
  ) %>% 
  arrange(signalname) %>% 
  select(signalname, pubdate, sampend, sampstart)


# merge and find OOS returns
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


## bootstrap mean distributions ====
# clustered by month

set.seed(6)

nboot = 500

# make wide dataset, use NA if not correct sample

bootfun = function(sampname){
  
  # make wide dataset, use NA if not correct sample
  wide_is = ret1 %>%
    filter(samptype == sampname) %>% 
    pivot_wider(
      c(date, "ret", signalname), names_from = "signalname", values_from = "ret"
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
rboot1 = bootfun('in-samp')
rboot2 = bootfun('out-of-samp')

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


ggsave(
  "results/MPrep.pdf",
  width = 12,
  height = 12
)



# Generate R2 replication figure ====

## Import ====

# All signals
doc =  fread("../data/SignalDoc.csv") %>%
  rename(signalname = Acronym) %>%
  select(-c(`Detailed Definition`, `Notes`))


# Chen-Zimmerman dataset
cz_all = fread(paste0("../data/",FILENAME)) %>%
  mutate(vint = 2022) %>%
  rbind(
    fread(paste0("../data/",FILENAME)) %>%
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

## Name Scatter cleaner ====

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

ggsave(
  "results/fitcomp.pdf",
  width = 12,
  height = 12
)





# Generate t too big figure ====

## settings ====

# for nboot = 1000, takes 20 seconds.  
# nboot = 10,000 takes 3 minutes.  This is about 20 million predictors
library(xtable) 

nboot = 10*1000
ndatemax = dim(emat)[1]
ndatemin = 240
tedge = c(seq(0,6,1), 20)

## Import ====

# All signals
doc =  fread("../data/SignalDoc.csv") %>%
  rename(signalname = Acronym) %>%
  select(-c(`Detailed Definition`, `Notes`))


# Chen-Zimmerman dataset
cz_all = fread(paste0("../data/PredictorPortsFull.csv"))  %>% 
  select(signalname, port, date, ret) %>%
  left_join(
    doc %>% select(signalname, SampleStartYear, SampleEndYear)
    , by = 'signalname'
  ) %>%
  mutate(
    insamp = (year(date) >= SampleStartYear) &  (year(date) <= SampleEndYear)
  )

# Performance Measures 
target_data = cz_all %>% 
  filter(insamp & port == 'LS') %>% 
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

## Bootstrap ====

# residuals
e = cz_all %>% 
  filter(insamp, port == 'LS') %>% 
  group_by(signalname) %>% 
  mutate(
    e = as.vector(scale(ret, center = T, scale = F))
  ) %>% 
  select(signalname, date, e) %>% 
  pivot_wider(names_from = signalname, values_from = e) %>% 
  select(-date) %>% 
  as.matrix()

# bootstrap
# fast binding based on https://stackoverflow.com/questions/19697700/how-to-speed-up-rbind
boot = data.table()

boot_once = function(){
  dateselect = sample(1:dim(emat)[1], ndatemax, replace = T)
  
  # draw returns, clustered by month
  eboot = e[dateselect, ]
  
  # summarize each predictor
  ebar = apply(eboot, 2, mean, na.rm=T)
  vol  = apply(eboot, 2, sd, na.rm=T)
  ndate  = apply(eboot, 2, function(x) sum(!is.na(x)))
  tstat = ebar/vol*sqrt(ndate)
  
  # remove if not enough observations
  tstat2 = tstat[ndate > ndatemin]
  
  # summarize across predictors
  hdat = hist(abs(tstat2), tedge, plot = F)
  temp = data.table(t_left = tedge[1:(length(tedge)-1)], p = hdat$counts/length(tstat2))
  return(temp)
}


# bootstrap!!

tic = Sys.time()
set.seed(1057)
bootdat = rbindlist(lapply(1:nboot, function(x) boot_once()))
toc = Sys.time()
toc - tic


# summarize across bootstraps
bootsum = bootdat %>% 
  group_by(t_left) %>% 
  summarize(
    p = mean(p)
  )


## Make Table ====

# empirical cdf
pemp = ecdf(fit_all$tstat)

# add bootstrap cdf??

# table settings
t_left  = tedge[1:(length(tedge)-1)]
t_right = tedge[2:length(tedge)]

# table
tab_too_big = data.frame(
  t_left 
  , t_right
  , N_emp = length(fit_all$tstat)*(pemp(t_right) - pemp(t_left))
  , prob_emp = pemp(t_right) - pemp(t_left)
  , prob_norm = 2*(pnorm(t_right) - pnorm(t_left))
) %>% 
  left_join(
    bootsum %>% rename(prob_boot = p)
    , by = 't_left'
  ) %>% 
  mutate(
    emp_to_norm = prob_emp / prob_norm
    , N_hack_norm = N_emp/prob_norm
  )



tab_too_big_wide = tab_too_big %>% t()

tab_too_big_wide

print(xtable(tab_too_big_wide, type = "latex"), file = "../results/tab-too-big.tex")
