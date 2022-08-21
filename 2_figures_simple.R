# 2022 07 Alec Erb on behalf of Andrew Chen: Generate Simple Figures for PBR

# SETUP ====

rm(list = ls())

## Globals ====
library(tidyverse)
library(data.table)
library(googledrive)
library(gridExtra)
library(nlme)
library(stargazer)
library(lubridate)


SUBDIR = 'Full Sets OP'; FILENAME = 'PredictorPortsFull.csv'

dir.create('../results/')

# Generate McLean-Pontiff bootstrapped mean returns figure ====
# read in PredictorPortsFull.csv
ret0 = fread('../data/PredictorPortsFull.csv') %>%                                           
  filter(!is.na(ret), port == 'LS')

# read in SignalDoc.csv and isolates to relevant date cols
signaldoc = fread('../data/SignalDoc.csv') %>% 
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
      , (date > sampend) & (date <= sampend %m+% months(36)) ~ 'out-of-samp'
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
    select(c(date, ret, signalname)) %>%
    pivot_wider(
      names_from = "signalname", values_from = "ret"
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

mean_insamp = ret1 %>% filter(samptype == 'in-samp') %>% 
  summarize(meanret = mean(ret)) %>% pull(meanret)

mean_oos = ret1 %>% filter(samptype == 'out-of-samp') %>% 
  summarize(meanret = mean(ret)) %>% pull(meanret)


# compile and plot
bootdat = data.frame(
  pooled_mean_ret = rboot1, samptype = 'in-samp' 
) %>% 
  rbind(
    data.frame(
      pooled_mean_ret = rboot2, samptype = 'out-of-samp' 
    )
  ) %>% 
  mutate(
    mean_ret_scaled = pooled_mean_ret/mean_insamp
  )



## main fig ----------------------------------------------------------------


bootdat %>% 
  ggplot(aes(x=pooled_mean_ret, fill=samptype)) +
  geom_histogram(alpha = 0.8
    , position = 'identity'
    , breaks = seq(0,1,0.025)
    , aes(y=..density..)
  ) +
  theme_minimal(
    base_size = 15,
  ) + 
  theme(
    text = element_text(size=30)
  ) + 
  scale_fill_manual(
    values = c('blue', 'gray'), name = "Sample Type"
  ) +
  labs(x='Pooled Mean Return (% monthly)', y='Desntiy') +
  geom_vline(xintercept = 0)

ggsave(
  "../results/MPrep.png",
  width = 12,
  height = 12
)



## alt figure --------------------------------------------------------------

  
bootdat %>% 
  filter(samptype == 'out-of-samp') %>% 
  ggplot(aes(x=mean_ret_scaled, fill=samptype)) +
  geom_histogram(alpha = 0.8
                 , position = 'identity'
                 , breaks = seq(0.2,1.3,0.025)
                 , aes(y=..density..)
  ) +
  theme_minimal(
    base_size = 15
  ) + 
  theme(
    text = element_text(size=30)
  ) +
  labs(x='Pooled Mean Return (% monthly)', y='Desntiy') +
  geom_vline(xintercept = mean_oos/mean_insamp) +
  scale_x_continuous(breaks = seq(0,1.25,0.2))



## alt figure 2 --------------------------------------------------------------


bootdat %>% 
  mutate(retplot = mean_ret_scaled*100) %>% 
  ggplot(aes(x=retplot, fill=samptype)) +
  geom_histogram(alpha = 0.8
                 , position = 'identity'
                 , breaks = seq(0,125, 5)
                 , aes(y=..density..)
  ) +
  theme_minimal(
    base_size = 15,
  ) + 
  theme(
    text = element_text(size=30)
  ) + 
  scale_fill_manual(
    values = c('blue', 'gray'), name = "Sample Type"
  ) +
  labs(x='Pooled Mean Return (bps monthly)', y='Density') +
  geom_vline(xintercept = 0) +
  scale_x_continuous(breaks = seq(0,125,25))+
  geom_vline(xintercept = mean_oos/mean_insamp*100)



ggsave(
  "../results/MPrep_scaled.pdf",
  width = 12,
  height = 8
)


# Generate R2 Replication Figure ----

## Import ====

# cz
doc =  fread("../data/SignalDoc.csv") %>% 
  rename(signalname = Acronym) %>% 
  select(-c(`Detailed Definition`, `Notes`))


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
fit_raw = target_data[
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

## Prepare signal doc ====

catname = c('raw','factor or char adjusted','nonstandard lag')

# select comparable t-stats
fit_OP = doc %>% 
  mutate(
    tstat_OP = abs(as.numeric(`T-Stat`))
  ) %>% 
  select(
    signalname, tstat_OP, `Predictability in OP`
    , `Signal Rep Quality`, `Test in OP`, `Evidence Summary`
    , SampleEndYear
  ) %>% 
  filter(
    `Signal Rep Quality` %in% c('1_good','2_fair')
    , grepl('port', `Test in OP`)
    , `Predictability in OP` != 'indirect'
  ) %>% 
  filter(!is.na(tstat_OP)) %>%   # some port sorts have only point estimates 
  mutate(
    adjusted = if_else(
      `Test in OP` %in% c('port sort', 'LS port') # these are raw long-shorts
      , catname[1]
      , catname[2]
    ) 
    , adjusted = if_else(
      grepl('nonstandard', `Evidence Summary`)
      , catname[3], adjusted
    )
    , adjusted = factor(
      adjusted
      , levels = catname
    )
  )

## Raw vs OP plot ====

# merge 
fitcomp = fit_raw %>% 
  filter(model == 'raw') %>% 
  rename(tstat_CZ = tstat) %>% 
  inner_join(
    fit_OP
    , by = 'signalname'
  ) 

# plot
tempname = 'Original Method'
fitcomp %>% 
  ggplot(aes(x=tstat_OP, y = tstat_CZ)) +
  geom_point(size=4, aes(shape =adjusted, fill = adjusted)) +
  theme_minimal(
    base_size = 15
  ) +
  theme(
    legend.position = c(.8, .25)
  ) +
  geom_abline(
    aes(slope = 1, intercept = 0)
  ) +
  scale_shape_manual(
    values = c(21, 22, 23), name = tempname
  ) +
  scale_fill_manual(
    values = c('blue', 'white', 'gray'), name = tempname
  ) +
  labs(x = 't-stat Original Paper (see legend)'
       , y = 't-stat Replicated (raw)')  +
  coord_trans(x='log10', y='log10', xlim = c(1.5, 17), ylim = c(1.0, 15)) +
  scale_x_continuous(breaks=c(2, 5, 10, 15)) +
  scale_y_continuous(breaks=c(2, 5, 10, 15))

ggsave(
  "../results/raw_op.pdf",
  width = 12,
  height = 12
)




# Generate t too big figure ====

## settings ====

# for nboot = 1000, takes 20 seconds.  
# nboot = 10,000 takes 3 minutes.  This is about 20 million predictors
library(xtable) 

nboot = 10*1000
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
  dateselect = sample(1:dim(e)[1], ndatemax, replace = T)
  
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
ndatemax = dim(e)[1]
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
tab_too_big = tab_too_big %>% 
  mutate(t_bound = paste(t_left, " - ", t_right), .before=t_left) %>%
  mutate(across(-t_bound, round, 2)) %>%
  as_tibble() %>%
  mutate_all(as.character())

tab_too_big
# mydf %>% mutate_at(vars(-vch1), funs(round(., 1)))
# mydf %>% mutate(across(starts_with("vnum"), round, 1))

tab_too_big_wide = tab_too_big %>% t()


tab_too_big_wide
stargazer(tab_too_big_wide, type = "latex", title = "Results", align=TRUE)

print(xtable(tab_too_big_wide, type = "latex"), file = "../results/tab-too-big.tex")



# Correlation Figures ------------------------------------------------------

# import data
doc =  fread("../data/SignalDoc.csv") %>%
  rename(signalname = Acronym) %>%
  select(-c(`Detailed Definition`, `Notes`))


cz_all = fread(paste0("../data/PredictorPortsFull.csv"))  %>% 
  select(signalname, port, date, ret) %>%
  left_join(
    doc %>% select(signalname, SampleStartYear, SampleEndYear)
    , by = 'signalname'
  ) %>%
  mutate(
    insamp = (year(date) >= SampleStartYear) &  (year(date) <= SampleEndYear)
  )

# use full sample for most overlap
retmat = cz_all %>% 
  filter(port == 'LS') %>% 
  select(signalname, date, ret) %>% 
  pivot_wider(names_from = signalname, values_from = ret) %>%
  select(-date) %>% 
  as.matrix()

## correlation dist ====
cormat = cor(retmat, use = 'pairwise.complete.obs')

corlong = cormat[lower.tri(cormat)]


ggplot(data.frame(cor = corlong), aes(x=cor)) +
  geom_histogram(alpha = .8, fill="blue") +
  theme_minimal(
    base_size = 15
  ) + 
  theme(
    text = element_text(size=30)
  ) + 
  labs(x = 'Pairwise Correlations'
       ,y = 'Density')

ggsave(
  "../results/correlation.pdf",
  width = 12,
  height = 12
)


## PCA ====
# some eigenvalues are negative, since we use available case correlations
# this is not a real correlation matrix
# let's just ignore this, it's not worth doing a "nearest PSD"

eigval = eigen(cormat)$values
pct_explained = cumsum(eigval/sum(eigval))*100

plotme = data.table(
  Num_PC = 1:length(pct_explained), pct_explained
)

ggplot(plotme, aes(color="blue", x=Num_PC, y = pct_explained)) + geom_line() +
  coord_cartesian(
    xlim = c(0,80)
  ) + 
  geom_line(size = 1.0) + 
  theme_minimal(
    base_size = 15
  ) + 
  theme(
    text = element_text(size=30)
  ) + 
  labs(x="Number of Pairwise Correlations", y="% Varaince Explained") +
  scale_color_manual(values="blue") +
  theme(legend.position="none")


ggsave(
  "../results/PCA.pdf",
  width = 12,
  height = 12
)





# Structural Figures ----
# Based off of qml-pub-bias checkset script

## environment ====
rm(list = ls())
source('functions.r')
tic_all = Sys.time()

## estimation settings  ====
set.check = list(
  opt_method = 'two-stage' # 'two-stage' or 'crs-only' or 'pif-grid'
  , opt_list1 = opts.crs(maxeval = 200)
  , opt_list2 = opts.qa(xtol_rel = 1e-3)  
  , model_fam = data.frame(
    mufam   = 'lognormraw' 
    , pif     = c(NA, 0.01, 0.99)
    , mua     = c(NA,    0, 2) # mua = 0 => median = 1.0
    , siga    = c(NA, 0.05, 1) # siga > 1 leads to crazy variances
    , pubfam  = 'stair' # 'stair' or 'piecelin' or 'trunc'  
    , pubpar1 = c(NA,  1/3, 2/3)
    , row.names  = c('base','lb','ub')  
  )  
)

## simulation settings ====

piflist = seq(0.1,0.9,length.out = 3)
rholist = seq(0.1,0.9,length.out = 3)
nportlist = c(200)
n_miniboot = 100


## Point Estimate ----------------------------------------------------------

tic = Sys.time()

cz_filt_tabs = import_cz(dl=F) %>% filter(tabs >= 1.96) %>% pull(tabs)
t_sample_mean = mean(cz_filt_tabs)

# Build and solve truncated normal dist. for Lambda (x) [see Wikipedia - Truncated Nomrla]
trunc_normal <- function(x) ((dnorm(2/x, mean = 0, sd = 1)) / (1-pnorm(2/x, mean =0, sd = 1))) * x - t_sample_mean
lambda <- uniroot(trunc_normal, interval = c(.5, 10), tol = .00001)

# Shrinkage formula (RAPS)
shrink <- (1/(lambda$root**2))


est.point = estimate(
  est.set = set.check
  , tabs = cz_filt_tabs
  , par.guess = random_guess(set.check,1,1243)
  , print_level = 3
)
temp = make_stats_pub(cz_filt_tabs,est.point$par) # for now assume signs don't matter
est.point$bias = temp$bias
est.point$fdrloc = temp$fdr_tabs

est.point$par %>% print()
est.point$negloglike %>% print()

## merge bias onto cz data with signalnames ====
# this isn't pretty but it'll do for now

temp = data.table(
  tabs = cz_filt_tabs, bias = est.point$bias
)

bias_dat = import_cz(dl = F) %>% 
  left_join(temp, by = 'tabs')

# read in PredictorPortsFull.csv
ret0 = fread('../data/PredictorPortsFull.csv') %>%                                           
  filter(!is.na(ret), port == 'LS')

# read in SignalDoc.csv and isolates to relevant date cols
signaldoc = fread('../data/SignalDoc.csv') %>% 
  mutate(
    signalname = Acronym
    , pubdate = as.Date(paste0(Year, '-12-31'))
    , sampend = as.Date(paste0(SampleEndYear, '-12-31'))
    , sampstart = as.Date(paste0(SampleStartYear, '-01-01'))
  ) %>% 
  arrange(signalname) %>% 
  select(signalname, pubdate, sampend, sampstart)

# merge 
ret1 = ret0 %>%                                                                 ## merging, on signalname, returns and date (sampend, sampstart, pubdate) data
  left_join(signaldoc) %>% 
  mutate(
    samptype = case_when(
      (date >= sampstart) & (date <= sampend) ~ 'insamp'
      , (date > sampend) & (date <= pubdate) ~ 'between'
      , (date > pubdate) ~ 'postpub'
      , TRUE ~ NA_character_
    )
  ) %>% 
  filter(!is.na(samptype), !is.na(ret))

# find mean returns by sample, merge and find muhat
retsum = ret1 %>% 
  group_by(signalname, samptype) %>%  
  summarize(
    rbar = mean(ret)
  ) %>% 
  pivot_wider(names_from = samptype, values_from = rbar) %>% 
  left_join(
    bias_dat %>% select(signalname, bias)
    , by = 'signalname'
  ) %>% 
  mutate(
    muhat = insamp * (1-shrink)
  )


# plot sketch
plotme = retsum %>% select(signalname, insamp, between, muhat) %>% 
  pivot_longer(cols = c(between,muhat), names_to = 'group', values_to = 'y')

ggplot(plotme, aes(color = "blue", x=insamp, y=y, group = group)) +
  geom_point(aes(color = group)) +
  theme_minimal(
    base_size = 15) +
  theme(
    text = element_text(size=30)
  ) +
  scale_color_manual(
    values = c('gray', 'blue'), name = "Sample Type"
  ) +
  labs(x = "In Sample Return", y = "Predicted Return")
  

ggsave(
  "../results/structural_figure.pdf",
  width = 12,
  height = 12
)

# ggplot(plotme, aes(color="blue", x=Num_PC, y = pct_explained)) + geom_line() +
#   coord_cartesian(
#     xlim = c(0,80)
#   ) + 
#   geom_line(size = 1.0) + 
#   theme_minimal(
#     base_size = 15
#   ) + 
#   theme(
#     text = element_text(size=30)
#   ) + 
#   labs(x="Number of Pairwise Correlations", y="% Varaince Explained") +
#   scale_color_manual(values="blue") +
#   theme(legend.position="none")


