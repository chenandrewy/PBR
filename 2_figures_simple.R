# 2022 07 Alec Erb on behalf of Andrew Chen: Generate Simple Figures for PBR

# SETUP ====

rm(list = ls())
source('functions.r')
tic = Sys.time()
dir.create('../results/')

## Globals ====
library(tidyverse)
library(data.table)
library(googledrive)
library(gridExtra)
library(nlme)
library(stargazer)
library(lubridate)
library(latex2exp)
library(Cairo)
library(xtable) 


# graph defaults
MATBLUE = rgb(0,0.4470,0.7410)
MATRED = rgb(0.8500, 0.3250, 0.0980)
MATYELLOW = rgb(0.9290, 0.6940, 0.1250)

groupdat = list(
  group = c('null', 'fit')
  , color = c(MATRED, MATBLUE, MATYELLOW)
  , labels = c(
    TeX('Null ($\\sigma_\\mu=0$)')
    , TeX('Fit ($\\sigma_\\mu=3$)')
  )
  , linetype = c('dashed','solid')
)

theme_set(
  theme_minimal() +
    theme(
      text = element_text(family = "Palatino Linotype.ttf")
    )
)




# Import Data ====

cz_all = fread("../data/PredictorPortsFull.csv")

signaldoc = fread('../data/SignalDoc.csv')



# Generate McLean-Pontiff bootstrapped mean returns figure ====
# read in PredictorPortsFull.csv

# read in SignalDoc.csv and isolates to relevant date cols
signaldoc_mp = signaldoc %>% 
  mutate(
    signalname = Acronym
    , pubdate = as.Date(paste0(Year, '-12-31'))
    , sampend = as.Date(paste0(SampleEndYear, '-12-31'))
    , sampstart = as.Date(paste0(SampleStartYear, '-01-01'))
  ) %>% 
  arrange(signalname) %>% 
  select(signalname, pubdate, sampend, sampstart)


# merge and find OOS returns
ret1 = cz_all %>%                                         
  filter(!is.na(ret), port == 'LS') %>%                                                           
  left_join(signaldoc_mp) %>% 
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
    text = element_text(size=40, family = "Palatino Linotype"),
    legend.title=element_text(size=40)
  ) + 
  scale_fill_manual(
    values = c(groupdat$color[2], 'grey'), name = "Sample Type"
  ) +
  labs(x='Pooled Mean Return (bps monthly)', y='Density') +
  geom_vline(xintercept = 0) +
  scale_x_continuous(breaks = seq(0,125,25))+
  geom_vline(xintercept = mean_oos/mean_insamp*100)

ggsave(
  "../results/MPrep_scaled.pdf",
  width = 12,
  height = 8,
  device = cairo_pdf
)


# Generate R2 Replication Figure ----

## Import ====
signaldoc_r2 =  signaldoc %>% 
  rename(signalname = Acronym) %>% 
  select(-c(`Detailed Definition`, `Notes`))


cz_r2 = cz_all %>% 
  mutate(vint = 2022) %>% 
  rbind(
    cz_all %>% 
      mutate(vint = 2021)
  ) %>% 
  select(vint, signalname, port, date, ret) %>% 
  left_join(
    signaldoc_r2 %>% select(signalname, SampleStartYear, SampleEndYear)
    , by = 'signalname'
  ) %>% 
  mutate(
    insamp = (year(date) >= SampleStartYear) &  (year(date) <= SampleEndYear)
  ) 

## Performance Measures ====
target_data = cz_r2 %>%
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
fit_OP = signaldoc_r2 %>% 
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
fit_OP = signaldoc_r2 %>% 
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
    legend.position = c(.8, .25),
    text = element_text(size=30, family = "Palatino Linotype")
  ) +
  geom_abline(
    aes(slope = 1, intercept = 0)
  ) +
  scale_shape_manual(
    values = c(21, 22, 23), name = tempname
  ) +
  scale_fill_manual(
    values = c(groupdat$color), name = tempname
  ) +
  labs(x = 't-stat Original Paper (see legend)'
       , y = 't-stat Replicated (raw)')  +
  coord_trans(x='log10', y='log10', xlim = c(1.5, 17), ylim = c(1.0, 15)) +
  scale_x_continuous(breaks=c(2, 5, 10, 15)) +
  scale_y_continuous(breaks=c(2, 5, 10, 15))

ggsave(
  "../results/raw_op.pdf",
  width = 12,
  height = 8,
  device = cairo_pdf
)




# Generate t too big figure ====

## settings ====
# for nboot = 1000, takes 20 seconds.  
# nboot = 10,000 takes 3 minutes.  This is about 20 million predictors
nboot = 10*1000
ndatemin = 240
tedge = c(seq(0,9, 1), 20)


## Import ====
# All signals
signaldoc_ttb = signaldoc %>%
  rename(signalname = Acronym) %>%
  select(-c(`Detailed Definition`, `Notes`))


# Chen-Zimmerman dataset
cz_ttb = cz_all  %>% 
  select(signalname, port, date, ret) %>%
  left_join(
    signaldoc_ttb %>% select(signalname, SampleStartYear, SampleEndYear)
    , by = 'signalname'
  ) %>%
  mutate(
    insamp = (year(date) >= SampleStartYear) &  (year(date) <= SampleEndYear)
  )

# Performance Measures 
target_data = cz_ttb %>% 
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


## bootstrap ====
# residuals
cz_ttb = cz_ttb %>% 
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
  dateselect = sample(1:dim(cz_ttb)[1], ndatemax, replace = T)
  
  # draw returns, clustered by month
  eboot = cz_ttb[dateselect, ]
  
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
ndatemax = dim(cz_ttb)[1]
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
  ) %>% 
  mutate(t_bound = paste(t_left, " - ", t_right), .before=t_left) %>%
  mutate(across(-t_bound, round, 6)) %>%
  as_tibble() %>%
  mutate_all(as.character())


tab_too_big_wide = tab_too_big %>% t()

# make table with Stargazer and export to LaTeX
stargazer(tab_too_big_wide, type = "latex", title = "Results", align=TRUE)
print(xtable(tab_too_big_wide, type = "latex"), file = "../results/tab-too-big.tex")




# Correlation Figures ------------------------------------------------------

# import data
signaldoc_corr =  signaldoc %>%
  rename(signalname = Acronym) %>%
  select(-c(`Detailed Definition`, `Notes`))

cz_corr = cz_all  %>% 
  select(signalname, port, date, ret) %>%
  left_join(
    signaldoc_corr %>% select(signalname, SampleStartYear, SampleEndYear)
    , by = 'signalname'
  ) %>%
  mutate(
    insamp = (year(date) >= SampleStartYear) &  (year(date) <= SampleEndYear)
  )

# use full sample for most overlap
retmat = cz_corr %>% 
  filter(port == 'LS') %>% 
  select(signalname, date, ret) %>% 
  pivot_wider(names_from = signalname, values_from = ret) %>%
  select(-date) %>% 
  as.matrix()


## correlation dist ====
cormat = cor(retmat, use = 'pairwise.complete.obs')
corlong = cormat[lower.tri(cormat)]

ggplot(data.frame(cor = corlong), aes(x=cor)) +
  geom_histogram(alpha = .8, fill=groupdat$color[2]) +
  theme_minimal(
    base_size = 15
  ) + 
  theme(
    text = element_text(size=40, family = "Palatino Linotype")
  ) + 
  labs(x = 'Pairwise Correlations'
       ,y = 'Count')

ggsave(
  "../results/correlation.pdf",
  width = 12,
  height = 12,
  device = cairo_pdf
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
    text = element_text(size=40, family = "Palatino Linotype")
  ) + 
  labs(x="Number of Principal Components", y="% Varaince Explained") +
  scale_color_manual(values=groupdat$color[2]) +
  theme(legend.position="none")

ggsave(
  "../results/PCA.pdf",
  width = 12,
  height = 12,
  device = cairo_pdf
)


# Filling the Gap ====


## environment -------------------------------------------------------------
# read returns
cz_gap = cz_all %>%                                           
  filter(!is.na(ret), port == 'LS')

# read in SignalDoc.csv and isolates to relevant date cols
signaldoc_gap = signaldoc %>% 
  mutate(
    signalname = Acronym
    , pubdate = as.Date(paste0(Year, '-12-31'))
    , sampend = as.Date(paste0(SampleEndYear, '-12-31'))
    , sampstart = as.Date(paste0(SampleStartYear, '-01-01'))
  ) %>% 
  arrange(signalname) %>% 
  select(signalname, pubdate, sampend, sampstart)

# merge 
ret1 = cz_gap %>%                                                                
  left_join(signaldoc_gap) %>% 
  mutate(
    samptype = case_when(
      (date >= sampstart) & (date <= sampend) ~ 'in-samp'
      , (date > sampend) & (date <= pubdate) ~ 'out-of-samp'
      , (date > pubdate) ~ 'post-pub'
      , TRUE ~ NA_character_
    )
  ) %>% 
  filter(!is.na(samptype))

# find empirical t-stats
t_emp = ret1 %>% filter(samptype == 'in-samp', port == 'LS') %>% 
  group_by(signalname) %>% 
  summarize(
    tstat = mean(ret)/sd(ret)*sqrt(dplyr::n())
  )  %>% 
  pull(tstat)


## make plotting data ----------------------------------------------------------------
# set up 
edge = seq(0,20,0.5)
t_left = edge[1:(length(edge)-1)]
t_right = edge[2:length(edge)]
mid = t_left + diff(edge)/2

edge2  = seq(0,20,0.1)
t_left2 = edge2[1:(length(edge2)-1)]
t_right2 = edge2[2:length(edge2)]
mid2 = t_left2 + diff(edge2)/2


# empirical
F_emp = ecdf(t_emp)
dat_emp = data.frame(
  t_mid = mid
  , prob = (F_emp(t_right) - F_emp(t_left))
  , group = 'emp'
)

# null
rescalefac = mean(diff(edge))/mean(diff(edge2))
dat_null =data.frame(
  t_mid = mid2
  , prob = (pnorm(t_right2) - pnorm(t_left2))/(1-pnorm(2))*rescalefac*(1-F_emp(2))
  , group = 'null'
) %>% 
  mutate(prob = if_else(t_mid > 2, prob, 0))

# fit
sigfit = sqrt(1+3^2)
dat_fit =data.frame(
  t_mid = mid2
  , prob = (pnorm(t_right2, sd = sigfit) - pnorm(t_left2, sd = sigfit))/
    (1-pnorm(2, sd = sigfit))*rescalefac*(1-F_emp(2))
  , group = 'fit'
) %>% 
  mutate(prob = if_else(t_mid > 2, prob, 0))

# merge and factor group
dat_all = rbind(dat_emp,dat_null,dat_fit) %>% 
  mutate(
    group = factor(group, levels = c('emp','null','fit'))
  )


## plot --------------------------------------------------------------------
groupdat = tibble(
  group = c('null', 'fit')
  , color = c(MATRED, MATBLUE)
  , labels = c(
    TeX('Null ($\\sigma_\\mu=0$)')
    , TeX('Fit ($\\sigma_\\mu=3$)')
  )
  , linetype = c('dashed','solid')
)


ggplot(dat_all %>%  filter(group == 'emp'), aes(x=t_mid,y=prob)) +
  geom_bar(stat='identity',position='identity',alpha=0.6, aes(fill = group)) +
  scale_fill_manual(
    values = 'gray', labels = 'Published', name = NULL
  ) +    
  geom_line(
    data = dat_all %>% filter(group != 'emp'), aes(color = group, linetype = group)
  ) +
  xlab('t-statistic') +
  ylab('Frequency') +
  scale_color_manual(
    values = groupdat$color, labels = groupdat$labels, name = NULL
  ) + 
  scale_linetype_manual(
    values = groupdat$linetype, labels = groupdat$labels, name = NULL
  ) +
  theme(
    legend.position = c(75,75)/100
    , legend.margin = margin(t = -15, r = 20, b = 0, l = 5)
    , text = element_text(family = "Palatino Linotype")
    
  )  +
  coord_cartesian(
    xlim = c(0,10), ylim = c(0,0.25)
  )  

ggsave('../results/filling-the-gap.pdf', width = 12, height = 8, device = cairo_pdf)




# Shrinkage Figure ----
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


## point estimate ----------------------------------------------------------
tic = Sys.time()
t_sample_mean = mean(t_emp[t_emp > 1.96])

# build and solve truncated normal dist. for Lambda (x) [see Wikipedia - Truncated Nomrla]
trunc_normal <- function(x) ((dnorm(2/x, mean = 0, sd = 1)) / (1-pnorm(2/x, mean =0, sd = 1))) * x - t_sample_mean
lambda <- uniroot(trunc_normal, interval = c(.5, 10), tol = .00001)

# shrinkage formula (RAPS)
shrink <- (1/(lambda$root**2))

est.point = estimate(
  est.set = set.check
  , tabs = t_emp
  , par.guess = random_guess(set.check,1,1243)
  , print_level = 3
)

temp = make_stats_pub(t_emp,est.point$par) # for now assume signs don't matter
est.point$bias = temp$bias
est.point$fdrloc = temp$fdr_tabs

est.point$par %>% print()

## merge bias onto cz data with signalnames ====
# this isn't pretty but it'll do for now
temp = data.table(
  tabs = t_emp, bias = est.point$bias
)

bias_dat = fread("output/pubcross.csv") %>% 
  left_join(temp, by = 'tabs')

# read in PredictorPortsFull.csv
cz_bias = cz_all %>%                                           
  filter(!is.na(ret), port == 'LS')

# read in SignalDoc.csv and isolates to relevant date cols
signaldoc_bias = signaldoc %>% 
  mutate(
    signalname = Acronym
    , pubdate = as.Date(paste0(Year, '-12-31'))
    , sampend = as.Date(paste0(SampleEndYear, '-12-31'))
    , sampstart = as.Date(paste0(SampleStartYear, '-01-01'))
  ) %>% 
  arrange(signalname) %>% 
  select(signalname, pubdate, sampend, sampstart)

# merge 
ret1 = cz_bias %>%                                                                 ## merging, on signalname, returns and date (sampend, sampstart, pubdate) data
  left_join(signaldoc_bias) %>% 
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
  height = 12,
  device = cairo_pdf
)


