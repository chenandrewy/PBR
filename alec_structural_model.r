# Based off of qml-pub-bias checkset script

# ENVIRONMENT ====
rm(list = ls())
source('functions.r')
tic_all = Sys.time()

library(mosaicCalc) # SOLVER
library(mosaic)
library(pracma) # access to ERF function

# run this if you get an error on import_cz below
# import_cz(dl = T) 

# SETTINGS ====

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


# Point Estimate ----------------------------------------------------------

tic = Sys.time()

cz_filt_tabs = import_cz(dl=F) %>% filter(tabs >= 1.96)


## Solver ====
t_sample_mean = mean(cz_filt_tabs$tabs)

trunc_normal = makeFun( ((dnorm(2/x, mean=0, sd=1)) / (1 - pnorm(2/x, mean = 0, sd = 1))) * x ~ x )
lambda = findZeros( trunc_normal(x) - t_sample_mean~x, xlim = range(2, 10000)) %>% pull()
shrink = (1/(lambda**2))


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



# Plot EStimation Overview -----------------------------------------------------------


p_point = hist_emp_vs_par(cz_filt_tabs,est.point$par)
p_bias = est.point$bias %>% as_tibble() %>% ggplot(aes(x=value)) + 
  geom_histogram() + xlab('bias')
p_fdrloc = est.point$fdrloc %>% as_tibble() %>% ggplot(aes(x=value)) + 
  geom_histogram() + xlab('Prob False')

grid.arrange(p_point,p_bias,p_fdrloc, nrow = 1)

toc = Sys.time()
toc - tic


# Predictor Scatter -------------------------------------------------------


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

ggplot(plotme, aes(x=insamp, y=y, group = group)) +
  geom_point(aes(color = group))
