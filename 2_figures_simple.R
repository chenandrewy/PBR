# SETUP ====
rm(list = ls())
tic = Sys.time()
dir.create('../results/')
source('functions.r')


## globals ====
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

set.seed(1057)

## graph defaults ====
MATBLUE = rgb(0,0.4470,0.7410)
MATRED = rgb(0.8500, 0.3250, 0.0980)
MATYELLOW = rgb(0.9290, 0.6940, 0.1250)

chen_theme = theme_minimal() +
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

# Import/Prepare Data ====
cz_all = fread("../data/PredictorPortsFull.csv")
signaldoc = fread('../data/SignalDoc.csv')

signaldoc = signaldoc %>% rename(signalname = Acronym)

## CZRET ====
signaldoc_dates = signaldoc %>% 
  mutate(
    pubdate = as.Date(paste0(Year, '-12-31'))
    , sampend = as.Date(paste0(SampleEndYear, '-12-31'))
    , sampstart = as.Date(paste0(SampleStartYear, '-01-01'))
  ) %>% 
  arrange(signalname) %>% 
  select(signalname, pubdate, sampend, sampstart)

czret = cz_all %>%                                         
  filter(!is.na(ret), port == 'LS') %>%                                                           
  left_join(signaldoc_dates) %>% 
  mutate(
    samptype = case_when(
      (date >= sampstart) & (date <= sampend) ~ 'in-samp'
      , (date > sampend) & (date <= sampend %m+% months(36)) ~ 'out-of-samp'
      , (date > pubdate) ~ 'post-pub'
      , TRUE ~ NA_character_
    )
  )


## CZSUM ====
czsum = cz_all %>%
  filter(port == 'LS', !is.na(ret)) %>%
  select(signalname, port, date, ret) %>% 
  filter(!is.na(ret), port == 'LS') %>%
  left_join(signaldoc %>% 
              select(signalname, SampleStartYear, SampleEndYear)
            , by = 'signalname'
  ) %>% 
  mutate(
    insamp = (year(date) >= SampleStartYear) &  (year(date) <= SampleEndYear)
  ) %>%
  group_by(signalname) %>% 
  summarize(
    tstat = mean(ret)/sd(ret)*sqrt(dplyr::n())
  ) 

## CZRETMAT ====
# used in correlations and one of the bootstraps

# use full sample for most overlap
czretmat = czret %>% 
  filter(port == 'LS') %>% 
  select(signalname, date, ret) %>% 
  pivot_wider(names_from = signalname, values_from = ret) %>%
  select(-date) %>% 
  as.matrix()



# Generate McLean-Pontiff bootstrapped mean returns figure ====
## bootstrap mean distributions ====
set.seed(6)
nboot = 50

bootfun = function(sampname){
  # make wide dataset, use NA if not correct sample
  wide_is = czret %>%
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
}

# bootstrap for each sample type
rboot1 = bootfun('in-samp')
rboot2 = bootfun('out-of-samp')

mean_insamp = czret %>% filter(samptype == 'in-samp') %>% 
  summarize(meanret = mean(ret)) %>% pull(meanret)

mean_oos = czret %>% filter(samptype == 'out-of-samp') %>% 
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

## plot --------------------------------------------------------------
bootdat %>% 
  mutate(retplot = mean_ret_scaled*100) %>% 
  ggplot(aes(x=retplot, fill=samptype)) +
  geom_histogram(alpha = 0.8
                 , position = 'identity'
                 , breaks = seq(0,125, 5)
                 , aes(y=..density..)
  ) +
  chen_theme + 
  theme(
    legend.position = c(.25, .8)
  ) + 
  scale_fill_manual(
    labels=c('Original Sample', 'Years 1-3 Post Sample'),
    values = c(MATBLUE, 'grey')
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
## performance measures ====
# t_insamp = czsum$tstat

# select comparable t-stats
fit_OP = signaldoc %>% 
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


## prepare signal doc ====
catname = c('raw','factor or char adjusted','nonstandard lag')

# select comparable t-stats
fit_OP = signaldoc %>% 
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


## plot ====
# merge 
ablines = tibble(slope = 1, 
                 intercept = 0,
                 group = factor(x = c('45 degree line'),
                                levels = c('45 degree line')))

fitcomp = czsum %>% 
  rename(tstat_CZ = tstat) %>% 
  inner_join(
    fit_OP
    , by = 'signalname'
  ) 

# plot
legname = 'Original Method'
fitcomp %>% 
  ggplot(aes(x=tstat_OP, y = tstat_CZ)) +
  geom_point(size=4, aes(shape =adjusted, fill = adjusted)) +
  chen_theme +
  geom_abline(
    aes(slope = 1, intercept = 0)
  ) +
  scale_shape_manual(
    values = c(21, 22, 23), name = legname
  ) +
  scale_fill_manual(
    values = c(MATRED, MATBLUE, MATYELLOW), name = legname
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

## bootstrap ====
# residuals
czbs = czret %>% 
  filter(samptype == 'in-samp', port == 'LS') %>% 
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
  dateselect = sample(1:dim(czbs)[1], ndatemax, replace = T)
  
  # draw returns, clustered by month
  eboot = czbs[dateselect, ]
  
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

# bootstrap
ndatemax = dim(czbs)[1]
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
pemp = ecdf(czsum$tstat)

# table settings
t_left  = tedge[1:(length(tedge)-1)]
t_right = tedge[2:length(tedge)]

# generate table
tab_too_big = data.frame(
  t_left 
  , t_right
  , N_emp = length(czsum$tstat)*(pemp(t_right) - pemp(t_left))
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
  mutate_all(as.character()) %>%
  t()

# make table with Stargazer and export to LaTeX
stargazer(tab_too_big, type = "latex", title = "Results", align=TRUE)
print(xtable(tab_too_big, type = "latex"), file = "../results/tab-too-big.tex")


# Correlation Figures ------------------------------------------------------
## correlation dist ====
cormat = cor(czretmat, use = 'pairwise.complete.obs')
corlong = cormat[lower.tri(cormat)]

ggplot(data.frame(cor = corlong), aes(x=cor)) +
  geom_histogram(alpha = .8, fill=MATBLUE) +
  chen_theme + 
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
  chen_theme +
  geom_line(size = 1.0) +
  # theme_minimal(
  #   base_size = 15
  # ) + 
  # theme(
  #   text = element_text(size=40, family = "Palatino Linotype")
  # ) + 
  labs(x="Number of Principal Components", y="% Varaince Explained") +
  scale_color_manual(values=MATBLUE) +
  theme(legend.position="none") +
  scale_y_continuous(breaks = seq(0, 100, 20))

ggsave(
  "../results/PCA.pdf",
  width = 12,
  height = 12,
  device = cairo_pdf
)


# Filling the Gap ====
# find empirical t-stats


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
F_emp = ecdf(czsum$tstat)
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
dat_fit = data.frame(
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
  
  scale_color_manual(
    values = groupdat$color, labels = groupdat$labels, name = NULL
  ) +
  
  scale_linetype_manual(
    values = groupdat$linetype, labels = groupdat$labels, name = NULL
  ) +
  
  geom_line(
    data = dat_all %>% filter(group != 'emp'), aes(color = group, linetype = group)
  ) +
  chen_theme +
  theme(
    legend.position = c(75,50)/100
    # , legend.margin = margin(t = -15, r = 20, b = 0, l = 5)
  )  +
  xlab('t-statistic') +
  ylab('Frequency') +
  coord_cartesian(
    xlim = c(0,10), ylim = c(0,0.9)
  )  

ggsave('../results/filling-the-gap.pdf', width = 12, height = 8, device = cairo_pdf)


# Shrinkage Figure ----
## prep data (Mostly Bootstrap) ---------------------------------------------------------------

t_cut = -Inf

nboot = 100
ndatemin = 15*12

boot_once = function(){
  dateselect = sample(1:dim(czretmat)[1], ndatemax, replace = T)
  
  # draw returns, clustered by month
  rboot = czretmat[dateselect, ]
  
  # summarize each predictor
  rbar = apply(rboot, 2, mean, na.rm=T)
  vol  = apply(rboot, 2, sd, na.rm=T)
  ndate  = apply(rboot, 2, function(x) sum(!is.na(x)))
  tstat = rbar/vol*sqrt(ndate)
  
  # remove if not enough observations
  tstat2 = tstat[ndate > ndatemin]
  
  dat = data.table(mean_t_trunc = mean(tstat2[tstat2 > t_cut]))
  return(dat)
}

# bootstrap!!
ndatemax = dim(czretmat)[1]
tic = Sys.time()
set.seed(1057)
bootdat = rbindlist(lapply(1:nboot, function(x) boot_once()))
toc = Sys.time()
toc - tic

# compile
dat_mean_t = bootdat %>% mutate(group = 'boot') %>% 
  rbind(
    data.table(
      mean_t_trunc = mean(czsum$tstat[czsum$tstat > t_cut]), group = 'emp'
    )
  )

## estimate ----------------------------------------------------------
tic = Sys.time()


for (i in 1:(nboot+1)){
  
  mean_t_curr = dat_mean_t$mean_t_trunc[i]

  # build and solve truncated normal dist. for Lambda (x) [see Wikipedia - Truncated Nomrla]
  trunc_normal <- function(x) 
    ((dnorm(2/x, mean = 0, sd = 1)) / (1-pnorm(2/x, mean =0, sd = 1))) * x - mean_t_curr
  lambda <- uniroot(trunc_normal, interval = c(.5, 10), tol = .00001)
  
  # shrinkage formula (RAPS)
  shrink <- (1/(lambda$root**2))
  
  # store
  dat_mean_t$shrink[i] = shrink

} # end for i in 1:nboot

dat_mean_t

# find point estimate and standard errors
estsum = list(
  point = dat_mean_t %>% filter(group == 'emp') %>% pull(shrink)
  , se  = dat_mean_t %>% filter(group == 'boot') %>% summarize(se = sd(shrink)) %>% pull(se)
)


## plot ----
plotme = czret %>% 
  group_by(samptype, signalname) %>% 
  summarise(mean = mean(ret)) %>%
  pivot_wider(names_from = samptype, values_from = mean) %>%
  mutate(
    muhat = `in-samp` * (1-shrink)
  ) %>%
  select(signalname, `out-of-samp`, `in-samp`, muhat) %>%
  pivot_longer(cols = c(`out-of-samp`, muhat), names_to = 'group', values_to = 'y')

# get regression coefficients for muhat data to draw abline
coefficients = plotme[plotme$group == "muhat", c('in-samp', 'y')]
coefficients = summary(lm(y ~ `in-samp`, data=coefficients))
muhat_slope = coefficients$coefficients['`in-samp`', 'Estimate']

ggplot(plotme %>% filter(group == "out-of-samp"), aes(x=`in-samp`, y=y)) +
  chen_theme +
  theme(
    legend.position = c(.7, .3)
  ) +
  geom_point(size = 2, color = 'gray', aes(alpha = "OOS")) +
  geom_abline(size = 1, color = MATBLUE,
              aes(alpha = "Bias adjusted \nin-samp", slope = muhat_slope, intercept = 0, color="Muhat")
  ) +
  scale_alpha_manual(name = NULL,
                       values = c(1, 1),
                       breaks = c("OOS", "Bias adjusted \nin-samp"),
                       guide = guide_legend(override.aes = list(linetype = c(0, 1),
                                                                shape = c(16, NA),
                                                                color = c('gray', MATBLUE) ) ) ) +

  labs(x="In Sample Returns", y="Out of Sample Returns") +
  xlim(-1, 3) + 
  ylim(-2, 3)

ggsave(
  "../results/shrinkage_figure.pdf",
  width = 12,
  height = 8,
  device = cairo_pdf
)


