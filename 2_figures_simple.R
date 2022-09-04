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
library(haven)

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

## Import/Prepare Data ====
cz_all = fread("../data/PredictorPortsFull.csv")
signaldoc = fread('../data/SignalDoc.csv')


signaldoc = signaldoc %>% rename(signalname = Acronym)

# czret (monthly returns)
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
      , TRUE ~ 'NA_character_'
    )
  ) %>% 
  select(signalname, date, ret, samptype)


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


# Replication Figure ------------------------------------------
## performance measures ====

# make data on hand collection
fit_OP = signaldoc %>% 
  mutate(
    tstat_OP = abs(as.numeric(`T-Stat`))
  ) %>% 
  select(
    signalname, tstat_OP, `Predictability in OP`, `Signal Rep Quality`
    , `Test in OP`, `Evidence Summary`
  ) %>% 
  filter(
    `Signal Rep Quality` %in% c('1_good','2_fair')
    , grepl('port', `Test in OP`)
    , `Predictability in OP` != 'indirect'
  ) %>% 
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
  ) %>% 
  select(signalname, starts_with('tstat'), adjusted)

fitcomp = czsum %>% 
  rename(tstat_CZ = tstat) %>% 
  filter(samptype == "in-samp") %>%
  inner_join(
    fit_OP
    , by = 'signalname'
  ) 


## plot ====


# plot
legname = 'Original Method'
breaks = c(2,4, 6, 8, 12,  16)
fitcomp %>% 
  filter(adjusted == 'raw') %>% 
  ggplot(aes(x=tstat_OP, y = tstat_CZ)) +
  geom_point(size=4, color = MATBLUE) +
  chen_theme +
  theme(
    legend.position = c(.75, .25)
  ) + 
  geom_abline(
    aes(slope = 1, intercept = 0)
  )+
  labs(x = 't-stat Original Paper'
       , y = 't-stat Replicated')  +
  coord_trans(x='log10', y='log10', xlim = c(1.5, 17), ylim = c(1.5, 15)) +
  scale_x_continuous(breaks=breaks) +
  scale_y_continuous(breaks=breaks) 

ggsave(
  "../results/raw_op.pdf",
  width = 12,
  height = 8,
  device = cairo_pdf
)



# Generate McLean-Pontiff bootstrapped mean returns figure ====
## bootstrap mean distributions ====
set.seed(6)
nboot = 500 

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

mean_insamp = mean(czsum$rbar[czsum$samptype == "in-samp"])

mean_oos = mean(czsum$rbar[czsum$samptype == "out-of-samp"])


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





# T Too Big Table ---------------------------------------------------------
## settings ====
t_left = c(seq(2,9,1), Inf)
sig_pct = pnorm(-t_left)*200
nboot = 100

## bootstrap ====
ndatemin = 240

# residuals
czresid = czret %>% 
  group_by(signalname) %>% 
  mutate(
    e = as.vector(scale(ret, center = T, scale = F))
  ) %>% 
  select(signalname, date, e) %>% 
  pivot_wider(names_from = signalname, values_from = e) %>% 
  select(-date) %>% 
  as.matrix()

boot_once = function(){
  dateselect = sample(1:dim(czresid)[1], ndatemax, replace = T)
  
  # draw returns, clustered by month
  eboot = czresid[dateselect, ]
  
  # summarize each predictor
  ebar = apply(eboot, 2, mean, na.rm=T)
  vol  = apply(eboot, 2, sd, na.rm=T)
  ndate  = apply(eboot, 2, function(x) sum(!is.na(x)))
  tstat = ebar/vol*sqrt(ndate)
  
  # remove if not enough observations
  tstat2 = tstat[ndate > ndatemin]
  
  # summarize across predictors
  F_emp = ecdf(abs(tstat2))
  
  temp = data.table(t_left, pval = 1-F_emp(t_left))
  return(temp)
}


# bootstrap!!
ndatemax = dim(czresid)[1]
## should this be: dim(czresid)[1]
tic = Sys.time()
set.seed(1057)
bootdat = rbindlist(lapply(1:nboot, function(x) boot_once()))
toc = Sys.time()
toc - tic

# summarize across bootstraps
bootsum = bootdat %>% 
  group_by(t_left) %>% 
  summarize(
    pval = mean(pval)
  )

## make table ====
Fcz = ecdf(abs(czsum$tstat[czsum$samptype == "in-samp"]))
Fyz = ecdf(abs(yzsum$tstat))
Ncz = length(czsum$tstat[czsum$samptype == "in-samp"])
Nyz = length(yzsum$tstat)

tab_long = data.frame(
  t_left
  , Count_cz = Ncz*(1-Fcz(t_left))
  , Count_yz = Nyz*(1-Fyz(t_left))
  , pct_cz = (1-Fcz(t_left))*100
  , pct_yz = (1-Fyz(t_left))*100
  , sig_pct
  , sig_boot = bootsum$pval*100
)

tab_wide = tab_long %>% t()

stargazer(tab_wide, type = "latex", title = "Results", align=TRUE)






# Correlation Figures ------------------------------------------------------
## correlation dist ====
cormat = cor(czretmat, use = 'pairwise.complete.obs')
corlong = cormat[lower.tri(cormat)]

ggplot(data.frame(cor = corlong), aes(x=cor)) +
  geom_histogram(alpha = .8, fill=MATBLUE) +
  chen_theme + 
  labs(x = 'Pairwise Corr Between Monthly Returns'
       ,y = 'Count')

ggsave(
  "../results/correlation.pdf",
  width = 12,
  height = 12,
  device = cairo_pdf,
  scale = 0.7
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
  labs(x="Number of Principal Components", y="% Varaince Explained") +
  scale_color_manual(values=MATBLUE) +
  theme(legend.position="none") +
  scale_y_continuous(breaks = seq(0, 100, 20))

ggsave(
  "../results/PCA.pdf",
  width = 12,
  height = 12,
  device = cairo_pdf,
  scale = 0.7
)


# Filling the Gap ----------------------------------------------------------------
## make plotting data ----
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
F_emp = ecdf(czsum$tstat[czsum$samptype == "in-samp"])
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

# miss
sigfit = sqrt(1+1.5^2)
dat_miss = data.frame(
  t_mid = mid2
  , prob = (pnorm(t_right2, sd = sigfit) - pnorm(t_left2, sd = sigfit))/
    (1-pnorm(2, sd = sigfit))*rescalefac*(1-F_emp(2))
  , group = 'miss'
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
dat_all = rbind(dat_emp,dat_null,dat_miss,dat_fit) %>% 
  mutate(
    group = factor(group, levels = c('emp','null','miss','fit'))
  )


## plot --------------------------------------------------------------------
groupdat = tibble(
  group = c('null', 'miss', 'fit')
  , color = c(MATRED, MATYELLOW, MATBLUE)
  , labels = c(
    TeX('$\\sigma_\\mu=0$ (Null)')
    , TeX('$\\sigma_\\mu=1.5$')    
    , TeX('$\\sigma_\\mu=3$ (Fit)')
  )
  , linetype = c('dashed','dotdash','solid')
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
    data = dat_all %>% filter(group != 'emp')
    , aes(color = group, linetype = group)
    , size = 1
  ) +
  chen_theme +
  theme(
    legend.position = c(75,50)/100
    # , legend.margin = margin(t = -15, r = 20, b = 0, l = 5)
  )  +
  xlab('t-statistic') +
  ylab('Frequency') +
  coord_cartesian(
    xlim = c(0,10), ylim = c(0,0.5)
  )  

ggsave('../results/filling-the-gap.pdf', 
       width = 12,
       height = 8,
       device = cairo_pdf
      )

# Lit Comp Figure ----
## generate theta data -----
n = 1e4

# hlz
v  = runif(n) > 0.444
se = 1500/sqrt(12)/sqrt(240)
mu = rexp(n, 1/55.5)
mu[!v] = 0
theta_hlz = mu / se

dat.hlz = data.table(
  paper = 'hlz', theta = theta_hlz
)

# cz
mu = rt(n, 4)*45
se = 22

dat.cz = data.table(
  paper = 'cz', theta = mu/se
)


#jkp
# with pub bias ( p 44)
tauc_alt = 29
tauw = 21
sigmu = (tauc_alt^2 + tauw^2)^0.5

# footnote 36
sig = (1000^2/12)^(1/2)
T = 420
se = sig/sqrt(T)
sigtheta = sigmu/se

dat.jkp = data.table(
  paper = 'jkp', theta = rnorm(n,0,sigtheta)
)

# ez
dat.ez = data.table(
  paper = 'ez', theta = rnorm(n,0,3)
)

dat = rbind(dat.hlz,dat.cz, dat.jkp, dat.ez) %>% 
  mutate(
    paper = factor(paper, levels = c('ez','hlz','cz','jkp'))
  )


## plot -----

groupdat = tibble(
  group =  c('ez','hlz','cz','jkp')
  , color = c( MATBLUE, MATRED, MATYELLOW, 'gray')
  , labels = c(
    TeX("Simple $\\theta_i \\sim$ Normal")
    , "Harvey, Liu and Zhu 2016"
    , "Chen, Zimmerman 2020"
    , "Jensen, Kelly and Pederson Forth"
  )
  , linetype = c(1,2,3,4)
)



ggplot(dat, aes(x=theta, linetype = paper, color = paper)) +
  geom_density(
    position = 'identity', alpha = 0.6, adjust = 2, size = 1
  ) +
  coord_cartesian(
    xlim = c(-10,10), ylim = c(0, 0.4)
  ) +
  scale_color_manual(
    values = groupdat$color, labels = groupdat$labels, name = NULL
  ) +
  scale_linetype_manual(
    values = groupdat$linetype, labels = groupdat$labels, name = NULL
  ) +
  chen_theme +
  theme(
    legend.position = c(.25, .75)
  ) + 
  xlab(TeX("True Predictability $\\theta_i$"))  

ggsave(
  "../results/lit-comp.pdf",
  width = 12,
  height = 8,
  device = cairo_pdf
)


## plot -----
ggplot(dat, aes(x=theta, fill = paper)) +
  geom_density(
    position = 'identity', alpha = 0.6, breaks = seq(-10,10,0.5)
  ) +
  scale_fill_manual(
    values = c(hlz = MATBLUE,
               cz  = MATRED,
               jkp = "grey",
               ez = MATYELLOW)
    ,labels = 
      c(
        "Harvey, Liu and Zhu 2016"
        , "Chen, Zimmerman 2020"
        , "Jensen, Kelly and Pederson Forth"
        , TeX("Simple $\\theta_i \\sim$ Normal")
        )
  ) +
  coord_cartesian(
    xlim = c(-10,10)
  ) + 
  chen_theme +
  theme(
    legend.position = c(.25, .75)
  ) + 
  xlab(TeX("True Predictability $\\theta_i$")) 

ggsave(
  "../results/lit-comp.pdf",
  width = 12,
  height = 8,
  device = cairo_pdf
)


# Monte Carlo: EZ -----

## simulate ez ----------------------------------------------------------------
n = 1e6
set.seed(459)

# ez
theta_ez = rnorm(n,0,3)

# simulate
dat = data.table(
  Z = rnorm(n), theta = theta_ez
) %>% 
  mutate(
    t = theta + Z, v = theta > 0, tabs = abs(t)
  ) %>% 
  mutate(
    v = factor(
      v, levels = c(TRUE,FALSE), labels = c('True Predictor', 'False Predictor')
    )
  ) %>% 
  mutate(
    tselect = t
  ) 

# fit shrinkage
datsum = dat %>% 
  mutate(
    tgroup = ntile(tselect,100)
  ) %>% 
  group_by(tgroup) %>% 
  summarize(
    tselect_left = min(tselect), tselect_right = max(tselect)
    , Etselect = mean(tselect), Etheta = mean(theta), n = dplyr::n()
    , nfalse = sum(!v)
  )


# fit "cdf"
datsum = datsum %>% 
  arrange(-Etselect) %>% 
  mutate(
    nfalse_cum = cumsum(nfalse)
    , n_cum = cumsum(n)
    , fdr_tselect_left = nfalse_cum / n_cum * 100
  ) %>% 
  arrange(Etselect)

# find hurdles
hurdle_05 = min(datsum$tselect_left[which(datsum$fdr_tselect_left < 5)])
hurdle_01 = min(datsum$tselect_left[which(datsum$fdr_tselect_left < 1)])



## plot ez (top panel) --------------------------------------------------------------------

# settings for both panels here
xlimnum = c(0,7)

set.seed(89)
iselect = sample(1:n, 1000)


ggplot(
  dat[iselect,]
  , aes(x=tselect,y=theta)) +
  geom_point(aes(group = v, color = v, shape = v)) +
  scale_shape_manual(values = c(16, 1)) +
  geom_vline(xintercept = 1.96) +
  geom_line(
    data = datsum, aes(x=Etselect, y=Etheta, linetype = "Shrinkage Correction")
  ) +
  geom_abline(aes(slope = 1, intercept = 0, linetype = "Naive Estimate (45 deg)")) +
  scale_linetype_manual(values = c(2,1)) +
  scale_color_manual(values=c(MATRED, MATBLUE)) +
  coord_cartesian(
    xlim = xlimnum, ylim = c(-2,10)
  ) +
  annotate(geom="text",
           label="Classical Hurdle",
           x=1.95, y=-1.5, vjust=-1,
           family = "Palatino Linotype",
           angle = 90
  ) +
  chen_theme +
  theme(
    legend.position = c(.25, .75)
  ) +
  xlab("t-statistic") +
  ylab(TeX("True Predictability $\\theta_i$"))

ggsave('../results/monte-carlo-ez.pdf', 
       width = 12,
       height = 8,
       device = cairo_pdf
)


## numbers for text --------------------------------------------------------

dat %>% 
  filter(tselect>2) %>% 
  summarize(
    mean(tselect)
    , mean(theta)
    , mean(theta <= 0)    
    , mean(theta) / mean(tselect)
  )

# Monte Carlo: HLZ -----

## simulate hlz ----------------------------------------------------------------

n = 1e6

# hlz baseline
v  = runif(n) > 0.444
se = 1500/sqrt(12)/sqrt(240)
mu = rexp(n, 1/55.5*se)
null = rnorm(n, 0, 0.1)
mu[!v] = null[!v]
theta_hlz = mu

# simulate
dat = data.table(
  Z = rnorm(n), theta = theta_hlz, v = v
) %>% 
  mutate(
    t = theta + Z, tabs = abs(t)
  ) %>% 
  mutate(
    v = factor(
      v, levels = c(TRUE,FALSE), labels = c('True Predictor', 'False Predictor')
    )
  ) %>%   
  mutate(
    tselect = t
  )

# fit shrinkage
datsum = dat %>% 
  mutate(
    tgroup = ntile(tselect,100)
  ) %>% 
  group_by(tgroup) %>% 
  summarize(
    tselect_left = min(tselect), tselect_right = max(tselect)
    , Etselect = mean(tselect), Etheta = mean(theta), n = dplyr::n()
    , nfalse = sum(!v)
  )


# fit "cdf"
datsum = datsum %>% 
  arrange(-Etselect) %>% 
  mutate(
    nfalse_cum = cumsum(nfalse)
    , n_cum = cumsum(n)
    , fdr_tselect_left = nfalse_cum / n_cum * 100
  ) %>% 
  arrange(Etselect)

# find hurdles
hurdle_05 = min(datsum$tselect_left[which(datsum$fdr_tselect_left < 5)])
hurdle_01 = min(datsum$tselect_left[which(datsum$fdr_tselect_left < 1)])
hurdle_bonf05 = qnorm(1-0.05/300/2)

## plot hlz (bottom pannel) ----

# settings for both panels here
xlimnum = c(-2,6)
nplot = 1000
set.seed(1)
ggplot(
  dat[sample(1:n,nplot),]
  , aes(x=tselect,y=theta)
) +
  geom_point(aes(group = v, color = v, shape = v)) +
  scale_shape_manual(values = c(16, 1)) +
  geom_vline(xintercept = 1.96, size = .75) +
  geom_vline(xintercept = hurdle_bonf05, size = .75) +
  geom_abline(aes(slope = 1, intercept = 0, linetype = "Naive Estimate (45 deg)")) +
  geom_line(
    data = datsum, aes(x=Etselect, y=Etheta, linetype = "Shrinkage Correction")
  ) +
  scale_linetype_manual(values = c(2,1)) +
  scale_color_manual(values=c(MATRED, MATBLUE)) +
  coord_cartesian(
    xlim = xlimnum, ylim = c(-2,10)
  ) +
  annotate(geom="text", 
           label="Classical Hurdle", 
           x=1.95, y=8.5, vjust=-1, 
           family = "Palatino Linotype", 
           angle = 90
  ) +
  annotate(geom="text", 
           label="HLZ Bonferroni Hurdle", 
           x=3.7, y=8.5, vjust=-1, 
           family = "Palatino Linotype", 
           angle = 90
  ) +
  chen_theme +
  theme(
    legend.position = c(.25, .75)
  ) +
  xlab("t-statistic") +
  ylab(TeX("True Predictability $\\theta_i$"))

ggsave('../results/monte-carlo-hlz.pdf', 
       width = 12,
       height = 8,
       device = cairo_pdf
)

## summarize
dat %>% 
  filter(tselect>2) %>% 
  summarize(
    mean(tselect)
    , mean(theta)
    , mean(!v)    
    , mean(theta) / mean(tselect)
  )



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
      mean_t_trunc = mean(czsum$tstat[czsum$tstat > t_cut & czsum$samptype == "in-samp"]), group = 'emp'
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

fitslope = summary(lm(y ~ 0 + `in-samp`, data=plotme %>% filter(group == 'out-of-samp')))


plotme = rename(plotme, insamp = `in-samp`)
ggplot(data = plotme %>% filter(group == "out-of-samp")) + 
  geom_point(aes(x = insamp, y = y, colour = "grey"),
           stat = "Identity", fill = NA) +
  geom_abline(aes(slope = 1-shrink, intercept = 0,linetype = "shrink est"), colour = MATBLUE, size = 1) +
  geom_abline(aes(slope = 1-shrink + 2*estsum$se, intercept = 0,linetype = "error "), size = .7, colour = MATBLUE, size = 1) +
  geom_abline(slope = 1-shrink - 2*estsum$se, intercept = 0, linetype = "dotted", color = MATBLUE, size = .7) + 
  
  # geom_smooth(method = lm, x = plotmeaes(x=insamp, y=y)) + 
  geom_smooth(method = lm, aes(x = plotme$insamp[plotme$group == "out-of-samp"], y = y, linetype = "in samp"), colour=MATRED, size = 1) +
  
  scale_linetype_manual(labels = c("Shrinkage Est.", "2 SE C.I.", "In-Samp Mean"), values = c("dotted", "solid", "solid")) +
  scale_colour_manual(labels = c("Years 1-3 Post Sample"), values = "grey") +
  guides(colour = guide_legend(override.aes = list(colour = "grey", size = 1), order = 1),
         linetype = guide_legend(title = NULL, override.aes = list(linetype = c("solid", "dotted", "solid"),
                                                                  colour = c(MATBLUE, MATBLUE, MATRED),
                                                                  size = c(1, .5, .5)),
                                                                  order = 3)) + 
  chen_theme +
  theme(
    legend.position = c(.8, .3)
  ) +
  xlim(-1, 3) + 
  ylim(-2, 3) +
  labs(y = "Out of Sample / \nBias Adjusted Return (% Monthly)", x = "In Sample Return (% Monthly)")



ggsave(
  "../results/shrinkage_figure.pdf",
  width = 12,
  height = 8,
  device = cairo_pdf
)


# YZ Figure ====


yzsum = fread('../data/yz_sum.csv')

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
F_emp = ecdf(yzsum$tabs)
dat_emp = data.frame(
  t_mid = mid
  , prob = (F_emp(t_right) - F_emp(t_left))
  , group = 'emp'
)

# null
rescalefac = mean(diff(edge))/mean(diff(edge2))
dat_null =data.frame(
  t_mid = mid2
  , prob = (pnorm(t_right2) - pnorm(t_left2))*rescalefac*2
  , group = 'null'
) 


## plot --------------------------------------------------------------------
groupdat = tibble(
  group = c('null', 'fit')
  , color = c(MATRED, MATBLUE)
  , labels = c(
    TeX('Null (Normal(0,1))')
  )
  , linetype = c('dashed','solid')
)


ggplot(dat_emp, aes(x=t_mid,y=prob)) +
  geom_bar(stat='identity',position='identity',alpha=0.6, aes(fill = group)) +
  scale_fill_manual(
    values = 'gray', labels = 'Data-Mined', name = NULL
  ) +
  scale_color_manual(
    values = groupdat$color, labels = groupdat$labels, name = NULL
  ) +
  scale_linetype_manual(
    values = groupdat$linetype, labels = groupdat$labels, name = NULL
  ) +
  geom_line(
    data = dat_null, aes(color = group, linetype = group)
  ) +
  chen_theme +
  theme(
    legend.position = c(75,75)/100
    # , legend.margin = margin(t = -15, r = 20, b = 0, l = 5)
  )  +
  xlab('t-statistic') +
  ylab('Frequency') +
  coord_cartesian(
    xlim = c(0,7), ylim = c(0,0.4)
  )  

ggsave(
  "../results/yan-zheng.pdf",
  width = 12,
  height = 8,
  device = cairo_pdf
)



