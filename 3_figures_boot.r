# SETUP ====
tic = Sys.time()
rm(list = ls())
source('setup.r')
testrun = F # T to run quickly for testing


# Generate McLean-Pontiff bootstrapped mean returns figure ====
## bootstrap mean distributions ====
set.seed(6)
nboot = 2000

if (testrun) {nboot = 200}

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
                 , breaks = seq(0,125, 2.5)
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
nboot = 10*1000

if (testrun) {nboot = 200}

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

write.csv(tab_wide, '../results/t-too-big.csv')






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


toc = Sys.time()

tic - toc