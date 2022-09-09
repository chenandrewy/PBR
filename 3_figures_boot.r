# SETUP ====
rm(list = ls())
tic = Sys.time()
dir.create('../results/')
source('functions.r')


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

# ret from disc
cz_all = fread("../data/PredictorPortsFull.csv")
yzsum = fread('../data/yz_sum.csv')
signaldoc = fread('../data/SignalDoc.csv') %>% 
  rename(signalname = Acronym) %>% 
  mutate(
    pubdate = as.Date(paste0(Year, '-12-31'))
    , sampend = as.Date(paste0(SampleEndYear, '-12-31'))
    , sampstart = as.Date(paste0(SampleStartYear, '-01-01'))
  )


## CZRET ===
czret = cz_all %>%         
  filter(port == 'LS', !is.na(ret)) %>% 
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
 

## CZSUM ===
czsum = czret %>%
  group_by(signalname, samptype) %>%
  summarize(
    rbar = mean(ret)
    , vol = sd(ret)
    , tstat = mean(ret)/sd(ret)*sqrt(dplyr::n())
  )

# Generate McLean-Pontiff bootstrapped mean returns figure ====
## bootstrap mean distributions ====
set.seed(6)
nboot = 2000

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




