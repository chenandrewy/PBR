# Setup -------------------------------------------------------------------
rm(list = ls())
source('setup.r')
 


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
    , `Test in OP` %in% c('port sort', 'LS port') # these are raw long-shorts
  ) %>% 
  select(signalname, tstat_OP) %>% 
  filter(!is.na(tstat_OP))

fitcomp = czsum %>% 
  rename(tstat_CZ = tstat) %>% 
  filter(samptype == "in-samp") %>%
  inner_join(
    fit_OP
    , by = 'signalname'
  ) 


## plot ====


# plot
breaks = c(2,4, 6, 8, 12,  16)
fitcomp %>% 
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









# Correlation Figures ------------------------------------------------------
## correlation dist ====

czret %>% 
  group_by(signalname) %>% 
  summarize(
    
  )
  

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
    TeX('$\\sigma_\\theta$ = 0 (Null)')
    , TeX('$\\sigma_\\theta$ = $1.5$')    
    , TeX('$\\hat{\\sigma}_\\theta$ = $3$ (Fit)')
  )
  , linetype = c('dotted','longdash','solid')
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
    xlim = c(0,15), ylim = c(0,0.2)
  )  

ggsave('../results/filling-the-gap.pdf', 
       width = 12,
       height = 8,
       device = cairo_pdf
      )

# Lit Comp Figure ----
## generate theta data -----
n = 1e5

# ez (really, CZ 20 appendix)
dat.ez = data.table(
  paper = 'ez', theta = rnorm(n,0,3)
)


# hlz
v  = runif(n) > 0.444
se = 1500/sqrt(12)/sqrt(240)
mu = rexp(n, 1/55.5)
mu[!v] = 0
theta_hlz = mu / se

dat.hlz = data.table(
  paper = 'hlz', theta = theta_hlz
)

# cz: page 18 Table 3 All
mu = rt(n, 4)*45
se = exp(rnorm(n, -1.67, 0.52))*100

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
nmonth = 420
se = sig/sqrt(nmonth)
sigtheta = sigmu/se

dat.jkp = data.table(
  paper = 'jkp', theta = rnorm(n,0,sigtheta)
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
    TeX("Simple Normal")
    , "Harvey-Liu-Zhu"
    , "Chen-Zimmerman"
    , "Jensen-Kelly-Pederson"
  )
  , linetype = c(1,2,3,4)
)



ggplot(dat, aes(x=theta, linetype = paper, color = paper)) +
  geom_density(
    position = 'identity', alpha = 0.6, adjust = 1.2, size = 1.2
  ) +
  scale_color_manual(
    values = groupdat$color, labels = groupdat$labels, name = NULL
  ) +
  scale_linetype_manual(
    values = groupdat$linetype, labels = groupdat$labels, name = NULL
  ) +
  chen_theme +
  theme(legend.position = c(.2, .75)) + 
  coord_cartesian(xlim = c(-8,8), ylim = c(0, 0.3)) +  
  scale_x_continuous(breaks = seq(-20,20,2)) +
  xlab(TeX("Corrected t-statistic $\\theta_i$"))  +
  ylab('Density')

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
    t = theta + Z, v = theta > 0, tabs = abs(t), pub = t > 2
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
    , nfalse = sum(v == 'False Predictor')
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


# find pub stats
pubstat = dat %>% 
  filter(tselect>2) %>% 
  summarize(
    Et = mean(tselect)
    , Etheta = mean(theta)
    , FDR = mean(theta <= 0)    
    , shrink = 1-mean(theta) / mean(tselect)
  )

pubplot = tibble(
  tselect = seq(2,10,0.1), theta = pubstat$Et, group = 'naive' 
)  %>% rbind(
  tibble(
    tselect = seq(2,10,0.1), theta = pubstat$Etheta, group = 'shrinkage' 
  )   
)




## plot pub only bar --------------------------------------------------------------------

plotme = dat %>% 
  filter( pub ) %>% 
  select(theta, t) %>% 
  pivot_longer(cols  = c(theta,t)) 

plotme = tibble(
  value = czsum$tstat[czsum$samptype == 'in-samp' & czsum$tstat > 0]
                , name = 't') %>% rbind(
  tibble(value = dat$theta[dat$t>2], name = 'theta')
)

plotmeans = plotme %>% group_by(name) %>% summarize(mean = mean(value)) %>% 
  pivot_wider(names_from = name, values_from = mean)

color_theta = MATBLUE
color_emp = 'dimgrey'
textsize = 7

ggplot(plotme, aes(x = value, group = name))   +
  geom_histogram(
    aes(y = ..density.., fill = name), position = 'identity', alpha = 0.5
    , breaks = c(seq(-1,16,0.5)) , color = 'white'
  ) +
  scale_fill_manual(
    values = c('dimgrey', MATBLUE)
    , labels = c(TeX('Observed $(t_i|pub_i)$')
                 ,TeX('Corrected ($\\hat{\\theta}_i|pub_i$)'))
    , name = NULL
  ) +  
  # shrinkage
  geom_vline(xintercept = plotmeans$t, color = 'dimgrey', size = 1
             , linetype = 'dashed') +
  geom_vline(xintercept = plotmeans$theta, color = color_theta, size = 1.1
             , linetype = 'dashed') +
  annotate(geom="text"
           , label='<-- Shrinkage',
           x=51/10, y=0.32, vjust=-1,
           family = "Palatino Linotype",  angle = 0, size = textsize, color = 'black' ) +
  # fdr
  geom_vline(xintercept = 0.25/2, color = MATRED, size = 1) +
  annotate(geom="text"
           , label="<- False",
           x=-8/10, y=0.28, vjust=-1,
           family = "Palatino Linotype",  size = textsize, color = 'black'
  ) +
  annotate(geom="text"
           , label="True ->",
           x=10/10, y=0.28, vjust=-1,
           family = "Palatino Linotype",  size = textsize, color = 'black'
  ) +  
  chen_theme  +
  theme(legend.position = c(75,60)/100) +
  coord_cartesian(xlim = c(-1, 15), ylim = c(0, 0.35)) +
  scale_x_continuous(breaks = seq(-10,20,2)) +  
  xlab('Published t-statistic') + ylab('Density')


ggsave('../results/monte-carlo-ez-bar.pdf', 
       width = 14, height = 8, device = cairo_pdf, scale = 0.8)


## plot ez scatter  --------------------------------------------------------------------
# (not used)

set.seed(89)
iselect = sample(1:n, 1000)


ggplot(
  dat[iselect,]
  , aes(x=tselect,y=theta)) +
  geom_point(aes(group = v, color = v, shape = v), size = 2.5) +
  scale_shape_manual(values = c(16, 1)) +
  # geom_line(data = pubplot, aes(group = group)) +  
  geom_line(
    data = datsum, aes(x=Etselect, y=Etheta, linetype = "Shrinkage Correction")
  ) +
  geom_abline(aes(slope = 1, intercept = 0, linetype = "Naive Estimate (45 deg)")) +
  scale_linetype_manual(values = c(1,2)) +
  scale_color_manual(values=c(MATBLUE, MATRED)) +
  coord_cartesian(
    xlim = c(0,7), ylim = c(-2,8)
  ) +
  geom_vline(xintercept = 1.96) +  
  annotate(geom="text",
           label="Classical Hurdle",
           x=1.95, y=6.5, vjust=-1,
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'black'
  ) +
  chen_theme +
  theme(
    legend.position = c(.8, .25)
  ) +
  xlab(TeX("t-statistic $\\t_i$")) +
  ylab(TeX("Standardized Expected Return $\\theta_i$"))

ggsave('../results/monte-carlo-ez.pdf', 
       width = 12,
       height = 8,
       device = cairo_pdf
)





## plot pub only line ---------------------------------------------------------
# (not used)

# set up 
edge = seq(-1,20,0.5)
t_left = edge[1:(length(edge)-1)]
t_right = edge[2:length(edge)]
mid = t_left + diff(edge)/2

edge2  = seq(-1,20,0.5)
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

# theta
F_theta = ecdf(dat$theta[dat$pub])
rescalefac = mean(diff(edge))/mean(diff(edge2))
dat_theta =data.frame(
  t_mid = mid2
  , prob = (F_theta(t_right2) - F_theta(t_left2))*rescalefac
  , group = 'null'
) 


labelmod = TeX('$E(\\theta_i|pub_i)$')

ggplot(dat_emp, aes(x=t_mid,y=prob)) +
  geom_bar(stat='identity',position='identity',alpha=0.6, aes(fill = group)) +
  scale_fill_manual(
    values = 'gray', labels = 'Published', name = NULL
  ) +
  geom_line(
    data = dat_theta
    , aes(color = group, linetype = group)
    , size = 1.5
  ) +
  scale_color_manual(
    values = MATBLUE, labels = labelmod, name = NULL
  ) +
  scale_linetype_manual(
    values = 'solid', labels = labelmod, name = NULL
  ) +  
  chen_theme +
  theme(
    legend.position = c(75,50)/100
    # , legend.margin = margin(t = -15, r = 20, b = 0, l = 5)
  )  +
  geom_vline(
    xintercept = mean(czsum$tstat[czsum$samptype == "in-samp"]), color = 'darkgrey'
    , size = 1
    ) +
  geom_vline(xintercept = mean(dat$theta[dat$pub]), color = MATBLUE
             , size = 1) +
  geom_vline(xintercept = 0, color = MATRED
             , size = 1) +  
  xlab('t-statistic') +
  ylab('Frequency') +
  coord_cartesian(
    xlim = c(-2,15), ylim = c(0,0.2)
  )  +
  scale_x_continuous(breaks = seq(-10,20,2))



## numbers for text --------------------------------------------------------

dat %>% 
  filter(tselect>2) %>% 
  summarize(
    mean(tselect)
    , mean(theta)
    , mean(theta <= 0)    
    , mean(theta) / mean(tselect)
  )

# for fdr approx
1-pnorm(2, 0, sqrt(10))

0.025/(1-pnorm(2, 0, sqrt(10))) * 0.5

# Monte Carlo: HLZ -----

## simulate hlz ----------------------------------------------------------------

n = 1e6

# hlz baseline
v  = runif(n) > 0.444
se = 1500/sqrt(12)/sqrt(240)
mu = rexp(n, 1/55.5); mu[!v] = 0
theta = mu / se
theta_scatter = theta; theta_scatter[!v] = rnorm(sum(!v), 0, 0.1)
pubnoise = runif(n)
 

# simulate
dat = data.table(
  Z = rnorm(n), theta, v, pubnoise, theta_scatter
) %>% 
  mutate(
    t = theta + Z
    , tabs = abs(t)
    , pub = case_when(
      tabs <= 1.96 ~ FALSE
      , (1.96 < tabs) & (tabs <= 2.57) ~ pubnoise < 0.5
      , 2.57 < tabs  ~ TRUE
    )
  ) %>% 
  mutate(
    v = factor(
      v, levels = c(TRUE,FALSE), labels = c('True Predictor', 'False Predictor')
    )
  ) %>%   
  mutate(
    tselect = tabs
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
    , nfalse = sum(v == 'False Predictor') 
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


# find pub stats
pubstat = dat %>% 
  filter(tselect>2) %>% 
  summarize(
    Et = mean(tselect)
    , Etheta = mean(theta)
    , FDR = mean(theta <= 0)    
    , shrink = 1-mean(theta) / mean(tselect)
  )

pubplot = tibble(
  tselect = seq(2,10,0.1), theta = pubstat$Et, group = 'naive' 
)  %>% rbind(
  tibble(
    tselect = seq(2,10,0.1), theta = pubstat$Etheta, group = 'shrinkage' 
  )   
)



## plot pub only bar --------------------------------------------------------------------

plotme = tibble(
  value = czsum$tstat[czsum$samptype == 'in-samp' & czsum$tstat > 0]
  , name = 't') %>% rbind(
    tibble(value = dat$theta[dat$pub], name = 'theta')
  )

plotmeans = plotme %>% group_by(name) %>% summarize(mean = mean(value)) %>% 
  pivot_wider(names_from = name, values_from = mean)

color_theta = MATBLUE
color_emp = 'dimgrey'
textsize = 7

ggplot(plotme, aes(x = value, group = name))   +
  geom_histogram(
    aes(y = ..density.., fill = name), position = 'identity', alpha = 0.5
    , breaks = c(seq(-1,16,0.5)) , color = 'white'
  ) +
  scale_fill_manual(
    values = c('dimgrey', MATBLUE)
    , labels = c(TeX('Observed $(|t_i||pub_i)$')
                 ,TeX('Corrected ($\\hat{\\theta}_i|pub_i$)'))
    , name = NULL
  ) +  
  # shrinkage
  geom_vline(xintercept = plotmeans$t, color = 'dimgrey', size = 1
             , linetype = 'dashed') +
  geom_vline(xintercept = plotmeans$theta, color = color_theta, size = 1.1
             , linetype = 'dashed') +
  annotate(geom="text"
           , label='<-- Shrinkage',
           x=51.2/10, y=0.32, vjust=-1,
           family = "Palatino Linotype",  angle = 0, size = textsize, color = 'black' ) +
  # fdr
  geom_vline(xintercept = 0.25/2, color = MATRED, size = 1) +
  annotate(geom="text"
           , label="<- False",
           x=-8/10, y=0.28, vjust=-1,
           family = "Palatino Linotype",  size = textsize, color = 'black'
  ) +
  annotate(geom="text"
           , label="True ->",
           x=10/10, y=0.28, vjust=-1,
           family = "Palatino Linotype",  size = textsize, color = 'black'
  ) +  
  chen_theme  +
  theme(legend.position = c(75,60)/100) +
  coord_cartesian(xlim = c(-1, 15), ylim = c(0, 0.35)) +
  scale_x_continuous(breaks = seq(-10,20,2)) +  
  xlab('Published Absolute t-statistic') + ylab('Density')


ggsave('../results/monte-carlo-hlz-bar.pdf', 
       width = 14, height = 8, device = cairo_pdf, scale = 0.8)


## numbers for text --------------------------------------------------------

print('shrinkage')
1-plotmeans$theta/plotmeans$t


print('FDR')
plotme %>% filter(name == 'theta') %>% summarize(mean(value<=0)*100) 



## plot scatter  ----

# holm algo

ntotal = 300/mean(dat$pub)

set.seed(430)

holm_05 = dat[sample(1:n,ntotal), ] %>% select(tabs) %>% 
  arrange(desc(tabs)) %>% 
  mutate(
    pval = 2*pnorm(-tabs)
    , k = row_number()
    , signif = pval < 0.05/(ntotal + 1 - k)
  ) %>% 
  filter(signif == F) %>% 
  filter(row_number() == 1) %>% 
  pull(tabs)
  
  


# settings for both panels here
nplot = 1000
set.seed(5)

texty = 7
textsize = 7

linesize = 1.1

ggplot(dat[sample(1:n,nplot),], aes(x=tselect,y=theta_scatter)) +
  geom_point(aes(group = v, color = v, shape = v), size = 2.5) +
  scale_shape_manual(values = c(16, 1)) +
  scale_color_manual(values=c(MATBLUE, MATRED)) +
  # HURDLES
  geom_vline(xintercept = 1.96, size = linesize) +
  annotate(geom="text", label="Classical Hurdle", 
           x=1.95, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'black'
  ) +
  geom_vline(xintercept = hurdle_05, size = linesize, color = 'dimgrey', linetype = 'longdash') +    
  annotate(geom="text", 
           label="FDR = 5%", 
           x=24/10, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'dimgrey'
  ) +  
  geom_vline(xintercept = hurdle_01, size = linesize, color = MATRED, linetype = 'dotdash') +  
  annotate(geom="text", 
           label="FDR = 1%", 
           x=3, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = MATRED
  ) +
  geom_vline(xintercept = holm_05, size = linesize, color = 'darkorchid', linetype = 'dotted') +
  annotate(geom="text", 
           label=TeX("Holm 5\\%"), 
           x=holm_05, y=texty, vjust=-1, 
           family = "Palatino Linotype", angle = 90, size = textsize, color = 'darkorchid'
  ) +  
  # SHRINKAGE
  # geom_abline(aes(slope = 1, intercept = 0, linetype = "Naive (45 deg)")) +
  # geom_line(
  #   data = datsum, aes(x=Etselect, y=Etheta, linetype = "Shrinkage")
  # ) +
  # scale_linetype_manual(values = c(1,2)) +
  coord_cartesian(xlim = c(-0.1,10), ylim = c(-0.5,10)) +
  scale_x_continuous(breaks = seq(-10,20,2)) +
  scale_y_continuous(breaks = seq(-10,20,2)) +  
  chen_theme +
  theme(
    legend.position = c(.80, .15)
    , panel.grid.major = element_blank()
    , panel.grid.minor = element_blank()
  ) +
  xlab(TeX("Absolute t-statistic $|\\t_i|$")) +
  ylab(TeX("Corrected t-statistic $\\theta_i$"))

ggsave('../results/hlz-scatter.pdf', 
       width = 12,
       height = 8,
       device = cairo_pdf
)




## numbers for text --------------------------------------------------------




dat %>% 
  filter(tselect > 1.96) %>% 
  summarize(
    mean(tselect < hurdle_01)
  )


dat %>% 
  filter(tselect > 2,  tselect < hurdle_01) %>% 
  summarize(
    mean(t)
    , mean(theta)
    , mean(v == 'False Predictor')    
    , mean(theta) / mean(t)
  )








# Liquidity Screen bars -------------------------------------------------------

rbarbarbase = czsum %>% filter(samptype == 'in-samp') %>% 
  ungroup() %>% 
  summarize(rbarbar = mean(rbar)) %>% 
  pull(rbarbar)


group_signal_sum = cz_alt %>% 
  filter(samptype == 'in-samp') %>% 
  group_by(group, signalname) %>% 
  summarize(
    rbar = mean(ret)/rbarbarbase*100
  ) %>% 
  mutate(
    group = factor(group)
  ) %>% 
  rbind(
    czsum %>% 
      filter(samptype == 'in-samp') %>% 
      mutate(
        group = 'base', rbar = rbar / rbarbarbase * 100
      ) %>% select(group, signalname, rbar)
  )

groupsum = group_signal_sum %>% 
  group_by(group) %>% 
  summarize(
    rbarbar = mean(rbar)
    , sd = sd(rbar)
    , se = sd/sqrt(dplyr::n())
  ) %>% 
  filter(group != 'holdper06') %>%   
  mutate(
    group = factor(
      group
      , levels = c('base','holdper12','mescreen','VWforce')
      , labels = c(
        'Original Implementation', 'Annual Rebalancing','ME > NYSE 20 Pct','Value Weighted'
      )
    )
  )

# plot
groupsum %>% 
  ggplot(aes(x = group, y = rbarbar)) +
  geom_bar(stat = 'identity', fill = 'grey', color = 'black') + 
  geom_errorbar(
    aes(ymin = rbarbar-2*se, ymax = rbarbar + 2*se), width = 0.2
  ) +
  chen_theme +
  scale_y_continuous(
    breaks = seq(0,100,20)
    , sec.axis = sec_axis
    (~./1,  breaks = seq(0,100,20)
    )
  ) +
  coord_cartesian(ylim = c(5, 110)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12)) +
  ylab('Grand Mean Return \n (bps Monthly)') +
  theme(
    axis.title.x = element_blank()
  )


ggsave('../results/liq_screen.pdf', width = 8, height = 5, scale = 1.3,  device = cairo_pdf )



## Numbers for paper -------------------------------------------------------

signaldoc %>% 
  filter(Cat.Signal == 'Predictor') %>% 
  filter(!is.na(`Portfolio Period`)) %>% 
  pull(`Portfolio Period`) %>% 
  quantile()


# Sample Time Returns ----------------------------------------------------------------

rollmonths = 12*3

rescale = 0.67

# find rolling stats
tempret = czret %>% 
  select(signalname, date, ret, sampend) %>% 
  mutate(
    samp_time = year(date) + month(date)/12
    - (year(sampend) + month(sampend)/12)
  )
rbar_samp_time = tempret %>% 
  arrange(samp_time) %>% 
  group_by(samp_time) %>% 
  summarize(
    rbar = mean(ret)*100/rescale, nsignal = dplyr::n()
  ) %>% 
  mutate(
    roll_rbar = rollmean(rbar, k = rollmonths, fill = NA, align = 'right')
  )


# big pic stats
czgrand = czret %>% group_by(samptype) %>% summarize(rbarbar = mean(ret)*100/rescale)

guidedat = tibble(
  time = seq(-50,50,0.02)
) %>% 
  mutate(
    rbar = case_when(
      time < 0 ~ 100
      , time < 5 ~ 100 - 12
      , time >= 5 ~ czgrand %>% filter(samptype == 'post-pub') %>% pull(rbarbar)
    )
  )


ggplot(rbar_samp_time, aes(x = samp_time, y = roll_rbar)) +
  geom_vline(xintercept = 0) +
  geom_vline(xintercept = 3, linetype = 'dashed') +  
  geom_line() +
  coord_cartesian(
    xlim = c(-10, 20), ylim = c(-0, 120)
  ) +
  scale_y_continuous(breaks = seq(0,180,20)) +
  scale_x_continuous(breaks = seq(-50,25,5)) +  
  geom_hline(
    yintercept = czgrand %>% filter(samptype == 'in-samp') %>% pull(rbarbar) 
    , color = MATBLUE
  ) +
  geom_hline(
    yintercept = czgrand %>% filter(samptype == 'post-pub') %>% pull(rbarbar) 
    , color = MATRED, linestyle = 'dotted'
  ) +
  geom_hline(
    yintercept = 88
    , color = MATYELLOW
  ) +
  # geom_line(dat = guidedat, aes(x= time, y = rbar)) +
  chen_theme +
  annotate(geom="text",
           label=TeX("Publication Bias")
           , x=-8, y=87, vjust=-1, family = "Palatino Linotype"
  ) +      
  annotate('segment', x=-5, xend = -5, y = 98, yend = 90
           , arrow = arrow(length = unit(0.3,'cm')), size = 0.5) +  
  annotate(geom="text",
           label=TeX('\\Delta Expected Return'), x=-4.5, y=60, vjust=-1, 
           family = "Palatino Linotype"
  ) +
  annotate('segment', x=-1, xend = -1, y = 85, yend = 54
           , arrow = arrow(length = unit(0.3,'cm')), size = 0.5) +    
  ylab('Trailing 3 Years \n Mean Return (bps p.m.)') +
  xlab('Years Since Original Sample Ended') +
  theme(
    axis.title.y = element_text(size = 18)
    , axis.title.x = element_text(size = 18)
  ) 


ggsave('../results/roll_rbar.pdf', height = 4, width = 8, scale = 1, device = cairo_pdf)



## numbers for paper -------------------------------------------------------


czret %>% 
  group_by(signalname) %>% 
  summarize(
    acor1 = cor(ret, dplyr::lag(ret,1), use = 'pairwise.complete.obs')
    , acor2 = cor(ret, dplyr::lag(ret,2), use = 'pairwise.complete.obs')
    , acor3 = cor(ret, dplyr::lag(ret,3), use = 'pairwise.complete.obs')    
    , acor4 = cor(ret, dplyr::lag(ret,4), use = 'pairwise.complete.obs')        
  ) %>% 
  summarize(
    across(-'signalname', ~ mean(.x, na.rm=T))
  )

signaldoc %>% 
  mutate(between = pubdate - sampend) %>% 
  summarize(
    mean(between)/365
  )

# Gregorian Time Returns ----------------------------------------------------------------

rollmonths = 12*3

# find rolling stats
tempret = czret %>% 
  select(signalname, date, ret, sampend) %>% 
  mutate(
    samp_time = year(date) + month(date)/12
    - (year(sampend) + month(sampend)/12)
  )
rbar_time = tempret %>% 
  arrange(date) %>% 
  group_by(date) %>% 
  summarize(
    rbar = mean(ret)*100/.67, nsignal = dplyr::n()
  ) %>% 
  mutate(
    roll_rbar = rollmean(rbar, k = rollmonths, fill = NA, align = 'right')
  )


# big pic stats
mean_old = rbar_time %>% filter(date >= 2000)


ggplot(rbar_time, aes(x = date, y = roll_rbar)) +
  geom_line() 



# YZ Figure ====

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
  xlab(TeX('Absolute t-statistic $|t_i|$')) +
  ylab('Frequency') +
  coord_cartesian(
    xlim = c(0,7), ylim = c(0,0.4)
  )  

ggsave(
  "../results/yan-zheng.pdf",
  width = 12,
  height = 8,
  scale = 0.9,
  device = cairo_pdf
)

