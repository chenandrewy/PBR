library(data.table)


# Find ez stats by monte carlo --------------------------------------------------------------



## simulate ez ----------------------------------------------------------------

n = 1e6
set.seed(456)

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

  

## plot top panel --------------------------------------------------------------------

# settings for both panels here
xlimnum = c(0,6)

ggplot(
  dat[1:2000,]
  , aes(x=tselect,y=theta)) +
  geom_point(aes(group = v, color = v)) +
  geom_vline(xintercept = 2) +
  geom_line(
    data = datsum, aes(x=Etselect, y=Etheta)
  ) +
  geom_abline(slope = 1) +
  coord_cartesian(
    xlim = xlimnum, ylim = c(-2,10)
  ) +
  theme(
    legend.position = c(25,75)/100
  )


ggsave('../results/monte-carlo-ez.pdf')


## numbers for text --------------------------------------------------------

dat %>% 
  filter(tselect>2) %>% 
  summarize(
    mean(tselect)
    , mean(theta)
    , mean(theta <= 0)    
    , mean(theta) / mean(tselect)
  )




# Find HLZ stats by monte carlo --------------------------------------------------------------



## simulate hlz ----------------------------------------------------------------#

n = 1e6
set.seed(456)

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

## plot top panel ==== #

# settings for both panels here
xlimnum = c(-2,6)
nplot = 1000
  
set.seed(24)
ggplot(
  dat[sample(1:n,nplot),]
  , aes(x=tselect,y=theta)) +
  geom_point(aes(group = v, color = v)) +
  geom_vline(xintercept = 2) +
  geom_vline(xintercept = hurdle_bonf05) +
  geom_line(
    data = datsum, aes(x=Etselect, y=Etheta)
  ) +
  geom_abline(slope = 1) +
  coord_cartesian(
    xlim = xlimnum, ylim = c(-2,10)
  ) +
  theme(
    legend.position = c(25,75)/100
  )


ggsave('../results/monte-carlo-hlz.pdf')



## numbers for text --------------------------------------------------------

dat %>% 
  filter(tselect>2) %>% 
  summarize(
    mean(tselect)
    , mean(theta)
    , mean(!v)    
    , mean(theta) / mean(tselect)
  )

