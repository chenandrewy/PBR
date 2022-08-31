# Find MT stats by monte carlo --------------------------------------------------------------


library(data.table)



## simulate hlz ----------------------------------------------------------------

n = 1e6
set.seed(456)

# hlz baseline
v  = runif(n) > 0.444
mu = rexp(n, 1/55.5*se)
mu[!v] = 0
se = 1500/sqrt(12)/sqrt(240)
theta_hlz = mu

# simulate
dat = data.table(
  Z = rnorm(n), theta = theta_hlz
) %>% 
  mutate(
    t = theta + Z, v = theta > 0, tabs = abs(t)
  )

# fit shrinkage
datsum = dat %>% mutate(tgroup = ntile(tabs,100)) %>% 
  group_by(tgroup) %>% 
  summarize(
    tabs_left = min(tabs), tabs_right = max(tabs)
    , Etabs = mean(tabs), Etheta = mean(theta), n = dplyr::n(), nfalse = sum(!v)
  )

# fit "cdf"
datsum = datsum %>% 
  arrange(-Etabs) %>% 
  mutate(
    nfalse_cum = cumsum(nfalse)
    , n_cum = cumsum(n)
    , fdr_tabs_left = nfalse_cum / n_cum * 100
  ) %>% 
  arrange(Etabs)

# find hurdles
hurdle_05 = min(datsum$tabs_left[which(datsum$fdr_tabs_left < 5)])
hurdle_01 = min(datsum$tabs_left[which(datsum$fdr_tabs_left < 1)])
  

## plot top panel --------------------------------------------------------------------

# settings for both panels here
xlimnum = c(0,6)

ggplot(
  dat[1:2000,]
  , aes(x=tabs,y=theta)) +
  geom_point(aes(group = v, color = v)) +
  geom_vline(xintercept = 2) +
  geom_vline(xintercept = hurdle_05) +
  geom_vline(xintercept = hurdle_01) +
  geom_line(
    data = datsum, aes(x=Etabs, y=Etheta)
  ) +
  geom_abline(slope = 1) +
  coord_cartesian(
    xlim = xlimnum, ylim = c(-2,10)
  ) +
  theme(
    legend.position = c(25,75)/100
  )


ggsave('../results/monte-carlo-top.pdf')


## plot bottom panel -------------------------------------------------------

ggplot(datsum, aes(x=tabs_left, y=fdr_tabs_left)) +
  geom_line() +
  geom_vline(xintercept = 2) +
  geom_vline(xintercept = hurdle_05) +
  geom_vline(xintercept = hurdle_01) +
  coord_cartesian(
    xlim = xlimnum, ylim = c(0, 20)
  )

ggsave('../results/monte-carlo-bottom.pdf')

## numbers for text --------------------------------------------------------



dat %>% 
  filter(tabs>2) %>% 
  summarize(
    mean(tabs)
    , mean(theta)
    , mean(theta <= 0)    
    , mean(theta) / mean(tabs)
  )

dat %>% 
  filter(t>2) %>% 
  summarize(
    mean(t)
    , mean(theta)
    , mean(theta <= 0)    
    , mean(theta) / mean(t)
  )


# lm(theta ~ t, dat) %>% summary()


