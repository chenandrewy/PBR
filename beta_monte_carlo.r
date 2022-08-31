# Find MT stats by monte carlo --------------------------------------------------------------


library(data.table)

n = 1e6

# t-dist based on normal fit
sigtheta = 3
nu_g = 4
sig_g = sqrt((nu_g-2)/nu_g)*sigtheta
theta_t = rt(n,nu_g)*sig_g

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
    t = theta + Z, v = theta > 0
  )

# fit
datsum = dat %>% mutate(tgroup = ntile(t,100)) %>% 
  group_by(tgroup) %>% 
  summarize(
    Etheta = mean(theta), Et = mean(t)
  )



Var_t = var(dat$t)

dat %>% 
  filter(t>2) %>% 
  summarize(
    mean(t)
    , mean(theta)
    , mean(theta <= 0)    
    , mean(theta <= 1)
    , mean(theta)
    , mean(theta/t)
  )

lm(theta ~ t, dat) %>% summary()


set.seed(456)


ggplot(
  dat[1:2000,]
  , aes(x=t,y=theta)) +
  geom_point(aes(group = v, color = v)) +
  geom_vline(xintercept = 2) +
  geom_line(
    data = datsum, aes(x=Et, y=Etheta)
  ) +
  coord_cartesian(
    xlim = c(0,10), ylim = c(-5,12)
  )


ggsave('../results/monte-carlo.pdf')