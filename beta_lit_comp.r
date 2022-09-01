# theta comp across lit ---------------------------------------------------
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

dat = rbind(dat.hlz,dat.cz, dat.jkp, dat.ez)


ggplot(dat, aes(x=theta)) +
  geom_histogram(
    aes(group = paper, fill = paper), position = 'identity', alpha = 0.6
    , breaks = seq(-10,10,0.5)
  ) +
  coord_cartesian(
    xlim = c(-10,10)
  )


ggsave('../results/lit-comp.pdf')

dat %>% 
  filter(theta > 0) %>% 
  group_by(paper) %>% 
  summarize(mean(theta))