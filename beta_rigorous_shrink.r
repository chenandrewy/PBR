# Monte carlo shrinkage estimate for HLZ 


# Do stuff -------------------------------------------------------------------------
library(tidyverse)

N = 2e4

p0 = 0.444
lambda = 55.5
SE = 28

dat = tibble(
  v = runif(N) > p0
  , mu = rexp(N, rate = 1/lambda)
  , ep = rnorm(N, 0, SE)
) %>% 
  mutate(
    mu = if_else(v, mu, 0)
    , rbar = mu+ep
    , t = rbar / SE
    , tabs = abs(t)
    , pub = case_when(
                tabs < 1.96 ~ F
                , tabs <= 2.57 ~ runif(1) < 0.5
                , tabs > 2.57 ~ T
            )
  )



dat %>% 
  filter(pub) %>% 
  summarize(
    mean(mu)/mean(rbar)
    , mean(1-mu/rbar)
  )




# jkp stuff ---------------------------------------------------------------

sigsq = 1000^2/12

T = 420

tausq = 43^2

kappa = 1/(1+sigsq/(T*tausq))

kappa

(sigsq/T)^0.5


# with pub bias

tauc_alt = 29
tauw = 21

vtheta = ((tauc_alt^2 + tauw^2)/14^2)^0.5


1/(vtheta^2 + 1)
