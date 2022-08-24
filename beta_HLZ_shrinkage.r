# Monte carlo shrinkage estimate for HLZ 


# Do stuff -------------------------------------------------------------------------
library(tidyverse)

N = 1e4

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

