
# t too big table ---------------------------------------------------------

## Import CZ ====

# All signals
doc =  fread("../data/SignalDoc.csv") %>%
  rename(signalname = Acronym) %>%
  select(-c(`Detailed Definition`, `Notes`))


# Chen-Zimmerman dataset
czret = fread(paste0("../data/PredictorPortsFull.csv"))  %>% 
  select(signalname, port, date, ret) %>%
  left_join(
    doc %>% select(signalname, SampleStartYear, SampleEndYear)
    , by = 'signalname'
  ) %>%
  filter(
    (year(date) >= SampleStartYear) &  (year(date) <= SampleEndYear)
    , port == 'LS', !is.na(ret)
  )

czsum = czret %>% 
  group_by(signalname) %>% 
  summarize(
    tstat = mean(ret)/sd(ret)*sqrt(n()) %>% abs()
  )

## Import YZ ====
library(haven)
temp = read_sas('../data-YZ/Yan_Zheng_RFS_Data.sas7bdat')

yzsum = temp %>%
  mutate(
    signalname = paste(transformation, fsvariable, sep = '.')
  ) %>%
  transmute(
    signalname, date = DATE, ret = 100*ddiff_ew
  ) %>% 
  filter(!is.na(ret)) %>% 
  group_by(signalname) %>% 
  summarize(
    tstat = mean(ret)/sd(ret)*sqrt(n()) %>% abs()
  )


## Settings ====

# sig_pct = c(5, 1, 0.1, 0.01, 0.001, 0.0001, 0.00001)
# t_left = -1*qnorm(sig_level/2/100)

t_left = c(seq(2,9,1), Inf)
sig_pct = pnorm(-t_left)*200
nboot = 10*1000

## Bootstrap ====

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
  eboot = e[dateselect, ]
  
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
ndatemax = dim(e)[1]
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

## Make table ====
Fcz = ecdf(abs(czsum$tstat))
Fyz = ecdf(abs(yzsum$tstat))
Ncz = length(czsum$tstat)
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

write.csv(tab_wide, '../results/t-too-big-v2.csv', col.names = NULL)

tab_wide