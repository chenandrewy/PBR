# Make time-series figures #
# # # # # # # # # # # # # #  

# Load packages -----------------------------------------------------------

library(tidyverse)
source('setup.r')

pathExhibits = '../../Apps/Overleaf/Publication Bias in Asset Pricing/exhibits/' # for Tom
pathExhibits = '../results/' # for general purpose use

# Load data ---------------------------------------------------------------

df = readxl::read_xlsx('Data_ts/GoyalWelchZafirov2021.xlsx', sheet = 'Table2') %>%   # for OP t-stats and replicated t-stats
  mutate(
    t_rep_sign = t_rep * sign(t_OP)
    , t_OP_sign = t_OP * sign(t_OP)
  ) 
# TS predictors from Goyal and Welsh (2008)

# the dividend-price ratio (d/p), the dividend-yield (d/y), the earnings-
#   price ratio (e/p), the dividend-payout ratio (d/e), as in Campbell and Shiller (1988); stock
# volatility (svar), as in Guo (2006); book-market (b/m), as in Kothari and Shanken (1997) and
# Pontiff and Schall (1998); net issuing activity (ntis), as in Boudoukh, Michaely, Richardson, and
# Roberts (2007); equity issuing activity (eqis), as in Baker and Wurgler (2000); the T-Bill rate
# (tbl), as in Campbell (1987); the long-term yield (lty), the long-term bond rate of return (ltr),
# the term-spread (tms), the default yield (dfy), the default rate of return (dfr), as in Fama and
# French (1989), the inflation rate (infl), as in Fama and Schwert (1977); private investment (i/k),
# as in Cochrane (1991), and “cay, ” 

df_preds_m = readxl::read_xlsx('Data_ts/PredictorData2021.xlsx', sheet = 'Monthly') %>% 
  mutate_all(as.numeric) %>% 
  transmute(dp = log(D12) - log(Index)
            , dy = log(D12) - lag(log(Index), 1)
            , ep = log(E12) - log(Index)
            , de = log(D12) - log(E12)
            , svar
            , csp
            , bm = `b/m`
            , ntis
             #, eqis
            , tbl
            , lty
            , ltr
            , tms = lty - tbl
            , dfy = BAA - AAA
            , dfr = corpr - ltr
            , infl = -1*infl
            # ik
            # cay
  ) %>% 
  mutate(tbl = -1*tbl)

            
df_preds_q = readxl::read_xlsx('Data_ts/PredictorData2021.xlsx', sheet = 'Quarterly') %>% 
  mutate_all(as.numeric) %>% 
  transmute(dp = log(D12) - log(Index)
            , dy = log(D12) - lag(log(Index), 1)
            , ep = log(E12) - log(Index)
            , de = log(D12) - log(E12)
            , svar
            , csp 
            , bm = `b/m`
            , ntis
            #, eqis
            , tbl = -1*tbl
            , lty
            , ltr
            , tms = lty - tbl
            , dfy = BAA - AAA
            , dfr = corpr - ltr
            , infl = -1*infl
            , ik = -1*ik
            # cay
  ) %>% 
  mutate(tbl = -1*tbl)

df_preds_a = readxl::read_xlsx('Data_ts/PredictorData2021.xlsx', sheet = 'Annual') %>% 
  mutate_all(as.numeric) %>% 
  transmute(dp = log(D12) - log(Index)
            , dy = log(D12) - lag(log(Index), 1)
            , ep = log(E12) - log(Index)
            , de = log(D12) - log(E12)
            , svar
            , csp
            , bm = `b/m`
            , ntis
            , eqis = -1*eqis
            , tbl
            , lty
            , ltr
            , tms = lty - tbl
            , dfy = BAA - AAA
            , dfr = corpr - ltr
            , infl = -1*infl
            , ik = -1*ik
            , cay
  ) %>% 
  mutate(tbl = -1*tbl)

# Fact 1: OPs can be replicated -------------------------------------------

df %>% 
  ggplot(aes(x=abs(t_OP), y = abs(t_rep))) +
  chen_theme +
  geom_text(aes(label = Vrbl, family="Palatino Linotype"), size = 6) +
  theme(
    legend.position = c(.8, .25)
    # , text = element_text(size=30, family = "Palatino Linotype")
  ) +
  geom_abline(
    aes(slope = 1, intercept = 0, linetype = "45 deg. line")
  ) +
  labs(x = 't-stat Original Paper'
       , y = 't-stat Replicated')  +
  coord_trans(x='log10', y='log10', xlim = c(1.5, 6), ylim = c(1.0, 6)) +
  scale_x_continuous(breaks=c(1, 2, 3, 4, 5)) +
  scale_y_continuous(breaks=c(1, 2, 3, 4, 5))

ggsave(
  "../results/ts_raw_op.pdf",
  width = 12,
  height = 8,
  device = cairo_pdf
)


# Fact 3: Empirical t-stats are not very large ----------------------------

# set up 
edge = seq(-1,20,0.5)
t_left = edge[1:(length(edge)-1)]
t_right = edge[2:length(edge)]
mid = t_left + diff(edge)/2

# cross sectional frequency
F_emp = ecdf(czsum$tstat[czsum$samptype == "in-samp"])
cs_emp_dat = data.frame(
  t_mid = mid
  , prob = (F_emp(t_right) - F_emp(t_left))
  , group = 'cs'
)

# time series frequency
ts_emp = ecdf(df$t_rep_sign)
ts_emp_dat = data.frame(
  t_mid = mid,
  prob = (ts_emp(t_right) - ts_emp(t_left)),
  group = "ts"
)

# EB model fit
mean_t_pub <- function(sigma) {
  ((dnorm(2/sigma, mean = 0, sd = 1)) / (1-pnorm(2/sigma, mean =0, sd = 1))) * sigma 
}  
vartheta = 1.5^2
mean_t_pub(sqrt(vartheta + 1))
mean(df$t_rep_sign[df$t_rep_sign>2])
1/(vartheta + 1)


sigfit = sqrt(vartheta + 1)
dat_fit = data.frame(
  t_mid = mid
  , prob = (pnorm(t_right, sd = sigfit) - pnorm(t_left, sd = sigfit))/
    (1-pnorm(2, sd = sigfit))*(1-F_emp(2))
  , group = 'fit'
) %>% 
  mutate(prob = if_else(t_mid > 2, prob, 0))



# rbind cs and ts
dat_all_ts_cs = rbind(ts_emp_dat, cs_emp_dat, dat_fit) %>% 
  mutate(
    group = factor(group, levels = c('ts', 'cs', 'fit'))
  ) %>% 
  filter(group %in% c('ts','cs'))

ggplot(dat_all_ts_cs, aes(x=t_mid,y=prob, fill = as.factor(group))) +
  geom_bar(data = dat_all_ts_cs, stat = 'identity', position='identity',alpha=0.6) +
  scale_fill_manual(
    values = c(MATRED, 'grey', MATBLUE)
    , labels = c("Equity Premium", "Cross-Section", 'fit'), name = NULL
  ) +
  scale_x_continuous(limits=c(0, 15), breaks = seq(-10,20,2)) +
  chen_theme +
  xlab("Reported t-stat") +
  ylab("Frequency") 


ggsave(
  "../results/ts_t_stat.pdf",
  width = 12,
  height = 8,
  device = cairo_pdf
)


# Table
tmp = hist(abs(df$t_OP), breaks = c(1, 2, 3, 4, 5, 6), plot = F)

tbl = tibble(breaks = tmp$breaks + 1,
             NLarger = nrow(df) - c(cumsum(tmp$counts), nrow(df)),
             Share = (nrow(df) - c(cumsum(tmp$counts), nrow(df)))/nrow(df))




# Fact 3: with fit  ----------------------------------------------------------------
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
F_emp = ecdf(df$t_rep_sign)
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

 
# fit
sigfit = sqrt(1+1.5^2)
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
    group = factor(group, levels = c('emp','null','fit'))
  )


## plot --------------------------------------------------------------------
groupdat = tibble(
  group = c('null',  'fit')
  , color = c(MATRED, MATBLUE)
  , labels = c(
    TeX('$\\sigma_\\theta$ = 0 (Null)')
    , TeX('$\\hat{\\sigma}_\\theta$ = $1.5$ (Fit)')
  )
  , linetype = c('longdash','solid')
)

ggplot(dat_all %>%  filter(group == 'emp'), aes(x=t_mid,y=prob)) +
  geom_bar(stat = 'identity', position='identity',alpha=0.6, aes(fill = group)) +
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
  coord_cartesian(xlim = c(0,8), ylim = c(0,0.5))  +
  scale_x_continuous(breaks = seq(-10,20,2))

ggsave('../results/ts-filling-the-gap.pdf', 
       width = 12,
       height = 8,
       device = cairo_pdf
)


## numbers for paper -------------------------------------------------------
1/sigfit^2



# Fact 4: Pairwise correlations are high ----------------------------------

tmp = cor(as.matrix(df_preds_m), use = "pairwise.complete.obs")
cor_m =  tmp[lower.tri(tmp)] %>% c()

tmp = cor(as.matrix(df_preds_q), use = "pairwise.complete.obs")
cor_q =  tmp[lower.tri(tmp)] %>% c()

tmp = cor(as.matrix(df_preds_a), use = "pairwise.complete.obs")
cor_a =  tmp[lower.tri(tmp)] %>% c()

cors = tibble(cor = cor_m, freq = 'monthly') %>% 
  bind_rows(tibble(cor = cor_q, freq = 'quarterly')) %>% 
  bind_rows(tibble(cor = cor_a, freq = 'annually'))
    
cors %>% 
  ggplot(aes(x = cor, color = freq)) +
  geom_histogram(fill = MATBLUE, color = NA, show.legend = FALSE, size = 1) +
  facet_wrap(~freq) +
  chen_theme +
  theme(
    text = element_text(size = 24)
  ) +
  # scale_color_manual(MATBLUE) +
  labs(x = 'Correlation', y = 'Count')

ggsave(
  "../results/fig_ts_correlations.pdf",
  width = 12,
  height = 8,
  device = cairo_pdf
)

cors %>% 
  group_by(freq) %>% 
  summarise(meanCor = mean(abs(cor)))
    