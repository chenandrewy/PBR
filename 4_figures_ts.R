# Make time-series figures #
# # # # # # # # # # # # # #  

# Load packages -----------------------------------------------------------

library(tidyverse)

pathExhibits = '../../Apps/Overleaf/Publication Bias in Asset Pricing/exhibits/' # for Tom
pathExhibits = '../results/' # for general purpose use

# Load data ---------------------------------------------------------------

df = readxl::read_xlsx('Data_ts/GoyalWelshZafirov2021.xlsx', sheet = 'Table2')  # for OP t-stats and replicated t-stats


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
# ggplot(dat_all %>%  filter(group == 'emp'), aes(x=t_mid)) +
#   geom_histogram(position='identity',alpha=0.6,breaks = c(0, 1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6)) +
#   chen_theme

# cross sectional frequency
F_emp = ecdf(czsum$tstat[czsum$samptype == "in-samp"])
cs_emp_dat = data.frame(
  t_mid = mid
  , prob = (F_emp(t_right) - F_emp(t_left))
  , group = 'cs'
)

# time series emp
ts_emp = ecdf(df$t_OP)
ts_emp_dat = data.frame(
  t_mid = mid,
  prob = (ts_emp(t_right) - ts_emp(t_left)),
  group = "ts"
)

# dat_all_ts_cs
dat_all_ts_cs = rbind(ts_emp_dat, cs_emp_dat) %>% 
  mutate(
    group = factor(group, levels = c('cs', 'ts'))
  )

ggplot(dat_all_ts_cs, aes(x=t_mid,y=prob, fill = as.factor(group))) +
  geom_bar(data = subset(dat_all_ts_cs, group == "ts"), stat = 'identity', position='identity',alpha=0.6) +
  geom_bar(data = subset(dat_all_ts_cs, group == "cs"), stat = 'identity', position='identity',alpha=0.6) +
  scale_fill_manual(
    values = c("grey", MATBLUE), labels = c("Time Series", "Cross Sectional"), name = NULL
  ) +
  scale_x_continuous(limits=c(0, 15)) +
  chen_theme +
  xlab("Reported t-stat") +
  ylab("Frequency")

ggsave(paste0(pathExhibits, 'fig_ts_t_stat.png'), width = 8, height = 6)
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
    