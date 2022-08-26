# Make time-series figures #
# # # # # # # # # # # # # #  

# Load packages -----------------------------------------------------------

library(tidyverse)

pathExhibits = '../../Apps/Overleaf/Publication Bias in Asset Pricing/exhibits/'

# Load data ---------------------------------------------------------------

df = readxl::read_xlsx('Data/GoyalWelshZafirov2021.xlsx', sheet = 'Table2')  # for OP t-stats and replicated t-stats


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

df_preds_m = readxl::read_xlsx('Data/PredictorData2021.xlsx', sheet = 'Monthly') %>% 
  mutate_all(as.numeric) %>% 
  transmute(dp = log(D12) - log(Index)
            , dy = log(D12) - lag(log(Index), 1)
            , ep = log(E12) - log(Index)
            , de = log(D12) - log(E12)
            , svar
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

            
df_preds_q = readxl::read_xlsx('Data/PredictorData2021.xlsx', sheet = 'Quarterly') %>% 
  mutate_all(as.numeric) %>% 
  transmute(dp = log(D12) - log(Index)
            , dy = log(D12) - lag(log(Index), 1)
            , ep = log(E12) - log(Index)
            , de = log(D12) - log(E12)
            , svar
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

df_preds_a = readxl::read_xlsx('Data/PredictorData2021.xlsx', sheet = 'Annual') %>% 
  mutate_all(as.numeric) %>% 
  transmute(dp = log(D12) - log(Index)
            , dy = log(D12) - lag(log(Index), 1)
            , ep = log(E12) - log(Index)
            , de = log(D12) - log(E12)
            , svar
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
  geom_text(aes(label = Vrbl)) +
  theme_bw(
    base_size = 15
  ) +
  theme(
    legend.position = c(.8, .25)
    # , text = element_text(size=30, family = "Palatino Linotype")
  ) +
  geom_abline(
    aes(slope = 1, intercept = 0)
  ) +
  labs(x = 't-stat Original Paper'
       , y = 't-stat Replicated')  +
  coord_trans(x='log10', y='log10', xlim = c(1.5, 6), ylim = c(1.0, 6)) +
  scale_x_continuous(breaks=c(1, 2, 3, 4)) +
  scale_y_continuous(breaks=c(1, 2, 3, 4))

ggsave(paste0(pathExhibits, 'fig_ts_scatterOPvsRep.png'), width = 8, height = 6)


# Fact 3: Empirical t-stats are not very large ----------------------------

df %>% 
  ggplot(aes(x = abs(t_OP))) +
  geom_histogram(center = .5, breaks = c(1, 2, 2.5, 3, 3.5, 4, 4.5, 5, 5.5, 6)) +
  geom_vline(xintercept = c(quantile(abs(df$t_OP), c(.5, .9)))) + 
  theme_bw(base_size = 15) +
  labs(x = 'Reported t-stat', y = 'Count')

ggsave(paste0(pathExhibits, 'fig_ts_t_stat.png'), width = 8, height = 6)

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
  ggplot(aes(x = cor)) +
  geom_histogram() +
  facet_wrap(~freq) +
  theme_bw(base_size = 15)

ggsave(paste0(pathExhibits, 'fig_ts_correlations.png'), width = 10, height = 6)

cors %>% 
  group_by(freq) %>% 
  summarise(meanCor = mean(abs(cor)))
    