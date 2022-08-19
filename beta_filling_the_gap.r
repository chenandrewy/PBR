
# Environment -------------------------------------------------------------
source('functions.r')
library(latex2exp)

# plot globals
theme_set(
  theme_minimal() +
    theme(
      text = element_text(family = "Palatino Linotype")
    )
)


# read returns
ret0 = fread('../data/PredictorPortsFull.csv') %>%                                           
  filter(!is.na(ret), port == 'LS')

# read in SignalDoc.csv and isolates to relevant date cols
signaldoc = fread('../data/SignalDoc.csv') %>% 
  mutate(
    signalname = Acronym
    , pubdate = as.Date(paste0(Year, '-12-31'))
    , sampend = as.Date(paste0(SampleEndYear, '-12-31'))
    , sampstart = as.Date(paste0(SampleStartYear, '-01-01'))
  ) %>% 
  arrange(signalname) %>% 
  select(signalname, pubdate, sampend, sampstart)


# merge 
ret1 = ret0 %>%                                                                 ## merging, on signalname, returns and date (sampend, sampstart, pubdate) data
  left_join(signaldoc) %>% 
  mutate(
    samptype = case_when(
      (date >= sampstart) & (date <= sampend) ~ 'in-samp'
      , (date > sampend) & (date <= pubdate) ~ 'out-of-samp'
      , (date > pubdate) ~ 'post-pub'
      , TRUE ~ NA_character_
    )
  ) %>% 
  filter(!is.na(samptype))

# find empirical t-stats
t_emp = ret1 %>% filter(samptype == 'in-samp', port == 'LS') %>% 
  group_by(signalname) %>% 
  summarize(
    tstat = mean(ret)/sd(ret)*sqrt(dplyr::n())
  )  %>% 
  pull(tstat)



# Make plotting data ----------------------------------------------------------------

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
F_emp = ecdf(t_emp)
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
sigfit = sqrt(1+3^2)
dat_fit =data.frame(
  t_mid = mid2
  , prob = (pnorm(t_right2, sd = sigfit) - pnorm(t_left2, sd = sigfit))/
        (1-pnorm(2, sd = sigfit))*rescalefac*(1-F_emp(2))
  , group = 'fit'
) %>% 
  mutate(prob = if_else(t_mid > 2, prob, 0))

# merge and factor group
dat_all = rbind(dat_emp,dat_null,dat_fit) %>% 
  mutate(
    group = factor(group, levels = c('emp','null','fit'))
  )


# Plot --------------------------------------------------------------------


groupdat = tibble(
  group = c('null', 'fit')
  , color = c(MATRED, MATBLUE)
  , labels = c(
    TeX('Null ($\\sigma_\\mu=0$)')
    , TeX('Fit ($\\sigma_\\mu=3$)')
    )
  , linetype = c('dashed','solid')
)


ggplot(dat_all %>%  filter(group == 'emp'), aes(x=t_mid,y=prob)) +
  geom_bar(stat='identity',position='identity',alpha=0.6, aes(fill = group)) +
  scale_fill_manual(
    values = 'gray', labels = 'Published', name = NULL
  ) +    
  geom_line(
    data = dat_all %>% filter(group != 'emp'), aes(color = group, linetype = group)
  ) +
  xlab('t-statistic') +
  ylab('Frequency') +
  scale_color_manual(
    values = groupdat$color, labels = groupdat$labels, name = NULL
  ) + 
  scale_linetype_manual(
    values = groupdat$linetype, labels = groupdat$labels, name = NULL
  ) +
  theme(
    legend.position = c(75,75)/100
    , legend.margin = margin(t = -15, r = 20, b = 0, l = 5),
  )  +
  coord_cartesian(
    xlim = c(0,10), ylim = c(0,0.25)
  )  


ggsave('../results/filling-the-gap.pdf', scale = 1, height = 2.5, width = 5, device = cairo_pdf)
