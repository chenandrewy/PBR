n# SETUP ====
rm(list = ls())
tic = Sys.time()
dir.create('../results/')
source('functions.r')


set.seed(1057)

## graph defaults ====
MATBLUE = rgb(0,0.4470,0.7410)
MATRED = rgb(0.8500, 0.3250, 0.0980)
MATYELLOW = rgb(0.9290, 0.6940, 0.1250)

chen_theme = theme_minimal() +
  theme(
    text = element_text(family = "Palatino Linotype")
    , panel.border = element_rect(colour = "black", fill=NA, size=1)    
    
    # Font sizes
    , axis.title.x = element_text(size = 26),
    axis.title.y = element_text(size = 26),
    axis.text.x = element_text(size = 22),
    axis.text.y = element_text(size = 22),
    legend.text = element_text(size = 18),
    
    # Tweaking legend
    legend.position = c(0.7, 0.8),
    legend.text.align = 0,
    legend.background = element_rect(fill = "white", color = "black"),
    legend.margin = margin(t = 5, r = 20, b = 5, l = 5), 
    legend.key.size = unit(1.5, "cm")
    , legend.title = element_blank()    
  ) 

# Import/Prepare Data ====
cz_all = fread("../data/PredictorPortsFull.csv")
signaldoc = fread('../data/SignalDoc.csv') %>% 
  rename(signalname = Acronym) %>% 
  mutate(
    pubdate = as.Date(paste0(Year, '-12-31'))
    , sampend = as.Date(paste0(SampleEndYear, '-12-31'))
    , sampstart = as.Date(paste0(SampleStartYear, '-01-01'))
  )

cz_alt = rbind(
  fread('../data/PredictorAltPorts_HoldPer_12.csv') %>% 
    mutate(group = 'holdper12')
  ,   fread('../data/PredictorAltPorts_HoldPer_6.csv') %>% 
    mutate(group = 'holdper06')
  , fread('../data/PredictorAltPorts_LiqScreen_ME_gt_NYSE20pct.csv') %>% 
    mutate(group = 'mescreen')  
  ,   fread('../data/PredictorAltPorts_LiqScreen_VWforce.csv') %>% 
    mutate(group = 'VWforce')
) %>% 
  filter(port == 'LS', !is.na(ret)) %>% 
  select(group, signalname, date, ret)


## CZRET ====
czret = cz_all %>%         
  filter(port == 'LS', !is.na(ret)) %>% 
  left_join(signaldoc) %>% 
  mutate(
    samptype = case_when(
      (date >= sampstart) & (date <= sampend) ~ 'in-samp'
      , (date > sampend) & (date <= sampend %m+% months(36)) ~ 'out-of-samp'
      , (date > pubdate) ~ 'post-pub'
      , TRUE ~ 'NA_character_'
    )
  ) %>% 
  select(signalname, date, ret, samptype, sampstart, sampend)

cz_alt = cz_alt %>% 
  left_join(signaldoc)  %>% 
  mutate(
    samptype = case_when(
      (date >= sampstart) & (date <= sampend) ~ 'in-samp'
      , (date > sampend) & (date <= sampend %m+% months(36)) ~ 'out-of-samp'
      , (date > pubdate) ~ 'post-pub'
      , TRUE ~ 'NA_character_'
    )
  ) %>% 
  select(group, signalname, date, ret, samptype, sampstart, sampend)
  


## CZSUM ====
czsum = czret %>%
  group_by(signalname, samptype) %>%
  summarize(
    rbar = mean(ret)
    , vol = sd(ret)
    , tstat = mean(ret)/sd(ret)*sqrt(dplyr::n())
  )



# Liquidity Screen bars -------------------------------------------------------

rbarbarbase = czsum %>% filter(samptype == 'in-samp') %>% 
  ungroup() %>% 
  summarize(rbarbar = mean(rbar)) %>% 
  pull(rbarbar)


group_signal_sum = cz_alt %>% 
  filter(samptype == 'in-samp') %>% 
  group_by(group, signalname) %>% 
  summarize(
    rbar = mean(ret)/rbarbarbase*100
  ) %>% 
  mutate(
    group = factor(group)
  ) %>% 
  rbind(
    czsum %>% 
      filter(samptype == 'in-samp') %>% 
      mutate(
        group = 'base', rbar = rbar / rbarbarbase * 100
      ) %>% select(group, signalname, rbar)
  )

groupsum = group_signal_sum %>% 
  group_by(group) %>% 
  summarize(
    rbarbar = mean(rbar)
    , sd = sd(rbar)
    , se = sd/sqrt(dplyr::n())
  ) %>% 
  filter(group != 'holdper06') %>%   
  mutate(
    group = factor(
      group
      , levels = c('base','holdper12','mescreen','VWforce')
      , labels = c(
        'Original Implementation', 'Annual Rebalancing','ME > NYSE 20 Pct','Value Weighted'
      )
    )
  )
  
# plot
groupsum %>% 
  ggplot(aes(x = group, y = rbarbar)) +
  geom_bar(stat = 'identity', fill = 'grey', color = 'black') + 
  geom_errorbar(
    aes(ymin = rbarbar-2*se, ymax = rbarbar + 2*se), width = 0.2
  ) +
  chen_theme +
  scale_y_continuous(
    breaks = seq(0,100,20)
    , sec.axis = sec_axis
    (~./1,  breaks = seq(0,100,20)
      )
  ) +
  coord_cartesian(ylim = c(5, 110)) +
  scale_x_discrete(labels = function(x) stringr::str_wrap(x, width = 12)) +
  ylab('Grand Mean Return \n (bps Monthly)') +
  theme(
    axis.title.x = element_blank()
  )


ggsave('../results/liq_screen.pdf', width = 8, height = 5, scale = 1.3,  device = cairo_pdf )



# Sample Time Returns ----------------------------------------------------------------

rollmonths = 12*5

rescale = 0.67

# find rolling stats
tempret = czret %>% 
  select(signalname, date, ret, sampend) %>% 
  mutate(
    samp_time = year(date) + month(date)/12
      - (year(sampend) + month(sampend)/12)
  )
rbar_samp_time = tempret %>% 
  arrange(samp_time) %>% 
  group_by(samp_time) %>% 
  summarize(
    rbar = mean(ret)*100/rescale, nsignal = dplyr::n()
  ) %>% 
  mutate(
    roll_rbar = rollmean(rbar, k = rollmonths, fill = NA, align = 'right')
  )


# big pic stats
czgrand = czret %>% group_by(samptype) %>% summarize(rbarbar = mean(ret)*100/rescale)

guidedat = tibble(
  time = seq(-50,50,0.02)
) %>% 
  mutate(
    rbar = case_when(
      time < 0 ~ 100
      , time < 5 ~ 100 - 12
      , time >= 5 ~ czgrand %>% filter(samptype == 'post-pub') %>% pull(rbarbar)
    )
  )


ggplot(rbar_samp_time, aes(x = samp_time, y = roll_rbar)) +
  geom_line() +
  coord_cartesian(
    xlim = c(-10, 20), ylim = c(-0, 120)
  ) +
  scale_y_continuous(breaks = seq(0,180,20)) +
  scale_x_continuous(breaks = seq(-50,25,5)) +  
  geom_hline(yintercept = czgrand %>% filter(samptype == 'in-samp') %>% pull(rbarbar) ) +
  geom_hline(yintercept = czgrand %>% filter(samptype == 'post-pub') %>% pull(rbarbar) ) +
  geom_hline(yintercept = 88) +
  # geom_line(dat = guidedat, aes(x= time, y = rbar)) +
  geom_vline(xintercept = 0) +
  chen_theme +
  # annotate(geom="text",
  #          label=TeX("In-Sample Mean")
  #          , x=10, y=100, vjust=-1,
  #          family = "Palatino Linotype"
  # ) +
  annotate(geom="text",
           label=TeX("Publication Bias")
           , x=-5, y=87, vjust=-1, family = "Palatino Linotype"
  ) +      
  annotate('segment', x=-1, xend = -1, y = 100, yend = 90
           , arrow = arrow(length = unit(0.3,'cm')), size = 0.5) +  
  # annotate(geom="text",
  #          label=TeX("Post-Pub Mean")
  #          , x=10, y=40, vjust=-1,
  #          family = "Palatino Linotype"
  # ) +
  annotate(geom="text",
           label=TeX('\\Delta E(Return)'), x=2, y=60, vjust=-1, 
           family = "Palatino Linotype"
  ) +
  annotate('segment', x=5, xend = 5, y = 85, yend = 54
           , arrow = arrow(length = unit(0.3,'cm')), size = 0.5) +    
  ylab('Trailing 5 Years \n Mean Return (bps p.m.)') +
  xlab('Years Since Original Sample Ended') +
  theme(
    axis.title.y = element_text(size = 20)
    , axis.title.x = element_text(size = 20)
  )

  
ggsave('../results/roll_rbar.pdf', height = 4, width = 8, scale = 1, device = cairo_pdf)



# Gregorian Time Returns ----------------------------------------------------------------

rollmonths = 12*3

# find rolling stats
tempret = czret %>% 
  select(signalname, date, ret, sampend) %>% 
  mutate(
    samp_time = year(date) + month(date)/12
    - (year(sampend) + month(sampend)/12)
  )
rbar_time = tempret %>% 
  arrange(date) %>% 
  group_by(date) %>% 
  summarize(
    rbar = mean(ret)*100/.67, nsignal = dplyr::n()
  ) %>% 
  mutate(
    roll_rbar = rollmean(rbar, k = rollmonths, fill = NA, align = 'right')
  )


# big pic stats
mean_old = rbar_time %>% filter(date >= 2000)


ggplot(rbar_time, aes(x = date, y = roll_rbar)) +
  geom_line() 




