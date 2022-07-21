# 2022 04 Andrew: estimations on simulated data with correlations

# Based off experiments/qml5*.ret
# takes 20 minutes for n_miniboot = 10
# 3 hours for n_miniboot = 100

# ENVIRONMENT ====
rm(list = ls())
source('0_Settings.r')
tic_all = Sys.time()

# run this if you get an error on import_cz below
# import_cz(dl = T) 

# SETTINGS ====

## estimation settings  ====
set.check = list(
  opt_method = 'two-stage' # 'two-stage' or 'crs-only' or 'pif-grid'
  , opt_list1 = opts.crs(maxeval = 200)
  , opt_list2 = opts.qa(xtol_rel = 1e-3)  
  , model_fam = data.frame(
    mufam   = 'lognormraw' 
    , pif     = c(NA, 0.01, 0.99)
    , mua     = c(NA,    0, 2) # mua = 0 => median = 1.0
    , siga    = c(NA, 0.05, 1) # siga > 1 leads to crazy variances
    , pubfam  = 'stair' # 'stair' or 'piecelin' or 'trunc'  
    , pubpar1 = c(NA,  1/3, 2/3)
    , row.names  = c('base','lb','ub')  
  )  
)

## simulation settings ====

piflist = seq(0.1,0.9,length.out = 3)
rholist = seq(0.1,0.9,length.out = 3)
nportlist = c(200)
n_miniboot = 100


# TEST POINT ESTIMATE ====
tic = Sys.time()

cz_filt_tabs = import_cz(dl=F) %>% filter(tabs >= 1.96) %>% pull(tabs)

est.point = estimate(
  est.set = set.check
  , tabs = cz_filt_tabs
  , par.guess = random_guess(set.check,1,1243)
  , print_level = 3
)
temp = make_stats_pub(cz_filt_tabs,est.point$par) # for now assume signs don't matter
est.point$bias = temp$bias
est.point$fdrloc = temp$fdr_tabs

est.point$par %>% print()
est.point$negloglike %>% print()

p_point = hist_emp_vs_par(cz_filt_tabs,est.point$par)
p_bias = est.point$bias %>% as_tibble() %>% ggplot(aes(x=value)) + 
  geom_histogram() + xlab('bias')
p_fdrloc = est.point$fdrloc %>% as_tibble() %>% ggplot(aes(x=value)) + 
  geom_histogram() + xlab('Prob False')

grid.arrange(p_point,p_bias,p_fdrloc, nrow = 1)

toc = Sys.time()
toc - tic

