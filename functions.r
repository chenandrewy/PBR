# LIBRARIES AND PATHS ====
library(tidyverse)

# for data download
library(data.table)
library(googledrive)
library(readxl)
library(writexl)
library(nloptr)
library(distr) 
library(ggplot2)
library(gridExtra)
library(lubridate)

dir.create('intermediate/', showWarnings = F)
dir.create('output/', showWarnings = F)

# root of March 2022 
pathRelease = 'https://drive.google.com/drive/folders/1O18scg9iBTiBaDiQFhoGxdn4FdsbMqGo'

# ESTIMATION FUNCTIONS ====

## Optimization settings ====

opts.crs = function(
  print_level = 0, maxeval = 200, xtol_rel = 0.001, ranseed = 1142
){
  list(
    algorithm = 'NLOPT_GN_CRS2_LM' 
    , xtol_rel = xtol_rel
    , maxeval = maxeval
    , ranseed = ranseed
    , print_level = print_level    
  )
} # end opts.crs  


opts.qa = function(
  print_level = 3, maxeval = 1000, xtol_rel = 1e-3, ftol_rel = 1e-8
){
  list(
    algorithm = 'NLOPT_LN_BOBYQA' 
    , xtol_rel = xtol_rel
    , ftol_rel = ftol_rel
    , maxeval = maxeval
    , print_level = print_level    
  )
} # end opts.qa

## Notation ====

# function: translate par df to parvec vector form
par2parvec = function(par, namevec){
  parvec = par %>% select(all_of(namevec)) %>% as.numeric()
  return = parvec
} # end par2parvec

parvec2par = function(parvec,parbase,namevec){
  par = parbase
  par[ , namevec] = parvec %>% c() %>% t() %>% as.numeric
  return(par)
} # end parvec2par

## Functional Forms ====

# prob of publication
pub_prob = function(tt,par){
  
  if (par$pubfam == 'piecelin'){
    
    prob = 0.5 + par$pubpar2*(tt-par$pubpar1)
    prob[prob<0] = 0; prob[prob>1] = 1
    return = prob    
    
  } else if (par$pubfam == 'logistic'){
    
    # this family may be required if mutrue is student's t
    prob = 1/(1+exp(-par$pubpar2*(tt-par$pubpar1)))
    
  } else if (par$pubfam == 'stair'){
    prob = (tt > 1.96 & tt <= 2.58)*par$pubpar1 +(tt>2.58)*1
    
  } else if (par$pubfam == 'stair2'){
    prob = 
      (tt > 0.00 & tt <= 1.50)*par$pubpar3 +
      (tt > 1.50 & tt <= 1.96)*par$pubpar2 +
      (tt > 1.96 & tt <= 2.58)*par$pubpar1 + 
      (tt>2.58)*1
    
  } else if (par$pubfam == 'stair3'){
    prob = 
      (tt > 0.00 & tt <= 1.00)*par$pubpar3 +
      (tt > 1.00 & tt <= 1.96)*par$pubpar2 +
      (tt > 1.96 & tt <= 2.58)*par$pubpar1 + 
      (tt>2.58)*1    
    
  } else if (par$pubfam == 'trunc'){
    prob = 1*(tt > par$pubpar1)
    
  } # end if par$pubfam
  
  
} # end pub_prob


# make latent density object with distr
make_tabs_object = function(par){
  
  # first make latent t object
  #   mixture normal case
  if (par$mufam == 'mix-norm'){
    
    t.o = UnivarMixingDistribution(
      Norm(0,1), Norm(par$mua,sqrt(par$siga^2+1)), Norm(par$mub,sqrt(par$sigb^2+1))
      , mixCoeff = c(par$pif, (1-par$pif)*par$pia, (1-par$pif)*(1-par$pia))
    )
    
    #   generalized student's t
  } else if (par$mufam == 't'){
    
    tempscale = par$siga*sqrt((par$nua-2)/par$nua)
    mutrue.o = tempscale*Td(df = par$nua) + par$mua
    
    t.o = UnivarMixingDistribution(
      Norm(0,1), mutrue.o + Norm(0,1)
      , mixCoeff = c(par$pif, (1-par$pif))
    )
    
    #   lognormal
  } else if (par$mufam == 'lognorm'){
    
    #    translate using wikipedia formulas
    tempmu = log(par$mua^2/sqrt(par$siga^2+par$mua^2))
    tempsig = log(1 + par$siga^2/par$mua^2)    
    # hack
    tempsig = max(tempsig, 0.1)    
    mutrue.o = Lnorm(tempmu, tempsig)
    
    t.o = UnivarMixingDistribution(
      Norm(0,1), mutrue.o  + Norm(0,1)
      , mixCoeff = c(par$pif, (1-par$pif))
    )    
    
  } else if (par$mufam == 'lognormraw'){
    
    t.o = UnivarMixingDistribution(
      Norm(0,1), Lnorm(par$mua, par$siga)  + Norm(0,1)
      , mixCoeff = c(par$pif, (1-par$pif))
    )
    
  } else if (par$mufam == 'exp'){
    
    t.o = UnivarMixingDistribution(
      Norm(0,1), Exp(1/par$mua)  + Norm(0,1)
      , mixCoeff = c(par$pif, (1-par$pif))
    )   
    
    
  } # end if par$mufam
  
  # take absolute value
  tabs.o = abs(t.o)
  
} # end make_tabs_object

# make mutrue object with distr
#   should only be used in sim
#   separate from tabs object to improve accuracy (hopefully)
make_mutrue_object = function(par){
  
  
  # first make latent t object
  #   mixture normal case
  if (par$mufam == 'mix-norm'){
    
    mutrue.o = UnivarMixingDistribution(
      Norm(par$mua,par$siga), Norm(par$mub,par$sigb)
      , mixCoeff = c(par$pia, (1-par$pia))
    )
    
    #   generalized student's t
  } else if (par$mufam == 't'){
    
    tempscale = par$siga*sqrt((par$nua-2)/par$nua)
    mutrue.o = tempscale*Td(df = par$nua) + par$mua
    
    #   log normal
  } else if (par$mufam == 'lognorm'){
    
    #    translate using wikipedia formulas
    tempmu = log(par$mua^2/sqrt(par$siga^2+par$mu^2))
    tempsig = log(1 + par$siga^2/par$mua^2)
    
    # hack
    tempsig = max(tempsig, 0.1)
    
    mutrue.o = Lnorm(tempmu, tempsig)
    
    #   log normal raw
  } else if (par$mufam == 'lognormraw'){
    
    mutrue.o =   Lnorm(par$mua, par$siga)
    
    #   exp
  } else if (par$mufam == 'exp'){
    
    mutrue.o = Exp(1/par$mua)
    
  } # end if par$mufam
  
  return = mutrue.o
  
} # end make_tabs_object

negloglike = function(tabs,par){
  
  # individual likes
  #   latent density
  tabs.o = make_tabs_object(par)
  
  #   observed density
  #   integrate in pieces for accuracy
  denom1 = integrate(function (tt) d(tabs.o)(tt)*pub_prob(tt,par), 0, 1.96)$value
  denom2 = integrate(function (tt) d(tabs.o)(tt)*pub_prob(tt,par), 1.96, Inf)$value  
  denom = denom1 + denom2
  
  singlelikes = d(tabs.o)(tabs)*pub_prob(tabs,par)/denom
  
  # clean up and take mean    
  singlelikes[singlelikes<=0] = 1e-6
  return = -1*mean(log(singlelikes))
  
} # end f_obs

## Estimation  ====

random_guess = function(est.set, nguess, seed = NULL){
  
  set.seed(seed)
  
  # first replicate base n times
  par.guess.all = est.set$model_fam['base',]
  par.guess.all=  do.call("rbind", replicate(nguess, par.guess.all, simplify = FALSE))
  
  # then randomize est_names columns uniform between lb and ub
  par.lb = est.set$model_fam['lb',]  
  par.ub = est.set$model_fam['ub',]  
  
  est_names = est.set$model_fam['base', is.na(est.set$model_fam['base',])] %>% colnames
  for (name in est_names){
    if (par.lb %>% pull(name) %>% is.numeric()){
      par.guess.all[[name]] = runif(nguess, par.lb %>% pull(name), par.ub %>% pull(name))
    } else {
      par.guess.all[[name]] = par.lb %>% pull(name)
    }
  } # for name
  
  return = par.guess.all
  
} # end random_tuess


estimate = function(est.set, tabs, par.guess, print_level = 0){
  # est.set$opt_method= 'two-stage' or 'crs-only' or 'pif-grid'
  
  # setup
  est_names = est.set$model_fam['base', is.na(est.set$model_fam['base',])] %>% colnames
  
  # === optimization methods ===
  
  ## two-stage ====
  if (est.set$opt_method == 'two-stage'){
    
    # grab options
    if ('opts_list1' %in% names(est.set)){
      temp.opts.stage1 = est.set$opts_list1
    } else {
      temp.opts.stage1 = opts.crs()
    }
    temp.opts.stage1$print_level = print_level
    
    if ('opts_list2' %in% names(est.set)){
      temp.opts.stage2 = est.set$opts_list1
    } else {
      temp.opts.stage2 = opts.qa()
    }
    temp.opts.stage2$print_level = print_level    
    
    
    # likelihood that gets optimized
    minme = function(parvec){
      
      par = parvec2par(parvec,par.guess,est_names)
      return = negloglike(tabs, par)
      
    } # end minme
    
    # test likelihood at guess
    parvecguess = par2parvec(par.guess,est_names)
    minme(parvecguess)    
    
    # optimize
    opt = nloptr(
      x0 = par2parvec(par.guess, est_names)
      , lb = par2parvec(est.set$model_fam['lb',], est_names)
      , ub = par2parvec(est.set$model_fam['ub',], est_names)
      , eval_f = minme
      , opts = temp.opts.stage1
    )
    
    # optimize again 
    opt = nloptr(
      x0 = opt$solution
      , lb = par2parvec(est.set$model_fam['lb',], est_names)
      , ub = par2parvec(est.set$model_fam['ub',], est_names)
      , eval_f = minme
      , opts = temp.opts.stage2
    )
    
    
    # store
    parvec.hat = opt$solution
    parhat = parvec2par(parvec.hat,par.guess,est_names)
    
    
    ## crs-only ====    
  } else if (est.set$opt_method == 'crs-only'){
    
    # likelihood that gets optimized
    minme = function(parvec){
      
      par = parvec2par(parvec,par.guess,est_names)
      return = negloglike(tabs, par)
      
    } # end minme
    
    # test likelihood at guess
    parvecguess = par2parvec(par.guess,est_names)
    minme(parvecguess)    
    
    opt = nloptr(
      x0 = par2parvec(par.guess, est_names)
      , lb = par2parvec(est.set$model_fam['lb',], est_names)
      , ub = par2parvec(est.set$model_fam['ub',], est_names)
      , eval_f = minme
      , opts = opts.crs(print_level = print_level, maxeval = 1000)
    )        
    
    # store
    parvec.hat = opt$solution
    parhat = parvec2par(parvec.hat,par.guess,est_names)
    
    ### pif-grid ====
  } else if (est.set$opt_method == 'pif-grid'){
    
    # grab options
    if ('opt_pif_grid_base' %in% names(est.set)){
      pif_grid_base = est.set$opt_pif_grid_base
    } else {
      pif_grid_base = c(0,1,0.05)
    }
    if ('opts_list' %in% names(est.set)){
      temp.opts = est.set$opts_list
    } else {
      temp.opts = opts.qa(print_level = print_level)
    }
    temp.opts$print_level = print_level
    
    
    # create pif grid from bounds
    templb = est.set$model_fam[2,'pif']
    tempub = est.set$model_fam[3,'pif']
    pif_grid = c(pif_grid_base, templb, tempub) %>% sort()
    pif_grid = pif_grid[
      pif_grid >= templb & pif_grid <= tempub
    ]
    
    
    # loop
    parhat_grid = tibble()
    opt_grid = list()
    temp.par.guess = par.guess 
    for (i in 1:length(pif_grid)){
      
      # print(pif_grid[i])
      
      # fix pif
      temp.par.guess = temp.par.guess %>% mutate(pif = pif_grid[i])
      
      # remove pif from stuff that is optimized directly
      temp_est_names = est_names[!est_names == 'pif']
      
      # likelihood that gets optimized, with fixed pif
      minme = function(parvec){
        
        par = parvec2par(parvec,temp.par.guess,temp_est_names)
        return = negloglike(tabs, par)
        
      } # end minme      
      
      # test likelihood at guess
      parvecguess = par2parvec(temp.par.guess,temp_est_names)
      minme(parvecguess)    
      
      # optimize
      temp_opt = nloptr(
        x0 = par2parvec(temp.par.guess, temp_est_names)
        , lb = par2parvec(est.set$model_fam['lb',], temp_est_names)
        , ub = par2parvec(est.set$model_fam['ub',], temp_est_names)
        , eval_f = minme
        , opts = temp.opts
      )         
      
      # enforce bounds
      parvechat = pmax(temp_opt$solution, par2parvec(est.set$model_fam['lb',], temp_est_names)) 
      parvechat = pmin(parvechat, par2parvec(est.set$model_fam['ub',], temp_est_names))
      
      temp_parhat = parvec2par(parvechat,temp.par.guess,temp_est_names) %>% 
        mutate(negloglike = temp_opt$objective)
      
      # store
      parhat_grid = rbind(
        parhat_grid
        , temp_parhat
      )
      opt_grid[[i]] = temp_opt
      
      # update guess
      temp.par.guess = temp_parhat
      
    } # for i in pif_grid
    
    # find best
    ibest = which.min(parhat_grid$negloglike)
    parhat = parhat_grid[ibest,] %>% select(-negloglike)
    opt = opt_grid[[ibest]]
    
    # optimize again over all est_names
    
    # likelihood that gets optimized
    minme = function(parvec){
      par = parvec2par(parvec,par.guess,est_names)
      return = negloglike(tabs, par)
    } # end minme
    
    # test likelihood at guess
    parvecguess = par2parvec(par.guess,est_names)
    minme(parvecguess)    
    
    opt = nloptr(
      x0 = par2parvec(parhat, est_names)
      , lb = par2parvec(est.set$model_fam['lb',], est_names)
      , ub = par2parvec(est.set$model_fam['ub',], est_names)
      , eval_f = minme
      , opts = temp.opts
    )      
    
    # store
    parvec.hat = opt$solution
    parhat = parvec2par(parvec.hat,par.guess,est_names)
    
  } # end if est.set$opt_method
  
  
  return = list(
    opt = opt
    , par = parhat
    , negloglike = opt$objective
  )
  
} # end estimate

toc = Sys.time()

## Bootstrap ====

# bootstrap
bootstrap = function(tabs,set.boot,nboot,bootname = 'deleteme'){
  
  nsignal = length(tabs)
  
  # initialize bootstrap
  set.seed(1053)
  id_all = matrix(
    sample(1:nsignal, size = nsignal*nboot, replace = T), nrow = nsignal 
  )
  id_all[,1] = 1:nsignal # fix first column is always the empirical data
  
  # initialize random guesses
  par.guess.all = random_guess(set.boot, nboot, seed = 155)  
  
  start_time = Sys.time()
  bootpar = tibble()
  bootstat = tibble()
  for (booti in 1:nboot){
    print(paste0('booti = ', booti))    
    
    tic = Sys.time()
    
    # estimate one ===
    
    # if first round, use empirical data
    if (booti == 1){
      boottabs = tabs
      
      # if second round, do either simple boot  
    } else if (set.boot$boot_method == 'simple'){
      
      boottabs = tabs[id_all[,booti]]
      
      
      # or semi-parameteric boot
    } else if (set.boot$boot_method == 'semipar') {
      
      simcross = sim_pubcross(par = bootpar[1, ], nport = 5000, ndate = 350, seed = id_all[1, booti])
      tempn = min(dim(simcross)[1], length(tabs)) # use nport that is smaller of simcross or emp nport
      boottabs = simcross$tabs[1:tempn]
      
    } # end if booti == 1
    
    par.guess = par.guess.all[booti,]    
    est = estimate(set.boot,boottabs,par.guess)
    
    # make stats
    stat = make_stats(est$par)
    temp = make_stats_pub(boottabs, est$par)
    stat = stat %>% 
      mutate(
        bias_mean = mean(temp$bias)
        , bias_med = median(temp$bias)
        , fdrloc_mean = mean(temp$fdr_tabs)
        , fdrloc_med = median(temp$fdr_tabs)        
        , sbias_mean = mean(temp$bias_signed)
        , sbias_med = median(temp$bias_signed)        
      )
    
    
    # store === 
    bootpar = rbind(bootpar, est$par %>% mutate(booti = booti))
    bootstat = rbind(bootstat, stat  %>% mutate(booti = booti))
    
    
    if (booti%%10 == 0){
      filename = paste0('boot ', bootname, ' started ', start_time)
      filename = gsub('[:]', '-', filename)
      filename = substr(filename, 1, nchar(filename)-6)
      save(
        list = ls(all.names = T)
        , file = paste0('intermediate/', filename, '.Rdata')
        , envir = environment()
      )
    } # if booti
    
    # feedback
    toc = Sys.time()
    print(toc - tic)
    print(est$par)
    
    p.fit = hist_emp_vs_par(boottabs,est$par) +
      ggtitle(
        paste0('pif = ', round(est$par$pif,2), ', ', substr(est$opt$message, 7, 20))
      )
    p.pif = bootpar %>% 
      ggplot(aes(x=pif)) +
      geom_histogram(aes(y=stat(density)), breaks = seq(0,1,0.1)) +
      theme_minimal()    
    p.hurdle = bootstat %>% 
      ggplot(aes(x=h_fdr5)) + 
      geom_histogram(aes(y=stat(density)), breaks = seq(0,5,0.5)) + 
      theme_minimal() +
      xlab('t-hurdle 5%')
    p.bias = bootstat %>% 
      ggplot(aes(x=bias_mean)) + 
      geom_histogram(aes(y=stat(density)), breaks = seq(0,1,0.05)) + 
      theme_minimal() +
      xlab('mean bias')
    p.sbias = bootstat %>% 
      ggplot(aes(x=sbias_mean)) + 
      geom_histogram(aes(y=stat(density)), breaks = seq(0,1,0.05)) + 
      theme_minimal() +
      xlab('mean sbias')    
    p.fdrpub = bootstat %>% 
      ggplot(aes(x=fdrloc_mean)) + 
      geom_histogram(aes(y=stat(density)), breaks = seq(0,1,0.05)) + 
      theme_minimal() +
      xlab('fdr pub')    
    
    grid.arrange(p.fit,p.pif,p.hurdle,p.fdrpub, p.bias,p.sbias, nrow = 3)  
    
    
  } # for booti
  
  bootlist = list(par = bootpar, stat = bootstat)
  return = bootlist
  
} # end function bootstrap

# DATA GENERATION ====

import_cz = function(dl = F){
  
  # download if desired
  if (dl) {
    
    # download signal documentation 
    target_dribble = pathRelease %>% drive_ls() %>% 
      filter(name=='SignalDoc.csv')
    
    drive_download(target_dribble, path = 'intermediate/deleteme.csv', overwrite = T)
    
    signaldoc = fread('intermediate/deleteme.csv') %>% 
      mutate(
        signalname = Acronym
        , pubdate = as.Date(paste0(Year, '-12-31'))
        , sampend = as.Date(paste0(SampleEndYear, '-12-31'))
        , sampstart = as.Date(paste0(SampleStartYear, '-01-01'))
      ) %>% 
      arrange(signalname) %>% 
      select(signalname, pubdate, sampend, sampstart)
    
    # download returns and add samptype
    target_dribble = pathRelease %>% drive_ls() %>%
      filter(name=='Portfolios') %>% drive_ls() %>%
      filter(name=='Full Sets OP') %>% drive_ls() %>%
      filter(name=='PredictorLSretWide.csv')
    
    drive_download(target_dribble[1,], path = 'intermediate/deleteme2.csv', overwrite = T)
    
    ret = fread('intermediate/deleteme2.csv') %>% 
      pivot_longer(
        -date, names_to = 'signalname', values_to = 'ret'
      ) %>% 
      filter(!is.na(ret)) %>% 
      left_join(signaldoc, by = 'signalname') %>% 
      mutate(
        samptype = case_when(
          (date >= sampstart) & (date <= sampend) ~ 'in-samp'
          , (date > sampend) & (date <= pubdate) ~ 'out-of-samp'
          , (date > pubdate) ~ 'post-pub'
          , TRUE ~ NA_character_
        )
      ) %>% 
      select(-pubdate,-sampend,-sampstart)
    
    # summarize into pubcross
    pubcross = ret %>% 
      filter(samptype == 'in-samp') %>% 
      group_by(signalname) %>% 
      summarize(
        rbar = mean(ret), vol = sd(ret), ndate = dplyr::n()
      ) %>% 
      mutate(
        tstat = rbar/vol*sqrt(ndate), tabs = abs(tstat)
      )
    
    write.csv(x = pubcross, file = 'output/pubcross.csv', row.names = F)
    write.csv(x = ret, file = 'output/pubpanel.csv', row.names = F)
    
  } # if dl 
  
  pubcross = fread('output/pubcross.csv') # loads pubcross
  
  return = pubcross
  
} # end download_cz

sim_cz_residuals = function(nport, ndate, seed = NULL){
  
  set.seed(seed)
  
  # create matrix of returns ====
  ret = fread('output/pubpanel.csv') %>% 
    filter(
      !is.na(ret), year(date) >= 1963 # focus on more observed years
    ) 
  
  # demean
  meanret = ret %>% 
    group_by(signalname) %>% 
    summarize(rbar = mean(ret, na.rm=T))
  ret = ret %>% 
    left_join(meanret, by = 'signalname') %>% 
    mutate(resid = ret - rbar)
  
  # bootstrap signals
  signallist = ret$signalname %>% unique
  bootsignals = tibble(
    portid = 1:nport, signalname = sample(signallist, nport, replace = T)
  )
  
  # bootstrap dates
  datelist = ret$date %>%  unique()
  bootdates = sample(datelist, ndate, replace = T) %>% sort()
  
  # make bootstrapped dataset
  bootret = expand.grid(
    portid = 1:nport, date = bootdates
  ) %>% 
    as_tibble() %>% 
    left_join(bootsignals, by = c('portid') ) %>% 
    left_join(ret, by = c('signalname','date') )
  
  return = bootret
  
} # end sim_cz_residuals



# function for simulating truth, but only cross-sectional
sim_pubcross = function(par, nport, ndate, eptype = 'ar1', rho = 0.5, seed = NULL){
  set.seed(seed)
  
  # cross-section of mus  
  type = runif(nport) > par$pif
  
  #   start with everything false
  mu = numeric(nport) 
  
  #   then add true stuff
  mutrue.o = make_mutrue_object(par)
  mu[type] = r(mutrue.o)(sum(type))
  
  # cross-section of observables
  #   ar1 noise
  if (eptype  == 'ar1'){
    ep = numeric(nport)
    for (i in 1:nport){
      if (i==1){
        ep[i] = rnorm(1,0,1)
      } else{
        ep[i] = rho*ep[i-1] + sqrt(1-rho^2)*rnorm(1,0,1)
      }
    } # end for i
    
    #   bootstrapped from cz
  } else if (eptype == 'boot'){
    
    bootret = sim_cz_residuals(nport = nport, ndate = ndate)
    
    bootsum = bootret %>% 
      group_by(portid) %>% 
      filter(!is.na(resid)) %>% 
      summarize(
        m = mean(resid), vol = sd(resid), n = sum(!is.na(resid)), ep = m/vol*sqrt(n)
      )  
    
    ep = bootsum$ep
    
    # check correlations
    # tempmat = bootret %>% 
    #   filter(portid <= 200) %>% 
    #   group_by(portid) %>% 
    #   arrange(portid,date) %>% 
    #   mutate(dateid = row_number()) %>% 
    #   select(dateid, portid, resid) %>% 
    #   pivot_wider(names_from = portid, values_from = resid) %>% 
    #   select(-dateid) 
    # 
    # tempcor = cor(tempmat, use = 'pairwise.complete.obs')
    # hist(tempcor)
    
    
  } # if eptype
  
  # cross-sectional stats
  cross = tibble(
    type = type
    , mu = mu
    , tstat = mu + ep
    , tabs = abs(tstat)
    , pub = runif(nport) <= pub_prob(tabs, par)
  )
  
  # published subset
  pubcross = cross %>% filter(pub)
  
  return = pubcross
  
} # end sim_pubcross

sim_hlz = function(nsim, nfac, pif, mua, pubpar1, rho){
  # this is its own function bc the exchangeable
  # correlations are awkward (even though they don't really
  # change anything)
  
  # cross-section of mus  
  type = runif(nsim*nfac) > pif
  
  #   start with everything false
  mu = numeric(nsim*nfac) 
  
  #   then add true stuff
  mu[type] = rexp(sum(type), rate = 1/mua)
  
  # simulate ep
  # use model u = c + v, Cov(c,v) = 0, so rho = sigv^2/(sigu*sigv) = sigv/sigu
  # so sigv = sigu*rho = 1*rho 
  # and then sigu^2 = sigc^2 + sigv^2 so sigc = sqrt(1-rho^2)
  sigv = rho
  sigc = sqrt(1-rho^2)
  
  temp = list(nsim)
  for (simi in 1:nsim){
    
    v = rnorm(nfac, 0, sigv)
    c = rnorm(nfac, 0, sigc)
    u = c+v
    
    temp[[simi]] = data.frame(
      simi = simi, ep = u
    )
    
  } # for simi
  noise_dat = do.call(rbind, temp)  %>% as_tibble() 
  
  temp = runif(nsim*nfac)
  hlz_dat = noise_dat %>% 
    mutate(
      type = type
      , mu = mu, t = mu + ep, tabs = abs(t)
      , pub = (tabs > 2.57) | (tabs > 1.96 & tabs < 2.57 & temp <= pubpar1 )
    )
  
  return = hlz_dat
  
} # end sim_hlz

# STATS AND PLOTTING ====

shrink_one = function(tabs,par,mutrue.o=NULL, take_abs = T){
  
  if (is.null(mutrue.o)){
    mutrue.o = make_mutrue_object(par)
  }
  
  # kernel_pos is f(tabs,mu|mutrue).  it gets reused in the numerator and denom.
  # note kernel_pos depends on observed tabs
  # the dirac term disappears in the numerator, and becomes dnorm(tabs,0,1) in the denom
  # take_abs = T => additional dnorm terms
  kernel_pos = function(mm) (dnorm(tabs,mm,1)+take_abs*dnorm(-tabs,mm,1))*d(mutrue.o)(mm) 
  numer = (1-par$pif)*integrate(function (mm) mm*kernel_pos(mm), -Inf, Inf)$value
  denom = (1+take_abs)*par$pif*dnorm(tabs,0,1) + (1-par$pif)*integrate(kernel_pos, -Inf, Inf)$value
  
  Emu_tabs = numer/denom 
  Eep_tabs = tabs-Emu_tabs
  bias = Eep_tabs/tabs
  
  return = tibble(Emu_tabs = Emu_tabs, Eep_tabs = Eep_tabs, bias = bias)
  
} # end shrink_one


make_stats_pub = function(t, par){
  
  mutrue.o = make_mutrue_object(par)
  tabs.o = make_tabs_object(par)
  
  bias = numeric(length(t))
  fdr_tabs = numeric(length(t))
  bias_signed = numeric(length(t))
  for (i in 1:length(t)){
    bias[i] = shrink_one(t[i], par, mutrue.o, take_abs = T)$bias
    bias_signed[i] = shrink_one(t[i], par, mutrue.o, take_abs = F)$bias
    fdr_tabs[i] = 2*dnorm(t[i])*par$pif/d(tabs.o)(t[i]) # times 2 bc abs val
  } # for i
  
  
  stats_pub = list(bias = bias, fdr_tabs = fdr_tabs, bias_signed = bias_signed)
  
  return = stats_pub
  
} # make_stats_pub


make_stats = function(par, fdrlist = c(10, 5, 1)){
  
  # latent density 
  tabs.o = make_tabs_object(par)
  
  # find fdr
  tlist = c(1.96, seq(0,5,0.1)) %>% sort()
  pr_discovery = NA*numeric(length(tlist))
  fdr = NA*numeric(length(tlist))
  for (ti in 1:length(tlist)){
    tt = tlist[ti]
    temp1 = integrate(function(ttt) d(tabs.o)(ttt), tt, 10)$value
    temp2 = integrate(function(ttt) d(tabs.o)(ttt), 10, Inf)$value
    pr_discovery[ti] = temp1 + temp2
    fdr[ti] = 2*pnorm(-tt)/pr_discovery[ti]*par$pif*100
    if (fdr[ti]<=min(fdrlist) & tlist[ti] >= 1.96){
      break
    }
  } # end for ti
  
  # find hurdles
  hurdlelist = NA*numeric(length(fdrlist))
  for (fdri in 1:length(fdrlist)){
    hurdlelist[fdri] = min(tlist[which(fdr<fdrlist[fdri])])
  }
  
  # save into stat
  stat = tibble(
    stat = paste0('h_fdr', fdrlist), value = hurdlelist
  ) %>%
    rbind(
      tibble(
        stat = 'pr_tgt_2', value = pr_discovery[tlist == 1.96]
      )
    ) %>% 
    pivot_wider(names_from = stat, values_from = value)
  
  return = stat
  
} # end make_stats


hist_emp_vs_par = function(tabs,par){
  
  # latent density
  tabs.o = make_tabs_object(par)
  
  # observed density
  #   find denominator for conditional density
  #   integrate in pieces for accuracy
  denom1 = integrate(function (tt) d(tabs.o)(tt)*pub_prob(tt,par), 0, 1.96)$value
  denom2 = integrate(function (tt) d(tabs.o)(tt)*pub_prob(tt,par), 1.96, 6)$value  
  denom3 = integrate(function (tt) d(tabs.o)(tt)*pub_prob(tt,par), 6, Inf)$value  
  denom = denom1 + denom2 + denom3
  
  f_obs = function(tt) d(tabs.o)(tt)*pub_prob(tt,par)/denom
  
  # setup
  edge = c(seq(0,15,1), 40, 100)
  freq = numeric(length(edge)-1)
  for (i in 2:length(edge)){
    a = edge[i-1]; b = edge[i]
    freq[i-1] = integrate(f_obs, a, b)$value
  }
  
  hemp = hist(tabs, edge, plot = F)
  
  plotme = tibble(
    group = 'emp', tabs = hemp$mids, freq = hemp$density
  ) %>% rbind(
    tibble(
      group = 'mod', tabs = hemp$mids, freq = freq    
    )
  ) %>% 
    filter(!is.na(freq))
  
  xmin = 0; xmax = 15
  plotme %>% 
    filter(tabs >= xmin, tabs <= xmax) %>% 
    ggplot(aes(x=tabs,y=freq,group=group,color=group)) +
    geom_line() +
    geom_point(size=3) +
    theme_minimal() +
    scale_color_manual(values = c('blue','red')) + 
    theme(legend.position = c(0.8,0.8)) +
    xlim(xmin,xmax)  
  
} # end hist_emp_vs_par




# FIGURES ====

library(latex2exp)
library(extrafont)



MATBLUE = rgb(0,0.4470,0.7410)
MATRED = rgb(0.8500, 0.3250, 0.0980)
MATYELLOW = rgb(0.9290, 0.6940, 0.1250)

NICEBLUE = "#619CFF"
NICEGREEN = "#00BA38"
NICERED = "#F8766D"

chen_theme =   theme_classic() +
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
