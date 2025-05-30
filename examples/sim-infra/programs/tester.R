library(bcf)
library(tidyverse)
library(glue)
library(pROC)
library(mvtnorm)
library(modi)
library(corpcor)
library(readxl)
library(rstan)

if (running_local) {
  source('programs/jennifer-dgp-local.R')  
  source('programs/test_acic.R')  
} else {
  source('/home/bcf-sim-infra/programs/jennifer-dgp.R')
  source('/home/bcf-sim-infra/programs/test_acic.R')  
}

#I need test_bcf and test_stan to have identical inputs, 
#so I don't get arguments for the other function passed to the DGP as part of ...
test_bcf <- function(...,
                     seed=02139, log=NULL, treedir=NULL,
                     ibcf=TRUE,
                     ubcf=FALSE,
                     obcf=FALSE,
                     oldbcf=FALSE,
                     isstan=FALSE,
                     w_as_covar=FALSE,
                     dgp_use_res = TRUE,
                     trt_eff_scenario = 'het',
                     hardcode_sigma_u     = FALSE,
                     hardcode_sigma_u_val = 0,
                     hardcode_sigma_v     = FALSE,
                     hardcode_sigma_v_val = 0,
                     hardcode_rho         = FALSE,
                     hardcode_rho_val     = 0,
                     block_b0_b1 = TRUE,
                     block_batch_size=100, 
                     #force_pos_b1=FALSE,
                     #include_tau0=FALSE,
                     #taux_stz=FALSE,
                     #tau0_sd=NULL, tau0_mean=0,
                     simplified_return=TRUE,
                     ret_fit=FALSE,
                     sd_moderate = NULL,
                     ate_prior_sd=NULL,
                     base_moderate=.25,
                     savedir='tests/results', savename=NULL, s3_save=FALSE, 
                     just_run = FALSE,
                     nburn=1000, nsim=1000, nthin=1,
                     n_chains=4, n_cores=n_chains,
                     model_name='linear_stan_model',
                     sigu_hyperprior = NULL,
                     sigu_pcthyperprior = NULL,
                     control=NULL) {
  #Sanitize input - AWS passes pars as strings
  seed <- as.numeric(seed)
  dgp_use_res <- as.logical(dgp_use_res)
  ibcf <- as.logical(ibcf)
  ubcf <- as.logical(ubcf)
  obcf <- as.logical(obcf)
  oldbcf <- as.logical(oldbcf)
  isstan <- as.logical(isstan)
  block_batch_size <- as.numeric(block_batch_size)
  simplified_return <- as.logical(simplified_return)
  ret_fit <- as.logical(ret_fit)
  s3_save <- as.logical(s3_save)
  nburn <- as.numeric(nburn)
  nsim <- as.numeric(nsim)
  nthin <- as.numeric(nthin)
  n_chains <- as.numeric(n_chains)
  n_cores <- as.numeric(n_cores)
  ate_prior_sd <- as.numeric(ate_prior_sd)
  
  set.seed(seed)
  if (dgp_use_res) {
    data <- genData(trt_eff_scenario=trt_eff_scenario,
                    ...)
  } else {
    data <- genData(sig_u=0,
                    sig_v=0,
                    trt_eff_scenario=trt_eff_scenario,
                    ...)
  }
  
  if (ubcf & obcf) {
    stop('Can\'t do both u and o flavors of ibcf')
  }
  if (!ibcf & (ubcf | obcf)) {
    stop('Can\'t do u or i with wbcf')
  }
  if (ibcf & oldbcf) {
    stop('Can\'t do ibcf with oldbcf')
  }
  if (isstan) {
    stop('I am not stan')
  }
  
  if (w_as_covar) {
    xmat <- as.matrix(data[,c('x1','x2','x3','x4','x5','w')])
  } else {
    xmat <- as.matrix(data[,c('x1','x2','x3','x4','x5')])  
  }
  
  pihat <- data$pi
  data$sigma <- data$sigma_y
  
  if (is.null(treedir)) {
    treedir <- tempdir()
  }
  if (is.null(log)) {
    log <- tempfile(fileext='.log')
  }
  if (n_cores==1) {
    sink(log, type='output')
  }
  
  #BCF default is SD(y) so we'll use that
  if (is.null(sd_moderate)) {
    sd_moderate <- sqrt(Hmisc::wtd.var(data$y, data$w))
  }
  
  model_name <- case_when(ibcf & ubcf ~ 'uBCF',
                          ibcf & obcf ~ 'oBCF',
                          ibcf ~ 'iBCF',
                          oldbcf ~ 'oldBCF',
                          TRUE ~ 'wBCF')
  weights_touse <- data$w
  
  if (model_name=='uBCF') {
    hardcode_sigma_v <- TRUE
    hardcode_sigma_v_val <- 0
    hardcode_rho <- TRUE
    hardcode_rho_val <- 0
  } else if (model_name=='oBCF') {
    hardcode_sigma_u <- TRUE
    hardcode_sigma_u_val <- data$sigma_u[1]
    hardcode_sigma_v <- TRUE
    hardcode_sigma_v_val <- data$sigma_v[1]
    hardcode_rho <- TRUE
    hardcode_rho_val <- data$rho[1]
  } else if (model_name=='oldBCF') {
    weights_touse <- rep(1, nrow(data))
  }
  
  if (is.null(sigu_hyperprior) & !is.null(sigu_pcthyperprior)) {
    #internally if sigu hyperprior is specified it gets divided by sdy then by .674
    sigu_hyperprior <- sigu_pcthyperprior * sqrt(Hmisc::wtd.var(data$y, weights_touse)) * .674
  } else if (!is.null(sigu_hyperprior)) {
    sigu_pcthyperprior <- sigu_hyperprior/sqrt(Hmisc::wtd.var(data$y, weights_touse))
  } else {
    sigu_pcthyperprior <- 2/3
  }
  
  #Reset seed so fit is separately reproducible
  set.seed(seed)
  tictoc::tic()
  fit <- bcf(y          = data$y,
             z          = data$z,
             x_control  = xmat,
             x_moderate = xmat,
             pihat      = pihat,
             w          = weights_touse,
             ate_prior_sd  = ate_prior_sd,
             sd_moderate   = sd_moderate,
             base_moderate = base_moderate,
             n_chains      = n_chains,
             n_cores       = n_cores,
             n_threads     = 1,
             nburn         = nburn,
             nsim          = nsim,
             nthin         = nthin,
             include_random_effects = ibcf,
             hardcode_sigma_u       = hardcode_sigma_u,
             hardcode_sigma_u_val   = hardcode_sigma_u_val,
             hardcode_sigma_v       = hardcode_sigma_v,
             hardcode_sigma_v_val   = hardcode_sigma_v_val,
             hardcode_rho           = hardcode_rho,
             hardcode_rho_val       = hardcode_rho_val,
             block_v_rho            = !hardcode_rho & !hardcode_sigma_v,
             block_batch_size       = block_batch_size,
             block_b0_b1            = block_b0_b1,
             sigu_hyperprior        = sigu_hyperprior,
             #force_pos_b1=force_pos_b1,
             #include_tau0=include_tau0,
             #taux_stz=taux_stz,
             #tau0_sd=tau0_sd,
             #tau0_mean=tau0_mean,
             save_tree_directory = treedir,
             log_file            = log,
             simplified_return   = simplified_return,
             verbose             = 0)
  timing <- tictoc::toc()
  
  if (n_cores==1) {
    sink(NULL)
  }
  
  stanchains <- stanify(fit$raw_chains,data$w)
  mixing <- rstan::monitor(stanchains,warmup=0,print=FALSE,probs=c(.025,.05,.1,.5,.9,.95,.975)) %>%
    as_tibble(rownames='par') %>%
    mutate(real = case_when(par=='sigma'      ~ data$sigma[1],
                            par=='sigma_y'    ~ data$sigma_y[1],
                            par=='sigma_u'    ~ data$sigma_u[1],
                            par=='sigma_v'    ~ data$sigma_v[1],
                            par=='sigma2'     ~ data$sigma[1]^2,
                            par=='sigma2_y'   ~ data$sigma_y[1]^2,
                            par=='sigma2_u'   ~ data$sigma_u[1]^2,
                            par=='sigma2_v'   ~ data$sigma_v[1]^2,
                            par=='rho'        ~ data$rho[1],
                            par=='sigv_delta' ~ data$sigma_v[1]^2 + 2*data$rho[1]*data$sigma_u[1]*data$sigma_v[1],
                            par=='tau_bar'    ~ weighted.mean(data$tau,data$w),
                            par=='mu_bar'     ~ weighted.mean(data$mu,data$w),
                            par=='yhat_bar'   ~ weighted.mean(data$y,data$w)),
           cover80 = real>=`10%` & real<=`90%`,
           cover90 = real>=`5%` & real<=`95%`,
           cover95 = real>=`2.5%` & real<=`97.5%`,
           width80 = `90%`-`10%`,
           width90 = `95%`-`5%`,
           width95 = `97.5%`-`2.5%`) %>%
    select(par, real, mean, sd,
           cover80, cover90, cover95, 
           width80, width90, width95, 
           p025=`2.5%`, p5 = `5%`, p10=`10%`, p50=`50%`,
           p90=`90%`, p95=`95%`, p975=`97.5%`,
           Rhat, n_eff, Bulk_ESS, Tail_ESS)
  
  if(just_run) {
    ret <- list(data=data,mixing=mixing)
    if (ret_fit) {
      ret$fit <- fit
    }
    return(ret)
  }
  
  acceptance <- fit$acceptance
  
  scalecorr <- get_corr(stanchains,c('b0','b1','tau_scale','mu_scale','delta_mu','real_mu_scale'))
  
  taus <- lapply(fit$raw_chains,`[[`,'tau') %>% do.call(what=rbind)
  
  indiv <- summarize_indiv(z=data$z, id=data$id, real=data$tau, ests=taus, name='tau', draws=nsim, chains=n_chains) %>%
    mutate(taux=data$taux,
           u=data$u,
           v=data$v,
           w=data$w)
  
  yhats <- lapply(fit$raw_chains,`[[`,'yhat') %>% do.call(what=rbind)
  
  exemplar <- list(tau  = exemplar_ptiles(data=data, real=data$tau, ests=taus,                                                       name='tau'),
                   yhat = exemplar_ptiles(data=data, real=data$y,   ests=yhats, name='yhat'))
  
  exemplar_summy <- list(tau  = exemplarsummy(exemplar=exemplar$tau,  name='tau'),
                         yhat = exemplarsummy(exemplar=exemplar$yhat, name='yhat'))
  
  real_resid <- data$u + data$z*data$v
  
  if (ibcf) {
    vs <- lapply(fit$raw_chains,`[[`,'v') %>% do.call(what=rbind)
    exemplar$v <- exemplar_ptiles(data=data, real=data$v, ests=vs, name='v')
    exemplar_summy$v <- exemplarsummy(exemplar=exemplar$v, name='v')
    
    resid_hat <- lapply(fit$raw_chains,`[[`,'u') %>% do.call(what=rbind) + sweep(vs, STATS = data$z, MARGIN = 2, FUN='*')
  } else {
    resid_hat <- -1 * sweep(yhats, STATS = data$y, MARGIN=2, FUN='-')
  }
  
  exemplar$resid <- exemplar_ptiles(data=data, real=real_resid, ests=resid_hat, name='resid')
  exemplar_summy$resid <- exemplarsummy(exemplar=exemplar$resid, name='resid')
  
  #Also add indiv resid
  indiv_resid <- summarize_indiv(z=data$z, id=data$id, real=real_resid, ests=resid_hat, name='resid', draws=nsim, chains=n_chains) %>%
    mutate(taux=data$taux,
           u=data$u,
           v=data$v,
           w=data$w)
  
  zidx <- data$z==1
  SATT <- weighted.mean(data$tau[zidx], data$w[zidx])
  SATT_draws <- apply(taus[,zidx],1,weighted.mean, data$w[zidx])
  SATT_hat <- mean(SATT_draws)
  lb80 <- quantile(SATT_draws,.1)
  ub80 <- quantile(SATT_draws,.9)
  lb90 <- quantile(SATT_draws,.05)
  ub90 <- quantile(SATT_draws,.95)
  lb95 <- quantile(SATT_draws,.025)
  ub95 <- quantile(SATT_draws,.975)
  
  zsmallidx <- data$z==1 & data$w <= quantile(data$w,.25)
  SATTsmall <- weighted.mean(data$tau[zsmallidx], data$w[zsmallidx])
  SATTsmall_draws <- apply(taus[,zsmallidx],1,weighted.mean, data$w[zsmallidx])
  SATTsmall_hat <- mean(SATTsmall_draws)
  lb80small <- quantile(SATTsmall_draws,.1)
  ub80small <- quantile(SATTsmall_draws,.9)
  lb90small <- quantile(SATTsmall_draws,.05)
  ub90small <- quantile(SATTsmall_draws,.95)
  lb95small <- quantile(SATTsmall_draws,.025)
  ub95small <- quantile(SATTsmall_draws,.975)
  
  zmid1idx <- data$z==1 & data$w > quantile(data$w,.25) & data$w <= quantile(data$w,.5)
  SATTmid1 <- weighted.mean(data$tau[zmid1idx], data$w[zmid1idx])
  SATTmid1_draws <- apply(taus[,zmid1idx],1,weighted.mean, data$w[zmid1idx])
  SATTmid1_hat <- mean(SATTmid1_draws)
  lb80mid1 <- quantile(SATTmid1_draws,.1)
  ub80mid1 <- quantile(SATTmid1_draws,.9)
  lb90mid1 <- quantile(SATTmid1_draws,.05)
  ub90mid1 <- quantile(SATTmid1_draws,.95)
  lb95mid1 <- quantile(SATTmid1_draws,.025)
  ub95mid1 <- quantile(SATTmid1_draws,.975)
  
  zmid2idx <- data$z==1 & data$w > quantile(data$w,.5) & data$w <= quantile(data$w,.75)
  SATTmid2 <- weighted.mean(data$tau[zmid2idx], data$w[zmid2idx])
  SATTmid2_draws <- apply(taus[,zmid2idx],1,weighted.mean, data$w[zmid2idx])
  SATTmid2_hat <- mean(SATTmid2_draws)
  lb80mid2 <- quantile(SATTmid2_draws,.1)
  ub80mid2 <- quantile(SATTmid2_draws,.9)
  lb90mid2 <- quantile(SATTmid2_draws,.05)
  ub90mid2 <- quantile(SATTmid2_draws,.95)
  lb95mid2 <- quantile(SATTmid2_draws,.025)
  ub95mid2 <- quantile(SATTmid2_draws,.975)
  
  zlargeidx <- data$z==1 & data$w > quantile(data$w,.75)
  SATTlarge <- weighted.mean(data$tau[zlargeidx], data$w[zlargeidx])
  SATTlarge_draws <- apply(taus[,zlargeidx],1,weighted.mean, data$w[zlargeidx])
  SATTlarge_hat <- mean(SATTlarge_draws)
  lb80large <- quantile(SATTlarge_draws,.1)
  ub80large <- quantile(SATTlarge_draws,.9)
  lb90large <- quantile(SATTlarge_draws,.05)
  ub90large <- quantile(SATTlarge_draws,.95)
  lb95large <- quantile(SATTlarge_draws,.025)
  ub95large <- quantile(SATTlarge_draws,.975)
  
  overall <- tibble(SATT=SATT,
                    SATT_hat=SATT_hat,
                    SATTcover80 = SATT>=lb80 & SATT<=ub80,
                    SATTcover90 = SATT>=lb90 & SATT<=ub90,
                    SATTcover95 = SATT>=lb95 & SATT<=ub95,
                    SATTwidth80 = ub80-lb80,
                    SATTwidth90 = ub90-lb90,
                    SATTwidth95 = ub95-lb95,
                    SATTbias = SATT_hat - SATT,
                    PEHE = sqrt(mean((indiv$mean - indiv$truth)^2)),
                    PEHT = sqrt(mean((indiv$mean[indiv$z==1] - indiv$truth[indiv$z==1])^2)),
                    resid_RMSE = sqrt(mean((indiv_resid$mean - indiv_resid$truth)^2)),
                    residT_RMSET = sqrt(mean((indiv_resid$mean[indiv_resid$z==1] - indiv_resid$truth[indiv_resid$z==1])^2)),
                    tau_var = mean(apply(taus,1,Hmisc::wtd.var,weights=data$w)),
                    tau_var_t = mean(apply(taus[,data$z==1],1,Hmisc::wtd.var,weights=data$w[data$z==1])),
                    CATEcover80 = mean(indiv$cover80),
                    CATEcover90 = mean(indiv$cover90),
                    CATEcover95 = mean(indiv$cover95),
                    CATTcover80 = mean(indiv$cover80[indiv$z==1]),
                    CATTcover90 = mean(indiv$cover90[indiv$z==1]),
                    CATTcover95 = mean(indiv$cover95[indiv$z==1]),
                    residcover80 = mean(indiv_resid$cover80),
                    residcover90 = mean(indiv_resid$cover90),
                    residcover95 = mean(indiv_resid$cover95),
                    residTcover80 = mean(indiv_resid$cover80[indiv_resid$z==1]),
                    residTcover90 = mean(indiv_resid$cover90[indiv_resid$z==1]),
                    residTcover95 = mean(indiv_resid$cover95[indiv_resid$z==1]),
                    SATTp025=lb95,
                    SATTp5=lb90,
                    SATTp10=lb80,
                    SATTp90=ub80,
                    SATTp95=ub90,
                    SATTp975=ub95,
                    sdy = sqrt(Hmisc::wtd.var(data$y, data$w)),
                    tau_bar = weighted.mean(data$tau,data$w),
                    sd_tau = sqrt(Hmisc::wtd.var(data$tau, data$w)),
                    sd_taux = sqrt(Hmisc::wtd.var(data$taux, data$w)),
                    sd_mu = sqrt(Hmisc::wtd.var(data$mu, data$w)),
                    sd_mux = sqrt(Hmisc::wtd.var(data$mux, data$w)),
                    timing = timing$toc - timing$tic) %>%
    #Also grab CATE/CATT widths from indiv output
    bind_cols(indiv %>%
                summarize(CATEwidth80 = mean(width80),
                          CATEwidth90 = mean(width90),
                          CATEwidth95 = mean(width95),
                          CATTwidth80 = mean(ifelse(z==1,width80,NA_real_), na.rm=TRUE),
                          CATTwidth90 = mean(ifelse(z==1,width90,NA_real_), na.rm=TRUE),
                          CATTwidth95 = mean(ifelse(z==1,width95,NA_real_), na.rm=TRUE),
                          CATUwidth80 = mean(ifelse(z==0,width80,NA_real_), na.rm=TRUE),
                          CATUwidth90 = mean(ifelse(z==0,width90,NA_real_), na.rm=TRUE),
                          CATUwidth95 = mean(ifelse(z==0,width95,NA_real_), na.rm=TRUE),
                          neff_lt10   = mean(n_eff < 10),
                          neff_lt100  = mean(n_eff < 100),
                          neff_q10    = quantile(n_eff, .1),
                          neff_q25    = quantile(n_eff, .25),
                          neff_q50    = quantile(n_eff, .5))) %>%
    bind_cols(indiv_resid %>%
                summarize(residwidth80 = mean(width80),
                          residwidth90 = mean(width90),
                          residwidth95 = mean(width95),
                          residTwidth80 = mean(ifelse(z==1,width80,NA_real_), na.rm=TRUE),
                          residTwidth90 = mean(ifelse(z==1,width90,NA_real_), na.rm=TRUE),
                          residTwidth95 = mean(ifelse(z==1,width95,NA_real_), na.rm=TRUE))) %>%
    bind_cols(indiv %>%
                filter(z==1 & w <= quantile(data$w,.25)) %>%
                summarize(SATTsmall = weighted.mean(truth,w),
                          SATTsmall_hat = weighted.mean(mean,w),
                          CATEsmallRMSE = sqrt(mean((truth-mean)^2)),
                          CATEsmallcover80 = mean(cover80),
                          CATEsmallcover90 = mean(cover90),
                          CATEsmallcover95 = mean(cover95),
                          CATEsmallwidth80 = mean(width80),
                          CATEsmallwidth90 = mean(width90),
                          CATEsmallwidth95 = mean(width95)) %>%
                mutate(SATTsmallcover80 = SATTsmall>=lb80small & SATTsmall<=ub80small,
                       SATTsmallcover90 = SATTsmall>=lb90small & SATTsmall<=ub90small,
                       SATTsmallcover95 = SATTsmall>=lb95small & SATTsmall<=ub95small,
                       SATTsmallwidth80 = ub80small-lb80small,
                       SATTsmallwidth90 = ub90small-lb90small,
                       SATTsmallwidth95 = ub95small-lb95small))%>%
    bind_cols(indiv %>%
                filter(z==1 & w > quantile(data$w, .25) & w <= quantile(data$w,.5)) %>%
                summarize(SATTmid1 = weighted.mean(truth,w),
                          SATTmid1_hat = weighted.mean(mean,w),
                          CATEmid1RMSE = sqrt(mean((truth-mean)^2)),
                          CATEmid1cover80 = mean(cover80),
                          CATEmid1cover90 = mean(cover90),
                          CATEmid1cover95 = mean(cover95),
                          CATEmid1width80 = mean(width80),
                          CATEmid1width90 = mean(width90),
                          CATEmid1width95 = mean(width95)) %>%
                mutate(SATTmid1cover80 = SATTmid1>=lb80mid1 & SATTmid1<=ub80mid1,
                       SATTmid1cover90 = SATTmid1>=lb90mid1 & SATTmid1<=ub90mid1,
                       SATTmid1cover95 = SATTmid1>=lb95mid1 & SATTmid1<=ub95mid1,
                       SATTmid1width80 = ub80mid1-lb80mid1,
                       SATTmid1width90 = ub90mid1-lb90mid1,
                       SATTmid1width95 = ub95mid1-lb95mid1))%>%
    bind_cols(indiv %>%
                filter(z==1 & w > quantile(data$w, .5) & w <= quantile(data$w,.75)) %>%
                summarize(SATTmid2 = weighted.mean(truth,w),
                          SATTmid2_hat = weighted.mean(mean,w),
                          CATEmid2RMSE = sqrt(mean((truth-mean)^2)),
                          CATEmid2cover80 = mean(cover80),
                          CATEmid2cover90 = mean(cover90),
                          CATEmid2cover95 = mean(cover95),
                          CATEmid2width80 = mean(width80),
                          CATEmid2width90 = mean(width90),
                          CATEmid2width95 = mean(width95)) %>%
                mutate(SATTmid2cover80 = SATTmid2>=lb80mid2 & SATTmid2<=ub80mid2,
                       SATTmid2cover90 = SATTmid2>=lb90mid2 & SATTmid2<=ub90mid2,
                       SATTmid2cover95 = SATTmid2>=lb95mid2 & SATTmid2<=ub95mid2,
                       SATTmid2width80 = ub80mid2-lb80mid2,
                       SATTmid2width90 = ub90mid2-lb90mid2,
                       SATTmid2width95 = ub95mid2-lb95mid2))%>%
    bind_cols(indiv %>%
                filter(z==1 & w > quantile(data$w,.75)) %>%
                summarize(SATTlarge = weighted.mean(truth,w),
                          SATTlarge_hat = weighted.mean(mean,w),
                          CATElargeRMSE = sqrt(mean((truth-mean)^2)),
                          CATElargecover80 = mean(cover80),
                          CATElargecover90 = mean(cover90),
                          CATElargecover95 = mean(cover95),
                          CATElargewidth80 = mean(width80),
                          CATElargewidth90 = mean(width90),
                          CATElargewidth95 = mean(width95)) %>%
                mutate(SATTlargecover80 = SATTlarge>=lb80large & SATTlarge<=ub80large,
                       SATTlargecover90 = SATTlarge>=lb90large & SATTlarge<=ub90large,
                       SATTlargecover95 = SATTlarge>=lb95large & SATTlarge<=ub95large,
                       SATTlargewidth80 = ub80large-lb80large,
                       SATTlargewidth90 = ub90large-lb90large,
                       SATTlargewidth95 = ub95large-lb95large))
  
  inputs = tibble(trt_eff_scenario=trt_eff_scenario,
                  ...,
                  n_chains=n_chains,
                  nburn=nburn,
                  nsim=nsim,
                  nthin=nthin,
                  ibcf=ibcf,
                  ubcf=ubcf,
                  obcf=obcf,
                  model=model_name,
                  seed=seed,
                  w_as_covar=w_as_covar,
                  sigu_hyperprior=sigu_pcthyperprior)
  
  calib <- calibcheck(data, fit$raw_chains)
  
  if (ibcf) {
    tauxs <- taus - vs
    #uv coverage
    indiv_uv <- bind_rows(summarize_indiv(z=data$z, id=data$id, real=data$u, ests=lapply(fit$raw_chains,`[[`,'u') %>% do.call(what=rbind), name='u', draws=nsim, chains=n_chains),
                          summarize_indiv(z=data$z, id=data$id, real=data$v, ests=vs, name='v', draws=nsim, chains=n_chains),
                          summarize_indiv(z=data$z, id=data$id, real=data$taux, ests=tauxs, name='taux', draws=nsim, chains=n_chains))
    
    #Also grab posterior diagnostics
    postcorr <- get_corr(stanchains,c('sigma_y','sigma_u','sigma_v','rho','sigv_delta'))
    
    ret <- list(overall=overall,
                inputs=inputs,
                mixing=mixing,
                indiv=indiv,
                indiv_resid=indiv_resid,
                indiv_uv=indiv_uv,
                exemplar=exemplar,
                exemplar_summy=exemplar_summy,
                postcorr=postcorr,
                scalecorr=scalecorr,
                acceptance=acceptance,
                calib=calib,
                data=data)
  } else {
    ret <- list(overall=overall,
                inputs=inputs,
                mixing=mixing,
                indiv=indiv,
                indiv_resid=indiv_resid,
                exemplar=exemplar,
                exemplar_summy=exemplar_summy,
                scalecorr=scalecorr,
                acceptance=acceptance,
                calib=calib,
                data=data)
  }
  
  if (ret_fit) {
    ret$fit <- fit
  }
  
  fit <- NULL
  gc()
  
  if (!is.null(savename)) {
    if (s3_save) {
      td <- tempdir()
      saveRDS(ret,glue('{td}/return.RDS'))
      system(glue('aws s3 cp {td}/return.RDS s3://{savedir}/{savename}.RDS'))
    } else {
      saveRDS(ret,paste0(savedir,'/',savename,'.RDS'))  
    }
  }
  return(ret)
}

compile_model <- function(model, model_folder=ifelse(running_local,'programs','/home/bcf-sim-infra/programs')) {
  #Recompile if needed
  if (!file.exists(sprintf('%s/%s.RDS', model_folder, model))) {
    comp_model <- stan_model(sprintf('%s/%s.stan', model_folder, model))
    saveRDS(comp_model, sprintf('%s/%s.RDS', model_folder, model))
  } else {
    comp_model <- readRDS(sprintf('%s/%s.RDS', model_folder, model))
    if (!identical(comp_model@model_code[1],  sprintf('%s/%s.stan', model_folder, model) %>% readLines %>% paste(collapse='\n'))) {
      comp_model <- stan_model(sprintf('%s/%s.stan', model_folder, model))
      saveRDS(comp_model, sprintf('%s/%s.RDS', model_folder, model))
    }
  }
  return(comp_model)
}

test_stan <- function(...,
                      seed=02139, log=NULL, treedir=NULL,
                      ibcf=TRUE,
                      ubcf=FALSE,
                      obcf=FALSE,
                      isstan=TRUE,
                      dgp_use_res = TRUE,
                      trt_eff_scenario = 'het',
                      hardcode_sigma_u     = FALSE,
                      hardcode_sigma_u_val = 0,
                      hardcode_sigma_v     = FALSE,
                      hardcode_sigma_v_val = 0,
                      hardcode_rho         = FALSE,
                      hardcode_rho_val     = 0,
                      block_b0_b1 = TRUE,
                      block_batch_size=100,
                      simplified_return=TRUE,
                      ret_fit=FALSE,
                      sd_moderate=NULL,
                      base_moderate=NULL,
                      just_run=FALSE,
                      ate_prior_sd=NULL,
                      savedir='tests/results', savename=NULL, s3_save=FALSE, 
                      nburn=1000, nsim=1000,nthin=1,
                      n_chains=4, n_cores=n_chains,
                      model_name='linear_stan_model',
                      sigu_hyperprior = NULL,
                      control=NULL) {
  #Sanitize input - AWS passes pars as strings
  seed <- as.numeric(seed)
  isstan <- as.logical(isstan)
  dgp_use_res <- as.logical(dgp_use_res)
  ret_fit <- as.logical(ret_fit)
  s3_save <- as.logical(s3_save)
  nburn <- as.numeric(nburn)
  nsim <- as.numeric(nsim)
  n_chains <- as.numeric(n_chains)
  n_cores <- as.numeric(n_cores)
  sigu_hyperprior <- as.numeric(sigu_hyperprior)
  ate_prior_sd <- as.numeric(ate_prior_sd)
  
  if (!isstan) {
    stop('But I am stan!')
  }
  
  set.seed(seed)
  if (dgp_use_res) {
    data <- genData(trt_eff_scenario=trt_eff_scenario,
                    ...)
  } else {
    data <- genData(sig_u=0,
                    sig_v=0,
                    trt_eff_scenario=trt_eff_scenario,
                    ...)
  }
  
  model <- compile_model(model_name)
  
  data$sigma <- data$sigma_y
  
  #standardize the Xs to be exact so we don't have to do it in stan
  xmat <- as.matrix(data[,c('x1','x2','x3','x4','x5')])
  xmat <- xmat %>%
    sweep(MARGIN=2, STATS=apply({.}, 2, mean), FUN='-') %>%
    sweep(MARGIN=2, STATS=apply({.}, 2, sd), FUN='/')
  
  sdy <- sqrt(modi::weighted.var(data$y, data$w))
  
  if (is.null(sigu_hyperprior) || identical(sigu_hyperprior, numeric(0))) {
    sigu_hyperprior <- 2*sdy/3
  }
  
  stan_data <- list(N=nrow(data),
                    sdy = sdy,
                    y=(data$y-weighted.mean(data$y)) / sdy,
                    z=data$z,
                    w=data$w/mean(data$w),
                    sigu_hyperprior=sigu_hyperprior/sdy,
                    ate_prior_sd=ate_prior_sd/sdy,
                    X=xmat,
                    K=ncol(xmat))
  
  #Reset seed so fit is separately reproducible
  tictoc::tic()
  fit <- sampling(model,
                  stan_data,
                  pars=c('sigma_y', 'sigma_u','tau_bar','u','rho','sigma_v'),
                  seed=seed,
                  iter = nsim+nburn,
                  warmup=nburn,
                  chains=n_chains,
                  cores=n_cores,
                  control=control,
                  refresh=0)
  timing <- tictoc::toc()
  
  stanchains <-  extract(fit, c('sigma_y', 'sigma_u','tau_bar','rho','sigma_v'), permuted=FALSE)
  mixing <- rstan::monitor(stanchains,warmup=0,print=FALSE,probs=c(.025,.05,.1,.5,.9,.95,.975)) %>%
    as_tibble(rownames='par') %>%
    mutate(real = case_when(par=='sigma'      ~ data$sigma[1],
                            par=='sigma_y'    ~ data$sigma_y[1],
                            par=='sigma_u'    ~ data$sigma_u[1],
                            par=='sigma_v'    ~ data$sigma_v[1],
                            par=='sigma2'     ~ data$sigma[1]^2,
                            par=='sigma2_y'   ~ data$sigma_y[1]^2,
                            par=='sigma2_u'   ~ data$sigma_u[1]^2,
                            par=='sigma2_v'   ~ data$sigma_v[1]^2,
                            par=='rho'        ~ data$rho[1],
                            par=='sigv_delta' ~ data$sigma_v[1]^2 + 2*data$rho[1]*data$sigma_u[1]*data$sigma_v[1],
                            par=='tau_bar'    ~ weighted.mean(data$tau,data$w),
                            par=='mu_bar'     ~ weighted.mean(data$mu,data$w),
                            par=='yhat_bar'   ~ weighted.mean(data$y,data$w)),
           cover80 = real>=`10%` & real<=`90%`,
           cover90 = real>=`5%` & real<=`95%`,
           cover95 = real>=`2.5%` & real<=`97.5%`,
           width80 = `90%`-`10%`,
           width90 = `95%`-`5%`,
           width95 = `97.5%`-`2.5%`) %>%
    select(par, real, mean, sd,
           cover80, cover90, cover95, 
           width80, width90, width95, 
           p025=`2.5%`, p5 = `5%`, p10=`10%`, p50=`50%`,
           p90=`90%`, p95=`95%`, p975=`97.5%`,
           Rhat, n_eff, Bulk_ESS, Tail_ESS)
  
  real_resid <- data$u + data$z*data$v
  est_resid <- extract(fit,'u')$u
  exemplar <- list(resid  = exemplar_ptiles(data=data, real=real_resid, ests=est_resid, name='resid'))
  
  exemplar_summy <- list(resid  = exemplarsummy(exemplar=exemplar$resid, name='resid'))
  
  indiv_resid <- summarize_indiv(z=data$z, id=data$id, real=real_resid, ests=est_resid, name='resid', draws=nsim, chains=n_chains) %>%
    mutate(taux=data$taux,
           u=data$u,
           v=data$v,
           w=data$w)
  
  SATT <- weighted.mean(data$tau[data$z==1], data$w[data$z==1])
  SATT_draws <- extract(fit,'tau_bar')$tau_bar
  SATT_hat <- mean(SATT_draws)
  lb80 <- quantile(SATT_draws,.1)
  ub80 <- quantile(SATT_draws,.9)
  lb90 <- quantile(SATT_draws,.05)
  ub90 <- quantile(SATT_draws,.95)
  lb95 <- quantile(SATT_draws,.025)
  ub95 <- quantile(SATT_draws,.975)
  
  overall <- tibble(SATT=SATT,
                    SATT_hat=SATT_hat,
                    SATTcover80 = SATT>=lb80 & SATT<=ub80,
                    SATTcover90 = SATT>=lb90 & SATT<=ub90,
                    SATTcover95 = SATT>=lb95 & SATT<=ub95,
                    SATTwidth80 = ub80-lb80,
                    SATTwidth90 = ub90-lb90,
                    SATTwidth95 = ub95-lb95,
                    SATTbias = SATT_hat - SATT,
                    SATTp025=lb95,
                    SATTp5=lb90,
                    SATTp10=lb80,
                    SATTp90=ub80,
                    SATTp95=ub90,
                    SATTp975=ub95,
                    resid_RMSE = sqrt(mean((indiv_resid$mean - indiv_resid$truth)^2)),
                    residT_RMSET = sqrt(mean((indiv_resid$mean[indiv_resid$z==1] - indiv_resid$truth[indiv_resid$z==1])^2)),
                    residcover80 = mean(indiv_resid$cover80),
                    residcover90 = mean(indiv_resid$cover90),
                    residcover95 = mean(indiv_resid$cover95),
                    residTcover80 = mean(indiv_resid$cover80[indiv_resid$z==1]),
                    residTcover90 = mean(indiv_resid$cover90[indiv_resid$z==1]),
                    residTcover95 = mean(indiv_resid$cover95[indiv_resid$z==1]),
                    sdy = sdy,
                    timing = timing$toc - timing$tic) %>%
    bind_cols(indiv_resid %>%
                summarize(residwidth80 = mean(width80),
                          residwidth90 = mean(width90),
                          residwidth95 = mean(width95),
                          residTwidth80 = mean(ifelse(z==1,width80,NA_real_), na.rm=TRUE),
                          residTwidth90 = mean(ifelse(z==1,width90,NA_real_), na.rm=TRUE),
                          residTwidth95 = mean(ifelse(z==1,width95,NA_real_), na.rm=TRUE)))
  
  inputs = tibble(trt_eff_scenario=trt_eff_scenario,
                  ...,
                  n_chains=n_chains,
                  nburn=nburn,
                  nsim=nsim,
                  model='stan',
                  seed=seed)
  
  ret <- list(overall=overall,
              inputs=inputs,
              mixing=mixing,
              indiv_resid=indiv_resid,
              exemplar=exemplar,
              exemplar_summy=exemplar_summy,
              data=data)
  
  if (ret_fit) {
    ret$fit <- fit
  }
  
  fit <- NULL
  gc()
  
  if (!is.null(savename)) {
    if (s3_save) {
      td <- tempdir()
      saveRDS(ret,glue('{td}/return.RDS'))
      system(glue('aws s3 cp {td}/return.RDS s3://{savedir}/{savename}.RDS'))
    } else {
      saveRDS(ret,paste0(savedir,'/',savename,'.RDS'))  
    }
  }
  return(ret)
}

stanify <- function(chains, w, pars=NULL, derived_pars=NULL) {
  ibcf <- 'sigma_y' %in% names(chains[[1]])
  
  if (is.null(pars)) {
    pars <- c('sigma', 'sigma_y','sigma_u','sigma_v','rho','tau0','tau_bar','mu_bar','yhat_bar','tau_scale','mu_scale','b0','b1','delta_mu')
  }
  pars <- intersect(pars, names(chains[[1]]))
  
  if (is.null(derived_pars)) {
    if (ibcf) {
      derived_pars <- c('tau_het', 'mu_het', 'real_mu_scale','sigv_delta', 'sigma2_y', 'sigma2_u', 'sigma2_v')
    } else {
      derived_pars <- c('tau_het', 'mu_het', 'real_mu_scale', 'sigma2')
    }
  }
  
  pars <- c(pars, derived_pars)
  
  npar = length(pars)
  nchain = length(chains)
  ndraw = nrow(chains[[1]]$mu)
  
  stanobj <- array(0,dim=c(ndraw, nchain, npar))
  dimnames(stanobj)[[3]] <- pars
  for(i in 1:nchain) {
    j <- 0
    for (par in pars) {
      j <- j+1
      if (par=='real_mu_scale') {
        stanobj[,i,j] <- chains[[i]]$mu_scale / sqrt(chains[[i]]$delta_mu)
      } else if (par %in% c('sigv_delta','trt_var_delta')) {
        stanobj[,i,j] <- chains[[i]]$sigma_v^2 + 2*chains[[i]]$rho*chains[[i]]$sigma_v*chains[[i]]$sigma_u
        #NB: heterogeneity here is inclusive of u/v. Maybe want a taux_het and mux_het for ibcf?
      } else if (str_detect(par,'het')) {
        opar <- str_remove(par,'_het')
        stanobj[,i,j] <- sqrt(apply(chains[[i]][[opar]], 1, modi::weighted.var, w))
      } else if (str_detect(par,'sigma2')) {
        opar <- str_replace(par,'sigma2','sigma')
        stanobj[,i,j] <- chains[[i]][[opar]]^2
      } else if (!str_detect(par,'bar')) {
        stanobj[,i,j] <- chains[[i]][[par]]
      } else {
        opar <- str_remove(par,'_bar')
        stanobj[,i,j] <- matrixStats::rowWeightedMeans(chains[[i]][[opar]], w)
      }
    }
  }
  return(stanobj)
}

summarize_indiv <- function(z, id, real, ests, name, draws, chains) {
  tibble(id=id,
         z=z,
         par=name,
         truth=real,
         mean = apply(ests,2,mean),
         sd = apply(ests,2,sd),
         #Data comes in from rbinding the chain values together. So rows = chains * iters, cols=obs
         #For stan ess we need to turn that into iters*chains*obs
         n_eff = apply(array(ests, dim=c(draws, chains, length(id))),3,rstan:::ess_rfun),
         pr_pos = apply(ests, 2, function(x) mean(x>0)),
         pr_gt_ate = apply(ests, 2, function(x) mean(x>mean(real))),
         lb80 = apply(ests,2,quantile,.1),
         ub80 = apply(ests,2,quantile,.9),
         lb90 = apply(ests,2,quantile,.05),
         ub90 = apply(ests,2,quantile,.95),
         lb95 = apply(ests,2,quantile,.025),
         ub95 = apply(ests,2,quantile,.975)) %>%
    mutate(cover80 = truth>=lb80 & truth<=ub80,
           cover90 = truth>=lb90 & truth<=ub90,
           cover95 = truth>=lb95 & truth<=ub95,
           width80 = ub80-lb80,
           width90 = ub90-lb90,
           width95 = ub95-lb95) %>%
    rename(p025=lb95,
           p5=lb90,
           p10=lb80,
           p90=ub80,
           p95=ub90,
           p975=ub95) %>%
    select(id, z, par, truth, mean, sd, n_eff, pr_pos, pr_gt_ate, matches('cover'), matches('width'), p025, p5, p10, p90, p95, p975)
}

exemplar_ptiles <- function(data, real, ests, name) {
  idx <- data$z==1
  est_ptiles <- t(apply(ests[,idx], 1, ptile))
  
  #We're going to do exemplars overall, but also by size (smallest quartile vs others)
  sizecut <- quantile(data$w[idx], .25)
  
  tibble(id                   = data$id[idx],
         z                    = data$z[idx],
         size_cat             = ifelse(data$w[idx] < sizecut, 'small','not_small'),
         xxx                  = real[idx],
         xxx_pos              = as.numeric(xxx > 0),
         xxx_gt_5             = as.numeric(xxx > 5),
         xxx_gt_10            = as.numeric(xxx > 10),
         xxx_gt_20            = as.numeric(xxx > 20),
         xxx_gt_30            = as.numeric(xxx > 30),
         xxx_gt_45            = as.numeric(xxx > 45),
         xxx_gt_60            = as.numeric(xxx > 60),
         xxx_gt_ate           = as.numeric(xxx > mean(real)),
         xxx_ptile            = ptile(xxx),
         xxx_top10pct         = as.numeric(xxx_ptile >  .9),
         xxx_top20pct         = as.numeric(xxx_ptile >  .8),
         xxx_top25pct         = as.numeric(xxx_ptile >  .75),
         xxx_top50pct         = as.numeric(xxx_ptile >  .5),
         xxx_mid50pct         = as.numeric(xxx_ptile >  .25 & xxx_ptile<=.75),
         xxx_bot50pct         = as.numeric(xxx_ptile <= .5),
         xxx_bot25pct         = as.numeric(xxx_ptile <= .25),
         xxx_bot20pct         = as.numeric(xxx_ptile <= .2),
         xxx_bot10pct         = as.numeric(xxx_ptile <= .1),
         xxx_hat              = apply(ests[,idx],2,mean),
         xxx_hat_prob_pos     = colMeans(ests[,idx] > 0),
         xxx_hat_prob_gt_5    = colMeans(ests[,idx] > 5),
         xxx_hat_prob_gt_10   = colMeans(ests[,idx] > 10),
         xxx_hat_prob_gt_20   = colMeans(ests[,idx] > 20),
         xxx_hat_prob_gt_30   = colMeans(ests[,idx] > 30),
         xxx_hat_prob_gt_45   = colMeans(ests[,idx] > 45),
         xxx_hat_prob_gt_60   = colMeans(ests[,idx] > 60),
         xxx_hat_prob_gt_ate  = colMeans(ests[,idx] > mean(real)),
         xxx_hat_ptile        = colMeans(est_ptiles),
         xxx_hat_probtop10pct = colMeans(est_ptiles >  .9),
         xxx_hat_probtop20pct = colMeans(est_ptiles >  .8),
         xxx_hat_probtop25pct = colMeans(est_ptiles >  .75),
         xxx_hat_probtop50pct = colMeans(est_ptiles >  .5),
         xxx_hat_probmid50pct = colMeans(est_ptiles >  .25 & est_ptiles <= .75),
         xxx_hat_probbot50pct = colMeans(est_ptiles <= .5),
         xxx_hat_probbot25pct = colMeans(est_ptiles <= .25),
         xxx_hat_probbot20pct = colMeans(est_ptiles <= .2),
         xxx_hat_probbot10pct = colMeans(est_ptiles <= .1)) %>%
    arrange(id) %>%
    rename_all(str_replace,'xxx',name)
}

gen_rank <- function(data, sortvar, ntop) {
  ranks <- nrow(data) + 1 - rank(data[[sortvar]], ties.method='first')
  ranks[ranks>ntop] <- NA_real_
  return(ranks)
}
# JS's utility function to return percentiles of a variable.
ptile = function(x){
  temp_ranks = rank(x, ties.method='first')
  ptile = cut(temp_ranks, quantile(temp_ranks, probs=0:100/100), include.lowest=TRUE, labels=(0:100/100)[-1])
  #Ptile is a factor, so as.numeric pulls out 1:100; as.character uses labels, so then we convert those
  return(as.numeric(as.character(ptile)))
}

exemplarsummy_ptile = function(exemplar, name, cut, size) {
  #input requires truth (true indicator), prob (estimated prob of indicator) and exemplar (exemplar indicator)
  
  auc = tibble(par=name,
               cutoff = cut, 
               size=size,
               auc = safe_auc(exemplar$truth, exemplar$prob))
  
  rmse <- exemplar %>%
    filter(exemplar==1) %>%
    mutate(diff = xxx - xxx_hat,
           diff_ptile = xxx_ptile - xxx_hat_ptile) %>%
    summarize(par=name,
              cutoff = cut,
              size=size,
              group='Exemplars',
              n=n(),
              avg_prob = mean(prob),
              share_small = mean(size_cat=='small'),
              rmse = sqrt(mean(diff^2)),
              maqe = mean(abs(diff_ptile)))
  
  rmse_misclass <- exemplar %>%
    filter(exemplar==1 & truth==0) %>%
    mutate(diff = xxx - xxx_hat,
           diff_ptile = xxx_ptile - xxx_hat_ptile) %>%
    summarize(par=name,
              cutoff = cut,
              size=size,
              group='Misclassified exemplars',
              n=n(),
              avg_prob = mean(prob),
              share_small = mean(size_cat=='small'),
              rmse = sqrt(mean(diff^2)),
              maqe = mean(abs(diff_ptile)))
  
  confus <- exemplar %>%
    filter(exemplar==1) %>%
    select(xxx_pos, xxx_gt_5, xxx_gt_10, xxx_gt_20, xxx_gt_30, xxx_gt_45, xxx_gt_60,
           xxx_top10pct, xxx_top20pct, xxx_top25pct, xxx_top50pct, 
           xxx_mid50pct, 
           xxx_bot50pct, xxx_bot25pct, xxx_bot20pct, xxx_bot10pct) %>%
    summarize_all(mean) %>%
    rename_all(str_remove_all,'xxx_|pct') %>%
    mutate(par=name,
           cutoff=cut,
           size=size) %>%
    select(par, cutoff, size, everything())
  
  return(list(auc=auc,
              rmse = bind_rows(rmse, rmse_misclass),
              confus=confus))
}

ex_cut_wrapper <- function(exname, size, exemplar, name) {
  if (size!='all') {
    exemplar <- filter(exemplar, size_cat==size)
  }
  
  #How many practices do we want to flag? 10%, but # depends on how many Ts there actually are
  #New: assume we're doing the top 5% only
  ntop <- round(.05*nrow(exemplar))
  
  if (exname=='pos' | str_detect(exname,'^gt_[0-9+]')) {
    exemplar$truth    <- exemplar[[paste0('xxx_', exname)]]
    exemplar$prob     <- exemplar[[paste0('xxx_hat_prob_', exname)]]
    cut <- exname
  } else {
    exemplar$truth    <- exemplar[[paste0('xxx_', exname, 'pct')]]
    exemplar$prob     <- exemplar[[paste0('xxx_hat_prob', exname, 'pct')]]
    cut <- paste0(exname,'pct')
  }
  
  #recalculate exemplars among the size sample only
  exemplar$exemplar <- as.numeric(!is.na(gen_rank(exemplar, 'prob', ntop)))
  
  exemplarsummy_ptile(exemplar, name=name, cut=cut, size=size)
}

exemplarsummy = function(exemplar, name){
  #First genericize the names
  colnames(exemplar) <- colnames(exemplar) %>% 
    str_replace('^(v|tau|yhat|resid)$','xxx') %>%
    str_replace('^(v|tau|yhat|resid)_','xxx_') %>%
    str_replace('exemplar_(v|tau|yhat|resid)_','exemplar_xxx_')
  
  ex_cuts <- expand_grid(exname=c('pos', 'gt_5','gt_10','gt_20','gt_30','gt_45','gt_60',
                                  'top10','top20','top25','top50',
                                  'mid50',
                                  'bot50','bot25','bot20','bot10'),
                         size=c('all','small','not_small'))
  
  ex_results <- purrr::pmap(ex_cuts, ex_cut_wrapper, exemplar, name)
  
  #Keep the pr pos stuff separate
  return(list(auc = lapply(ex_results, `[[`,'auc') %>% bind_rows(),
              rmse = lapply(ex_results, `[[`,'rmse') %>% bind_rows(),
              confus = lapply(ex_results, `[[`,'confus') %>% bind_rows()))
}

#AUC errors out if there's only one level of the truth (e.g. taus when homogenous and sig_v=0, or vs when sig_v=0)
#So write a version that returns NA if that's the case
safe_auc <- function(truth, pred) {
  if (length(unique(truth))==1) {
    return(NA_real_)
  } else {
    return(as.numeric(auc(truth, pred)))
  }
}

get_corr <- function(stanchains, pars){
  corr <- pars %>% lapply(function(x) as.vector(stanchains[,,x])) %>%
    do.call(what=cbind) %>%
    cor
  
  rownames(corr) <- colnames(corr) <- pars
  return(corr)
}

get_sunscreen <- function(data, chains) {
  u <- chains %>% lapply(`[[`,'u') %>% do.call(what=rbind)
  v <- chains %>% lapply(`[[`,'u') %>% do.call(what=rbind)
  tibble(id    = 1:nrow(data),
         treat = data$z,
         realu = data$u,
         realv = data$v,
         meanu = apply(u,2,mean),
         meanv = apply(v,2,mean),
         sdu   = apply(u,2,sd),
         sdv   = apply(v,2,sd),
         p10u  = apply(u,2,quantile,.1),
         p10v  = apply(v,2,quantile,.1),
         p90u  = apply(u,2,quantile,.9),
         p90v  = apply(v,2,quantile,.9)) %>%
    pivot_longer(cols=matches('u|v'),names_pattern = '(.+)([uv])',names_to = c('.value','par')) %>%
    mutate(cover80 = real >= p10 & real<=p90)
}

calibcheck = function(data, chains){
  # Calculates posterior mean estimates for the following.
  # - r2a:     The r-squared for residuals for treated observations.
  #              r2a = 1 - RSS / SST where
  #              RSS = sum_{i=1}^{n}((ui + vi + ei)^2), and
  #              SST = sum_{i=1}^{n}((y_i-ybar)^2)
  # -r2uv:     The ratio of treatment-related idiosyncratic variance to total idiosyncratic variance.
  #              r2uv = sig2v / (sig2v + sig2u)
  # -r2trt:    The proportion of treatment effect variance explained by covariates.
  #              r2trt = 1 - sig2v / (sig2v + var(taux))
  # -vartaux:  The variance of taux only, without v, for treated obs.
  # -vartau:   The variance of (taux + v), for treated obs.
  
  ibcf = 'u' %in% names(chains[[1]])
  
  isT <- data$z==1
  
  tau <- lapply(chains,`[[`,'tau') %>% do.call(what=rbind)
  mu  <- lapply(chains,`[[`,'mu') %>% do.call(what=rbind)
  if (ibcf) {
    taux <- tau - lapply(chains,`[[`,'v') %>% do.call(what=rbind)
    mux <- mu - lapply(chains,`[[`,'u') %>% do.call(what=rbind)
  } else {
    taux <- tau
    mux <- mu
  }
  
  #Save matrix of x-fitted values for treatment
  yxhat <- mux
  yxhat[,isT] <- mux[,isT] + taux[,isT]
  
  #Want y-yhat, so take -(yhat-y). Doesn't matter since we square anyways
  resids <- -sweep(yxhat,2,data$y,'-')
  var_tot <- weighted.var(data$y[isT],data$w[isT])
  var_resid_draws <- apply(resids[,isT], 1, weighted.var, data$w[isT])
  r2a_draws <- 1 - var_resid_draws / var_tot
  r2a_hat = mean(r2a_draws)
  
  # R2uv -------------------------------------------------------------------------------------------
  
  if(ibcf){
    # Posterior draws of sigma_u and sigma_v.
    sigu_draws = lapply(chains,`[[`,'sigma_u') %>% unlist
    sigv_draws = lapply(chains,`[[`,'sigma_v') %>% unlist
    
    r2uv_draws = sigv_draws^2 / (sigu_draws^2 + sigv_draws^2)
    r2uv_hat = mean(r2uv_draws)
  } else{
    r2uv_hat = NA
  }
  
  # R2trt ------------------------------------------------------------------------------------------
  
  if(ibcf){
    # Posterior draws of var(taux) and var(tau).
    var_taux_draws = apply(taux[,isT], 1, weighted.var, data$w[isT])
    var_tau_draws  = apply(tau[,isT],  1, weighted.var, data$w[isT])
    
    r2trt_draws = var_taux_draws / var_tau_draws
    r2trt_hat = mean(r2trt_draws)
  } else{
    r2trt_hat = 1
  }
  
  if(ibcf){
    # Posterior draws of var(taux) and var(tau).
    var_taux = mean(var_taux_draws)
    var_tau  = mean(var_tau_draws)
    
    yhat     = lapply(chains,`[[`,'yhat') %>% do.call(what=rbind)
    var_yhat = yhat %>% apply(1,weighted.var, data$w) %>% mean
    var_yhat_trt = yhat[,isT] %>% apply(1,weighted.var, data$w[isT]) %>% mean
  } else{
    var_tau  = tau[,isT] %>% apply(1, weighted.var, data$w[isT]) %>% mean
    var_taux = var_tau
    #yhatx already defined, and for vanilla BCF that is yhat
    var_yhat = yxhat %>% apply(1,weighted.var, data$w) %>% mean
    var_yhat_trt = yxhat[,isT] %>% apply(1,weighted.var, data$w[isT]) %>% mean
  }
  
  #-------------------------------------------------------------------------------------------------
  # Assemble and return dataframe of calibration metric estimates.
  out = tibble(r2a_input    = data$r2a_input[1],
               r2a_check    = data$r2a_check[1],
               r2a_hat      = r2a_hat,
               r2uv_input   = data$r2uv_input[1],
               r2uv_check   = data$r2uv_check[1],
               r2uv_hat     = r2uv_hat,
               r2trt_input  = data$r2trt_input[1],
               r2trt_check  = data$r2trt_check[1],
               r2trt_hat    = r2trt_hat,
               var_tau_trt  = var_tau,
               var_taux_trt = var_taux,
               var_yhat     = var_yhat,
               var_yhat_trt = var_yhat_trt)
  return(out)
}

make_jan_2023_control_file <- function() {
  if (running_local) {
    seeds <- c(readLines('data/10k_random_ints_20211202.txt'),
               readLines('data/10k_random_ints_20211203.txt'))
  } else {
    seeds <- c(readLines('/home/bcf-sim-infra/data/10k_random_ints_20211202.txt'),
               readLines('/home/bcf-sim-infra/data/10k_random_ints_20211203.txt'))
  }
  
  seeds <- seeds %>%
    as.numeric %>%
    unique
  
  base_su <- 1
  base_sv <- 1
  base_rho <- 0
  base_nt = 1000
  
  coarse_su <- c(.5,1,2)
  coarse_sv <- c(0,.5,1,2)
  coarse_rho <- c(-.5,0,.5)
  coarse_nt = c(500,1000)
  
  fine_su <- c(.25,.5,2/3,1,1.5,2,4)
  fine_sv <- c(0,.25,.5,2/3,1,1.5,2,4)
  fine_rho <- seq(-.75,.75,.25)
  fine_nt = c(250,500,1000,1500,2000)
  
  scenarios <- bind_rows(expand_grid(sigu_multiplier = coarse_su,
                                     sigv_multiplier = coarse_sv,
                                     rho             = coarse_rho,
                                     nT              = coarse_nt),
                         expand_grid(sigu_multiplier = fine_su,
                                     sigv_multiplier = base_sv,
                                     rho             = base_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = fine_sv,
                                     rho             = base_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = base_sv,
                                     rho             = fine_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = base_sv,
                                     rho             = base_rho,
                                     nT              = fine_nt)) %>%
    distinct() %>%
    mutate(world = 'medicare',
           ate_prior_sd = 20/sqrt(2/pi),
           nC = nT*2,
           scenario = glue('{world}_sigu{sigu_multiplier}_sigv{sigv_multiplier}_rho{str_replace(rho,"-","m")}_nT{nT}'))
  
  stopifnot(nrow(scenarios) == n_distinct(scenarios$scenario))
  
  #scenarios <- filter(scenarios, sigu_multiplier==1 & sigv_multiplier==1 & rho==0)
  
  test_list <- expand_grid(scenarios,
                           sim=1:50) %>%
    mutate(nburn   = 1000,
           nsim    = 1000,
           nthin   = 1,
           n_cores = 1,
           s3_path = 'bcf-sim-study/jan-2023-runs/1kburn-1ksim')
  
  #Add seeds before we duplicate for ibcf vs wbcf
  stopifnot(length(seeds) >= nrow(test_list))
  test_list$seed <- seeds[1:nrow(test_list)]
  
  test_list <- test_list %>%
    expand_grid(ibcf=c(TRUE,FALSE)) %>%
    mutate(fname = case_when(ibcf ~ paste0(scenario, '_ibcf_jan_2023_sim', str_pad(sim,width = 3, side='left',pad = '0')),
                             !ibcf ~ paste0(scenario, '_wbcf_jan_2023_sim', str_pad(sim,width = 3, side='left',pad = '0')))) %>%
    select(-scenario, -sim)
}

make_feb_2023_control_file <- function() {
  if (running_local) {
    seeds <- c(readLines('data/10k_random_ints_20211202.txt'),
               readLines('data/10k_random_ints_20211203.txt'))
  } else {
    seeds <- c(readLines('/home/bcf-sim-infra/data/10k_random_ints_20211202.txt'),
               readLines('/home/bcf-sim-infra/data/10k_random_ints_20211203.txt'))
  }
  
  seeds <- seeds %>%
    as.numeric %>%
    unique
  
  base_su <- 1
  base_sv <- 1
  base_rho <- 0
  base_nt = 1000
  
  coarse_su <- c(.5,1,2)
  coarse_sv <- c(0,.5,1,2)
  coarse_rho <- c(-.5,0,.5)
  coarse_nt = c(500,1000)
  
  fine_su <- c(.25,.5,2/3,1,1.5,2,4)
  fine_sv <- c(0,.25,.5,2/3,1,1.5,2,4)
  fine_rho <- seq(-.75,.75,.25)
  fine_nt = c(250,500,1000,1500,2000)
  
  scenarios <- bind_rows(expand_grid(sigu_multiplier = coarse_su,
                                     sigv_multiplier = coarse_sv,
                                     rho             = coarse_rho,
                                     nT              = coarse_nt),
                         expand_grid(sigu_multiplier = fine_su,
                                     sigv_multiplier = base_sv,
                                     rho             = base_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = fine_sv,
                                     rho             = base_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = base_sv,
                                     rho             = fine_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = base_sv,
                                     rho             = base_rho,
                                     nT              = fine_nt)) %>%
    distinct() %>%
    mutate(world = 'medicare',
           ate_prior_sd = 20/sqrt(2/pi),
           nC = nT*2,
           scenario = glue('{world}_sigu{sigu_multiplier}_sigv{sigv_multiplier}_rho{str_replace(rho,"-","m")}_nT{nT}'))
  
  stopifnot(nrow(scenarios) == n_distinct(scenarios$scenario))
  
  #scenarios <- filter(scenarios, sigu_multiplier==1 & sigv_multiplier==1 & rho==0)
  
  test_list <- expand_grid(scenarios,
                           sim=1:200) %>%
    mutate(nburn   = 1000,
           nsim    = 1000,
           nthin   = 1,
           n_cores = 1,
           s3_path = 'bcf-sim-study/feb-2023-runs/1kburn-1ksim')
  
  #Add seeds before we duplicate for ibcf vs wbcf
  stopifnot(length(seeds) >= nrow(test_list))
  test_list$seed <- seeds[1:nrow(test_list)]
  
  test_list <- test_list %>%
    expand_grid(model=c('ibcf','ubcf','obcf','wbcf')) %>%
    mutate(ibcf = model!='wbcf',
           ubcf = model=='ubcf',
           obcf = model=='obcf',
           fname = glue('{scenario}_{model}_feb_2023_sim{str_pad(sim,width = 3, side="left",pad = "0")}'),
           core = sigu_multiplier %in% coarse_su & sigv_multiplier %in% coarse_sv & rho %in% coarse_rho & nT %in% coarse_nt,
           base = sigu_multiplier == base_su & sigv_multiplier == base_sv & rho == base_rho & nT == base_nt) %>%
    #Don't run u/o bcf outside of the inner grid, and don't run more than 50 sims for the outer grid
    filter((model %in% c('ibcf','wbcf') & core) | (model %in% c('ibcf','wbcf') & sim <= 50) | (model %in% c('ubcf','obcf') & base) | (model %in% c('ubcf','obcf') & core & sim<=50)) %>%
    select(-scenario, -sim, -model, -core, -base)
}

make_small_feb_2023_control_file <- function() {
  make_feb_2023_control_file() %>%
    filter(sigu_multiplier == 1,
           sigv_multiplier == 1,
           rho             == 0,
           nT              == 1000)
}

first_50_feb_2023_control_file <- function() {
  make_feb_2023_control_file() %>%
    filter(sigu_multiplier %in% c(.5,1,2),
           sigv_multiplier %in% c(0,.5,1,2),
           rho             %in% c(-.5,0,.5),
           nT              %in% c(500,1000),
           as.numeric(str_match(fname,'sim([0-9]+)')[,2]) <= 50)
}

make_mar_2023_small_uiw_control_file <- function() {
  if (running_local) {
    seeds <- c(readLines('data/10k_random_ints_20211202.txt'),
               readLines('data/10k_random_ints_20211203.txt'))
  } else {
    seeds <- c(readLines('/home/bcf-sim-infra/data/10k_random_ints_20211202.txt'),
               readLines('/home/bcf-sim-infra/data/10k_random_ints_20211203.txt'))
  }
  
  seeds <- seeds %>%
    as.numeric %>%
    unique
  
  base_su <- 1
  base_sv <- 1
  base_rho <- 0
  base_nt = 1000
  
  coarse_su <- base_su
  coarse_sv <- base_sv
  coarse_rho <- base_rho
  coarse_nt = base_nt
  
  fine_su <- c(0,.25,.5,2/3,1,1.5,2,4)
  fine_sv <- c(0,.25,.5,2/3,1,1.5,2,4)
  fine_rho <- base_rho
  fine_nt = base_nt
  
  scenarios <- bind_rows(expand_grid(sigu_multiplier = coarse_su,
                                     sigv_multiplier = coarse_sv,
                                     rho             = coarse_rho,
                                     nT              = coarse_nt),
                         expand_grid(sigu_multiplier = fine_su,
                                     sigv_multiplier = base_sv,
                                     rho             = base_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = fine_sv,
                                     rho             = base_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = base_sv,
                                     rho             = fine_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = base_sv,
                                     rho             = base_rho,
                                     nT              = fine_nt)) %>%
    distinct() %>%
    mutate(world = 'medicare',
           ate_prior_sd = 20/sqrt(2/base::pi),
           nC = nT*2,
           scenario = glue('{world}_sigu{sigu_multiplier}_sigv{sigv_multiplier}_rho{str_replace(rho,"-","m")}_nT{nT}'))
  
  stopifnot(nrow(scenarios) == n_distinct(scenarios$scenario))
  
  #scenarios <- filter(scenarios, sigu_multiplier==1 & sigv_multiplier==1 & rho==0)
  
  test_list <- expand_grid(scenarios,
                           sim=1:200) %>%
    mutate(nburn   = 1000,
           nsim    = 1000,
           nthin   = 1,
           n_cores = 1,
           s3_path = 'bcf-sim-study/mar-2023-uiw/1kburn-1ksim')
  
  #Add seeds before we duplicate for ibcf vs wbcf
  stopifnot(length(seeds) >= nrow(test_list))
  test_list$seed <- seeds[1:nrow(test_list)]
  
  test_list <- test_list %>%
    expand_grid(model=c('ibcf','ubcf','wbcf')) %>%
    mutate(ibcf = model!='wbcf',
           ubcf = model=='ubcf',
           obcf = model=='obcf',
           fname = glue('{scenario}_{model}_mar_2023_uiw_sim{str_pad(sim,width = 3, side="left",pad = "0")}'),
           core = sigu_multiplier %in% coarse_su & sigv_multiplier %in% coarse_sv & rho %in% coarse_rho & nT %in% coarse_nt,
           base = sigu_multiplier == base_su & sigv_multiplier == base_sv & rho == base_rho & nT == base_nt) %>%
    select(-scenario, -sim, -model, -core, -base)
}

make_mar_2023_corrd_w_control_file <- function() {
  if (running_local) {
    seeds <- c(readLines('data/10k_random_ints_20211202.txt'),
               readLines('data/10k_random_ints_20211203.txt'),
               readLines('data/10k_random_ints_20211204.txt'))
  } else {
    seeds <- c(readLines('/home/bcf-sim-infra/data/10k_random_ints_20211202.txt'),
               readLines('/home/bcf-sim-infra/data/10k_random_ints_20211203.txt'),
               readLines('/home/bcf-sim-infra/data/10k_random_ints_20211204.txt'))
  }
  
  seeds <- seeds %>%
    as.numeric %>%
    unique
  
  base_su <- 1
  base_sv <- 1
  base_rho <- 0
  base_nt = 1000
  
  coarse_su <- base_su
  coarse_sv <- base_sv
  coarse_rho <- base_rho
  coarse_nt <- base_nt
  
  fine_su <- c(0,1,2,4)
  fine_sv <- base_sv
  fine_rho <- base_rho
  fine_nt <- base_nt
  
  scenarios <- bind_rows(expand_grid(sigu_multiplier = coarse_su,
                                     sigv_multiplier = coarse_sv,
                                     rho             = coarse_rho,
                                     nT              = coarse_nt),
                         expand_grid(sigu_multiplier = fine_su,
                                     sigv_multiplier = base_sv,
                                     rho             = base_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = fine_sv,
                                     rho             = base_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = base_sv,
                                     rho             = fine_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = base_sv,
                                     rho             = base_rho,
                                     nT              = fine_nt)) %>%
    distinct() %>%
    expand_grid(r2_w_pi=c(0,.25,.5),
                r2_w_mu=c(0,.25,.5),
                r2_w_tau=c(0,.25,.5),
                w_as_covar=c(TRUE,FALSE)) %>%
    mutate(world = 'medicare',
           ate_prior_sd = 20/sqrt(2/base::pi),
           nC = nT*2,
           scenario = glue('{world}_sigu{sigu_multiplier}_sigv{sigv_multiplier}_rho{str_replace(rho,"-","m")}_nT{nT}_r2wpi{r2_w_pi}_r2wmu{r2_w_mu}_r2wtau{r2_w_tau}_wcovar{w_as_covar}'))
  
  stopifnot(nrow(scenarios) == n_distinct(scenarios$scenario))
  
  #scenarios <- filter(scenarios, sigu_multiplier==1 & sigv_multiplier==1 & rho==0)
  
  test_list <- expand_grid(scenarios,
                           sim=1:50) %>%
    mutate(nburn   = 1000,
           nsim    = 1000,
           nthin   = 1,
           n_cores = 1,
           s3_path = 'bcf-sim-study/mar-2023-wcorr/1kburn-1ksim')
  
  #Add seeds before we duplicate for ibcf vs wbcf
  stopifnot(length(seeds) >= nrow(test_list))
  test_list$seed <- seeds[1:nrow(test_list)]
  
  test_list <- test_list %>%
    expand_grid(model=c('ubcf','wbcf')) %>%
    mutate(ibcf = model!='wbcf',
           ubcf = model=='ubcf',
           obcf = model=='obcf',
           fname = glue('{scenario}_{model}_mar_2023_wcorr_sim{str_pad(sim,width = 3, side="left",pad = "0")}'),
           core = sigu_multiplier %in% coarse_su & sigv_multiplier %in% coarse_sv & rho %in% coarse_rho & nT %in% coarse_nt,
           base = sigu_multiplier == base_su & sigv_multiplier == base_sv & rho == base_rho & nT == base_nt) %>%
    select(-scenario, -sim, -model, -core, -base)
}

make_apr_2023_r2e_control_file <- function() {
  if (running_local) {
    seeds <- c(readLines('data/10k_random_ints_20211202.txt'),
               readLines('data/10k_random_ints_20211203.txt'),
               readLines('data/10k_random_ints_20211204.txt'))
  } else {
    seeds <- c(readLines('/home/bcf-sim-infra/data/10k_random_ints_20211202.txt'),
               readLines('/home/bcf-sim-infra/data/10k_random_ints_20211203.txt'),
               readLines('/home/bcf-sim-infra/data/10k_random_ints_20211204.txt'))
  }
  
  seeds <- seeds %>%
    as.numeric %>%
    unique
  
  base_su <- 1
  base_sv <- 1
  base_rho <- 0
  base_nt = 1000
  
  coarse_su <- base_su
  coarse_sv <- base_sv
  coarse_rho <- base_rho
  coarse_nt <- base_nt
  
  fine_su <- 1
  fine_sv <- 1
  fine_rho <- base_rho
  fine_nt <- base_nt
  
  scenarios <- bind_rows(expand_grid(sigu_multiplier = coarse_su,
                                     sigv_multiplier = coarse_sv,
                                     rho             = coarse_rho,
                                     nT              = coarse_nt),
                         expand_grid(sigu_multiplier = fine_su,
                                     sigv_multiplier = base_sv,
                                     rho             = base_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = fine_sv,
                                     rho             = base_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = base_sv,
                                     rho             = fine_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = base_sv,
                                     rho             = base_rho,
                                     nT              = fine_nt)) %>%
    distinct() %>%
    expand_grid(r2_error=seq(0,1,.1)) %>%
    mutate(world = 'medicare',
           ate_prior_sd = 20/sqrt(2/base::pi),
           nC = nT*2,
           scenario = glue('{world}_r2e{r2_error}_sigu{sigu_multiplier}_sigv{sigv_multiplier}'))
  
  stopifnot(nrow(scenarios) == n_distinct(scenarios$scenario))
  
  #scenarios <- filter(scenarios, sigu_multiplier==1 & sigv_multiplier==1 & rho==0)
  
  test_list <- expand_grid(scenarios,
                           sim=1:250) %>%
    mutate(nburn   = 1000,
           nsim    = 1000,
           nthin   = 1,
           n_cores = 1,
           s3_path = 'bcf-sim-study/apr-2023-r2e/1kburn-1ksim')
  
  #Add seeds before we duplicate for ibcf vs wbcf
  stopifnot(length(seeds) >= nrow(test_list))
  test_list$seed <- seeds[1:nrow(test_list)]
  
  test_list <- test_list %>%
    expand_grid(model=c('ubcf','wbcf')) %>%
    mutate(ibcf = model!='wbcf',
           ubcf = model=='ubcf',
           obcf = model=='obcf',
           fname = glue('{scenario}_{model}_apr_2023_r2e_sim{str_pad(sim,width = 3, side="left",pad = "0")}'),
           core = sigu_multiplier %in% coarse_su & sigv_multiplier %in% coarse_sv & rho %in% coarse_rho & nT %in% coarse_nt,
           base = sigu_multiplier == base_su & sigv_multiplier == base_sv & rho == base_rho & nT == base_nt) %>%
    select(-scenario, -sim, -model, -core, -base)
}

make_stan_comp_file <- function() {
  if (running_local) {
    seeds <- c(readLines('data/10k_random_ints_20211202.txt'),
               readLines('data/10k_random_ints_20211203.txt'))
  } else {
    seeds <- c(readLines('/home/bcf-sim-infra/data/10k_random_ints_20211202.txt'),
               readLines('/home/bcf-sim-infra/data/10k_random_ints_20211203.txt'))
  }
  
  seeds <- seeds %>%
    as.numeric %>%
    unique
  
  base_su <- 1
  base_sv <- 1
  base_rho <- 0
  base_nt = 1000
  
  scenarios <- expand_grid(sigu_multiplier = base_su,
                           sigv_multiplier = base_sv,
                           rho             = base_rho,
                           nT              = base_nt,
                           nonlin_frac     = c(0, .1, 1/3, 2/3)) %>%
    distinct() %>%
    mutate(world = 'medicare',
           nC = nT*2,
           scenario = glue('{world}_sigu{sigu_multiplier}_sigv{sigv_multiplier}_rho{str_replace(rho,"-","m")}_nT{nT}_nonlin{round(nonlin_frac,2)}'))
  
  stopifnot(nrow(scenarios) == n_distinct(scenarios$scenario))
  
  test_list <- expand_grid(scenarios,
                           sim=1:200) %>%
    mutate(nburn   = 1000,
           nsim    = 1000,
           nthin   = 1,
           n_cores = 1,
           s3_path = 'bcf-sim-study/feb-2023-nonlin/1kburn-1ksim')
  
  #Add seeds before we duplicate for ibcf vs wbcf
  stopifnot(length(seeds) >= nrow(test_list))
  test_list$seed <- seeds[1:nrow(test_list)]
  
  test_list <- test_list %>%
    expand_grid(model=c('ibcf','wbcf','stan')) %>%
    mutate(ibcf = model=='ibcf',
           ubcf = model=='ubcf',
           obcf = model=='obcf',
           isstan = model=='stan',
           fname = glue('{scenario}_{model}_feb_2023_sim{str_pad(sim,width = 3, side="left",pad = "0")}'),
           ate_prior_sd = ifelse(model=='stan',20, 20/sqrt(2/pi))) %>%
    select(-scenario, -sim, -model)
}

make_hardcoded_rho_control_file <- function() {
  if (running_local) {
    seeds <- readLines('data/10k_random_ints_20211204.txt')
  } else {
    seeds <- readLines('/home/bcf-sim-infra/data/10k_random_ints_20211204.txt')
  }
  
  seeds <- seeds %>%
    as.numeric %>%
    unique
  
  base_su <- 1
  base_sv <- 1
  base_rho <- 0
  base_nt = 1000
  
  scenarios <- expand_grid(sigu_multiplier = base_su,
                           sigv_multiplier = base_sv,
                           rho             = base_rho,
                           nT              = base_nt,
                           hardcode_rho=TRUE) %>%
    distinct() %>%
    mutate(world = 'medicare',
           ate_prior_sd = 20/sqrt(2/pi),
           nC = nT*2,
           scenario = glue('{world}_sigu{sigu_multiplier}_sigv{sigv_multiplier}_rhohard0_nT{nT}'))
  
  stopifnot(nrow(scenarios) == n_distinct(scenarios$scenario))
  
  #scenarios <- filter(scenarios, sigu_multiplier==1 & sigv_multiplier==1 & rho==0)
  
  test_list <- expand_grid(scenarios,
                           sim=1:200) %>%
    mutate(nburn   = 1000,
           nsim    = 1000,
           nthin   = 1,
           n_cores = 1,
           s3_path = 'bcf-sim-study/jan-2023-runs/1kburn-1ksim-rhohard0')
  
  #Add seeds before we duplicate for ibcf vs wbcf
  stopifnot(length(seeds) >= nrow(test_list))
  test_list$seed <- seeds[1:nrow(test_list)]
  
  test_list <- test_list %>%
    expand_grid(ibcf=c(TRUE)) %>%
    mutate(fname = case_when(ibcf ~ paste0(scenario, '_ibcf_jan_2023_sim', str_pad(sim,width = 3, side='left',pad = '0')),
                             !ibcf ~ paste0(scenario, '_wbcf_jan_2023_sim', str_pad(sim,width = 3, side='left',pad = '0')))) %>%
    select(-scenario, -sim)
}

make_aug_2023_simple_control_file <- function() {
  if (running_local) {
    seeds <- c(readLines('data/10k_random_ints_20211202.txt'),
               readLines('data/10k_random_ints_20211203.txt'))
  } else {
    seeds <- c(readLines('/home/bcf-sim-infra/data/10k_random_ints_20211202.txt'),
               readLines('/home/bcf-sim-infra/data/10k_random_ints_20211203.txt'))
  }
  
  seeds <- seeds %>%
    as.numeric %>%
    unique
  
  base_su <- 1
  base_sv <- 1
  base_rho <- 0
  base_nt = 1000
  
  coarse_su <- base_su
  coarse_sv <- base_sv
  coarse_rho <- base_rho
  coarse_nt = base_nt
  
  fine_su <- c(0,.25,.5,2/3,1,1.5,2,4)
  fine_sv <- c(0,.25,.5,2/3,1,1.5,2,4)
  fine_rho <- base_rho
  fine_nt = base_nt
  
  scenarios <- bind_rows(expand_grid(sigu_multiplier = coarse_su,
                                     sigv_multiplier = coarse_sv,
                                     rho             = coarse_rho,
                                     nT              = coarse_nt),
                         expand_grid(sigu_multiplier = fine_su,
                                     sigv_multiplier = base_sv,
                                     rho             = base_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = fine_sv,
                                     rho             = base_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = base_sv,
                                     rho             = fine_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = base_sv,
                                     rho             = base_rho,
                                     nT              = fine_nt)) %>%
    distinct() %>%
    mutate(world = 'medicare',
           ate_prior_sd = 20/sqrt(2/base::pi),
           nC = nT*2,
           scenario = glue('{world}_sigu{sigu_multiplier}_sigv{sigv_multiplier}_rho{str_replace(rho,"-","m")}_nT{nT}'))
  
  stopifnot(nrow(scenarios) == n_distinct(scenarios$scenario))
  
  #scenarios <- filter(scenarios, sigu_multiplier==1 & sigv_multiplier==1 & rho==0)
  
  test_list <- expand_grid(scenarios,
                           sim=1:200) %>%
    mutate(nburn   = 1000,
           nsim    = 1000,
           nthin   = 1,
           n_cores = 1,
           s3_path = 'bcf-sim-study/aug-2023-uold/1kburn-1ksim')
  
  #Add seeds before we duplicate for ibcf vs wbcf
  stopifnot(length(seeds) >= nrow(test_list))
  test_list$seed <- seeds[1:nrow(test_list)]
  
  test_list <- test_list %>%
    expand_grid(model=c('ubcf','oldbcf')) %>%
    filter((sigu_multiplier==base_su & sigv_multiplier==base_sv) | model=='ubcf') %>%
    mutate(ibcf = model=='ubcf',
           ubcf = model=='ubcf',
           obcf = FALSE,
           oldbcf = model=='oldbcf',
           fname = glue('{scenario}_{model}_aug-2023-uold_sim{str_pad(sim,width = 3, side="left",pad = "0")}'),
           core = sigu_multiplier %in% coarse_su & sigv_multiplier %in% coarse_sv & rho %in% coarse_rho & nT %in% coarse_nt,
           base = sigu_multiplier == base_su & sigv_multiplier == base_sv & rho == base_rho & nT == base_nt) %>%
    select(-scenario, -sim, -model, -core, -base)
}

make_sep_2023_hyperprior_control_file <- function() {
  if (running_local) {
    seeds <- c(readLines('data/10k_random_ints_20211202.txt'),
               readLines('data/10k_random_ints_20211203.txt'))
  } else {
    seeds <- c(readLines('/home/bcf-sim-infra/data/10k_random_ints_20211202.txt'),
               readLines('/home/bcf-sim-infra/data/10k_random_ints_20211203.txt'))
  }
  
  seeds <- seeds %>%
    as.numeric %>%
    unique
  
  base_su <- 1
  base_sv <- 1
  base_rho <- 0
  base_nt = 1000
  
  coarse_su <- base_su
  coarse_sv <- base_sv
  coarse_rho <- base_rho
  coarse_nt = base_nt
  
  fine_su <- base_su
  fine_sv <- base_sv
  fine_rho <- base_rho
  fine_nt = base_nt
  
  scenarios <- bind_rows(expand_grid(sigu_multiplier = coarse_su,
                                     sigv_multiplier = coarse_sv,
                                     rho             = coarse_rho,
                                     nT              = coarse_nt),
                         expand_grid(sigu_multiplier = fine_su,
                                     sigv_multiplier = base_sv,
                                     rho             = base_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = fine_sv,
                                     rho             = base_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = base_sv,
                                     rho             = fine_rho,
                                     nT              = base_nt),
                         expand_grid(sigu_multiplier = base_su,
                                     sigv_multiplier = base_sv,
                                     rho             = base_rho,
                                     nT              = fine_nt)) %>%
    distinct() %>%
    expand_grid(sigu_pcthyperprior=2/3*c(.25,.5,2,4)) %>%
    mutate(world = 'medicare',
           ate_prior_sd = 20/sqrt(2/base::pi),
           nC = nT*2,
           scenario = glue('{world}_sigu{sigu_multiplier}_sigv{sigv_multiplier}_rho{str_replace(rho,"-","m")}_nT{nT}_siguhp{sigu_pcthyperprior}'))
  
  stopifnot(nrow(scenarios) == n_distinct(scenarios$scenario))
  
  #scenarios <- filter(scenarios, sigu_multiplier==1 & sigv_multiplier==1 & rho==0)
  
  test_list <- expand_grid(scenarios,
                           sim=1:200) %>%
    mutate(nburn   = 1000,
           nsim    = 1000,
           nthin   = 1,
           n_cores = 1,
           s3_path = 'bcf-sim-study/sep-2023-siguhp/1kburn-1ksim')
  
  #Add seeds before we duplicate for ibcf vs wbcf
  stopifnot(length(seeds) >= nrow(test_list))
  test_list$seed <- seeds[1:nrow(test_list)]
  
  test_list <- test_list %>%
    expand_grid(model=c('ubcf')) %>%
    filter((sigu_multiplier==base_su & sigv_multiplier==base_sv) | model=='ubcf') %>%
    mutate(ibcf = model=='ubcf',
           ubcf = model=='ubcf',
           obcf = FALSE,
           oldbcf = model=='oldbcf',
           fname = glue('{scenario}_{model}_sep-2023-siguhp_sim{str_pad(sim,width = 3, side="left",pad = "0")}'),
           core = sigu_multiplier %in% coarse_su & sigv_multiplier %in% coarse_sv & rho %in% coarse_rho & nT %in% coarse_nt,
           base = sigu_multiplier == base_su & sigv_multiplier == base_sv & rho == base_rho & nT == base_nt) %>%
    select(-scenario, -sim, -model, -core, -base)
}

make_oct2023_paper_uhp_ctrl <- function() {
  make_aug_2023_simple_control_file() %>%
    filter(sigu_multiplier==1 & sigv_multiplier==1 & ubcf) %>%
    expand_grid(sigu_pcthyperprior=2/3*c(.25,.5,2,4)) %>%
    mutate(model = ifelse(oldbcf,'oldbcf','ubcf'),
           sim = row_number(),
           scenario = glue('{world}_sigu{sigu_multiplier}_sigv{sigv_multiplier}_rho{str_replace(rho,"-","m")}_nT{nT}_siguhp{sigu_pcthyperprior}'),
           s3_path = 'bcf-sim-study/oct-2023-paper-uhp/1kburn-1ksim',
           fname = glue('{scenario}_{model}_oct-2023-paper-uhp_sim{str_pad(sim,width = 3, side="left",pad = "0")}')) %>%
    select(-scenario, -sim, -model)
}

make_oct2023_paper_icc_ctrl <- function() {
  make_aug_2023_simple_control_file() %>%
    filter(sigu_multiplier==1 & sigv_multiplier==1) %>%
    expand_grid(r2_error=seq(0,1,.1)) %>%
    mutate(model = ifelse(oldbcf,'oldbcf','ubcf'),
           sim = row_number(),
           scenario = glue('{world}_sigu{sigu_multiplier}_sigv{sigv_multiplier}_rho{str_replace(rho,"-","m")}_nT{nT}_r2e{r2_error}'),
           s3_path = 'bcf-sim-study/oct-2023-paper-icc/1kburn-1ksim',
           fname = glue('{scenario}_{model}_oct-2023-paper-icc_sim{str_pad(sim,width = 3, side="left",pad = "0")}')) %>%
    select(-scenario, -sim, -model)
}

make_oct2023_paper_ibcf <- function() {
  make_aug_2023_simple_control_file() %>%
    filter(sigu_multiplier==1 & sigv_multiplier==1 & ubcf) %>%
    mutate(model = 'ibcf',
           ubcf = FALSE,
           sim = row_number(),
           scenario = glue('{world}_sigu{sigu_multiplier}_sigv{sigv_multiplier}_rho{str_replace(rho,"-","m")}_nT{nT}'),
           s3_path = 'bcf-sim-study/oct-2023-paper-ibcf/1kburn-1ksim',
           fname = glue('{scenario}_{model}_oct-2023-paper-ibcf_sim{str_pad(sim,width = 3, side="left",pad = "0")}')) %>%
    select(-scenario, -sim, -model)
}

simple_AWS_run <- function(s3_path, fname, ...) {
  result <- test_bcf(...)
  saveRDS(result,'/home/results.RDS')
  system(glue('aws s3 cp /home/results.RDS s3://{s3_path}/{fname}.RDS'))
  return(NULL)
}

AWS_run_failsafe <- function(s3_path, fname, isstan=FALSE, ...) {
  tryCatch({
    if (!isstan) {
      result <- test_bcf(...)  
    } else {
      result <- test_stan(...)
    }
    
    temp <- tempfile(fileext='.RDS')
    saveRDS(result,temp)
    system(glue('aws s3 cp {temp} s3://{s3_path}/{fname}.RDS'))
    #Now save small version
    temp <- tempfile(fileext='.RDS')
    if (!isstan) {
      saveRDS(result[c('overall','inputs','mixing','exemplar_summy','postcorr','scalecorr','acceptance','calib')], temp)  
    } else {
      saveRDS(result[c('overall','inputs','mixing','exemplar_summy')], temp)  
    }
    
    system(glue('aws s3 cp {temp} s3://{s3_path}/{fname}_small.RDS'))
  }, error=function(e){
    temp <- tempfile(fileext='.txt')
    writeLines(as.character(e),temp) 
    system(glue('aws s3 cp {temp} s3://{s3_path}/{fname}.txt'))  
  })
  return(NULL)
}

aws_run_multiple <- function(processes, control_file) {
  processes <- as.numeric(processes)
  
  #Download control file
  system(glue('aws s3 cp {control_file} /home/control_file.RDS'))
  pars <- readRDS('/home/control_file.RDS')
  
  library(furrr)
  plan(multicore, workers=processes)
  future_pmap(pars, AWS_run_failsafe)
  return(NULL)
}

aws_run_multiple_make_ctrl <- function(processes, control_fn, start, end) {
  processes <- as.numeric(processes)
  start     <- as.numeric(start)
  end       <- as.numeric(end)
  
  pars <- get(control_fn)()[start:end,]
  
  #Precompile stan model if desired
  if ('isstan' %in% colnames(pars)) {
    if (any(pars$isstan)) {
      #model <- compile_model('linear_stan_model')
    }
  }
  
  library(furrr)
  plan(multicore, workers=processes)
  future_pmap(pars, AWS_run_failsafe)
  return(NULL)
}

even_simpler_test <- function(s3_path,fname) {
  writeLines(c('abc','123'),'/home/results.txt')
  system(glue('aws s3 cp /home/results.txt s3://{s3_path}/{fname}.txt'))
  system(glue('aws s3 cp /opt/ml/input/config/hyperparameters.json s3://{s3_path}/hyperpars.json'))
  return(NULL)
}
