AWS_run_failsafe_acic <- function(year, par, sim, s3_path, ...) {
  tryCatch({
    test_acic(year=year, par=par, sim=sim, s3_path=s3_path, ...)
  }, error=function(e){
    temp <- tempfile(fileext='.txt')
    writeLines(as.character(e),temp) 
    system(glue('aws s3 cp {temp} s3://{s3_path}/acic{year}_p{str_pad(par,width=2,side="left",pad="0")}_s{str_pad(sim,width=3,side="left",pad="0")}.txt'))  
  })
  return(NULL)
}

aws_run_multiple_acic <- function(processes, control_file) {
  processes <- as.numeric(processes)
  
  #Download control file
  system(glue('aws s3 cp {control_file} /home/control_file.RDS'))
  pars <- readRDS('/home/control_file.RDS')
  
  library(furrr)
  plan(multicore, workers=processes)
  future_pmap(pars, AWS_run_failsafe_acic)
  return(NULL)
}

test_acic <- function(year, par, sim, s3_path=NULL, ...) {
  year <- as.numeric(year)
  par <- as.numeric(par)
  sim <- as.numeric(sim)
  
  if (year==2016) {
    x <- aciccomp2016::input_2016
    #z is treat, y is observed, y.0/y.1 are potential outcomes, and e is pscore
    #No idea what mu.0/mu.1 are - maybe true mu so denoised? But then why does mu vary with t/c?
    yz <- as.data.frame(aciccomp2016::dgp_2016(x, par, sim))
    yz$tau <- yz$y.1 - yz$y.0
    yz$mu <- yz$y.0
  } else {
    x <- aciccomp2017::input_2017
    #Just has z/y/alpha
    #I think alpha is treatment effect
    yz <- as.data.frame(aciccomp2017::dgp_2017(par, sim))
    yz$tau <- yz$alpha
    yz$mu <- yz$y - yz$z*yz$tau
    #Don't have true pscore so just blank
    yz$e <- NA_real_
  }
  
  fml <- as.formula(glue::glue('z ~ {paste(colnames(x),collapse="+")}'))
  xmat <- model.matrix(fml,cbind(x,yz))
  
  #generate a pscore - good news is we have real underlying propensity, so we can check how bad we are at pscore vs fit
  #logit <- glm(fml, data=cbind(x,yz), family=binomial)
  #yz$pihat <- logit$fitted.values
  bartfit <- dbarts::bart(xmat, yz$z)
  yz$pihat <- bartfit$yhat.train %>% pnorm %>% apply(2,mean)
  
  info <- tibble(year=year,
                 par=par,
                 sim=sim,
                 pscore_rmse = sqrt(mean(yz$pihat - yz$e)^2),
                 pscore_corr = cor(yz$pihat, yz$e),
                 pscore_rmse_z = sqrt(mean(yz$pihat - yz$z)^2),
                 pscore_corr_z = cor(yz$pihat, yz$z))
  
  ibcf <- acic_fit_bcf(yz, xmat, ibcf=TRUE, ...)
  gc()
  
  vbcf <- acic_fit_bcf(yz, xmat, ibcf=FALSE, ...)
  gc()
  
  ret <- list(info=info,
              ibcf=ibcf,
              vbcf=vbcf)
  
  if (!is.null(s3_path)) {
    temp <- tempfile(fileext='.RDS')
    saveRDS(ret,temp)
    system(glue('aws s3 cp {temp} s3://{s3_path}/acic{year}_p{str_pad(par,width=2,side="left",pad="0")}_s{str_pad(sim,width=3,side="left",pad="0")}.RDS'))
  }
  
  return(ret)
}

acic_fit_bcf <- function(yz, xmat, ibcf=FALSE, nburn=1000, nsim=1000, nthin=1, n_chains=4, n_cores=n_chains) {
  tictoc::tic()
  #Only hyperpar change is including pscore in both mod and con, since Jared did that for his ACIC runs
  fit <- bcf(y           = yz$y,
             z           = yz$z,
             x_control   = xmat,
             x_moderate  = xmat,
             pihat       = yz$pihat,
             include_pi  = 'both',
             nburn       = as.numeric(nburn),
             nsim        = as.numeric(nsim),
             nthin       = as.numeric(nthin),
             n_chains    = as.numeric(n_chains),
             n_cores     = as.numeric(n_cores),
             n_threads   = 1,
             block_v_rho = TRUE,
             include_random_effects=ibcf,
             simplified_return = TRUE,
             verbose=0)
  timing <- tictoc::toc()
  
  stanchains <- stanify(fit$raw_chains, rep(1, nrow(xmat)), re=ibcf)
  mixing <- rstan::monitor(stanchains,warmup=0,print=FALSE,probs=c(.025,.05,.1,.9,.95,.975)) %>%
    as_tibble(rownames='par') %>%
    mutate(real = case_when(par=='tau_bar'    ~ mean(yz$tau),
                            par=='mu_bar'     ~ mean(yz$y.0),
                            par=='yhat_bar'   ~ mean(yz$y)),
           cover80 = real>=`10%` & real<=`90%`,
           cover90 = real>=`5%` & real<=`95%`,
           cover95 = real>=`2.5%` & real<=`97.5%`,
           width80 = `90%`-`10%`,
           width90 = `95%`-`5%`,
           width95 = `97.5%`-`2.5%`) %>%
    select(par, real, mean, sd,
           cover80, cover90, cover95, 
           width80, width90, width95, 
           p025=`2.5%`, p5 = `5%`, p10=`10%`,
           p90=`90%`, p95=`95%`, p975=`97.5%`,
           n_eff, Rhat)
  
  #What do we care about? ATE, ATT?
  # BCF paper talks about coverage and interval length, but I'm not clear on what (SATT I think?), 
  # avg/sd of bias, avg/sd of |bias|, and average RMSE of "CATE for each unit" - I assume that means RMSE on tau for each indiv?
  taus <- lapply(fit$raw_chains,`[[`,'tau') %>% do.call(what=rbind)
  SATT <- mean(yz$tau[yz$z==1])
  SATT_draws <- apply(taus[,yz$z==1],1,mean)
  SATT_hat <- mean(SATT_draws)
  CATE_hat <- apply(taus,2,mean)
  lb80 <- quantile(SATT_draws,.1)
  ub80 <- quantile(SATT_draws,.9)
  lb90 <- quantile(SATT_draws,.05)
  ub90 <- quantile(SATT_draws,.95)
  lb95 <- quantile(SATT_draws,.025)
  ub95 <- quantile(SATT_draws,.975)
  
  overall <- tibble(SATT=SATT,
                    SATT_hat=SATT_hat,
                    cover80 = SATT>=lb80 & SATT<=ub80,
                    cover90 = SATT>=lb90 & SATT<=ub90,
                    cover95 = SATT>=lb95 & SATT<=ub95,
                    width80 = ub80-lb80,
                    width90 = ub90-lb90,
                    width95 = ub95-lb95,
                    bias = SATT_hat - SATT,
                    PEHE = sqrt(mean((CATE_hat - yz$tau)^2)),
                    PEHT = sqrt(mean((CATE_hat[yz$z==1] - yz$tau[yz$z==1])^2)),
                    p025=lb95,
                    p5=lb90,
                    p10=lb80,
                    p90=ub80,
                    p95=ub90,
                    p975=ub95,
                    timing = timing$toc -timing$tic)
  
  indiv <- tibble(z=yz$z,
                  tau=yz$tau,
                  tau_hat = CATE_hat,
                  lb80 = apply(taus,2,quantile,.1),
                  ub80 = apply(taus,2,quantile,.9),
                  lb90 = apply(taus,2,quantile,.05),
                  ub90 = apply(taus,2,quantile,.95),
                  lb95 = apply(taus,2,quantile,.025),
                  ub95 = apply(taus,2,quantile,.975)) %>%
    mutate(cover80 = tau>=lb80 & tau<=ub80,
           cover90 = tau>=lb90 & tau<=ub90,
           cover95 = tau>=lb95 & tau<=ub95,
           width80 = ub80-lb80,
           width90 = ub90-lb90,
           width95 = ub95-lb95) %>%
    rename(p025=lb95,
           p5=lb90,
           p10=lb80,
           p90=ub80,
           p95=ub90,
           p975=ub95) %>%
    select(z, tau, tau_hat, matches('cover'), matches('width'), p025, p5, p10, p90, p95, p975)
  
  return(list(mixing=mixing,
              overall=overall,
              indiv=indiv))
}
