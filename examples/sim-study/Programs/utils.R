library(tidyverse)
library(glue)

pal=list(green='#17A673',
         light_green='#6AB790',
         teal='#189394',
         light_red='#D9654A',
         red='#D02B27',
         gold='#F1B51C',
         purple='#5C4377',
         chocolate='#753726',
         gray='#5B6771',
         light_gray='#9FA2A9',
         beige='#E0D4B5',
         navy='#0B2949',
         tim_gray='#C8C9CC')

drop_identifiers <- function(data, all=TRUE) {
  if (all) {
    pardrop <- c('trt_eff_scenario', 'sig_u', 'sig_v', 'rho', 'weights', 'uv_dist', 'N', 'uv_dist_df', 'n_chains', 'nburn', 'nsim', 'nthin', 'seed')
    todrop <- intersect(pardrop, colnames(data))
    data %>% 
      select(-all_of(todrop)) %>% 
      return
  } else {
    pardrop <- c('nC', 'trt_eff_scenario', 'weights', 'uv_dist', 'uv_dist_df', 'n_chains', 'nburn', 'nsim', 'nthin', 'seed')
    todrop <- intersect(pardrop, colnames(data))
    data %>% 
      select(-all_of(todrop)) %>% 
      return
  }
}

create_group <- function(data, drop=TRUE, all=TRUE) {
  data <- data %>%
    mutate(scenario = case_when(set %in% c('main_het','main_homog') ~ glue('{set}_sigv{sig_v}_rho{rho}'),
                                set == 'no_res' ~ glue('{set}_{trt_eff_scenario}'),
                                set == 'no_weights' ~ glue('{set}_{trt_eff_scenario}_{ifelse(sig_v==0, "no","")}res'),
                                set == 't_dist' ~ glue('{set}_{trt_eff_scenario}_sigv{sig_v}')) %>%
             as.character)
  if (drop) {
    data <- drop_identifiers(data, all)
  }
  return(data)
}

combine_combined_output <- function(prefixes, dirs, small=TRUE, big=TRUE, outfix, outdir) {
  stopifnot(length(prefixes) == length(dirs))
  n <- length(prefixes)
  
  if (small) {
    print(glue('Combining smalls'))
    lapply(1:n, function(i) readRDS(glue('{dirs[i]}/{prefixes[i]}-overall.RDS'))) %>%
      bind_rows %>%
      saveRDS(glue('{outdir}/{outfix}-overall.RDS'))
    
    lapply(1:n, function(i) readRDS(glue('{dirs[i]}/{prefixes[i]}-mixing.RDS'))) %>%
      bind_rows %>%
      saveRDS(glue('{outdir}/{outfix}-mixing.RDS'))
    
    lapply(1:n, function(i) readRDS(glue('{dirs[i]}/{prefixes[i]}-misc.RDS'))) %>%
      unlist(recursive = FALSE) %>%
      saveRDS(glue('{outdir}/{outfix}-misc.RDS'))
    
    lapply(1:n, function(i) readRDS(glue('{dirs[i]}/{prefixes[i]}-exemplar-summy.RDS'))) %>%
      bind_ex_summys() %>%
      saveRDS(glue('{outdir}/{outfix}-exemplar-summy.RDS'))
    
    lapply(1:n, function(i) readRDS(glue('{dirs[i]}/{prefixes[i]}-calib.RDS'))) %>%
      bind_rows %>%
      saveRDS(glue('{outdir}/{outfix}-calib.RDS'))
    gc()
  }
  
  if (big) {
    print(glue('Combining datas'))
    lapply(1:n, function(i) readRDS(glue('{dirs[i]}/{prefixes[i]}-data.RDS'))) %>%
      bind_rows %>%
      saveRDS(glue('{outdir}/{outfix}-data.RDS'))
    gc()
    
    print(glue('Combining indivs'))
    lapply(1:n, function(i) readRDS(glue('{dirs[i]}/{prefixes[i]}-indiv.RDS'))) %>%
      bind_rows %>%
      saveRDS(glue('{outdir}/{outfix}-indiv.RDS'))
    gc()
    
    print(glue('Combining indiv resids'))
    lapply(1:n, function(i) readRDS(glue('{dirs[i]}/{prefixes[i]}-indiv_resid.RDS'))) %>%
      bind_rows %>%
      saveRDS(glue('{outdir}/{outfix}-indiv_resid.RDS'))
    gc()
    
    print(glue('Combining indiv uvs'))
    lapply(1:n, function(i) readRDS(glue('{dirs[i]}/{prefixes[i]}-indiv_uv.RDS'))) %>%
      bind_rows %>%
      saveRDS(glue('{outdir}/{outfix}-indiv_uv.RDS'))
    gc()
    
    print(glue('Combining exemplars'))
    lapply(1:n, function(i) readRDS(glue('{dirs[i]}/{prefixes[i]}-exemplar.RDS'))) %>%
      bind_exemplars() %>%
      saveRDS(glue('{outdir}/{outfix}-exemplar.RDS'))
    gc()
  }
}

combine_output_chunks <- function(files, dir, save_small=TRUE, save_big=TRUE, max_size=1000, print_iter=100, skipfirst=0) {
  sets <- split(files, ceiling(1:nrow(files) / max_size))
  idir <- paste0(dir,'/interim')
  if (!dir.exists(idir)) dir.create(idir)
  
  for (i in 1:length(sets)) {
    if (i <= skipfirst) {
      print(glue('Already did {i}'))
    } else {
      combine_output(sets[[i]], prefix=paste0('set',i), dir=idir, save_small=save_small, save_big=save_big, print_iter=print_iter)  
    }
  }
  
  dirs <- rep(paste0(dir,'/interim'),length(sets))
  for (d in dirs) {
    if (!dir.exists(d)) dir.create(d)
  }
  
  combine_combined_output(prefixes = paste0('set',1:length(sets)), 
                          dirs=dirs,
                          outfix='all',
                          outdir=dir,
                          small=save_small, big=save_big)
}

combine_output <- function(files, prefix, dir, save_small=TRUE, save_big=TRUE, print_iter=100) {
  print(glue('{prefix}'))
  env <- environment()
  #Initialize output holders with NULLs, since they play nice with bind_rows
  miscs <- overalls <- mixings <- indivs <- indiv_uvs <- indiv_resids <- datas <- exemplars <- ex_summys <- calibs <- vector('list',nrow(files))
  
  for (i in 1:nrow(files)) {
    if (i==1 | i %% print_iter==0) print(glue('  {i}'))
    add_sim_to_output(files$file[i], files$set[i], files$scenario[i], files$sim[i], files$method[i], i, env, do_big=save_big)
  }
  
  if (save_small) {
    print(glue('  saving smalls'))
    saveRDS(miscs, glue('{dir}/{prefix}-misc.RDS'))
    overalls  %>% bind_rows %>% saveRDS(glue('{dir}/{prefix}-overall.RDS'))
    mixings   %>% bind_rows %>% saveRDS(glue('{dir}/{prefix}-mixing.RDS'))
    ex_summys %>% bind_ex_summys %>% saveRDS(glue('{dir}/{prefix}-exemplar-summy.RDS'))
    calibs    %>% bind_rows %>% saveRDS(glue('{dir}/{prefix}-calib.RDS'))
  }
  if (save_big) {
    print(glue('  saving datas'))
    datas     %>% bind_rows %>% saveRDS(glue('{dir}/{prefix}-data.RDS'))
    
    print(glue('  saving indivs'))
    indivs    %>% bind_rows %>% saveRDS(glue('{dir}/{prefix}-indiv.RDS'))
    
    print(glue('  saving indiv resids'))
    indiv_resids %>% bind_rows %>% saveRDS(glue('{dir}/{prefix}-indiv_resid.RDS'))
    
    print(glue('  saving indiv uvs'))
    indiv_uvs %>% bind_rows %>% saveRDS(glue('{dir}/{prefix}-indiv_uv.RDS'))
    
    print(glue('  saving exemplars'))
    exemplars %>% bind_exemplars %>% saveRDS(glue('{dir}/{prefix}-exemplar.RDS'))
  }
  
  miscs <- overalls <- mixings <- indivs <- indiv_uvs <- indiv_resids <- datas <- exemplars <- ex_summys <- NULL
  gc()
  
  return(NULL)
}

add_sim_to_output <- function(file, set, scenario, sim, method, i, e, do_big=TRUE) {
  fit <- readRDS(file)
  
  fit$inputs$set      <- set
  fit$inputs$scenario <- scenario
  fit$inputs$sim      <- sim
  fit$inputs$method   <- method
  
  fit$inputs <- drop_identifiers(fit$inputs, all=FALSE)
  
  #Add CATE/CATT/CATU widths to overall if not present
  if (is.null(fit$overall$CATEwidth80) & do_big) {
    widths <- fit$indiv %>%
      summarize(CATEwidth80 = mean(width80),
                CATEwidth90 = mean(width90),
                CATEwidth95 = mean(width95),
                CATTwidth80 = mean(ifelse(z==1,width80,NA_real_), na.rm=TRUE),
                CATTwidth90 = mean(ifelse(z==1,width90,NA_real_), na.rm=TRUE),
                CATTwidth95 = mean(ifelse(z==1,width95,NA_real_), na.rm=TRUE),
                CATUwidth80 = mean(ifelse(z==0,width80,NA_real_), na.rm=TRUE),
                CATUwidth90 = mean(ifelse(z==0,width90,NA_real_), na.rm=TRUE),
                CATUwidth95 = mean(ifelse(z==0,width95,NA_real_), na.rm=TRUE))
    
    fit$overall <- bind_cols(fit$overall, widths) 
  }
  if (is.null(fit$overall$Residwidth80) & do_big) {
    widths <- fit$indiv_resid %>%
      summarize(Residcover80 = mean(cover80),
                Residcover90 = mean(cover90),
                Residcover95 = mean(cover95),
                Residwidth80 = mean(width80),
                Residwidth90 = mean(width90),
                Residwidth95 = mean(width95),
                Residwidth80 = mean(ifelse(z==1,width80,NA_real_), na.rm=TRUE),
                Residwidth90 = mean(ifelse(z==1,width90,NA_real_), na.rm=TRUE),
                Residwidth95 = mean(ifelse(z==1,width95,NA_real_), na.rm=TRUE),
                Residwidth80 = mean(ifelse(z==0,width80,NA_real_), na.rm=TRUE),
                Residwidth90 = mean(ifelse(z==0,width90,NA_real_), na.rm=TRUE),
                Residwidth95 = mean(ifelse(z==0,width95,NA_real_), na.rm=TRUE))
    
    fit$overall <- bind_cols(fit$overall, widths) 
  }
  
  #Add taux/v to indiv
  if (do_big) {
    if (is.null(fit$indiv$taux)) {
      stopifnot(identical(fit$indiv$id, fit$data$id))
      fit$indiv$taux <- fit$data$taux
      fit$indiv$v    <- fit$data$v  
    }
  }
  
  e$overalls[[i]]      <- bind_cols(fit$inputs, fit$overall)
  e$mixings[[i]]       <- bind_cols(fit$inputs, fit$mixing)
  e$calibs[[i]]        <- bind_cols(fit$inputs, fit$calib)
  if (do_big) {
    e$indivs[[i]]        <- bind_cols(fit$inputs, fit$indiv)
    e$indiv_resids[[i]]  <- bind_cols(fit$inputs, fit$indiv_resid)
    e$datas[[i]]         <- bind_cols(fit$inputs, fit$data %>% select(-z_str, -sigma_u, -sigma_v, -rho, -sigma_y, -sigma)) 
    e$exemplars[[i]]     <- list(tau   = bind_cols(fit$inputs, fit$exemplar$tau),
                                 yhat  = bind_cols(fit$inputs, fit$exemplar$yhat),
                                 resid = bind_cols(fit$inputs, fit$exemplar$resid),
                                 v     = bind_cols(fit$inputs, fit$exemplar$v))
  }
  
  if (method %in% c('ibcf','ubcf','obcf')) {
    if (do_big) {
      e$indiv_uvs[[i]]   <- bind_cols(fit$inputs, fit$indiv_uv) 
    }
    e$miscs[[i]]  <- fit[c('inputs','postcorr','scalecorr','acceptance')]
  } else {
    e$miscs[[i]]  <- fit[c('inputs','scalecorr','acceptance')]
  }
  
  e$ex_summys[[i]]     <- add_input_to_ex_summy(fit$exemplar_summy, fit$inputs)
  
  fit <- NULL
  gc()
  
  return(NULL)
}

combine_exemplars <- function(orig, x, inputs) {
  orig$tau <- bind_rows(orig$tau, bind_cols(inputs, x$tau))
  orig$yhat <- bind_rows(orig$yhat, bind_cols(inputs, x$yhat))
  orig$resid <- bind_rows(orig$resid, bind_cols(inputs, x$resid))
  if (is.null(inputs) || inputs$ibcf) {
    orig$v <- bind_rows(orig$v, bind_cols(inputs, x$v))
  }
  
  return(orig)
}

add_input_to_ex_summy <- function(ex_summy, inputs) {
  ex_summy$tau$auc                <- bind_cols(inputs, ex_summy$tau$auc)
  ex_summy$tau$rmse               <- bind_cols(inputs, ex_summy$tau$rmse)
  ex_summy$tau$confus             <- bind_cols(inputs, ex_summy$tau$confus)
  
  ex_summy$yhat$auc                <- bind_cols(inputs, ex_summy$yhat$auc)
  ex_summy$yhat$rmse               <- bind_cols(inputs, ex_summy$yhat$rmse)
  ex_summy$yhat$confus             <- bind_cols(inputs, ex_summy$yhat$confus)
  
  ex_summy$resid$auc                <- bind_cols(inputs, ex_summy$resid$auc)
  ex_summy$resid$rmse               <- bind_cols(inputs, ex_summy$resid$rmse)
  ex_summy$resid$confus             <- bind_cols(inputs, ex_summy$resid$confus)
  
  ex_summy$v$auc                <- bind_cols(inputs, ex_summy$v$auc)
  ex_summy$v$rmse               <- bind_cols(inputs, ex_summy$v$rmse)
  ex_summy$v$confus             <- bind_cols(inputs, ex_summy$v$confus)
  return(ex_summy)
}
combine_ex_summys <- function(orig, x, inputs) {
  orig$tau = list(auc    = bind_rows(orig$tau$auc,    bind_cols(inputs, x$tau$auc)),
                  rmse   = bind_rows(orig$tau$rmse,   bind_cols(inputs, x$tau$rmse)),
                  confus = bind_rows(orig$tau$confus, bind_cols(inputs, x$tau$confus)))
  orig$yhat = list(auc    = bind_rows(orig$yhat$auc,    bind_cols(inputs, x$yhat$auc)),
                   rmse   = bind_rows(orig$yhat$rmse,   bind_cols(inputs, x$yhat$rmse)),
                   confus = bind_rows(orig$yhat$confus, bind_cols(inputs, x$yhat$confus)))
  orig$resid = list(auc    = bind_rows(orig$resid$auc,    bind_cols(inputs, x$resid$auc)),
                    rmse   = bind_rows(orig$resid$rmse,   bind_cols(inputs, x$resid$rmse)),
                    confus = bind_rows(orig$resid$confus, bind_cols(inputs, x$resid$confus)))
  if (is.null(inputs) || inputs$ibcf) {
    orig$v = list(auc    = bind_rows(orig$v$auc,    bind_cols(inputs, x$v$auc)),
                  rmse   = bind_rows(orig$v$rmse,   bind_cols(inputs, x$v$rmse)),
                  confus = bind_rows(orig$v$confus, bind_cols(inputs, x$v$confus)))
  }
  
  return(orig)
}

bind_exemplars <- function(exemplars) {
  exemplar <- list()
  exemplar$tau <- lapply(exemplars, `[[`, 'tau') %>% bind_rows
  exemplar$yhat <- lapply(exemplars, `[[`, 'yhat') %>% bind_rows
  exemplar$resid <- lapply(exemplars, `[[`, 'resid') %>% bind_rows
  exemplar$v   <- lapply(exemplars, `[[`, 'v')   %>% bind_rows
  
  return(exemplar)
}

bind_ex_summys <- function(ex_summys) {
  ex_summy <- list()
  ex_summy$tau$auc                <- lapply(ex_summys, function(x) x$tau$auc)                %>% bind_rows
  ex_summy$tau$rmse               <- lapply(ex_summys, function(x) x$tau$rmse)               %>% bind_rows
  ex_summy$tau$confus             <- lapply(ex_summys, function(x) x$tau$confus)             %>% bind_rows
  
  ex_summy$yhat$auc                <- lapply(ex_summys, function(x) x$yhat$auc)                %>% bind_rows
  ex_summy$yhat$rmse               <- lapply(ex_summys, function(x) x$yhat$rmse)               %>% bind_rows
  ex_summy$yhat$confus             <- lapply(ex_summys, function(x) x$yhat$confus)             %>% bind_rows
  
  ex_summy$resid$auc                <- lapply(ex_summys, function(x) x$resid$auc)                %>% bind_rows
  ex_summy$resid$rmse               <- lapply(ex_summys, function(x) x$resid$rmse)               %>% bind_rows
  ex_summy$resid$confus             <- lapply(ex_summys, function(x) x$resid$confus)             %>% bind_rows
  
  ex_summy$v$auc                  <- lapply(ex_summys, function(x) x$v$auc)                  %>% bind_rows
  ex_summy$v$rmse                 <- lapply(ex_summys, function(x) x$v$rmse)                 %>% bind_rows
  ex_summy$v$confus               <- lapply(ex_summys, function(x) x$v$confus)               %>% bind_rows
  
  return(ex_summy)
}

danfig1 <- function(overall, varying='sig_v', sigv_const=1, sigu_const=1, rho_const=0, nT_const=1000, estimand='SATT', ci=FALSE, const_cap=FALSE) {
  touse <- overall
  
  if (varying!='sig_v') {
    touse <- filter(touse, sigv_multiplier==sigv_const) 
  } else {
    xmult <- 8.33
    touse$xvar <- factor(round(xmult*touse$sigv_multiplier))
    xlabel <- expression(sigma[v])
  }
  if (varying!='sig_u') {
    touse <- filter(touse, sigu_multiplier==sigu_const) 
  } else {
    xmult <- 61.3
    touse$xvar <- factor(round(xmult*touse$sigu_multiplier))
    xlabel <- expression(sigma[u])
  }
  if (varying!='rho') {
    touse <- filter(touse, rho==rho_const) 
  } else {
    xmult <- 1
    touse$xvar <- factor(touse$rho)
    xlabel <- expression(rho)
  }
  if (varying!='nT') {
    touse <- filter(touse, nT==nT_const) 
  } else {
    xmult <- 1
    touse$xvar <- factor(touse$nT)
    xlabel <- '# treatment practices'
  }
  
  if (const_cap & varying=='sig_u') {
    cap <- bquote(sigma[v]~'='~.(sigv_const)*','~rho~'='~.(rho_const)*', nT ='~.(nT_const))
  } else if (const_cap & varying=='sig_v') {
    cap <- bquote(sigma[u]~'='~.(sigu_const)*','~rho~'='~.(rho_const)*', nT ='~.(nT_const))
  } else if (const_cap & varying=='rho') {
    cap <- bquote(sigma[u]~'='~.(sigu_const)*','~sigma[v]~'='~.(sigv_const)*', nT ='~.(nT_const))
  } else if (const_cap & varying=='nT') {
    cap <- bquote(sigma[u]~'='~.(sigu_const)*','~sigma[v]~'='~.(sigv_const)*','~rho~'='~.(rho_const))
  } else {
    cap=''
  }
  
  stopifnot(nrow(touse) == touse %>% select(method, xvar, sim) %>% distinct %>% nrow)
  
  if (estimand=='SATT') {
    smry <- touse %>%
      group_by(method, xvar) %>%
      summarize(cover_est = mean(SATTcover90),
                cover_lb  = qbinom(.05, n(), cover_est)/n(),
                cover_ub  = qbinom(.95, n(), cover_est)/n(),
                RMSE_est  = sqrt(mean(SATTbias^2)),
                #Does this make sense?
                RMSE_lb   = sqrt(quantile(SATTbias^2, .05)),
                RMSE_ub   = sqrt(quantile(SATTbias^2, .95)),
                width_est = mean(SATTwidth90),
                width_lb  = quantile(SATTwidth90, .05),
                width_ub  = quantile(SATTwidth90, .95))
  } else if (estimand=='CATT') {
    smry <- touse %>%
      group_by(method, xvar) %>%
      summarize(cover_est = mean(CATTcover90),
                cover_lb  = qbinom(.05, n(), cover_est)/n(),
                cover_ub  = qbinom(.95, n(), cover_est)/n(),
                #PEHT is CATT RMSE, so square it to get to MSE, then mean that, then root, so it's really RMMSE; 
                #if Ns are equal across what we're summarizing we should be fine, and I think they are if we're keeping nt fixed
                RMSE_est  = sqrt(mean(PEHT^2)),
                RMSE_lb   = sqrt(quantile(PEHT^2, .05)),
                RMSE_ub   = sqrt(quantile(PEHT^2, .95)),
                width_est = mean(CATTwidth90),
                width_lb  = quantile(CATTwidth90, .05),
                width_ub  = quantile(CATTwidth90, .95))
  } else if (estimand=='resid') {
    smry <- touse %>%
      group_by(method, xvar) %>%
      summarize(cover_est = mean(Residcover90),
                cover_lb  = qbinom(.05, n(), cover_est)/n(),
                cover_ub  = qbinom(.95, n(), cover_est)/n(),
                #PEHT is CATT RMSE, so square it to get to MSE, then mean that, then root, so it's really RMMSE; 
                #if Ns are equal across what we're summarizing we should be fine, and I think they are if we're keeping nt fixed
                RMSE_est  = sqrt(mean(RMSET_resid^2)),
                RMSE_lb   = sqrt(quantile(RMSET_resid^2, .05)),
                RMSE_ub   = sqrt(quantile(RMSET_resid^2, .95)),
                width_est = mean(Residwidth90),
                width_lb  = quantile(Residwidth90, .05),
                width_ub  = quantile(Residwidth90, .95))
  }
  
  smry <- smry %>%
    pivot_longer(matches('est|[ul]b'), names_sep = '_', names_to=c('metric','.value')) %>%
    mutate(truth = ifelse(metric=='cover',.9,NA_real_))
  
  plot <- ggplot(smry) + 
    geom_line(aes(x=xvar, y=est, color=method, group=method)) +
    geom_point(aes(x=xvar, y=est, color=method)) +
    geom_hline(aes(yintercept=truth)) 
  
  if (ci) {
    plot <- plot + 
      #geom_pointrange(aes(x=xvar, y=est, ymin=lb, ymax=ub, color=ibcf))
      geom_ribbon(aes(x=xvar, ymin=lb, ymax=ub, fill=method, group=method), alpha=.1)
  }
  
  plot + 
    facet_wrap(~metric, scales='free') +
    theme(axis.title.y = element_blank()) +
    labs(subtitle=estimand,
         x=xlabel,
         caption=cap)
}

fig1_contingency_wrapper <- function(overall, covarying='sig_u',  ..., const_cap=TRUE) {
  if (covarying =='sig_u') {
    figs <- expand_grid(estimand=c('SATT','CATT'),
                        sigu_const=c(.5,1,2)) %>%
      pmap(danfig1, overall=overall, ..., const_cap=const_cap)
    names(figs) <- outer(c('sig_u_0.5','sig_u_1','sig_u_2'), c('SATT', 'CATT'), paste, sep='_') %>% as.vector
  }
  if (covarying =='sig_v') {
    figs <- expand_grid(estimand=c('SATT','CATT'),
                        sigv_const=c(0,.5,1,2)) %>%
      pmap(danfig1, overall=overall, ..., const_cap=const_cap)
    names(figs) <- outer(c('sig_v_0','sig_v_0.5','sig_v_1','sig_v_2'), c('SATT', 'CATT'), paste, sep='_') %>% as.vector
  }
  if (covarying =='rho') {
    figs <- expand_grid(estimand=c('SATT','CATT'),
                        rho_const=c(-0.5,0,.5)) %>%
      pmap(danfig1, overall=overall, ..., const_cap=const_cap)
    names(figs) <- outer(c('rho_-0.5','rho_0','rho_0.5'), c('SATT', 'CATT'), paste, sep='_') %>% as.vector
  }
  if (covarying =='nT') {
    figs <- expand_grid(estimand=c('SATT','CATT'),
                        nT_const=c(500,1000)) %>%
      pmap(danfig1, overall=overall, ..., const_cap=const_cap)
    names(figs) <- outer(c('nT_500','nT_1000'), c('SATT', 'CATT'), paste, sep='_') %>% as.vector
  }
  
  cowplot::plot_grid(plotlist=lapply(figs, function(x) x + theme(legend.position='none')))
}

danfig2 <- function(mixing, varying='sig_v', sigv_const=1, sigu_const=1, rho_const=0, nT_const=1000, metric='cover', ci=FALSE) {
  touse <- mixing
  
  if (varying!='sig_v') {
    touse <- filter(touse, sigv_multiplier==sigv_const) 
  } else {
    xmult <- 8.33
    touse$xvar <- factor(round(xmult*touse$sigv_multiplier))
    xlabel <- expression(sigma[v])
  }
  if (varying!='sig_u') {
    touse <- filter(touse, sigu_multiplier==sigu_const) 
  } else {
    xmult <- 61.3
    touse$xvar <- factor(round(xmult*touse$sigu_multiplier))
    xlabel <- expression(sigma[u])
  }
  if (varying!='rho') {
    touse <- filter(touse, rho==rho_const) 
  } else {
    xmult <- 1
    touse$xvar <- factor(touse$rho)
    xlabel <- expression(rho)
  }
  if (varying!='nT') {
    touse <- filter(touse, nT==nT_const) 
  } else {
    xmult <- 1
    touse$xvar <- factor(touse$nT)
    xlabel <- '# treatment practices'
  }
  
  stopifnot(nrow(touse) == touse %>% select(method, xvar, sim, par) %>% distinct %>% nrow)
  
  smry <- touse %>%
    mutate(par = case_when(method=='wbcf' & par=='sigma' ~ 'sigma_y',
                           method=='wbcf' & par=='sigma2' ~ 'sigma2_y',
                           TRUE ~ par)) %>%
    filter(par %in% c('sigma_y','sigma_u','sigma_v','rho')) %>%
    group_by(par, method, xvar) %>%
    summarize(cover = mean(cover90),
              coverlb  = qbinom(.05, n(), cover)/n(),
              coverub  = qbinom(.95, n(), cover)/n(),
              rmse = sqrt(mean((real-mean)^2)),
              bias = mean(mean-real),
              est = mean(mean),
              estlb = quantile(mean,.05),
              estub = quantile(mean,.95),
              real = mean(real)) %>%
    mutate(truth=.9,
           grp = paste(par, method))
  if (metric=='cover') {
    plot <- ggplot(smry) + 
      geom_line(aes(x=xvar, y=cover, group=grp, color=method)) +
      geom_point(aes(x=xvar, y=cover, color=method)) +
      geom_hline(aes(yintercept=truth)) 
  } else if (metric=='rmse') {
    plot <- ggplot(smry) + 
      geom_line(aes(x=xvar, y=rmse, group=grp, color=method)) +
      geom_point(aes(x=xvar, y=rmse, color=method))
  } else if (metric=='bias') {
    plot <- ggplot(smry) + 
      geom_line(aes(x=xvar, y=bias, group=grp, color=method)) +
      geom_point(aes(x=xvar, y=bias, color=method))
  } else if (metric=='est') {
    plot <- ggplot(smry) + 
      geom_line(aes(x=xvar, y=est, group=grp, color=method)) +
      geom_point(aes(x=xvar, y=est, color=method)) + 
      geom_line(aes(x=xvar, y=real, group=par), linetype='dashed')
  }
  
  if (ci & metric=='cover') {
    plot <- plot + 
      #geom_pointrange(aes(x=xvar, y=est, ymin=lb, ymax=ub, color=ibcf))
      geom_ribbon(aes(x=xvar, ymin=coverlb, ymax=coverub, group=grp, fill=method), alpha=.1)
  } else if (ci & metric=='est') {
    plot <- plot + 
      #geom_pointrange(aes(x=xvar, y=est, ymin=lb, ymax=ub, color=ibcf))
      geom_ribbon(aes(x=xvar, ymin=estlb, ymax=estub, group=grp, fill=method), alpha=.1)
  }
  
  plot + 
    facet_wrap(~par, scales='free') +
    labs(x=xlabel)
}

danfig3 <- function(ex_smy, metric='tau', sigv_const=1, sigu_const=1, rho_const=0, nT_const=1000, show='diff') {
  touse <- ex_smy[[metric]]$confus %>%
    filter(sigv_multiplier   == sigv_const &
             sigu_multiplier == sigu_const &
             rho             == rho_const &
             nT              == nT_const &
             !(cutoff %in% c('top20pct','bot20pct'))) %>%
    select(sim, ibcf, predicted=cutoff, matches('top|bot|mid'), -top20, -bot20) %>%
    mutate(top1025 = top25-top10,
           top2550 = top50-top25,
           bot2550 = bot50-bot25,
           bot1025 = bot25-bot10) %>%
    pivot_longer(matches('top|bot|mid'), names_to='truth', values_to='pct')
  
  stopifnot(nrow(touse) == touse %>% select(ibcf, predicted, truth, sim) %>% distinct %>% nrow)
  
  smry <- touse %>%
    group_by(ibcf, predicted, truth) %>% 
    summarize(pct = mean(pct)) %>%
    ungroup 
  
  if (metric=='v') {
    smry <- filter(smry, ibcf)
  }
  smry <- smry %>%
    pivot_wider(c(predicted, truth), names_from=ibcf, values_from='pct') %>%
    rename(iBCF=`TRUE`)
  if (metric=='v') {
    #wbcf doesn't have v, but we can treat it as guessing at random
    smry <- mutate(smry, wBCF=case_when(truth %in% c('bot10', 'top10') ~ .1,
                                        truth %in% c('bot25', 'top25') ~ .25,
                                        truth %in% c('bot50', 'top50', 'mid50') ~ .5,
                                        truth %in% c('bot1025', 'top1025') ~ .15,
                                        truth %in% c('bot2550', 'top2550') ~ .25,))
  } else {
    smry <- rename(smry, wBCF=`FALSE`)
  }
  smry <- smry %>%
    mutate(diff = iBCF-wBCF,
           truth = factor(truth, levels=c('bot10', 'bot1025', 'bot25', 'bot2550', 'bot50', 'mid50', 'top50', 'top2550', 'top25', 'top1025', 'top10')),
           predicted = factor(str_remove(predicted,'pct'), levels=c('bot10', 'bot25', 'bot50', 'mid50', 'top50', 'top25', 'top10')),
           labelibcf = round(iBCF*100),
           labelwbcf = round(wBCF*100),
           labeldiff = round(diff*100))
  
  cumlset <- c('bot10', 'bot25', 'bot50', 'mid50', 'top50', 'top25', 'top10')
  exclset <- c('bot10', 'bot1025', 'bot2550', 'top2550', 'top1025', 'top10')
  
  #The truths are being in the actual top/bottom X%
  #The predictions are the 10% most likely to be in the predicted group
  #so predicted=top50 = the 10% of practices most likely to be in the top half
  #It's _not_ the 50% of practices with the highest taus.
  #Which is why we can have 40% of the predicted top50 group in the actual top10, since that's only .1*.4 = 4% of total practices listed as top50 and in top10
  #_not .5*.4 = 20% of total practices listed in top 50 and in top 10.
  
  if (show=='ibcf') {
    ggplot(smry %>% filter(truth %in% exclset)) + 
      geom_tile(aes(x=truth, y=predicted, fill=iBCF)) + 
      geom_text(aes(x=truth, y=predicted, label=labelibcf), color='white')  
  } else if (show=='wbcf') {
    ggplot(smry %>% filter(truth %in% exclset)) + 
      geom_tile(aes(x=truth, y=predicted, fill=wBCF)) + 
      geom_text(aes(x=truth, y=predicted, label=labelwbcf), color='white')
  } else if (show=='diff') {
    ggplot(smry %>% filter(truth %in% exclset)) + 
      geom_tile(aes(x=truth, y=predicted, fill=diff)) + 
      geom_text(aes(x=truth, y=predicted, label=labeldiff), color='white')    
  } else {
    stop('show what now?')
  }
}
