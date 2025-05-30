source('programs/utils.R')
overall <- readRDS('data/aug-2023-uold/all-overall.RDS')
mixing <- readRDS('data/aug-2023-uold/all-mixing.RDS')
ex_smy  <- readRDS('data/aug-2023-uold/all-exemplar-summy.RDS')

icc_overall <- readRDS('data/oct-2023-paper-icc/all-overall.RDS')
uhp_overall <- readRDS('data/oct-2023-paper-uhp/all-overall.RDS')
uhp_mix     <- readRDS('data/oct-2023-paper-uhp/all-mixing.RDS')
ibcf_overall <- readRDS('data/oct-2023-paper-ibcf/all-overall.RDS')

feb_overall <- readRDS('data/feb-2023-runs/all-overall.RDS')
feb_mixing  <- readRDS('data/feb-2023-runs/all-mixing.RDS')
mar_overall <- readRDS('data/mar-2023-uiw/all-overall.RDS')
mar_mixing  <- readRDS('data/mar-2023-uiw/all-mixing.RDS')

overall %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1) %>%
  group_by(model) %>%
  summarize(SATT_RMSE = sqrt(mean((SATT_hat-SATT)^2)),
            SATT_cover90 = mean(SATTcover90),
            SATT_width90 = mean(SATTwidth90),
            UTE_RMSE = sqrt(mean(PEHT^2)),
            UTE_cover90 = mean(CATTcover90),
            UTE_width90 = mean(CATTwidth90)) %>%
  pivot_longer(-model, names_sep='_', names_to=c('estimand','metric')) %>%
  mutate(nominal = ifelse(metric=='cover90',.9,NA_real_),
         metric = factor(metric, levels=c('RMSE','width90','cover90'), labels=c('RMSE','90% interval width','90% interval coverage')),
         model = factor(model, levels=c('oldBCF','uBCF'), labels=c('BCF','aBCF'))) %>%
  ggplot() + 
  geom_col(aes(x=model, y=value, fill=model)) + 
  facet_grid(metric ~ estimand, scales='free_y') +
  geom_hline(aes(yintercept=nominal), linetype='dashed') +
  labs(x='',
       y='') + 
  theme(legend.position='none')
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/BCF paper 2023/outline_figs/fig1_rmsesatt.png', height=4.5, width=8)

overall %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1) %>%
  mutate(SATT_absbias = abs(SATTbias),
         SATT_SqEr = SATTbias^2,
         UTE_SqEr = PEHT^2) %>%
  select(sim, model, SATT_absbias, SATT_SqEr, SATT_cover90=SATTcover90, SATT_width90=SATTwidth90, 
         UTE_RMSE=PEHT, UTE_SqEr, UTE_cover90=CATTcover90, UTE_width90=CATTwidth90) %>%
  pivot_longer(cols=-c(sim, model)) %>%
  nest(.by=name) %>%
  mutate(fit = lapply(data, lm, formula = value ~ 0 + as.factor(sim) + model),
         diff = lapply(fit, function(m) coefficients(m)['modeluBCF']),
         se   = lapply(fit, function(m) sqrt(diag(vcov(m)))['modeluBCF'])) %>%
  unnest(cols=c(data, diff, se)) %>%
  group_by(name, model, diff, se) %>%
  summarize(mean = mean(value)) %>%
  ungroup %>%
  pivot_wider(names_from=model, values_from=mean) %>%
  mutate(p = 2*(1-pnorm(abs(diff/se)))) %>%
  select(metric=name, oldBCF, uBCF, diff, se, p) %>%
  write_csv('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/BCF paper 2023/outline_figs/table1_rmsecvg.csv')

overall %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1) %>%
  mutate(SATTabsbias = abs(SATTbias)) %>%
  select(sim, model, SATTbias, SATTabsbias, SATTcover90, SATTwidth90, PEHT, CATTcover90, CATTwidth90) %>%
  pivot_longer(cols=-c(sim, model)) %>%
  split({.}$name) %>%
  lapply(function(x) {
    x %>% lm(value ~ model + as.factor(sim), data=.) %>%
      summary %>%
      `$`('coefficients') %>%
      as_tibble(rownames='par') %>%
      filter(str_detect(par,'uBCF')) %>%
      mutate(old_mean = mean(x %>% filter(model=='oldBCF') %>% pull(value)))
  }) %>%
  bind_rows(.id='measure')
fixest::feols(value ~ model:name | sim^name) %>%
  summary

overall %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1) %>%
  group_by(model) %>%
  summarize(SATT_RMSE = sqrt(mean((SATT_hat-SATT)^2)),
            SATT_cover90 = mean(SATTcover90),
            SATT_width90 = mean(SATTwidth90),
            UTE_RMSE = sqrt(mean(PEHT^2)),
            UTE_cover90 = mean(CATTcover90),
            UTE_width90 = mean(CATTwidth90),
            resid_RMSE=sqrt(mean(residT_RMSET^2)),
            resid_cover90=mean(residTcover90),
            resid_width90=mean(residTwidth90)) %>%
  pivot_longer(-model, names_sep='_', names_to=c('estimand','metric')) %>%
  mutate(nominal = ifelse(metric=='cover90',.9,NA_real_)) %>%
  ggplot() + 
  geom_col(aes(x=model, y=value, fill=model)) + 
  facet_grid(metric ~ estimand, scales='free_y') +
  geom_hline(aes(yintercept=nominal), linetype='dashed') +
  labs(title='fig 1alt')

overall %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1) %>%
  group_by(model) %>%
  summarize(SATT_RMSE_est = sqrt(mean((SATT_hat-SATT)^2)),
            SATT_cover90_est = mean(SATTcover90),
            SATT_cover90_lb  = qbinom(p=.05, prob=SATT_cover90_est,size=n())/n(),
            SATT_cover90_ub  = qbinom(p=.95, prob=SATT_cover90_est,size=n())/n(),
            SATT_width90_est = mean(SATTwidth90),
            SATT_width90_lb  = quantile(SATTwidth90,.05),
            SATT_width90_ub  = quantile(SATTwidth90,.95),
            UTE_RMSE_est = sqrt(mean(PEHT^2)),
            UTE_RMSE_lb  = quantile(PEHT,.05),
            UTE_RMSE_ub  = quantile(PEHT,.95),
            UTE_cover90_est = mean(CATTcover90),
            UTE_cover90_lb  = quantile(CATTcover90,.05),
            UTE_cover90_ub  = quantile(CATTcover90,.95),
            UTE_width90_est = mean(CATTwidth90),
            UTE_width90_lb  = quantile(CATTwidth90,.05),
            UTE_width90_ub  = quantile(CATTwidth90,.95),
            resid_RMSE_est  = sqrt(mean(residT_RMSET^2)),
            resid_RMSE_lb   = quantile(residT_RMSET,.05),
            resid_RMSE_ub   = quantile(residT_RMSET,.95),
            resid_cover90_est = mean(residTcover90),
            resid_cover90_lb  = quantile(residTcover90, .05),
            resid_cover90_ub  = quantile(residTcover90, .95),
            resid_width90_est = mean(residTwidth90),
            resid_width90_lb  = quantile(residTwidth90,.05),
            resid_width90_ub  = quantile(residTwidth90,.95)) %>%
  pivot_longer(-model, names_sep='_', names_to=c('estimand','metric','.value')) %>%
  mutate(nominal = ifelse(metric=='cover90',.9,NA_real_)) %>%
  ggplot() + 
  geom_col(aes(x=model, y=est, fill=model)) + 
  geom_errorbar(aes(x=model, ymin=lb, ymax=ub, fill=model)) +
  facet_grid(metric ~ estimand, scales='free_y') +
  geom_hline(aes(yintercept=nominal), linetype='dashed') +
  labs(title='fig 1alt')

#but those CIs ignore covariance
wide_overall <- overall %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1) %>%
  select(sim, model, SATTbias, CATT_rmse=PEHT, residt_RMSE=residT_RMSET, matches('^(SATT|CATT|residT)(cover|width)90')) %>%
  pivot_longer(-c(sim, model)) %>%
  pivot_wider(c(sim, name), names_from=model, values_from=value)

ggplot(wide_overall) + 
  geom_point(aes(x=oldBCF, y=uBCF)) + 
  geom_abline() + 
  facet_wrap(~name, scales='free')

wide_overall %>%
  pivot_longer(c(oldBCF, uBCF), names_to='model') %>%
  split({.}$name) %>%
  lapply(function(x){
    as_tibble(summary(fixest::feols(value ~ model | sim, data=x))$coeftable)
  }) %>%
  bind_rows(.id='metric') %>%
  mutate(star = case_when(`Pr(>|t|)`<=.01 ~ '***',
                          `Pr(>|t|)`<=.05 ~ '**',
                          `Pr(>|t|)`<=.1 ~ '*',
                          TRUE ~ ''))

wide_overall %>% mutate(d=uBCF-oldBCF) %>%
  group_by(name) %>%
  summarize(mean = mean(d),
            pos  = mean(d>0),
            sd   = sd(d),
            lb   = quantile(d,.05),
            ub   = quantile(d,.95))

overall %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1) %>%
  group_by(model) %>%
  summarize(resid_RMSE=sqrt(mean(residT_RMSET^2)),
            resid_cover90 = mean(residTcover90),
            resid_width90 = mean(residTwidth90),) %>%
  pivot_longer(-model) %>%
  mutate(nominal = ifelse(name=='resid_cover90',.9,NA_real_),
         metric = factor(name, 
                         levels=c('resid_RMSE','resid_width90','resid_cover90'),
                         labels=c('RMSE','90% interval width','90% interval coverage')),
         model = factor(model, levels=c('oldBCF','uBCF'), labels=c('BCF','aBCF'))) %>%
  ggplot() + 
  geom_col(aes(x=model, y=value, fill=model)) + 
  facet_wrap(~metric, scales='free_y') +
  geom_hline(aes(yintercept=nominal), linetype='dashed') +
  theme_bw() +
  theme(legend.position='none',
        strip.background =element_rect(fill="white")) +
  labs(x='',
       y='')
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/BCF paper 2023/outline_figs/fig2_rmseresid.png', height=4.5, width=8)

ex_smy$resid$confus %>% 
  filter(sigu_multiplier==1 & sigv_multiplier==1 & cutoff=='gt_45' & size=='all') %>%
  group_by(model) %>%
  summarize(value = mean(gt_45)) %>%
  ggplot() + 
  geom_col(aes(x=model, y=value, fill=model)) + 
  labs(title='fig 3', y='Share of estimated exemplars that are true exemplars')
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/BCF paper 2023/outline_figs/fig3_exemplarid.png')

overall %>%
  mutate(sigu_hyperprior=ifelse(model=='oldBCF','BCF','0.67')) %>%
  bind_rows(uhp_overall %>% mutate(sigu_hyperprior = as.character(round(sigu_hyperprior,2)))) %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1) %>%
  group_by(sigu_hyperprior, model) %>%
  summarize(SATT_RMSE = sqrt(mean((SATT_hat-SATT)^2)),
            SATT_cover90 = mean(SATTcover90),
            SATT_width90 = mean(SATTwidth90),
            UTE_RMSE = sqrt(mean(PEHT^2)),
            UTE_cover90 = mean(CATTcover90),
            UTE_width90 = mean(CATTwidth90),
            resid_RMSE=sqrt(mean(residT_RMSET^2)),
            resid_cover90=mean(residTcover90),
            resid_width90=mean(residTwidth90)) %>%
  pivot_longer(-c(sigu_hyperprior, model), names_sep='_', names_to=c('estimand','metric')) %>%
  mutate(nominal = ifelse(metric=='cover90',.9,NA_real_),
         model = factor(model, levels=c('oldBCF','uBCF'), labels=c('BCF','aBCF'))) %>%
  filter(estimand!='resid') %>%
  ggplot() + 
  geom_col(aes(x=sigu_hyperprior, y=value, fill=model)) + 
  facet_grid(metric ~ estimand, scales='free_y') +
  geom_hline(aes(yintercept=nominal), linetype='dashed') +
  labs(title='I DON\'T PLAN TO ACTUALLY SHOW THIS',
       x=bquote('Hyperprior for '*sigma[u]*', in terms of '*sigma[y]),
       y='',
       fill='')
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/BCF paper 2023/outline_figs/figuhp_sens.png', height=4.5, width=8)

mixing %>%
  mutate(sigu_hyperprior=ifelse(model=='oldBCF',model,'0.67')) %>%
  bind_rows(uhp_mix %>% mutate(sigu_hyperprior = as.character(round(sigu_hyperprior,2)))) %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1 & par %in% c('sigma_y','sigma_u') & method=='ubcf') %>%
  group_by(sigu_hyperprior, par) %>%
  summarize(rmse = sqrt(mean((mean-real)^2)),
            width = mean(width90),
            cover = mean(cover90)) %>%
  pivot_longer(c(rmse, width, cover), names_to='metric') %>%
  mutate(nominal = ifelse(metric=='cover',.9,NA_real_)) %>%
  ggplot() + 
  geom_col(aes(x=sigu_hyperprior, y=value)) +
  geom_hline(aes(yintercept=nominal), linetype='dashed') +
  facet_grid(metric ~ par, scales='free_y')

mixing %>%
  mutate(sigu_hyperprior=ifelse(model=='oldBCF',model,'0.67')) %>%
  bind_rows(uhp_mix %>% mutate(sigu_hyperprior = as.character(round(sigu_hyperprior,2)))) %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1 & par %in% c('sigma_y','sigma_u') & method=='ubcf') %>%
  mutate(SqEr = (mean-real)^2,
         sigu_hyperprior = factor(sigu_hyperprior, levels=c('0.67','0.17','0.33','1.33','2.67'))) %>%
  select(sim, sigu_hyperprior, par, SqEr, cover90, width90) %>%
  pivot_longer(cols=-c(sim, sigu_hyperprior, par)) %>%
  nest(.by=c(par, name)) %>%
  mutate(fit = lapply(data, lm, formula = value ~ 0 + as.factor(sim) + sigu_hyperprior),
         diff17 = lapply(fit, function(m) coefficients(m)['sigu_hyperprior0.17']),
         se17   = lapply(fit, function(m) sqrt(diag(vcov(m)))['sigu_hyperprior0.17']),
         diff33 = lapply(fit, function(m) coefficients(m)['sigu_hyperprior0.33']),
         se33   = lapply(fit, function(m) sqrt(diag(vcov(m)))['sigu_hyperprior0.33']),
         diff133 = lapply(fit, function(m) coefficients(m)['sigu_hyperprior1.33']),
         se133   = lapply(fit, function(m) sqrt(diag(vcov(m)))['sigu_hyperprior1.33']),
         diff267 = lapply(fit, function(m) coefficients(m)['sigu_hyperprior2.67']),
         se267   = lapply(fit, function(m) sqrt(diag(vcov(m)))['sigu_hyperprior2.67'])) %>%
  unnest(cols=c(data, matches('diff|se'))) %>%
  mutate(diff = case_when(sigu_hyperprior=='0.17' ~ diff17,
                          sigu_hyperprior=='0.33' ~ diff33,
                          sigu_hyperprior=='1.33' ~ diff133,
                          sigu_hyperprior=='2.67' ~ diff267),
         se = case_when(sigu_hyperprior=='0.17' ~ se17,
                        sigu_hyperprior=='0.33' ~ se33,
                        sigu_hyperprior=='1.33' ~ se133,
                        sigu_hyperprior=='2.67' ~ se267)) %>%
  group_by(name, par, sigu_hyperprior, diff, se) %>%
  summarize(mean = mean(value)) %>%
  mutate(p = 2*(1-pnorm(abs(diff/se))),
         sig = ifelse(p<.05,'*','')) %>%
  ungroup %>%
  select(-se, -sig) %>%
  pivot_wider(names_from=sigu_hyperprior, values_from=c(mean, diff, p)) %>%
  select(-diff_0.67, -p_0.67) %>%
  rename_all(str_replace,'mean','m') %>%
  rename_all(str_replace,'diff','d') %>%
  rename_all(str_replace,'sig','s') %>%
  arrange(par, name)

mixing %>%
  mutate(sigu_hyperprior=ifelse(model=='oldBCF',model,'0.67')) %>%
  bind_rows(uhp_mix %>% mutate(sigu_hyperprior = as.character(round(sigu_hyperprior,2)))) %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1 & par %in% c('sigma_y','sigma_u') & method=='ubcf') %>%
  group_by(par, sigu_hyperprior) %>%
  summarize(rmse = sqrt(mean((mean-real)^2)),
            width = mean(width90),
            cover=mean(cover90))

overall %>%
  mutate(sigu_hyperprior=ifelse(model=='oldBCF',model,'0.67')) %>%
  bind_rows(uhp_overall %>% mutate(sigu_hyperprior = as.character(round(sigu_hyperprior,2)))) %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1 & method=='ubcf') %>%
  mutate(sigu_hyperprior = factor(sigu_hyperprior, levels=c('0.67','0.17','0.33','1.33','2.67')),
         SATT_absbias = abs(SATTbias),
         SATT_SqEr = SATTbias^2,
         UTE_SqEr = PEHT^2) %>%
  select(sim, sigu_hyperprior, SATT_absbias, SATT_SqEr, SATT_cover90=SATTcover90, SATT_width90=SATTwidth90, 
         UTE_RMSE=PEHT, UTE_SqEr, UTE_cover90=CATTcover90, UTE_width90=CATTwidth90) %>%
  pivot_longer(cols=-c(sim, sigu_hyperprior)) %>%
  nest(.by=name) %>%
  mutate(fit = lapply(data, lm, formula = value ~ 0 + as.factor(sim) + sigu_hyperprior),
         diff17 = lapply(fit, function(m) coefficients(m)['sigu_hyperprior0.17']),
         se17   = lapply(fit, function(m) sqrt(diag(vcov(m)))['sigu_hyperprior0.17']),
         diff33 = lapply(fit, function(m) coefficients(m)['sigu_hyperprior0.33']),
         se33   = lapply(fit, function(m) sqrt(diag(vcov(m)))['sigu_hyperprior0.33']),
         diff133 = lapply(fit, function(m) coefficients(m)['sigu_hyperprior1.33']),
         se133   = lapply(fit, function(m) sqrt(diag(vcov(m)))['sigu_hyperprior1.33']),
         diff267 = lapply(fit, function(m) coefficients(m)['sigu_hyperprior2.67']),
         se267   = lapply(fit, function(m) sqrt(diag(vcov(m)))['sigu_hyperprior2.67'])) %>%
  unnest(cols=c(data, matches('diff|se'))) %>%
  mutate(diff = case_when(sigu_hyperprior=='0.17' ~ diff17,
                          sigu_hyperprior=='0.33' ~ diff33,
                          sigu_hyperprior=='1.33' ~ diff133,
                          sigu_hyperprior=='2.67' ~ diff267),
         se = case_when(sigu_hyperprior=='0.17' ~ se17,
                        sigu_hyperprior=='0.33' ~ se33,
                        sigu_hyperprior=='1.33' ~ se133,
                        sigu_hyperprior=='2.67' ~ se267)) %>%
  group_by(name, sigu_hyperprior, diff, se) %>%
  summarize(mean = mean(value)) %>%
  mutate(p = 2*(1-pnorm(abs(diff/se))),
         sig = ifelse(p<.05,'*','')) %>%
  ungroup %>%
  select(-se, -sig) %>%
  pivot_wider(names_from=sigu_hyperprior, values_from=c(mean, diff, p)) %>%
  select(-diff_0.67, -p_0.67) %>%
  rename_all(str_replace,'mean','m') %>%
  rename_all(str_replace,'diff','d') %>%
  rename_all(str_replace,'sig','s')

#Bah these don't actually have the same seeds across hyperpriors. Damn
wide_sens <- overall %>%
  mutate(sigu_hyperprior=ifelse(model=='oldBCF',model,'0.67')) %>%
  bind_rows(uprior_sens %>% mutate(sigu_hyperprior = as.character(round(sigu_hyperprior,2)))) %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1) %>%
  select(sim, model, sigu_hyperprior, SATTbias, CATT_rmse=PEHT, residt_RMSE=residT_RMSET, matches('^(SATT|CATT|residT)(cover|width)90')) %>%
  pivot_longer(-c(sim, model, sigu_hyperprior)) %>%
  pivot_wider(c(sim, name), names_from=sigu_hyperprior, values_from=value)

wide_sens %>%
  pivot_longer(-c(sim, name), names_to='model') %>%
  split({.}$name) %>%
  lapply(function(x){
    x %>%
      mutate(model = factor(model, levels=c('oldBCF','0.17','0.33','0.67','1.33','2.67'))) %>%
      fixest::feols(value ~ model | sim, data=.) %>%
      summary() %>%
      `$`('coeftable') %>%
      as_tibble(rownames='par')
  }) %>%
  bind_rows(.id='metric') %>%
  mutate(star = case_when(`Pr(>|t|)`<=.01 ~ '***',
                          `Pr(>|t|)`<=.05 ~ '**',
                          `Pr(>|t|)`<=.1 ~ '*',
                          TRUE ~ '')) %>%
  print(n=Inf)


wide_sens %>%
  pivot_longer(-c(sim, name), names_to='model') %>%
  split({.}$name) %>%
  lapply(function(x){
    x %>%
      mutate(model = factor(model, levels=c('0.67','oldBCF','0.17','0.33','1.33','2.67'))) %>%
      fixest::feols(value ~ model | sim, data=.) %>%
      summary() %>%
      `$`('coeftable') %>%
      as_tibble(rownames='par')
  }) %>%
  bind_rows(.id='metric') %>%
  mutate(star = case_when(`Pr(>|t|)`<=.01 ~ '***',
                          `Pr(>|t|)`<=.05 ~ '**',
                          `Pr(>|t|)`<=.1 ~ '*',
                          TRUE ~ '')) %>%
  print(n=Inf)


#ICC figs
r2es <- seq(0,1,.1)
iccs <- c(0,.0002,.0004,.0007,.0011,.0016,.0024,.0038,.0065,.0146,1)
options(scipen=999)
lbls <- paste0(round(r2es,2), '\n(', iccs, ')')

icc_overall %>%
  group_by(model, r2_error) %>%
  summarize(SATT_RMSE = sqrt(mean((SATT_hat-SATT)^2)),
            SATT_cover90 = 100*mean(SATTcover90),
            SATT_width90 = mean(SATTwidth90),
            UTE_RMSE = sqrt(mean(PEHT^2)),
            UTE_cover90 = 100*mean(CATTcover90),
            UTE_width90 = mean(CATTwidth90)) %>%
  pivot_longer(-c(model, r2_error), names_sep='_', names_to=c('estimand','metric')) %>%
  mutate(nominal = ifelse(metric=='cover90',90,NA_real_),
         metric = factor(metric, levels=c('RMSE','width90','cover90'), labels=c('RMSE','90% interval width','90% interval coverage')),
         model = factor(model, levels=c('oldBCF','uBCF'), labels=c('BCF','aBCF'))) %>%
  mutate(min = 0, bad=NA_real_) %>%
  ggplot() + 
  geom_line(aes(x=r2_error, y=value, color=model)) + 
  geom_point(aes(x=r2_error, y=value, color=model)) + 
  geom_point(aes(x=bad, y=min)) + 
  facet_wrap(~ metric + estimand, scales='free_y', ncol=2, labeller=labeller(.multi_line = FALSE)) +
  geom_hline(aes(yintercept=nominal), linetype='dashed') +
  scale_x_continuous(breaks=r2es, labels=r2es) +
  ggh4x::facetted_pos_scales(y=list(NULL, scale_y_continuous(limits=c(10,20)),
                                    scale_y_continuous(limits=c(10,25)), scale_y_continuous(limits=c(30,50)),
                                    scale_y_continuous(limits=c(70,100)), scale_y_continuous(limits=c(70,100)))) +
  theme_bw() +
  theme(legend.position='bottom',
        strip.background =element_rect(fill="white")) +
  labs(x='Residual variance share',
       y='',
       color='')
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/BCF paper 2023/outline_figs/figicc_rmsesatt.png', width=6.5, height=5.5)

icc_overall %>%
  mutate(SATT_SqEr = SATTbias^2,
         UTE_SqEr = PEHT^2) %>%
  select(sim, r2_error, model, SATT_SqEr, SATT_cover90=SATTcover90, SATT_width90=SATTwidth90, 
         UTE_RMSE=PEHT, UTE_SqEr, UTE_cover90=CATTcover90, UTE_width90=CATTwidth90) %>%
  pivot_longer(cols=-c(sim, r2_error, model)) %>%
  nest(.by=c(name,r2_error)) %>%
  mutate(fit = lapply(data, lm, formula = value ~ 0 + as.factor(sim) + model),
         diff = lapply(fit, function(m) coefficients(m)['modeluBCF']),
         se   = lapply(fit, function(m) sqrt(diag(vcov(m)))['modeluBCF'])) %>%
  unnest(cols=c(data, matches('diff|se'))) %>%
  select(name, r2_error, diff, se) %>%
  distinct() %>%
  mutate(sig = abs(diff/se) > qnorm(.975)) %>%
  arrange(name, r2_error) %>%
  print(n=Inf)

#compare iBCF and uBCF
overall %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1 & model=='uBCF') %>%
  bind_rows(ibcf_overall) %>%
  mutate(SATT_absbias = abs(SATTbias),
         SATT_SqEr = SATTbias^2,
         UTE_SqEr = PEHT^2,
         model=factor(model, levels=c('uBCF','iBCF'))) %>%
  select(sim, model, SATT_absbias, SATT_SqEr, SATT_cover90=SATTcover90, SATT_width90=SATTwidth90, 
         UTE_RMSE=PEHT, UTE_SqEr, UTE_cover90=CATTcover90, UTE_width90=CATTwidth90) %>%
  pivot_longer(cols=-c(sim, model)) %>%
  nest(.by=name) %>%
  mutate(fit = lapply(data, lm, formula = value ~ 0 + as.factor(sim) + model),
         diff = lapply(fit, function(m) coefficients(m)['modeliBCF']),
         se   = lapply(fit, function(m) sqrt(diag(vcov(m)))['modeliBCF'])) %>%
  unnest(cols=c(data, diff, se)) %>%
  group_by(name, model, diff, se) %>%
  summarize(mean = mean(value)) %>%
  ungroup %>%
  pivot_wider(names_from=model, values_from=mean) %>%
  mutate(p = 2*(1-pnorm(abs(diff/se)))) %>%
  select(metric=name, uBCF, iBCF, diff, se, p)

overall %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1) %>%
  bind_rows(ibcf_overall) %>%
  group_by(model) %>%
  summarize(SATT_RMSE = sqrt(mean((SATT_hat-SATT)^2)),
            SATT_cover90 = mean(SATTcover90),
            SATT_width90 = mean(SATTwidth90),
            UTE_RMSE = sqrt(mean(PEHT^2)),
            UTE_cover90 = mean(CATTcover90),
            UTE_width90 = mean(CATTwidth90)) %>%
  pivot_longer(-model, names_sep='_', names_to=c('estimand','metric')) %>%
  mutate(nominal = ifelse(metric=='cover90',.9,NA_real_),
         metric = factor(metric, levels=c('RMSE','width90','cover90'), labels=c('RMSE','90% interval width','90% interval coverage')),
         model = factor(model, levels=c('oldBCF','uBCF','iBCF'), labels=c('BCF','aBCF','iBCF'))) %>%
  ggplot() + 
  geom_col(aes(x=model, y=value, fill=model)) + 
  facet_grid(metric ~ estimand, scales='free_y') +
  geom_hline(aes(yintercept=nominal), linetype='dashed') +
  scale_fill_manual(values=c(BCF=scales::hue_pal()(2)[1], aBCF=scales::hue_pal()(2)[2], iBCF=pal$purple)) +
  labs(x='',
       y='') + 
  theme(legend.position='none')
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/BCF paper 2023/outline_figs/figapdx_rmsesatt.png', height=4.5, width=8)

#Need a comparison of sigma v doing worse over values
sigv_est <- mar_mixing %>%
  filter(par=='sigma_v' & sigu_multiplier==1 & method=='ibcf') %>%
  group_by(sigv_multiplier) %>%
  summarize(real= mean(real),
            est = mean(mean),
            lb  = quantile(mean,.05),
            ub  = quantile(mean,.95)) %>%
  ungroup %>%
  mutate(sigv_multiplier = factor(round(8.33*sigv_multiplier)))

sigvplot <- ggplot(sigv_est, aes(x=real)) + 
  geom_line(aes(y=real), linetype='dashed') +
  geom_line(aes(y=est)) +
  geom_ribbon(aes(ymin=lb, ymax=ub), alpha=.1) +
  geom_point(aes(y=est)) +
  labs(y=bquote('Estimate of'~sigma[v]),
       x=bquote('True'~sigma[v])) +
  scale_x_continuous(breaks=sort(unique(sigv_est$real)), labels=as.numeric(levels(sigv_est$sigv_multiplier))) + 
  scale_y_continuous(breaks=sort(unique(sigv_est$real)), labels=as.numeric(levels(sigv_est$sigv_multiplier))) + 
  theme(panel.grid.minor = element_blank()) + 
  theme_bw()
print(sigvplot)
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/BCF paper 2023/outline_figs/sigv_ests.png', height=4.5, width=8)

rho_est <- feb_mixing %>%
  filter(par=='rho' & sigu_multiplier==1 & sigv_multiplier==1 & method=='ibcf') %>%
  group_by(rho) %>%
  summarize(real= mean(real),
            est = mean(mean),
            lb  = quantile(mean,.05),
            ub  = quantile(mean,.95)) %>%
  ungroup

rhoplot <- ggplot(rho_est, aes(x=real)) + 
  geom_line(aes(y=real), linetype='dashed') +
  geom_line(aes(y=est)) +
  geom_ribbon(aes(ymin=lb, ymax=ub), alpha=.1) +
  geom_point(aes(y=est)) +
  labs(y=bquote('Estimate of'~rho),
       x=bquote('True'~rho)) +
  scale_x_continuous(breaks=sort(unique(rho_est$real))) + 
  theme(panel.grid.minor = element_blank()) +
  theme_bw()
print(rhoplot)
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/BCF paper 2023/outline_figs/rho_ests.png', height=4.5, width=8)

cowplot::plot_grid(plotlist=list(sigvplot, rhoplot))
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/BCF paper 2023/outline_figs/sigvrho_ests.png', height=4.5, width=8)


