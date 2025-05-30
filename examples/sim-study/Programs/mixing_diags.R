source('programs/utils.R')
source('programs/load_data.R')

mlong <- mixing %>%
  mutate(par = ifelse(par=='sigma','sigma_y',par)) %>%
  pivot_longer(cols=c(-world, -r2uv_multiplier, -r2trt_multiplier, -set, -scenario, -sim, -method, -par), names_to='metric', values_to='value') %>%
  pivot_wider(id_cols=c(world, r2uv_multiplier, r2trt_multiplier, set, scenario, sim, par, metric), values_from=value, names_from=method)

msmry <- mlong %>%
  group_by(world, set, scenario, par, metric) %>%
  summarize_all(mean) %>%
  ungroup 

#Do we cover?
mixing %>%
  filter(method=='ibcf' & set=='main' & par %in% c('sigma_u','sigma_v','sigma_y','sigv_delta')) %>%
  ggplot() + 
  geom_point(aes(x=real, y=mean), alpha=.1) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~world + par, labeller=labeller(.multi_line = FALSE), scales='free')

mixing %>%
  filter(method=='ibcf' & set=='norho' & par %in% c('sigma_u','sigma_v','sigma_y','sigv_delta')) %>%
  ggplot() + 
  geom_point(aes(x=real, y=mean), alpha=.1) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~world + par, labeller=labeller(.multi_line = FALSE), scales='free')

mixing %>%
  filter(method=='ibcf' & set %in% c('main','norho') & par %in% c('sigma_u','sigma_v','sigma_y','rho','sigv_delta')) %>%
  group_by(par, world, set) %>%
  summarize(truth = mean(real),
            avg_est = mean(mean),
            avg_bias = mean(abs(real-mean)),
            avg_bias_pct = mean(abs(real-mean)/abs(real)),
            cover80 = mean(cover80),
            mean_sd = mean(sd),
            sd_mean = sd(mean)) %>%
  print(n=Inf)

mixing %>%
  filter(method=='ibcf' & set %in% c('main','norho') & par %in% c('sigma_u','sigma_v','sigma_y','sigv_delta')) %>%
  group_by(par, world, set, r2uv_multiplier) %>%
  summarize(truth = mean(real),
            avg_est = mean(mean),
            avg_bias = mean(abs(real-mean)/abs(real)),
            cover80 = mean(cover80),
            mean_sd = mean(sd),
            sd_mean = sd(mean)) %>%
  print(n=Inf)

ggplot(data=filter(msmry, metric=='n_eff' & par %in% c('tau_scale', 'real_mu_scale', 'sigma_y','tau_bar','mu_bar'))) + 
  geom_point(aes(x=wbcf, y=ibcf, color=set)) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~par, scales='free')

ggplot(data=filter(msmry, metric=='n_eff' & par %in% c('tau_scale', 'real_mu_scale', 'sigma_y','sigma_u','sigma_v','rho','tau_bar','mu_bar'))) + 
  geom_point(aes(x=ibcf, y=hbcf, color=set)) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~par, scales='free')

ggplot(data=filter(msmry, metric=='Rhat' & par %in% c('tau_scale', 'real_mu_scale', 'sigma_y','tau_bar','mu_bar'))) + 
  geom_point(aes(x=wbcf, y=ibcf, color=set)) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~par, scales='free')

ggplot(data=filter(mlong, metric=='n_eff' & par %in% c('tau_scale', 'real_mu_scale', 'sigma_y','tau_bar','mu_bar'))) + 
  geom_point(aes(x=wbcf, y=ibcf, color=set), alpha=0.1) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~par, scales='free')

ggplot(data=filter(mlong, metric=='n_eff' & par %in% c('tau_scale', 'real_mu_scale', 'sigma_y','tau_bar','mu_bar'))) + 
  geom_point(aes(x=wbcf, y=ibcf, color=world), alpha=0.1) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~par, scales='free')

#Why is neff so bad for medicare but not edu?
#Seems like it strongly correlates with r2trt
filter(mlong, par=='tau_bar' & metric=='n_eff') %>%
  mutate(abs_bad = ibcf<10,
         rel_bad = ibcf<wbcf) %>% 
  group_by(world, set, r2uv_multiplier, r2trt_multiplier) %>% 
  summarize(abs_bad = mean(abs_bad),
            rel_bad = mean(rel_bad)) %>% print(n=Inf)
#But only sometimes?
ggplot(data=filter(mlong, metric=='n_eff' & par=='tau_bar' & set %in% c('main','norho'))) + 
  geom_point(aes(x=wbcf, y=ibcf, color=set), alpha=0.1) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~world, scales='free')

ggplot(data=filter(mlong, metric=='n_eff' & par=='tau_bar' & set %in% c('main','norho'))) + 
  geom_point(aes(x=wbcf, y=ibcf, color=set), alpha=0.1) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~world, scales='free')

#distin of neff by world/setting
ggplot(data=filter(mixing, par=='tau_bar' & set %in% c('main','norho')) %>%
         mutate(sim_pars = as.factor(glue('r2uv: {r2uv_multiplier}, r2trt: {r2trt_multiplier}')))) + 
  geom_density(aes(x=n_eff, color=sim_pars)) + 
  facet_wrap(~method + world + set, scales='free_y', labeller=labeller(.multi_line = FALSE), ncol=2)

#distn of bad CATEs by world/setting
#NB: left good, right bad - opposite of above
ggplot(data=filter(overall, set %in% c('main','norho')) %>%
         mutate(sim_pars = as.factor(glue('r2uv: {r2uv_multiplier}, r2trt: {r2trt_multiplier}')))) + 
  geom_density(aes(x=neff_lt100, color=sim_pars)) + 
  facet_wrap(~method + world + set, scales='free_y', labeller=labeller(.multi_line = FALSE), ncol=2)

#And really bad CATEs
ggplot(data=filter(overall, set %in% c('main','norho')) %>%
         mutate(sim_pars = as.factor(glue('r2uv: {r2uv_multiplier}, r2trt: {r2trt_multiplier}')))) + 
  geom_density(aes(x=neff_lt10, color=sim_pars)) + 
  facet_wrap(~method + world + set, scales='free_y', labeller=labeller(.multi_line = FALSE), ncol=2)

#distin of rhats
ggplot(data=filter(mlong, metric=='Rhat' & par=='tau_bar' & set %in% c('main','norho'))) + 
  geom_point(aes(x=wbcf, y=ibcf, color=set), alpha=0.1) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~world, scales='free')

filter(mixing, par=='tau_bar' & set %in% c('main','norho') & method=='ibcf') %>%
  mutate(neff_bad = n_eff<100,
         rhat_bad = Rhat>1.1) %>% 
  group_by(world, set, method, r2trt_multiplier, r2uv_multiplier) %>% 
  summarize(neff_bad = mean(neff_bad),
            rhat_bad = mean(rhat_bad)) %>% print(n=Inf)

#Same true for individual taus?
overall %>% filter(set %in% c('main','norho')) %>%
  group_by(world, set, method) %>%
  summarize(across(c(neff_lt10, neff_lt100, neff_q10, neff_q25, neff_q50), mean))

lm(n_eff ~ set*r2trt_multiplier*r2uv_multiplier, data=filter(mixing, par=='tau_bar' & set %in% c('main','norho') & method=='ibcf' & world=='edu')) %>% summary
lm(n_eff ~ set*r2trt_multiplier*r2uv_multiplier, data=filter(mixing, par=='tau_bar' & set %in% c('main','norho') & method=='ibcf' & world=='medicare')) %>% summary

#Neff for main models, by world
ggplot(data=filter(mlong, set=='main' & metric=='n_eff' & par %in% c('tau_scale', 'real_mu_scale', 'sigma_y','tau_bar','mu_bar'))) + 
  geom_point(aes(x=wbcf, y=ibcf, color=world), alpha=0.1) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~par, scales='free')

ggplot(data=filter(mlong, metric=='n_eff' & par %in% c('tau_scale', 'real_mu_scale', 'sigma_y','tau_bar','mu_bar'))) + 
  geom_density(aes(x=wbcf, color=set, linetype='wbcf')) + 
  geom_density(aes(x=ibcf, color=set, linetype='ibcf')) + 
  facet_wrap(~par, scales='free')

#Coverage of terms
mixing %>%
  filter(set=='main' & !is.na(cover90)) %>%
  group_by(world, method, par) %>%
  summarize(cover90 = mean(cover90)) %>%
  ungroup %>%
  pivot_wider(c(world, par), names_from=method, values_from=cover90)

#density plot of sigmas
mixing %>%
  filter(set=='main' & method=='ibcf' & str_detect(par,'rho|sig')) %>%
  mutate(mults = glue('uv{r2uv_multiplier}_trt{r2trt_multiplier}')) %>%
  ggplot() + 
  geom_density(aes(x=mean, color=mults)) + 
  geom_vline(aes(xintercept=real, color=mults)) +
  facet_wrap(~ par + world, scales='free', ncol=2)

#How are we doing on sigv delta?
mixing %>% 
  filter(set=='main' & method=='ibcf' & par=='sigv_delta') %>%
  group_by(world, r2uv_multiplier, r2trt_multiplier) %>%
  summarize(truth = mean(real),
            est = mean(mean),
            rmse = sqrt(mean((real-mean)^2)),
            mae = mean(abs(real-mean)),
            avg_sd = mean(sd),
            sd_mean = sd(mean),
            cover90 = mean(cover90))

overall %>% 
  pivot_wider(id_cols=c(set, world, scenario, sim), names_from=method, values_from=timing) %>% 
  ggplot() + 
  geom_point(aes(x=wbcf, y=ibcf)) + geom_abline(slope=1, intercept=0)

time <- overall %>% 
  pivot_wider(id_cols=c(set, world, scenario, sim), names_from=method, values_from=timing) %>%
  group_by(scenario, set, world) %>%
  summarize_all(mean)

ggplot(time) + 
  geom_point(aes(x=wbcf, y=ibcf, color=set)) + 
  geom_abline(slope=1, intercept=0)

overall %>% 
  select(-SATT, -SATT_hat, -matches('^SATTp[0-9]')) %>%
  mutate(abs_bias = abs(SATTbias)) %>%
  pivot_longer(cols=c(-set, -scenario, -sim, -ibcf), names_to='metric',values_to='value') %>%
  pivot_wider(id_cols=c(set, scenario, sim, metric), names_from=ibcf, values_from=value) %>%
  rename(ibcf=`TRUE`, wbcf=`FALSE`) %>%
  group_by(set, metric) %>%
  summarize(ibcf=mean(ibcf),
            wbcf=mean(wbcf)) %>%
  ungroup %>%
  pivot_wider(id_cols=metric, values_from=c(ibcf, wbcf), names_from=set) %>%
  select(metric, matches('main_het'), matches('main_homog'), matches('no_res'), matches('no_weights'), matches('t_dist')) %>%
  print(n=Inf)

#How different are CATE widths?
indiv <- bind_rows(readRDS('Data/main_het-indiv.RDS'),
                   readRDS('Data/main_homog-indiv.RDS'),
                   readRDS('Data/no_weights-indiv.RDS'),
                   readRDS('Data/no_res-indiv.RDS'),
                   readRDS('Data/t_dist-indiv.RDS'))
saveRDS(indiv,'Data/all-indiv.RDS')

indiv %>% 
  select(set, scenario, ibcf, matches('width')) %>%
  group_by(set, scenario, ibcf) %>%
  summarize_all(mean)

avgs <- indiv %>% 
  group_by(scenario, ibcf, sim) %>% 
  summarize(width95 = mean(width95)) %>% 
  ungroup %>%
  mutate(type = str_extract(scenario,'het|homog'),
         sigv = str_match(scenario,'sigv([0-9.]+)')[,2],
         rho = str_match(scenario,'rho([0-9.-]+)')[,2],
         ibcf = as.numeric(ibcf))
lm(width95 ~ scenario + type:ibcf + sigv:ibcf + rho:ibcf, data=avgs) %>% summary

ex_smy$tau$auc %>%
  group_by(set, ibcf, cutoff) %>%
  summarize(auc = mean(auc)) %>%
  ungroup %>%
  pivot_wider(id_cols=c(set, cutoff), values_from=auc, names_from=ibcf) %>%
  rename(ibcf=`TRUE`,wbcf=`FALSE`) %>%
  filter(str_detect(cutoff,'top')) %>%
  arrange(cutoff, set)

#Some Nans e.g. if no one is misclassified
#Seems fair to ignore, but may mean measures are noisier than we think
ex_smy$tau$rmse %>%
  group_by(set, ibcf, cutoff, group) %>%
  summarize(rmse = mean(rmse, na.rm=TRUE),
            maqe = mean(maqe, na.rm=TRUE)) %>%
  ungroup %>%
  pivot_wider(id_cols=c(set, cutoff, group), values_from=c(rmse, maqe), names_from=ibcf) %>%
  filter(str_detect(cutoff,'top')) %>%
  arrange(cutoff, group, set) %>%
  print(n=Inf)

ex_smy$tau$confus %>%
  group_by(set, ibcf, cutoff) %>%
  select(-scenario, -par) %>%
  summarize_all(mean) %>%
  ungroup %>%
  filter(str_detect(cutoff,'top')) %>%
  arrange(cutoff, set, ibcf) %>%
  print(n=Inf)

#How does rho look?
mixing %>%
  left_join(scenarios %>% mutate(sig_v = paste0('sigv=',sig_v)) %>% select(scenario, te = trt_eff_scenario, sig_v, rho),
            by='scenario') %>%
  filter(str_detect(set,'main') & method %in% c('ibcf','hbcf') & par %in% c('rho','sigma_v','sigma_y','sigma_u','sigv_delta')) %>%
  group_by(scenario, te, sig_v, rho, method, par) %>%
  summarize(mean = mean(mean),
            truth = mean(real)) %>%
  ungroup %>%
  ggplot() + 
  geom_point(aes(x=rho, y=mean, color=te, shape=method)) + 
  geom_line(aes(x=rho, y=truth)) + 
  facet_wrap(~par + sig_v, labeller = labeller(.multi_line = FALSE), scales='free') + 
  labs(y='Mean posterior mean')

mixing %>%
  left_join(scenarios %>% mutate(sig_v = paste0('sigv=',sig_v)) %>% select(scenario, te = trt_eff_scenario, sig_v, rho),
            by='scenario') %>%
  filter(str_detect(set,'main') & method %in% c('ibcf','hbcf') & par %in% c('rho','sigma_v','sigma_y','sigma_u','sigv_delta')) %>%
  group_by(scenario, te, sig_v, rho, method, par) %>%
  summarize(mean = mean(mean),
            sd = mean(sd),
            truth = mean(real)) %>%
  ungroup %>%
  ggplot() + 
  geom_point(aes(x=rho, y=sd, color=te, shape=method)) + 
  facet_wrap(~par + sig_v, labeller = labeller(.multi_line = FALSE), scales='free') + 
  labs(y='Mean posterior sd')

mixing %>%
  left_join(scenarios %>% mutate(sig_v = paste0('sigv=',sig_v)) %>% select(scenario, te = trt_eff_scenario, sig_v, rho),
            by='scenario') %>%
  filter(str_detect(set,'main') & ibcf & par %in% c('rho','sigma_v','sigma_y','sigma_u','sigv_delta')) %>%
  group_by(scenario, te, sig_v, rho, par) %>%
  summarize(cover80 = mean(cover80),
            truth = mean(real)) %>%
  ungroup %>%
  ggplot() + 
  geom_point(aes(x=rho, y=cover80, color=te)) + 
  facet_wrap(~par + sig_v, labeller = labeller(.multi_line = FALSE))

mixing %>%
  left_join(scenarios %>% mutate(sig_v = paste0('sigv=',sig_v)) %>% select(scenario, te = trt_eff_scenario, sig_v, rho),
            by='scenario') %>%
  filter(str_detect(set,'main') & ibcf & par=='sigv_delta') %>%
  group_by(te, sig_v, rho) %>%
  select(te, sig_v, rho, real, mean, sd, matches('cover')) %>%
  summarize_all(mean) %>%
  ungroup %>%
  ggplot() + 
  geom_point(aes(x=rho, y=sd, color=te)) + 
  facet_wrap(~sig_v, labeller = labeller(.multi_line = FALSE))

overall %>% 
  filter(str_detect(set,'main')) %>%
  left_join(scenarios, by=c('scenario','set')) %>%
  group_by(method, trt_eff_scenario, sig_v, rho) %>%
  summarize(CATTwidth95 = mean(CATTwidth95)) %>%
  ungroup %>%
  ggplot() + 
  geom_point(aes(x=rho, y=CATTwidth95, color=method)) + 
  facet_wrap(~trt_eff_scenario + sig_v)

#What looks bad when tau_neff is bad?
mixing %>%
  filter(world=='medicare' & set %in% c('norho','main') & method=='ibcf') %>%
  group_by(world, set, scenario, sim) %>%
  mutate(tb_neff = max(case_when(par=='tau_bar' & n_eff<10 ~ 'tb_lt10',
                                 par=='tau_bar' & n_eff<100 ~ 'tb_lt100',
                                 par=='tau_bar' ~ 'tb_ge100'), na.rm = TRUE)) %>%
  group_by(world, par, tb_neff) %>%
  summarize(num = n(),
            x = mean(n_eff<10),
            b = mean(n_eff>=10 & n_eff<100),
            n = mean(n_eff)) %>%
  pivot_wider(id_cols=c(world, par), names_from=tb_neff, values_from=c(n, x, b)) %>%
  select(world, par, n_tb_lt10, n_tb_lt100, n_tb_ge100, x_tb_lt10, x_tb_lt100, x_tb_ge100, b_tb_lt10, b_tb_lt100, b_tb_ge100) %>%
  arrange(world, par) %>%
  print(n=Inf)

#Do we get worse performance when tau_bar is low?
overall %>%
  select(world, set, scenario, sim, method, SATT, SATT_hat, matches('SATTcover')) %>%
  inner_join(mixing %>% 
               filter(par=='tau_bar') %>% 
               mutate(tb_neff = case_when(par=='tau_bar' & n_eff< -10 ~ 1,
                                              par=='tau_bar' & n_eff<100 ~ 2,
                                              par=='tau_bar' ~ 3) %>%
                        factor(levels=1:3, labels=c('tb_lt10','tb_lt100','tb_ge100'))) %>%
               select(world, set, scenario, sim, method, tb_neff),
             by=c('world', 'set', 'scenario', 'sim', 'method')) %>%
  filter(world=='medicare' & set %in% c('main','norho') & method=='ibcf') %>%
  group_by(world, set, tb_neff) %>%
  summarize(n=n(),
            avg_bias = mean(SATT_hat - SATT),
            avb_abs_bias = mean(abs(SATT_hat - SATT)),
            rmse = sqrt(mean((SATT_hat - SATT)^2)),
            across(matches('cover'),mean))

tb_df <- overall %>%
  select(world, set, scenario, sim, method, SATT, SATT_hat, matches('SATTcover')) %>%
  inner_join(mixing %>% 
               filter(par=='tau_bar') %>% 
               mutate(tb_neff = n_eff) %>%
               select(world, set, scenario, sim, method, tb_neff),
             by=c('world', 'set', 'scenario', 'sim', 'method')) %>%
  filter(world=='medicare' & set %in% c('main','norho') & method=='ibcf')

ggplot(data=tb_df) + geom_point(aes(x=tb_neff, y=SATT_hat - SATT))

mixing %>% 
  filter(method=='ibcf' & set %in% c('main','norho') & par %in% c('tau_bar','mu_bar','tau_scale','real_mu_scale')) %>%
  ggplot() + 
  geom_density(aes(x=n_eff, color=world)) + 
  facet_wrap(~par, scales='free')

mixing %>% 
  filter(method=='ibcf' & set %in% c('main','norho') & par %in% c('tau_bar','mu_bar','tau_scale','real_mu_scale')) %>%
  pivot_wider(id_cols=c(world, set, scenario, sim), names_from=par, values_from=n_eff) %>%
  ggplot() +
  geom_point(aes(x=tau_scale, y=tau_bar)) + 
  facet_wrap(~world + set)
