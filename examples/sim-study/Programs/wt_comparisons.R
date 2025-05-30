source('programs/utils.R')
source('programs/load_data.R')

filter(scenarios, !weights)

wcomp_long <- overall %>%
  left_join(scenarios %>% select(scenario, set, sig_u, sig_v, rho, trt_eff_scenario, uv_dist, weights), by=c('scenario','set')) %>%
  filter(sig_u %in% c(0,1) & rho==0 & sig_v %in% c(0,1.5) & uv_dist=='normal' & method!='hbcf') %>%
  mutate(res = ifelse(sig_u==0,'no REs','has REs'),
         weights = ifelse(weights,'wtd','unwtd')) %>%
  group_by(method, trt_eff_scenario, res, weights) %>% 
  summarize(across(c(matches('width|cover'), PEHE, PEHT), mean),
            SATTrmse = sqrt(mean(SATTbias^2)),
            SATTmae = mean(abs(SATTbias))) %>%
  ungroup %>%
  pivot_longer(cols=c(-method, -trt_eff_scenario, -res, -weights), values_to='value',names_to='metric') %>%
  filter(!str_detect(metric,'CATU')) %>%
  mutate(combo = paste(trt_eff_scenario, res, sep=', '))

wcomp_wide <- wcomp_long %>%
  pivot_wider(id_cols=c(trt_eff_scenario, res, metric), values_from=value, names_from=c(weights, method), names_sep='_')

wcomp_semiwide <- wcomp_long %>%
  pivot_wider(id_cols=c(trt_eff_scenario, res, metric, combo, method), values_from=value, names_from=weights, names_sep='_')

ggplot(wcomp_semiwide) + 
  geom_point(aes(x=wtd, y=unwtd, color=method, shape=combo)) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~metric, scales='free')

#How does t(3) quantile compare?
rt3 <- rt(1e6,3) * sqrt(1/3)
mean(rt3 < qnorm(.9) & rt3 > qnorm(.1))
mean(rt3 < qnorm(.95) & rt3 > qnorm(.05))
mean(rt3 < qnorm(.975) & rt3 > qnorm(.025))

ggplot() + 
  geom_density(aes(x=rnorm(1e5),color='n')) +
  geom_density(aes(x=rt(1e5,3),color='rt3')) +
  geom_density(aes(x=rt(1e5,3)*sqrt(1/3),color='st3')) +
  scale_x_continuous(limits=c(-5,5))

wmix <- mixing %>%
  left_join(scenarios %>% select(scenario, set, sig_u, sig_v, rho, trt_eff_scenario, uv_dist, weights), by=c('scenario','set')) %>%
  filter(sig_u %in% c(0,1) & rho==0 & sig_v %in% c(0,1.5) & uv_dist=='normal' & method!='hbcf') %>%
  mutate(par = ifelse(par=='sigma','sigma_y',par),
         res = ifelse(sig_u==0,'no REs','has REs'),
         weights = ifelse(weights,'wtd','unwtd')) %>%
  group_by(method, trt_eff_scenario, res, weights, par) %>% 
  summarize(across(c(mean, sd, matches('width|cover'), n_eff), ~mean(.x)),
            rmse = sqrt(mean((mean-real)^2)),
            mae = mean(abs(mean-real))) %>%
  ungroup %>%
  pivot_longer(cols=c(-method, -trt_eff_scenario, -res, -weights, -par), values_to='value',names_to='metric') %>%
  mutate(combo = paste(trt_eff_scenario, res, sep=', ')) %>%
  pivot_wider(id_cols=c(trt_eff_scenario, res, par, combo, method, metric), values_from=value, names_from=weights, names_sep='_')

ggplot(wmix %>% filter(metric=='mean')) + 
  geom_point(aes(x=wtd, y=unwtd, color=method, shape=combo)) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~par, scales='free')

ggplot(wmix %>% filter(metric=='sd')) + 
  geom_point(aes(x=wtd, y=unwtd, color=method, shape=combo)) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~par, scales='free')

ggplot(wmix %>% filter(metric=='width95')) + 
  geom_point(aes(x=wtd, y=unwtd, color=method, shape=combo)) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~par, scales='free')

ggplot(wmix %>% filter(metric=='cover95')) + 
  geom_point(aes(x=wtd, y=unwtd, color=method, shape=combo)) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~par)

ggplot(wmix %>% filter(metric=='n_eff')) + 
  geom_point(aes(x=wtd, y=unwtd, color=method, shape=combo)) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~par, scales='free')
