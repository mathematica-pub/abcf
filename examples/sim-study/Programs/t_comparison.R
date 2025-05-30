source('programs/utils.R')
source('programs/load_data.R')

filter(scenarios, uv_dist=='t')

tcomp_long <- overall %>%
  left_join(scenarios %>% select(scenario, set, sig_v, rho, trt_eff_scenario, uv_dist), by=c('scenario','set')) %>%
  filter(set %in% c('main_het','main_homog','t_dist') & rho==0 & sig_v!=2 & method!='hbcf') %>%
  group_by(method, trt_eff_scenario, sig_v, uv_dist) %>% 
  summarize(across(c(matches('width|cover'), PEHE, PEHT), mean),
            SATTrmse = sqrt(mean(SATTbias^2)),
            SATTmae = mean(abs(SATTbias))) %>%
  ungroup %>%
  pivot_longer(cols=c(-method, -trt_eff_scenario, -sig_v, -uv_dist), values_to='value',names_to='metric') %>%
  filter(!str_detect(metric,'CATU')) %>%
  mutate(uv_dist = substr(uv_dist,1,1),
         combo = paste(uv_dist, method, sep='_'),
         sig_v = as.factor(sig_v))

tcomp_wide <- tcomp_long %>%
  pivot_wider(id_cols=c(trt_eff_scenario, sig_v, metric), values_from=value, names_from=c(method, uv_dist), names_sep='_')

ggplot(tcomp_long) + 
  geom_point(aes(x=sig_v, y=value, color=combo)) +
  geom_line(aes(x=sig_v, y=value, color=combo)) + 
  facet_wrap(~trt_eff_scenario + metric, scales='free', labeller=labeller(.multi_line = FALSE))

ggplot(tcomp_wide) + 
  geom_point(aes(x=ibcf_n, y=ibcf_t, color=trt_eff_scenario,shape=sig_v)) + 
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

tmix <- mixing %>%
  left_join(scenarios %>% select(scenario, set, sig_u, sig_v, rho, trt_eff_scenario, uv_dist, weights), by=c('scenario','set')) %>%
  filter(set %in% c('main_het','main_homog','t_dist') & rho==0 & sig_v!=2 & method!='hbcf') %>%
  mutate(par = ifelse(par=='sigma','sigma_y',par),
         uv_dist = substr(uv_dist,1,1),
         sig_v = as.factor(sig_v)) %>%
  group_by(method, trt_eff_scenario, uv_dist, sig_v, par) %>% 
  summarize(across(c(mean, sd, matches('width|cover'), n_eff), ~mean(.x)),
            rmse = sqrt(mean((mean-real)^2)),
            mae = mean(abs(mean-real))) %>%
  ungroup %>%
  pivot_longer(cols=c(-method, -trt_eff_scenario, -uv_dist, -sig_v, -par), values_to='value',names_to='metric') %>%
  pivot_wider(id_cols=c(trt_eff_scenario, sig_v, par, method, metric), values_from=value, names_from=uv_dist, names_sep='_')

ggplot(tmix %>% filter(metric=='mean')) + 
  geom_point(aes(x=n, y=t, color=method, shape=sig_v)) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~par, scales='free')

ggplot(tmix %>% filter(metric=='sd')) + 
  geom_point(aes(x=n, y=t, color=method, shape=sig_v)) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~par, scales='free')

ggplot(tmix %>% filter(metric=='cover95')) + 
  geom_point(aes(x=n, y=t, color=method, shape=sig_v)) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~par, scales='free')

ggplot(tmix %>% filter(metric=='n_eff')) + 
  geom_point(aes(x=n, y=t, color=method, shape=sig_v)) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~par, scales='free')
