source('programs/utils.R')
source('programs/load_data.R')

indiv <- readRDS('data/all-indiv.RDS') %>% select(set, scenario, sim, method, id, z, tau=truth, taux, v, tau_hat=mean, tau_sd=sd)
indiv_v <- readRDS('data/all-indiv_uv.RDS') %>% filter(par=='v') %>% select(set, scenario, sim, method, id, v_hat=mean, v_sd=sd)

indiv_taux <- left_join(indiv, indiv_v, by=c('set', 'scenario', 'sim', 'method', 'id')) %>%
  mutate(v_hat = ifelse(is.na(v_hat), 0, v_hat),
         taux_hat = tau_hat - v_hat)

sim_taux <- indiv_taux %>%
  filter(str_detect(set,'main')) %>%
  group_by(set, scenario, method, sim) %>%
  summarize(rmse_taux = sqrt(mean((taux_hat - taux)^2)),
            rmse_tau  = sqrt(mean((tau_hat  - tau)^2)),
            rmse_v    = sqrt(mean((v_hat    - v)^2))) %>%
  ungroup %>%
  left_join(scenarios %>% select(set, scenario, sig_v, rho, trt_eff_scenario),
            by=c('set','scenario'))

#Do I want mean of sim rmses, or rmse across all folks
scn_rmse <- sim_taux %>%
  group_by(set, scenario, method, sig_v, rho, trt_eff_scenario) %>%
  summarize_all(mean) %>%
  ungroup %>%
  pivot_longer(cols=matches('rmse'), names_to='par', values_to='rmse')

ggplot(scn_rmse) + 
  geom_point(aes(x=rho, y=rmse, color=method)) + 
  facet_wrap(~sig_v + trt_eff_scenario + par, labeller = labeller(.multi_line = FALSE), scales='free')

ggplot(scn_rmse %>% mutate(group = paste0(trt_eff_scenario, method))) + 
  geom_point(aes(x=rho, y=rmse, color=method, shape=trt_eff_scenario)) + 
  geom_line(aes(x=rho, y=rmse, color=method, group=group)) + 
  facet_wrap(~par + sig_v, labeller = labeller(.multi_line = FALSE), scales='free', nrow=3)

lm(taux_hat ~ taux + v, data=filter(indiv_taux, str_detect(set,'main') & method=='ibcf')) %>% summary
lm(taux_hat ~ taux + v, data=filter(indiv_taux, set=='main_het' & method=='ibcf')) %>% summary
lm(taux_hat ~ taux + v, data=filter(indiv_taux, set=='main_homog' & method=='ibcf')) %>% summary
