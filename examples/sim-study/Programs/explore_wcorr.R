source('programs/utils.R')

overall <- readRDS('data/mar-2023-wcorr/all-overall.RDS')

smry <- overall %>%
  group_by(model, sigu_multiplier, r2_w_pi, r2_w_mu, r2_w_tau) %>%
  summarize(est_RMSE = sqrt(mean(SATTbias^2)),
            est_cover = mean(SATTcover90),
            est_width = mean(SATTwidth90),
            lb_width = quantile(SATTwidth90,.05),
            ub_width = quantile(SATTwidth90,.95)) %>%
  ungroup %>%
  mutate(lb_cover = qbinom(.05,50,est_cover)/50,
         ub_cover = qbinom(.95,50,est_cover)/50) %>%
  pivot_longer(matches('est|lb|ub'), names_sep='_', names_to=c('.value', 'metric')) %>%
  mutate(pi_class = case_when(r2_w_pi==0   ~ 'none',
                              r2_w_pi==.25 ~ 'low pi',
                              r2_w_pi==.5  ~ 'high pi') %>%
           factor(levels=c('none','low pi','high pi')),
         mu_class = case_when(r2_w_mu==0   ~ 'none',
                              r2_w_mu==.25 ~ 'low mu',
                              r2_w_mu==.5  ~ 'high mu') %>%
           factor(levels=c('none','low mu','high mu')),
         tau_class = case_when(r2_w_tau==0   ~ 'none',
                               r2_w_tau==.25 ~ 'low tau',
                               r2_w_tau==.5  ~ 'high tau') %>%
           factor(levels=c('none','low tau','high tau')),
         class = case_when(r2_w_pi==0   & r2_w_mu==0   & r2_w_tau==0 ~ 'none',
                           r2_w_pi==0   & r2_w_mu==0   & r2_w_tau!=0 ~ as.character(tau_class),
                           r2_w_pi==0   & r2_w_mu!=0   & r2_w_tau==0 ~ as.character(mu_class),
                           r2_w_pi!=0   & r2_w_mu==0   & r2_w_tau==0 ~ as.character(pi_class),
                           r2_w_pi==0   & r2_w_mu!=0   & r2_w_tau!=0 ~ paste(as.character(mu_class), as.character(tau_class), sep=', '),
                           r2_w_pi!=0   & r2_w_mu==0   & r2_w_tau!=0 ~ paste(as.character(pi_class), as.character(tau_class), sep=', '),
                           r2_w_pi!=0   & r2_w_mu!=0   & r2_w_tau==0 ~ paste(as.character(pi_class), as.character(mu_class), sep=', '),
                           r2_w_pi!=0   & r2_w_mu!=0   & r2_w_tau!=0 ~ paste(as.character(pi_class), as.character(mu_class), as.character(tau_class), sep=', ')))

ggplot(smry, aes(x=sigu_multiplier, y=est, color=model)) + 
  geom_line() + 
  geom_point() + 
  facet_grid(metric~class, scales='free_y')

smry %>%
  filter(tau_class=='none' & metric=='RMSE') %>%
  ggplot(aes(x=sigu_multiplier, y=est, color=model)) + 
  geom_line() + 
  geom_point() + 
  facet_grid(pi_class~mu_class)

smry %>%
  filter(pi_class=='none' & metric=='RMSE') %>%
  ggplot(aes(x=sigu_multiplier, y=est, color=model)) + 
  geom_line() + 
  geom_point() + 
  facet_grid(tau_class~mu_class)
