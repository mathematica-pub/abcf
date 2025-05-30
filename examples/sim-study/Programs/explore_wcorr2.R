source('programs/utils.R')

overall <- readRDS('data/mar-2023-wcorr2/all-overall.RDS')

smry <- overall %>%
  expand_grid(level=c('all','small')) %>%
  mutate(SATTbias    = ifelse(level=='all',SATTbias,    SATTsmall_hat-SATTsmall),
         SATTcover90 = ifelse(level=='all',SATTcover90, SATTsmallcover90),
         SATTwidth90 = ifelse(level=='all',SATTwidth90, SATTsmallwidth90),
         PEHT        = ifelse(level=='all', PEHT, CATEsmallRMSE)) %>%
  group_by(model, level, sigu_multiplier, r2_w_pi, r2_w_mu, r2_w_tau, w_as_covar) %>%
  summarize(est_RMSE = sqrt(mean(SATTbias^2)),
            est_cover = mean(SATTcover90),
            est_width = mean(SATTwidth90),
            est_PEHT = sqrt(mean(PEHT^2)),
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

oldsmry <- filter(smry, level=='all' & w_as_covar==FALSE)

ggplot(oldsmry, aes(x=sigu_multiplier, y=est, color=model)) + 
  geom_line() + 
  geom_point() + 
  facet_grid(metric~class, scales='free_y')

smry %>%
  filter(tau_class=='none' & metric=='RMSE' & level=='all') %>%
  ggplot(aes(x=sigu_multiplier, y=est, color=model, linetype=w_as_covar)) + 
  geom_line() + 
  geom_point() + 
  facet_grid(pi_class~mu_class)

smry %>%
  filter(pi_class=='none' & metric=='RMSE' & level=='all') %>%
  ggplot(aes(x=sigu_multiplier, y=est, color=model, linetype=w_as_covar)) + 
  geom_line() + 
  geom_point() + 
  facet_grid(tau_class~mu_class)

smry %>%
  filter(tau_class=='none' & metric=='PEHT' & level=='all') %>%
  ggplot(aes(x=sigu_multiplier, y=est, color=model, linetype=w_as_covar)) + 
  geom_line() + 
  geom_point() + 
  facet_grid(pi_class~mu_class)

smry %>%
  filter(tau_class=='none' & metric=='PEHT' & w_as_covar) %>%
  ggplot(aes(x=sigu_multiplier, y=est, color=model, linetype=level)) + 
  geom_line() + 
  geom_point() + 
  facet_grid(pi_class~mu_class)

#So as we get more correlation with tau, there's more and more problems fitting small practices, with either model
#But also, we get more of a u vs w difference for small practices at sigu=2 than for large pracs
smry %>%
  filter(pi_class=='none' & metric=='RMSE' & w_as_covar) %>%
  ggplot(aes(x=sigu_multiplier, y=est, color=model, linetype=level)) + 
  geom_line() + 
  geom_point() + 
  facet_grid(tau_class~mu_class)

smry %>%
  filter(mu_class=='none' & metric=='RMSE' & w_as_covar) %>%
  ggplot(aes(x=sigu_multiplier, y=est, color=model, linetype=level)) + 
  geom_line() + 
  geom_point() + 
  facet_grid(tau_class~pi_class)

smry %>%
  filter(pi_class=='none' & metric=='PEHT' & w_as_covar) %>%
  ggplot(aes(x=sigu_multiplier, y=est, color=model, linetype=level)) + 
  geom_line() + 
  geom_point() + 
  facet_grid(tau_class~mu_class)

smry %>%
  filter(pi_class=='none' & metric=='PEHT') %>%
  ggplot(aes(x=sigu_multiplier, y=est, linetype=model, color=paste(level, w_as_covar))) + 
  geom_line() + 
  geom_point() + 
  facet_grid(tau_class~mu_class)

smry %>%
  filter(tau_class=='none' & metric=='PEHT') %>%
  ggplot(aes(x=sigu_multiplier, y=est, linetype=model, color=paste(level, w_as_covar))) + 
  geom_line() + 
  geom_point() + 
  facet_grid(pi_class~mu_class)

smry %>%
  filter(tau_class=='none' & mu_class=='none' & metric=='PEHT') %>%
  ggplot(aes(x=sigu_multiplier, y=est, color=model)) + 
  geom_line() + 
  geom_point() + 
  facet_grid(pi_class~paste(level, w_as_covar))
