source('programs/utils.R')
overall <- readRDS('Data/feb-2023-nonlin/all-overall.RDS')
mixing <- readRDS('Data/feb-2023-nonlin/all-mixing.RDS')
ex_smy <- readRDS('Data/feb-2023-nonlin/all-exemplar-summy.RDS')
msc_list <- readRDS('Data/feb-2023-nonlin/all-misc.RDS')
calib <- readRDS('Data/feb-2023-nonlin/all-calib.RDS')
scenarios <- readRDS('Data/feb-2023-nonlin/_scenarios.RDS')

overall %>%
  group_by(model, nonlin_frac) %>%
  summarize(SATTRMSE = sqrt(mean(SATTbias^2)),
            SATTwidth90 = mean(SATTwidth90),
            SATTcover90 = mean(SATTcover90))

overall %>%
  group_by(model, nonlin_frac) %>%
  summarize(CATTRMSE = sqrt(mean(PEHT^2)),
            CATTwidth90 = mean(CATTwidth90),
            CATTcover90 = mean(CATTcover90))

ex_smy$resid$rmse %>%
  filter(group=='Exemplars') %>%
  group_by(model, nonlin_frac, cutoff, size) %>%
  summarize(RMSE = sqrt(mean(rmse^2))) %>%
  filter(size=='all' & cutoff=='gt_45')

ex_smy$resid$confus %>%
  group_by(model, nonlin_frac, cutoff, size) %>%
  summarize(across(c(pos, matches('gt')), mean)) %>%
  filter(size=='small' & cutoff=='pos')

ex_smy$resid$auc %>%
  group_by(model, nonlin_frac, cutoff, size) %>%
  summarize(auc = mean(auc)) %>%
  filter(size=='all' & cutoff=='gt_45')

#Gah, I forgot to multiply sigma_y for bcf through by root mean weight sqrt(608)
mixing %>%
  mutate(par = ifelse(method=='wbcf' & par=='sigma','sigma_y',par)) %>%
  filter(str_detect(par,'sigma_|tau_bar')) %>%
  group_by(method, par) %>%
  summarize(across(c(real, mean, sd, cover90), ~mean(.x))) %>%
  arrange(par, method)

mixing %>%
  mutate(par = ifelse(method=='wbcf' & par=='sigma','sigma_y',par)) %>%
  filter(str_detect(par,'sigma_|tau_bar')) %>%
  group_by(method, par, nonlin_frac) %>%
  summarize(across(c(real, mean, sd, cover90), ~mean(.x))) %>%
  arrange(par, method, nonlin_frac) %>%
  print(n=Inf)

mixing %>%
  mutate(par = ifelse(method=='wbcf' & par=='sigma','sigma_y',par)) %>%
  filter(str_detect(par,'sigma_[uy]')) %>%
  group_by(par, method, nonlin_frac) %>%
  summarize(across(c(real, mean, sd, cover90), ~mean(.x))) %>%
  rename(truth=real, post_mean = mean, post_sd=sd) %>%
  mutate(post_mean = ifelse(method=='stan' & par=='sigma_y', post_mean*sqrt(608), post_mean),
         post_sd   = ifelse(method=='stan' & par=='sigma_y', post_sd  *sqrt(608), post_sd),
         cover90   = ifelse(method=='stan' & par=='sigma_y', NA_real_, cover90)) %>%
  arrange(par, method, nonlin_frac) %>%
  print(n=Inf)
