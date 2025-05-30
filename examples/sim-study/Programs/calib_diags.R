source('programs/utils.R')
source('programs/load_data.R')

calib_long <- calib %>%
  select(-var_tau_trt, -var_taux_trt, -var_yhat, -var_yhat_trt) %>%
  pivot_longer(matches('hat|check|input'), names_sep='_', names_to=c('term','.value'))

ggplot(calib_long) + 
  geom_point(aes(x=input, y=check, color=method), alpha=.1) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~ world + term, labeller=labeller(.multi_line = FALSE), scales='free')

ggplot(calib_long) + 
  geom_point(aes(x=input, y=hat, color=method), alpha=.1) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~ world + term, labeller=labeller(.multi_line = FALSE), scales='free')

ggplot(calib_long %>% filter(method=='ibcf')) + 
  geom_point(aes(x=input, y=hat, color=set), alpha=.1) + 
  geom_abline(slope=1, intercept=0) + 
  facet_wrap(~ world + term, labeller=labeller(.multi_line = FALSE), scales='free')

for (x in c('r2a','r2trt','r2uv')) {
  print(x)
  lm(hat ~ world*input, data=filter(calib_long, method=='ibcf' & term==x & set=='main')) %>%
    summary %>%
    `$`('coefficients') %>%
    as_tibble(rownames='par') %>%
    print
}

for (x in c('r2a','r2trt','r2uv')) {
  print(x)
  lm(check ~ world*input, data=filter(calib_long, method=='ibcf' & term==x & set=='main')) %>%
    summary %>%
    `$`('coefficients') %>%
    as_tibble(rownames='par') %>%
    print
}
