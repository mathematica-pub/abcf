source('programs/utils.R')
source('programs/load_data.R')

#Overall performance across main sims
overall %>%
  filter(set=='main') %>%
  group_by(world, method) %>%
  summarize(rmse = sqrt(mean(SATTbias^2)),
            mae = mean(abs(SATTbias)),
            PEHE = mean(PEHE),
            across(matches('timing|(SATT|CATE)cover'), mean)) %>%
  pivot_longer(cols=c(-world, -method), names_to='metric') %>%
  pivot_wider(id_cols=c(world, metric), names_from=method, values_from=value)

#How do we do covering the sigmas with different Ns?
mixing %>% 
  #drop pars without truths
  filter(ibcf & !is.na(cover80)) %>%
  mutate(bias = mean-real) %>%
  group_by(nT, par) %>%
  summarize(across(matches('mean|bias|cover'),~mean(.x)),
            rmse = sqrt(mean(bias^2))) %>%
  arrange(par, nT)

#Are we mucking up sigmas? yes.
mixing %>% 
  filter(ibcf & str_detect(par,'sigma2?')) %>%
  ggplot() + 
  geom_density(aes(x=mean, color=par)) + 
  geom_vline(aes(xintercept=real)) + 
  facet_wrap(~par, scales='free')

mixing %>% 
  filter(ibcf) %>%
  ggplot() + 
  geom_density(aes(x=mean, color=as.factor(nT))) + 
  geom_vline(aes(xintercept=real)) + 
  facet_wrap(~par, scales='free')

mixing %>% 
  filter(ibcf) %>%
  pivot_wider(c(sim, scenario, nT),names_from=par, values_from=mean) %>%
  ggplot() + 
  geom_point(aes(x=sigma2_u, y=sigma2_v, color=as.factor(nT)))

#A ha! So our prior in a rho-less vacuum is centered correctly, but the problem is our prior integrated over rho is not
mixing %>% 
  filter(ibcf) %>%
  pivot_wider(c(sim, scenario, nT),names_from=par, values_from=mean) %>%
  ggplot() + 
  geom_point(aes(x=sigma2_v/8.33^2, y=rho, color=as.factor(nT)))

mixing %>%
  filter(ibcf & sigu_multiplier==1 & sigv_multiplier==1 & rho==0 & par %in% c('sigma2_v','rho','sigma_v','sigv_delta')) %>%
  pivot_wider(c(nT, scenario, sim), names_from=par, values_from=c(mean, real)) %>%
  ggplot() + 
  geom_point(aes(x=mean_rho, y=mean_sigma2_v)) + 
  geom_hline(aes(yintercept=real_sigma2_v))

#Compare ibcf across rhos
overall %>%
  filter(r2uv_multiplier==1 & r2trt_multiplier==1 & set %in% c('main','norho','posrho')) %>%
  group_by(world, method, set) %>%
  summarize(rmse = sqrt(mean(SATTbias^2)),
            mae = mean(abs(SATTbias)),
            PEHE = mean(PEHE),
            across(matches('timing|(SATT|CATE)cover'), mean)) %>%
  pivot_longer(cols=c(-world, -set, -method), names_to='metric') %>%
  pivot_wider(id_cols=c(world, metric, method), names_from=set, values_from=value) %>%
  pivot_wider(id_cols=c(world, metric), names_from=method, values_from=c(main, posrho, norho))

#Compare across r2uv
overall %>%
  filter(r2trt_multiplier==1 & set == 'main') %>%
  group_by(world, method, r2uv_multiplier) %>%
  summarize(rmse = sqrt(mean(SATTbias^2)),
            mae = mean(abs(SATTbias)),
            PEHE = mean(PEHE),
            across(matches('timing|(SATT|CATE)cover'), mean)) %>%
  pivot_longer(cols=c(-world, -r2uv_multiplier, -method), names_to='metric') %>%
  pivot_wider(id_cols=c(world, metric, method), names_from=r2uv_multiplier, values_from=value) %>%
  pivot_wider(id_cols=c(world, metric), names_from=method, values_from=c(`0.75`, `1`, `1.25`))

#and across r2trt
overall %>%
  filter(r2uv_multiplier==1 & set == 'main') %>%
  group_by(world, method, r2trt_multiplier) %>%
  summarize(rmse = sqrt(mean(SATTbias^2)),
            mae = mean(abs(SATTbias)),
            PEHE = mean(PEHE),
            across(matches('timing|(SATT|CATE)cover'), mean)) %>%
  pivot_longer(cols=c(-world, -r2trt_multiplier, -method), names_to='metric') %>%
  pivot_wider(id_cols=c(world, metric, method), names_from=r2trt_multiplier, values_from=value) %>%
  pivot_wider(id_cols=c(world, metric), names_from=method, values_from=c(`0.75`, `1`, `1.25`))

#and across r2trt
overall %>%
  filter(set == 'main') %>%
  group_by(world, method, r2uv_multiplier, r2trt_multiplier) %>%
  summarize(across(matches('(SATT|CATE)cover'), mean)) %>%
  pivot_longer(cols=c(-world, -r2uv_multiplier, -r2trt_multiplier, -method), names_to='metric') %>%
  ungroup %>%
  mutate(level = as.numeric(str_extract(metric,'[0-9]+')),
         metric = str_extract(metric,'SATT|CATE'),
         thing = paste0(method, level)) %>%
  ggplot() + 
  geom_line(aes(y=value, x=r2trt_multiplier, color=method, group=thing)) + 
  geom_point(aes(y=value, x=r2trt_multiplier, color=method)) + 
  geom_hline(yintercept=0.8, linetype='dashed') +
  geom_hline(yintercept=0.9, linetype='dashed') +
  geom_hline(yintercept=0.95, linetype='dashed') +
  facet_wrap(~ metric + world + r2uv_multiplier, labeller=labeller(.multi_line = FALSE), ncol=3)


pivot_wider(id_cols=c(world, metric, method), names_from=r2trt_multiplier, values_from=value) %>%
  pivot_wider(id_cols=c(world, metric), names_from=method, values_from=c(`0.75`, `1`, `1.25`))

#coverage by sim pars
cvg <- overall %>% 
  filter(set %in% c('main','norho')) %>%
  mutate(rho = ifelse(set=='main','-','0')) %>%
  select(world, method, rho, r2uv = r2uv_multiplier, r2trt = r2trt_multiplier, matches('(SATT|CATE)cover')) %>%
  group_by(world, method, rho, r2uv, r2trt) %>%
  summarize_all(mean) %>%
  pivot_longer(matches('cover'), values_to='cover') %>%
  mutate(type = str_extract(name,'SATT|CATE'),
         nominal = as.numeric(str_extract(name,'[0-9]+')) / 100,
         grp = paste0(nominal, method))

#When rho is -ve, r2trt matters a tiny bit (higher = higher coverage) and r2uv matters a bunch (higher=lower cvg)
#When rho is negative, neither matter.
#That pattern seems to hold across ed and medicare
ggplot(filter(cvg, type=='CATE')) + 
  geom_line(aes(x=r2trt, y=cover, group=grp, color=method)) + 
  geom_line(aes(x=r2trt, y=nominal, group=grp), color='black',linetype='dashed') +
  facet_wrap(~world + rho + r2uv, labeller=labeller(.multi_line = FALSE), ncol=3)

ggplot(filter(cvg, type=='SATT')) + 
  geom_line(aes(x=r2trt, y=cover, group=grp, color=method)) + 
  geom_line(aes(x=r2trt, y=nominal, group=grp), color='black',linetype='dashed') +
  facet_wrap(~world + rho + r2uv, labeller=labeller(.multi_line = FALSE), ncol=3)

overall %>% 
  filter(set %in% c('main','norho') & world=='medicare') %>%
  mutate(rho = ifelse(set=='main','-','0')) %>%
  select(world, method, rho, r2uv = r2uv_multiplier, r2trt = r2trt_multiplier, matches('(SATT|CATE)cover')) %>%
  mutate(ibcf = method=='ibcf') %>%
  lm(SATTcover90 ~ rho*ibcf*(r2uv + r2trt), data=.) %>%
  summary %>% `$`(coefficients) %>%
  as_tibble(rownames='par') %>%
    mutate(Estimate = ifelse(abs(Estimate) < 1e-10,0,Estimate))

overall %>% 
  filter(set %in% c('main','norho') & world=='medicare') %>%
  mutate(rho = ifelse(set=='main','-','0')) %>%
  select(world, method, rho, r2uv = r2uv_multiplier, r2trt = r2trt_multiplier, matches('(SATT|CATE)cover')) %>%
  mutate(ibcf = method=='ibcf') %>%
  lm(CATEcover90 ~ rho*ibcf*(r2uv + r2trt), data=.) %>%
  summary %>% `$`(coefficients) %>%
  as_tibble(rownames='par') %>%
  mutate(Estimate = ifelse(abs(Estimate) < 1e-10,0,Estimate))

#Yep, this is the story - r2s don't matter except when rho is negative
#And there what mostly matters is uv, which makes sense since negative rho and large sigma v is what give the problems
overall %>% 
  filter(set %in% c('main','norho') & world=='medicare' & method=='ibcf') %>%
  mutate(nrho = set=='main') %>%
  select(world, method, nrho, r2uv = r2uv_multiplier, r2trt = r2trt_multiplier, matches('(SATT|CATE)cover')) %>%
  lm(CATEcover90 ~ nrho*(r2uv + r2trt), data=.) %>%
  summary %>% `$`(coefficients) %>%
  as_tibble(rownames='par') %>%
  mutate(Estimate = ifelse(abs(Estimate) < 1e-10,0,Estimate))

#Average CATE width
widths <- overall %>% 
  filter(set %in% c('main','norho')) %>%
  mutate(rho = ifelse(set=='main','-','0')) %>%
  select(world, method, rho, r2uv = r2uv_multiplier, r2trt = r2trt_multiplier, matches('(SATT|CATE)width')) %>%
  group_by(world, method, rho, r2uv, r2trt) %>%
  summarize_all(mean) %>%
  pivot_longer(matches('width'), values_to='width') %>%
  mutate(type = str_extract(name,'SATT|CATE'),
         nominal = as.numeric(str_extract(name,'[0-9]+')) / 100,
         grp = paste0(nominal, method),
         #Varies too much by ed vs medicare, but free scales is misleading, so scale ed to be ~ medicare
         width = ifelse(world=='edu',width*1000,width))

#So for medicare when rho is -ve we don't change CATE widths, but we do when rho is 0
ggplot(filter(widths, type=='CATE')) + 
  geom_line(aes(x=r2trt, y=width, group=grp, color=method)) + 
  facet_wrap(~world + rho + r2uv, labeller=labeller(.multi_line = FALSE), ncol=3)

#Weird - r2trt effects medicare but not edu
ggplot(filter(widths, type=='SATT')) + 
  geom_line(aes(x=r2trt, y=width, group=grp, color=method)) + 
  facet_wrap(~world + rho + r2uv, labeller=labeller(.multi_line = FALSE), ncol=3)
#yuuuup
overall %>% 
  filter(set %in% c('main','norho') & world=='medicare' & method=='ibcf') %>%
  mutate(nrho = set=='main') %>%
  select(world, method, nrho, r2uv = r2uv_multiplier, r2trt = r2trt_multiplier, matches('(SATT|CATE)width')) %>%
  lm(CATEwidth95 ~ nrho*(r2uv + r2trt), data=.) %>%
  summary %>% `$`(coefficients) %>%
  as_tibble(rownames='par') %>%
  mutate(Estimate = ifelse(abs(Estimate) < 1e-10,0,Estimate))
