source('programs/utils.R')
overall <- readRDS('data/jan-2023-runs/all-overall.RDS')

fig1s <- expand_grid(varying=c('sig_v','sig_u','rho','nT'),
                     estimand=c('SATT','CATT')) %>%
  pmap(danfig1, overall=overall, ci=TRUE)

fig2s <- expand_grid(varying=c('sig_v','sig_u','rho','nT'),
                     metric=c('cover','rmse')) %>%
  pmap(danfig2, mixing=mixing, ci=TRUE)

fig3s <- expand_grid(sigv_const=c(0,.25, .5,2/3,1,1.5,2,4),
                     metric=c('tau','v')) %>%
  pmap(danfig3, ex_smy=ex_smy, show='diff')

names(fig1s) <- outer(c('SATT', 'CATT'), c('sig_v','sig_u','rho','nT'), paste, sep='_') %>% as.vector
names(fig2s) <- outer(c('cover','rmse'), c('sig_v','sig_u','rho','nT'), paste, sep='_') %>% as.vector
names(fig3s) <- outer(c('tau_diff','v_diff'), round(c(0,.25, .5,2/3,1,1.5,2,4),2), paste, sep='_') %>% as.vector

cowplot::plot_grid(plotlist=lapply(fig1s, function(x) x + theme(legend.position='none')))
cowplot::plot_grid(plotlist=lapply(fig2s, function(x) x + theme(legend.position='none')))


readRDS('data/jan-2023-runs/rhohard0/all-mixing.RDS') %>% 
  filter(ibcf) %>%
  ggplot() + 
  geom_density(aes(x=mean)) + 
  geom_vline(aes(xintercept=real)) + 
  facet_wrap(~par, scales='free')

readRDS('data/jan-2023-runs/rhohard0/all-mixing.RDS') %>% 
  group_by(par) %>%
  summarize(real=mean(real),
            avg = mean(mean),
            sd = sd(mean),
            q25 = quantile(mean,.25),
            q75 = quantile(mean,.75))

readRDS('data/jan-2023-runs/rhohard0/all-mixing.RDS') %>% 
  group_by(par) %>%
  summarize(min_neff = min(n_eff),
            q025_neff = quantile(n_eff, .025, na.rm=TRUE),
            max_rhat = max(Rhat),
            q975_Rhat = quantile(Rhat, .975, na.rm=TRUE))


#Why so bimodal? Could also be a dumb computer use thing with the parallelization.
overall %>%
  filter(rho==0 & sigu_multiplier==1 & sigv_multiplier==1) %>%
  ggplot() + 
  geom_density(aes(x=timing, color=ibcf)) + 
  facet_wrap(~nT)

overall %>% select(scenario, sim, ibcf, nT, timing) %>%
  pivot_wider(c(scenario, sim, nT), values_from=timing, names_from=ibcf) %>%
  ggplot() + 
  geom_point(aes(x=`FALSE`, y=`TRUE`, color=as.factor(nT))) + 
  geom_abline()

#Not due to tauscale
overall %>%
  filter(rho==0 & sigu_multiplier==1 & sigv_multiplier==1) %>%
  inner_join(mixing %>%
               filter(par=='tau_scale') %>%
               select(scenario, sim, ibcf, mean, n_eff, Rhat),
  by=c('scenario', 'sim', 'ibcf')) %>%
  ggplot() +
  geom_point(aes(x=Rhat, y=timing, color=ibcf)) + 
  facet_wrap(~nT)

overall %>%
  filter(rho==0 & sigu_multiplier==1 & sigv_multiplier==1) %>%
  inner_join(mixing %>%
               filter(par=='mu_scale') %>%
               select(scenario, sim, ibcf, mean, n_eff, Rhat),
             by=c('scenario', 'sim', 'ibcf')) %>%
  ggplot() +
  geom_point(aes(x=mean, y=timing, color=ibcf)) + 
  facet_wrap(~nT)

make_tree_plot <- function (x, y, cp=.1, pfx='runtime', ...) {
  tree_data <- data.frame(x, yhat = y)
  tree_fit <- rpart::rpart(yhat ~ .,data = tree_data) # Want to determine drivers of cont. propensities
  print(tree_fit$cptable)
  
  tree_fit_pruned <- if(is.null(cp)){tree_fit} else {rpart::prune(tree_fit, cp = cp)}
  rpart.plot::prp(tree_fit_pruned, faclen = 0, clip.facs = TRUE, 
                  box.palette = "YlGnBl", 
                  compress = TRUE, ycompress = TRUE, extra = 101, prefix = paste0(pfx,': '), 
                  shadow.col = "gray", type = 4, varlen = 0)
}

xxx <- overall %>%
  filter(rho==0 & sigu_multiplier==1 & sigv_multiplier==1 & nT==2000) %>%
  inner_join(mixing %>%
               pivot_wider(c(scenario, sim, ibcf), names_from=par, values_from=c(mean, Rhat)),
             by=c('scenario', 'sim', 'ibcf'))

xxxi <- filter(xxx, ibcf)

make_tree_plot(xxxi %>% select(matches('mean|Rhat')), xxxi$timing)

overall %>%
  filter(rho==0 & sigu_multiplier==1 & sigv_multiplier==1) %>%
  group_by(nT, ibcf) %>%
  summarize(mean = mean(timing),
            q5 = quantile(timing,.05),
            q95 = quantile(timing,.95)) %>%
  ggplot() + 
  geom_point(aes(x=nT, y=mean, color=ibcf)) + 
  geom_line(aes(x=nT, y=mean, color=ibcf, group=ibcf)) + 
  geom_ribbon(aes(x=nT, ymin=q5, ymax=q95, fill=ibcf), alpha=.1)
