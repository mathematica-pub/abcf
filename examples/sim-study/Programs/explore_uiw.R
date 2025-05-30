source('programs/utils.R')
overall <- readRDS('data/mar-2023-uiw/all-overall.RDS')
mixing  <- readRDS('data/mar-2023-uiw/all-mixing.RDS')
ex_smy  <- readRDS('data/mar-2023-uiw/all-exemplar-summy.RDS')

pw = 9
plotdir = 'data/mar-2023-uiw'
plotsave <- function(name, w=pw, h=w*9/16, d=plotdir) {
  ggsave(filename=glue('{d}/{name}.png'), height=h, width=w)
}

pal$wbcf <- pal$teal
pal$ibcf <- pal$purple
pal$ubcf <- pal$gold

fig1s <- list(sigv_satt=danfig1(overall, varying='sig_v', estimand='SATT', ci=TRUE),
              sigv_catt=danfig1(overall, varying='sig_v', estimand='CATT', ci=TRUE),
              sigv_resid=danfig1(overall, varying='sig_v', estimand='resid', ci=TRUE),
              sigu_satt=danfig1(overall, varying='sig_u', estimand='SATT', ci=TRUE),
              sigu_catt=danfig1(overall, varying='sig_u', estimand='CATT', ci=TRUE),
              sigu_resid=danfig1(overall, varying='sig_u', estimand='resid', ci=TRUE))

cowplot::plot_grid(plotlist=fig1s[c(1,4)], nrow = 2)              
plotsave('fig1_satt')
cowplot::plot_grid(plotlist=fig1s[c(2,5)], nrow = 2)              
plotsave('fig1_catt')
cowplot::plot_grid(plotlist=fig1s[c(3,6)], nrow = 2)              
plotsave('fig1_resid')

confus <- ex_smy$resid$confus %>%
  filter(sigv_multiplier==1 & cutoff=='gt_45' & size %in% c('all','small')) %>%
  group_by(method, sigu_multiplier, size) %>%
  summarize(pos_est = mean(pos),
            pos_lb=quantile(pos,.05),
            pos_ub=quantile(pos,.95),
            gt45_est = mean(gt_45),
            gt45_lb=quantile(gt_45,.05),
            gt45_ub=quantile(gt_45,.95),) %>%
  arrange(sigu_multiplier, size, method) %>%
  mutate(sigu_multiplier = factor(round(61.3*sigu_multiplier))) %>%
  pivot_longer(matches('pos|gt45'), names_sep='_', names_to=c('metric','.value')) %>%
  mutate(real = ifelse(metric=='pos',.5,1-pnorm(45, sd=sqrt(8.33^2 + as.numeric(as.character(sigu_multiplier))^2))),
         all=1)

ggplot(confus) + 
  geom_point(aes(x=sigu_multiplier, y=est, color=method, group=method)) + 
  geom_line(aes(x=sigu_multiplier, y=est, color=method, group=method)) + 
  geom_line(aes(x=sigu_multiplier, y=real, group=all), linetype='dashed', color='black') + 
  geom_ribbon(aes(x=sigu_multiplier, ymin=lb, ymax=ub, fill=method, group=method), alpha=.1) + 
  facet_grid(metric~size)
plotsave('fig4_exemplar_gt45')

confus %>%
  filter(metric=='gt45' & size=='all') %>%
  ggplot() + 
  geom_point(aes(x=sigu_multiplier, y=est, color=method, group=method)) + 
  geom_line(aes(x=sigu_multiplier, y=real, group=all,color='coin'), linetype='dashed') + 
  scale_color_manual(values=c('coin'='transparent','wbcf'='transparent', 'ibcf'='transparent', 'ubcf'='transparent')) +
  theme(legend.text=element_text(color='transparent'))+
  labs(color=NULL, fill=NULL,
       y='Percent of exemplars actually over 45',
       x=expression(sigma[u]))
plotsave('fig4_exemplar_gt45_00', w=6)

confus %>%
  filter(metric=='gt45' & size=='all') %>%
  ggplot() + 
  geom_point(aes(x=sigu_multiplier, y=est, color=method, group=method)) + 
  geom_line(aes(x=sigu_multiplier, y=real, group=all,color='coin'), linetype='dashed') + 
  scale_color_manual(values=c('coin'='black','wbcf'='transparent', 'ibcf'='transparent', 'ubcf'='transparent')) +
  labs(color=NULL, fill=NULL,
       y='Percent of exemplars actually over 45',
       x=expression(sigma[u]))
plotsave('fig4_exemplar_gt45_0', w=6)

confus %>%
  filter(metric=='gt45' & size=='all') %>%
  ggplot() + 
  geom_point(aes(x=sigu_multiplier, y=est, color=method, group=method)) + 
  geom_line(aes(x=sigu_multiplier, y=est, color=method, group=method)) + 
  geom_line(aes(x=sigu_multiplier, y=real, group=all), linetype='dashed', color='black') + 
  geom_ribbon(aes(x=sigu_multiplier, ymin=lb, ymax=ub, fill=method, group=method, alpha=method)) +
  scale_color_manual(values=c('coin'='black','wbcf'=pal$wbcf, 'ibcf'='transparent', 'ubcf'='transparent')) +
  scale_fill_manual(values=c('coin'='black','wbcf'=pal$wbcf, 'ibcf'='transparent', 'ubcf'='transparent')) +
  scale_alpha_manual(values=c('coin'=.1,'wbcf'=.1, 'ibcf'=0, 'ubcf'=0)) +
  labs(color=NULL, fill=NULL, alpha=NULL,
       y='Percent of exemplars actually over 45',
       x=expression(sigma[u]))
plotsave('fig4_exemplar_gt45_1', w=6)

confus %>%
  filter(metric=='gt45' & size=='all') %>%
  ggplot() + 
  geom_point(aes(x=sigu_multiplier, y=est, color=method, group=method)) + 
  geom_line(aes(x=sigu_multiplier, y=est, color=method, group=method)) + 
  geom_line(aes(x=sigu_multiplier, y=real, group=all), linetype='dashed', color='black') + 
  geom_ribbon(aes(x=sigu_multiplier, ymin=lb, ymax=ub, fill=method, group=method, alpha=method)) +
  scale_color_manual(values=c('coin'='black','wbcf'=pal$wbcf, 'ibcf'=pal$ibcf, 'ubcf'='transparent')) +
  scale_fill_manual(values=c('coin'='black','wbcf'=pal$wbcf, 'ibcf'=pal$ibcf, 'ubcf'='transparent')) +
  scale_alpha_manual(values=c('coin'=.1,'wbcf'=.1, 'ibcf'=.1, 'ubcf'=0)) +
  labs(color=NULL, fill=NULL,alpha=NULL,
       y='Percent of exemplars actually over 45',
       x=expression(sigma[u]))
plotsave('fig4_exemplar_gt45_2', w=6)

confus %>%
  filter(metric=='gt45' & size=='all') %>%
  ggplot() + 
  geom_point(aes(x=sigu_multiplier, y=est, color=method, group=method)) + 
  geom_line(aes(x=sigu_multiplier, y=est, color=method, group=method)) + 
  geom_line(aes(x=sigu_multiplier, y=real, group=all), linetype='dashed', color='black') + 
  geom_ribbon(aes(x=sigu_multiplier, ymin=lb, ymax=ub, fill=method, group=method), alpha=.1) +
  scale_color_manual(values=c('coin'='black','wbcf'=pal$wbcf, 'ibcf'=pal$ibcf, 'ubcf'=pal$ubcf)) +
  scale_fill_manual(values=c('coin'='black','wbcf'=pal$wbcf, 'ibcf'=pal$ibcf, 'ubcf'=pal$ubcf)) +
  labs(color=NULL, fill=NULL,
       y='Percent of exemplars actually over 45',
       x=expression(sigma[u]))
plotsave('fig4_exemplar_gt45_3', w=6)

confus %>%
  filter(metric=='gt45' & size=='small') %>%
  ggplot() + 
  geom_point(aes(x=sigu_multiplier, y=est, color=method, group=method)) + 
  geom_line(aes(x=sigu_multiplier, y=est, color=method, group=method)) + 
  geom_line(aes(x=sigu_multiplier, y=real, group=all), linetype='dashed', color='black') + 
  geom_ribbon(aes(x=sigu_multiplier, ymin=lb, ymax=ub, fill=method, group=method), alpha=.1) +
  scale_color_manual(values=c('coin'='black','wbcf'=pal$wbcf, 'ibcf'=pal$ibcf, 'ubcf'=pal$ubcf)) +
  scale_fill_manual(values=c('coin'='black','wbcf'=pal$wbcf, 'ibcf'=pal$ibcf, 'ubcf'=pal$ubcf)) +
  labs(color=NULL, fill=NULL,
       y='Percent of small exemplars actually over 45',
       x=expression(sigma[u]))
plotsave('fig4_exemplar_gt45_small', w=6)

sigv_est <- mixing %>%
  filter(par=='sigma_v' & sigu_multiplier==1) %>%
  expand_grid(temp=1:2) %>%
  filter(method=='ubcf' | temp==1) %>%
  mutate(method = ifelse(temp==2 & method=='ubcf', 'wbcf',method)) %>%
  group_by(method, sigv_multiplier) %>%
  summarize(real= mean(real),
            est = mean(mean),
            lb  = quantile(mean,.05),
            ub  = quantile(mean,.95)) %>%
  ungroup %>%
  mutate(sigv_multiplier = factor(round(8.33*sigv_multiplier)))

ggplot(sigv_est, aes(x=real)) + 
  geom_line(aes(y=real, color='truth', group='truth'), linetype='dashed') +
  geom_line(aes(y=est, group=method, color=method)) +
  geom_ribbon(aes(ymin=lb, ymax=ub, group=method, fill=method), alpha=.1) +
  geom_point(aes(y=est, group=method, color=method)) +
  scale_color_manual(values=c('truth'='black','wbcf'='transparent', 'ibcf'='transparent', 'ubcf'='transparent')) +
  scale_fill_manual(values=c('truth'='black','wbcf'='transparent', 'ibcf'='transparent', 'ubcf'='transparent')) +
  labs(color=NULL, fill=NULL,
       y=bquote('Estimate of'~sigma[v]),
       x=bquote('True'~sigma[v])) +
  scale_x_continuous(breaks=sort(unique(sigv_est$real)), labels=as.numeric(levels(sigv_est$sigv_multiplier))) + 
  scale_y_continuous(breaks=sort(unique(sigv_est$real)), labels=as.numeric(levels(sigv_est$sigv_multiplier))) + 
  theme(panel.grid.minor = element_blank())
plotsave('fig5_sigv_0', w=6)

ggplot(sigv_est, aes(x=real)) + 
  geom_line(aes(y=real, color='truth', group='truth'), linetype='dashed') +
  geom_line(aes(y=est, group=method, color=method)) +
  geom_ribbon(aes(ymin=lb, ymax=ub, group=method, fill=method), alpha=.1) +
  geom_point(aes(y=est, group=method, color=method)) +
  scale_color_manual(values=c('truth'='black','wbcf'=pal$wbcf, 'ibcf'='transparent', 'ubcf'='transparent')) +
  scale_fill_manual(values=c('truth'='black','wbcf'=pal$wbcf, 'ibcf'='transparent', 'ubcf'='transparent')) +
  labs(color=NULL, fill=NULL,
       y=bquote('Estimate of'~sigma[v]),
       x=bquote('True'~sigma[v])) +
  scale_x_continuous(breaks=sort(unique(sigv_est$real)), labels=as.numeric(levels(sigv_est$sigv_multiplier))) + 
  scale_y_continuous(breaks=sort(unique(sigv_est$real)), labels=as.numeric(levels(sigv_est$sigv_multiplier))) + 
  theme(panel.grid.minor = element_blank())
plotsave('fig5_sigv_1', w=6)

ggplot(sigv_est, aes(x=real)) + 
  geom_line(aes(y=real, color='truth', group='truth'), linetype='dashed') +
  geom_line(aes(y=est, group=method, color=method)) +
  geom_ribbon(aes(ymin=lb, ymax=ub, group=method, fill=method), alpha=.1) +
  geom_point(aes(y=est, group=method, color=method)) +
  scale_color_manual(values=c('truth'='black','wbcf'=pal$wbcf, 'ibcf'='transparent', 'ubcf'=pal$ubcf)) +
  scale_fill_manual(values=c('truth'='black','wbcf'=pal$wbcf, 'ibcf'='transparent', 'ubcf'=pal$ubcf)) +
  labs(color=NULL, fill=NULL,
       y=bquote('Estimate of'~sigma[v]),
       x=bquote('True'~sigma[v])) +
  scale_x_continuous(breaks=sort(unique(sigv_est$real)), labels=as.numeric(levels(sigv_est$sigv_multiplier))) + 
  scale_y_continuous(breaks=sort(unique(sigv_est$real)), labels=as.numeric(levels(sigv_est$sigv_multiplier))) + 
  theme(panel.grid.minor = element_blank())
plotsave('fig5_sigv_2', w=6)

ggplot(sigv_est %>% 
         filter(method!='ibcf' | sigv_multiplier==8) %>% 
         mutate(lb = ifelse(method=='ibcf',lb,NA_real_),
                ub = ifelse(method=='ibcf',ub,NA_real_)), 
       aes(x=real)) + 
  geom_line(aes(y=real, color='truth', group='truth'), linetype='dashed') +
  geom_line(aes(y=est, group=method, color=method)) +
  geom_errorbar(aes(ymin=lb, ymax=ub, group=method, color=method)) +
  geom_point(aes(y=est, group=method, color=method)) +
  scale_color_manual(values=c('truth'='black','wbcf'=pal$wbcf, 'ibcf'=pal$ibcf, 'ubcf'=pal$ubcf)) +
  scale_fill_manual(values=c('truth'='black','wbcf'=pal$wbcf, 'ibcf'=pal$ibcf, 'ubcf'=pal$ubcf)) +
  labs(color=NULL, fill=NULL,
       y=bquote('Estimate of'~sigma[v]),
       x=bquote('True'~sigma[v])) +
  scale_x_continuous(breaks=sort(unique(sigv_est$real)), labels=as.numeric(levels(sigv_est$sigv_multiplier)), limits=c(0, max(sigv_est$real))) + 
  scale_y_continuous(breaks=sort(unique(sigv_est$real)), labels=as.numeric(levels(sigv_est$sigv_multiplier))) + 
  theme(panel.grid.minor = element_blank())
plotsave('fig5_sigv_3', w=6)

ggplot(sigv_est, aes(x=real)) + 
  geom_line(aes(y=real, color='truth', group='truth'), linetype='dashed') +
  geom_line(aes(y=est, group=method, color=method)) +
  geom_ribbon(aes(ymin=lb, ymax=ub, group=method, fill=method), alpha=.25) +
  geom_point(aes(y=est, group=method, color=method)) +
  scale_color_manual(values=c('truth'='black','wbcf'=pal$wbcf, 'ibcf'=pal$ibcf, 'ubcf'=pal$ubcf)) +
  scale_fill_manual(values=c('truth'='black','wbcf'=pal$wbcf, 'ibcf'=pal$ibcf, 'ubcf'=pal$ubcf)) +
  labs(color=NULL, fill=NULL,
       y=bquote('Estimate of'~sigma[v]),
       x=bquote('True'~sigma[v])) +
  scale_x_continuous(breaks=sort(unique(sigv_est$real)), labels=as.numeric(levels(sigv_est$sigv_multiplier))) + 
  scale_y_continuous(breaks=sort(unique(sigv_est$real)), labels=as.numeric(levels(sigv_est$sigv_multiplier))) + 
  theme(panel.grid.minor = element_blank())
plotsave('fig5_sigv_4', w=6)

sigu_est <- mixing %>%
  filter(par=='sigma_u' & sigv_multiplier==1) %>%
  expand_grid(temp=1:2) %>%
  filter(method=='ubcf' | temp==1) %>%
  mutate(method = ifelse(temp==2 & method=='ubcf', 'wbcf',method),
         mean = ifelse(method=='wbcf',0,mean)) %>%
  group_by(method, sigu_multiplier) %>%
  summarize(real= mean(real),
            est = mean(mean),
            lb  = quantile(mean,.05),
            ub  = quantile(mean,.95)) %>%
  ungroup %>%
  mutate(sigu_multiplier = factor(round(61.3*sigu_multiplier)))

ggplot(sigu_est, aes(x=real)) + 
  geom_line(aes(y=real, color='truth', group='truth'), linetype='dashed') +
  geom_line(aes(y=est, group=method, color=method)) +
  geom_ribbon(aes(ymin=lb, ymax=ub, group=method, fill=method), alpha=0) +
  geom_point(aes(y=est, group=method, color=method)) +
  scale_color_manual(values=c('truth'='black','wbcf'='transparent', 'ibcf'='transparent', 'ubcf'='transparent')) +
  scale_fill_manual(values=c('truth'='black','wbcf'='transparent', 'ibcf'='transparent', 'ubcf'='transparent')) +
  labs(color=NULL, fill=NULL,
       y=bquote('Estimate of'~sigma[u]),
       x=bquote('True'~sigma[u])) +
  scale_x_continuous(breaks=sort(unique(sigu_est$real)), labels=as.numeric(levels(sigu_est$sigu_multiplier))) + 
  scale_y_continuous(breaks=sort(unique(sigu_est$real)), labels=as.numeric(levels(sigu_est$sigu_multiplier))) +
  theme(panel.grid.minor = element_blank())
plotsave('fig6_sigu_0', w=6)

ggplot(sigu_est, aes(x=real)) + 
  geom_line(aes(y=real, color='truth', group='truth'), linetype='dashed') +
  geom_line(aes(y=est, group=method, color=method)) +
  geom_ribbon(aes(ymin=lb, ymax=ub, group=method, fill=method), alpha=0) +
  geom_point(aes(y=est, group=method, color=method)) +
  scale_color_manual(values=c('truth'='black','wbcf'=pal$wbcf, 'ibcf'='transparent', 'ubcf'='transparent')) +
  scale_fill_manual(values=c('truth'='black','wbcf'=pal$wbcf, 'ibcf'='transparent', 'ubcf'='transparent')) +
  labs(color=NULL, fill=NULL,
       y=bquote('Estimate of'~sigma[u]),
       x=bquote('True'~sigma[u])) +
  scale_x_continuous(breaks=sort(unique(sigu_est$real)), labels=as.numeric(levels(sigu_est$sigu_multiplier))) + 
  scale_y_continuous(breaks=sort(unique(sigu_est$real)), labels=as.numeric(levels(sigu_est$sigu_multiplier))) +
  theme(panel.grid.minor = element_blank())
plotsave('fig6_sigu_1', w=6)

ggplot(sigu_est %>% 
         filter(method!='ibcf' | sigu_multiplier==61) %>%
         mutate(lb = ifelse(method=='ibcf',lb,NA_real_),
                ub = ifelse(method=='ibcf',ub,NA_real_)), 
       aes(x=real)) + 
  geom_line(aes(y=real, color='truth', group='truth'), linetype='dashed') +
  geom_line(aes(y=est, group=method, color=method)) +
  geom_errorbar(aes(ymin=lb, ymax=ub, group=method, color=method)) +
  geom_point(aes(y=est, group=method, color=method)) +
  scale_color_manual(values=c('truth'='black','wbcf'=pal$wbcf, 'ibcf'=pal$ibcf, 'ubcf'='transparent')) +
  scale_fill_manual(values=c('truth'='black','wbcf'=pal$wbcf, 'ibcf'=pal$ibcf, 'ubcf'='transparent')) +
  labs(color=NULL, fill=NULL,
       y=bquote('Estimate of'~sigma[u]),
       x=bquote('True'~sigma[u])) +
  scale_x_continuous(breaks=sort(unique(sigu_est$real)), labels=as.numeric(levels(sigu_est$sigu_multiplier)), limits=c(0, max(sigu_est$real))) + 
  scale_y_continuous(breaks=sort(unique(sigu_est$real)), labels=as.numeric(levels(sigu_est$sigu_multiplier)), limits=c(0, max(sigu_est$ub))) +
  theme(panel.grid.minor = element_blank())
plotsave('fig6_sigu_2', w=6)

ggplot(sigu_est, aes(x=real)) + 
  geom_line(aes(y=real, color='truth', group='truth'), linetype='dashed') +
  geom_line(aes(y=est, group=method, color=method)) +
  geom_ribbon(aes(ymin=lb, ymax=ub, group=method, fill=method), alpha=0.1) +
  geom_point(aes(y=est, group=method, color=method)) +
  scale_color_manual(values=c('truth'='black','wbcf'=pal$wbcf, 'ibcf'=pal$ibcf, 'ubcf'='transparent')) +
  scale_fill_manual(values=c('truth'='black','wbcf'=pal$wbcf, 'ibcf'=pal$ibcf, 'ubcf'='transparent')) +
  labs(color=NULL, fill=NULL,
       y=bquote('Estimate of'~sigma[u]),
       x=bquote('True'~sigma[u])) +
  scale_x_continuous(breaks=sort(unique(sigu_est$real)), labels=as.numeric(levels(sigu_est$sigu_multiplier))) + 
  scale_y_continuous(breaks=sort(unique(sigu_est$real)), labels=as.numeric(levels(sigu_est$sigu_multiplier))) +
  theme(panel.grid.minor = element_blank())
plotsave('fig6_sigu_3', w=6)

ggplot(sigu_est, aes(x=real)) + 
  geom_line(aes(y=real, color='truth', group='truth'), linetype='dashed') +
  geom_line(aes(y=est, group=method, color=method)) +
  geom_ribbon(aes(ymin=lb, ymax=ub, group=method, fill=method), alpha=.1) +
  geom_point(aes(y=est, group=method, color=method)) +
  scale_color_manual(values=c('truth'='black','wbcf'=pal$wbcf, 'ibcf'=pal$ibcf, 'ubcf'=pal$ubcf)) +
  scale_fill_manual(values=c('truth'='black','wbcf'=pal$wbcf, 'ibcf'=pal$ibcf, 'ubcf'=pal$ubcf)) +
  labs(color=NULL, fill=NULL,
       y=bquote('Estimate of'~sigma[u]),
       x=bquote('True'~sigma[u])) +
  scale_x_continuous(breaks=sort(unique(sigu_est$real)), labels=as.numeric(levels(sigu_est$sigu_multiplier))) + 
  scale_y_continuous(breaks=sort(unique(sigu_est$real)), labels=as.numeric(levels(sigu_est$sigu_multiplier))) +
  theme(panel.grid.minor = element_blank())
plotsave('fig6_sigu_4', w=6)

fig2s <- list(sigu_est   = danfig2(mixing, varying='sig_u', ci=TRUE, metric='est'),
              sigu_cover = danfig2(mixing, varying='sig_u', ci=TRUE, metric='cover'),
              sigv_est   = danfig2(mixing, varying='sig_v', ci=TRUE, metric='est'),
              sigv_cover = danfig2(mixing, varying='sig_v', ci=TRUE, metric='cover'))
cowplot::plot_grid(plotlist=fig2s[1:2])
plotsave('fig2_sigu')
cowplot::plot_grid(plotlist=fig2s[3:4])
plotsave('fig2_sigv')



#Is the overall RMSE actually different? Within seed comparison
matched <- overall %>%
  filter(method!='ibcf' & sigv_multiplier==1) %>%
  pivot_wider(id_cols=c(sigu_multiplier, sim), names_from=method, values_from=c(SATT, SATT_hat, SATTbias)) %T>%
  {mutate(., diff = abs(SATT_ubcf-SATT_wbcf)) %>% pull(diff) %>% summary %>% print} %>%
  filter(SATT_ubcf==SATT_wbcf) %>%
  mutate(SATT=SATT_ubcf) %>%
  select(-SATT_ubcf, -SATT_wbcf) 

matched %>%
  split({.}$sigu_multiplier) %>%
  lapply(function(x){
    summary(lm(SATT ~ SATT_hat_ubcf + SATT_hat_wbcf, data=x))$coefficients
  })

matched %>% 
  mutate(d = abs(SATTbias_wbcf) - abs(SATTbias_ubcf)) %>%
  group_by(sigu_multiplier) %>%
  summarize(meand = mean(d),
            q25d = quantile(d,.25),
            q50d = quantile(d,.50),
            q75d = quantile(d,.75))

matched %>% 
  mutate(d = abs(SATTbias_wbcf) - abs(SATTbias_ubcf)) %>%
  ggplot() + 
  geom_density(aes(x=d)) + 
  geom_vline(xintercept=0) +
  facet_wrap(~sigu_multiplier, scales='free')

indiv_resid  <- readRDS('data/mar-2023-uiw/all-indiv_resid.RDS')
ir_base <- filter(indiv_resid, sigu_multiplier==1 & sigv_multiplier)
ir_wide <- indiv_resid %>%
  select(sigu_multiplier, sigv_multiplier, sim, id, z, truth, mean, method) %>%
  pivot_wider(names_from=method, values_from=mean)
ir_base_wide <- ir_base %>%
  select(sim, id, z, truth, mean, method) %>%
  pivot_wider(names_from=method, values_from=mean)

indiv_uv  <- readRDS('data/mar-2023-uiw/all-indiv_uv.RDS')
iv_base <- filter(indiv_uv, sigu_multiplier==1 & sigv_multiplier & par=='v' & z==1)
iv_wide <- indiv_uv %>%
  filter(par=='v' & z==1) %>%
  select(sigu_multiplier, sigv_multiplier, sim, id, z, truth, mean, method) %>%
  pivot_wider(names_from=method, values_from=mean)
iv_base_wide <- iv_base %>%
  select(sim, id, z, truth, mean, method) %>%
  pivot_wider(names_from=method, values_from=mean)

cor(ir_base$truth[ir_base$method=='wbcf'], ir_base$mean[ir_base$method=='wbcf'])
cor(ir_base$truth[ir_base$method=='ubcf'], ir_base$mean[ir_base$method=='ubcf'])
cor(ir_base$truth[ir_base$method=='ibcf'], ir_base$mean[ir_base$method=='ibcf'])
cor(ir_base$truth[ir_base$method=='ubcf' & ir_base$z==1], ir_base$mean[ir_base$method=='ubcf' & ir_base$z==1])
cor(ir_base$truth[ir_base$method=='ibcf' & ir_base$z==1], ir_base$mean[ir_base$method=='ibcf' & ir_base$z==1])
cor(ir_base$truth[ir_base$method=='ubcf' & ir_base$z==0], ir_base$mean[ir_base$method=='ubcf' & ir_base$z==0])
cor(ir_base$truth[ir_base$method=='ibcf' & ir_base$z==0], ir_base$mean[ir_base$method=='ibcf' & ir_base$z==0])


summary(ir_base_wide)
lm(truth ~ ibcf, data=ir_base_wide)
lm(truth ~ ubcf, data=ir_base_wide)
lm(truth ~ wbcf, data=ir_base_wide)

lm(ibcf ~ truth, data=ir_base_wide)
lm(ubcf ~ truth, data=ir_base_wide)
lm(wbcf ~ truth, data=ir_base_wide)

lm(truth ~ ubcf + ibcf, data=ir_base_wide)
lm(truth-ubcf ~ ibcf, data=ir_base_wide)
lm(truth-ibcf ~ ubcf, data=ir_base_wide)
lm(ibcf ~ ubcf, data=ir_base_wide)
lm(ubcf ~ ibcf, data=ir_base_wide)

summary(iv_base_wide)
lm(truth ~ ibcf, data=iv_base_wide)
lm(ibcf ~ truth, data=iv_base_wide)

ir_wide

ir_wide %>%
  filter(sigv_multiplier==1) %>%
  mutate(sigu = as.factor(round(61.3*sigu_multiplier))) %>%
  lm(truth ~ sigu + sigu:ubcf, data=.) %>%
  summary

ir_wide %>%
  filter(sigv_multiplier==1) %>%
  mutate(sigu = as.factor(round(61.3*sigu_multiplier))) %>%
  lm(ubcf ~ sigu + sigu:truth, data=.) %>%
  summary

iv_wide %>%
  filter(sigu_multiplier==1) %>%
  mutate(sigv = as.factor(round(8.33*sigv_multiplier))) %>%
  lm(truth ~ sigv + sigv:ibcf, data=.) %>%
  summary

iv_wide %>%
  filter(sigu_multiplier==1) %>%
  mutate(sigv = as.factor(round(8.33*sigv_multiplier))) %>%
  lm(ibcf ~ sigv + sigv:truth, data=.) %>%
  summary
