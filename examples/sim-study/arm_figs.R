source('programs/utils.R')
overall <- readRDS('data/aug-2023-uold/all-overall.RDS')
mixing <- readRDS('data/aug-2023-uold/all-mixing.RDS')
ex_smy  <- readRDS('data/aug-2023-uold/all-exemplar-summy.RDS')

icc_overall <- readRDS('data/oct-2023-paper-icc/all-overall.RDS')
uhp_overall <- readRDS('data/oct-2023-paper-uhp/all-overall.RDS')
uhp_mix     <- readRDS('data/oct-2023-paper-uhp/all-mixing.RDS')
ibcf_overall <- readRDS('data/oct-2023-paper-ibcf/all-overall.RDS')

feb_overall <- readRDS('data/feb-2023-runs/all-overall.RDS')
feb_mixing  <- readRDS('data/feb-2023-runs/all-mixing.RDS')
mar_overall <- readRDS('data/mar-2023-uiw/all-overall.RDS')
mar_mixing  <- readRDS('data/mar-2023-uiw/all-mixing.RDS')

figh = 4.5
figw = 8

pal = list(bcf = scales::hue_pal()(2)[1],
           abcf = scales::hue_pal()(2)[2])

pal = list(bcf = '#046B5C',
           abcf = '#5C4377')

pal = list(bcf = '#65A7A9',
           abcf = '#414B66')

pal = list(bcf = '#0B2949',
           abcf = '#046B5C')


#######################################
#######################################
#     ARM
#######################################
#######################################

psatt <- overall %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1) %>%
  group_by(model) %>%
  summarize(SATT_RMSE = sqrt(mean((SATT_hat-SATT)^2)),
            SATT_cover90 = mean(SATTcover90),
            SATT_width90 = mean(SATTwidth90),
            UTE_RMSE = sqrt(mean(PEHT^2)),
            UTE_cover90 = mean(CATTcover90),
            UTE_width90 = mean(CATTwidth90)) %>%
  pivot_longer(-model, names_sep='_', names_to=c('estimand','metric')) %>%
  mutate(nominal = ifelse(metric=='cover90',.9,NA_real_),
         metric = factor(metric, levels=c('RMSE','width90','cover90'), labels=c('RMSE','90% interval width','90% interval coverage')),
         model = factor(model, levels=c('oldBCF','uBCF'), labels=c('BCF','aBCF'))) %>%
  filter(estimand=='SATT') %>%
  ggplot() + 
  geom_col(aes(x=model, y=value, fill=model)) + 
  facet_wrap( ~ metric, scales='free_y') +
  geom_hline(aes(yintercept=nominal), linetype='dashed') +
  ggh4x::facetted_pos_scales(y=list(scale_y_continuous(expand=expansion(mult=c(0, .05)), labels=scales::dollar),
                                    scale_y_continuous(expand=expansion(mult=c(0, .05)), labels=scales::dollar),
                                    scale_y_continuous(expand=expansion(mult=c(0, .05)), labels=scales::percent, breaks = c(0, .25, .5, .75, .9)))) +
  scale_fill_manual(values = c(BCF=pal$bcf, aBCF=pal$abcf)) +
  labs(x='',
       y='') + 
  theme_bw() +
  theme(legend.position='none',
        strip.background =element_rect(fill="white"))

psatt + scale_fill_manual(values = c(BCF=pal$bcf, aBCF='transparent'))
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/ARM 2025 aBCF/figs/satt_abcf_vs_bcf_intro.png', height=figh, width=figw)

psatt
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/ARM 2025 aBCF/figs/satt_abcf_vs_bcf.png', height=figh, width=figw)

pute <- overall %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1) %>%
  group_by(model) %>%
  summarize(SATT_RMSE = sqrt(mean((SATT_hat-SATT)^2)),
            SATT_cover90 = mean(SATTcover90),
            SATT_width90 = mean(SATTwidth90),
            UTE_RMSE = sqrt(mean(PEHT^2)),
            UTE_cover90 = mean(CATTcover90),
            UTE_width90 = mean(CATTwidth90)) %>%
  pivot_longer(-model, names_sep='_', names_to=c('estimand','metric')) %>%
  mutate(nominal = ifelse(metric=='cover90',.9,NA_real_),
         metric = factor(metric, levels=c('RMSE','width90','cover90'), labels=c('RMSE','90% interval width','90% interval coverage')),
         model = factor(model, levels=c('oldBCF','uBCF'), labels=c('BCF','aBCF'))) %>%
  filter(estimand=='UTE') %>%
  ggplot() + 
  geom_col(aes(x=model, y=value, fill=model)) + 
  facet_wrap( ~ metric, scales='free_y') +
  geom_hline(aes(yintercept=nominal), linetype='dashed') +
  ggh4x::facetted_pos_scales(y=list(scale_y_continuous(expand=expansion(mult=c(0, .05)), labels=scales::dollar),
                                    scale_y_continuous(expand=expansion(mult=c(0, .05)), labels=scales::dollar),
                                    scale_y_continuous(expand=expansion(mult=c(0, .05)), labels=scales::percent, breaks = c(0, .25, .5, .75, .9)))) +
  scale_fill_manual(values = c(BCF=pal$bcf, aBCF=pal$abcf)) +
  labs(x='',
       y='') + 
  theme_bw() +
  theme(legend.position='none',
        strip.background =element_rect(fill="white"))

pute
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/ARM 2025 aBCF/figs/ute_abcf_vs_bcf.png', height=figh, width=figw)

sigu_lvls <- c('BCF',.17,.33,.67,1.33,2.67)
psigu_satt <- overall %>%
  mutate(sigu_hyperprior=ifelse(model=='oldBCF','BCF','0.67')) %>%
  bind_rows(uhp_overall %>% mutate(sigu_hyperprior = as.character(round(sigu_hyperprior,2)))) %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1) %>%
  group_by(sigu_hyperprior, model) %>%
  summarize(SATT_RMSE = sqrt(mean((SATT_hat-SATT)^2)),
            SATT_cover90 = mean(SATTcover90),
            SATT_width90 = mean(SATTwidth90),
            UTE_RMSE = sqrt(mean(PEHT^2)),
            UTE_cover90 = mean(CATTcover90),
            UTE_width90 = mean(CATTwidth90),
            resid_RMSE=sqrt(mean(residT_RMSET^2)),
            resid_cover90=mean(residTcover90),
            resid_width90=mean(residTwidth90)) %>%
  pivot_longer(-c(sigu_hyperprior, model), names_sep='_', names_to=c('estimand','metric')) %>%
  mutate(nominal = ifelse(metric=='cover90',.9,NA_real_),
         metric = factor(metric, levels=c('RMSE','width90','cover90'), labels=c('RMSE','90% interval width','90% interval coverage')),
         model = factor(model, levels=c('oldBCF','uBCF'), labels=c('BCF','aBCF')),
         sigu_hyperprior = factor(sigu_hyperprior, levels=sigu_lvls)) %>%
  filter(estimand=='SATT') %>%
  ggplot() + 
  geom_col(aes(x=sigu_hyperprior, y=value, fill=model, alpha=sigu_hyperprior)) + 
  facet_wrap(~ metric, scales='free_y') +
  geom_hline(aes(yintercept=nominal), linetype='dashed') +
  scale_alpha_manual(values=c('BCF'=1, '0.67'=1, '0.17'=.5, '0.33'=.5, '1.33'=.5, '2.67'=.5)) +
  ggh4x::facetted_pos_scales(y=list(scale_y_continuous(expand=expansion(mult=c(0, .05)), labels=scales::dollar),
                                    scale_y_continuous(expand=expansion(mult=c(0, .05)), labels=scales::dollar),
                                    scale_y_continuous(expand=expansion(mult=c(0, .05)), labels=scales::percent, breaks = c(0, .25, .5, .75, .9)))) +
  scale_fill_manual(values = c(BCF=pal$bcf, aBCF=pal$abcf)) +
  guides(alpha='none') +
  labs(x=bquote('Hyperprior for '*sigma[u]*', in terms of '*sd(y)),
       y='',
       fill='') +
  theme_bw() +
  theme(legend.position='bottom',
        strip.background =element_rect(fill="white"))

psigu_satt + scale_alpha_manual(values=c('BCF'=1, '0.67'=1, '0.17'=0, '0.33'=0, '1.33'=0, '2.67'=0))
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/ARM 2025 aBCF/figs/satt_siguhp_sens_intro.png', height=figh, width=figw)

psigu_satt
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/ARM 2025 aBCF/figs/satt_siguhp_sens.png', height=figh, width=figw)

psigu_ute <- overall %>%
  mutate(sigu_hyperprior=ifelse(model=='oldBCF','BCF','0.67')) %>%
  bind_rows(uhp_overall %>% mutate(sigu_hyperprior = as.character(round(sigu_hyperprior,2)))) %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1) %>%
  group_by(sigu_hyperprior, model) %>%
  summarize(SATT_RMSE = sqrt(mean((SATT_hat-SATT)^2)),
            SATT_cover90 = mean(SATTcover90),
            SATT_width90 = mean(SATTwidth90),
            UTE_RMSE = sqrt(mean(PEHT^2)),
            UTE_cover90 = mean(CATTcover90),
            UTE_width90 = mean(CATTwidth90),
            resid_RMSE=sqrt(mean(residT_RMSET^2)),
            resid_cover90=mean(residTcover90),
            resid_width90=mean(residTwidth90)) %>%
  pivot_longer(-c(sigu_hyperprior, model), names_sep='_', names_to=c('estimand','metric')) %>%
  mutate(nominal = ifelse(metric=='cover90',.9,NA_real_),
         metric = factor(metric, levels=c('RMSE','width90','cover90'), labels=c('RMSE','90% interval width','90% interval coverage')),
         model = factor(model, levels=c('oldBCF','uBCF'), labels=c('BCF','aBCF')),
         sigu_hyperprior = factor(sigu_hyperprior, levels=sigu_lvls)) %>%
  filter(estimand=='UTE') %>%
  ggplot() + 
  geom_col(aes(x=sigu_hyperprior, y=value, fill=model, alpha=sigu_hyperprior)) + 
  facet_wrap(~ metric, scales='free_y') +
  geom_hline(aes(yintercept=nominal), linetype='dashed') +
  scale_alpha_manual(values=c('BCF'=1, '0.67'=1, '0.17'=.5, '0.33'=.5, '1.33'=.5, '2.67'=.5)) +
  ggh4x::facetted_pos_scales(y=list(scale_y_continuous(expand=expansion(mult=c(0, .05)), labels=scales::dollar),
                                    scale_y_continuous(expand=expansion(mult=c(0, .05)), labels=scales::dollar),
                                    scale_y_continuous(expand=expansion(mult=c(0, .05)), labels=scales::percent, breaks = c(0, .25, .5, .75, .9)))) +
  scale_fill_manual(values = c(BCF=pal$bcf, aBCF=pal$abcf)) +
  guides(alpha='none') +
  labs(x=bquote('Hyperprior for '*sigma[u]*', in terms of '*sd(y)),
       y='',
       fill='') +
  theme_bw() +
  theme(legend.position='bottom',
        strip.background =element_rect(fill="white"))

psigu_ute + scale_alpha_manual(values=c('BCF'=1, '0.67'=1, '0.17'=0, '0.33'=0, '1.33'=0, '2.67'=0))
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/ARM 2025 aBCF/figs/ute_siguhp_sens_intro.png', height=figh, width=figw)

psigu_ute
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/ARM 2025 aBCF/figs/ute_siguhp_sens.png', height=figh, width=figw)

psigu_sigu <- mixing %>%
  mutate(sigu_hyperprior=ifelse(model=='oldBCF','BCF','0.67')) %>%
  bind_rows(uhp_mix %>% mutate(sigu_hyperprior = as.character(round(sigu_hyperprior,2)))) %>%
  filter(sigu_multiplier==1 & sigv_multiplier==1 & par=='sigma_u') %>%
  group_by(sigu_hyperprior, model) %>% 
  summarize(su_RMSE = sqrt(mean((mean-real)^2)),
            su_cover90 = mean(cover90),
            su_width90 = mean(width90)) %>%
  pivot_longer(-c(sigu_hyperprior, model), names_sep='_', names_to=c('estimand','metric')) %>%
  mutate(value = ifelse(model=='oldBCF', NA, value),
         nominal = ifelse(metric=='cover90',.9,NA_real_),
         metric = factor(metric, levels=c('RMSE','width90','cover90'), labels=c('RMSE','90% interval width','90% interval coverage')),
         model = factor(model, levels=c('oldBCF','uBCF'), labels=c('BCF','aBCF')),
         sigu_hyperprior = factor(sigu_hyperprior, levels=sigu_lvls)) %>%
  ggplot() +
  geom_col(aes(x=sigu_hyperprior, y=value, fill=model, alpha=sigu_hyperprior)) + 
  facet_wrap(~ metric, scales='free_y') +
  geom_hline(aes(yintercept=nominal), linetype='dashed') +
  scale_alpha_manual(values=c('BCF'=1, '0.67'=1, '0.17'=.5, '0.33'=.5, '1.33'=.5, '2.67'=.5)) +
  ggh4x::facetted_pos_scales(y=list(scale_y_continuous(expand=expansion(mult=c(0, .05)), labels=scales::dollar),
                                    scale_y_continuous(expand=expansion(mult=c(0, .05)), labels=scales::dollar),
                                    scale_y_continuous(expand=expansion(mult=c(0, .05)), labels=scales::percent, breaks = c(0, .25, .5, .75, .9)))) +
  scale_fill_manual(values = c(BCF=pal$bcf, aBCF=pal$abcf)) +
  guides(alpha='none') +
  labs(x=bquote('Hyperprior for '*sigma[u]*', in terms of '*sd(y)),
       y='',
       fill='') +
  theme_bw() +
  theme(legend.position='none',
        strip.background =element_rect(fill="white"))

psigu_sigu + scale_alpha_manual(values=c('BCF'=1, '0.67'=1, '0.17'=0, '0.33'=0, '1.33'=0, '2.67'=0))
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/ARM 2025 aBCF/figs/sigu_siguhp_sens_intro.png', height=figh, width=figw)

psigu_sigu + scale_alpha_manual(values=c('BCF'=1, '0.67'=1, '0.17'=0, '0.33'=0, '1.33'=.5, '2.67'=.5))
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/ARM 2025 aBCF/figs/sigu_siguhp_sens_mid.png', height=figh, width=figw)

psigu_sigu
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/ARM 2025 aBCF/figs/sigu_siguhp_sens.png', height=figh, width=figw)

#ICC figs
r2es <- seq(0,1,.1)
iccs <- c(0,.0002,.0004,.0007,.0011,.0016,.0024,.0038,.0065,.0146,1)
options(scipen=999)
lbls <- paste0(round(r2es,2), '\n(', iccs, ')')
base_r2e <- 61.3^2 / (61.3^2 + 2557^2/608)

picc_satt <- icc_overall %>%
  mutate(sens='sens') %>%
  bind_rows(overall %>% filter(sigu_multiplier==1 & sigv_multiplier==1) %>% mutate(r2_error = base_r2e, sens='main')) %>%
  group_by(model, r2_error, sens) %>%
  summarize(SATT_RMSE = sqrt(mean((SATT_hat-SATT)^2)),
            SATT_cover90 = mean(SATTcover90),
            SATT_width90 = mean(SATTwidth90),
            UTE_RMSE = sqrt(mean(PEHT^2)),
            UTE_cover90 = mean(CATTcover90),
            UTE_width90 = mean(CATTwidth90)) %>%
  pivot_longer(-c(model, r2_error, sens), names_sep='_', names_to=c('estimand','metric')) %>%
  mutate(nominal = ifelse(metric=='cover90',.90,NA_real_),
         metric = factor(metric, levels=c('RMSE','width90','cover90'), labels=c('RMSE','90% interval width','90% interval coverage')),
         model = factor(model, levels=c('oldBCF','uBCF'), labels=c('BCF','aBCF')),
         dumb = 'dumb') %>%
  filter(estimand=='SATT') %>%
  ggplot() + 
  geom_line(aes(x=r2_error, y=value, color=model, linetype=dumb)) + 
  geom_point(aes(x=r2_error, y=value, color=model, alpha=sens)) + 
  facet_wrap(~ metric, scales='free_y', ncol=3, labeller=labeller(.multi_line = FALSE)) +
  geom_hline(aes(yintercept=nominal), linetype='dashed') +
  scale_x_continuous(breaks=c(0, base_r2e, 1), labels=c(0, 'Sim\ndefault', 1), minor_breaks=r2es) +
  ggh4x::facetted_pos_scales(y=list(scale_y_continuous(limits=c(0, 10), breaks = seq(0, 10, 2), labels=scales::dollar(seq(0, 10, 2)), expand = expansion(mult=c(0, 0))),
                                    scale_y_continuous(limits=c(10,30), expand = expansion(mult=c(0, 0)), labels=scales::dollar),
                                    scale_y_continuous(limits=c(.7,1), expand = expansion(mult=c(0, 0)), labels=scales::percent))) +
  scale_color_manual(values = c(BCF=pal$bcf, aBCF=pal$abcf)) +
  scale_alpha_manual(values=c(main=1, sens=1)) +
  theme_bw() +
  theme(legend.position='bottom',
        strip.background =element_rect(fill="white"),
        panel.grid.major.x=element_blank()) +
  guides(alpha='none',
         linetype='none') +
  labs(x='ICC',
       y='',
       color='')

picc_satt + scale_alpha_manual(values=c(main=1, sens=0)) + scale_linetype_manual(values=c(dumb='blank'))
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/ARM 2025 aBCF/figs/icc_sens_satt_intro.png', width=figw, height=figh)

picc_satt
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/ARM 2025 aBCF/figs/icc_sens_satt.png', width=figw, height=figh)

picc_ute <- icc_overall %>%
  mutate(sens='sens') %>%
  bind_rows(overall %>% filter(sigu_multiplier==1 & sigv_multiplier==1) %>% mutate(r2_error = 61.3^2 / (61.3^2 + 2557^2/608), sens='main')) %>%
  group_by(model, r2_error, sens) %>%
  summarize(SATT_RMSE = sqrt(mean((SATT_hat-SATT)^2)),
            SATT_cover90 = mean(SATTcover90),
            SATT_width90 = mean(SATTwidth90),
            UTE_RMSE = sqrt(mean(PEHT^2)),
            UTE_cover90 = mean(CATTcover90),
            UTE_width90 = mean(CATTwidth90)) %>%
  pivot_longer(-c(model, r2_error, sens), names_sep='_', names_to=c('estimand','metric')) %>%
  mutate(nominal = ifelse(metric=='cover90',.9,NA_real_),
         metric = factor(metric, levels=c('RMSE','width90','cover90'), labels=c('RMSE','90% interval width','90% interval coverage')),
         model = factor(model, levels=c('oldBCF','uBCF'), labels=c('BCF','aBCF')),
         dumb = 'dumb') %>%
  filter(estimand=='UTE') %>%
  ggplot() + 
  geom_line(aes(x=r2_error, y=value, color=model, linetype=dumb)) + 
  geom_point(aes(x=r2_error, y=value, color=model, alpha=sens)) + 
  facet_wrap(~ metric, scales='free_y', ncol=3, labeller=labeller(.multi_line = FALSE)) +
  geom_hline(aes(yintercept=nominal), linetype='dashed') +
  scale_x_continuous(breaks=c(0, base_r2e, 1), labels=c(0, 'Sim\ndefault', 1), minor_breaks = r2es) +
  ggh4x::facetted_pos_scales(y=list(scale_y_continuous(limits=c(0, 20), expand = expansion(mult=c(0, 0)), labels=scales::dollar),
                                    scale_y_continuous(limits=c(30,50), expand = expansion(mult=c(0, 0)), labels=scales::dollar),
                                    scale_y_continuous(limits=c(.7,1), expand = expansion(mult=c(0, 0)), labels=scales::percent))) +
  scale_color_manual(values = c(BCF=pal$bcf, aBCF=pal$abcf)) +
  scale_alpha_manual(values=c(main=1, sens=1)) +
  theme_bw() +
  theme(legend.position='bottom',
        strip.background =element_rect(fill="white"),
        panel.grid.major.x=element_blank()) +
  guides(alpha='none',
         linetype='none') +
  labs(x='Residual variance share',
       y='',
       color='')

picc_ute + scale_alpha_manual(values=c(main=1, sens=0)) + scale_linetype_manual(values=c(dumb='blank'))
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/ARM 2025 aBCF/figs/icc_sens_UTE_intro.png', width=figw, height=figh)

picc_ute
ggsave('C:/Users/dthal/OneDrive - Mathematica/Desktop/Papers and presentations/ARM 2025 aBCF/figs/icc_sens_UTE.png', width=figw, height=figh)
