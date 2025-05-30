source('programs/utils.R')

from <- 's3://bcf-sim-study/dec-2021-runs/4000burn-10ksim'
to <- 'C:/Projects/bcf-sim-study/Data/dec-2021-runs/4000burn-10ksim'
#system(glue('bash -c "aws s3 sync {from} {to} --profile bcf-sims --region us-east-1 --exclude "*control-files/*" "'))

#Too big to load all at once, try separating into chunks
files <- tibble(file = dir(to, pattern = '*.RDS', recursive = FALSE, full.names = TRUE)) %>%
  mutate(name = file %>% str_remove('C:/Projects/bcf-sim-study/Data/dec-2021-runs/4000burn-10ksim/') %>% str_remove('\\.RDS'),
         set = case_when(str_detect(file,'weightsF') ~ 'no_weights',
                         str_detect(file,'uvdistt') ~ 't_dist',
                         str_detect(file,'sigu0') ~ 'no_res',
                         str_detect(file,'het') ~ 'main_het',
                         str_detect(file,'homog') ~ 'main_homog'),
         sim = as.numeric(str_match(file,'sim([0-9][0-9][0-9])')[,2]),
         method = str_extract(name,'ibcf|wbcf'))

combine_output_chunks(files, dir='C:/Projects/bcf-sim-study/Data/dec-2021-runs')

combine_combined_output(prefixes=c('all','all'),
                        dirs=c('C:/Projects/bcf-sim-study/Data/dec-2021-runs',
                               'C:/Projects/bcf-sim-study/Data/dec-2021-hcrho-runs'),
                        small=TRUE, big=TRUE,
                        outfix='all',
                        outdir='C:/Projects/bcf-sim-study/Data')

sets <- split(files,files$set)

save_small <- TRUE
save_big <- TRUE

for (i in 1:length(sets)) {
  set <- names(sets)[i]
  print(glue('{set}'))
  fs <- sets[[i]]$file
  fits <- overalls <- mixings <- indivs <- indiv_uvs <- exemplars <- ex_summys <- miscs <- datas <- vector('list',length(fs))
  for (j in 1:length(fs)) {
    if (j==1 | j%%100 == 0) print(glue('  {j}'))
    fit <- readRDS(fs[j])
    #Add sim and set to input object
    fit$inputs$set <- set
    fit$inputs$sim <- sets[[i]]$sim[j]
    fit$inputs <- create_group(fit$inputs, drop=FALSE)
    
    fits[[j]] <- fit
    #Now that that's saved, drop the long set of identifiers
    fit$inputs <- select(fit$inputs, set, scenario, ibcf, sim)
    
    overalls[[j]] <- bind_cols(fit$inputs, fit$overall)
    mixings[[j]] <- bind_cols(fit$inputs, fit$mixing)
    indivs[[j]] <- bind_cols(fit$inputs, fit$indiv)
    exemplars[[j]] <- lapply(fit$exemplar, function(x) bind_cols(fit$inputs, x))
    ex_summys[[j]] <- fit[c('inputs', 'exemplar_summy')]
    datas[[j]] <- bind_cols(fit$inputs, fit$data)
    if (fit$inputs$ibcf) {
      indiv_uvs[[j]] <- bind_cols(fit$inputs, fit$indiv_uv)
      miscs[[j]] <- fit[c('inputs','postcorr','scalecorr','acceptance')]
    } else {
      miscs[[j]] <- fit[c('inputs','scalecorr','acceptance')]
    }
  }
  
  if (save_small) {
    overalls %>% bind_rows %>% saveRDS(glue('Data/{set}-overall.RDS'))
    mixings %>% bind_rows %>% saveRDS(glue('Data/{set}-mixing.RDS'))
    #Bind together ex_summys so constituent pieces are bind_rows type data sets
    ex_summys <- list(tau = list(auc    = ex_summys %>% lapply(function(x) bind_cols(x$inputs, x$exemplar_summy$tau$auc)) %>% bind_rows,
                                 rmse   = ex_summys %>% lapply(function(x) bind_cols(x$inputs, x$exemplar_summy$tau$rmse)) %>% bind_rows,
                                 confus = ex_summys %>% lapply(function(x) bind_cols(x$inputs, x$exemplar_summy$tau$confus)) %>% bind_rows,
                                 prob_pos = list(auc    = ex_summys %>% lapply(function(x) bind_cols(x$inputs, x$exemplar_summy$tau$prob_pos$auc)) %>% bind_rows,
                                                 rmse   = ex_summys %>% lapply(function(x) bind_cols(x$inputs, x$exemplar_summy$tau$prob_pos$rmse)) %>% bind_rows,
                                                 confus = ex_summys %>% lapply(function(x) bind_cols(x$inputs, x$exemplar_summy$tau$prob_pos$confus)) %>% bind_rows),
                                 prob_gt_ate = list(auc    = ex_summys %>% lapply(function(x) bind_cols(x$inputs, x$exemplar_summy$tau$prob_gt_ate$auc)) %>% bind_rows,
                                                    rmse   = ex_summys %>% lapply(function(x) bind_cols(x$inputs, x$exemplar_summy$tau$prob_gt_ate$rmse)) %>% bind_rows,
                                                    confus = ex_summys %>% lapply(function(x) bind_cols(x$inputs, x$exemplar_summy$tau$prob_gt_ate$confus)) %>% bind_rows)),
                      v   = list(auc    = ex_summys %>% lapply(function(x) bind_cols(x$inputs, x$exemplar_summy$v$auc)) %>% bind_rows %>% filter(ibcf),
                                 rmse   = ex_summys %>% lapply(function(x) bind_cols(x$inputs, x$exemplar_summy$v$rmse)) %>% bind_rows %>% filter(ibcf),
                                 confus = ex_summys %>% lapply(function(x) bind_cols(x$inputs, x$exemplar_summy$v$confus)) %>% bind_rows %>% filter(ibcf),
                                 prob_pos = list(auc    = ex_summys %>% lapply(function(x) bind_cols(x$inputs, x$exemplar_summy$v$prob_pos$auc)) %>% bind_rows %>% filter(ibcf),
                                                 rmse   = ex_summys %>% lapply(function(x) bind_cols(x$inputs, x$exemplar_summy$v$prob_pos$rmse)) %>% bind_rows %>% filter(ibcf),
                                                 confus = ex_summys %>% lapply(function(x) bind_cols(x$inputs, x$exemplar_summy$v$prob_pos$confus)) %>% bind_rows %>% filter(ibcf)),
                                 prob_gt_ate = list(auc    = ex_summys %>% lapply(function(x) bind_cols(x$inputs, x$exemplar_summy$v$prob_gt_ate$auc)) %>% bind_rows %>% filter(ibcf),
                                                    rmse   = ex_summys %>% lapply(function(x) bind_cols(x$inputs, x$exemplar_summy$v$prob_gt_ate$rmse)) %>% bind_rows %>% filter(ibcf),
                                                    confus = ex_summys %>% lapply(function(x) bind_cols(x$inputs, x$exemplar_summy$v$prob_gt_ate$confus)) %>% bind_rows %>% filter(ibcf))))
    saveRDS(ex_summys, glue('Data/dec-2021-runs/{set}-exemplar-summy.RDS'))
    saveRDS(miscs,     glue('Data/dec-2021-runs/{set}-misc.RDS'))
  }
  
  if (save_big) {
    indivs %>% bind_rows %>% saveRDS(glue('Data/dec-2021-runs/{set}-indiv.RDS'))
    indiv_uvs %>% bind_rows %>% saveRDS(glue('Data/dec-2021-runs/{set}-indiv_uv.RDS'))
    #Exemplars is a set of lists of tau and v, bind together separately so we still have different tau/v lists
    list(tau = exemplars %>% lapply(`[[`,'tau') %>% bind_rows,
         v   = exemplars %>% lapply(`[[`,'v') %>% bind_rows) %>%
      saveRDS(glue('Data/dec-2021-runs/{set}-exemplar.RDS'))
    
    #Don't bind fits or data together
    saveRDS(fits,  glue('Data/dec-2021-runs/{set}-fits.RDS'))
    saveRDS(datas, glue('Data/dec-2021-runs/{set}-data.RDS'))  
  }
  
  fits <- overalls <- mixings <- indivs <- ex_summys <- miscs <- datas <- NULL
  gc()
}

#combined versions of small files
if (save_small) {
  bind_rows(readRDS('Data/dec-2021-runs/main_het-overall.RDS'),
            readRDS('Data/dec-2021-runs/main_homog-overall.RDS'),
            readRDS('Data/dec-2021-runs/no_res-overall.RDS'),
            readRDS('Data/dec-2021-runs/no_weights-overall.RDS'),
            readRDS('Data/dec-2021-runs/t_dist-overall.RDS')) %>%
    saveRDS('Data/dec-2021-runs/all-overall.RDS')
  
  bind_rows(readRDS('Data/dec-2021-runs/main_het-mixing.RDS'),
            readRDS('Data/dec-2021-runs/main_homog-mixing.RDS'),
            readRDS('Data/dec-2021-runs/no_res-mixing.RDS'),
            readRDS('Data/dec-2021-runs/no_weights-mixing.RDS'),
            readRDS('Data/dec-2021-runs/t_dist-mixing.RDS')) %>%
    saveRDS('Data/dec-2021-runs/all-mixing.RDS')
  
  #Fuck this
  ex_summys <- list(a=readRDS('Data/dec-2021-runs/main_het-exemplar-summy.RDS'),
                    b=readRDS('Data/dec-2021-runs/main_homog-exemplar-summy.RDS'),
                    c=readRDS('Data/dec-2021-runs/no_res-exemplar-summy.RDS'),
                    d=readRDS('Data/dec-2021-runs/no_weights-exemplar-summy.RDS'),
                    f=readRDS('Data/dec-2021-runs/t_dist-exemplar-summy.RDS'))
  ex_smy <- list(tau = list(auc    = ex_summys %>% lapply(function(x) x$tau$auc) %>% bind_rows,
                            rmse   = ex_summys %>% lapply(function(x) x$tau$rmse) %>% bind_rows,
                            confus = ex_summys %>% lapply(function(x) x$tau$confus) %>% bind_rows,
                            prob_pos = list(auc    = ex_summys %>% lapply(function(x) x$tau$prob_pos$auc) %>% bind_rows,
                                            rmse   = ex_summys %>% lapply(function(x) x$tau$prob_pos$rmse) %>% bind_rows,
                                            confus = ex_summys %>% lapply(function(x) x$tau$prob_pos$confus) %>% bind_rows),
                            prob_gt_ate = list(auc    = ex_summys %>% lapply(function(x) x$tau$prob_gt_ate$auc) %>% bind_rows,
                                               rmse   = ex_summys %>% lapply(function(x) x$tau$prob_gt_ate$rmse) %>% bind_rows,
                                               confus = ex_summys %>% lapply(function(x) x$tau$prob_gt_ate$confus) %>% bind_rows)),
                 v   = list(auc    = ex_summys %>% lapply(function(x) x$v$auc) %>% bind_rows,
                            rmse   = ex_summys %>% lapply(function(x) x$v$rmse) %>% bind_rows,
                            confus = ex_summys %>% lapply(function(x) x$v$confus) %>% bind_rows,
                            prob_pos = list(auc    = ex_summys %>% lapply(function(x) x$v$prob_pos$auc) %>% bind_rows,
                                            rmse   = ex_summys %>% lapply(function(x) x$v$prob_pos$rmse) %>% bind_rows,
                                            confus = ex_summys %>% lapply(function(x) x$v$prob_pos$confus) %>% bind_rows),
                            prob_gt_ate = list(auc    = ex_summys %>% lapply(function(x) x$v$prob_gt_ate$auc) %>% bind_rows,
                                               rmse   = ex_summys %>% lapply(function(x) x$v$prob_gt_ate$rmse) %>% bind_rows,
                                               confus = ex_summys %>% lapply(function(x) x$v$prob_gt_ate$confus) %>% bind_rows)))
  
  saveRDS(ex_smy, 'Data/dec-2021-runs/all-exemplar-summy.RDS')
  
  c(readRDS('Data/dec-2021-runs/main_het-misc.RDS'),
    readRDS('Data/dec-2021-runs/main_homog-misc.RDS'),
    readRDS('Data/dec-2021-runs/no_res-misc.RDS'),
    readRDS('Data/dec-2021-runs/no_weights-misc.RDS'),
    readRDS('Data/dec-2021-runs/t_dist-misc.RDS')) %>%
    saveRDS('Data/dec-2021-runs/all-misc.RDS')
}

if (save_big) {
  indiv <- bind_rows(readRDS('Data/dec-2021-runs/main_het-indiv.RDS'),
                     readRDS('Data/dec-2021-runs/main_homog-indiv.RDS'),
                     readRDS('Data/dec-2021-runs/no_weights-indiv.RDS'),
                     readRDS('Data/dec-2021-runs/no_res-indiv.RDS'),
                     readRDS('Data/dec-2021-runs/t_dist-indiv.RDS'))
  saveRDS(indiv,'Data/dec-2021-runs/all-indiv.RDS')
  #Also save cate widths
  indiv %>%
    group_by(set, scenario, ibcf, sim) %>%
    summarize(CATEwidth80 = mean(width80),
              CATEwidth90 = mean(width90),
              CATEwidth95 = mean(width95),
              CATTwidth80 = mean(ifelse(z==1,width80,NA_real_), na.rm=TRUE),
              CATTwidth90 = mean(ifelse(z==1,width90,NA_real_), na.rm=TRUE),
              CATTwidth95 = mean(ifelse(z==1,width95,NA_real_), na.rm=TRUE),
              CATUwidth80 = mean(ifelse(z==0,width80,NA_real_), na.rm=TRUE),
              CATUwidth90 = mean(ifelse(z==0,width90,NA_real_), na.rm=TRUE),
              CATUwidth95 = mean(ifelse(z==0,width95,NA_real_), na.rm=TRUE)) %>%
    ungroup %>%
    saveRDS('Data/dec-2021-runs/all-cate_width.RDS')
  
  bind_rows(readRDS('Data/dec-2021-runs/main_het-indiv_uv.RDS'),
            readRDS('Data/dec-2021-runs/main_homog-indiv_uv.RDS'),
            readRDS('Data/dec-2021-runs/no_weights-indiv_uv.RDS'),
            readRDS('Data/dec-2021-runs/no_res-indiv_uv.RDS'),
            readRDS('Data/dec-2021-runs/t_dist-indiv_uv.RDS')) %>%
    saveRDS('Data/dec-2021-runs/all-indiv_uv.RDS')
  
  #For size, don't save columns we can get elsewhere
  bind_rows(readRDS('Data/dec-2021-runs/main_het-data.RDS') %>% bind_rows %>% select(-z_str, -sigma_u, -sigma_v, -rho, -sigma_y, -sigma),
            readRDS('Data/dec-2021-runs/main_homog-data.RDS') %>% bind_rows %>% select(-z_str, -sigma_u, -sigma_v, -rho, -sigma_y, -sigma),
            readRDS('Data/dec-2021-runs/no_weights-data.RDS') %>% bind_rows %>% select(-z_str, -sigma_u, -sigma_v, -rho, -sigma_y, -sigma),
            readRDS('Data/dec-2021-runs/no_res-data.RDS') %>% bind_rows %>% select(-z_str, -sigma_u, -sigma_v, -rho, -sigma_y, -sigma),
            readRDS('Data/dec-2021-runs/t_dist-data.RDS') %>% bind_rows %>% select(-z_str, -sigma_u, -sigma_v, -rho, -sigma_y, -sigma)) %>%
    saveRDS('Data/dec-2021-runs/all-data.RDS')
  
  exemp <- list(readRDS('Data/dec-2021-runs/main_het-exemplar.RDS'),
                readRDS('Data/dec-2021-runs/main_homog-exemplar.RDS'),
                readRDS('Data/dec-2021-runs/no_weights-exemplar.RDS'),
                readRDS('Data/dec-2021-runs/no_res-exemplar.RDS'),
                readRDS('Data/dec-2021-runs/t_dist-exemplar.RDS'))
  
  list(tau = lapply(exemp,`[[`,'tau') %>% bind_rows,
       v   = lapply(exemp,`[[`,'v')   %>% bind_rows) %>%
    saveRDS('Data/dec-2021-runs/all-exemplar.RDS')
}
#Save full scenario file
read_csv('Data/dec-2021-runs/sim_study_scenarios.csv') %>%
  mutate(set = case_when(str_detect(scenario_desc,'weightsF') ~ 'no_weights',
                         str_detect(scenario_desc,'uvdistt') ~ 't_dist',
                         str_detect(scenario_desc,'sigu0') ~ 'no_res',
                         str_detect(scenario_desc,'het') ~ 'main_het',
                         str_detect(scenario_desc,'homog') ~ 'main_homog')) %>%
  create_group(drop=FALSE) %>%
  saveRDS('Data/dec-2021-runs/_scenarios.RDS')
