source('programs/utils.R')

acceptance <- array(NA_real_, dim=c(4,4,64,2,200))
dimnames(acceptance) <- list(rownames(msc_list[[1]]$acceptance),
                             colnames(msc_list[[1]]$acceptance),
                             scenarios$scenario,
                             c('iBCF','wBCF'),
                             NULL)
for (i in 1:length(msc_list)) {
  scn_idx <- match(msc_list[[i]]$inputs$scenario, scenarios$scenario)
  sim_idx <- msc_list[[i]]$inputs$sim
  model_idx <- ifelse(msc_list[[i]]$inputs$ibcf,1,2)
  acceptance[,,scn_idx,model_idx,sim_idx] <- msc_list[[i]]$acceptance
}
avg_acc <- apply(acceptance,c(2,3,4),mean)
avg_acc[,,1] %>%
  t %>% as_tibble(rownames='scenario')
