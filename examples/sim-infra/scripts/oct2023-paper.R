source('programs/tester.R')  
source('scripts/launch_utils.R')

test_list <- make_oct2023_paper_uhp_ctrl()
n_simul <- 80
depth <- 5

n_machine <- ceiling(nrow(test_list)/n_simul/depth)
print(n_machine)
(nrow(test_list) - (n_simul*depth*(n_machine-1))) / n_simul/depth

sdir <- 'oct-2023-paper-uhp'
dir.create(glue('scripts/{sdir}'))
saveRDS(test_list,glue('scripts/{sdir}/test_list.RDS'))

cmds <- make_ctrl_cmds(test_list,
                       processes=n_simul,
                       each=n_simul*depth,
                       fun='aws_run_multiple_make_ctrl',
                       control_fn = 'make_oct2023_paper_uhp_ctrl',
                       instance='ml.m5.24xlarge',
                       script_path=glue('scripts/{sdir}'),
                       name_prefix=glue('20231010-bigm5-{sdir}-001'),
                       s3_folder=glue('bcf-sim-study/{sdir}/1kburn-1ksim'),
                       volume_size=30)

for (i in 1:length(cmds)) {
  print(cmds[i])
  system(cmds[i])
  Sys.sleep(1)
}

test_list <- make_oct2023_paper_icc_ctrl()
n_simul <- 80
depth <- 5

n_machine <- ceiling(nrow(test_list)/n_simul/depth)
print(n_machine)
(nrow(test_list) - (n_simul*depth*(n_machine-1))) / n_simul/depth

sdir <- 'oct-2023-paper-icc'
dir.create(glue('scripts/{sdir}'))
saveRDS(test_list,glue('scripts/{sdir}/test_list.RDS'))

cmds <- make_ctrl_cmds(test_list,
                       processes=n_simul,
                       each=n_simul*depth,
                       fun='aws_run_multiple_make_ctrl',
                       control_fn = 'make_oct2023_paper_icc_ctrl',
                       instance='ml.m5.24xlarge',
                       script_path=glue('scripts/{sdir}'),
                       name_prefix=glue('20231010-bigm5-{sdir}-001'),
                       s3_folder=glue('bcf-sim-study/{sdir}/1kburn-1ksim'),
                       volume_size=30)

for (i in 1:length(cmds)) {
  print(cmds[i])
  system(cmds[i])
  Sys.sleep(1)
}

test_list <- make_oct2023_paper_ibcf()
n_simul <- 80
depth <- 5

n_machine <- ceiling(nrow(test_list)/n_simul/depth)
print(n_machine)
(nrow(test_list) - (n_simul*depth*(n_machine-1))) / n_simul/depth

sdir <- 'oct-2023-paper-ibcf'
dir.create(glue('scripts/{sdir}'))
saveRDS(test_list,glue('scripts/{sdir}/test_list.RDS'))

cmds <- make_ctrl_cmds(test_list,
                       processes=n_simul,
                       each=n_simul*depth,
                       fun='aws_run_multiple_make_ctrl',
                       control_fn = 'make_oct2023_paper_ibcf',
                       instance='ml.m5.24xlarge',
                       script_path=glue('scripts/{sdir}'),
                       name_prefix=glue('20231010-bigm5-{sdir}-001'),
                       s3_folder=glue('bcf-sim-study/{sdir}/1kburn-1ksim'),
                       volume_size=30)

for (i in 1:length(cmds)) {
  print(cmds[i])
  system(cmds[i])
  Sys.sleep(1)
}
