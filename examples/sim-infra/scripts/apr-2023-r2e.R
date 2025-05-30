source('programs/tester.R')  
source('scripts/launch_utils.R')

test_list <- make_apr_2023_r2e_control_file()
n_simul <- 92
depth <- 10

n_machine <- ceiling(nrow(test_list)/n_simul/depth)
print(n_machine)
(nrow(test_list) - (n_simul*depth*(n_machine-1))) / n_simul/depth

#dir.create('scripts/apr-2023-r2e')
saveRDS(test_list,'scripts/apr-2023-r2e/test_list.RDS')

cmds <- make_ctrl_cmds(test_list,
                       processes=n_simul,
                       each=n_simul*depth,
                       fun='aws_run_multiple_make_ctrl',
                       control_fn = 'make_apr_2023_r2e_control_file',
                       instance='ml.m5.24xlarge',
                       script_path='scripts/apr-2023-r2e',
                       name_prefix='20230224-bigm5-apr-2023-r2e-002',
                       s3_folder='bcf-sim-study/apr-2023-r2e/1kburn-1ksim',
                       volume_size=30)

for (i in 1:length(cmds)) {
  print(cmds[i])
  system(cmds[i])
  Sys.sleep(1)
}
