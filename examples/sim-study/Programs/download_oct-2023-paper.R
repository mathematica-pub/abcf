source('programs/utils.R')

sdir <- 'oct-2023-paper-uhp'
dir.create(glue('data/{sdir}'))
file.copy('c:/projects/bcf-sim-infra/data/sim_study_config_medicare_and_edu.xlsx',
          glue('data/{sdir}/sim_study_config_medicare_and_edu.xlsx'),
          overwrite=TRUE)

file.copy(glue('c:/projects/bcf-sim-infra/scripts/{sdir}/test_list.RDS'),
          glue('data/{sdir}/test_list.RDS'),
          overwrite=TRUE)

scenarios <- readRDS(glue('data/{sdir}/test_list.RDS')) %>%
  mutate(set = case_when(nT==500 ~ 'small',
                         nT==1000 ~ 'large'),
         scenario = glue('{world}_sigu{sigu_multiplier}_sigv{sigv_multiplier}_rho{str_replace(rho,"-","m")}_nT{nT}_siguhp{sigu_pcthyperprior}')) %>% 
  select(-seed, -ibcf, -ubcf, -obcf, -oldbcf, -fname, -ate_prior_sd) %>%
  distinct %T>%
  saveRDS(glue('data/{sdir}/_scenarios.RDS')) %>%
  select(scenario, world, set)

from <- glue('s3://bcf-sim-study/{sdir}/1kburn-1ksim')
to <- glue('C:/Projects/bcf-sim-study/Data/{sdir}/1kburn-1ksim')
if (!dir.exists(to)) dir.create(to)
#system(glue('bash -c "aws s3 sync {from} {to} --profile bcf-sims --region us-east-1"'))
#This works in bash itself put not via system(); idk why.
#system(glue('bash -c "aws s3 sync {from} {to} --profile bcf-sims --region us-east-1 --exclude "*" --include "*small*" "'))

#Too big to load all at once, try separating into chunks
files <- tibble(file = dir(to, pattern = '*.RDS', recursive = FALSE, full.names = TRUE)) %>%
  mutate(name = file %>% str_remove(glue('C:/Projects/bcf-sim-study/Data/{sdir}/1kburn-1ksim/')) %>% str_remove('\\.RDS'),
         world = str_extract(file,'medicare|edu'),
         scenario = str_match(file,'((medicare|edu)_(.+))_(([iwou]|old)bcf|stan)')[,2],
         method = str_extract(file,'([iwou]|old)bcf|stan'),
         sim = as.numeric(str_match(file,'sim([0-9][0-9][0-9])')[,2])) %>%
  left_join(scenarios, by=c('world','scenario')) %>%
  #I created the sim number all wonky; fix
  group_by(scenario, method) %>%
  arrange(sim) %>%
  mutate(sim = row_number()) %>%
  ungroup

smalls <- filter(files, str_detect(file,'small'))
bigs <- filter(files, !str_detect(file,'small'))

combine_output_chunks(smalls, dir=glue('C:/Projects/bcf-sim-study/Data/{sdir}'), save_big=FALSE)
#combine_output_chunks(bigs, dir=glue('C:/Projects/bcf-sim-study/Data/{sdir}'), save_big=TRUE)

sdir <- 'oct-2023-paper-icc'
dir.create(glue('data/{sdir}'))
file.copy('c:/projects/bcf-sim-infra/data/sim_study_config_medicare_and_edu.xlsx',
          glue('data/{sdir}/sim_study_config_medicare_and_edu.xlsx'),
          overwrite=TRUE)

file.copy(glue('c:/projects/bcf-sim-infra/scripts/{sdir}/test_list.RDS'),
          glue('data/{sdir}/test_list.RDS'),
          overwrite=TRUE)

scenarios <- readRDS(glue('data/{sdir}/test_list.RDS')) %>%
  mutate(set = case_when(nT==500 ~ 'small',
                         nT==1000 ~ 'large'),
         scenario = glue('{world}_sigu{sigu_multiplier}_sigv{sigv_multiplier}_rho{str_replace(rho,"-","m")}_nT{nT}_r2e{r2_error}')) %>% 
  select(-seed, -ibcf, -ubcf, -obcf, -oldbcf, -fname, -ate_prior_sd) %>%
  distinct %T>%
  saveRDS(glue('data/{sdir}/_scenarios.RDS')) %>%
  select(scenario, world, set)

from <- glue('s3://bcf-sim-study/{sdir}/1kburn-1ksim')
to <- glue('C:/Projects/bcf-sim-study/Data/{sdir}/1kburn-1ksim')
if (!dir.exists(to)) dir.create(to)
#system(glue('bash -c "aws s3 sync {from} {to} --profile bcf-sims --region us-east-1"'))
#This works in bash itself put not via system(); idk why.
#system(glue('bash -c "aws s3 sync {from} {to} --profile bcf-sims --region us-east-1 --exclude "*" --include "*small*" "'))

#Too big to load all at once, try separating into chunks
files <- tibble(file = dir(to, pattern = '*.RDS', recursive = FALSE, full.names = TRUE)) %>%
  mutate(name = file %>% str_remove(glue('C:/Projects/bcf-sim-study/Data/{sdir}/1kburn-1ksim/')) %>% str_remove('\\.RDS'),
         world = str_extract(file,'medicare|edu'),
         scenario = str_match(file,'((medicare|edu)_(.+))_(([iwou]|old)bcf|stan)')[,2],
         method = str_extract(file,'([iwou]|old)bcf|stan'),
         sim = as.numeric(str_match(file,'sim([0-9][0-9][0-9][0-9]?)')[,2])) %>%
  left_join(scenarios, by=c('world','scenario')) %>%
  #I created the sim number all wonky; fix
  group_by(scenario, method) %>%
  arrange(sim) %>%
  mutate(sim = row_number()) %>%
  ungroup

smalls <- filter(files, str_detect(file,'small'))
bigs <- filter(files, !str_detect(file,'small'))

combine_output_chunks(smalls, dir=glue('C:/Projects/bcf-sim-study/Data/{sdir}'), save_big=FALSE)

sdir <- 'oct-2023-paper-ibcf'
dir.create(glue('data/{sdir}'))
file.copy('c:/projects/bcf-sim-infra/data/sim_study_config_medicare_and_edu.xlsx',
          glue('data/{sdir}/sim_study_config_medicare_and_edu.xlsx'),
          overwrite=TRUE)

file.copy(glue('c:/projects/bcf-sim-infra/scripts/{sdir}/test_list.RDS'),
          glue('data/{sdir}/test_list.RDS'),
          overwrite=TRUE)

scenarios <- readRDS(glue('data/{sdir}/test_list.RDS')) %>%
  mutate(set = case_when(nT==500 ~ 'small',
                         nT==1000 ~ 'large'),
         scenario = glue('{world}_sigu{sigu_multiplier}_sigv{sigv_multiplier}_rho{str_replace(rho,"-","m")}_nT{nT}')) %>% 
  select(-seed, -ibcf, -ubcf, -obcf, -oldbcf, -fname, -ate_prior_sd) %>%
  distinct %T>%
  saveRDS(glue('data/{sdir}/_scenarios.RDS')) %>%
  select(scenario, world, set)

from <- glue('s3://bcf-sim-study/{sdir}/1kburn-1ksim')
to <- glue('C:/Projects/bcf-sim-study/Data/{sdir}/1kburn-1ksim')
if (!dir.exists(to)) dir.create(to)
#system(glue('bash -c "aws s3 sync {from} {to} --profile bcf-sims --region us-east-1"'))
#This works in bash itself put not via system(); idk why.
#system(glue('bash -c "aws s3 sync {from} {to} --profile bcf-sims --region us-east-1 --exclude "*" --include "*small*" "'))

#Too big to load all at once, try separating into chunks
files <- tibble(file = dir(to, pattern = '*.RDS', recursive = FALSE, full.names = TRUE)) %>%
  mutate(name = file %>% str_remove(glue('C:/Projects/bcf-sim-study/Data/{sdir}/1kburn-1ksim/')) %>% str_remove('\\.RDS'),
         world = str_extract(file,'medicare|edu'),
         scenario = str_match(file,'((medicare|edu)_(.+))_(([iwou]|old)bcf|stan)')[,2],
         method = str_extract(file,'([iwou]|old)bcf|stan'),
         sim = as.numeric(str_match(file,'sim([0-9][0-9][0-9])')[,2])) %>%
  left_join(scenarios, by=c('world','scenario')) %>%
  #I created the sim number all wonky; fix
  group_by(scenario, method) %>%
  arrange(sim) %>%
  mutate(sim = row_number()) %>%
  ungroup

smalls <- filter(files, str_detect(file,'small'))
bigs <- filter(files, !str_detect(file,'small'))

combine_output_chunks(smalls, dir=glue('C:/Projects/bcf-sim-study/Data/{sdir}'), save_big=FALSE)
