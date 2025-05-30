source('programs/utils.R')

file.copy('c:/projects/bcf-sim-infra/data/sim_study_config_medicare_and_edu.xlsx',
          'data/mar-2023-wcorr/sim_study_config_medicare_and_edu.xlsx',
          overwrite=TRUE)

file.copy('c:/projects/bcf-sim-infra/scripts/mar-2023-wcorr/test_list.RDS',
          'data/mar-2023-wcorr/test_list.RDS',
          overwrite=TRUE)

scenarios <- readRDS('data/mar-2023-wcorr/test_list.RDS') %>%
  bind_rows(readRDS('data/mar-2023-wcorr/test_list_old.RDS')) %>%
  mutate(set = case_when(nT==500 ~ 'small',
                         nT==1000 ~ 'large'),
         scenario = glue('{world}_sigu{sigu_multiplier}_sigv{sigv_multiplier}_rho{str_replace(rho,"-","m")}_nT{nT}_r2wpi{r2_w_pi}_r2wmu{r2_w_mu}_r2wtau{r2_w_tau}')) %>% 
  select(-seed, -ibcf, -ubcf, -obcf, -fname, -ate_prior_sd) %>%
  distinct %T>%
  saveRDS('data/mar-2023-wcorr/_scenarios.RDS') %>%
  select(scenario, world, set)

from <- 's3://bcf-sim-study/mar-2023-wcorr/1kburn-1ksim'
to <- 'C:/Projects/bcf-sim-study/Data/mar-2023-wcorr/1kburn-1ksim'
if (!dir.exists(to)) dir.create(to)
#system(glue('bash -c "aws s3 sync {from} {to} --profile bcf-sims --region us-east-1"'))
#This works in bash itself put not via system(); idk why.
#system(glue('bash -c "aws s3 sync {from} {to} --profile bcf-sims --region us-east-1 --exclude "*" --include "*small*" "'))

#Too big to load all at once, try separating into chunks
files <- tibble(file = dir(to, pattern = '*.RDS', recursive = FALSE, full.names = TRUE)) %>%
  mutate(name = file %>% str_remove('C:/Projects/bcf-sim-study/Data/mar-2023-wcorr/1kburn-1ksim/') %>% str_remove('\\.RDS'),
         world = str_extract(file,'medicare|edu'),
         scenario = str_match(file,'((medicare|edu)_(.+))_([iwou]bcf|stan)')[,2],
         method = str_extract(file,'[iwou]bcf|stan'),
         sim = as.numeric(str_match(file,'sim([0-9][0-9][0-9])')[,2])) %>%
  left_join(scenarios, by=c('world','scenario'))

smalls <- filter(files, str_detect(file,'small'))
bigs <- filter(files, !str_detect(file,'small'))

combine_output_chunks(smalls, dir='C:/Projects/bcf-sim-study/Data/mar-2023-wcorr', save_big=FALSE, skipfirst=35)
#combine_output_chunks(bigs, dir='C:/Projects/bcf-sim-study/Data/mar-2023-wcorr', save_big=TRUE)
