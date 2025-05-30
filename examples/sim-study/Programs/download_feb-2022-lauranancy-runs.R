source('programs/utils.R')

scenarios <- readxl::read_xlsx('data/feb-2022-lauranancy-runs/sim_study_config_medicare_and_edu.xlsx', sheet='scenarios') %>%
  mutate(set = case_when(trt_eff_scenario=='homog' ~ 'homog',
                         rho==0 ~ 'norho',
                         rho>0  ~ 'posrho',
                         weights %in% c('F','FALSE') ~ 'unweighted',
                         uv_dist=='t' ~ 'tdist',
                         TRUE ~ 'main'),
         scenario=scenario_num) %>%
  select(scenario, world, set, everything(), -matches('adj'), -scenario_num, -scenario_desc) %T>%
  saveRDS('Data/jan-2022-runs/_scenarios.RDS') %>%
  select(scenario, world, set)

from <- 's3://bcf-sim-study/feb-2022-lauranancy-runs/4kburn-2ksim'
to <- 'C:/Projects/bcf-sim-study/Data/feb-2022-lauranancy-runs/4kburn-2ksim'
#system(glue('bash -c "aws s3 sync {from} {to} --profile bcf-sims --region us-east-1 --exclude "*control-files/*" "'))
#system(glue('bash -c "aws s3 sync {from} {to} --profile bcf-sims --region us-east-1 --exclude "*" --include "*small*" "'))

#Too big to load all at once, try separating into chunks
files <- tibble(file = dir(to, pattern = '*.RDS', recursive = FALSE, full.names = TRUE)) %>%
  mutate(name = file %>% str_remove('C:/Projects/bcf-sim-study/Data/dec-2021-runs/4000burn-10ksim/') %>% str_remove('\\.RDS'),
         world = str_extract(file,'medicare|edu'),
         scenario = as.numeric(str_match(file,'(medicare|edu)_([0-9]+)_(ibcf|wbcf)')[,3]),
         method = str_extract(file,'ibcf|wbcf'),
         sim = as.numeric(str_match(file,'sim([0-9][0-9][0-9])')[,2])) %>%
  left_join(scenarios, by=c('world','scenario'))

smalls <- filter(files, str_detect(file,'small'))
bigs <- filter(files, !str_detect(file,'small'))

#combine_output_chunks(bigs, dir='C:/Projects/bcf-sim-study/Data/feb-2022-lauranancy-runs', save_big=TRUE)
combine_output_chunks(smalls, dir='C:/Projects/bcf-sim-study/Data/feb-2022-lauranancy-runs', save_big=FALSE)
