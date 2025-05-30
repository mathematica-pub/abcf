source('programs/utils.R')

file.copy('c:/projects/bcf-sim-infra/data/sim_study_config_medicare_and_edu.xlsx',
          'data/feb-2023-runs/sim_study_config_medicare_and_edu.xlsx',
          overwrite=TRUE)

file.copy('c:/projects/bcf-sim-infra/scripts/feb-2023-small/test_list.RDS',
          'data/feb-2023-runs/test_list.RDS',
          overwrite=TRUE)

scenarios <- readRDS('data/feb-2023-runs/test_list.RDS') %>%
  mutate(set = case_when(nT==500 ~ 'small',
                              nT==1000 ~ 'large'),
         scenario = glue('{world}_sigu{sigu_multiplier}_sigv{sigv_multiplier}_rho{str_replace(rho,"-","m")}_nT{nT}')) %>% 
  select(-seed, -ibcf, -ubcf, -obcf, -fname) %>%
  distinct %T>%
  saveRDS('data/feb-2023-runs/_scenarios.RDS') %>%
  select(scenario, world, set)

from <- 's3://bcf-sim-study/feb-2023-runs/1kburn-1ksim'
to <- 'C:/Projects/bcf-sim-study/Data/feb-2023-runs/1kburn-1ksim'
if (!dir.exists(to)) dir.create(to)
#system(glue('bash -c "aws s3 sync {from} {to} --profile bcf-sims --region us-east-1"'))
#This works in bash itself put not via system(); idk why.
#system(glue('bash -c "aws s3 sync {from} {to} --profile bcf-sims --region us-east-1 --exclude "*" --include "*small*" "'))

#Too big to load all at once, try separating into chunks
files <- tibble(file = dir(to, pattern = '*.RDS', recursive = FALSE, full.names = TRUE)) %>%
  mutate(name = file %>% str_remove('C:/Projects/bcf-sim-study/Data/feb-2023-runs/1kburn-1ksim/') %>% str_remove('\\.RDS'),
         world = str_extract(file,'medicare|edu'),
         scenario = str_match(file,'((medicare|edu)_(.+))_([iwou]bcf)')[,2],
         method = str_extract(file,'[iwou]bcf'),
         sim = as.numeric(str_match(file,'sim([0-9][0-9][0-9])')[,2])) %>%
  left_join(scenarios, by=c('world','scenario'))

smalls <- filter(files, str_detect(file,'small'))
bigs <- filter(files, !str_detect(file,'small'))

combine_output_chunks(smalls, dir='C:/Projects/bcf-sim-study/Data/feb-2023-runs', save_big=FALSE)
