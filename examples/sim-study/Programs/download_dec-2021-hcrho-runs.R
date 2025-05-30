source('programs/utils.R')

from <- 's3://bcf-sim-study/dec-2021-hcrho-runs/4000burn-10ksim'
to <- 'C:/Projects/bcf-sim-study/Data/dec-2021-hcrho-runs/4000burn-10ksim'
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
         method = 'hbcf')

combine_output_chunks(files, dir='C:/Projects/bcf-sim-study/Data/dec-2021-hcrho-runs')
