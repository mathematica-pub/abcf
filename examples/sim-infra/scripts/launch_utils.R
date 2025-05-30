library(tidyverse)
library(glue)

role_arn <- 'fill this in'
train_image <- 'fill this in/bcf-simmer:latest'

split_run_cmds <- function(test_list, 
                           processes,
                           each=processes, 
                           instance='ml.c5.xlarge', 
                           script_path,
                           name_prefix,
                           fun,
                           s3_folder,
                           upload_control=TRUE,
                           volume_size=10) {
  
  n_cmds <- ceiling(nrow(test_list)/each)
  cmds <- lapply(1:n_cmds, function(i) {
    #First save and upload the control file for this batch
    start <- (i-1)*each + 1
    end <- min((i-1)*each + each, nrow(test_list))
    subset <- test_list[start:end,]
    if (upload_control) {
      saveRDS(subset,glue('{script_path}/control_file{i}.RDS'))
      system(glue('bash -c "aws s3 cp {script_path}/control_file{i}.RDS s3://{s3_folder}/control_file{i}.RDS --profile bcf-sims --region us-east-1"'))  
      #file.remove(glue('{script_path}/control_file{i}.RDS'))
    }
    
    #Next save a sagemaker cli args json
    json <- list(TrainingJobName = paste0(name_prefix,i),
                 RoleArn = role_arn,
                 ResourceConfig=list(InstanceType=instance,
                                     InstanceCount=1,
                                     VolumeSizeInGB=volume_size),
                 EnableNetworkIsolation=FALSE,
                 StoppingCondition=list(MaxRuntimeInSeconds=432000),
                 AlgorithmSpecification=list(TrainingImage=train_image, 
                                             TrainingInputMode='File'),
                 InputDataConfig=list(list(ChannelName='train', DataSource=list(S3DataSource=list(S3DataType='S3Prefix',
                                                                                                  S3Uri='s3://bcf-sim-study/testing/useless.txt',
                                                                                                  S3DataDistributionType='FullyReplicated')))),
                 OutputDataConfig=list(S3OutputPath = glue('s3://{s3_folder}')),
                 HyperParameters=list(fun=fun,
                                      processes=as.character(processes),
                                      control_file=glue('s3://{s3_folder}/control_file{i}.RDS')))
    
    json %>%jsonlite::toJSON(auto_unbox=TRUE, pretty=TRUE) %>% writeLines(glue('{script_path}/autojson{i}.json'))
    
    return(glue('bash -c "aws sagemaker create-training-job --cli-input-json file://{script_path}/autojson{i}.json --region us-east-1 --profile bcf-sims"'))
  })
  unlist(cmds)
}

make_ctrl_cmds <- function(test_list, 
                           processes,
                           each=processes, 
                           control_fn,
                           instance='ml.c5.xlarge', 
                           script_path,
                           clear_scripts=TRUE,
                           name_prefix,
                           fun,
                           s3_folder,
                           volume_size=30) {

  if (!dir.exists(script_path)) {
    dir.create(script_path)  
  }
  
  if (clear_scripts) {
    todel <- dir(script_path, pattern = 'autojson', full.names = TRUE)
    lapply(todel, file.remove)
  }
  
  s3_exist <- system(glue('bash -c "aws s3 ls {s3_folder} --region us-east-1 --profile bcf-sims"'))==0
  if (!s3_exist) {
    parent <- str_split_fixed(s3_folder,'/',2)[1,1]
    child <- str_split_fixed(s3_folder,'/',2)[1,2]
    system(glue('bash -c "aws s3api put-object --bucket {parent} --key {child}/ --region us-east-1 --profile bcf-sims"'))
  }
  
  n_cmds <- ceiling(nrow(test_list)/each)
  cmds <- lapply(1:n_cmds, function(i) {
    #First save and upload the control file for this batch
    start <- (i-1)*each + 1
    end <- min((i-1)*each + each, nrow(test_list))
    subset <- test_list[start:end,]
    
    #Next save a sagemaker cli args json
    json <- list(TrainingJobName = paste0(name_prefix,i),
                 RoleArn = role_arn,
                 ResourceConfig=list(InstanceType=instance,
                                     InstanceCount=1,
                                     VolumeSizeInGB=volume_size),
                 EnableNetworkIsolation=FALSE,
                 StoppingCondition=list(MaxRuntimeInSeconds=432000),
                 AlgorithmSpecification=list(TrainingImage=train_image, 
                                             TrainingInputMode='File'),
                 InputDataConfig=list(list(ChannelName='train', DataSource=list(S3DataSource=list(S3DataType='S3Prefix',
                                                                                                  S3Uri='s3://bcf-sim-study/testing/useless.txt',
                                                                                                  S3DataDistributionType='FullyReplicated')))),
                 OutputDataConfig=list(S3OutputPath = glue('s3://{s3_folder}')),
                 HyperParameters=list(fun=fun,
                                      processes  = as.character(processes),
                                      control_fn = as.character(control_fn), 
                                      start      = as.character(start),
                                      end        = as.character(end)
                                      ))
    
    json %>% jsonlite::toJSON(auto_unbox=TRUE, pretty=TRUE) %>% writeLines(glue('{script_path}/autojson{i}.json'))
    
    return(glue('bash -c "aws sagemaker create-training-job --cli-input-json file://{script_path}/autojson{i}.json --region us-east-1 --profile bcf-sims"'))
  })
  unlist(cmds)
}
