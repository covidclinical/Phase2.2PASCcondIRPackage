#' Run analysis
#'
#' @return NULL. Result files are written out to the `conditional_testing_results` directory.
#'
#' @keywords 4CE
#' @export
#' @param dir.data path of input data (string)
#' @param dir.repo path of output data (string)
#' @import tidyverse
#' @import data.table
#' @import caret
#' @import randomForest
#' @import e1071
#' @import glmnet
#' @import stats
#' @import metafor
#' @import poolr
#'

runAnalysis <- function(dir.data, dir.repo){

  # create output result folder
  dir.create(paste0(dir.repo,siteid,"_conditional_testing_results"))

  # read the data
  obs = fread(paste0(dir.data,"Phase22all_LocalPatientObservations.csv"),stringsAsFactors = F)
  summary = fread(paste0(dir.data,"Phase22all_LocalPatientSummary.csv"),stringsAsFactors = F)

  ### phecode mapping
  #load(sysdata)
  ### Construct phenotype data (takes 15-30 minutes)
  construct_conditional_matrix(dir.repo,
                               siteid,
                               obs,
                               summary)

  ### Load newly constructed phenotype data
  load(paste0(dir.repo, siteid, "_conditional_testing_data_phase22.Rdata"))

  comorbid=expand.grid(0:1,0:1,0:1)
  colnames(comorbid)=c("250.2","278.1","401")
  output = as.list(NULL)

  prevalence_main(comorbid,
                  summary.dcrt,
                  siteid,
                  dir.repo)

  phecode.pass=prescreen(comorbid,
                         summary.dcrt,
                         siteid,
                         dir.repo)

  phecode.pass.dCRT=conditional_testing_dCRT(comorbid,
                                             summary.dcrt,
                                             siteid,
                                             dir.repo,
                                             phecode.pass)

  conditional_testing_DML(comorbid,
                          summary.dcrt,
                          siteid,
                          dir.repo,
                          phecode.pass.dCRT)


}
