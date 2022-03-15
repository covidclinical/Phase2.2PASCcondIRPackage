#' Run analysis
#'
#' @return NULL. Result files are written out to the `conditional_testing_results` directory.
#'
#' @keywords 4CE
#' @export
#' @param dir.data path of input data (string)
#' @param dir.repo path of output data (string)
#' @param siteid your 4CE site id (string)
#' @param run.DML whether or not to run GBM model, TRUE by default
#' @import data.table
#' @import dplyr
#' @import caret
#' @import glmnet
#' @import metafor
#' @import poolr
#' @import e1071
#' @import gbm
#' @import nnet

# files=list.files(paste0("~/Documents/GitHub/Phase2.2PASCcondIRPackage/R/"))
# files=files[grepl(".R",files)]
# invisible(sapply(paste0("~/Documents/GitHub/Phase2.2PASCcondIRPackage/R/",files), function(x) tryCatch(source(x),error=function(e) NA)))


runAnalysis <- function(dir.data, dir.repo, siteid, run.DML=T){

  # create output result folder
   dir.create(paste0(dir.repo,siteid,"_conditional_testing_results"))
   dir.create(paste0(dir.repo,siteid,"_conditional_testing_results/",siteid,"_DML_null_distribution"))

  # read the data
   obs = fread(paste0(dir.data,"LocalPatientObservations.csv"),stringsAsFactors = F)
   summary = fread(paste0(dir.data,"LocalPatientSummary.csv"),stringsAsFactors = F)

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

  prevalence_main(comorbid,
                  summary.dcrt,
                  siteid,
                  dir.repo,
                  res.out.90.final,
                  res.out.180.final,
                  res.conf.final)

  phecode.pass=prescreen(comorbid,
                         summary.dcrt,
                         siteid,
                         dir.repo,
                         res.out.90.final,
                         res.out.180.final,
                         res.conf.final)

  phecode.pass.dCRT=conditional_testing_dCRT(comorbid,
                                             summary.dcrt,
                                             siteid,
                                             dir.repo,
                                             phecode.pass,
                                             res.out.90.final,
                                             res.out.180.final,
                                             res.conf.final)
  if(run.DML==TRUE){
    conditional_testing_DML(comorbid,
                            summary.dcrt,
                            siteid,
                            dir.repo,
                            phecode.pass.dCRT,
                            res.out.90.final,
                            res.out.180.final,
                            res.conf.final)

  }

  res.conf.final=NULL
  res.out.180.final=NULL
  res.out.90.final=NULL
  summary.dcrt=NULL
  save(res.conf.final,
       res.out.180.final,
       res.out.90.final,
       summary.dcrt,
       file=paste0(dir.repo,siteid,"_conditional_testing_data_phase22.Rdata"))


}
