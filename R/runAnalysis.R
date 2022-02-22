#' Run analysis
#'
#' @return NULL. Result files are written out to the `conditional_testing_results` directory.
#'
#' @keywords 4CE
#' @export
#' @import data.table
#' @import caret
#' @import randomForest
#' @import e1071
#' @import gbm
#' @import plyr
#' @import nnet
#' @import stats
#'
runAnalysis <- function(){
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

  ### 90 days testing
  for(tt in 1:3){
    for(aa in 1:3){
      for(cc in 1:nrow(comorbid)){
        print(paste0(siteid,"_tt_",tt,"_aa_",aa,"_cc_",cc,"_90"))
        output[[paste0(siteid,"_tt_",tt,"_aa_",aa,"_cc_",cc,"_90")]]=conditional_testing(summary.dcrt,
                                                                                         tt,
                                                                                         aa,
                                                                                         cc,
                                                                                         time.period=90,
                                                                                         comorbid)
      }
    }
  }

  ### 180 days testing
  for(tt in 1:3){
    for(aa in 1:3){
      for(cc in 1:nrow(comorbid)){
        print(paste0(siteid,"_tt_",tt,"_aa_",aa,"_cc_",cc,"_180"))
        output[[paste0(siteid,"_tt_",tt,"_aa_",aa,"_cc_",cc,"_180")]]=conditional_testing(summary.dcrt,
                                                                                          tt,
                                                                                          aa,
                                                                                          cc,
                                                                                          time.period=180,
                                                                                          comorbid)
      }
    }
  }


  dir.create(paste0(dir.repo,Sys.Date(),"_",siteid,"_conditional_testing_results"))

  save(output,
       file=paste0(dir.repo,Sys.Date(),"_",siteid,"_conditional_testing_results/",siteid,"_conditional_testing_output.Rdata"))
}
