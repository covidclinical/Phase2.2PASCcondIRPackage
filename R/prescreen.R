#' @import data.table
#' @import dplyr
#' @import caret
#' @import glmnet
#' @import metafor
#' @import poolr
#' @import e1071
#' @import gbm
#' @import nnet

prescreen=function(comorbid,
                   summary.dcrt,
                   siteid,
                   dir.repo){

  ### Logistic prescreening
  res.logistic.prescreen=NULL
  for(tt in 1:3){
    for(aa in 1:3){
      for(cc in 1:nrow(comorbid)){
        for(post.period in c(90,180)){
          for(hosp in c(0,1)){
            tryCatch({
              print(paste0(siteid,"_tt_",tt,"_aa_",aa,"_cc_",cc,"_",post.period))
              res.logistic.prescreen=rbind.data.frame(res.logistic.prescreen,
                                                      logistic_prescreen(summary.dcrt,
                                                                         tt,
                                                                         aa,
                                                                         cc,
                                                                         post.period,
                                                                         comorbid,
                                                                         siteid,
                                                                         hosp,
                                                                         res.out.90.final,
                                                                         res.out.180.final,
                                                                         res.conf.final))
            },error=function(e){NA})

          }
        }
      }
    }
  }

  save(res.logistic.prescreen,
       file=paste0(dir.repo,siteid,"_conditional_testing_results/",siteid,"_logistic_prescreen.Rdata"))

  res.logistic.prescreen.meta=logistic_meta(res.logistic.prescreen)

  save(res.logistic.prescreen.meta,
       file=paste0(dir.repo,siteid,"_conditional_testing_results/",siteid,"_logistic_prescreen_meta.Rdata"))

  phecode.pass=res.logistic.prescreen.meta %>%
    dplyr::filter(PVAL_SIG==TRUE) %>%
    dplyr::select(phecode,post_period,beta,pval_adjust_BH,PVAL_SIG)

  save(phecode.pass,
       file=paste0(dir.repo,siteid,"_conditional_testing_results/",siteid,"_phecode_pass.Rdata"))

  return(phecode.pass)


}
