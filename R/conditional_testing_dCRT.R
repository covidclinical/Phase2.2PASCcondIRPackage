conditional_testing_dCRT = function(comorbid,
                                    summary.dcrt,
                                    siteid,
                                    dir.repo,
                                    phecode.pass){
  
  ### Logistic prescreening
  res.dCRT=NULL
  for(tt in 1:3){
    for(aa in 1:3){
      for(cc in 1:nrow(comorbid)){
        for(post.period in c(90,180)){
          for(hosp in c(0,1)){
            tryCatch({
              print(paste0(siteid,"_tt_",tt,"_aa_",aa,"_cc_",cc,"_",post.period))
              res.dCRT=rbind.data.frame(res.dCRT,
                                                      analysis_dCRT(summary.dcrt,
                                                                         tt,
                                                                         aa,
                                                                         cc,
                                                                         post.period,
                                                                         comorbid,
                                                                         siteid,
                                                                         hosp,
                                                                    as.character(phecode.pass$phecode[phecode.pass$post_period==post.period])))
            },error=function(e){NA})
            
          }
        }
      }
    }
    save(res.dCRT,
         file=paste0(dir.repo,siteid,"_conditional_testing_results/",siteid,"_dCRT_intermediate_","tt_",tt,".Rdata"))
    
  }
  
  save(res.dCRT,
       file=paste0(dir.repo,siteid,"_conditional_testing_results/",siteid,"_dCRT.Rdata"))
  tryCatch({
  res.dCRT.meta=dCRT_meta(res.dCRT)
  
  save(res.dCRT.meta,
       file=paste0(dir.repo,siteid,"_conditional_testing_results/",siteid,"_dCRT_meta.Rdata"))
},error=function(e){NA})
  phecode.pass=res.dCRT.meta %>%
    filter(PVAL_SIG==TRUE) %>%
    select(phecode,post_period,beta,pval_adjust_BH,PVAL_SIG)
  save(phecode.pass,
       file=paste0(dir.repo,siteid,"_conditional_testing_results/",siteid,"_phecode_pass_dCRT",".Rdata"))
  
  return(phecode.pass)
  
  
  
  
  
}