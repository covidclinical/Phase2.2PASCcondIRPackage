conditional_testing_DML = function(comorbid,
                                    summary.dcrt,
                                    siteid,
                                    dir.repo,
                                    phecode.pass){


  res.DML=NULL
  for(tt in 1:3){
    for(aa in 1:3){
      for(cc in 1:nrow(comorbid)){
        for(post.period in c(90,180)){
          for(hosp in c(0,1)){
            tryCatch({
              print(paste0(siteid,"_tt_",tt,"_aa_",aa,"_cc_",cc,"_",post.period))
              res.DML=rbind.data.frame(res.DML,
                                       analysis_DML(summary.dcrt,
                                                   tt,
                                                   aa,
                                                   cc,
                                                   post.period,
                                                   comorbid,
                                                   siteid,
                                                   hosp,
                                                  code.kp=as.character(phecode.pass$phecode[phecode.pass$post_period==post.period]),
                                                  res.out.90.final,
                                                  res.out.180.final,
                                                  res.conf.final))
            },error=function(e){NA})

          }
        }
      }
      save(res.DML,
           file=paste0(dir.repo,siteid,"_conditional_testing_results/",siteid,"_DML_intermediate_","tt_",tt,"_aa_",aa,".Rdata"))
    }
      }

  save(res.DML,
       file=paste0(dir.repo,siteid,"_conditional_testing_results/",siteid,"_DML.Rdata"))

}



