prevalence_main=function(comorbid,
                   summary.dcrt,
                   siteid,
                   dir.repo){

  res.prev=NULL
  for(tt in 1:3){
    for(aa in 1:3){
      for(cc in 1:nrow(comorbid)){
        for(post.period in c(90,180)){
          for(hosp in c(0,1)){
            tryCatch({
              # print(paste0(siteid,"_tt_",tt,"_aa_",aa,"_cc_",cc,"_",post.period))
              res.prev=rbind.data.frame(res.prev,
                                                      prevalence(summary.dcrt,
                                                                         tt,
                                                                         aa,
                                                                         cc,
                                                                         post.period,
                                                                         comorbid,
                                                                         siteid,
                                                                         hosp))
            },error=function(e){NA})

          }
        }
      }
    }
  }

  save(res.prev,
       file=paste0(dir.repo,siteid,"_conditional_testing_results/",siteid,"_prevalence.Rdata"))



}
