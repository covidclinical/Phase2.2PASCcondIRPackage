#' @import dplyr
prevalence = function(summary.dcrt,
                      tt,
                      aa,
                      cc,
                      post.period,
                      comorbid,
                      siteid,
                      hosp){

  if(post.period==90){
    res.out.final=res.out.90.final
  }else if(post.period==180){
    res.out.final=res.out.180.final
  }

  if(aa==1){
    summary.tmp=filter(summary.dcrt,
                       age>18,
                       age<=49,
                       period==tt,
                       hospital_flag==hosp)
    age="18to49"
  }else if(aa==2){
    summary.tmp=filter(summary.dcrt,
                       age>49,
                       age<=69,
                       period==tt,
                       hospital_flag==hosp)
    age="49to69"
  } else{
    summary.tmp=filter(summary.dcrt,
                       age>69,
                       period==tt,
                       hospital_flag==hosp)
    age="69plus"
  }
  summary.tmp=as.matrix(summary.tmp)
  rownames(summary.tmp)=summary.tmp[,"patient_num"]

  # select comorbid combo\
  ######### BUG HERE
  # pat1=rownames(res.conf.final)[res.conf.final[,colnames(comorbid)[1]]==comorbid[cc,1]]
  # pat2=rownames(res.conf.final)[res.conf.final[,colnames(comorbid)[2]]==comorbid[cc,2]]
  # pat3=rownames(res.conf.final)[res.conf.final[,colnames(comorbid)[3]]==comorbid[cc,3]]

  # pat.keep=as.character(intersect(intersect(intersect(pat1,pat2),pat3),summary.tmp[,"patient_num"]))
  # pat.keep=as.character(intersect(pat.keep,rownames(res.out.final)))
  # pat.keep=as.character(intersect(pat.keep,rownames(res.conf.final)))
  # print(paste0("strata_size: ",length(pat.keep)))
  #########
  pat1=tryCatch(rownames(res.conf.final)[res.conf.final[,colnames(comorbid)[1]]==comorbid[cc,1]],error=function(e){NA})
  pat2=tryCatch(rownames(res.conf.final)[res.conf.final[,colnames(comorbid)[2]]==comorbid[cc,2]],error=function(e){NA})
  pat3=tryCatch(rownames(res.conf.final)[res.conf.final[,colnames(comorbid)[3]]==comorbid[cc,3]],error=function(e){NA})

  list.pat=Filter(Negate(anyNA),list(pat1,pat2,pat3))
  pat.keep=as.character(intersect(Reduce(intersect, list.pat),summary.tmp[,"patient_num"]))
  pat.keep=as.character(intersect(pat.keep,rownames(res.out.final)))
  pat.keep=as.character(intersect(pat.keep,rownames(res.conf.final)))
  print(paste0("strata_size: ",length(pat.keep)))

  ######## BUG HERE else?
  if(length(pat.keep)>200 & (0.02*length(pat.keep)<=sum(as.numeric(summary.tmp[pat.keep,"exposure"])))){

    summary.tmp=summary.tmp[pat.keep,]
    res.out.tmp=res.out.final[pat.keep,]
    res.conf.tmp=res.conf.final[pat.keep,]

    X = (res.conf.tmp)
    Z = as.matrix(res.out.tmp)


    prev_Z=apply(Z,MARGIN = 2,mean)
    index.keep.Z = which(prev_Z > 0.01)
    Z=Z[,index.keep.Z]

    X=X[,-which(colnames(X)=="util")]
    prev_X <- colMeans(ifelse(X > 0, 1, 0))
    index.keep.X <- which(prev_X > 0.01)
    X=X[,index.keep.X]

    prev_Z=apply(Z,MARGIN = 2,mean)
    prev_X <- colMeans(ifelse(X > 0, 1, 0))

    prev_X=data.frame("phecode"=names(prev_X),
                      "cohort"="all",
                      "setting"="baseline",
                      "age"=age,
                      "period"=tt,
                      "post_period"=post.period,
                      "hospital_flag"=hosp,
                      "comorbid"=paste0("T2D_",comorbid[cc,1],"_obesity_",comorbid[cc,2],"_hyp_",comorbid[cc,3]),
                      "prev"=prev_X,
                      "n"=nrow(Z))
    prev_Z=data.frame("phecode"=names(prev_Z),
                      "cohort"="all",
                      "setting"="post",
                      "age"=age,
                      "period"=tt,
                      "post_period"=post.period,
                      "hospital_flag"=hosp,
                      "comorbid"=paste0("T2D_",comorbid[cc,1],"_obesity_",comorbid[cc,2],"_hyp_",comorbid[cc,3]),
                      "prev"=prev_Z,
                      "n"=nrow(Z))
    id.covid=as.character(summary.tmp[,"patient_num"][summary.tmp[,"exposure"]==1])
    prev_Z_covid=apply(Z[id.covid,],MARGIN = 2,mean)
    prev_X_covid =apply(X[id.covid,],MARGIN=2,mean)


    prev_X_covid=data.frame("phecode"=names(prev_X_covid),
                            "cohort"="covid",
                            "setting"="baseline",
                            "age"=age,
                            "period"=tt,
                            "post_period"=post.period,
                            "hospital_flag"=hosp,
                            "comorbid"=paste0("T2D_",comorbid[cc,1],"_obesity_",comorbid[cc,2],"_hyp_",comorbid[cc,3]),
                            "prev"=prev_X_covid,
                            "n"=nrow(Z))
    prev_Z_covid=data.frame("phecode"=names(prev_Z_covid),
                            "cohort"="covid",
                            "setting"="post",
                            "age"=age,
                            "period"=tt,
                            "post_period"=post.period,
                            "hospital_flag"=hosp,
                            "comorbid"=paste0("T2D_",comorbid[cc,1],"_obesity_",comorbid[cc,2],"_hyp_",comorbid[cc,3]),
                            "prev"=prev_Z_covid,
                            "n"=nrow(Z[id.covid,]))
    prev=rbind.data.frame(prev_X,prev_Z,prev_X_covid,prev_Z_covid)
    prev=cbind.data.frame("siteid"=siteid,
                          prev)

    return(prev)
  }else{print("Sample size is too small")}

}
