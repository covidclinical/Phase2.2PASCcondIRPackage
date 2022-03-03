#' @import dplyr
#' @import stats
logistic_prescreen = function(summary.dcrt,
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

  # select comorbid combo
  pat1=rownames(res.conf.final)[res.conf.final[,colnames(comorbid)[1]]==comorbid[cc,1]]
  pat2=rownames(res.conf.final)[res.conf.final[,colnames(comorbid)[2]]==comorbid[cc,2]]
  pat3=rownames(res.conf.final)[res.conf.final[,colnames(comorbid)[3]]==comorbid[cc,3]]

  pat.keep=as.character(intersect(intersect(intersect(pat1,pat2),pat3),summary.tmp[,"patient_num"]))
  pat.keep=as.character(intersect(pat.keep,rownames(res.out.final)))
  pat.keep=as.character(intersect(pat.keep,rownames(res.conf.final)))
  #print(paste0("strata_size: ",length(pat.keep)))

  if(length(pat.keep)>200 & (0.02*length(pat.keep)<=sum(as.numeric(summary.tmp[pat.keep,"exposure"])))){

    summary.tmp=summary.tmp[pat.keep,]
    res.out.tmp=res.out.final[pat.keep,]
    res.conf.tmp=res.conf.final[pat.keep,]

    #X = as.matrix(res.conf.tmp)
    Z = as.matrix(res.out.tmp)

    prev_Z=apply(Z,MARGIN = 2,mean)
    index.keep.Z = names(prev_Z)[prev_Z>0.01]

    ### Preliminary Screening using Logistic model

    #print("starting logistic screening")
    formula=as.formula(paste0("Y ~A_junk",sep=""))
    A_junk=as.numeric(summary.tmp[,"exposure"])
    index.keep.Z.filter=NULL

    for(zz in 1:length(index.keep.Z)){
      tryCatch({
        index = index.keep.Z[zz]
        index = as.character(index)
        Y=Z[,index]
        junk=cbind.data.frame(Y,A_junk)
        fit.glm=glm(formula,family="binomial",maxit=25,data=junk)
        fit.glm=summary(fit.glm);

        index.keep.Z.filter=rbind.data.frame(index.keep.Z.filter,
                                             cbind.data.frame("siteid"=siteid,
                                                              "phecode"=index,
                                                              "zz"=zz,
                                                              "method"="logistic_marginal",
                                                              "age"=age,
                                                              "period"=tt,
                                                              "hospital_flag"=hosp,
                                                              "post_period"=post.period,
                                                              "comorbid"=paste0("T2D_",comorbid[cc,1],"_obesity_",comorbid[cc,2],"_hyp_",comorbid[cc,3]),
                                                              "beta"=fit.glm$coefficients["A_junk",1],
                                                              "se"=fit.glm$coefficients["A_junk",2],
                                                              "pval"=fit.glm$coefficients["A_junk",4]))

      },error=function(e){NA})
    }


    # save(testing.output,
    #      file=paste0(dir.repo,siteid,"_conditional_testing_results/",
    #                  siteid,"_tt_",tt,"_aa_",aa,"_cc_",cc,"_",post.period,".Rdata"))
    #
    return(index.keep.Z.filter)
  }
} # end of function
