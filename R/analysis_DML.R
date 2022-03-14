#' @import data.table
#' @import dplyr
#' @import caret
#' @import glmnet
#' @import metafor
#' @import poolr
#' @import e1071
#' @import gbm
#' @import nnet
analysis_DML = function(summary.dcrt,
                        tt,
                        aa,
                        cc,
                        post.period,
                        comorbid,
                        siteid,
                        hosp,
                        code.kp,
                        res.out.90.final,
                        res.out.180.final,
                        res.conf.final){
  res.ML=NULL
  if(post.period==90){
    res.out.final=res.out.90.final
  }else if(post.period==180){
    res.out.final=res.out.180.final
  }

  if(aa==1){
    summary.tmp=dplyr::filter(summary.dcrt,
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
  pat1=tryCatch(rownames(res.conf.final)[res.conf.final[,colnames(comorbid)[1]]==comorbid[cc,1]],error=function(e){NA})
  pat2=tryCatch(rownames(res.conf.final)[res.conf.final[,colnames(comorbid)[2]]==comorbid[cc,2]],error=function(e){NA})
  pat3=tryCatch(rownames(res.conf.final)[res.conf.final[,colnames(comorbid)[3]]==comorbid[cc,3]],error=function(e){NA})

  list.pat=Filter(Negate(anyNA),list(pat1,pat2,pat3))
  pat.keep=as.character(intersect(Reduce(intersect, list.pat),summary.tmp[,"patient_num"]))
  pat.keep=as.character(intersect(pat.keep,rownames(res.out.final)))
  pat.keep=as.character(intersect(pat.keep,rownames(res.conf.final)))

  if(length(pat.keep)>200 & (0.02*length(pat.keep)<=sum(as.numeric(summary.tmp[pat.keep,"exposure"])))){

    print(paste0("strata_size: ",length(pat.keep)))

    summary.tmp=summary.tmp[pat.keep,]
    res.out.tmp=res.out.final[pat.keep,]
    res.conf.tmp=res.conf.final[pat.keep,]

    X = (res.conf.tmp)
    Z = as.matrix(res.out.tmp)

    prev_Z=apply(Z,MARGIN = 2,mean)
    index.keep.Z = names(prev_Z)[prev_Z>0.01]
    index.keep.Z=index.keep.Z[index.keep.Z %in% code.kp]

    prev_X <- colMeans(ifelse(X > 0, 1, 0))
    index.keep.X <- which(prev_X > 0.025)
    X <- X[,index.keep.X]
    X=as.matrix(X)

    for(zz in 1:length(index.keep.Z)){

      tryCatch({
        index = index.keep.Z[zz]
        index = as.character(index)

        if(index %in% colnames(X)){
          keep.id=rownames(X)[X[,index]==0]
          X.tmp=X[keep.id,]
          A = as.numeric(summary.tmp[keep.id,"exposure"])
          names(A)=keep.id
          Z.tmp=Z[keep.id,]
        }else{
          X.tmp=X
          A=as.numeric(summary.tmp[,"exposure"])
          names(A)=summary.tmp[,"patient_num"]
          Z.tmp=Z
        }

        ## Third filtering: down-sample Y=0 (1:10)
        id.1 = rownames(Z.tmp)[Z.tmp[,index]==1]
        set.seed(2022)
        id.2 = sample(rownames(Z.tmp)[Z.tmp[,index]==0],min(length(id.1)*5,nrow(Z.tmp)-length(id.1)))
        X.tmp=X.tmp[c(id.1,id.2),]
        Z.tmp=Z.tmp[c(id.1,id.2),]
        A=A[c(id.1,id.2)]

        # Use gbm to estimate the nuisance models:
        K = 2
        t = proc.time()
        method = tryCatch(DML(Z.tmp[,index],as.numeric(A),data.frame(X.tmp),K,givenEstimator = "gbm"), error = function(e){0})
        t = (proc.time() - t)[3]

        meth = method
        if(typeof(meth) == "list"){

          # point estimate of beta:
          betas = tryCatch(Estimate(Z.tmp[,index],as.numeric(A),meth),error= function(e){0})

          # 95% confidence interval of beta by the boostrap:
          bb = Boostrap(Z.tmp[,index],as.numeric(A),meth,B=2000)

          if(betas != 0){
            res.ML=rbind.data.frame(res.ML,
                                    cbind.data.frame("siteid"=siteid,
                                                     "method"="DML",
                                                     "phecode"=index,
                                                     "age"=age,
                                                     "period"=tt,
                                                     "post_period"=post.period,
                                                     "hospital_flag"=hosp,
                                                     "comorbid"=paste0("T2D_",comorbid[cc,1],"_obesity_",comorbid[cc,2],"_hyp_",comorbid[cc,3]),
                                                     "model"="gbm",
                                                     "k"=K,
                                                     "pval"=bb["pval"],
                                                     "beta"=betas,
                                                     "sd"=bb["sd"],
                                                     "ci.95L"=bb["2.5%"],
                                                     "ci.95U"=bb["97.5%"],
                                                     "n"=nrow(Z.tmp)))
          }
        }



      },error=function(e){NA})
    }

    return(res.ML)

  }else(return(NULL)) # if statement


} # end of function




