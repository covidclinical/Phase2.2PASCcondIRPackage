conditional_testing = function(summary.dcrt,
                               tt,
                               aa,
                               cc,
                               time.period,
                               comorbid){
res=NULL
res.ML=NULL
  
if(time.period==90){
  res.out.final=res.out.90.final
}else if(time.period==180){
  res.out.final=res.out.180.final
}

if(aa==1){
  summary.tmp=filter(summary.dcrt,
                     age>18,
                     age<=49,
                     period==tt)
  age="18to49"
}else if(aa==2){
  summary.tmp=filter(summary.dcrt,
                     age>49,
                     age<=69,
                     period==tt)
  age="49to69"
} else{
  summary.tmp=filter(summary.dcrt,
                     age>69,
                     period==tt)
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
print(paste0("strata_size: ",length(pat.keep)))

if(length(pat.keep)>250 & (0.02*length(pat.keep)<=sum(as.numeric(summary.tmp[pat.keep,"exposure"])))){

summary.tmp=summary.tmp[pat.keep,]
res.out.tmp=res.out.final[pat.keep,]
res.conf.tmp=res.conf.final[pat.keep,]
  
X = as.matrix(res.conf.tmp)
Z = as.matrix(res.out.tmp)

prev_Z=apply(Z,MARGIN = 2,mean)
index.keep.Z = names(prev_Z)[prev_Z>0.01]

prev_X <- colMeans(ifelse(X > 0, 1, 0))
index.keep.X <- which(prev_X > 0.025)
X <- X[,index.keep.X]

### Preliminary Screening using Logistic model
print("starting logistic screening")
#formula=as.formula(paste0(" Y ~A_junk+",paste(colnames(junk)[3:ncol(junk)],collapse="+"),sep=""))
formula=as.formula(paste0(" Y ~A_junk",sep=""))
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
if(fit.glm$coefficients["A_junk",4]<0.05){
  index.keep.Z.filter=rbind.data.frame(index.keep.Z.filter,
                                       cbind.data.frame("phecode"=index,
                                                        "zz"=zz,
                                                        "method"="logistic_marginal",
                                                        "age"=age,
                                                        "period"=tt,
                                                        "comorbid"=paste0("T2D_",comorbid[cc,1],"_obesity_",comorbid[cc,2],"_hyp_",comorbid[cc,3]),
                                                        "beta"=fit.glm$coefficients["A_junk",1],
                                                        "se"=fit.glm$coefficients["A_junk",2],
                                                        "pval"=fit.glm$coefficients["A_junk",4]))

      }
  },error=function(e){NA})
}

print('starting confounding analysis')
for(zz in index.keep.Z.filter$zz){
 
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
id.2 = sample(rownames(Z.tmp)[Z.tmp[,index]==0],min(length(id.1)*10,nrow(Z.tmp)-length(id.1)))
X.tmp=X.tmp[c(id.1,id.2),]
Z.tmp=Z.tmp[c(id.1,id.2),]
A=A[c(id.1,id.2)]

tryCatch({
# Fit conditional model for A (getting COVID-19):
# X baseline covariates
Cond_A <- fit_cond_Z(X.tmp, A)
d0CRT_result <- dCRT(A = Z.tmp[,index], Z = A, X = X.tmp, mean_Z = Cond_A, model = 'Binomial_lasso',
                     k = 0, M = 7500,
                     RF.num.trees = c(100,30), MC_free = F, Gen_Z = example_Gen_Z)

res=rbind.data.frame(res,
                     cbind.data.frame("method"="dCRT",
                                      "phecode"=index,
                                      "age"=age,
                                      "period"=tt,
                                      "time.period"=time.period,
                                      "comorbid"=paste0("T2D_",comorbid[cc,1],"_obesity_",comorbid[cc,2],"_hyp_",comorbid[cc,3]),
                                      "model"="Binomial_lasso",
                                      "k"=0,
                                      "pval"=d0CRT_result$pvl))
},error=function(e){NA})
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
  bb = Boostrap(Z.tmp[,index],as.numeric(A),meth,B=3000)
  
  if(betas != 0){
    res.ML=rbind.data.frame(res.ML,
                         cbind.data.frame("method"="DML",
                                          "phecode"=index,
                                          "age"=age,
                                          "period"=tt,
                                          "time.period"=time.period,
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


prev_X=data.frame("phecode"=names(prev_X),
                  "cohort"="all",
                  "setting"="baseline",
                  "age"=age,
                  "period"=tt,
                  "time.period"=time.period,
                  "comorbid"=paste0("T2D_",comorbid[cc,1],"_obesity_",comorbid[cc,2],"_hyp_",comorbid[cc,3]),
                  "prev"=prev_X,
                  "n"=nrow(Z))
prev_Z=data.frame("phecode"=names(prev_Z),
                  "cohort"="all",
                  "setting"="post",
                  "age"=age,
                  "period"=tt,
                  "time.period"=time.period,
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
                        "time.period"=time.period,
                        "comorbid"=paste0("T2D_",comorbid[cc,1],"_obesity_",comorbid[cc,2],"_hyp_",comorbid[cc,3]),
                        "prev"=prev_X_covid,
                        "n"=nrow(Z))
prev_Z_covid=data.frame("phecode"=names(prev_Z_covid),
                        "cohort"="covid",
                        "setting"="post",
                        "age"=age,
                        "period"=tt,
                        "time.period"=time.period,
                        "comorbid"=paste0("T2D_",comorbid[cc,1],"_obesity_",comorbid[cc,2],"_hyp_",comorbid[cc,3]),
                        "prev"=prev_Z_covid,
                        "n"=nrow(Z[id.covid,]))
prev=rbind.data.frame(prev_X,prev_Z,prev_X_covid,prev_Z_covid)

            } # if statement

testing.output=as.list(rep(0,4))
tryCatch({
testing.output[[1]]=index.keep.Z.filter
testing.output[[2]]=res
testing.output[[3]]=res.ML
testing.output[[4]]=prev
},error=function(e){NA})
names(testing.output)=c("logistic","dCRT","DML","prev")
  
return(testing.output)

      } # end of function








