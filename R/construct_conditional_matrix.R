#' @import data.table
#' @import dplyr
#' @import caret
#' @import glmnet
#' @import metafor
#' @import poolr
#' @import e1071
#' @import gbm
#' @import nnet
construct_conditional_matrix = function(dir.repo,
                                        siteid,
                                        obs,
                                        summary){

  ### PheCode mapping
  dat.icd=dplyr::filter(obs,concept_type=="DIAG-ICD10")
  dat.icd=left_join(dat.icd, icd10.phecode.map[,c("concept_code","phecode")],
                    by="concept_code")
  dat.icd=dplyr::filter(dat.icd,!is.na(dat.icd$phecode))
  pat.rm=table(dat.icd$patient_num)
  pat.rm=pat.rm[pat.rm<=10]
  dat.icd=dplyr::filter(dat.icd,!(patient_num %in% names(pat.rm)))
  patients.keep=as.character(unique(dat.icd$patient_num))
  # pat.rm=table(dat.icd$patient_num[dat.icd$days_since_admission<0])
  # pat.rm=pat.rm[pat.rm<=5]
  # dat.icd=filter(dat.icd,!(patient_num %in% names(pat.rm)))

  ### Create data matrices

  ### Confounding matrix X
  icd.tmp=dplyr::filter(dat.icd,days_since_admission<0)
  icd.tmp=dplyr::select(icd.tmp,patient_num,phecode)
  icd.tmp=icd.tmp[!duplicated(icd.tmp),]
  icd.tmp$patient_num=as.character(icd.tmp$patient_num)
  icd.tmp$phecode=as.character(icd.tmp$phecode)
  icd.tmp=table(icd.tmp)
  icd.tmp=data.frame(icd.tmp)
  icd.tmp=arrange(icd.tmp,phecode,patient_num)

  num.patients=length(unique(icd.tmp$patient_num))
  num.patients
  sum(icd.tmp$patient_num[1:num.patients]==icd.tmp$patient_num[(num.patients+1):(num.patients*2)])

  phecodes.conf=as.character(unique(icd.tmp$phecode))
  res.conf=as.list(phecodes.conf)
  counter=0
  for(pp in 1:length(phecodes.conf)){
    res.conf[[pp]]=c(icd.tmp$Freq[((pp-1)*num.patients+1):(pp*num.patients)])
  }
  res.conf.final=do.call("cbind",res.conf)
  colnames(res.conf.final)=phecodes.conf
  rownames(res.conf.final)=icd.tmp$patient_num[1:num.patients]
  pat.add=patients.keep[!(patients.keep%in%rownames(res.conf.final))]
  junk=matrix(rep(0,length(pat.add)*ncol(res.conf.final)),nrow=length(pat.add),ncol=ncol(res.conf.final))
  rownames(junk)=pat.add
  colnames(junk)=colnames(res.conf.final)
  res.conf.final=rbind(res.conf.final,junk)
  health_util = obs %>%
    dplyr::filter(days_since_admission<0) %>%
    dplyr::select(patient_num) %>%
    dplyr::mutate("count"=as.numeric(1)) %>%
    dplyr::group_by(patient_num) %>%
    dplyr::summarize("util"=sum(count))
  xx=as.numeric(log(1+health_util$util))
  health_util[,"util"]=xx
  ##### BUG HERE
  # health_util=rbind.data.frame(health_util,
  #                              cbind.data.frame("patient_num"=rownames(res.conf.final)[!(rownames(res.conf.final) %in% health_util$patient_num)],
  #                                             "util"=0))
  #####
  if(all((rownames(res.conf.final) %in% health_util$patient_num)>0)!=T){
    health_util=rbind.data.frame(health_util,
                                 cbind.data.frame("patient_num"=rownames(res.conf.final)[!(rownames(res.conf.final) %in% health_util$patient_num)],
                                                  "util"=0))

  }

  health_util=as.data.frame(health_util)
  rownames(health_util)=health_util$patient_num
  health_util=health_util[rownames(res.conf.final),]
  health_util=as.data.frame(health_util)
  res.conf.final=as.data.frame(res.conf.final)
  res.conf.final=bind_cols(health_util,
                           res.conf.final)
  res.conf.final=res.conf.final[,-which(colnames(res.conf.final)=="patient_num")]
  res.conf.final[,"util"]=as.numeric(res.conf.final[,"util"])
  #res.conf.final=as.matrix(res.conf.final)

  ### Outcomes Matrix Z 90 days
  icd.tmp=dplyr::filter(dat.icd,days_since_admission>=90)
  icd.tmp=dplyr::select(icd.tmp,patient_num,phecode)
  icd.tmp=icd.tmp[!duplicated(icd.tmp),]
  icd.tmp$patient_num=as.character(icd.tmp$patient_num)
  icd.tmp$phecode=as.character(icd.tmp$phecode)
  icd.tmp=table(icd.tmp)
  icd.tmp=data.frame(icd.tmp)
  icd.tmp=arrange(icd.tmp,phecode,patient_num)

  num.patients=length(unique(icd.tmp$patient_num))
  num.patients
  sum(icd.tmp$patient_num[1:num.patients]==icd.tmp$patient_num[(num.patients+1):(num.patients*2)])

  phecodes.conf=as.character(unique(icd.tmp$phecode))
  res.out.90=as.list(phecodes.conf)
  counter=0
  for(pp in 1:length(phecodes.conf)){
    res.out.90[[pp]]=c(icd.tmp$Freq[((pp-1)*num.patients+1):(pp*num.patients)])
  }
  res.out.90.final=do.call("cbind",res.out.90)
  colnames(res.out.90.final)=phecodes.conf
  rownames(res.out.90.final)=icd.tmp$patient_num[1:num.patients]


  ### Outcomes Matrix Z 180 days
  icd.tmp=dplyr::filter(dat.icd,days_since_admission>=180)
  icd.tmp=dplyr::select(icd.tmp,patient_num,phecode)
  icd.tmp=icd.tmp[!duplicated(icd.tmp),]
  icd.tmp$patient_num=as.character(icd.tmp$patient_num)
  icd.tmp$phecode=as.character(icd.tmp$phecode)
  icd.tmp=table(icd.tmp)
  icd.tmp=data.frame(icd.tmp)
  icd.tmp=arrange(icd.tmp,phecode,patient_num)

  num.patients=length(unique(icd.tmp$patient_num))
  num.patients
  sum(icd.tmp$patient_num[1:num.patients]==icd.tmp$patient_num[(num.patients+1):(num.patients*2)])

  phecodes.conf=as.character(unique(icd.tmp$phecode))
  res.out.180=as.list(phecodes.conf)
  counter=0
  for(pp in 1:length(phecodes.conf)){
    res.out.180[[pp]]=c(icd.tmp$Freq[((pp-1)*num.patients+1):(pp*num.patients)])
  }
  res.out.180.final=do.call("cbind",res.out.180)
  colnames(res.out.180.final)=phecodes.conf
  rownames(res.out.180.final)=icd.tmp$patient_num[1:num.patients]

  ### Outcomes Matrix Z 30 days
  # icd.tmp=filter(dat.icd,days_since_admission>=30)
  # icd.tmp=select(icd.tmp,patient_num,phecode)
  # icd.tmp=icd.tmp[!duplicated(icd.tmp),]
  # icd.tmp$patient_num=as.character(icd.tmp$patient_num)
  # icd.tmp$phecode=as.character(icd.tmp$phecode)
  # icd.tmp=table(icd.tmp)
  # icd.tmp=data.frame(icd.tmp)
  # icd.tmp=arrange(icd.tmp,phecode,patient_num)
  #
  # num.patients=length(unique(icd.tmp$patient_num))
  # num.patients
  # sum(icd.tmp$patient_num[1:num.patients]==icd.tmp$patient_num[(num.patients+1):(num.patients*2)])
  #
  # phecodes.conf=as.character(unique(icd.tmp$phecode))
  # res.out.30=as.list(phecodes.conf)
  # counter=0
  # for(pp in 1:length(phecodes.conf)){
  #   res.out.30[[pp]]=c(icd.tmp$Freq[((pp-1)*num.patients+1):(pp*num.patients)])
  # }
  # res.out.30.final=do.call("cbind",res.out.30)
  # colnames(res.out.30.final)=phecodes.conf
  # rownames(res.out.30.final)=icd.tmp$patient_num[1:num.patients]



  ### Exposure vector A
  junk=obs %>%
    dplyr::filter(concept_type=="COVID-TEST",
           concept_code=="covidpos") %>%
    dplyr::select(patient_num,concept_type,concept_code)
  patients.keep=unique(junk$patient_num)
  exposure.tmp=ifelse(summary$patient_num %in% patients.keep,1,0)
  summary$exposure=exposure.tmp
  summary=dplyr::filter(summary,
                 !(grepl("U07",cohort)&exposure==0))
  summary$admission_date=as.Date(summary$admission_date)
  time.period=unlist(lapply(summary$admission_date,FUN=function(ll){
    if(as.Date(ll)<=as.Date("2020-08-31")){
      return(1)
    }else if(as.Date(ll)>as.Date("2020-08-31") & as.Date(ll)<as.Date("2021-05-30")){
      return(2)
    }else if(as.Date(ll)>=as.Date("2021-05-30") & as.Date(ll)<as.Date("2022-10-31")){
      return(3)
    }
  }))
  summary$period=time.period
  hosp.flag=ifelse(grepl("NotAdm",summary$cohort),0,1)
  summary$hospital_flag=hosp.flag
  summary.dcrt=summary

  ### Save Data
  save(res.conf.final,
       res.out.180.final,
       res.out.90.final,
       summary.dcrt,
       file=paste0(dir.repo,siteid,"_conditional_testing_data_phase22.Rdata"))

}


