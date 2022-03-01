### Libraries needed: dplyr, metafor, poolr

logistic_meta=function(res){
  
  phecode.keep=table(res$phecode)
  phecode.keep=names(phecode.keep)[phecode.keep>1]
  res=data.frame(res)
  res = res %>%
    filter(phecode %in% phecode.keep) 
  res.meta=NULL
  for(pp in phecode.keep){
    tryCatch({
    junk=filter(res,phecode==pp,post_period==90)
    sigmap=fisher(junk$pval)
    res.rma=rma(yi=junk$beta,
                sei=junk$se,
                method="DL")
    res.meta=rbind.data.frame(res.meta,
                                       cbind.data.frame("phecode"=as.character(pp),
                                                        "beta"=res.rma$b,
                                                        "se"=res.rma$se,
                                                        "ci.95L"=res.rma$ci.lb,
                                                        "ci.95U"=res.rma$ci.ub,
                                                        "pval_RE_meta"=res.rma$pval,
                                                        "pval_fisher"=sigmap$p,
                                                        "post_period"=90))
    },error=function(e){NA})
     tryCatch({
    junk1=filter(res,phecode==pp,post_period==180)
    sigmap1=fisher(junk1$pval)
    res.rma1=rma(yi=junk1$beta,
                sei=junk1$se,
                method="DL")
    res.meta=rbind.data.frame(res.meta,
                              cbind.data.frame("phecode"=as.character(pp),
                                               "beta"=res.rma1$b,
                                               "se"=res.rma1$se,
                                               "ci.95L"=res.rma1$ci.lb,
                                               "ci.95U"=res.rma1$ci.ub,
                                               "pval_RE_meta"=res.rma1$pval,
                                               "pval_fisher"=sigmap1$p,
                                               "post_period"=180))
    
    
    },error=function(e){NA})
    
  }
  res.meta$pval_adjust_BH=NA
  res.meta$pval_adjust_BH[res.meta$post_period==90]=p.adjust(res.meta$pval_RE_meta[res.meta$post_period==90],method="BH")
  res.meta$pval_adjust_BH[res.meta$post_period==180]=p.adjust(res.meta$pval_RE_meta[res.meta$post_period==180],method="BH")
  res.meta=res.meta[order(res.meta$pval_adjust_BH),]
  res.meta$PVAL_SIG=res.meta$pval_adjust_BH<0.20
  
  return(res.meta)
}
