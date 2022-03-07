#' @import dplyr
#' @import metafor
#' @import poolr
dCRT_meta=function(res){

  phecode.keep=table(res$phecode)
  phecode.keep=names(phecode.keep)[phecode.keep>1]
  res=data.frame(res)
  res = res %>%
    dplyr::filter(phecode %in% phecode.keep)
  res.meta=NULL
  for(pp in phecode.keep){
    tryCatch({
      junk=dplyr::filter(res,phecode==pp,post_period==90)
      sigmap=fisher(junk$pval)
      res.meta=rbind.data.frame(res.meta,
                                cbind.data.frame("phecode"=as.character(pp),
                                                 "beta"=NA,
                                                 "se"=NA,
                                                 "ci.95L"=NA,
                                                 "ci.95U"=NA,
                                                 "pval_RE_meta"=NA,
                                                 "pval_fisher"=sigmap$p,
                                                 "post_period"=90))
    },error=function(e){NA})
    tryCatch({
      junk1=dplyr::filter(res,phecode==pp,post_period==180)
      sigmap1=fisher(junk1$pval)
      res.meta=rbind.data.frame(res.meta,
                                cbind.data.frame("phecode"=as.character(pp),
                                                 "beta"=NA,
                                                 "se"=NA,
                                                 "ci.95L"=NA,
                                                 "ci.95U"=NA,
                                                 "pval_RE_meta"=NA,
                                                 "pval_fisher"=sigmap1$p,
                                                 "post_period"=180))


    },error=function(e){NA})

  }
  res.meta$pval_adjust_BH=NA
  res.meta$pval_adjust_BH[res.meta$post_period==90]=p.adjust(res.meta$pval_fisher[res.meta$post_period==90],method="BH")
  res.meta$pval_adjust_BH[res.meta$post_period==180]=p.adjust(res.meta$pval_fisher[res.meta$post_period==180],method="BH")
  res.meta=res.meta[order(res.meta$pval_adjust_BH),]
  res.meta$PVAL_SIG=res.meta$pval_adjust_BH<0.20

  return(res.meta)
}
