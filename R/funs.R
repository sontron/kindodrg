#' fn_ADRGRulesLst
#' 
#' transform ADRGRulesDf to ageRuleDf,bornWtRuleDf,PDXRuleDf,PPROCRuleDf objs.
fn_ADRGRulesLst<-function(adrgrules,rulevars=c('age','bornWt','statusOut','LOS','PDX','PPROC','CC')){
  require('stringi')
  #require('parallel') # 并行计算
  
  ADRGRulesLst<-list()
  ADRGRulesLst[["ADRGRulesDf"]]<-adrgrules
  
  for(i in c('age','bornWt','statusOut','LOS','PDX','PPROC','CC')){
    
    adrgrules[,i]->ruleRaw
    stri_replace_all_fixed(ruleRaw,";","|")->ruleRegex
    stri_replace_all_regex(ruleRegex,"[\\(\\)\\[\\]]","")->ruleelements
    unique(unlist(stri_split_fixed(ruleelements,"|")))->RuleElements
    which(RuleElements%in%c('NULL','NtAvb'))->ind
    c(RuleElements[-ind],'NtAvb')->NameRule
    NameRule[nchar(NameRule)>=2]->NameRule
    
    ind<-numeric(length(NameRule))
    for(j in 1:length(NameRule)){
      as.vector(unlist(sapply(ruleRegex,function(k){
        if(k=='[NULL]'){
          T->res
        } else {
          if(stri_detect_fixed(k,"(")){
            stri_detect_regex(NameRule[j],k)->res
          }
          if(stri_detect_fixed(k,"[")){
            stri_replace_all_fixed(k,"[","(")->k
            stri_replace_all_fixed(k,"]",")")->k
            !stri_detect_regex(NameRule[j],k)->res
          }
        }
        return(res)
      }
      )))->res
      paste(which(res),collapse="|")->ind[j]
    }
    
    dt<-data.frame(NameRule,ind,stringsAsFactors = F)
    names(dt)<-c(i,paste(i,'IndADRGRule',sep=''))
    
    ADRGRulesLst[[paste(i,"RuleDf",sep='')]]<-dt
    
  }
  return(ADRGRulesLst)
}


fn_null<-function(){
  #fn_ADRGRulesLst(adrgrules)->res
  res$PDXRuleDf[-which(res$PDXRuleDf=='NtAvb'),]->res$PDXRuleDf
  res$LOSRuleDf[2,1]<-'5-Inf'
  res$ageRuleDf[4,1]<-'17-Inf'
  #res$ageRuleDf[4,1]<-'17-Inf'
  res->CNDRGRulesLst20170619
  save(CNDRGRulesLst20170619,file='d:/github/kindodrg/data/CNDRGRulesLst20170619.RData')
  
}



#' fn_multiVars
#' 
#' this is a inner function in ADRGGrouper for multiple decisions
fn_multiVars<-function(x,y){
  if(y=='[NULL]'){
    T->res
  } else {
    if(grepl("(",y,fixed=T)){
      res<-any(stri_detect_regex(x,y))
    }
    if(grepl("[",y,fixed=T)){
      stri_replace_all_fixed(y,'[','(')->y
      stri_replace_all_fixed(y,']',')')->y
      res<-all(!stri_detect_regex(x,y))
    }
  }
  return(res)
}




#' fn_Trans
#' 
#' transform ICD10 and ICD9 code to match DRGRulesLst objs
#' 
#' @export

fn_Trans<-function(data,id='id',drgruleslst=CNDRGRulesLst20170619){
  require(stringi)
  require(stronmisc)
  as.data.frame(data,stringsAsFactors=F)->data
  data[data=='']<-NA
  drgruleslst[['PDXRuleDf']][,1]->PDX
  drgruleslst[['PPROCRuleDf']][,1]->PPROC
  
  unique(data$PDX)->pdx
  unique(data$PROC1)->pproc
  
  data.frame(pdx=pdx,stringsAsFactors = F)->dtpdx
  data.frame(pproc=pproc,stringsAsFactors = F)->dtpproc
  
  data.frame(PDX=PDX,stringsAsFactors = F)->dtPDX
  data.frame(PPROC=PPROC,stringsAsFactors = F)->dtPPROC
  
  fn_trans(dataFrom=dtpdx,dataTo=dtPDX,byFrom='pdx',byTo='PDX',no.res=1)->distpdx
  fn_trans(dataFrom=dtpproc,dataTo=dtPPROC,byFrom='pproc',byTo='PPROC',no.res=1)->distpproc
  
  distpdx[,c(1,3)]->distpdx
  distpproc[,c(1,3)]->distpproc
  names(distpdx)[2]<-'matchedPDX'
  names(distpproc)[2]<-'matchedPPROC'
  
  merge(data,distpdx,by.x='PDX',by.y='pdx',all.x=T)->data
  data$PDX<-data$matchedPDX
  
  merge(data,distpproc,by.x='PROC1',by.y='pproc',all.x=T)->data
  data$PROC1<-data$matchedPPROC
  
  
  #data[,c('PDX',paste("ADX",1:15,sep=''))]->dx_mat
  #unique(as.vector(as.matrix(dx_mat)))->dx
  
  #data[,paste("PROC",1:8,sep='')]->proc_mat
  #unique(as.vector(as.matrix(proc_mat)))->proc
  data[,-which(names(data)%in%c('matchedPDX','matchedPPROC'))]->data
  return(data)
}


#' fn_ccexcl
#' 
#' test
#' 
#' @export
fn_ccexcl<-function(data=NA,id='id',dx='dx'){
  as.data.frame(data,stringsAsFactors=F)->data
  CCEcl->cc_excl
  require("stringi")
  dx_mat<-as.matrix(data[,c('PDX',paste('ADX',1:15,sep=''))])
  pdx<-as.vector(as.character(data$PDX))
  adx<-as.matrix(data[,c(paste('ADX',1:15,sep=''))])
  sapply(1:length(pdx),function(i){
    which(stri_detect_fixed(pdx[i],cc_excl[,1]))[1]->ind.cc
    if(length(ind.cc)==0) {
      return(as.vector(dx_mat[i,]))} else {
        cc_excl[ind.cc,2]->adx_excl
        adxx<-as.vector(adx[i,])
        ifelse(stri_detect_regex(adxx,adx_excl),NA,adxx)->adxx
        return(c(pdx[i],adxx))
      }
    
  })->dx_correct
  data[,c('PDX',paste('ADX',1:15,sep=''))]<-t(dx_correct)
  return(data)
  
}




#' fn_cc
#' 
#' test
#' 
#' @export
fn_cc<-function(data,id='id'){
  require('stringi')
  tryCatch(id<-data[,id],error=function(e)id<<-NA)
  if(all(is.na(id))) id<-1:nrow(data) else id
  
  
  adx<-as.matrix(data[,c(paste('ADX',1:15,sep=''))])
  CCStatus->cc_mcc
  cc_mcc$cc->CC
  cc_mcc$adx->ADX
  names(data)->nameData
  
  
  
  for(i in 1:15){
    substr(data[,paste('ADX',i,sep='')],1,5)->data[,paste('subADX',i,sep='')]
    merge(data,cc_mcc,by.x=paste('subADX',i,sep=''),by.y='adx',all.x=T,sort=F)->data
    which(names(data)=='cc')->ind
    names(data)[ind]<-paste('cc',i,sep='')
  }
  
  apply(data[,paste('cc',1:15,sep='')],1,function(i){
    ifelse(any(i=='MCC',na.rm=T),'MCC',ifelse(any(i=='CC',na.rm=T),'CC','noCC'))->CC
    return(CC)
  })->data$CC
  
  return(data[,c(nameData,'CC')])
}



#' CNDRGGrouper
#' 
#' CNDRGGrouper takes data(which is transformed by dataTrans into standard data) and return ADRG results containing MDC code, MDC name, ADRG code, ADRG name.
CNDRGGrouper<-function(data,ADRGRulesLst,groupVars=c('age','bornWt','statusOut','LOS','PDX','PPROC','ADX','APROC','CC'),weight=c(1,2,4,8,16,32,64,128,256)){
  require('data.table')
  require('stringi')
  
  ADRGRulesLst[["ADRGRulesDf"]]->ADRGRulesDf
  weightMt<-matrix(nc=length(weight),nr=nrow(ADRGRulesDf))
  for(i in 1:ncol(weightMt)){
    weightMt[,i]<-weight[i]
  }
  
  names(data)->namesData
  
  ifelse(ADRGRulesDf$age=='[NULL]',0,weightMt[,1])->weightMt[,1]
  ifelse(ADRGRulesDf$bornWt=='[NULL]',0,weightMt[,2])->weightMt[,2]
  
  ifelse(ADRGRulesDf$statusOut=='[NULL]',0,weightMt[,3])->weightMt[,3]
  ifelse(ADRGRulesDf$LOS=='[NULL]',0,weightMt[,4])->weightMt[,4]
  ifelse(ADRGRulesDf$PDX=='[NULL]',0,weightMt[,5])->weightMt[,5]
  ifelse(ADRGRulesDf$PPROC=='[NULL]',0,weightMt[,6])->weightMt[,6]
  ifelse(ADRGRulesDf$ADX=='[NULL]',0,weightMt[,7])->weightMt[,7]
  ifelse(ADRGRulesDf$APROC=='[NULL]',0,weightMt[,8])->weightMt[,8]
  ifelse(ADRGRulesDf$CC=='[NULL]',0,weightMt[,9])->weightMt[,9]
  
  ifelse(data$age<=28/365,'0-28d',ifelse(data$age<=1,'29d-1',ifelse(data$age<=17,'1-17','17-Inf')))->data$ageGroup
  
  data$bornWtGroup<-ifelse(data$bornWt<1500,'0-1500',
                           ifelse(data$bornWt<2500,'1500-2500','2500-Inf'))
  
  data$LOSGroup<-ifelse(data$LOS<5,'0-5','5-Inf')
  
  data[is.na(data)]<-'NtAvb'
  as.data.table(data)->data
  
  as.data.table(ADRGRulesLst[['ageRuleDf']])->ageRuleDf
  as.data.table(ADRGRulesLst[['bornWtRuleDf']])->bornWtRuleDf
  as.data.table(ADRGRulesLst[['LOSRuleDf']])->LOSRuleDf
  as.data.table(ADRGRulesLst[['PDXRuleDf']])->PDXRuleDf
  as.data.table(ADRGRulesLst[['PPROCRuleDf']])->PPROCRuleDf
  as.data.table(ADRGRulesLst[['CCRuleDf']])->CCRuleDf
  as.data.table(ADRGRulesLst[['statusOutRuleDf']])->statusOutRuleDf
  
  
  as.data.frame(PPROCRuleDf,stringsAsFactors=F)->PPROCRuleDf
  unique(PPROCRuleDf[,'PPROC'])->procs
  sapply(procs,function(i)which(PPROCRuleDf[,'PPROC']==i)[1])->ind
  as.data.table(PPROCRuleDf[ind,])->PPROCRuleDf
  
  
  merge(data,ageRuleDf,by.x='ageGroup',by.y='age',all.x=T)->data
  merge(data,bornWtRuleDf,by.x='bornWtGroup',by.y='bornWt',all.x=T)->data
  merge(data,LOSRuleDf,by.x='LOSGroup',by.y='LOS',all.x=T)->data
  merge(data,PDXRuleDf,by.x='PDX',by.y='PDX',all.x=T)->data
  merge(data,PPROCRuleDf,by.x='PROC1',by.y='PPROC',all.x=T)->data
  merge(data,CCRuleDf,by.x='CC',by.y='CC',all.x=T)->data
  merge(data,statusOutRuleDf,by.x='statusOut',by.y='statusOut',all.x=T)->data
  data[is.na(data)]<-'NtAvb'
  noRules<-nrow(ADRGRulesDf)
  ifelse(data$PDXIndADRGRule=='NtAvb'&data$PPROCIndADRGRule=='NtAvb','PDXErr&PROCErr',
         ifelse(data$PDXIndADRGRule=='NtAvb','PDXErr',
                ifelse(data$PPROCIndADRGRule=='NtAvb','PPROCErr','NoErr')))->data$Err
  unique(data[,c("ageIndADRGRule","bornWtIndADRGRule","LOSIndADRGRule","PDXIndADRGRule","PPROCIndADRGRule","CCIndADRGRule","statusOutIndADRGRule")])->dataMerge
  
  
  unlist(lapply(1:nrow(dataMerge),function(i){
    is.element(as.character(1:noRules),unlist(stri_split_fixed(dataMerge$ageIndADRGRule[i],"|")))&
      is.element(as.character(1:noRules),unlist(stri_split_fixed(dataMerge$bornWtIndADRGRule[i],"|")))&
      is.element(as.character(1:noRules),unlist(stri_split_fixed(dataMerge$LOSIndADRGRule[i],"|")))&
      is.element(as.character(1:noRules),unlist(stri_split_fixed(dataMerge$PDXIndADRGRule[i],"|")))&
      is.element(as.character(1:noRules),unlist(stri_split_fixed(dataMerge$PPROCIndADRGRule[i],"|")))&
      is.element(as.character(1:noRules),unlist(stri_split_fixed(dataMerge$CCIndADRGRule[i],"|")))&
      is.element(as.character(1:noRules),unlist(stri_split_fixed(dataMerge$statusOutIndADRGRule[i],"|")))->ind
    which(ind)->intSet
    if(length(intSet)==0) NA else {
      paste(intSet,collapse='|')}
  }))->dataMerge$joinRule
  
  
  as.data.table(dataMerge)->dataMerge
  merge(data,dataMerge,by=c("ageIndADRGRule","bornWtIndADRGRule","LOSIndADRGRule","PDXIndADRGRule","PPROCIndADRGRule","CCIndADRGRule","statusOutIndADRGRule"),all.x=T)->data
  
  ifelse(is.na(data$joinRule)&data$Err=='NoErr','joinErr',data$Err)->data$Err
  
  data[is.na(data)]<-'NA'
  IND<-NA
  
  paste("ADX",1:15,sep='')->ADXData
  paste('PROC',2:8,sep='')->APROCData
  
  as.data.frame(data)->data
  
  apply(data,1,function(x){
    if(x['joinRule']=='NA') NA else {
      as.numeric(unlist(strsplit(x['joinRule'],'|',fixed=T)))->ind
      lapply(ind,function(j){
        fn_multiVars(unlist(x[ADXData]),ADRGRulesDf[j,"ADX"])&fn_multiVars(unlist(x[APROCData]),ADRGRulesDf[j,"APROC"])
      })->res
      ind[which(unlist(res))]->indd
      if(length(indd)==1) {
        sum(weightMt[indd,],na.rm=T)->sumWeight
      } else {
        apply(weightMt[indd,],1,sum,na.rm=T)->sumWeight
      }
      indd[which.max(sumWeight)][1]
    }
    
  })->IND
  
  data$indRule<-IND
  if(all(is.na(IND))) {
    ADRGRulesDf[IND,c("MDCCode","DRGCode","DRGName")]->ADRGRes
    ADRGRes[1:nrow(data),]->ADRGRes
  } else {
    ADRGRulesDf[IND,c("MDCCode","DRGCode","DRGName")]->ADRGRes
  }
  data.frame(data,ADRGRes)->dataADRGRes
  return(dataADRGRes[,c(namesData,"Err",names(ADRGRes))])
}


#' CNDRGGrouperAll
#' 
#' CNDRGGrouper all takes fn_Trans fn_excl fn_cc and CNDRGGrouper together for grouping
#' 
#' @export
CNDRGGrouperAll<-function(data,adrgruleslst=CNDRGRulesLst20170619){
  require(dplyr)
  fn_Trans(data=data,id='id') %>% fn_ccexcl(data=.) %>% fn_cc(data=.) %>% CNDRGGrouper(data=.,ADRGRulesLst=adrgruleslst)->res
  return(res)
}
