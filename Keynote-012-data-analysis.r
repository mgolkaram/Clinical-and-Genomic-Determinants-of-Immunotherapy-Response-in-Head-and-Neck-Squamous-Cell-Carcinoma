df<-read.csv('PATH2FILE/Summary.csv')

clin.fac<-df
clin.fac$TMB<-ifelse(clin.fac$TMB>89.99,'high','low')
clin.fac$TP53<-ifelse(clin.fac$TP53>0,'mut','wt')
clin.fac$PDL1.del<-ifelse(clin.fac$LRR_threshold< 0,1,0)
clin.fac$purity<-ifelse(clin.fac$purity>0.54,'high','low')

clin.fac$subtype[clin.fac$TP53=="mut" & clin.fac$TMB == "low"]<-'subtype1'
clin.fac$subtype[clin.fac$TP53=="mut" & clin.fac$TMB == "high" & (clin.fac$PDL1.del==1)]<-'subtype2'
clin.fac$subtype[clin.fac$TP53=="mut" & clin.fac$TMB == "high" & !(clin.fac$PDL1.del==1)]<-'subtype3'
clin.fac$subtype[clin.fac$TP53=="wt" & clin.fac$TMB == "high"]<-'subtype4'
clin.fac$subtype[clin.fac$TP53=="wt" & clin.fac$TMB == "low" & clin.fac$purity == "low"]<-'subtype5'
clin.fac$subtype[clin.fac$TP53=="wt" & clin.fac$TMB == "low" & clin.fac$purity == "high"]<-'subtype6'
clin.fac<-clin.fac[!clin.fac$tumorName=='P177T',] ### POLE mutant
table(clin.fac$subtype[clin.fac$Response==1])
table(clin.fac$subtype[clin.fac$Response==0])


par(mfrow=c(1,2))
R<-c()
for(i in 1:6){
  R<-c(R,(mean(clin.fac$Response[clin.fac$subtype==paste0('subtype',i)],na.rm = T)))
}
barplot(R,ylab='response rate',names.arg = paste0('subtype',1:6,' \n(N=',table(clin.fac$subtype),')'),
        col='orange',main='Keynote-012')

clin.fac1<-clin.fac[!clin.fac$tumorName=='P177T',]
clin.fac1$Response[clin.fac1$tumorName=='P210T']<-0
R<-c()
for(i in 1:6){
  R<-c(R,(mean(clin.fac1$Response[clin.fac1$subtype==paste0('subtype',i)],na.rm = T)))
}
barplot(R,ylab='response rate',names.arg = paste0('subtype',1:6,' \n(N=',table(clin.fac$subtype),')'),
        col='orange',main='Keynote-012 \nafter removing the POLE mutant sample (P177T)')
clin.fac$risk<-NA
clin.fac$risk[clin.fac$subtype=='subtype1' | clin.fac$subtype=='subtype2' | clin.fac$subtype=='subtype6']<-'high risk'
clin.fac$risk[clin.fac$subtype=='subtype3' | clin.fac$subtype=='subtype4' | clin.fac$subtype=='subtype5']<-'low risk'

par(mfrow=c(1,2))
barplot(c(mean(clin.fac$Response[clin.fac$risk=='low risk'],na.rm = T),mean(clin.fac$Response[clin.fac$risk=='high risk'],na.rm = T)),ylab='response rate',
        col='orange',main='Keynote-012',names.arg = c('low risk','high risk'),ylim=c(0,0.4))
barplot(c(mean(clin.fac$Response[clin.fac$risk=='low risk'],na.rm = T),mean(clin.fac$Response[clin.fac$risk=='high risk' & clin.fac$tumorName!='P177T'],na.rm = T)),ylab='response rate',
        col='orange',main='Keynote-012\nafter removing the POLE mutant sample (P177T)',ylim=c(0,0.4),names.arg = c('low risk','high risk'))



