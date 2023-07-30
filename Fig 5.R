##### Figure 5 ######
#### Update the path to the location of resource files: ~/Desktop/publications/H&N
#### 70/30 CV ####
library(randomForest)
library(ROCR)
require(sva)
library(ggfortify)
library(M3C)
library(pROC)
require(randomForestSRC)
library("survcomp")
library(glmnet)
require(gbm)
require(matrixStats)

WES<-read.csv('~/Desktop/publications/H&N/WES-HNSC-features.csv')
# features<-c(
#     "TP53","PIK3CA","TERT.promoter","PDL1","del9q34.3","MET.amp","HLALOH","Smoking_Signature","APOBEC_Signature",
#     "nonsyn.TMB","Indel_load",
#     "HED","purity","ploidy",
#     "GENDER","Age_at_IO_start","ECOG","ALCOHOL_NEVER_EVER",
#     "AnyVirus","ALLERGIES","PRIMARY_TUMOR_SITE","STEROIDS","ANTIBIOTICS","TOXICITY",
#     "METASTATIC_VS_RECURRENT","Number_sites_mets","Neut.Mono.Lymp","Platelets","AUTOIMMUNE_DISEASES",
#     "HGB","Albumin","PFS","Progression"
# )
# features<-c(
#   "TP53","PIK3CA","PDL1","del9q34.3","MET.amp","HLALOH","Smoking_Signature","APOBEC_Signature",
#   "nonsyn.TMB",
#   "purity","Age_at_IO_start","ECOG",
#   "AnyVirus","PRIMARY_TUMOR_SITE","TOXICITY",
#   "METASTATIC_VS_RECURRENT","Number_sites_mets","Neut.Mono.Lymp","Platelets",
#   "HGB","Albumin","PFS","Progression"
# )

features<-c("TP53","PIK3CA","PDL1","del9q34.3","MET.amp","Smoking_Signature","APOBEC_Signature",
            "nonsyn.TMB","Indel_load","Age_at_IO_start","GENDER","ECOG","AUTOIMMUNE_DISEASES",
            "METASTATIC_VS_RECURRENT","Neut.Mono.Lymp","Platelets","Albumin","INFECTION_DURING_IO",
            "ANTIBIOTICS","PRIMARY_TUMOR_SITE","HLALOH","del9p","AnyVirus",
            "OS","Overall_Survival_Event")
features<-c("TP53","PIK3CA","PDL1","del9q34.3","MET.amp","Smoking_Signature","APOBEC_Signature",
            "nonsyn.TMB","Indel_load","Age_at_IO_start","GENDER","ECOG","AUTOIMMUNE_DISEASES",
            "METASTATIC_VS_RECURRENT","Neut.Mono.Lymp","Platelets","Albumin","INFECTION_DURING_IO",
            "ANTIBIOTICS","PRIMARY_TUMOR_SITE","HLALOH","del9p","AnyVirus",
            "PFS","Progression")

WES$Number_sites_mets[WES$Number_sites_mets=="N/A"]<-0
WES$Number_sites_mets<-as.numeric(WES$Number_sites_mets)
WES$METASTATIC_VS_RECURRENT<-ifelse(WES$METASTATIC_VS_RECURRENT=="Metastatic","Metastatic","Non-metastatic")
WES$HED<-(WES$tumor.HLA.A+WES$tumor.HLA.B+WES$tumor.HLA.C)/3

train.idx<-read.csv('~/Desktop/publications/H&N/training.csv')$x
test.idx<-read.csv('~/Desktop/publications/H&N/test.csv')$x
#test.idx<-train.idx
clin.fac<-WES[WES$Patient.ID%in%c(train.idx,test.idx),]
#clin.fac$PFS<-clin.fac$OS; clin.fac$Progression<-clin.fac$Overall_Survival_Event
df<- clin.fac[,colnames(clin.fac)%in%features]
rownames(df)<-clin.fac$sample
df<-df[complete.cases(df), ]
df$PDL1[df$PDL1==1]<-'amp'
df$PDL1[df$PDL1==-1]<-'del'
df$PDL1[df$PDL1==0]<-'wt'
df$PRIMARY_TUMOR_SITE<-as.factor(ifelse(df$PRIMARY_TUMOR_SITE=="Oral Cavity", "Oral Cavity", "Non-oral cavity"))
df$APOBEC_Signature<-as.factor(ifelse(df$APOBEC_Signature=='yes','Present','Absent'))
df$Smoking_Signature<-as.factor(ifelse(df$Smoking_Signature=='yes','Present','Absent'))
df$AUTOIMMUNE_DISEASES<-as.factor(ifelse(df$AUTOIMMUNE_DISEASES=='Yes','Present','Absent'))
#df$ALLERGIES<-as.factor(ifelse(df$ALLERGIES=='Yes','Present','Absent'))
df$MET.amp<-as.factor(df$MET.amp)

df$ECOG<-(ifelse(df$ECOG==0,"0",'1-2')) 

for(i in 1:ncol(df)){
  if(colnames(df)[i]%in%c('METASTATIC_VS_RECURRENT','AnyVirus','ALCOHOL_NEVER_EVER','GENDER','TP53','PIK3CA',
                          'Smoking_Signature','APOBEC_Signature','PDL1','TERT.promoter','HLALOH','del9q34.3',
                          "ANTIBIOTICS",'STEROIDS','TOXICITY','ECOG','INFECTION_DURING_IO')){
    df[,i]<-as.factor(df[,i])
  } else if(colnames(df)[i]%in%c('HED','OS','PFS','Albumin','HGB','Platelets','Neut.Mono.Lymp','Age_at_IO_start','ploidy',
                                 'purity','Indel_load','nonsyn.TMB')) {
    df[,i]<-as.numeric(df[,i])
  }
}

# for(i in 1:1000){
#   P=0; IDX<-c()
#   idx<-sort(sample(seq(nrow(df)),size = round(nrow(df)*0.6)))
#   penalty<-0
#   for(j in 1:ncol(df)){
#     if(class(df[,j])=='character'){
#       penalty<-penalty+fisher.test(rbind(as.numeric(table(df[idx,j])),as.numeric(table(df[-idx,j]))))$p.val
#     } else {
#       penalty<-penalty+wilcox.test(as.numeric(df[idx,j]),as.numeric(df[-idx,j]))$p.val
#     }
#   }
#   if(penalty>P){P<-penalty; IDX<-idx}
#   
# }
# for(j in 1:ncol(df)){
#   if(class(df[,j])=='character'){
#     print(rbind(table(df[idx,j])/sum(table(df[idx,j])),table(df[-idx,j])/sum(table(df[-idx,j]))))
#   } else {
#     print(wilcox.test(as.numeric(df[idx,j]),as.numeric(df[-idx,j]))$p.val)
#   }
# }
#   


y<-c()

x<-c()



if(sum(colnames(df)%in%'PFS')>0){df$PFS<-df$PFS*30}
if(sum(colnames(df)%in%'OS')>0){df$OS<-df$OS*30}

train.idx<-sort(which(rownames(df)%in%(WES$sample[WES$Patient.ID%in%train.idx])))
test.idx<-sort(which(rownames(df)%in%(WES$sample[WES$Patient.ID%in%test.idx])))



K<-which(colnames(df)%in%c('Overall_Survival_Event','OS'))
K<-which(colnames(df)%in%c('Progression','PFS'))

training.df<-as.data.frame(df[train.idx,])
validation1.df<-as.data.frame(df[test.idx,-K])

for(count in 1:100){
  rf_classifier  <- rfsrc(Surv(PFS,Progression) ~ ., data = training.df, importance = "permute",na.action = "na.impute",
                          ntree = 100, block.size = 1)
  # rf_classifier  <- rfsrc(Surv(OS,Overall_Survival_Event) ~ ., data = training.df, importance = "permute",na.action = "na.impute",
  #                         ntree = 100, block.size = 1)

  rf_classifier
  y<-cbind(y,(rf_classifier$importance)[order(names(rf_classifier$importance))])
}
data<-y
rownames(data)<-c('Age','Albumin','Antibiotics','Viral status','APOBEC signature','Autoimmune disease','9p loss',
                  '9q34.3 loss','ECOG','Sex','HLA LOH','Indel load','Infection','MET gain','Metastatic vs Non-metastatic',
                  'SIRI','nonsyn TMB','PDL1 SCNA','PIK3CA status','Platelets','Tumor site','Smoking signature','TP53 status')
data<-t(data)
dev.off()
par(mar=c(7.1, 12.1, 4.1, 2.1))
boxplot(data[,order(colSums(data))],las =2,col='lightblue',
        horizontal=TRUE,cex.axis=0.8,outline=F,xlab='Feature contribution (permutation)')
#prediction_for_roc_curve <- predict(rf_classifier,validation1.df[,-K])

prediction_for_roc_curve <- predict(rf_classifier,as.data.frame(df[test.idx,-K]))

V<-prediction_for_roc_curve
x<-apply(X = V$survival, MARGIN = 1, FUN = function(W) V$time.interest[which.min(abs(W-0.5))]) 
#x<-apply(X = V$survival, MARGIN = 1, FUN = function(W) sum(W)) 
#x<-prediction_for_roc_curve
###
clin.fac<-WES[WES$sample%in%rownames(df),]
clin.fac$OS<-clin.fac$OS*30; clin.fac$PFS<-clin.fac$PFS*30; clin.fac$DSS<-clin.fac$DSS*30
clin.fac$RF23[clin.fac$sample%in%rownames(df[test.idx,])]<-x


#### RF11 ####
features<-c((rownames(y)[order(rowMedians(y),decreasing = T)])[1:11],"Overall_Survival_Event","OS")
features<-c((rownames(y)[order(rowMedians(y),decreasing = T)])[1:14],"Progression","PFS")

df<- df[,colnames(df)%in%features]
x<-y<-c()


K<-which(colnames(df)%in%c('Overall_Survival_Event','OS'))
K<-which(colnames(df)%in%c('Progression','PFS'))

training.df<-as.data.frame(df[train.idx,])
validation1.df<-as.data.frame(df[test.idx,-K])

rf_classifier  <- rfsrc(Surv(OS,Overall_Survival_Event) ~ ., data = training.df, importance = "permute",na.action = "na.impute",
                        ntree = 100, block.size = 1)
rf_classifier  <- rfsrc(Surv(PFS,Progression) ~ ., data = training.df, importance = "permute",na.action = "na.impute",
                        ntree = 100, block.size = 1)
rf_classifier
y<-cbind(y,(rf_classifier$importance)[order(names(rf_classifier$importance))])

prediction_for_roc_curve <- predict(rf_classifier,validation1.df)
V<-prediction_for_roc_curve
x<-apply(X = V$survival, MARGIN = 1, FUN = function(W) V$time.interest[which.min(abs(W-0.5))]) 
#}
clin.fac$RF11[clin.fac$sample%in%rownames(df[test.idx,])]<-x

###

clin.fac1<-clin.fac[!is.na(clin.fac$RF11),]
my_palette <- colorRampPalette(c("yellow","purple"))(n = 4)
dev.off()

clin.fac1$BEST_RESPONSE<-ifelse(clin.fac1$BEST_RESPONSE=='CR' | clin.fac1$BEST_RESPONSE=='PR', 1,0 )
RRF23<-as.numeric(roc( clin.fac1$BEST_RESPONSE,clin.fac1$RF23,smoothed = TRUE)$auc)
RRF11<-as.numeric(roc( clin.fac1$BEST_RESPONSE,clin.fac1$RF11,smoothed = TRUE)$auc)
RTMB<-as.numeric(roc( clin.fac1$BEST_RESPONSE,clin.fac1$nonsyn.TMB,smoothed = TRUE)$auc)

par(mfrow=c(1,2))

C.idx.24<-concordance.index(x=-clin.fac1$RF23, surv.time=clin.fac1$PFS, surv.event=clin.fac1$Progression, method="noether")$c.index
C.idx.12<-concordance.index(x=-clin.fac1$RF11, surv.time=clin.fac1$PFS, surv.event=clin.fac1$Progression, method="noether")$c.index
C.idx.TMB<-concordance.index(x=-clin.fac1$nonsyn.TMB, surv.time=clin.fac1$PFS, surv.event=clin.fac1$Progression, method="noether")$c.index
barplot(c(C.idx.24,C.idx.12,C.idx.TMB), ylim=c(0,1), col=c('red','blue','green'), names=c('RF23','RF11','nonsyn TMB'), ylab='C-index',main = 'PFS', las=2, cex.names = 0.8)

C.idx.24<-concordance.index(x=-clin.fac1$RF23, surv.time=clin.fac1$OS, surv.event=clin.fac1$Overall_Survival_Event, method="noether")$c.index
C.idx.12<-concordance.index(x=-clin.fac1$RF11, surv.time=clin.fac1$OS, surv.event=clin.fac1$Overall_Survival_Event, method="noether")$c.index
C.idx.TMB<-concordance.index(x=-clin.fac1$nonsyn.TMB, surv.time=clin.fac1$OS, surv.event=clin.fac1$Overall_Survival_Event, method="noether")$c.index
barplot(c(C.idx.24,C.idx.12,C.idx.TMB), ylim=c(0,1), col=c('red','blue','green'), names=c('RF23','RF11','nonsyn TMB'), ylab='C-index',main = 'OS', las=2, cex.names = 0.8)

#### ROC ####
## Note the MASS package masks select()!
library(tidyverse)
## https://github.com/tidyverse/tibble/issues/395
options(crayon.enabled = FALSE)
## Used for the dataset.
library(survival)
## Used for visualizaiton.
library(survminer)
## Load the Ovarian Cancer Survival Data

FP.idx<-seq(0.01,0.99,length.out = 1000)

library(timeROC)
library(survivalROC)

## Define a helper functio nto evaluate at various t
dev.off()
#### PFS ROC  ####
par(mfrow=c(1,3))
t=30*6; my.title = paste('PFS at',t/30,'months')

coxph1 <- coxph(Surv(PFS, Progression) ~ clin.fac1$RF23, data = clin.fac1)
HR24<-summary(coxph1)$coefficients[2]
P24<-summary(coxph1)$coefficients[5]

clin.fac1$lp <- predict(coxph1, type = "lp")
V<-timeROC(T=clin.fac1$PFS,delta=clin.fac1$Progression,
           marker=clin.fac1$lp,cause=1,
           weighting="marginal",
           times=c(t),ROC=TRUE)
RF23.TP<-predict(loess(V$TP[,2]~V$FP[,2]), newdata = FP.idx)
AUC1<-round(V$AUC[2],2)
RF23.AUC<-AUC1
plot(V,time=t,col="red",lwd=2,add=F, title=FALSE)
title(my.title)

# V<-survivalROC(Stime        = clin.fac1$PFS,
#             status       = clin.fac1$Progression,
#             marker       = clin.fac1$lp,
#             predict.time = t,
#             method       = "NNE",
#             span = 0.25 * nrow(clin.fac1)^(-0.20))
# AUC1<-round(V$AUC,2)
# plot(V$FP,V$TP,col="red",lwd=2,type='l')

coxph1 <- coxph(Surv(PFS, Progression) ~ clin.fac1$RF11, data = clin.fac1)
HR12<-summary(coxph1)$coefficients[2]
P12<-summary(coxph1)$coefficients[5]
print(summary(coxph1)$coefficients[5])
clin.fac1$lp <- predict(coxph1, type = "lp")
V<-timeROC(T=clin.fac1$PFS,delta=clin.fac1$Progression,
           marker=clin.fac1$lp,cause=1,
           weighting="marginal",
           times=c(t),ROC=TRUE)
RF11.TP<-predict(loess(V$TP[,2]~V$FP[,2]), newdata = FP.idx)

AUC2<-round(V$AUC[2],2)
RF11.AUC<-AUC2

plot(V,time=t,col="blue",add=TRUE,title=FALSE,lwd=2)


coxph1 <- coxph(Surv(PFS, Progression) ~ clin.fac1$nonsyn.TMB, data = clin.fac1)
HRTMB<-summary(coxph1)$coefficients[2]
PTMB<-summary(coxph1)$coefficients[5]
print(summary(coxph1)$coefficients[5])
clin.fac1$lp <- predict(coxph1, type = "lp")
V<-timeROC(T=clin.fac1$PFS,delta=clin.fac1$Progression,
           marker=clin.fac1$lp,cause=1,
           weighting="marginal",
           times=c(t),ROC=TRUE)
TMB.TP<-predict(loess(V$TP[,2]~V$FP[,2]), newdata = FP.idx)

AUC3<-round(V$AUC[2],2)
TMB.AUC<-AUC3
plot(V,time=t,col="green",add=TRUE,title=FALSE,lwd=2)


legend(0.5,0.5,col=c('red','blue','green'),lwd=2,y.intersp=1.5,c(paste0('RF23\n (AUC=',AUC1,')'),
                                                                 paste0('RF11\n (AUC=',AUC2,')'),
                                                                 paste0('nonsyn TMB\n (AUC=',AUC3,')')),bty = 'n',cex=1.2)

#### OS ROC  ####
t=30*12; my.title = paste('OS at',t/30,'months')

coxph1 <- coxph(Surv(OS, Overall_Survival_Event) ~ clin.fac1$RF23, data = clin.fac1)
HR24<-summary(coxph1)$coefficients[2]
P24<-summary(coxph1)$coefficients[5]

clin.fac1$lp <- predict(coxph1, type = "lp")
V<-timeROC(T=clin.fac1$OS,delta=clin.fac1$Overall_Survival_Event,
           marker=clin.fac1$lp,cause=1,
           weighting="marginal",
           times=c(t),ROC=TRUE)
RF23.TP<-predict(loess(V$TP[,2]~V$FP[,2]), newdata = FP.idx)
AUC1<-round(V$AUC[2],2)
RF23.AUC<-AUC1
plot(V,time=t,col="red",lwd=2,add=F, title=FALSE)
title(my.title)

# V<-survivalROC(Stime        = clin.fac1$PFS,
#             status       = clin.fac1$Progression,
#             marker       = clin.fac1$lp,
#             predict.time = t,
#             method       = "NNE",
#             span = 0.25 * nrow(clin.fac1)^(-0.20))
# AUC1<-round(V$AUC,2)
# plot(V$FP,V$TP,col="red",lwd=2,type='l')

coxph1 <- coxph(Surv(OS, Overall_Survival_Event) ~ clin.fac1$RF11, data = clin.fac1)
HR12<-summary(coxph1)$coefficients[2]
P12<-summary(coxph1)$coefficients[5]
print(summary(coxph1)$coefficients[5])
clin.fac1$lp <- predict(coxph1, type = "lp")
V<-timeROC(T=clin.fac1$OS,delta=clin.fac1$Overall_Survival_Event,
           marker=clin.fac1$lp,cause=1,
           weighting="marginal",
           times=c(t),ROC=TRUE)
RF11.TP<-predict(loess(V$TP[,2]~V$FP[,2]), newdata = FP.idx)

AUC2<-round(V$AUC[2],2)
RF11.AUC<-AUC2

plot(V,time=t,col="blue",add=TRUE,title=FALSE,lwd=2)


coxph1 <- coxph(Surv(OS, Overall_Survival_Event) ~ clin.fac1$nonsyn.TMB, data = clin.fac1)
HRTMB<-summary(coxph1)$coefficients[2]
PTMB<-summary(coxph1)$coefficients[5]
print(summary(coxph1)$coefficients[5])
clin.fac1$lp <- predict(coxph1, type = "lp")
V<-timeROC(T=clin.fac1$OS,delta=clin.fac1$Overall_Survival_Event,
           marker=clin.fac1$lp,cause=1,
           weighting="marginal",
           times=c(t),ROC=TRUE)
TMB.TP<-predict(loess(V$TP[,2]~V$FP[,2]), newdata = FP.idx)

AUC3<-round(V$AUC[2],2)
TMB.AUC<-AUC3
plot(V,time=t,col="green",add=TRUE,title=FALSE,lwd=2)


legend(0.5,0.5,col=c('red','blue','green'),lwd=2,y.intersp=1.5,c(paste0('RF23\n (AUC=',AUC1,')'),
                                                                 paste0('RF11\n (AUC=',AUC2,')'),
                                                                 paste0('nonsyn TMB\n (AUC=',AUC3,')')),bty = 'n',cex=1.2)
#### response ROC  ####
res.RF23<-roc(clin.fac1$BEST_RESPONSE,clin.fac1$RF23,smoothed = TRUE)
res.RF11<-roc(clin.fac1$BEST_RESPONSE,clin.fac1$RF11,smoothed = TRUE)
res.TMB<-roc(clin.fac1$BEST_RESPONSE,clin.fac1$nonsyn.TMB,smoothed = TRUE)
plot(1-res.RF23$specificities,res.RF23$sensitivities,col='red',type='l',lwd=2,xlab='1-Specificy',ylab='Sensitivity',main='Response')
lines(1-res.RF11$specificities,res.RF11$sensitivities,col='blue',lwd=2)
lines(1-res.TMB$specificities,res.TMB$sensitivities,col='green',lwd=2)
lines(c(-1,2),c(-1,2),lty=2)
AUC1=round(res.RF23$auc,2)
AUC2=round(res.RF11$auc,2)
AUC3=round(res.TMB$auc,2)
legend(0.5,0.5,col=c('red','blue','green'),lwd=2,y.intersp=1.5,c(paste0('RF23\n (AUC=',AUC1,')'),
                                                                 paste0('RF11\n (AUC=',AUC2,')'),
                                                                 paste0('nonsyn TMB\n (AUC=',AUC3,')')),bty = 'n',cex=1.2)
####
clin.fac1$RF11=as.factor(ifelse(clin.fac1$RF11>median(clin.fac1$RF11),'Predicted R','Predicted NR'))
clin.fac1$RF23=as.factor(ifelse(clin.fac1$RF23>median(clin.fac1$RF23),'Predicted R','Predicted NR'))

#### RF23 ####
HR.HPV<-round(summary(coxph(Surv(PFS, Progression) ~ RF23, data = clin.fac1))$coefficients[2],2)
fit.HPV <- survfit(Surv(PFS, Progression) ~ RF23, data = clin.fac1)
ggsurv.HPV<-ggsurvplot(fit.HPV, data = clin.fac1, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (years)', xlim = c(0, 365*4), xscale = "d_y",  break.x.by = 365,  palette = c("#E7B800", "#2E9FDF"))
ggsurv.HPV$plot <- ggsurv.HPV$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV), size = 5)
ggsurv.HPV
HR.HPV<-round(summary(coxph(Surv(OS, Overall_Survival_Event) ~ RF23, data = clin.fac1))$coefficients[2],2)
fit.HPV <- survfit(Surv(OS, Overall_Survival_Event) ~ RF23, data = clin.fac1)
ggsurv.HPV<-ggsurvplot(fit.HPV, data = clin.fac1, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'OS', xlab = 'Time (years)', xlim = c(0, 365*4), xscale = "d_y",  break.x.by = 365,  palette = c("#E7B800", "#2E9FDF"))
ggsurv.HPV$plot <- ggsurv.HPV$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV), size = 5)
ggsurv.HPV
#### RF11 ####
HR.HPV<-round(summary(coxph(Surv(PFS, Progression) ~ RF11, data = clin.fac1))$coefficients[2],2)
fit.HPV <- survfit(Surv(PFS, Progression) ~ RF11, data = clin.fac1)
ggsurv.HPV<-ggsurvplot(fit.HPV, data = clin.fac1, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (years)', xlim = c(0, 365*4), xscale = "d_y",  break.x.by = 365,  palette = c("#E7B800", "#2E9FDF"))
ggsurv.HPV$plot <- ggsurv.HPV$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV), size = 5)
ggsurv.HPV
HR.HPV<-round(summary(coxph(Surv(OS, Overall_Survival_Event) ~ RF11, data = clin.fac1))$coefficients[2],2)
fit.HPV <- survfit(Surv(OS, Overall_Survival_Event) ~ RF11, data = clin.fac1)
ggsurv.HPV<-ggsurvplot(fit.HPV, data = clin.fac1, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'OS', xlab = 'Time (years)', xlim = c(0, 365*4), xscale = "d_y",  break.x.by = 365,  palette = c("#E7B800", "#2E9FDF"))
ggsurv.HPV$plot <- ggsurv.HPV$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV), size = 5)
ggsurv.HPV

#### TMB ####
clin.fac1$nonsyn.TMB<-as.factor(ifelse(clin.fac1$nonsyn.TMB<3.34,"High","Low"))
levels(clin.fac1$nonsyn.TMB)<-c('Low','High')
HR.HPV<- round(summary(coxph(Surv(PFS, Progression) ~ nonsyn.TMB, data = clin.fac1))$coefficients[2],2)
fit.HPV <- survfit(Surv(PFS, Progression) ~ nonsyn.TMB, data = clin.fac1)
ggsurv.HPV<-ggsurvplot(fit.HPV, data = clin.fac1, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (years)', xlim = c(0, 365*4), xscale = "d_y",  break.x.by = 365,  palette = c("#E7B800","#2E9FDF"))
ggsurv.HPV$plot <- ggsurv.HPV$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV), size = 5)
ggsurv.HPV
HR.HPV<- round(summary(coxph(Surv(OS, Overall_Survival_Event) ~ nonsyn.TMB, data = clin.fac1))$coefficients[2],2)
fit.HPV <- survfit(Surv(OS, Overall_Survival_Event) ~ nonsyn.TMB, data = clin.fac1)
ggsurv.HPV<-ggsurvplot(fit.HPV, data = clin.fac1, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'OS', xlab = 'Time (years)', xlim = c(0, 365*4), xscale = "d_y",  break.x.by = 365,  palette = c("#E7B800", "#2E9FDF"))
ggsurv.HPV$plot <- ggsurv.HPV$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV), size = 5)
ggsurv.HPV


#### IMPACT ####
require(randomForestSRC)
require(survcomp)
require(pROC)
WES<-read.csv('~/Desktop/publications/H&N/WES-HNSC-features.csv')
features<-c("TP53","PIK3CA","PDL1","del9q34.3","Smoking_Signature",
            "nonsyn.TMB","Age_at_IO_start","GENDER","ECOG","AUTOIMMUNE_DISEASES",
            "METASTATIC_VS_RECURRENT","Neut.Mono.Lymp","Platelets","Albumin","INFECTION_DURING_IO",
            "ANTIBIOTICS","PRIMARY_TUMOR_SITE","HLALOH","del9p","AnyVirus",
            "OS","Overall_Survival_Event")
features<-c("TP53","PIK3CA","PDL1","del9q34.3","Smoking_Signature",
            "nonsyn.TMB","Age_at_IO_start","GENDER","ECOG","AUTOIMMUNE_DISEASES",
            "METASTATIC_VS_RECURRENT","Neut.Mono.Lymp","Platelets","Albumin","INFECTION_DURING_IO",
            "ANTIBIOTICS","PRIMARY_TUMOR_SITE","HLALOH","del9p","AnyVirus",
            "PFS","Progression")

features<-c("TP53","PIK3CA","PDL1","del9q34.3","Smoking_Signature",
            "nonsyn.TMB","Age_at_IO_start","ECOG","AUTOIMMUNE_DISEASES",
            "Neut.Mono.Lymp","Platelets",
            "del9p","AnyVirus","METASTATIC_VS_RECURRENT",
            "PFS","Progression")
features<-c("TP53","PIK3CA","PDL1","del9q34.3","Smoking_Signature",
            "nonsyn.TMB","Age_at_IO_start","ECOG","AUTOIMMUNE_DISEASES",
            "Neut.Mono.Lymp","Platelets",'Indel_load','APOBEC_Signature',
            "del9p",
            "PFS","Progression")
WES$Number_sites_mets[WES$Number_sites_mets=="N/A"]<-0
WES$Number_sites_mets<-as.numeric(WES$Number_sites_mets)
WES$METASTATIC_VS_RECURRENT<-ifelse(WES$METASTATIC_VS_RECURRENT=="Metastatic","Metastatic","Non-metastatic")
WES$HED<-(WES$tumor.HLA.A+WES$tumor.HLA.B+WES$tumor.HLA.C)/3

train.idx<-read.csv('~/Desktop/publications/H&N/training.csv')$x
test.idx<-c()
clin.fac<-WES[WES$Patient.ID%in%c(train.idx,test.idx),]
df<- clin.fac[,colnames(clin.fac)%in%features]
rownames(df)<-clin.fac$sample
df<-df[complete.cases(df), ]
df$PDL1[df$PDL1==1]<-'amp'
df$PDL1[df$PDL1==-1]<-'del'
df$PDL1[df$PDL1==0]<-'wt'
df$PRIMARY_TUMOR_SITE<-as.factor(ifelse(df$PRIMARY_TUMOR_SITE=="Oral Cavity", "Oral Cavity", "Non-oral cavity"))
df$Smoking_Signature<-as.factor(ifelse(df$Smoking_Signature=='yes','Present','Absent'))
df$AUTOIMMUNE_DISEASES<-as.factor(ifelse(df$AUTOIMMUNE_DISEASES=='Yes','Present','Absent'))
df$ECOG<-(ifelse(df$ECOG==0,"0",'1-2')) 

for(i in 1:ncol(df)){
  if(colnames(df)[i]%in%c('METASTATIC_VS_RECURRENT','AnyVirus','ALCOHOL_NEVER_EVER','GENDER','TP53','PIK3CA',
                          'Smoking_Signature','APOBEC_Signature','PDL1','TERT.promoter','HLALOH','del9q34.3',
                          "ANTIBIOTICS",'STEROIDS','TOXICITY','ECOG','INFECTION_DURING_IO')){
    df[,i]<-as.factor(df[,i])
  } else if(colnames(df)[i]%in%c('HED','PFS','OS','Albumin','HGB','Platelets','Neut.Mono.Lymp','Age_at_IO_start','ploidy',
                                 'purity','Indel_load','nonsyn.TMB')) {
    df[,i]<-as.numeric(df[,i])
  }
}

y<-c()
x<-c()
if(sum(colnames(df)%in%'PFS')>0){df$PFS<-df$PFS*30}
if(sum(colnames(df)%in%'OS')>0){df$OS<-df$OS*30}

train.idx<-sort(which(rownames(df)%in%(WES$sample[WES$Patient.ID%in%train.idx])))

#K<-which(colnames(df)%in%c('Overall_Survival_Event','OS'))
K<-which(colnames(df)%in%c('Progression','PFS'))
training.df<-as.data.frame(df[train.idx,])

for(count in 1:100){
  rf_classifier  <- rfsrc(Surv(PFS,Progression) ~ ., data = training.df, importance = "permute",na.action = "na.impute",
                          ntree = 100, block.size = 1)
  # rf_classifier  <- rfsrc(Surv(OS,Overall_Survival_Event) ~ ., data = training.df, importance = "permute",na.action = "na.impute",
  #                         ntree = 100, block.size = 1)
  # 
}
test.set<-read.csv('~/Desktop/publications/H&N/IMPACT-validation.csv',stringsAsFactors = F)

test.set$Indel_load<-NA
test.set$APOBEC_Signature<-NA
test.set<-test.set[,colnames(test.set)%in%colnames(df)]
tmp<-test.set
for(k in 1:ncol(tmp)){tmp[,k]<-test.set[,colnames(test.set)%in%colnames(df)[k]]}
test.set<-tmp
colnames(test.set)<-colnames(df)
tmp<-rbind(df,test.set)
for(i in 1:ncol(tmp)){
  if(colnames(tmp)[i]%in%c('METASTATIC_VS_RECURRENT','AnyVirus','ALCOHOL_NEVER_EVER','GENDER','TP53','PIK3CA',
                          'Smoking_Signature','APOBEC_Signature','PDL1','TERT.promoter','HLALOH','del9q34.3',
                          "ANTIBIOTICS",'STEROIDS','TOXICITY','ECOG','INFECTION_DURING_IO')){
    tmp[,i]<-as.factor(tmp[,i])
  } else if(colnames(tmp)[i]%in%c('HED','OS','PFS','Albumin','HGB','Platelets','Neut.Mono.Lymp','Age_at_IO_start','ploidy',
                                 'purity','Indel_load','nonsyn.TMB')) {
    tmp[,i]<-as.numeric(tmp[,i])
  }
}
test.set<-tmp[-c(1:nrow(df)),]
prediction_for_roc_curve <- predict(rf_classifier,as.data.frame(test.set[,-K]),na.action="na.impute")

V<-prediction_for_roc_curve
x<-apply(X = V$survival, MARGIN = 1, FUN = function(W) V$time.interest[which.min(abs(W-0.5))]) 
#x<-apply(X = V$survival, MARGIN = 1, FUN = function(W) sum(W)) 
#x<-prediction_for_roc_curve
###
clin.fac1<-test.set
clin.fac1$RF14<-x
clin.fac1$BEST_RESPONSE<-as.factor(ifelse(read.csv('~/Desktop/publications/H&N/IMPACT-validation.csv')$BEST_RESPONSE=='CR' | read.csv('~/Desktop/publications/H&N/IMPACT-validation.csv')$BEST_RESPONSE=='PR', 1,0 ))
RRF14<-as.numeric(roc( clin.fac1$BEST_RESPONSE,clin.fac1$RF14,smoothed = TRUE)$auc)
RTMB<-as.numeric(roc( clin.fac1$BEST_RESPONSE,clin.fac1$nonsyn.TMB,smoothed = TRUE)$auc)

clin.fac1$Overall_Survival_Event<-read.csv('~/Desktop/publications/H&N/IMPACT-validation.csv')$Overall_Survival_Event
clin.fac1$OS<-read.csv('~/Desktop/publications/H&N/IMPACT-validation.csv')$OS*30
clin.fac1$Progression<-read.csv('~/Desktop/publications/H&N/IMPACT-validation.csv')$Progression
clin.fac1$PFS<-read.csv('~/Desktop/publications/H&N/IMPACT-validation.csv')$PFS*30

par(mfrow=c(1,2))

C.idx.24<-concordance.index(x=-clin.fac1$RF14, surv.time=clin.fac1$PFS, surv.event=clin.fac1$Progression, method="noether")$c.index
C.idx.TMB<-concordance.index(x=-clin.fac1$nonsyn.TMB, surv.time=clin.fac1$PFS, surv.event=clin.fac1$Progression, method="noether")$c.index
barplot(c(C.idx.24,C.idx.TMB), ylim=c(0,1), col=c('red','green'), names=c('RF14','nonsyn TMB'), ylab='C-index',main = 'PFS', las=2, cex.names = 0.8)

C.idx.24<-concordance.index(x=-clin.fac1$RF14, surv.time=clin.fac1$OS, surv.event=clin.fac1$Overall_Survival_Event, method="noether")$c.index
C.idx.TMB<-concordance.index(x=-clin.fac1$nonsyn.TMB, surv.time=clin.fac1$OS, surv.event=clin.fac1$Overall_Survival_Event, method="noether")$c.index
barplot(c(C.idx.24,C.idx.TMB), ylim=c(0,1), col=c('red','green'), names=c('RF14','nonsyn TMB'), ylab='C-index',main = 'OS', las=2, cex.names = 0.8)

#### ROC ####
## Note the MASS package masks select()!
library(tidyverse)
## https://github.com/tidyverse/tibble/issues/395
options(crayon.enabled = FALSE)
## Used for the dataset.
library(survival)
## Used for visualizaiton.
library(survminer)
## Load the Ovarian Cancer Survival Data

FP.idx<-seq(0.01,0.99,length.out = 1000)

library(timeROC)
library(survivalROC)

## Define a helper functio nto evaluate at various t
dev.off()
#### PFS ROC  ####
par(mfrow=c(1,3))
t=30*6; my.title = paste('PFS at',t/30,'months')

coxph1 <- coxph(Surv(PFS, Progression) ~ clin.fac1$RF14, data = clin.fac1)
HR24<-summary(coxph1)$coefficients[2]
P24<-summary(coxph1)$coefficients[5]

clin.fac1$lp <- predict(coxph1, type = "lp")
V<-timeROC(T=clin.fac1$PFS,delta=clin.fac1$Progression,
           marker=clin.fac1$lp,cause=1,
           weighting="marginal",
           times=c(t),ROC=TRUE)
RF14.TP<-predict(loess(V$TP[,2]~V$FP[,2]), newdata = FP.idx)
AUC1<-round(V$AUC[2],2)
RF14.AUC<-AUC1
plot(V,time=t,col="red",lwd=2,add=F, title=FALSE)
title(my.title)

# V<-survivalROC(Stime        = clin.fac1$PFS,
#             status       = clin.fac1$Progression,
#             marker       = clin.fac1$lp,
#             predict.time = t,
#             method       = "NNE",
#             span = 0.25 * nrow(clin.fac1)^(-0.20))
# AUC1<-round(V$AUC,2)
# plot(V$FP,V$TP,col="red",lwd=2,type='l')

coxph1 <- coxph(Surv(PFS, Progression) ~ clin.fac1$nonsyn.TMB, data = clin.fac1)
HRTMB<-summary(coxph1)$coefficients[2]
PTMB<-summary(coxph1)$coefficients[5]
print(summary(coxph1)$coefficients[5])
clin.fac1$lp <- predict(coxph1, type = "lp")
V<-timeROC(T=clin.fac1$PFS,delta=clin.fac1$Progression,
           marker=clin.fac1$lp,cause=1,
           weighting="marginal",
           times=c(t),ROC=TRUE)
TMB.TP<-predict(loess(V$TP[,2]~V$FP[,2]), newdata = FP.idx)

AUC3<-round(V$AUC[2],2)
TMB.AUC<-AUC3
plot(V,time=t,col="green",add=TRUE,title=FALSE,lwd=2)


legend(0.5,0.5,col=c('red','green'),lwd=2,y.intersp=1.5,c(paste0('RF14\n (AUC=',AUC1,')'),
                                                                 paste0('nonsyn TMB\n (AUC=',AUC3,')')),bty = 'n',cex=1.2)

#### OS ROC  ####
t=30*12; my.title = paste('OS at',t/30,'months')

coxph1 <- coxph(Surv(OS, Overall_Survival_Event) ~ clin.fac1$RF14, data = clin.fac1)
HR24<-summary(coxph1)$coefficients[2]
P24<-summary(coxph1)$coefficients[5]

clin.fac1$lp <- predict(coxph1, type = "lp")
V<-timeROC(T=clin.fac1$OS,delta=clin.fac1$Overall_Survival_Event,
           marker=clin.fac1$lp,cause=1,
           weighting="marginal",
           times=c(t),ROC=TRUE)
RF14.TP<-predict(loess(V$TP[,2]~V$FP[,2]), newdata = FP.idx)
AUC1<-round(V$AUC[2],2)
RF14.AUC<-AUC1
plot(V,time=t,col="red",lwd=2,add=F, title=FALSE)
title(my.title)

# V<-survivalROC(Stime        = clin.fac1$PFS,
#             status       = clin.fac1$Progression,
#             marker       = clin.fac1$lp,
#             predict.time = t,
#             method       = "NNE",
#             span = 0.25 * nrow(clin.fac1)^(-0.20))
# AUC1<-round(V$AUC,2)
# plot(V$FP,V$TP,col="red",lwd=2,type='l')

coxph1 <- coxph(Surv(OS, Overall_Survival_Event) ~ clin.fac1$nonsyn.TMB, data = clin.fac1)
HRTMB<-summary(coxph1)$coefficients[2]
PTMB<-summary(coxph1)$coefficients[5]
print(summary(coxph1)$coefficients[5])
clin.fac1$lp <- predict(coxph1, type = "lp")
V<-timeROC(T=clin.fac1$OS,delta=clin.fac1$Overall_Survival_Event,
           marker=clin.fac1$lp,cause=1,
           weighting="marginal",
           times=c(t),ROC=TRUE)
TMB.TP<-predict(loess(V$TP[,2]~V$FP[,2]), newdata = FP.idx)

AUC3<-round(V$AUC[2],2)
TMB.AUC<-AUC3
plot(V,time=t,col="green",add=TRUE,title=FALSE,lwd=2)


legend(0.5,0.5,col=c('red','green'),lwd=2,y.intersp=1.5,c(paste0('RF14\n (AUC=',AUC1,')'),
                                                                 paste0('nonsyn TMB\n (AUC=',AUC3,')')),bty = 'n',cex=1.2)
#### response ROC  ####
require(pROC)
res.RF14<-roc(clin.fac1$BEST_RESPONSE,clin.fac1$RF14,smoothed = TRUE)
res.TMB<-roc(clin.fac1$BEST_RESPONSE,clin.fac1$nonsyn.TMB,smoothed = TRUE)
plot(1-res.RF14$specificities,res.RF14$sensitivities,col='red',type='l',lwd=2,xlab='1-Specificy',ylab='Sensitivity',main='Response')
lines(1-res.TMB$specificities,res.TMB$sensitivities,col='green',lwd=2)
lines(c(-1,2),c(-1,2),lty=2)
AUC1=round(res.RF14$auc,2)
AUC3=round(res.TMB$auc,2)
legend(0.5,0.5,col=c('red','green'),lwd=2,y.intersp=1.5,c(paste0('RF14\n (AUC=',AUC1,')'),
                                                                 paste0('nonsyn TMB\n (AUC=',AUC3,')')),bty = 'n',cex=1.2)
####
clin.fac1$RF14=as.factor(ifelse(clin.fac1$RF14>median(clin.fac1$RF14),'Predicted R','Predicted NR'))

#### RF14 ####
HR.HPV<-round(summary(coxph(Surv(PFS, Progression) ~ RF14, data = clin.fac1))$coefficients[2],2)
fit.HPV <- survfit(Surv(PFS, Progression) ~ RF14, data = clin.fac1)
ggsurv.HPV<-ggsurvplot(fit.HPV, data = clin.fac1, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (years)', xlim = c(0, 365*4), xscale = "d_y",  break.x.by = 365,  palette = c("#E7B800", "#2E9FDF"))
ggsurv.HPV$plot <- ggsurv.HPV$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV), size = 5)
ggsurv.HPV
HR.HPV<-round(summary(coxph(Surv(OS, Overall_Survival_Event) ~ RF14, data = clin.fac1))$coefficients[2],2)
fit.HPV <- survfit(Surv(OS, Overall_Survival_Event) ~ RF14, data = clin.fac1)
ggsurv.HPV<-ggsurvplot(fit.HPV, data = clin.fac1, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'OS', xlab = 'Time (years)', xlim = c(0, 365*4), xscale = "d_y",  break.x.by = 365,  palette = c("#E7B800", "#2E9FDF"))
ggsurv.HPV$plot <- ggsurv.HPV$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV), size = 5)
ggsurv.HPV

#### TMB ####
clin.fac1$nonsyn.TMB<-as.factor(ifelse(clin.fac1$nonsyn.TMB<3.34,"High","Low"))
levels(clin.fac1$nonsyn.TMB)<-c('Low','High')
HR.HPV<- round(summary(coxph(Surv(PFS, Progression) ~ nonsyn.TMB, data = clin.fac1))$coefficients[2],2)
fit.HPV <- survfit(Surv(PFS, Progression) ~ nonsyn.TMB, data = clin.fac1)
ggsurv.HPV<-ggsurvplot(fit.HPV, data = clin.fac1, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (years)', xlim = c(0, 365*4), xscale = "d_y",  break.x.by = 365,  palette = c("#E7B800","#2E9FDF"))
ggsurv.HPV$plot <- ggsurv.HPV$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV), size = 5)
ggsurv.HPV
HR.HPV<- round(summary(coxph(Surv(OS, Overall_Survival_Event) ~ nonsyn.TMB, data = clin.fac1))$coefficients[2],2)
fit.HPV <- survfit(Surv(OS, Overall_Survival_Event) ~ nonsyn.TMB, data = clin.fac1)
ggsurv.HPV<-ggsurvplot(fit.HPV, data = clin.fac1, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'OS', xlab = 'Time (years)', xlim = c(0, 365*4), xscale = "d_y",  break.x.by = 365,  palette = c("#E7B800", "#2E9FDF"))
ggsurv.HPV$plot <- ggsurv.HPV$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV), size = 5)
ggsurv.HPV



##### IMPACT non-IO #####
require(randomForestSRC)
require(survcomp)
require(pROC)
WES<-read.csv('~/Desktop/publications/H&N/WES-HNSC-features.csv')
features<-c("TP53","PIK3CA","PDL1","del9q34.3","Smoking_Signature",
            "nonsyn.TMB","Age_at_IO_start","GENDER","ECOG","AUTOIMMUNE_DISEASES",
            "METASTATIC_VS_RECURRENT","Neut.Mono.Lymp","Platelets","Albumin","INFECTION_DURING_IO",
            "ANTIBIOTICS","PRIMARY_TUMOR_SITE","HLALOH","del9p","AnyVirus",
            "OS","Overall_Survival_Event")
features<-c("TP53","PIK3CA","PDL1","del9q34.3","Smoking_Signature",
            "nonsyn.TMB","Age_at_IO_start","GENDER","ECOG","AUTOIMMUNE_DISEASES",
            "METASTATIC_VS_RECURRENT","Neut.Mono.Lymp","Platelets","Albumin","INFECTION_DURING_IO",
            "ANTIBIOTICS","PRIMARY_TUMOR_SITE","HLALOH","del9p","AnyVirus",
            "PFS","Progression")

features<-c("TP53","PIK3CA","PDL1","del9q34.3","Smoking_Signature",
            "nonsyn.TMB","Age_at_IO_start","ECOG","AUTOIMMUNE_DISEASES",
            "Neut.Mono.Lymp","Platelets",
            "del9p","AnyVirus","METASTATIC_VS_RECURRENT",
            "PFS","Progression")
features<-c("TP53","PIK3CA","PDL1","del9q34.3","Smoking_Signature",
            "nonsyn.TMB","Age_at_IO_start","ECOG","AUTOIMMUNE_DISEASES",
            "Neut.Mono.Lymp","Platelets",'Indel_load','APOBEC_Signature',
            "del9p",
            "PFS","Progression")
WES$Number_sites_mets[WES$Number_sites_mets=="N/A"]<-0
WES$Number_sites_mets<-as.numeric(WES$Number_sites_mets)
WES$METASTATIC_VS_RECURRENT<-ifelse(WES$METASTATIC_VS_RECURRENT=="Metastatic","Metastatic","Non-metastatic")
WES$HED<-(WES$tumor.HLA.A+WES$tumor.HLA.B+WES$tumor.HLA.C)/3

train.idx<-read.csv('~/Desktop/publications/H&N/training.csv')$x
test.idx<-c()
clin.fac<-WES[WES$Patient.ID%in%c(train.idx,test.idx),]
df<- clin.fac[,colnames(clin.fac)%in%features]
rownames(df)<-clin.fac$sample
df<-df[complete.cases(df), ]
df$PDL1[df$PDL1==1]<-'amp'
df$PDL1[df$PDL1==-1]<-'del'
df$PDL1[df$PDL1==0]<-'wt'
#df$PRIMARY_TUMOR_SITE<-as.factor(ifelse(df$PRIMARY_TUMOR_SITE=="Oral Cavity", "Oral Cavity", "Non-oral cavity"))
df$Smoking_Signature<-as.factor(ifelse(df$Smoking_Signature=='yes','Present','Absent'))
df$AUTOIMMUNE_DISEASES<-as.factor(ifelse(df$AUTOIMMUNE_DISEASES=='Yes','Present','Absent'))
df$ECOG<-(ifelse(df$ECOG==0,"0",'1-2')) 

for(i in 1:ncol(df)){
  if(colnames(df)[i]%in%c('METASTATIC_VS_RECURRENT','AnyVirus','ALCOHOL_NEVER_EVER','GENDER','TP53','PIK3CA',
                          'Smoking_Signature','APOBEC_Signature','PDL1','TERT.promoter','HLALOH','del9q34.3',
                          "ANTIBIOTICS",'STEROIDS','TOXICITY','ECOG','INFECTION_DURING_IO')){
    df[,i]<-as.factor(df[,i])
  } else if(colnames(df)[i]%in%c('HED','PFS','OS','Albumin','HGB','Platelets','Neut.Mono.Lymp','Age_at_IO_start','ploidy',
                                 'purity','Indel_load','nonsyn.TMB')) {
    df[,i]<-as.numeric(df[,i])
  }
}

y<-c()
x<-c()
if(sum(colnames(df)%in%'PFS')>0){df$PFS<-df$PFS*30}
if(sum(colnames(df)%in%'OS')>0){df$OS<-df$OS*30}

train.idx<-sort(which(rownames(df)%in%(WES$sample[WES$Patient.ID%in%train.idx])))

#K<-which(colnames(df)%in%c('Overall_Survival_Event','OS'))
K<-which(colnames(df)%in%c('Progression','PFS'))
training.df<-as.data.frame(df[train.idx,])

for(count in 1:100){
  rf_classifier  <- rfsrc(Surv(PFS,Progression) ~ ., data = training.df, importance = "permute",na.action = "na.impute",
                          ntree = 100, block.size = 1)
  # rf_classifier  <- rfsrc(Surv(OS,Overall_Survival_Event) ~ ., data = training.df, importance = "permute",na.action = "na.impute",
  #                         ntree = 100, block.size = 1)
  # 
}
test.set<-read.csv('~/Desktop/publications/H&N/IMPACT-nonIO.csv',stringsAsFactors = F)

test.set$Indel_load<-NA
test.set$APOBEC_Signature<-NA
test.set<-test.set[,colnames(test.set)%in%colnames(df)]
tmp<-test.set
for(k in 1:ncol(tmp)){tmp[,k]<-test.set[,colnames(test.set)%in%colnames(df)[k]]}
test.set<-tmp
colnames(test.set)<-colnames(df)
tmp<-rbind(df,test.set)
for(i in 1:ncol(tmp)){
  if(colnames(tmp)[i]%in%c('METASTATIC_VS_RECURRENT','AnyVirus','ALCOHOL_NEVER_EVER','GENDER','TP53','PIK3CA',
                           'Smoking_Signature','APOBEC_Signature','PDL1','TERT.promoter','HLALOH','del9q34.3',
                           "ANTIBIOTICS",'STEROIDS','TOXICITY','ECOG','INFECTION_DURING_IO')){
    tmp[,i]<-as.factor(tmp[,i])
  } else if(colnames(tmp)[i]%in%c('HED','OS','PFS','Albumin','HGB','Platelets','Neut.Mono.Lymp','Age_at_IO_start','ploidy',
                                  'purity','Indel_load','nonsyn.TMB')) {
    tmp[,i]<-as.numeric(tmp[,i])
  }
}
test.set<-tmp[-c(1:nrow(df)),]
prediction_for_roc_curve <- predict(rf_classifier,as.data.frame(test.set[,-K]),na.action="na.impute")

V<-prediction_for_roc_curve
x<-apply(X = V$survival, MARGIN = 1, FUN = function(W) V$time.interest[which.min(abs(W-0.5))]) 
#x<-apply(X = V$survival, MARGIN = 1, FUN = function(W) sum(W)) 
#x<-prediction_for_roc_curve
###
clin.fac1<-test.set
clin.fac1$RF14<-x
clin.fac1$BEST_RESPONSE<-as.factor(ifelse(read.csv('~/Desktop/publications/H&N/IMPACT-nonIO.csv')$BEST_RESPONSE=='CR' | read.csv('~/Desktop/publications/H&N/IMPACT-nonIO.csv')$BEST_RESPONSE=='PR', 1,0 ))
# RRF14<-as.numeric(roc( clin.fac1$BEST_RESPONSE,clin.fac1$RF14,smoothed = TRUE)$auc)
# RTMB<-as.numeric(roc( clin.fac1$BEST_RESPONSE,clin.fac1$nonsyn.TMB,smoothed = TRUE)$auc)

clin.fac1$Overall_Survival_Event<-read.csv('~/Desktop/publications/H&N/IMPACT-nonIO.csv')$Overall_Survival_Event
clin.fac1$OS<-read.csv('~/Desktop/publications/H&N/IMPACT-nonIO.csv')$OS*30
clin.fac1$Progression<-read.csv('~/Desktop/publications/H&N/IMPACT-nonIO.csv')$Progression
clin.fac1$PFS<-read.csv('~/Desktop/publications/H&N/IMPACT-nonIO.csv')$PFS*30

#par(mfrow=c(1,2))

# C.idx.24<-concordance.index(x=-clin.fac1$RF14, surv.time=clin.fac1$PFS, surv.event=clin.fac1$Progression, method="noether")$c.index
# C.idx.TMB<-concordance.index(x=-clin.fac1$nonsyn.TMB, surv.time=clin.fac1$PFS, surv.event=clin.fac1$Progression, method="noether")$c.index
# barplot(c(C.idx.24,C.idx.TMB), ylim=c(0,1), col=c('red','green'), names=c('RF14','nonsyn TMB'), ylab='C-index',main = 'PFS', las=2, cex.names = 0.8)

C.idx.24<-concordance.index(x=-clin.fac1$RF14, surv.time=clin.fac1$OS, surv.event=clin.fac1$Overall_Survival_Event, method="noether")$c.index
C.idx.TMB<-concordance.index(x=-clin.fac1$nonsyn.TMB, surv.time=clin.fac1$OS, surv.event=clin.fac1$Overall_Survival_Event, method="noether")$c.index
barplot(c(C.idx.24,C.idx.TMB), ylim=c(0,1), col=c('red','green'), names=c('RF14','nonsyn TMB'), ylab='C-index',main = 'OS', las=2, cex.names = 0.8)

#### ROC ####
## Note the MASS package masks select()!
library(tidyverse)
## https://github.com/tidyverse/tibble/issues/395
options(crayon.enabled = FALSE)
## Used for the dataset.
library(survival)
## Used for visualizaiton.
library(survminer)
## Load the Ovarian Cancer Survival Data

FP.idx<-seq(0.01,0.99,length.out = 1000)

library(timeROC)
library(survivalROC)

## Define a helper functio nto evaluate at various t
#dev.off()

#### OS ROC  ####
t=30*12; my.title = paste('OS at',t/30,'months')

coxph1 <- coxph(Surv(OS, Overall_Survival_Event) ~ clin.fac1$RF14, data = clin.fac1)
HR24<-summary(coxph1)$coefficients[2]
P24<-summary(coxph1)$coefficients[5]

clin.fac1$lp <- predict(coxph1, type = "lp")
V<-timeROC(T=clin.fac1$OS,delta=clin.fac1$Overall_Survival_Event,
           marker=clin.fac1$lp,cause=1,
           weighting="marginal",
           times=c(t),ROC=TRUE)
RF14.TP<-predict(loess(V$TP[,2]~V$FP[,2]), newdata = FP.idx)
AUC1<-round(V$AUC[2],2)
RF14.AUC<-AUC1
plot(V,time=t,col="red",lwd=2,add=F, title=FALSE)
title(my.title)

# V<-survivalROC(Stime        = clin.fac1$PFS,
#             status       = clin.fac1$Progression,
#             marker       = clin.fac1$lp,
#             predict.time = t,
#             method       = "NNE",
#             span = 0.25 * nrow(clin.fac1)^(-0.20))
# AUC1<-round(V$AUC,2)
# plot(V$FP,V$TP,col="red",lwd=2,type='l')

coxph1 <- coxph(Surv(OS, Overall_Survival_Event) ~ clin.fac1$nonsyn.TMB, data = clin.fac1)
HRTMB<-summary(coxph1)$coefficients[2]
PTMB<-summary(coxph1)$coefficients[5]
print(summary(coxph1)$coefficients[5])
clin.fac1$lp <- -1*predict(coxph1, type = "lp")
V<-timeROC(T=clin.fac1$OS,delta=clin.fac1$Overall_Survival_Event,
           marker=clin.fac1$lp,cause=1,
           weighting="marginal",
           times=c(t),ROC=TRUE)
TMB.TP<-predict(loess(V$TP[,2]~V$FP[,2]), newdata = FP.idx)

AUC3<-round(V$AUC[2],2)
TMB.AUC<-AUC3
plot(V,time=t,col="green",add=TRUE,title=FALSE,lwd=2)


legend('topleft',col=c('red','green'),lwd=2,y.intersp=1.5,c(paste0('RF14\n (AUC=',AUC1,')'),
                                                          paste0('nonsyn TMB\n (AUC=',AUC3,')')),bty = 'n',cex=0.8)

#### RF14 ####
clin.fac1$RF14=as.factor(ifelse(clin.fac1$RF14>median(clin.fac1$RF14),'Predicted R','Predicted NR'))
HR.HPV<-round(summary(coxph(Surv(OS, Overall_Survival_Event) ~ RF14, data = clin.fac1))$coefficients[2],2)
fit.HPV <- survfit(Surv(OS, Overall_Survival_Event) ~ RF14, data = clin.fac1)
ggsurv.HPV<-ggsurvplot(fit.HPV, data = clin.fac1, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(150,0.4), ylab = 'OS', xlab = 'Time (years)', xlim = c(0, 365*4), xscale = "d_y",  break.x.by = 365,  palette = c("#E7B800", "#2E9FDF"))
ggsurv.HPV$plot <- ggsurv.HPV$plot + ggplot2::annotate("text", x = 530, y = 0.3, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV), size = 5)
ggsurv.HPV

#### TMB ####
clin.fac1$nonsyn.TMB<-as.factor(ifelse(clin.fac1$nonsyn.TMB<3.34,"High","Low"))
levels(clin.fac1$nonsyn.TMB)<-c('Low','High')

HR.HPV<- round(summary(coxph(Surv(OS, Overall_Survival_Event) ~ nonsyn.TMB, data = clin.fac1))$coefficients[2],2)
fit.HPV <- survfit(Surv(OS, Overall_Survival_Event) ~ nonsyn.TMB, data = clin.fac1)
ggsurv.HPV<-ggsurvplot(fit.HPV, data = clin.fac1, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(150,0.4), ylab = 'OS', xlab = 'Time (years)', xlim = c(0, 365*4), xscale = "d_y",  break.x.by = 365,  palette = c("#E7B800", "#2E9FDF"))
ggsurv.HPV$plot <- ggsurv.HPV$plot + ggplot2::annotate("text", x = 530, y = 0.3, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV), size = 5)
ggsurv.HPV





####
require(rpart)
for(i in 1:nrow(df)){
  df$smoker[i]<-WES$SMOKE_NEVER_EVER[WES$sample%in%rownames(df)[i]]
}
pfit2  <- rpart(Surv(PFS,Progression) ~  nonsyn.TMB + Neut.Mono.Lymp + Platelets + Smoking_Signature + ECOG , 
                data = df)
par(mar = rep(3.5, 4))

pfit2 <- prune(pfit2, cp = 0.05)
plot(pfit2, uniform = TRUE, compress = TRUE)
text(pfit2, use.n = TRUE,cex=0.6)
temp <- snip.rpart(pfit2, 3)
km <- survfit(Surv(PFS,Progression) ~ temp$where, df)
plot(km, lty = 1:3, mark.time = FALSE,
       xlab = "days", ylab = "Progression")
legend(10, 0.3, lty = 1:3)
