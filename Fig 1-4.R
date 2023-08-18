#### Analysis of MSK-ILMN collaboration: WES of Head and Neck cancer patients treated with IO ########
Copy right reserved by Mahdi Golkaram, PhD.


library(ggfortify)
library(survival)
library(survminer)
library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)
library(matrixStats)
library(gplots)
library(finalfit)
library(boot)
require(heatmap.plus)
require(matrixTests)
require(maftools)
require(survival)
require(ggsurv)
library(gridExtra)
require(GenVisR)

#### read data ####
WES<-read.csv('~/PATH2FILE/H&N/WES-HNSC-features.csv')
WES$EBV[is.na(WES$EBV)]<-'not tested'
WES$EBV<-ifelse(WES$EBV=='positive' & WES$HPV_Final_2==0,'EBV+',
                ifelse(WES$HPV_Final_2==1 & WES$EBV=='positive','HPV+/EBV+',
                       ifelse(WES$HPV_Final_2==1 & WES$EBV!='positive','HPV+','negative')))
WES$Number_sites_mets[WES$Number_sites_mets=="N/A"]<-0
WES$Number_sites_mets<-as.numeric(WES$Number_sites_mets)
WES$METASTATIC_VS_RECURRENT<-ifelse(WES$METASTATIC_VS_RECURRENT=="Metastatic","Metastatic","Non-metastatic")
clinicalData<-WES
clinicalData$HLALOH<-ifelse(clinicalData$HLALOH==1,'LOH','no LOH')
colnames(clinicalData)[1]<-'Tumor_Sample_Barcode'
laml.mutsig <- read.delim('~/PATH2FILE/H&N/mutsig.txt')
TCGA<-read.delim('~/PATH2FILE/H&N/IntOGen-DriverGenes_HNSC_TCGA.tsv')
clinicalData$TERT_promoter<-ifelse(clinicalData$TERT.promoter==1,'mut','wt')
clinicalData$HPV_status<-ifelse(clinicalData$HPV_Final_2==1,'positive','negative')
clinicalData$Virus_status<-clinicalData$EBV
clinicalData$PRIMARY_TUMOR_SITE<-clinicalData$PRIMARY_TUMOR_SITE
laml<-read.maf(maf = '~/PATH2FILE/H&N/MAF/allsamples_all.vep-complete.maf',
               #cnTable = cn,
               clinicalData = clinicalData,
               verbose = FALSE)

##### Figure 1 #####
subset.samples<-as.character(getSampleSummary(laml)$Tumor_Sample_Barcode)[as.character(getSampleSummary(laml)$Tumor_Sample_Barcode)%in%WES$sample]
laml<-subsetMaf(laml,tsb = subset.samples, mafObj = T)
#HPV.pos.samples<-as.character(getSampleSummary(laml)$Tumor_Sample_Barcode)[as.character(getSampleSummary(laml)$Tumor_Sample_Barcode)%in%(WES$sample[WES$HPV.status=='Positive'])]
HPV.pos.samples<-as.character(getSampleSummary(laml)$Tumor_Sample_Barcode)[as.character(getSampleSummary(laml)$Tumor_Sample_Barcode)%in%(WES$sample[WES$AnyVirus==1])]
#HPV.neg.samples<-as.character(getSampleSummary(laml)$Tumor_Sample_Barcode)[as.character(getSampleSummary(laml)$Tumor_Sample_Barcode)%in%(WES$sample[WES$HPV.status!='Positive'])]
HPV.neg.samples<-as.character(getSampleSummary(laml)$Tumor_Sample_Barcode)[as.character(getSampleSummary(laml)$Tumor_Sample_Barcode)%in%(WES$sample[WES$AnyVirus!=1])]

laml.HPV.pos<-subsetMaf(laml,tsb = HPV.pos.samples, mafObj = T)
laml.HPV.neg<-subsetMaf(laml,tsb = HPV.neg.samples, mafObj = T)
laml.mutsig[,c(19)]<- -10*log10(laml.mutsig[,c(19)])
q=data.table::data.table(laml.mutsig[,c(2,19)])
M<-mutCountMatrix(maf = laml)
M<-M[,order(colnames(M))]
#feq.HPV.pos<-rowSums(M[,WES$HPV.status=='Positive'])
feq.HPV.pos<-rowSums(M[,WES$AnyVirus==1])
#feq.HPV.neg<-rowSums(M[,WES$HPV.status!='Positive'])
feq.HPV.neg<-rowSums(M[,WES$HPV.status!=0])

library(BSgenome.Hsapiens.UCSC.hg19, quietly = TRUE)
library('NMF')
laml.tnm = trinucleotideMatrix(maf = laml,add = F, ignoreChr = " chr", ref_genome = "BSgenome.Hsapiens.UCSC.hg19")
#laml.sign = estimateSignatures(mat = laml.tnm, nTry = 6)
#plotCophenetic(res = laml.sign)
laml.sig = extractSignatures(mat = laml.tnm, n = 4) 
#maftools::plotSignatures(nmfRes = laml.sig, title_size = 0.8, sig_db = "SBS")
#barplot(laml.sig$contributions[,WES$HPV.status=='Positive'], col = 2:7, main = 'signature contribtuion (HPV+)',las =2, cex.names = 0.4)
barplot(laml.sig$contributions[,WES$AnyVirus==1], col = 2:7, main = 'signature contribtuion (V+)',las =2, cex.names = 0.4)
#barplot(laml.sig$contributions[,WES$HPV.status!='Positive'], col = 2:7, main = 'signature contribtuion (V-)',las =2, cex.names = 0.4)
barplot(laml.sig$contributions[,WES$AnyVirus!=1], col = 2:7, main = 'signature contribtuion (V-)',las =2, cex.names = 0.4)
legend('topright', legend = c('APOBEC (SBS13)','APOBEC (SBS2)',
                              'spontaneous or enzymatic deamination of 5-methylcytosine (Age) (SBS1)',
                              'exposure to tobacco (smoking) mutagens (SBS4)'),
       col = 2:7,pch = 15)



##### Figure 1 #####
my_palette <- colorRampPalette(c("gray","darkslategrey"))(n = 1000)
barplot(t(cbind(TCGA[1:10,4],as.numeric(feq.HPV.pos[TCGA$Symbol[1:10]]),as.numeric(feq.HPV.neg[TCGA$Symbol[1:10]])))/100,beside = T, ylim=c(0,1), 
        names=TCGA$Symbol[1:10],las =2, ylab='Mutation frequency')
legend('topleft',c('TCGA','This study: V+', 'This study: V-'), bty='n', pch = 15, col = c(my_palette)[c(911,310,1)])
par(mfrow = c(1, 3))
barplot(cbind(mean(WES$APOBEC_Signature[WES$AnyVirus==1]=='yes'),
              mean(WES$APOBEC_Signature[WES$AnyVirus==0]=='yes')),
        names = c('V+', 'V-'), ylab = 'APOBEC signature frequency',
        ylim=c(0,1),las =2)
legend('top',bty = 'n',legend = paste0('P=', round(fisher.test(matrix(c(as.numeric(table((WES$APOBEC_Signature[WES$AnyVirus==0]=='yes'))),
                                                                        as.numeric(table((WES$APOBEC_Signature[WES$AnyVirus==1]=='yes')))),2,2))$p.val,3)))
barplot(cbind(mean(WES$Smoking_Signature[WES$AnyVirus==1]=='yes'),
              mean(WES$Smoking_Signature[WES$AnyVirus==0]=='yes')),
        names = c('V+', 'V-'), ylab = 'Smoking signature frequency',
        ylim=c(0,1),las =2)
legend('top',bty = 'n',legend = paste0('P=', round(fisher.test(matrix(c(as.numeric(table((WES$Smoking_Signature[WES$AnyVirus==0]=='yes'))),
                                                                        as.numeric(table((WES$Smoking_Signature[WES$AnyVirus==1]=='yes')))),2,2))$p.val,5)))

barplot(cbind(mean(WES$BEST_RESPONSE[WES$AnyVirus==1]=='PR' | WES$BEST_RESPONSE[WES$AnyVirus==1]=='CR'),
              mean(WES$BEST_RESPONSE[WES$AnyVirus==0]=='PR' | WES$BEST_RESPONSE[WES$AnyVirus==0]=='CR')),
        names = c('V+', 'V-'), ylab = 'Response rate',
        ylim=c(0,1),las =2)
legend('top',bty = 'n',legend = paste0('P=', round(fisher.test(matrix(c(as.numeric(table((WES$BEST_RESPONSE[WES$AnyVirus==0]=='PR' | WES$BEST_RESPONSE[WES$AnyVirus==0]=='CR'))),
                                                                        as.numeric(table((WES$BEST_RESPONSE[WES$AnyVirus==1]=='PR' | WES$BEST_RESPONSE[WES$AnyVirus==1]=='CR')))),2,2))$p.val,2)))

dev.off()

p=data.table::data.table(TCGA[,1:2])
p$Mutations<-1
colnames(p)<-c('genes','Driver mutation')

HLAcolors=c("#D53E4F", "#F46D43")
names(HLAcolors)=c("LOH","no LOH")

sigcolors = c("#FDAE61", "#FEE08B")
names(sigcolors)=c("yes","no")

MSIcolors=c("green", "brown")
names(MSIcolors)=c("MSH","MSS")

Responsecolors=c("#3288BD","darkgreen","orange", "red")
names(Responsecolors)=c("CR","PR","SD","PD")

CNLcolors=c("purple", "pink")
names(CNLcolors)=c(1,0)

CNGcolors=c("cyan", "blue")
names(CNGcolors)=c(1,0)

CNVcolors=c("purple", "green",'pink')
names(CNVcolors)=c(-1,0,1)

TERTcolors=c("cyan", "blue")
names(TERTcolors)=c('mut','wt')

Vcolors=c("purple", "pink")
names(Vcolors)=c('HPV+','EBV+')


primarycolors=c('gray','maroon','darkcyan','burlywood1','black','aquamarine1','magenta','#640707')
names(primarycolors)=names(table(WES$PRIMARY_TUMOR_SITE))
#This is the final colored list. Names of the list elements should match those in clinicalFeatures arguments 
anno_cols = list(HLALOH = HLAcolors,
                 Smoking_Signature = sigcolors,
                 APOBEC_Signature = sigcolors,
                 BEST_RESPONSE = Responsecolors,
                 MSI = MSIcolors,
                 HPV = Vcolors,
                 EBV = Vcolors,
                 CDKN2A = CNLcolors,
                 PTEN = CNLcolors,
                 q34.3 = CNLcolors,
                 CCND1 = CNGcolors,
                 PRIMARY_TUMOR_SITE=primarycolors,
                 HPV_status=Vcolors,
                 Virus_status=Vcolors,
                 #PDL1 = CNVcolors
                 TERT_promoter = TERTcolors
                 #  ALCOHOL_NEVER_EVER = ALCOHOL_NEVER_EVERcolors
)
topBarData = as.data.frame(cbind(laml.HPV.neg@clinical.data$Tumor_Sample_Barcode,
                                 laml.HPV.neg@clinical.data$nonsyn.TMB))
#topBarData$V2[which.max(topBarData$V2)]<-20
colnames(topBarData)<-c('Tumor_Sample_Barcode','TMB')
A<-oncoplot(showTumorSampleBarcodes = F, rightBarLims = c(0, 400),maf = laml.HPV.neg,rightBarData =q,genes = laml.mutsig$gene[1:20], fill = TRUE,removeNonMutated = F,sortByMutation = T,
         clinicalFeatures=c('TERT_promoter',
                            'HLALOH','MSI',
                            'APOBEC_Signature','Smoking_Signature','PRIMARY_TUMOR_SITE',
                            'BEST_RESPONSE'),anno_height=5,legendFontSize = 0.55,
         annotationFontSize = 0.55,leftBarData = p,leftBarLims = c(0,100),annotationColor=anno_cols,topBarData=topBarData)
topBarData = as.data.frame(cbind(laml.HPV.pos@clinical.data$Tumor_Sample_Barcode,
                                 laml.HPV.pos@clinical.data$nonsyn.TMB))
topBarData$V2[which.max(topBarData$V2)]<-20
colnames(topBarData)<-c('Tumor_Sample_Barcode','TMB')
B<-oncoplot(rightBarLims = c(0, 400),maf = laml.HPV.pos,rightBarData =q,genes = laml.mutsig$gene[1:20], fill = TRUE,removeNonMutated = F,sortByMutation = T,
         clinicalFeatures=c('TERT_promoter',
                            'HLALOH','MSI',
                            'APOBEC_Signature','Smoking_Signature',
                            'Virus_status','PRIMARY_TUMOR_SITE',
                            'BEST_RESPONSE'),anno_height=5,legendFontSize = 0.55,
         annotationFontSize = 0.55,leftBarData = p,leftBarLims = c(0,100),annotationColor=anno_cols,topBarData=topBarData)

#### reorder 1B ####
dev.off()
par(mfrow=c(2,1))
barplot((laml.sig$contributions[,WES$AnyVirus==1])[,B], col = 2:7, main = 'signature contribtuion (V+)',las =2, cex.names = 0.4)
#barplot(laml.sig$contributions[,WES$HPV.status!='Positive'], col = 2:7, main = 'signature contribtuion (V-)',las =2, cex.names = 0.4)
barplot((laml.sig$contributions[,WES$AnyVirus!=1])[,A], col = 2:7, main = 'signature contribtuion (V-)',las =2, cex.names = 0.4)
legend('topright', legend = c('APOBEC (SBS13)','APOBEC (SBS2)',
                              'spontaneous or enzymatic deamination of 5-methylcytosine (Age) (SBS1)',
                              'exposure to tobacco (smoking) mutagens (SBS4)'),
       col = 2:7,pch = 15)

####



G<-read.csv('~/PATH2FILE/H&N/CNV_impact.csv')
G$ID<-gsub('-T.*','',G$ID)
G<-G[G$ID%in%WES$Patient.ID,]
G$loc<-paste0('chr',G$chrom,':',G$loc.start,'-',G$loc.end)
g<-G[,c(2,3,4,6,1)]
colnames(g)<-c("chromosome", "start", "end", "segmean","sample")
g$segmean<-2^(g$segmean+1)
g.neg<-g[g$sample%in%WES$Patient.ID[ WES$AnyVirus == 0 ],]
g.pos<-g[g$sample%in%WES$Patient.ID[ WES$AnyVirus == 1 ],]
cnFreq(g.neg, genome="hg19",CN_low_cutoff = 1.86,CN_high_cutoff = 2.14,facet_lab_size=2)
cnFreq(g.pos, genome="hg19",CN_low_cutoff = 1.86,CN_high_cutoff = 2.14,facet_lab_size=2)


my_palette <- colorRampPalette(c("yellow","purple"))(n = 2)

par(mfrow = c(1, 7))
p<-round(wilcox.test(WES$nonsyn.TMB[WES$AnyVirus==1],WES$nonsyn.TMB[WES$AnyVirus==0])$p.value,2)
boxplot(WES$nonsyn.TMB~WES$AnyVirus, col = my_palette, outline = F, names=c('V-', 'V+'), ylab='Nonsyn mutational burden [mut/Mbp]', xlab ='',main=paste0('wilcox P=',p))
stripchart(WES$nonsyn.TMB~WES$AnyVirus,
           vertical = TRUE, method = "jitter",
           pch = 21, col = "blue", bg = "bisque",
           add = TRUE, cex = 0.5)

p<-round(wilcox.test(WES$Indel_load[WES$AnyVirus==1],WES$Indel_load[WES$AnyVirus==0])$p.value,2)
boxplot(WES$Indel_load/30~WES$AnyVirus, col = my_palette, outline = F, names=c('V-', 'V+'), ylab='Indel mutational burden [mut/Mbp]', xlab ='',main=paste0('wilcox P=',p))
stripchart(WES$Indel_load/30~WES$AnyVirus,
           vertical = TRUE, method = "jitter",
           pch = 21, col = "blue", bg = "bisque",
           add = TRUE, cex = 0.5)

p<-round(wilcox.test(WES$Clonal_Mutational_Load[WES$AnyVirus==1],WES$Clonal_Mutational_Load[WES$AnyVirus==0])$p.value,2)
boxplot(WES$Clonal_Mutational_Load/30~WES$AnyVirus, col = my_palette, outline = F, names=c('V-', 'V+'), ylab='Clonal mutational burden [mut/Mbp]', xlab ='',main=paste0('wilcox P=',p))
stripchart(WES$Clonal_Mutational_Load/30~WES$AnyVirus,
           vertical = TRUE, method = "jitter",
           pch = 21, col = "blue", bg = "bisque",
           add = TRUE, cex = 0.5)

p<-round(wilcox.test(WES$dN.dS.ON.OFF[WES$AnyVirus==1],WES$dN.dS.ON.OFF[WES$AnyVirus==0])$p.value,2)
boxplot(WES$dN.dS.ON.OFF~WES$AnyVirus, col = my_palette, outline = F, names=c('V-', 'V+'), ylab='dN/dS (ON/OFF)', xlab ='',main=paste0('wilcox P=',p))
stripchart(WES$dN.dS.ON.OFF~WES$AnyVirus,
           vertical = TRUE, method = "jitter",
           pch = 21, col = "blue", bg = "bisque",
           add = TRUE, cex = 0.5)

p<-round(wilcox.test(WES$ITH[WES$AnyVirus==1],WES$ITH[WES$AnyVirus==0])$p.value,2)
boxplot(WES$ITH~WES$AnyVirus, col = my_palette, outline = F, names=c('V-', 'V+'), ylab='ITH', xlab ='',main=paste0('wilcox P=',p))
stripchart(WES$ITH~WES$AnyVirus,
           vertical = TRUE, method = "jitter",
           pch = 21, col = "blue", bg = "bisque",
           add = TRUE, cex = 0.5)

p<-round(wilcox.test(WES$ploidy[WES$AnyVirus==1],WES$ploidy[WES$AnyVirus==0])$p.value,4)
boxplot(WES$ploidy~WES$AnyVirus, col = my_palette, outline = F, names=c('V-', 'V+'), ylab='Tumor ploidy', xlab ='',main=paste0('wilcox P=',p))
stripchart(WES$ploidy~WES$AnyVirus,
           vertical = TRUE, method = "jitter",
           pch = 21, col = "blue", bg = "bisque",
           add = TRUE, cex = 0.5)

p<-round(wilcox.test(WES$purity[WES$AnyVirus==1],WES$purity[WES$AnyVirus==0])$p.value,4)
boxplot(WES$purity~WES$AnyVirus, col = my_palette, outline = F, names=c('V-', 'V+'), ylab='Tumor purity', xlab ='',main=paste0('wilcox P=',p))
stripchart(WES$purity~WES$AnyVirus,
           vertical = TRUE, method = "jitter",
           pch = 21, col = "blue", bg = "bisque",
           add = TRUE, cex = 0.5)

dev.off()


clin.fac<-WES
clin.fac$HED<-(clin.fac$tumor.HLA.A+clin.fac$tumor.HLA.B+clin.fac$tumor.HLA.C)/3
clin.fac$OS<-clin.fac$OS*30; clin.fac$PFS<-clin.fac$PFS*30; clin.fac$DSS<-clin.fac$DSS*30
clin.fac$AnyVirus=ifelse(WES$AnyVirus==1,'Positive','Negative')
clin.fac$ALCOHOL_NEVER_EVER=ifelse(WES$ALCOHOL_NEVER_EVER==1,'Yes','No')
clin.fac$ECOG=as.factor(clin.fac$ECOG)
clin.fac$del9q34.3=as.factor(ifelse(WES$del9q34.3==1,'del','wt'))
clin.fac$MET=as.factor(ifelse(WES$MET.amp==1,'amp','wt'))
clin.fac$TERT.promoter=as.factor(ifelse(WES$TERT.promoter==1,'mut','wt'))
clin.fac$HLALOH=as.factor(ifelse(WES$HLALOH==1,'LOH','No LOH'))
#### V ####
HR.HPV<-round(summary(coxph(Surv(PFS, Progression) ~ AnyVirus, data = clin.fac))$coefficients[2],2)
fit.HPV <- survfit(Surv(PFS, Progression) ~ AnyVirus, data = clin.fac)
ggsurv.HPV<-ggsurvplot(fit.HPV, data = clin.fac, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (years)', xlim = c(0, 365*4), xscale = "d_y",  break.x.by = 365,  palette = c("#E7B800", "#2E9FDF"))
ggsurv.HPV$plot <- ggsurv.HPV$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV), size = 5)
ggsurv.HPV
HR.HPV<-round(summary(coxph(Surv(OS, Overall_Survival_Event) ~ AnyVirus, data = clin.fac))$coefficients[2],2)
fit.HPV <- survfit(Surv(OS, Overall_Survival_Event) ~ AnyVirus, data = clin.fac)
ggsurv.HPV<-ggsurvplot(fit.HPV, data = clin.fac, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'OS', xlab = 'Time (years)', xlim = c(0, 365*4), xscale = "d_y",  break.x.by = 365,  palette = c("#E7B800", "#2E9FDF"))
ggsurv.HPV$plot <- ggsurv.HPV$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV), size = 5)
ggsurv.HPV

#### V #####
clin.fac$V=clin.fac$EBV
HR.HPV<-round(summary(coxph(Surv(PFS, Progression) ~ V, data = clin.fac))$coefficients[2],2)
fit.HPV <- survfit(Surv(PFS, Progression) ~ V, data = clin.fac)
ggsurv.HPV<-ggsurvplot(fit.HPV, data = clin.fac, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (years)', xlim = c(0, 365*4), xscale = "d_y",  break.x.by = 365)
ggsurv.HPV$plot <- ggsurv.HPV$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV), size = 5)
ggsurv.HPV
HR.HPV<-round(summary(coxph(Surv(OS, Overall_Survival_Event) ~ V, data = clin.fac))$coefficients[2],2)
fit.HPV <- survfit(Surv(OS, Overall_Survival_Event) ~ V, data = clin.fac)
ggsurv.HPV<-ggsurvplot(fit.HPV, data = clin.fac, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'OS', xlab = 'Time (years)', xlim = c(0, 365*4), xscale = "d_y",  break.x.by = 365)
ggsurv.HPV$plot <- ggsurv.HPV$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV), size = 5)
ggsurv.HPV

my_palette <- colorRampPalette(c("yellow","purple"))(n = 3)
round(wilcox.test(WES$nonsyn.TMB[WES$EBV=='EBV+'],WES$nonsyn.TMB[WES$EBV=='HPV+'])$p.value,2)
round(wilcox.test(WES$nonsyn.TMB[WES$EBV=='EBV+'],WES$nonsyn.TMB[WES$EBV=='negative'])$p.value,2)
boxplot(WES$nonsyn.TMB~WES$EBV, col = my_palette, outline = F, names=c('EBV+','HPV+', 'V-'), ylab='Nonsyn mutational burden [mut/Mbp]', xlab ='',ylim=c(0,15))
stripchart(WES$nonsyn.TMB~WES$EBV,
           vertical = TRUE, method = "jitter",
           pch = 21, col = "blue", bg = "bisque",
           add = TRUE, cex = 0.5)
barplot(cbind(mean(WES$BEST_RESPONSE[WES$EBV=='EBV+']=='PR' | WES$BEST_RESPONSE[WES$EBV=='EBV+']=='CR'),
              mean(WES$BEST_RESPONSE[WES$EBV=='HPV+']=='PR' | WES$BEST_RESPONSE[WES$EBV=='HPV+']=='CR'),
              mean(WES$BEST_RESPONSE[WES$EBV=='negative']=='PR' | WES$BEST_RESPONSE[WES$EBV=='negative']=='CR')),
        names = c('EBV+','HPV+','V-'), ylab = 'Response rate',
        ylim=c(0,1),las =2)
legend('top',bty = 'n',legend = paste0('P=', round(fisher.test(matrix(c(as.numeric(table((WES$BEST_RESPONSE[WES$AnyVirus==0]=='PR' | WES$BEST_RESPONSE[WES$AnyVirus==0]=='CR'))),
                                                                        as.numeric(table((WES$BEST_RESPONSE[WES$AnyVirus==1]=='PR' | WES$BEST_RESPONSE[WES$AnyVirus==1]=='CR')))),2,2))$p.val,2)))





#### threshold setting ####
a<-as.numeric(quantile(WES$purity,probs = 0.1,na.rm = T));b<-as.numeric(quantile(WES$purity,probs = 0.9,na.rm = T))

for(i in  seq(a,b,length.out = 50)){
  clin.fac.HPV.pos<-clin.fac[clin.fac$AnyVirus=='Positive',]
  clin.fac.HPV.neg<-clin.fac[clin.fac$AnyVirus=='Negative',]
  clin.fac.HPV.pos$purity<-as.factor(ifelse(clin.fac.HPV.pos$purity>i,'high','low'))
  clin.fac.HPV.neg$purity<-as.factor(ifelse(clin.fac.HPV.neg$purity>i,'high','low'))
  print(paste('i=',i))
  print(round(summary(coxph(Surv(PFS, Progression) ~ clin.fac.HPV.pos$purity, data = clin.fac.HPV.pos))$coefficients[5],4))
  print(round(summary(coxph(Surv(PFS, Progression) ~ clin.fac.HPV.neg$purity, data = clin.fac.HPV.neg))$coefficients[5],4))
}
###
clin.fac.HPV.pos<-clin.fac[clin.fac$AnyVirus=='Positive',]
clin.fac.HPV.neg<-clin.fac[clin.fac$AnyVirus=='Negative',]

clin.fac.HPV.pos$BMI<-as.factor(ifelse(clin.fac.HPV.pos$BMI>25,'>25','<25'))
clin.fac.HPV.neg$BMI<-as.factor(ifelse(clin.fac.HPV.neg$BMI>25,'>25','<25'))

clin.fac.HPV.pos$Albumin<-as.factor(ifelse(clin.fac.HPV.pos$Albumin>3.75,'>3.75','<3.75'))
clin.fac.HPV.neg$Albumin<-as.factor(ifelse(clin.fac.HPV.neg$Albumin>3.75,'>3.75','<3.75'))

clin.fac.HPV.pos$HGB<-as.factor(ifelse(clin.fac.HPV.pos$HGB>9.6,'>9.6','<9.6'))
clin.fac.HPV.neg$HGB<-as.factor(ifelse(clin.fac.HPV.neg$HGB>9.6,'>9.6','<9.6'))

clin.fac.HPV.pos$Platelets<-as.factor(ifelse(clin.fac.HPV.pos$Platelets>200,'>200','<200'))
clin.fac.HPV.neg$Platelets<-as.factor(ifelse(clin.fac.HPV.neg$Platelets>200,'>200','<200'))

clin.fac.HPV.pos$Neut.Mono.Lymp<-as.factor(ifelse(clin.fac.HPV.pos$Neut.Mono.Lymp>8.3,'>8.3','<8.3'))
clin.fac.HPV.neg$Neut.Mono.Lymp<-as.factor(ifelse(clin.fac.HPV.neg$Neut.Mono.Lymp>8.3,'>8.3','<8.3'))

clin.fac.HPV.pos$Age_at_IO_start<-as.factor(ifelse(clin.fac.HPV.pos$Age_at_IO_start>60,'>60','<60'))
clin.fac.HPV.neg$Age_at_IO_start<-as.factor(ifelse(clin.fac.HPV.neg$Age_at_IO_start>60,'>60','<60'))

clin.fac.HPV.pos$Number_sites_mets<-as.factor(ifelse(clin.fac.HPV.pos$Number_sites_mets>2,'>2','<2'))
clin.fac.HPV.neg$Number_sites_mets<-as.factor(ifelse(clin.fac.HPV.neg$Number_sites_mets>2,'>2','<2'))

clin.fac.HPV.pos$nonsyn.TMB<-as.factor(ifelse(clin.fac.HPV.pos$nonsyn.TMB>3.34,'high','low'))
clin.fac.HPV.neg$nonsyn.TMB<-as.factor(ifelse(clin.fac.HPV.neg$nonsyn.TMB>3.34,'high','low'))

clin.fac.HPV.pos$Indel_load<-as.factor(ifelse(clin.fac.HPV.pos$Indel_load>5,'high','low'))
clin.fac.HPV.neg$Indel_load<-as.factor(ifelse(clin.fac.HPV.neg$Indel_load>5,'high','low'))

clin.fac.HPV.pos$Clonal_Mutational_Load<-as.factor(ifelse(clin.fac.HPV.pos$Clonal_Mutational_Load>103.7,'high','low'))
clin.fac.HPV.neg$Clonal_Mutational_Load<-as.factor(ifelse(clin.fac.HPV.neg$Clonal_Mutational_Load>103.7,'high','low'))

clin.fac.HPV.pos$dN.dS.ON.OFF<-as.factor(ifelse(clin.fac.HPV.pos$dN.dS.ON.OFF>1.6,'high','low'))
clin.fac.HPV.neg$dN.dS.ON.OFF<-as.factor(ifelse(clin.fac.HPV.neg$dN.dS.ON.OFF>1.6,'high','low'))

clin.fac.HPV.pos$tumor.HLA.A<-as.factor(ifelse(clin.fac.HPV.pos$tumor.HLA.A>15,'high','low'))
clin.fac.HPV.neg$tumor.HLA.A<-as.factor(ifelse(clin.fac.HPV.neg$tumor.HLA.A>15,'high','low'))

clin.fac.HPV.pos$tumor.HLA.B<-as.factor(ifelse(clin.fac.HPV.pos$tumor.HLA.B>15,'high','low'))
clin.fac.HPV.neg$tumor.HLA.B<-as.factor(ifelse(clin.fac.HPV.neg$tumor.HLA.B>15,'high','low'))

clin.fac.HPV.pos$tumor.HLA.C<-as.factor(ifelse(clin.fac.HPV.pos$tumor.HLA.C>15,'high','low'))
clin.fac.HPV.neg$tumor.HLA.C<-as.factor(ifelse(clin.fac.HPV.neg$tumor.HLA.C>15,'high','low'))

clin.fac.HPV.pos$ploidy<-as.factor(ifelse(clin.fac.HPV.pos$ploidy>2.6,'high','low'))
clin.fac.HPV.neg$ploidy<-as.factor(ifelse(clin.fac.HPV.neg$ploidy>2.6,'high','low'))

clin.fac.HPV.pos$purity<-as.factor(ifelse(clin.fac.HPV.pos$purity>0.54,'high','low'))
clin.fac.HPV.neg$purity<-as.factor(ifelse(clin.fac.HPV.neg$purity>0.54,'high','low'))

clin.fac.HPV.pos$ITH<-as.factor(ifelse(clin.fac.HPV.pos$ITH>25,'high','low'))
clin.fac.HPV.neg$ITH<-as.factor(ifelse(clin.fac.HPV.neg$ITH>25,'high','low'))

#### HLALOH ####
HR.HPV.pos<-round(summary(coxph(Surv(PFS, Progression) ~ HLALOH, data = clin.fac.HPV.pos))$coefficients[2],2)
HR.HPV.neg<-round(summary(coxph(Surv(PFS, Progression) ~ HLALOH, data = clin.fac.HPV.neg))$coefficients[2],2)
fit.HPV.pos <- survfit(Surv(PFS, Progression) ~ HLALOH, data = clin.fac.HPV.pos)
fit.HPV.neg <- survfit(Surv(PFS, Progression) ~ HLALOH, data = clin.fac.HPV.neg)
ggsurv.HPV.pos<-ggsurvplot(fit.HPV.pos, data = clin.fac.HPV.pos, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.neg<-ggsurvplot(fit.HPV.neg, data = clin.fac.HPV.neg, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.pos$plot <- ggsurv.HPV.pos$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.pos), size = 5)
ggsurv.HPV.neg$plot <- ggsurv.HPV.neg$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.neg), size = 5)
ggsurv.HPV.pos;ggsurv.HPV.neg

#### TP53 ####
HR.HPV.pos<-round(summary(coxph(Surv(PFS, Progression) ~ TP53, data = clin.fac.HPV.pos))$coefficients[2],2)
HR.HPV.neg<-round(summary(coxph(Surv(PFS, Progression) ~ TP53, data = clin.fac.HPV.neg))$coefficients[2],2)
fit.HPV.pos <- survfit(Surv(PFS, Progression) ~ TP53, data = clin.fac.HPV.pos)
fit.HPV.neg <- survfit(Surv(PFS, Progression) ~ TP53, data = clin.fac.HPV.neg)
ggsurv.HPV.pos<-ggsurvplot(fit.HPV.pos, data = clin.fac.HPV.pos, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.neg<-ggsurvplot(fit.HPV.neg, data = clin.fac.HPV.neg, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.pos$plot <- ggsurv.HPV.pos$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.pos), size = 5)
ggsurv.HPV.neg$plot <- ggsurv.HPV.neg$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.neg), size = 5)
ggsurv.HPV.pos;ggsurv.HPV.neg

#### PIK3CA ####
HR.HPV.pos<-round(summary(coxph(Surv(PFS, Progression) ~ PIK3CA, data = clin.fac.HPV.pos))$coefficients[2],2)
HR.HPV.neg<-round(summary(coxph(Surv(PFS, Progression) ~ PIK3CA, data = clin.fac.HPV.neg))$coefficients[2],2)
fit.HPV.pos <- survfit(Surv(PFS, Progression) ~ PIK3CA, data = clin.fac.HPV.pos)
fit.HPV.neg <- survfit(Surv(PFS, Progression) ~ PIK3CA, data = clin.fac.HPV.neg)
ggsurv.HPV.pos<-ggsurvplot(fit.HPV.pos, data = clin.fac.HPV.pos, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.neg<-ggsurvplot(fit.HPV.neg, data = clin.fac.HPV.neg, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.pos$plot <- ggsurv.HPV.pos$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.pos), size = 5)
ggsurv.HPV.neg$plot <- ggsurv.HPV.neg$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.neg), size = 5)
ggsurv.HPV.pos;ggsurv.HPV.neg

#### Smoking_Signature ####
HR.HPV.pos<-round(summary(coxph(Surv(PFS, Progression) ~ Smoking_Signature, data = clin.fac.HPV.pos))$coefficients[2],2)
HR.HPV.neg<-round(summary(coxph(Surv(PFS, Progression) ~ Smoking_Signature, data = clin.fac.HPV.neg))$coefficients[2],2)
fit.HPV.pos <- survfit(Surv(PFS, Progression) ~ Smoking_Signature, data = clin.fac.HPV.pos)
fit.HPV.neg <- survfit(Surv(PFS, Progression) ~ Smoking_Signature, data = clin.fac.HPV.neg)
ggsurv.HPV.pos<-ggsurvplot(fit.HPV.pos, data = clin.fac.HPV.pos, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.neg<-ggsurvplot(fit.HPV.neg, data = clin.fac.HPV.neg, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.pos$plot <- ggsurv.HPV.pos$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.pos), size = 5)
ggsurv.HPV.neg$plot <- ggsurv.HPV.neg$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.neg), size = 5)
ggsurv.HPV.pos;ggsurv.HPV.neg

#### APOBEC_Signature ####
HR.HPV.pos<-round(summary(coxph(Surv(PFS, Progression) ~ APOBEC_Signature, data = clin.fac.HPV.pos))$coefficients[2],2)
HR.HPV.neg<-round(summary(coxph(Surv(PFS, Progression) ~ APOBEC_Signature, data = clin.fac.HPV.neg))$coefficients[2],2)
fit.HPV.pos <- survfit(Surv(PFS, Progression) ~ APOBEC_Signature, data = clin.fac.HPV.pos)
fit.HPV.neg <- survfit(Surv(PFS, Progression) ~ APOBEC_Signature, data = clin.fac.HPV.neg)
ggsurv.HPV.pos<-ggsurvplot(fit.HPV.pos, data = clin.fac.HPV.pos, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.neg<-ggsurvplot(fit.HPV.neg, data = clin.fac.HPV.neg, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.pos$plot <- ggsurv.HPV.pos$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.pos), size = 5)
ggsurv.HPV.neg$plot <- ggsurv.HPV.neg$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.neg), size = 5)
ggsurv.HPV.pos;ggsurv.HPV.neg

#### nonsyn TMB ####
HR.HPV.pos<-round(summary(coxph(Surv(PFS, Progression) ~ nonsyn.TMB, data = clin.fac.HPV.pos))$coefficients[2],2)
HR.HPV.neg<-round(summary(coxph(Surv(PFS, Progression) ~ nonsyn.TMB, data = clin.fac.HPV.neg))$coefficients[2],2)
fit.HPV.pos <- survfit(Surv(PFS, Progression) ~ nonsyn.TMB, data = clin.fac.HPV.pos)
fit.HPV.neg <- survfit(Surv(PFS, Progression) ~ nonsyn.TMB, data = clin.fac.HPV.neg)
ggsurv.HPV.pos<-ggsurvplot(fit.HPV.pos, data = clin.fac.HPV.pos, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.neg<-ggsurvplot(fit.HPV.neg, data = clin.fac.HPV.neg, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.pos$plot <- ggsurv.HPV.pos$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.pos), size = 5)
ggsurv.HPV.neg$plot <- ggsurv.HPV.neg$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.neg), size = 5)
ggsurv.HPV.pos;ggsurv.HPV.neg


#### Indel load ####
HR.HPV.pos<-round(summary(coxph(Surv(PFS, Progression) ~ Indel_load, data = clin.fac.HPV.pos))$coefficients[2],2)
HR.HPV.neg<-round(summary(coxph(Surv(PFS, Progression) ~ Indel_load, data = clin.fac.HPV.neg))$coefficients[2],2)
fit.HPV.pos <- survfit(Surv(PFS, Progression) ~ Indel_load, data = clin.fac.HPV.pos)
fit.HPV.neg <- survfit(Surv(PFS, Progression) ~ Indel_load, data = clin.fac.HPV.neg)
ggsurv.HPV.pos<-ggsurvplot(fit.HPV.pos, data = clin.fac.HPV.pos, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.neg<-ggsurvplot(fit.HPV.neg, data = clin.fac.HPV.neg, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.pos$plot <- ggsurv.HPV.pos$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.pos), size = 5)
ggsurv.HPV.neg$plot <- ggsurv.HPV.neg$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.neg), size = 5)
ggsurv.HPV.pos;ggsurv.HPV.neg

#### Clonal mutational load ####
HR.HPV.pos<-round(summary(coxph(Surv(PFS, Progression) ~ Clonal_Mutational_Load, data = clin.fac.HPV.pos))$coefficients[2],2)
HR.HPV.neg<-round(summary(coxph(Surv(PFS, Progression) ~ Clonal_Mutational_Load, data = clin.fac.HPV.neg))$coefficients[2],2)
fit.HPV.pos <- survfit(Surv(PFS, Progression) ~ Clonal_Mutational_Load, data = clin.fac.HPV.pos)
fit.HPV.neg <- survfit(Surv(PFS, Progression) ~ Clonal_Mutational_Load, data = clin.fac.HPV.neg)
ggsurv.HPV.pos<-ggsurvplot(fit.HPV.pos, data = clin.fac.HPV.pos, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.neg<-ggsurvplot(fit.HPV.neg, data = clin.fac.HPV.neg, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.pos$plot <- ggsurv.HPV.pos$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.pos), size = 5)
ggsurv.HPV.neg$plot <- ggsurv.HPV.neg$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.neg), size = 5)
ggsurv.HPV.pos;ggsurv.HPV.neg
#### dN/dS (ON/OFF) ####
HR.HPV.pos<-round(summary(coxph(Surv(PFS, Progression) ~ dN.dS.ON.OFF, data = clin.fac.HPV.pos))$coefficients[2],2)
HR.HPV.neg<-round(summary(coxph(Surv(PFS, Progression) ~ dN.dS.ON.OFF, data = clin.fac.HPV.neg))$coefficients[2],2)
fit.HPV.pos <- survfit(Surv(PFS, Progression) ~ dN.dS.ON.OFF, data = clin.fac.HPV.pos)
fit.HPV.neg <- survfit(Surv(PFS, Progression) ~ dN.dS.ON.OFF, data = clin.fac.HPV.neg)
ggsurv.HPV.pos<-ggsurvplot(fit.HPV.pos, data = clin.fac.HPV.pos, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.neg<-ggsurvplot(fit.HPV.neg, data = clin.fac.HPV.neg, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.pos$plot <- ggsurv.HPV.pos$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.pos), size = 5)
ggsurv.HPV.neg$plot <- ggsurv.HPV.neg$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.neg), size = 5)
ggsurv.HPV.pos;ggsurv.HPV.neg

#### HLA A diversity ####
HR.HPV.pos<-round(summary(coxph(Surv(PFS, Progression) ~ normal.HLA.A, data = clin.fac.HPV.pos))$coefficients[2],2)
HR.HPV.neg<-round(summary(coxph(Surv(PFS, Progression) ~ normal.HLA.A, data = clin.fac.HPV.neg))$coefficients[2],2)
fit.HPV.pos <- survfit(Surv(PFS, Progression) ~ normal.HLA.A, data = clin.fac.HPV.pos)
fit.HPV.neg <- survfit(Surv(PFS, Progression) ~ normal.HLA.A, data = clin.fac.HPV.neg)
ggsurv.HPV.pos<-ggsurvplot(fit.HPV.pos, data = clin.fac.HPV.pos, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.neg<-ggsurvplot(fit.HPV.neg, data = clin.fac.HPV.neg, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.pos$plot <- ggsurv.HPV.pos$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.pos), size = 5)
ggsurv.HPV.neg$plot <- ggsurv.HPV.neg$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.neg), size = 5)
ggsurv.HPV.pos;ggsurv.HPV.neg

#### HLA B diversity ####
HR.HPV.pos<-round(summary(coxph(Surv(PFS, Progression) ~ normal.HLA.B, data = clin.fac.HPV.pos))$coefficients[2],2)
HR.HPV.neg<-round(summary(coxph(Surv(PFS, Progression) ~ normal.HLA.B, data = clin.fac.HPV.neg))$coefficients[2],2)
fit.HPV.pos <- survfit(Surv(PFS, Progression) ~ normal.HLA.B, data = clin.fac.HPV.pos)
fit.HPV.neg <- survfit(Surv(PFS, Progression) ~ normal.HLA.B, data = clin.fac.HPV.neg)
ggsurv.HPV.pos<-ggsurvplot(fit.HPV.pos, data = clin.fac.HPV.pos, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.neg<-ggsurvplot(fit.HPV.neg, data = clin.fac.HPV.neg, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.pos$plot <- ggsurv.HPV.pos$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.pos), size = 5)
ggsurv.HPV.neg$plot <- ggsurv.HPV.neg$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.neg), size = 5)
ggsurv.HPV.pos;ggsurv.HPV.neg

#### HLA C diversity ####
HR.HPV.pos<-round(summary(coxph(Surv(PFS, Progression) ~ normal.HLA.C, data = clin.fac.HPV.pos))$coefficients[2],2)
HR.HPV.neg<-round(summary(coxph(Surv(PFS, Progression) ~ normal.HLA.C, data = clin.fac.HPV.neg))$coefficients[2],2)
fit.HPV.pos <- survfit(Surv(PFS, Progression) ~ normal.HLA.C, data = clin.fac.HPV.pos)
fit.HPV.neg <- survfit(Surv(PFS, Progression) ~ normal.HLA.C, data = clin.fac.HPV.neg)
ggsurv.HPV.pos<-ggsurvplot(fit.HPV.pos, data = clin.fac.HPV.pos, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.neg<-ggsurvplot(fit.HPV.neg, data = clin.fac.HPV.neg, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.pos$plot <- ggsurv.HPV.pos$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.pos), size = 5)
ggsurv.HPV.neg$plot <- ggsurv.HPV.neg$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.neg), size = 5)
ggsurv.HPV.pos;ggsurv.HPV.neg

#### purity ####
HR.HPV.pos<-round(summary(coxph(Surv(PFS, Progression) ~ purity, data = clin.fac.HPV.pos))$coefficients[2],2)
HR.HPV.neg<-round(summary(coxph(Surv(PFS, Progression) ~ purity, data = clin.fac.HPV.neg))$coefficients[2],2)
fit.HPV.pos <- survfit(Surv(PFS, Progression) ~ purity, data = clin.fac.HPV.pos)
fit.HPV.neg <- survfit(Surv(PFS, Progression) ~ purity, data = clin.fac.HPV.neg)
ggsurv.HPV.pos<-ggsurvplot(fit.HPV.pos, data = clin.fac.HPV.pos, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.neg<-ggsurvplot(fit.HPV.neg, data = clin.fac.HPV.neg, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.pos$plot <- ggsurv.HPV.pos$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.pos), size = 5)
ggsurv.HPV.neg$plot <- ggsurv.HPV.neg$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.neg), size = 5)
ggsurv.HPV.pos;ggsurv.HPV.neg

#### ploidy ####
HR.HPV.pos<-round(summary(coxph(Surv(PFS, Progression) ~ ploidy, data = clin.fac.HPV.pos))$coefficients[2],2)
HR.HPV.neg<-round(summary(coxph(Surv(PFS, Progression) ~ ploidy, data = clin.fac.HPV.neg))$coefficients[2],2)
fit.HPV.pos <- survfit(Surv(PFS, Progression) ~ ploidy, data = clin.fac.HPV.pos)
fit.HPV.neg <- survfit(Surv(PFS, Progression) ~ ploidy, data = clin.fac.HPV.neg)
ggsurv.HPV.pos<-ggsurvplot(fit.HPV.pos, data = clin.fac.HPV.pos, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.neg<-ggsurvplot(fit.HPV.neg, data = clin.fac.HPV.neg, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.pos$plot <- ggsurv.HPV.pos$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.pos), size = 5)
ggsurv.HPV.neg$plot <- ggsurv.HPV.neg$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.neg), size = 5)
ggsurv.HPV.pos;ggsurv.HPV.neg

#### ITH ####
HR.HPV.pos<-round(summary(coxph(Surv(PFS, Progression) ~ ITH, data = clin.fac.HPV.pos))$coefficients[2],2)
HR.HPV.neg<-round(summary(coxph(Surv(PFS, Progression) ~ ITH, data = clin.fac.HPV.neg))$coefficients[2],2)
fit.HPV.pos <- survfit(Surv(PFS, Progression) ~ ITH, data = clin.fac.HPV.pos)
fit.HPV.neg <- survfit(Surv(PFS, Progression) ~ ITH, data = clin.fac.HPV.neg)
ggsurv.HPV.pos<-ggsurvplot(fit.HPV.pos, data = clin.fac.HPV.pos, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.neg<-ggsurvplot(fit.HPV.neg, data = clin.fac.HPV.neg, risk.table = TRUE,pval = T, conf.int = T, ggtheme = theme_bw(), pval.coord = c(750,0.75), ylab = 'PFS', xlab = 'Time (days)')
ggsurv.HPV.pos$plot <- ggsurv.HPV.pos$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.pos), size = 5)
ggsurv.HPV.neg$plot <- ggsurv.HPV.neg$plot + ggplot2::annotate("text", x = 930, y = 0.9, vjust = 1, hjust = 1, label = paste0("HR=",HR.HPV.neg), size = 5)
ggsurv.HPV.pos;ggsurv.HPV.neg


#### ####
clin.fac.HPV<-clin.fac
df<-clin.fac.HPV
df$PDL1[df$PDL1==1]<-'amp'
df$PDL1[df$PDL1==-1]<-'del'
df$PDL1[df$PDL1==0]<-'wt'
df$del9p<-as.factor(ifelse(df$del9p==1,'del','wt'))
df$PRIMARY_TUMOR_SITE<-as.factor(ifelse(df$PRIMARY_TUMOR_SITE=="Oral Cavity", "Oral Cavity", "Non-oral cavity"))
df$APOBEC_Signature<-as.factor(ifelse(df$APOBEC_Signature=='yes','Present','Absent'))
df$Smoking_Signature<-as.factor(ifelse(df$Smoking_Signature=='yes','Present','Absent'))
df$AUTOIMMUNE_DISEASES<-as.factor(ifelse(df$AUTOIMMUNE_DISEASES=='Yes','Present','Absent'))
df$ALLERGIES<-as.factor(ifelse(df$ALLERGIES=='Yes','Present','Absent'))

colnames(df)[colnames(df)=='TERT.promoter']<-'TERT promoter'
colnames(df)[colnames(df)=='Smoking_Signature']<-'Smoking signature'
colnames(df)[colnames(df)=='APOBEC_Signature']<-'APOBEC signature'
colnames(df)[colnames(df)=='Smoking_Signature']<-'Smoking signature'
colnames(df)[colnames(df)=='nonsyn.TMB']<-'nonsyn TMB'
colnames(df)[colnames(df)=='Indel_load']<-'Indel load'
colnames(df)[colnames(df)=='Clonal_Mutational_Load']<-'Clonal mutational load'
colnames(df)[colnames(df)=='dN.dS.ON.OFF']<-'dN/dS ON/OFF'
colnames(df)[colnames(df)=='purity']<-'Purity'
colnames(df)[colnames(df)=='ploidy']<-'Ploidy'
colnames(df)[colnames(df)=='SEX']<-'Sex'
colnames(df)[colnames(df)=='Age_at_IO_start']<-'Age at Immunotherapy start'
colnames(df)[colnames(df)=='ALCOHOL_NEVER_EVER']<-'Alcohol use'
colnames(df)[colnames(df)=='AUTOIMMUNE_DISEASES']<-'Autoimmune disease'
colnames(df)[colnames(df)=='ALLERGIES']<-'Allergy'
colnames(df)[colnames(df)=='METASTATIC_VS_RECURRENT']<-'Type of tumor'
colnames(df)[colnames(df)=='PRIMARY_TUMOR_SITE']<-'Tumor site'
colnames(df)[colnames(df)=='Neut.Mono.Lymp']<-'SIRI'
colnames(df)[colnames(df)=='INFECTION_DURING_IO']<-'Infection'
colnames(df)[colnames(df)=='STEROIDS']<-'Steroids'
colnames(df)[colnames(df)=='ANTIBIOTICS']<-'Antibiotics'
colnames(df)[colnames(df)=='HLALOH']<-'HLA LOH'
colnames(df)[colnames(df)=='AnyVirus']<-'Viral status'
colnames(df)[colnames(df)=='del9q34.3']<-'9q34.3'
colnames(df)[colnames(df)=='del9p']<-'9p'

df$ECOG<-(ifelse(df$ECOG==0,"0",'1-2'))


map(vars(TP53,PIK3CA,'TERT promoter',PDL1,'9q34.3',MET,'HLA LOH','Smoking signature','APOBEC signature',
         'nonsyn TMB','Indel load','Clonal mutational load','dN/dS ON/OFF',
         HED,'Purity','Ploidy',ITH,
         Sex,'Age at Immunotherapy start',ECOG,BMI,'Alcohol use',
         'Autoimmune disease','Allergy','Viral status',
         'Type of tumor','Tumor site','SIRI',Platelets,
         HGB,Albumin,'Steroids','Infection','Antibiotics'), function(by)
         {
           analyse_multivariate(df,
                                reference_level_dict=c(TP53='wt',PIK3CA='wt','TERT promoter'='wt',
                                                       '9q34.3'='wt', MET='wt',PDL1='wt','HLA LOH'='No LOH',ECOG="0",
                                                       'APOBEC signature'='Absent','Smoking signature'='Absent','Type of tumor'='Non-metastatic',
                                                       'Autoimmune disease'='Absent','Alcohol use'='No','Allergy'='Absent','Viral status'='Negative',
                                                       'Antibiotics'='No','Steroids'='No','Infection'='No','Tumor site'='Non-oral cavity'
                                ),
                                vars(OS,Overall_Survival_Event),
                                covariates = list(by))
         }) %>% forest_plot



df<-df[df$Patient.ID%in%c(read.csv('~/PATH2FILE/H&N/test.csv')$x,read.csv('~/PATH2FILE/H&N/training.csv')$x),]
df %>% analyse_multivariate(
  reference_level_dict=c(TP53='wt',PIK3CA='wt','TERT promoter'='wt','9p'='wt',
                         '9q34.3'='wt', MET='wt',PDL1='wt','HLA LOH'='No LOH',ECOG="0",
                         'APOBEC signature'='Absent','Smoking signature'='Absent','Type of tumor'='Non-metastatic',
                         'Autoimmune disease'='Absent','Alcohol use'='No','Allergy'='Absent','Viral status'='Negative',
                         'Antibiotics'='No','Steroids'='No','Infection'='No','Tumor site'='Non-oral cavity'
  ),
  vars(PFS, Progression),
  covariates = vars(TP53,PIK3CA,PDL1,'9q34.3',MET,'Smoking signature',
                    'nonsyn TMB',"9p",
                    ECOG,'Viral status','Autoimmune disease',Platelets,
                    'Tumor site','SIRI')) %>% forest_plot

df %>% analyse_multivariate(
  reference_level_dict=c(TP53='wt',PIK3CA='wt','TERT promoter'='wt',
                         del9q34.3='wt', MET='wt',PDL1='wt','HLA LOH'='No LOH',ECOG="0",
                         'APOBEC signature'='Absent','Smoking signature'='Absent','Type of tumor'='Non-metastatic',
                         'Autoimmune disease'='Absent','Alcohol use'='No','Allergy'='Absent','Viral status'='Negative',
                         'Antibiotics'='No','Steroids'='No','Infection'='No','Tumor site'='Non-oral cavity'
  ),
  vars(OS, Overall_Survival_Event),
  covariates = vars(TP53,PIK3CA,'TERT promoter',PDL1,del9q34.3,MET,'HLA LOH','Smoking signature','APOBEC signature',
                    'nonsyn TMB','Indel load','Clonal mutational load','dN/dS ON/OFF',
                    HED,'Purity','Ploidy',ITH,
                    Sex,'Age at Immunotherapy start',ECOG,BMI,'Alcohol use',
                    'Autoimmune disease','Allergy','Viral status',
                    'Type of tumor','Tumor site','SIRI',Platelets,
                    HGB,Albumin,'Steroids','Infection','Antibiotics')) %>% forest_plot




clin.fac<-WES

clin.fac$nonsyn.TMB<-(ifelse(clin.fac$nonsyn.TMB>3.34,'high','low'))
clin.fac$Indel_load<-(ifelse(clin.fac$Indel_load>5,'high','low'))
clin.fac$Clonal_Mutational_Load<-(ifelse(clin.fac$Clonal_Mutational_Load>103.7,'high','low'))
clin.fac$dN.dS.ON.OFF<-(ifelse(clin.fac$dN.dS.ON.OFF>1.6,'high','low'))
clin.fac$tumor.HLA.A<-(ifelse(clin.fac$tumor.HLA.A>15,'high','low'))
clin.fac$tumor.HLA.B<-(ifelse(clin.fac$tumor.HLA.B>15,'high','low'))
clin.fac$tumor.HLA.C<-(ifelse(clin.fac$tumor.HLA.C>15,'high','low'))
clin.fac$ploidy<-(ifelse(clin.fac$ploidy>2.6,'high','low'))
clin.fac$purity<-(ifelse(clin.fac$purity>0.54,'high','low'))
clin.fac$ITH<-(ifelse(clin.fac$ITH>25,'high','low'))

M<-clin.fac[,c("TP53","Smoking_Signature","PIK3CA","ITH","PDL1.amp","PDL1.del","APOBEC_Signature",
               "del9q34.3","Clonal_Mutational_Load","purity",
               "ploidy","nonsyn.TMB","MET.amp","AnyVirus")]
M<-clin.fac[,c("TP53","Smoking_Signature","ITH","del9p","APOBEC_Signature",
               "purity",
               "nonsyn.TMB","AnyVirus")]
M[M=='high']=1;M[M=='low']=0;M[M=="mut"]=1;M[M=='wt']=0;M[M=='yes']=1;M[M=='no']=0;M[M=='Positive']=1;M[M=='Negative']=0
mat=t(data.matrix(M))
mat[c(4,8),]<-mat[c(4,8),]+1

colnames(mat)<-clin.fac$sample
my_palette <- colorRampPalette(c('red','blue'))(n = 1000)
mat<-mat-1
H<-heatmap.2(mat,trace='none', cexRow = 0.4, cexCol = 0.5, col = my_palette)
C<-WES[H$colInd,c('BEST_RESPONSE','PRIMARY_TUMOR_SITE')]
cols<-RColorBrewer::brewer.pal(n = 4, name = "RdYlGn")
C$BEST_RESPONSE <- factor(C$BEST_RESPONSE , levels=c("CR", "PR", "SD", "PD"))

row_ha = columnAnnotation(Response = C$BEST_RESPONSE, 'Tumor site' = C$PRIMARY_TUMOR_SITE,
                       col = list( Response = c("CR" = cols[4], "PR" = cols[3], "SD" = cols[2], "PD" = cols[1]),
                                   'Tumor site' = c("Hypopharynx"=1, "Larynx"=2, "Nasopharynx"=3, "Oral Cavity"=4, "Oropharynx"=5, "Sinonasal"=6, "Multiple sites"=7, "Unknown"=8)))
Heatmap((mat)[,H$colInd],cluster_rows = T,cluster_columns = F, column_names_gp = gpar(fontsize = 5), name = "mat",bottom_annotation =  row_ha)



A<-cutree(hclust(dist(t(mat))), k=6)
for(i in 1:nrow(clin.fac)){if(sum(names(A)%in%clin.fac$sample[i])==0){clin.fac$subtype[i]<-NA; next}
  clin.fac$subtype[i]<-A[names(A)==clin.fac$sample[i]]}
clin.fac$subtype[clin.fac$subtype==6]<-'subtype1'
clin.fac$subtype[clin.fac$subtype==3]<-'subtype2'
clin.fac$subtype[clin.fac$subtype==4]<-'subtype3'
clin.fac$subtype[clin.fac$subtype==2]<-'subtype4'
clin.fac$subtype[clin.fac$subtype==1]<-'subtype5'
clin.fac$subtype[clin.fac$subtype==5]<-'subtype6'

clin.fac$subtype[clin.fac$TP53=="mut" & clin.fac$nonsyn.TMB == "low"]<-'subtype1'
clin.fac$subtype[clin.fac$TP53=="mut" & clin.fac$nonsyn.TMB == "high" & (clin.fac$PDL1.del==1)]<-'subtype2'
clin.fac$subtype[clin.fac$TP53=="mut" & clin.fac$nonsyn.TMB == "high" & !(clin.fac$PDL1.del==1)]<-'subtype3'
clin.fac$subtype[clin.fac$TP53=="wt" & clin.fac$nonsyn.TMB == "high"]<-'subtype4'
clin.fac$subtype[clin.fac$TP53=="wt" & clin.fac$nonsyn.TMB == "low" & clin.fac$purity == "low"]<-'subtype5'
clin.fac$subtype[clin.fac$TP53=="wt" & clin.fac$nonsyn.TMB == "low" & clin.fac$purity == "high"]<-'subtype6'

fit <- rpart(subtype ~ TP53 + nonsyn.TMB + PDL1.del + purity , data = clin.fac)
plot(fit)
text(fit, use.n = TRUE)
summary(fit)


clin.fac$subtype<-as.factor(clin.fac$subtype)
clin.fac$OS<-clin.fac$OS*30; clin.fac$PFS<-clin.fac$PFS*30; clin.fac$DSS<-clin.fac$DSS*30
clin.fac %>% analyse_survival(vars(PFS,Progression), clin.fac$subtype,cox_reference_level='subtype1', p_adjust_method="BH") -> OS.result
clin.fac %>% analyse_survival(vars(OS,Overall_Survival_Event), clin.fac$subtype,cox_reference_level='subtype1', p_adjust_method="BH") -> OS.result
#clin.fac %>% analyse_survival(vars(DSS,Disease_Specific_Survival_Event), clin.fac$subtype,cox_reference_level='subtype1', p_adjust_method="BH") -> OS.result
default_args <- list(xlab=c("PFS (year)"),
                     break.time.by="breakByYear",
                     hazard.ratio=T,
                     risk.table=TRUE,
                     table.layout="clean", 
                     tables.height = 0.15,xscale = "d_y",
                     ggtheme=ggplot2::theme_bw(10))

K<-kaplan_meier_plot(OS.result, default_args,hazard.position = "top",
                     mapped_plot_args=list(legend.title=c('subtype'),
                                           legend=c('right','right')
                     ))
K$plot$layers[[4]]$data$y<-0.85
K$plot$layers[[4]]$data$x<-600
K


#### risk ###
clin.fac$subtype<-ifelse(clin.fac$Clonal_Mutational_Load=='high', 'Low risk','High risk')
clin.fac$subtype<-ifelse(clin.fac$PDL1.del==0 , 'Low risk','High risk')

clin.fac$subtype<-as.factor(clin.fac$subtype)
#clin.fac$OS<-clin.fac$OS*30; clin.fac$PFS<-clin.fac$PFS*30; clin.fac$DSS<-clin.fac$DSS*30
df<-clin.fac[clin.fac$AnyVirus==0,]
df %>% analyse_survival(vars(PFS,Progression), df$subtype,cox_reference_level='subtype1', p_adjust_method="BH") -> OS.result
df %>% analyse_survival(vars(OS,Overall_Survival_Event), df$subtype,cox_reference_level='subtype1', p_adjust_method="BH") -> OS.result
df %>% analyse_survival(vars(DSS,Disease_specific_survival_event), df$subtype,cox_reference_level='subtype1', p_adjust_method="BH") -> OS.result
default_args <- list(xlab=c("PFS (months)"),
                     break.time.by="breakByYear",
                     hazard.ratio=T,
                     risk.table=TRUE,
                     table.layout="clean", 
                     ggtheme=ggplot2::theme_bw(10))

K<-kaplan_meier_plot(OS.result, default_args,hazard.position = "top",
                     mapped_plot_args=list(legend.title=c('subtype'),
                                           legend=c('right','right')
                     ))
K$plot$layers[[4]]$data$y<-0.9
K$plot$layers[[4]]$data$x<-400
K



par(mfrow=c(2,10))
barplot(aggregate(clin.fac$PDL1.del, list(clin.fac$subtype), FUN=mean)$x,
        names=aggregate(clin.fac$PDL1.del, list(clin.fac$subtype), FUN=mean)$Group.1,
        ylab='Proportion',ylim=c(0,1),main='PD-L1 del',las =2,
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$PDL1.amp, list(clin.fac$subtype), FUN=mean)$x,
        names=aggregate(clin.fac$PDL1.amp, list(clin.fac$subtype), FUN=mean)$Group.1,
        ylab='Proportion',ylim=c(0,1),main='PD-L1 amp',las =2,
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$del9p, list(clin.fac$subtype), FUN=mean)$x,
        names=aggregate(clin.fac$del9p, list(clin.fac$subtype), FUN=mean)$Group.1,
        ylab='Proportion',ylim=c(0,1),main='9p del',las =2,
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$CDKN2A, list(clin.fac$subtype), FUN=mean)$x,
        names=aggregate(clin.fac$CDKN2A, list(clin.fac$subtype), FUN=mean)$Group.1,
        ylab='Proportion',ylim=c(0,1),main='CDKN2A del',las =2,
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$del9q34.3, list(clin.fac$subtype), FUN=mean)$x,
        names=aggregate(clin.fac$del9q34.3, list(clin.fac$subtype), FUN=mean)$Group.1,
        ylab='Proportion',ylim=c(0,1),main='9q34.3 del',las =2,
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$PTEN.del, list(clin.fac$subtype), FUN=mean)$x,
        names=aggregate(clin.fac$PTEN.del, list(clin.fac$subtype), FUN=mean)$Group.1,
        ylab='Proportion',ylim=c(0,1),main='PTEN del',las =2,
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$MET.amp, list(clin.fac$subtype), FUN=mean)$x,
        names=aggregate(clin.fac$MET.amp, list(clin.fac$subtype), FUN=mean)$Group.1,
        ylab='Proportion',ylim=c(0,1),main='MET amp',las =2,
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$CCND1.amp, list(clin.fac$subtype), FUN=mean)$x,
        names=aggregate(clin.fac$CCND1.amp, list(clin.fac$subtype), FUN=mean)$Group.1,
        ylab='Proportion',ylim=c(0,1),main='CCND1 amp',las =2,
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$TERT.promoter, list(clin.fac$subtype), FUN=mean)$x,
        names=aggregate(clin.fac$TERT.promoter, list(clin.fac$subtype), FUN=mean)$Group.1,
        ylab='Proportion',ylim=c(0,1),main='TERT promoter mut',las =2,
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$TP53=="mut", list(clin.fac$subtype), FUN=mean)$x,
        names=aggregate(clin.fac$TP53, list(clin.fac$subtype), FUN=mean)$Group.1,
        ylab='Proportion',ylim=c(0,1),main='TP53',las =2,
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$PIK3CA=="mut", list(clin.fac$subtype), FUN=mean)$x,
        names=aggregate(clin.fac$PIK3CA, list(clin.fac$subtype), FUN=mean)$Group.1,
        ylab='Proportion',ylim=c(0,1),main='PIK3CA',las =2,
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$Smoking_Signature=='yes', list(clin.fac$subtype), FUN=mean)$x,
        names=aggregate(clin.fac$Smoking_Signature=='yes', list(clin.fac$subtype), FUN=mean)$Group.1,
        ylab='Proportion',ylim=c(0,1),main='Smoking signature',las =2,
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$APOBEC_Signature=='yes', list(clin.fac$subtype), FUN=mean)$x,
        names=aggregate(clin.fac$Smoking_Signature=='yes', list(clin.fac$subtype), FUN=mean)$Group.1,
        ylab='Proportion',ylim=c(0,1),main='APOBEC signature',las =2,
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$nonsyn.TMB=='high', list(clin.fac$subtype), FUN=mean)$x,
        names=aggregate(clin.fac$nonsyn.TMB=='high', list(clin.fac$subtype), FUN=mean)$Group.1,
        ylab='Proportion',main='High nonsyn TMB',las =2,ylim=c(0,1),
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$Clonal_Mutational_Load=='high', list(clin.fac$subtype), FUN=function(y) mean(y,na.rm = T))$x,
        names=aggregate(clin.fac$Clonal_Mutational_Load=='high', list(clin.fac$subtype), FUN=function(y) mean(y,na.rm = T))$Group.1,
        ylab='Proportion',main='High clonal TMB',las =2,ylim=c(0,1),
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$purity=='high', list(clin.fac$subtype), FUN=mean)$x,
        names=aggregate(clin.fac$purity=='high', list(clin.fac$subtype), FUN=mean)$Group.1,
        ylab='Proportion',main='High purity',las =2,ylim=c(0,1),
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$ploidy=='high', list(clin.fac$subtype), FUN=mean)$x,
        names=aggregate(clin.fac$ploidy=='high', list(clin.fac$subtype), FUN=mean)$Group.1,
        ylab='Proportion',main='High ploidy',las =2,ylim=c(0,1),
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$ITH=='high', list(clin.fac$subtype), FUN=function(y) mean(y,na.rm = T))$x,
        names=aggregate(clin.fac$ITH=='high', list(clin.fac$subtype), FUN=function(y) mean(y,na.rm = T))$Group.1,
        ylab='Proportion',main='High ITH',las =2,ylim=c(0,1),
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$AnyVirus, list(clin.fac$subtype), FUN=mean)$x,
        names=aggregate(clin.fac$AnyVirus, list(clin.fac$subtype), FUN=mean)$Group.1,
        ylab='Proportion',ylim=c(0,1),main='V+',las =2,
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
barplot(aggregate(clin.fac$Neut.Mono.Lymp, list(clin.fac$subtype), FUN=function(y) mean(y,na.rm = T))$x,
        names=aggregate(clin.fac$Neut.Mono.Lymp, list(clin.fac$subtype), FUN=function(y) mean(y,na.rm = T))$Group.1,
        ylab='Mean',main='SIRI',las =2,
        col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
# barplot(aggregate(clin.fac$Clinical_benefit.CR.PR.SD.6m., list(clin.fac$subtype), FUN=function(y) mean(y,na.rm = T))$x,
#         names=aggregate(clin.fac$Clinical_benefit.CR.PR.SD.6m., list(clin.fac$subtype), FUN=function(y) mean(y,na.rm = T))$Group.1,
#         ylab='Proportion',main='Clinical benefit',las =2,ylim=c(0,1),
#         col=c('red','darkgoldenrod3','green','aquamarine2','deepskyblue2','magenta'))
clin.fac1<-clin.fac
clin.fac1$subtype<-factor(clin.fac1$subtype,
                          levels = c("subtype4",
                                     "subtype3",
                                     "subtype5",
                                     "subtype1",
                                     "subtype6",
                                     "subtype2"))
par(mfrow=c(3,3))
for(i in 129:136){
  boxplot(clin.fac1[,i]~clin.fac1$subtype,col=2:7,outline=F,las =2,xlab='',ylab=colnames(clin.fac1)[i] )
  stripchart(clin.fac1[,i]~clin.fac1$subtype,
             vertical = TRUE, method = "jitter",
             pch = 21, col = "blue", bg = "bisque",
             add = TRUE, cex = 0.7)
  
}
WES$BEST_RESPONSE <- factor(WES$BEST_RESPONSE , levels=c("CR", "PR", "SD", "PD"))



####
WES<-read.csv('~/PATH2FILE/H&N/WES-HNSC-features.csv')
clin.fac<-WES[!is.na(WES$RiskGroup),]
clin.fac$OS<-clin.fac$OS*30; clin.fac$PFS<-clin.fac$PFS*30; clin.fac$DSS<-clin.fac$DSS*30
clin.fac$RiskGroup[clin.fac$RiskGroup==1]<-"High"
clin.fac$RiskGroup[clin.fac$RiskGroup==2]<-"Intermediate"
clin.fac$RiskGroup[clin.fac$RiskGroup==3]<-"Low"
clin.fac$RiskGroup<-as.factor(clin.fac$RiskGroup)

clin.fac %>% analyse_survival(vars(PFS,Progression), clin.fac$RiskGroup,cox_reference_level='Low', p_adjust_method="BH") -> OS.result
clin.fac %>% analyse_survival(vars(OS,Overall_Survival_Event), clin.fac$RiskGroup,cox_reference_level='Low', p_adjust_method="BH") -> OS.result

default_args <- list(xlab=c("OS (year)"),
                     break.time.by="breakByYear",
                     hazard.ratio=T,
                     risk.table=TRUE,
                     table.layout="clean", xlim = c(0, 365*4),
                     tables.height = 0.15,xscale = "d_y",
                     ggtheme=ggplot2::theme_bw(10))

K<-kaplan_meier_plot(OS.result, default_args,hazard.position = "top",
                     mapped_plot_args=list(legend.title=c('Risk Group'),
                                           legend=c('right','right')
                     ))
K$plot$layers[[4]]$data$y<-0.85
K$plot$layers[[4]]$data$x<-600
K

IMPACT<-read.csv('~/PATH2FILE/H&N/IMPACT-validation.csv')
IMPACT$RiskGroup<-"Low"
IMPACT$RiskGroup[(IMPACT$nonsyn.TMB<2 & IMPACT$Neut.Mono.Lymp <13) | (  IMPACT$Neut.Mono.Lymp <13 & IMPACT$Smoking_Signature=="Present")]<-"Intermediate"
IMPACT$RiskGroup[IMPACT$Neut.Mono.Lymp>13]<-"High"
IMPACT$RiskGroup[IMPACT$Neut.Mono.Lymp<13 & (IMPACT$Smoking_Signature!="Present" & IMPACT$nonsyn.TMB>2)]<-"Low"

clin.fac<-IMPACT
clin.fac$OS<-clin.fac$OS*30; clin.fac$PFS<-clin.fac$PFS*30; 
clin.fac$RiskGroup<-as.factor(clin.fac$RiskGroup)
clin.fac %>% analyse_survival(vars(PFS,Progression), clin.fac$RiskGroup,cox_reference_level='Low', p_adjust_method="BH") -> OS.result
clin.fac %>% analyse_survival(vars(OS,Overall_Survival_Event), clin.fac$RiskGroup,cox_reference_level='Low', p_adjust_method="BH") -> OS.result

default_args <- list(xlab=c("PFS (year)"),
                     break.time.by="breakByYear",
                     hazard.ratio=T,
                     risk.table=TRUE,
                     table.layout="clean", xlim = c(0, 365*4),
                     tables.height = 0.15,xscale = "d_y",
                     ggtheme=ggplot2::theme_bw(10))

K<-kaplan_meier_plot(OS.result, default_args,hazard.position = "top",
                     mapped_plot_args=list(legend.title=c('Risk Group'),
                                           legend=c('right','right')
                     ))
K$plot$layers[[4]]$data$y<-0.85
K$plot$layers[[4]]$data$x<-500
K
####


                  
