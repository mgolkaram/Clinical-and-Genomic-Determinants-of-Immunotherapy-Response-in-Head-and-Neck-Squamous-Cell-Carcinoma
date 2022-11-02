require(gplots)
WES<-read.csv('~/PATH2FILE/H&N/WES-HNSC-features.csv')
arm<-read.csv('~/PATH2FILE/H&N/Merck WES samples pipeline output/arm_result_for_Luc.csv')
del.9p<-arm$sample[arm$arm=='9p' & arm$arm_mean_tcn<2]

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
M[M=='high']=1;M[M=='low']=0;M[M=="mut"]=1;M[M=='wt']=0;M[M=='yes']=1;M[M=='no']=0;M[M=='Positive']=1;M[M=='Negative']=0
mat=t(data.matrix(M))
mat[c(5,6,8,13,14),]<-mat[c(5,6,8,13,14),]+1
colnames(mat)<-clin.fac$sample
my_palette <- colorRampPalette(c('red','blue'))(n = 1000)
mat<-mat-1
H<-heatmap.2(mat,trace='none', cexRow = 0.4, cexCol = 0.5, col = my_palette)
C<-WES[H$colInd,c('BEST_RESPONSE','PRIMARY_TUMOR_SITE')]
cols<-RColorBrewer::brewer.pal(n = 4, name = "RdYlGn")
C$BEST_RESPONSE <- factor(C$BEST_RESPONSE , levels=c("CR", "PR", "SD", "PD"))


A<-cutree(hclust(dist(t(mat))), k=6)
for(i in 1:nrow(clin.fac)){if(sum(names(A)%in%clin.fac$sample[i])==0){clin.fac$subtype[i]<-NA; next}
  clin.fac$subtype[i]<-A[names(A)==clin.fac$sample[i]]}
#clin.fac$subtype[clin.fac$subtype==1]<-'subtype1'
clin.fac$subtype[clin.fac$subtype==6]<-'subtype1'
clin.fac$subtype[clin.fac$subtype==3]<-'subtype2'
clin.fac$subtype[clin.fac$subtype==4]<-'subtype3'
clin.fac$subtype[clin.fac$subtype==2]<-'subtype4'
clin.fac$subtype[clin.fac$subtype==1]<-'subtype5'
clin.fac$subtype[clin.fac$subtype==5]<-'subtype6'

clin.fac1<-clin.fac[,colnames(clin.fac)%in%c("purity","ITH","PDL1.del","del9q34.3","AnyVirus","MET.amp","PD1.amp","PIK3CA",
                                             "Clonal_Mutational_Load","nonsyn.TMB","ploidy","TP53","Smoking_Signature",
                                             "APOBEC_Signature","subtype")]

fit <- rpart::rpart(subtype ~ purity + ITH + clin.fac1$PDL1.del + clin.fac1$del9q34.3  + clin.fac1$AnyVirus + 
                      clin.fac1$MET.amp + clin.fac1$PD1.amp +clin.fac1$PIK3CA + clin.fac1$Clonal_Mutational_Load + 
                      clin.fac1$nonsyn.TMB + clin.fac1$ploidy + clin.fac1$TP53 + clin.fac1$Smoking_Signature + 
                      clin.fac1$APOBEC_Signature , data = clin.fac1,cost = rep(0.05,14))
fit <- randomForest::randomForest(as.factor(subtype) ~ ., data = clin.fac1[,-ncol(clin.fac1)], importance=TRUE, na.action=na.omit)


fit <- rpart::rpart(subtype ~ purity + PDL1.del + AnyVirus + 
                      nonsyn.TMB + TP53, data = clin.fac1,cost = rep(0.05,5))

fit <- rpart::rpart(subtype ~ purity + PDL1.del +  
                      nonsyn.TMB + TP53, data = clin.fac1,cost = rep(0.05,4))

fit <- rpart::rpart(subtype ~ purity + PDL1.del +  PD1.amp +
                      nonsyn.TMB + TP53, data = clin.fac1,cost = rep(0.05,5))


plot(fit); text(fit, use.n = TRUE)

for(TMB.cutoff in 85:120){
  TMB.cutoff<-90
df<-read.csv('~/Box/H&N/Merck WES samples pipeline output/Summary.csv')
df$del.9p<-ifelse(df$tumorName%in%del.9p,1,0)
maf<-read.delim('~/PATH2FILE/H&N/MERCK-MAF.tsv') # MERCK MAF
maf<-maf[maf$Variant_Classification%in%unique(maf$Variant_Classification)[c(7,13,15)],]

#df$TMB<-as.numeric(table(maf$TUMOR_SAMPLE))
df$HPV_type<-ifelse(df$HPV_type=='negative',0,1)

colnames(df)[c(7,9,14)]<-c('AnyVirus','PDL1.del','nonsyn.TMB')
df$PDL1.amp<-ifelse(df$PDL1.del>0,1,0)
df$PDL1.del<-ifelse(df$PDL1.del<0,1,0)

df$TP53<-ifelse(df$TP53>0,'mut','wt')
df$purity<-ifelse(df$purity>0.54,'high','low')
df$nonsyn.TMB<-ifelse(df$nonsyn.TMB>TMB.cutoff,'high','low')
rownames(df)<-df$tumorName
# df<-df[,colnames(df)%in%c("purity","PDL1.del","PDL1.del",
#                           "nonsyn.TMB","TP53")]
df<-df[,colnames(df)%in%c("purity","PDL1.del","AnyVirus","PDL1.del",
                          "nonsyn.TMB","TP53")]
df<-df[rownames(df)!='P177T',]
df<-df[-which(df$TP53=='wt' & df$nonsyn.TMB =='low' & is.na(df$purity)),]
dg<-read.csv('~/Box/H&N/Merck WES samples pipeline output/Summary.csv')
#dg$TMB<-as.numeric(table(maf$TUMOR_SAMPLE))
dg<-dg[dg$tumorName !='P177T',]
dg<-dg[-which(dg$TP53<1 & dg$TMB < TMB.cutoff+0.001 & is.na(dg$purity)),]

mahdi<-paste0('subtype',apply(predict(fit,clin.fac1),FUN = which.max,MARGIN = 1))
clin.fac1$risk<-ifelse(clin.fac1$subtype=='subtype1' | clin.fac1$subtype=='subtype2' | clin.fac1$subtype=='subtype6','high risk','low risk')
tmp<-table(ifelse(mahdi=='subtype1' | mahdi=='subtype2' | mahdi=='subtype6', 'high risk', 'low risk' ),clin.fac1$risk)

1-(tmp[2,1]+tmp[1,2])/(tmp[1,1]+tmp[2,2])

mahdi<-paste0('subtype',apply(predict(fit,df),FUN = which.max,MARGIN = 1))
df$risk<-ifelse(mahdi=='subtype1' | mahdi=='subtype2' | mahdi=='subtype6',0,1)
tmp<-table(df$risk,dg$Response)
print(TMB.cutoff)
print(fisher.test(tmp))
}

1-(tmp[2,1]+tmp[1,2])/(tmp[1,1]+tmp[2,2])
