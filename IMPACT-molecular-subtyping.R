IMPACT<-read.csv('~/Box/H&N/IMPACT-validation.csv')
clin.fac<-IMPACT
clin.fac$nonsyn.TMB<-ifelse(clin.fac$nonsyn.TMB>3.34,'high','low')

clin.fac$purity<-ifelse(clin.fac$TUMOR_PURITY_IMPACT>0.54,'high','low')

clin.fac$subtype[clin.fac$TP53=="mut" & clin.fac$nonsyn.TMB == "low"]<-'subtype1'
clin.fac$subtype[clin.fac$TP53=="mut" & clin.fac$nonsyn.TMB == "high" & (clin.fac$PDL1.del==1)]<-'subtype2'
clin.fac$subtype[clin.fac$TP53=="mut" & clin.fac$nonsyn.TMB == "high" & !(clin.fac$PDL1.del==1)]<-'subtype3'
clin.fac$subtype[clin.fac$TP53=="wt" & clin.fac$nonsyn.TMB == "high"]<-'subtype4'
clin.fac$subtype[clin.fac$TP53=="wt" & clin.fac$nonsyn.TMB == "low" & clin.fac$purity == "low"]<-'subtype5'
clin.fac$subtype[clin.fac$TP53=="wt" & clin.fac$nonsyn.TMB == "low" & clin.fac$purity == "high"]<-'subtype6'

clin.fac$risk<-NA
clin.fac$risk[clin.fac$subtype=='subtype1' | clin.fac$subtype=='subtype2' | clin.fac$subtype=='subtype6']<-0
clin.fac$risk[clin.fac$subtype=='subtype3' | clin.fac$subtype=='subtype4' | clin.fac$subtype=='subtype5']<-1

table(clin.fac$subtype)
table(clin.fac$risk)

fisher.test(table(clin.fac$risk,clin.fac$Nonresponders_Responders))

clin.fac$OS<-clin.fac$OS*30; clin.fac$PFS<-clin.fac$PFS*30; 
clin.fac$risk<-as.factor(clin.fac$risk)

clin.fac %>% analyse_survival(vars(PFS,Progression), clin.fac$subtype,cox_reference_level='subtype1', p_adjust_method="BH") -> OS.result

default_args <- list(xlab=c("Time (months)"),
                     break.time.by="breakByYear",
                     hazard.ratio=T,
                     risk.table=TRUE,
                     table.layout="clean", xlim = c(0, 365*6),
                     tables.height = 0.15,xscale = "d_m",
                     ggtheme=ggplot2::theme_bw(10))

K<-kaplan_meier_plot(OS.result, default_args,hazard.position = "top",
                     mapped_plot_args=list(legend.title=c('subtype'),
                                           legend=c('right','right')
                     ))
K$plot$layers[[4]]$data$y<-0.85
K$plot$layers[[4]]$data$x<-500
K

