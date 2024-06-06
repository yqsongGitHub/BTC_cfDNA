source("./src/library.R")
source("./src/function.R")
# TP is True Positive fraction，sensitivity；
# FP is False Positive fraction，1-specificity
dat1 <- read.csv("./data/metadata.csv")
dat2 <- read.csv("./data/metadata_clinical.csv")
# cfDNA1
dat1.1 <- dat1[grepl("1",dat1$ID1),]
dat1.1$ID1 <- gsub("1","",dat1.1$ID1)
dat <- merge(dat2,dat1.1,by=c("ID1","cancer"))
dat <- dat[!(is.na(dat$time) & is.na(dat$state)),]
dat$state1 <- dat$state=="dead"

#### fig 4b ####
dat$group <- ifelse(dat$AFV > median(dat$AFV,na.rm = T),"high","low")
table(dat$group)
y <- Surv(dat$time,dat$state=="dead")
mx.km <- survfit(y~group, data = dat, type= "kaplan-meier")
pdf("./images/fig4b_AFV.pdf",width = 10,height = 8,onefile=F)
km_plot(mx.km,dat)
dev.off()

#### fig 4c ####
dat$group <- ifelse(dat$CA199 > median(dat$CA199,na.rm = T),"high","low")
table(dat$group)
y <- Surv(dat$time,dat$state=="dead")
mx.km <- survfit(y~group, data = dat, type= "kaplan-meier")
pdf("./images/fig4c_CA199.pdf",width = 10,height = 8,onefile=F)
km_plot(mx.km,dat)
dev.off()

#### fig 4d ####
dat$group <- ifelse(dat$CEA > median(dat$CEA,na.rm = T),"high","low")
table(dat$group)
y <- Surv(dat$time,dat$state=="dead")
mx.km <- survfit(y~group, data = dat, type= "kaplan-meier")
km_plot(mx.km,dat)
pdf("./images/fig4d_CEA.pdf",width = 10,height = 8,onefile=F)
km_plot(mx.km,dat)
dev.off()

#### fig 4e####
dat.Roc <- survivalROC(Stime=dat$time,  
                       status=dat$state=="dead",      
                       marker = dat$AFV,  
                       predict.time = 1500,
                       span = 0.25*NROW(dat)^(-0.20),
                       method="KM")
dat.Roc1 <- survivalROC(Stime=dat$time,  
                        status=dat$state=="dead",      
                        marker = dat$CA199,     
                        predict.time = 1500,
                        method="KM")
dat.Roc2 <- survivalROC(Stime=dat$time,  
                        status=dat$state=="dead",      
                        marker = dat$CEA,     
                        predict.time = 1500,
                        method="KM")

pdf("./images/fig4e_ROC.pdf",width=10,height = 10)
par(mar=c(6,6,6,6))
plot(dat.Roc$FP*100, dat.Roc$TP*100, 
     type="l",
     col="#b30000",
     xlim=c(0,100),ylim=c(0,100),   
     xlab="False positive rate (%)",
     ylab="True positive rate (%)",
     main="Time dependent ROC",
     lwd = 4,
     font.main=1,
     cex.main=2.5, cex.lab=2.5, cex.axis=2)
abline(0,1,col="gray",lty=2,lwd = 3)
lines(dat.Roc1$FP*100, dat.Roc1$TP*100, lwd = 4,type="l",col="#377eb8",xlim=c(0,100), ylim=c(0,100))
lines(dat.Roc2$FP*100, dat.Roc2$TP*100, lwd = 4,type="l",col="#a65628",xlim=c(0,100), ylim=c(0,100))
#lines(dat.Roc3$FP*100, dat.Roc3$TP*100, lwd = 4,type="l",col="#e41a1c",xlim=c(0,100), ylim=c(0,100))
# legend(55,15,
#        c(paste("AFV =",round(dat.Roc$AUC,3)),
#          paste("CA199 =",round(dat.Roc1$AUC,3)),
#          paste("CEA =",round(dat.Roc2$AUC,3))),
#        x.intersp=1, y.intersp=0.8,
#        lty= 1 ,lwd= 3,col=c("#b30000","#377eb8","#a65628"),
#        bty = "n",
#        seg.len=1,cex=2.0)
dev.off()

dat.Roc.list <- list(dat.Roc,dat.Roc1,dat.Roc2)
df0 <- data.frame()
for(i in 1:length(dat.Roc.list)){
  dat.Roc.s <- dat.Roc.list[[i]]
  Youden <- dat.Roc.s$TP - dat.Roc.s$FP
  cut.value <- dat.Roc.s$cut.values[which.max(Youden)]
  df <- data.frame(
    "sensitivity"=dat.Roc.s$TP[which.max(Youden)],
    "specificity"= 1- dat.Roc.s$FP[which.max(Youden)],
    "AUC"=dat.Roc.s$AUC
  )
  df0 <- rbind(df0,df)
}
rownames(df0) <- c("AFV","CA199","CEA")
write.csv(df0,"./data/sensitivity_specificity.csv")

#### fig 4f  ####
# chol_tcga_pan_can_atlas_2018
cli <- read.table("./data/cbioportal/chol_tcga_pan_can_atlas_2018/data_clinical_patient.txt",sep = "\t",header = T)
mut <- read.table("./data/cbioportal/chol_tcga_pan_can_atlas_2018/data_mutations.txt",sep = "\t",header = T)
mut$vaf_t <- mut$t_alt_count/(mut$t_alt_count+mut$t_ref_count)
mut$vaf_n <- mut$n_alt_count/(mut$n_alt_count+mut$n_ref_count)

df0 <- data.frame()
for (i in unique(mut$Tumor_Sample_Barcode)) {
  dat <- mut[mut$Tumor_Sample_Barcode==i,]
  AFVvalue <- AFV(dat$vaf_t,dat$vaf_n)
  df <- data.frame("sample"=i,
                   "AFV"=AFVvalue)
  df0 <- rbind(df0,df)
}
df0$sample <- gsub("-01.*","",unique(df0$sample))
cli1 <- cli[,c(1,30,31)]
dat <- merge(cli1,df0,by.x = "PATIENT_ID",by.y="sample")
dat$OS_STATUS <- gsub(":.*","",dat$OS_STATUS)
dat$OS_STATUS <- as.numeric(dat$OS_STATUS)
dat <- dat[!(dat$OS_STATUS==" "&is.na(dat$OS_MONTHS)),]
dat$OS_DAYS <- as.numeric(dat$OS_MONTHS)*30
res.cut <- surv_cutpoint(dat, time = "OS_DAYS", event = "OS_STATUS",
                         variables = c("AFV"))
dat$group <- ifelse(dat$AFV >res.cut$cutpoint$cutpoint,"high","low")
y <- Surv(dat$OS_DAYS,dat$OS_STATUS==1)
mx.km <- survfit(y~group, data = dat, type= "kaplan-meier")
p <- ggsurvplot(mx.km, data = dat, palette = c("#e41a1c", "#377eb8"), 
                conf.int = F,
                pval.size = 11,
                xlab = 'Time (Days)', 
                conf.int.style='step', legend.title='', title = NULL,
                legend =  c(0.8, 0.9),legend.labs = c("High", "Low"),
                pval = T, pval.method = F,
                risk.table = T, tables.height = 0.3, 
                fontsize = 6,
                risk.table.title = 'Number at risk', risk.table.y.text = F,
                censor = T, censor.shape = 4, censor.size = 8,
                ggtheme = theme(axis.title = element_text(size = 32, colour = 'black'),
                                plot.title = element_text(size = 24, hjust = 0.5, colour = 'black'),
                                axis.text = element_text(size = 24, colour = 'black'),
                                legend.text = element_text(size = 24, colour = 'black'),
                                legend.key = element_rect(fill='transparent'),
                                legend.background = element_blank(),
                                panel.background=element_rect(fill='transparent'),
                                axis.line = element_line(size = 1)))
pdf("./images/fig4f_TCGA.pdf",onefile = F)
print(p)
dev.off()

#### fig 4g #### 
# Gallbladder Cancer (MSK, 2022)
cli <- read.table("./data/cbioportal/gbc_mskcc_2022/data_clinical_patient.txt",sep = "\t",header = T)
mut <- read.table("./data/cbioportal/gbc_mskcc_2022/data_mutations.txt",sep = "\t",header = T)
mut$vaf_t <- mut$t_alt_count/(mut$t_alt_count+mut$t_ref_count)
mut$vaf_n <- mut$n_alt_count/(mut$n_alt_count+mut$n_ref_count)
unique(mut$Tumor_Sample_Barcode)
df0 <- data.frame()
for (i in unique(mut$Tumor_Sample_Barcode)) {
  dat <- mut[mut$Tumor_Sample_Barcode==i,]
  AFVvalue <- AFV(dat$vaf_t,dat$vaf_n)
  df <- data.frame("sample"=i,
                   "AFV"=AFVvalue)
  df0 <- rbind(df0,df)
}
df0$sample <- gsub("-T.*","",unique(df0$sample))
cli1 <- cli[,c(1,26,27)]
dat <- merge(cli1,df0,by.x = "PATIENT_ID",by.y="sample")
dat$OS_STATUS <- gsub(":.*","",dat$OS_STATUS)
dat$OS_STATUS <- as.numeric(dat$OS_STATUS)
dat <- dat[!(dat$OS_STATUS==" "&is.na(dat$OS_MONTHS)),]
dat$OS_DAYS <- as.numeric(dat$OS_MONTHS)*30
res.cut <- surv_cutpoint(dat, time = "OS_DAYS", event = "OS_STATUS",
                         variables = c("AFV"))
dat$group <- ifelse(dat$AFV >res.cut$cutpoint$cutpoint,"high","low")
y <- Surv(dat$OS_DAYS,dat$OS_STATUS==1)
mx.km <- survfit(y~group, data = dat, type= "kaplan-meier")
p <- ggsurvplot(mx.km, data = dat, palette = c("#e41a1c", "#377eb8"), 
                conf.int = F,
                pval.size = 11,
                xlab = 'Time (Days)', 
                conf.int.style='step', legend.title='', title = NULL,
                legend =  c(0.8, 0.9),legend.labs = c("High", "Low"),
                pval = T, pval.method = F,
                risk.table = T, tables.height = 0.3, 
                fontsize = 6,
                risk.table.title = 'Number at risk', risk.table.y.text = F,
                censor = T, censor.shape = 4, censor.size = 8,
                ggtheme = theme(axis.title = element_text(size = 32, colour = 'black'),
                                plot.title = element_text(size = 24, hjust = 0.5, colour = 'black'),
                                axis.text = element_text(size = 24, colour = 'black'),
                                legend.text = element_text(size = 24, colour = 'black'),
                                legend.key = element_rect(fill='transparent'),
                                legend.background = element_blank(),
                                panel.background=element_rect(fill='transparent'),
                                axis.line = element_line(size = 1)))
pdf("./images/fig4g_MSKCC.pdf",onefile = F)
print(p)
dev.off()