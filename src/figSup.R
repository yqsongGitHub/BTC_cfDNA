source("./src/library.R")
source("./src/function.R")

#### Fig S3a ####
meta <- read.csv("./data/metadata.csv")
all1 <- read.csv("./data/Gallbladder.csv")
all1$cancer <- "Gallbladder"
all2 <- read.csv("./data/Cholangiocarcinoma.csv")
all2$cancer <- "Cholangiocarcinoma"
all1 <- rbind(all1,all2)
all1 <- all1[,c("ID1","Gene_Name","Annotation","cancer")]
colnames(all1) <- c("ID1","GENE","VAR_TYPE_SX","Cancer")
all1 <- all1[!is.na(all1$ID1),]

meta <- read.csv("./data/metadata.csv")
keepsample <- read.csv("./data/keepsamples.csv")
keepsample$cfDNA <- NA
keepsample$AFP <- NA
keepsample$CEA <- NA
keepsample$CA199 <- NA
keepsample$CA125 <- NA
keepsample$AFV <- NA
keepsample <- keepsample[,c(1,3:9,2)]
keepsample1 <- keepsample[!keepsample$ID1 %in% meta$ID1,]
meta <- rbind(meta,keepsample1)

mut <- all1[grepl("2",all1$ID1),]
allsamples <- unique(mut$ID1)
mut <- mut[,1:3]
sample.keep <- keepsample$ID1[grepl("2",keepsample$ID1)]

pdf("./images/figS3a_cfDNA2.pdf",width = 12 ,height = 10)
landscape_cancer_cfDNA2 <- draw_landscape1(mut,sample.keep)
dev.off()
write.csv(landscape_cancer_cfDNA2,"./data/figS3a_cancer_cfDNA2.csv")

#### Fig S3b ####
dt1 <- read.csv("./data/Cholangiocarcinoma.csv")
dt2 <- read.csv("./data/Gallbladder.csv")
dt <- rbind(dt1,dt2)
allsamples <- unique(dt$ID1)
mut <- read.csv("./data/mutation.csv",row.names = 1)
mut <- data.frame(rowSums(mut))
colnames(mut) <- "TMB"
meta <- read.csv("./data/metadata.csv")
meta <- merge(meta,mut,by.x="ID1",by.y="row.names",all.x = T)
meta <- meta[meta$ID1 %in% allsamples,]

dat1 <- reshape2::melt(meta)
dat1$value <- log2(dat1$value+1)
dat.box <- dat1[dat1$variable=="TMB",]
dat.box <- na.omit(dat.box)
dat.box$type1 <- NA
dat.box$type1[grepl("1",dat.box$ID1)] <- "cfDNA1"
dat.box$type1[grepl("2",dat.box$ID1)] <- "cfDNA2"
dat.box <- dat.box[!is.na(dat.box$type1),]

p <- draw_jitter1(dat.box,ytitle="log2(TMB+1)",ymax = 6)
pdf("./images/figS3b_TMB_cfDNA1_cfDNA2.pdf",width=6, height = 8)
p 
dev.off()

#### fig S4a####
dat1 <- read.csv("./data/metadata.csv")
dat2 <- read.csv("./data/metadata_clinical.csv")
# cfDNA1
dat1.1 <- dat1[grepl("1",dat1$ID1),]
dat1.1$ID1 <- gsub("1","",dat1.1$ID1)
dat <- merge(dat2,dat1.1,by=c("ID1","cancer"))
dat <- dat[!(is.na(dat$time) & is.na(dat$state)),]
# CA199
dd <- na.omit(dat[,c(12,14)])
dd[,1] <- log2(dd[,1]+1)
dd[,2] <- log2(dd[,2]+1)

p <- ggscatter(dd, x = 'AFV', y = 'CA199', size = 3,palette = "jco",
               #add = "reg.line", conf.int = TRUE,
               add.params = list(color = "black", fill = "gray",
                                 size =1))+
  labs(y="log2(CA199+1)",
       x="log2(AFV+1)",
       fill = "")+
  # geom_abline(intercept = 0, slope = 1)+
  stat_cor(method = "spearman",size=7)+
  theme_linedraw()+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        panel.border = element_rect(fill=NA,color="black", linewidth=2, linetype="solid"),
        title = element_text(colour = 'black', angle = 0,size = 2),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5,size=2),
        strip.text.x = element_text(size = 32,color="black", vjust=0.5, hjust=0.5),
        axis.title.y=element_text(size=36, color="black", vjust=0.5, hjust=0.5),
        axis.title.x=element_text(size=36, color="black", vjust=0.5, hjust=0.5),
        axis.text.x= element_text(size=32, color="black", vjust=0.5, hjust=0.5,angle = 0),
        axis.text.y= element_text(size=32, color="black", vjust=0.5, hjust=0.5))+
  xlim(c(0,5.5))+ylim(c(0,15))
p
pdf("./images/figS4a_Correlation_AFV_CA199.pdf")
p
dev.off()

#### fig S4b ####
dd <- na.omit(dat[,c(11,14)])
dd[,1] <- log2(dd[,1]+1)
dd[,2] <- log2(dd[,2]+1)

p <- ggscatter(dd, x = 'AFV', y = 'CEA', size = 3,palette = "jco",
               #add = "reg.line", conf.int = TRUE,
               add.params = list(color = "black", fill = "gray",
                                 size =1))+
  labs(y="log2(CEA+1)",
       x="log2(AFV+1)",
       fill = "")+
  stat_cor(method = "spearman",size=7)+
  theme_linedraw()+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        panel.border = element_rect(fill=NA,color="black", linewidth=2, linetype="solid"),
        title = element_text(colour = 'black', angle = 0,size = 2),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5,size=2),
        strip.text.x = element_text(size = 32,color="black", vjust=0.5, hjust=0.5),
        axis.title.y=element_text(size=36, color="black", vjust=0.5, hjust=0.5),
        axis.title.x=element_text(size=36, color="black", vjust=0.5, hjust=0.5),
        axis.text.x= element_text(size=32, color="black", vjust=0.5, hjust=0.5,angle = 0),
        axis.text.y= element_text(size=32, color="black", vjust=0.5, hjust=0.5))+
  xlim(c(0,5.5))+ylim(c(0,8))
pdf("./images/figS4a_Correlation_AFV_CEA.pdf")
p
dev.off()

#### fig S5a ####
dat1 <- read.csv("./data/metadata.csv")
dat2 <- read.csv("./data/metadata_clinical.csv")

dat1.1 <- dat1[grepl("2",dat1$ID1),]
dat1.1$ID1 <- gsub("2","",dat1.1$ID1)
dat <- merge(dat2,dat1.1,by=c("ID1","cancer"))
dat <- dat[!(is.na(dat$time) & is.na(dat$state)),]
dat$state1 <- dat$state=="dead"
dat$group <- ifelse(dat$CA199 > median(dat$CA199,na.rm = T),"high","low")
table(dat$group)
y <- Surv(dat$time,dat$state=="dead")
mx.km <- survfit(y~group, data = dat, type= "kaplan-meier")
pdf("./images/figS5a_CA199_cfDNA2.pdf",width = 10,height = 8,onefile=F)
km_plot(mx.km,dat)
dev.off()

#### fig S5b ####
dat$group <- ifelse(dat$CEA > median(dat$CEA,na.rm = T),"high","low")
table(dat$group)
y <- Surv(dat$time,dat$state=="dead")
mx.km <- survfit(y~group, data = dat, type= "kaplan-meier")
km_plot(mx.km,dat)
pdf("./images/figS5b_CEA_cfDNA2.pdf",width = 10,height = 8,onefile=F)
km_plot(mx.km,dat)
dev.off()
