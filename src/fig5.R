
#### fig 5a ####
dat <- read.csv("./data/metadata.csv")
dat <- dat[,c(1,8)]
dat <- na.omit(dat)
colnames(dat)[1] <- "sample"
dat$AFV <- log2(dat$AFV+1)
dat$type <- "suspect"
dat$type[grepl("1",dat$sample)] <- "cfDNA1"
dat$type[grepl("2",dat$sample)] <- "cfDNA2"
dat.box <- dat[!dat$type %in% "suspect",]
colnames(dat.box) <- c("sample","value","type1")
dat.box <- na.omit(dat.box)
dat.box$pair <- gsub("1|2","",dat.box$sample)
dat.box1 <- reshape2::dcast(dat.box,pair~type1)
dat.box1 <- na.omit(dat.box1)
colnames(dat.box1) <- c("pair","Before","After")
dat.box1$diff <- dat.box1$After - dat.box1$Before
table(dat.box1$diff >0)

pdf("./images/fig5a_AFV_cfDNA1_cfDNA2.pdf")
ggpaired(dat.box1, cond1 = "Before", cond2 = "After",
         line.color = "gray",
         point.size = 3,
         line.size = 1,
         fill=c("#BD0031", "#A0C9F0"))+
  #stat_compare_means(paired = TRUE)+
  labs(y="AFV", x="", fill = "")+
  ylim(c(0,6))+
  theme_linedraw()+
  theme(panel.grid=element_blank(),
        panel.background=element_blank(),
        panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
        title = element_text(colour = 'black', angle = 0,size = 2),
        legend.position = 'none',
        plot.title = element_text(hjust = 0.5,size=2),
        axis.title.x=element_blank(),
        strip.text.x = element_text(size = 32,color="black", face= "bold", vjust=0.5, hjust=0.5),
        axis.title.y=element_text(size=46, color="black", face= "bold", vjust=0.5, hjust=0.5),
        axis.text.x= element_text(size=32, color="black", face= "bold", vjust=0.5, hjust=0.5,angle = 0),
        axis.text.y= element_text(size=32, color="black", face= "bold", vjust=0.5, hjust=0.5),
        #axis.line = element_line(colour = "black", size = 1, linetype = "solid")
  )
dev.off()

#### fig 5b ####
dat1 <- read.csv("./data/metadata.csv")
dat1.1 <- dat1[grepl("2",dat1$ID1),]
dat1.1$ID1 <- gsub("2","",dat1.1$ID1)
dat <- merge(dat2,dat1.1,by=c("ID1","cancer"))
dat <- dat[!(is.na(dat$time) & is.na(dat$state)),]
dat$state1 <- dat$state=="dead"
dat$group <- ifelse(dat$AFV > median(dat$AFV,na.rm = T),"high","low")
y <- Surv(dat$time,dat$state=="dead")
mx.km <- survfit(y~group, data = dat, type= "kaplan-meier")
pdf("./images/fig5b_cfDNA2.pdf",width = 10,height = 8,onefile=F)
km_plot(mx.km,dat)
dev.off()

#### fig 5c ####
dat1 <- read.csv("./data/metadata.csv")
dat2 <- read.csv("./data/metadata_clinical.csv")

dat1.1 <- dat1[grepl("1",dat1$ID1),]
dat1.1$ID1 <- gsub("1","",dat1.1$ID1)
dat <- merge(dat2,dat1.1,by=c("ID1","cancer"))
dat <- dat[!(is.na(dat$time) & is.na(dat$state)),]
dat$state1 <- dat$state=="dead"
dat.cfDNA1 <- dat

dat1.1 <- dat1[grepl("2",dat1$ID1),]
dat1.1$ID1 <- gsub("2","",dat1.1$ID1)
dat <- merge(dat2,dat1.1,by=c("ID1","cancer"))
dat <- dat[!(is.na(dat$time) & is.na(dat$state)),]
dat$state1 <- dat$state=="dead"
dat.cfDNA2 <- dat

dat.Roc1 <- survivalROC(Stime=dat.cfDNA1$time,  
                        status=dat.cfDNA1$state=="dead",      
                        marker = dat.cfDNA1$AFV,  
                        predict.time = 1500,
                        span = 0.25*NROW(dat.cfDNA1)^(-0.20),
                        method="KM")
dat.Roc2 <- survivalROC(Stime=dat.cfDNA2$time,  
                        status=dat.cfDNA2$state=="dead",      
                        marker = dat.cfDNA2$AFV,  
                        predict.time = 1500,
                        span = 0.25*NROW(dat.cfDNA2)^(-0.20),
                        method="KM")

dat.Roc.list <- list(dat.Roc1,dat.Roc2)
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
rownames(df0) <- c("cfDNA1","cfDNA2")
write.csv(df0,"./data/cfDNA1_cfDNA2_sensitivity_specificity.csv")

pdf("./images/fig5c_cfDNA1_cfDNA2.pdf",width=10,height = 10)
par(mar=c(6,6,6,6))
plot(dat.Roc1$FP*100, dat.Roc1$TP*100, 
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
lines(dat.Roc2$FP*100, dat.Roc2$TP*100, lwd = 4,type="l",col="#377eb8",xlim=c(0,100), ylim=c(0,100))
legend(55,10,
       c(paste("cfDNA1 =",round(dat.Roc1$AUC,3)),
         paste("cfDNA2 =",round(dat.Roc2$AUC,3))),
       x.intersp=1, y.intersp=0.8,
       lty= 1 ,lwd= 3,col=c("#b30000","#377eb8"),
       bty = "n",
       seg.len=1,cex=2.0)
dev.off()





