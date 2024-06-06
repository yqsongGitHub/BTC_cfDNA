source("./src/library.R")
source("./src/function.R")

# load data 
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

#### fig 2a ####
mut <- all1[grepl("1",all1$ID1),]
allsamples <- unique(mut$ID1)
mut <- mut[,1:3]
sample.keep <- keepsample$ID1[grepl("1",keepsample$ID1)]
pdf("./images/fig2a_cancer.pdf",width = 12 ,height = 10)
fig2a_landscape_cancer <- draw_landscape1(mut,sample.keep)
dev.off()
write.csv(fig2a_landscape_cancer,"./data/fig2a_cancer.csv")

#### fig 2b ####
mut <- all1[!grepl("1|2",all1$ID1),]
allsamples <- unique(mut$ID1)
mut <- mut[,1:3]
sample.keep <- keepsample$ID1[!grepl("1|2",keepsample$ID1)]
pdf("./images/fig2b_inflammation.pdf",width = 12 ,height = 10)
fig2a_landscape_suspect <- draw_landscape1(mut,sample.keep)
dev.off()
write.csv(fig2a_landscape_suspect,"./data/fig2a_inflammation.csv")

#### fig 2d ####
mut.gene <- unique(all1[grepl("1",all1$ID1),"GENE"])
sus.gene <- unique(all1[!grepl("1|2",all1$ID1),"GENE"])
pdf(file='./images/fig2d_cancer_inflammation_venn.pdf')
grid.draw(vennplot(mut.gene,sus.gene,setn=2,
                   cnames=c("cancer","inflammation"),
                   cols=c("#d01c8b","#66a61e"),
                   filen= './fig/tmp.tiff'))
dev.off()

#### fig 2c ####
dat <- data.frame("Type"=c("BTC","BTI"),
                  "Overlap"=c(24,24),
                  "Only"=c(27,8))
data1 <- reshape2::melt(dat)
p <- ggplot(data1, aes(x=Type,y=value,fill=variable))+
  geom_bar(position="fill", stat="identity",width = .8)+
  scale_fill_manual(#breaks=c("coding" ,"non_coding"),
  values=c("#d9d9d9","#969696")) +
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size=18, color="black", face= "bold", vjust=0.5, hjust=0.5),
        axis.title.y= element_text(size=22, color="black", face= "bold", vjust=0.5, hjust=0.5),
        axis.text.x= element_text(size=26, color="black", face= "bold", vjust=0.5, hjust=0.5, angle=90),
        axis.text.y= element_text(size=20, color="black", face= "bold", vjust=0.5, hjust=0.5),
        legend.text = element_text(size=20, color="black", face= "bold"),
        axis.title.x= element_blank(),
        plot.margin = margin(1, 0, 1, 0, "cm"))+
  labs(fill="",title = "",x="",y="")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+scale_y_continuous(labels = scales::percent)  
pdf("./images/fig2c_cancer_inflammation_intersect.pdf")
p
dev.off()

#### fig 2e ####
genes <- c("EGFR","BRAF","IDH1","ARID1A","BAP1","NF1","APC","TP53","LRP1B")
dat <- all1[!grepl("2",all1$ID1),]
colnames(dat)[2] <- "Gene_Name"
dat.fi <- fisher_func(dat)
dat.fi <- dat.fi[dat.fi$Gene %in% genes,]
dat <- dat.fi[,1:5]
data1 <- reshape2::melt(dat)
data1$type <- gsub(".*_","",data1$variable)
data1$variable <- gsub("_.*","",data1$variable)
data1$variable <- factor(data1$variable,levels = c("Not","Mut"))
data2 <- data1[data1$type=="cancer",]
ord <- data2[data2$variable=="Mut",]
ord <- ord[order(-ord$value),]
data2$Gene <- factor(data2$Gene,levels = ord$Gene)
data3 <- data1[data1$type=="normal",]
data3$Gene <- factor(data3$Gene,levels = ord$Gene)
p1 <- ggplot(data2, aes(x=Gene,y=value,fill=variable))+
  geom_bar(position="fill", stat="identity",width = .8)+
  scale_fill_manual(#breaks=c("coding" ,"non_coding"),
    #labels=c("Gene-harboring hotspots", "Non-gene hotspots"),
    values=c("#d9d9d9","#969696")) +
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size=18, color="black", face= "bold", vjust=0.5, hjust=0.5),
        axis.title.y= element_text(size=22, color="black", face= "bold", vjust=0.5, hjust=0.5),
        axis.text.x= element_text(size=22, color="black", face= "bold", vjust=0.5, hjust=0.5, angle=90),
        axis.text.y= element_text(size=20, color="black", face= "bold", vjust=0.5, hjust=0.5),
        legend.text = element_text(size=20, color="black", face= "bold"),
        axis.title.x= element_blank(),
        plot.margin = margin(1, 0, 1, 0, "cm"))+
  labs(fill="",title = "",x="",y="")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+scale_y_continuous(labels = scales::percent)+
  guides(fill="none")
p2 <- ggplot(data3, aes(x=Gene,y=value,fill=variable))+
  geom_bar(position="fill", stat="identity",width = .8)+
  scale_fill_manual(#breaks=c("coding" ,"non_coding"),
    #labels=c("Gene-harboring hotspots", "Non-gene hotspots"),
    values=c("#d9d9d9","#969696")) +
  theme(panel.grid.major =element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(size=18, color="black", face= "bold", vjust=0.5, hjust=0.5),
        axis.title.y= element_text(size=22, color="black", face= "bold", vjust=0.5, hjust=0.5),
        axis.text.x= element_text(size=22, color="black", face= "bold", vjust=0.5, hjust=0.5, angle=90),
        axis.text.y= element_text(size=20, color="black", face= "bold", vjust=0.5, hjust=0.5),
        legend.text = element_text(size=20, color="black", face= "bold"),
        axis.title.x= element_blank(),
        plot.margin = margin(1, 0, 1, 0, "cm"))+
  labs(fill="",title = "",x="",y="")+
  scale_y_continuous(expand = c(0, 0), limits = c(0, NA))+scale_y_continuous(labels = scales::percent)  
ggpubr::ggarrange(p1,p2,nrow = 1,widths = c(0.45,0.55))
ggsave("./images/fig2e_BTC_BTI_barplot.pdf",width = 10,height = 5)
