source("./src/library.R")
source("./src/function.R")

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
meta2 <- meta[!grepl("2",meta$ID1),]
meta2$type1 <- meta2$type
meta2$type1[meta2$type %in% c("cfDNA1","cfDNA2")] <- "cancer"
rownames(meta2) <- meta2$ID1
meta2 <- meta2[,-c(1,2,4)]
dat1 <- reshape2::melt(meta2)
dat1$value <- log2(dat1$value+1)
dat1 <- dat1[!is.na(dat1$type1),]

#### fig 3a;TMB ####
dat.box <- dat1[dat1$variable=="TMB",]
dat.box <- na.omit(dat.box)
p <- draw_jitter(dat.box,ytitle="log2(TMB+1)",ymax = 6)
pdf("./images/fig3a_TMB.pdf",width=6, height = 8)
p 
dev.off()

#### fig 3a;CA199 ####
dat.box <- dat1[dat1$variable=="CA199",]
dat.box <- na.omit(dat.box)
p <- draw_jitter(dat.box,ytitle="log2(CA199+1)",ymax = 15)
pdf("./images/fig3a_CA199.pdf",width=6, height = 8)
p 
dev.off()

#### fig 3a;CEA ####
dat.box <- dat1[dat1$variable=="CEA",]
dat.box <- na.omit(dat.box)
p <- draw_jitter(dat.box,ytitle="log2(CEA+1)",ymax = 7)
pdf("./images/fig3a_CEA.pdf",width=6, height = 8)
p 
dev.off()

#### fig 3b ####
dat1 <- read.csv("./data/metadata.csv")
dat2 <- read.csv("./data/metadata_clinical.csv")
dat1.1 <- dat1[grepl("1",dat1$ID1),]
dat1.1$ID1 <- gsub("1","",dat1.1$ID1)
mut <- mut[grepl("1",row.names(mut)),,drop=F]
rownames(mut) <- gsub("1","",rownames(mut))
dat <- merge(dat1.1,mut,by.x = "ID1",by.y = "row.names",all.y=T)
dd <- na.omit(dat[,c(10,6)])
dd[,1] <- log2(dd[,1]+1)
dd[,2] <- log2(dd[,2]+1)

p <- ggscatter(dd, x = 'TMB', y = 'CA199', size = 3,palette = "jco",
               add = "reg.line", conf.int = TRUE,
               add.params = list(color = "black", fill = "gray",
                                 size =1))+
  labs(
    x="log2(TMB+1)",
    y="log2(CA199+1)",
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
  xlim(c(0,4.5))+ylim(c(0,15))
pdf("./images/fig3b_Correlation_CA199_TMB.pdf")
p
dev.off()

dd <- na.omit(dat[,c(10,5)])
dd[,1] <- log2(dd[,1]+1)
dd[,2] <- log2(dd[,2]+1)
p <- ggscatter(dd, x = 'TMB', y = 'CEA', size = 3,palette = "jco",
               add = "reg.line", conf.int = TRUE,
               add.params = list(color = "black", fill = "gray",
                                 size =1))+
  labs(
    x="log2(TMB+1)",
    y="log2(CEA+1)",
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
  xlim(c(0,4.5))+ylim(c(0,6))
pdf("./images/fig3b_Correlation_CEA_TMB.pdf")
p
dev.off()

#### fig 3c ####
all1 <- read.csv("./data/Gallbladder.csv")
all1$cancer <- "Gallbladder"
all2 <- read.csv("./data/Cholangiocarcinoma.csv")
all2$cancer <- "Cholangiocarcinoma"
all0 <- rbind(all1,all2)
stat.all <- select_diff_gene(all0)
stat.all <- stat.all[stat.all$diff >0, ]
stat.all <- stat.all[1:10,]
stat.all <- stat.all[order(stat.all$FREQUENCY.x+stat.all$FREQUENCY.y,decreasing = T),]
panel.genelist <- as.character(stat.all$GENE)

# for colnames order, identical to fig2a
fig2a_landscape_cancer <- read.csv("./data/fig2a_cancer.csv",row.names = 1)
fig2a_landscape_suspect <- read.csv("./data/fig2a_inflammation.csv",row.names = 1)
fig2a_sample_order <- c(colnames(fig2a_landscape_cancer),colnames(fig2a_landscape_suspect))

all1 <- all0
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

mut <- all1[!grepl("2",all1$ID1),]
mut <- unique(mut[, c("ID1","GENE","VAR_TYPE_SX")])
mut <- data.table::dcast(mut, ID1~GENE, fun.aggregate = function(x){paste(x, collapse = ";")})
mut <- data.frame(mut)
row.names(mut) <- mut$ID1
mut <- base::t(mut[,-1])
rownames(mut) <- gsub("\\.","-",rownames(mut))
mut <- mut[panel.genelist, ]

annoInfo <- data.frame(meta[,c(1:3,6,5,9)])
rownames(annoInfo) <- annoInfo[,"ID1"]
annoInfo <- data.frame(annoInfo[,-1])
annoInfo[,2] <- log2(annoInfo[,2]+1)
annoInfo[,3] <- log2(annoInfo[,3]+1)
annoInfo[,4] <- log2(annoInfo[,4]+1)
annoInfo$type <- paste0(annoInfo$cancer,"_",annoInfo$type)
annoInfo$type <- gsub("cfDNA1|cfDNA2","cancer",annoInfo$type)
annoInfo <- annoInfo[match(fig2a_sample_order, rownames(annoInfo)),,drop=F]

mut <- data.frame(mut)
setdiffsample <- setdiff(fig2a_sample_order,colnames(mut))
mut[,setdiffsample] <- ""
mut <- mut[,match(fig2a_sample_order,colnames(mut))]
mut <- as.matrix(mut)
colsplit <- annoInfo[,1,drop=F]
colsplit$type <- gsub(".*_","",colsplit$type)
annotationsize <- 16
rownamessize <- 14

onco <- oncoPrint(mat = mut, alter_fun = alter_fun, 
                  get_type = function(x) unique(strsplit(x, ";")[[1]]),
                  col = oncocol,
                  column_names_gp =  gpar(fontsize = rownamessize, fontface='bold'),
                  row_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                  row_title_gp= gpar(fontsize = rownamessize, fontface='bold'),
                  show_pct = F, 
                  pct_gp = gpar(fontsize = annotationsize),
                  row_names_side = "right",
                  show_heatmap_legend=T,
                  gap = unit(0.1, 'mm'),
                  column_gap = unit(0.5, "mm"),
                  row_gap = unit(0.1, "mm"),
                  remove_empty_columns = F,
                  column_split = colsplit,
                  top_annotation = topAnno ,
                  right_annotation = rowAnnotation(
                    rbar = anno_oncoprint_barplot(
                      axis_param = list(
                        side = "top",
                        gp=gpar(fontsize=rownamessize+2)))),
                  column_order = fig2a_sample_order,
                  row_order = panel.genelist)
pdf("./images/fig3c_panel.pdf",width = 22 ,height = 9)
draw(onco)
dev.off()

#### fig 3d ####
#CA 199 levels (higher than 20.0 units/mL) 
#CEA levels (higher than 4.0 ng/mL)
pro <- read.xlsx("./data/protein_abbreviation.xlsx")
pro <- unique(pro[,3:4])
colnames(pro) <- c("Let3","Let1")
pro_change <- read.csv("./data/cbioportal_merge.csv")
pro_change$HGVSc <- gsub("\\..*:","",pro_change$HGVSc)
fig2a_landscape_cancer <- read.csv("./data/fig2a_cancer.csv",row.names = 1)
fig2a_landscape_suspect <- read.csv("./data/fig2a_inflammation.csv",row.names = 1)
fig2a_sample_order <- c(colnames(fig2a_landscape_cancer),colnames(fig2a_landscape_suspect))

meta <- read.csv("./data/metadata.csv")
meta <- meta[!is.na(meta$CA199),]
meta <- meta[!meta$type=="cfDNA2",]
meta <- meta[meta$CA199 > 20,]
samplelist <- fig2a_sample_order[fig2a_sample_order %in%meta$ID1]

pdf("./images/fig3d_cancer_CA199high.pdf",width = 22 ,height = 9)
draw_mutation(meta,samplelist)
dev.off()

#### fig 3e ####
meta <- meta <- read.csv("./data/metadata.csv")
meta <- meta[!is.na(meta$CA199),]
meta <- meta[!meta$type=="cfDNA2",]
meta <- meta[!meta$CA199 > 20,]
samplelist <- fig2a_sample_order[fig2a_sample_order %in%meta$ID1]
pdf("./images/fig3e_cancer_CA199low.pdf",width = 22 ,height = 9)
draw_mutation(meta,samplelist)
dev.off()

