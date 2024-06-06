convert_mut <- function(mut=mut, freq, all_sn, n = 20){
  mut <- unique(mut[, c("ID1","GENE","VAR_TYPE_SX")])
  mut <- data.table::dcast(mut, ID1~GENE, fun.aggregate = function(x){paste(x, collapse = ";")})
  mut <- data.frame(mut)
  row.names(mut) <- mut$ID1
  mut[setdiff(all_sn, mut$ID1), ] <- ""
  mut <- base::t(mut[,-1])
  freq <- head(freq, n = n)
  rownames(mut) <- gsub("\\.","-",rownames(mut))
  mut <- mut[rownames(mut) %in% freq$GENE, ]
  return(mut)
}
cal_freq <- function(mut, total){
  mut <- unique(mut[, c("ID1", "GENE")])
  mut <- as.data.frame(table(mut$GENE))
  colnames(mut) <- c("GENE", "MUTANT")
  mut$N <- total
  mut$FREQUENCY <- mut$MUTANT / mut$N
  mut$WT <- mut$N - mut$MUTANT
  mut <- mut[order(mut$FREQUENCY, decreasing = T), ]
  return(mut)
}
change_mut<- function(mut){
  mutN <- mut
  mutN[mutN==""] <- as.numeric(1)
  mutN[!mutN==1] <- as.numeric(0)
  mutN <- data.frame(apply(mutN,2,  as.numeric))
  rownames(mutN) <- row.names(mut)
  colnames(mutN) <- colnames(mut)
  mutN <- data.frame(t(mutN))
  return(mutN)
}
alter_fun <- function(x,y,w,h,v){
  n <- sum(v)
  h <-  h*0.8
  w <-  w*0.85
  grid.rect(x, y, w, h, gp = gpar(fill = "grey95", col = NA))
  if(n) grid.rect(x,y-h*0.5+1:n/n*h, w, 1/n*h, gp=gpar(fill=oncocol[names(which(v))], col=NA), just = "top")
}

varorder <- c('Missense variant', 'Splice',
              'Inframe deletion', 'Frameshift variant', 
              'Start lost','Stop gained')
oncocol <- c('Missense variant' = "#BD0031",  
             'Splice' = "#A0C9F0",           
             'Inframe deletion' = "#F38400",  
             'Frameshift variant' = "#865592",
             'Start lost' = "#F3C300",       
             'Stop gained' = "#999999")     

col_fun_cfDNA <- colorRamp2(c(0, 5), c( "white", "#e41a1c"))
col_fun_CA199 <- colorRamp2(c(0, 15), c("white", "#377eb8"))
col_fun_CEA <-  colorRamp2(c(0, 10), c("white", '#a65628'))

draw_landscape1 <- function(mut,sample.keep){
  freq <- cal_freq(mut, length(allsamples))
  mut <- convert_mut(mut,freq, allsamples, n = 20)
  mut <- data.frame(mut)
  setdiffsample <- setdiff(sample.keep,colnames(mut))
  mut[,setdiffsample] <- ""
  
  # annoInfo
  annoInfo <- data.frame(meta[,c(1:3,6,5,9)])
  rownames(annoInfo) <- annoInfo[,"ID1"]
  annoInfo <- data.frame(annoInfo[,-1])
  annoInfo[,2] <- log2(annoInfo[,2]+1)
  annoInfo[,3] <- log2(annoInfo[,3]+1)
  annoInfo[,4] <- log2(annoInfo[,4]+1)
  annoInfo$type <- paste0(annoInfo$cancer,"_",annoInfo$type)
  annoInfo$type <- gsub("cfDNA1|cfDNA2","cancer",annoInfo$type)
  
  mutN <- change_mut(mut)
  mutN <- merge(annoInfo,mutN,by="row.names",all.y = T)
  mutN$type <- gsub(".*_","",mutN$type)
  mutN <- mutN[do.call(order,c(mutN[,-c(1,3,4,5)],list(decreasing=F))),]
  
  mutN1 <- apply(mut,2,function(x){
    x1 <- strsplit(paste0(x,collapse= ";"),";")[[1]]
    return(sum(!x1==""))})
  mutN1 <- data.frame(sample=mutN1)
  mutN1 <-  merge(mutN1,mutN[,c(6,1)],by.x="row.names",by.y = "Row.names")
  mutN1 <- mutN1[,c(1,3,2)]
  mutN1 <- mutN1[do.call(order,c(mutN1[,-c(1)],list(decreasing=T))),,drop=F]
  
  mutN <- mutN[match(mutN1$Row.names,mutN$Row.names),]
  colnames(mutN)[1] <- "ID1"
  mut <- mut[, match(mutN$ID1, colnames(mut))]
  annoInfo <- annoInfo[match(mutN$ID1, rownames(annoInfo)),,drop=F]
  
  # annotation
  annotationsize <- 16
  rownamessize <- 14
  bottomAnno <- ComplexHeatmap::HeatmapAnnotation(
    Type = annoInfo[,1],
    cfDNA = annoInfo[,2],
    CA199= annoInfo[,3],
    CEA = annoInfo[,4],
    col = list(Type = c("Gallbladder_cancer" = "#E68FAC", 
                        "Cholangio_cancer" = "#B3446C", 
                        "Gallbladder_suspect" = "#8DB600",
                        "Cholangio_suspect" = "#008755"),
               cfDNA = col_fun_cfDNA,
               CA199 = col_fun_CA199,
               CEA = col_fun_CEA),
    annotation_name_gp = gpar(fontsize = annotationsize, fontface='bold'),
    show_legend = F
  )
  mut <- as.matrix(mut)
  onco <- oncoPrint(mat = mut, alter_fun = alter_fun, 
                    get_type = function(x) unique(strsplit(x, ";")[[1]]),
                    col = oncocol,
                    heatmap_legend_param = list(title='Alternations', 
                                                title_position=c('topleft'), 
                                                nrow=length(names(table(mut)))-1,
                                                title_gp=gpar(fontsize=rownamessize+2, fontface='bold'),
                                                labels_gp=gpar(fontsize=rownamessize+2)),
                    column_names_gp =  gpar(fontsize = rownamessize, fontface='bold'),
                    row_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                    row_title_gp= gpar(fontsize = rownamessize, fontface='bold'),
                    show_pct = T, 
                    pct_gp = gpar(fontsize = annotationsize),
                    row_names_side = "right",
                    show_heatmap_legend=F,
                    gap = unit(0.5, 'mm'),
                    column_gap = unit(0.5, "mm"), 
                    row_gap = unit(0.5,"mm"),
                    remove_empty_columns = F,
                    bottom_annotation = bottomAnno ,
                    right_annotation = rowAnnotation(
                      rbar = anno_oncoprint_barplot(
                        axis_param = list(
                          side = "top",
                          gp=gpar(fontsize=rownamessize+2)))),
                    top_annotation = columnAnnotation(
                      rbar = anno_oncoprint_barplot(
                        axis_param = list(
                          side = "left",
                          gp=gpar(fontsize=rownamessize+2)))),
                    column_order = mutN$ID1,
                    #row_order = all_genes_order$GENE
  )
  
  lgd_mut <- Legend(at = c("0","1","2","3","4","5"),
                    labels = c("Missense variant",
                               "Splice",
                               "Inframe deletion",
                               "Frameshift variant",
                               "Start lost",
                               'Stop gained'),
                    title = "Alternations",
                    legend_gp = gpar(fill = c("#BD0031",
                                              "#A0C9F0",
                                              "#F38400",
                                              "#865592",
                                              "#F3C300",
                                              "#C2B180"
                    )),
                    labels_gp = gpar( fontsize = annotationsize-2),
                    title_gp = gpar( fontsize = annotationsize-2, fontface="bold"))
  
  anno_legend_list = lapply(bottomAnno@anno_list,function(anno) color_mapping_legend(anno@color_mapping, plot = FALSE))
  anno_legend_list[["cfDNA"]]@grob$children[[1]]$gp$fontsize <- annotationsize-2
  anno_legend_list[["CA199"]]@grob$children[[1]]$gp$fontsize <- annotationsize-2
  anno_legend_list[["CEA"]]@grob$children[[1]]$gp$fontsize <- annotationsize-2
  anno_legend_list[["cfDNA"]]@name <- "log2(cfDNA+1)"
  anno_legend_list[["CA199"]]@name <- "log2(CA199+1)"
  anno_legend_list[["CEA"]]@name <- "log2(CEA+1)"
  # anno_legend_list[["Type"]]@grob$children[[2]]$children[[1]]$gp$fontsize <- annotationsize-2
  draw(onco, heatmap_legend_side='right',heatmap_legend_list=c(list(lgd_mut),anno_legend_list))
  return(mut)
}

vennplot <- function(set1,set2,set3=NULL,setn=2,
                     cnames=c("a","b"),
                     cols=c("#440154ff", '#21908dff', '#fde725ff'),
                     filen= 'venn.png'){
  #myCol <- brewer.pal(3, "Pastel2")
  if(is.null(set3)){
    data.list <- list(set1,set2)
  }else{
    data.list <- list(set1,set2,set3)
  }
  venn.diagram(
    x = data.list,
    category.names = cnames,
    #filename = filen,
    filename = NULL,
    output = T ,
    imagetype="tiff" ,
    height =800 , 
    width = 800 , 
    resolution = 1000,
    #compression = "lzw",
    lwd = 6,
    col=cols[1:setn],
    fill = c(alpha(cols[1],0.3), alpha(cols[2],0.3), alpha(cols[3],0.3))[1:setn],
    cex = 5,
    #main="Cholangiocarcinoma",
    #main.cex = 2.5
    #fontfamily = "sans",
    cat.cex = 0,
    cat.default.pos = "outer",
    cat.pos = c(-27, 27, 135)[1:setn],
    cat.dist = c(0.055, 0.055, 0.085)[1:setn],
    #cat.fontfamily = "sans",
    cat.col = cols[1:setn],
    #rotation = 1
  )
}
fisher_func <- function(dat){
  df0 <- data.frame()
  for(i in unique(dat$Gene_Name)){
    print(i)
    cancer <- unique(dat$ID1[grepl("1|2",dat$ID1)])
    normal <- unique(dat$ID1[!dat$ID1 %in% cancer])
    
    mutsam <- unique(dat$ID1[dat$Gene_Name %in% i])
    mut.can <- mutsam[mutsam %in% cancer]
    mut.nor <- mutsam[mutsam %in% normal]
    NN1 <- length(mut.can)
    NN2 <- length(cancer)- NN1
    NN3 <- length(mut.nor)
    NN4 <- length(normal)- NN3
    mt <- matrix(c(NN1,NN2,NN3,NN4),nrow = 2, 
                 dimnames = list(mut_or_not = c("Mut", "Not_mut"),
                                 can_or_nor = c("Can", "Nor")))
    re <- fisher.test(mt)
    df <- data.frame("Gene"= i,
                     "Mut_cancer" = NN1,
                     "Not_mut_cancer"= NN2,
                     "Mut_normal" = NN3,
                     "Not_mut_normal"= NN4,
                     "Pvalue"=re$p.value,
                     "OR"=re$estimate)
    df0 <- rbind(df0,df)
  }
  df0 <- df0[order(df0$Pvalue,decreasing = F),]
  return(df0)
}
draw_jitter <- function(dat.box,ytitle="log2(cfDNA)",ymax=5,ymin=0){
  p <- ggdotplot(dat.box, x = "type1", y = "value",
                 fill="type1",color = "type1", 
                 size = 0.8,
                 position = position_jitter(0.1),
                 #facet.by = c("variable"),
                 add.params = list(size=0.1),
                 add = "jitter"
                 #add = "mean"
  )+ 
    theme_linedraw()+
    scale_fill_manual(breaks=c("cancer", "suspect"),
                      values=c("#BD0031", "#A0C9F0"))+
    scale_color_manual(breaks=c("cancer", "suspect"),
                       values=c("#BD0031", "#A0C9F0"))+
    stat_compare_means(size=10,label.y = max(dat.box$value)+0.5) +# Add pairwise comparisons p-value
    labs(y=ytitle,fill = "")+
    ylim(c(ymin,ymax))+
    theme(panel.grid=element_blank(),
          panel.background=element_blank(),
          title = element_text(colour = 'black', angle = 0,size = 2),
          panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
          legend.position = 'none',
          axis.ticks.length=unit(0.4, "cm"),
          plot.title = element_text(hjust = 0.5,size=2),
          axis.title.x=element_blank(),
          strip.text.x = element_text(size = 32,color="black", vjust=0.5, hjust=0.5),
          axis.title.y=element_text(size=36, color="black", vjust=0.5, hjust=0.5),
          axis.text.x= element_text(size=32, color="black", vjust=0.5, hjust=0.5,angle = 0),
          axis.text.y= element_text(size=32, color="black", vjust=0.5, hjust=0.5),
          axis.line = element_line(colour = "black", size = 1, linetype = "solid"))+
    scale_x_discrete(labels = c("cancer" = paste0("Cancer","\n","( n = ",table(dat.box$type1)[1],")"),
                                "suspect" = paste0("Suspect","\n","( n = ",table(dat.box$type1)[2],")")))+ 
    stat_summary(fun.data = function(x) 
      data.frame(y=round(mean(x),2),
                 ymin=round(mean(x),2),
                 ymax=round(mean(x),2)),
      geom="crossbar",
      width=0.5,
      color=c("black", "black")
      #color=c("#BD0031", "#A0C9F0")
    )
  return(p)
}
draw_jitter1 <- function(dat.box,ytitle="log2(cfDNA)",ymax=5,ymin=0){
  p <- ggdotplot(dat.box, x = "type1", y = "value",
                 fill="type1",color = "type1", 
                 size = 0.8,
                 position = position_jitter(0.1),
                 #facet.by = c("variable"),
                 add.params = list(size=0.1),
                 add = "jitter"
                 #add = "mean"
  )+ 
    theme_linedraw()+
    scale_fill_manual(breaks=c("cfDNA1", "cfDNA2"),
                      values=c("#BD0031", "#A0C9F0"))+
    scale_color_manual(breaks=c("cfDNA1", "cfDNA2"),
                       values=c("#BD0031", "#A0C9F0"))+
    stat_compare_means(size=10,label.y = max(dat.box$value)+0.5) +# Add pairwise comparisons p-value
    labs(y=ytitle,fill = "")+
    ylim(c(ymin,ymax))+
    theme(panel.grid=element_blank(),
          panel.background=element_blank(),
          title = element_text(colour = 'black', angle = 0,size = 2),
          panel.border = element_rect(fill=NA,color="black", size=2, linetype="solid"),
          legend.position = 'none',
          axis.ticks.length=unit(0.4, "cm"),
          plot.title = element_text(hjust = 0.5,size=2),
          axis.title.x=element_blank(),
          strip.text.x = element_text(size = 32,color="black", vjust=0.5, hjust=0.5),
          axis.title.y=element_text(size=36, color="black", vjust=0.5, hjust=0.5),
          axis.text.x= element_text(size=32, color="black", vjust=0.5, hjust=0.5,angle = 0),
          axis.text.y= element_text(size=32, color="black", vjust=0.5, hjust=0.5),
          axis.line = element_line(colour = "black", size = 1, linetype = "solid"))+
    scale_x_discrete(labels = c("cfDNA1" = paste0("cfDNA1","\n","( n = ",table(dat.box$type1)[1],")"),
                                "cfDNA2" = paste0("cfDNA2","\n","( n = ",table(dat.box$type1)[2],")")))+ 
    stat_summary(fun.data = function(x) 
      data.frame(y=round(mean(x),2),
                 ymin=round(mean(x),2),
                 ymax=round(mean(x),2)),
      geom="crossbar",
      width=0.5,
      color=c("black", "black")
      #color=c("#BD0031", "#A0C9F0")
    )
  return(p)
}
select_diff_gene <- function(all0){
  all1 <- all0
  all1 <- all1[grepl("1",all1$ID1),]
  all1 <- all1[all1$ID1 %in% meta$ID1,]
  colnames(all1)[17] <- "GENE"
  mut <- all1
  mut <- unique(mut[, c("ID1", "GENE")])
  mut <- as.data.frame(table(mut$GENE))
  colnames(mut) <- c("GENE", "MUTANT")
  mut$N <- length(unique(all1$ID1))
  mut$FREQUENCY <- mut$MUTANT / mut$N
  mut <- mut[order(mut$FREQUENCY, decreasing = T), ]
  stat.cancer <- mut
  
  all1 <- all0
  all1 <- all1[!grepl("1|2",all1$ID1),]
  all1 <- all1[all1$ID1 %in% meta$ID1,]
  colnames(all1)[17] <- "GENE"
  mut <- all1
  mut <- unique(mut[, c("ID1", "GENE")])
  mut <- as.data.frame(table(mut$GENE))
  colnames(mut) <- c("GENE", "MUTANT")
  mut$N <- length(unique(all1$ID1))
  mut$FREQUENCY <- mut$MUTANT / mut$N
  mut <- mut[order(mut$FREQUENCY, decreasing = T), ]
  stat.sus <- mut
  
  stat <- merge(stat.cancer,stat.sus,by="GENE",all=T)
  stat$FREQUENCY.x[is.na(stat$FREQUENCY.x)] <- 0
  stat$FREQUENCY.y[is.na(stat$FREQUENCY.y)] <- 0
  stat$MUTANT.x[is.na(stat$MUTANT.x)] <- 0
  stat$MUTANT.y[is.na(stat$MUTANT.y)] <- 0
  stat$diff <- stat$FREQUENCY.x - stat$FREQUENCY.y
  stat$mut_sum <- stat$MUTANT.x + stat$MUTANT.y
  stat <- stat[order(stat$diff,decreasing = T),]
  return(stat)
}
topAnno <- ComplexHeatmap::HeatmapAnnotation(
  type= annoInfo[,1],
  cfDNA = annoInfo[,2],
  CA199= annoInfo[,3],
  CEA=annoInfo[,4],
  col = list(type = c("Gallbladder_cancer" = "#E68FAC", 
                      "Cholangio_cancer" = "#B3446C", 
                      "Gallbladder_suspect" = "#8DB600",
                      "Cholangio_suspect" = "#008755"),
             cfDNA=col_fun_cfDNA,
             CA199=col_fun_CA199,
             CEA=col_fun_CEA),
  annotation_name_gp =gpar(fontsize = annotationsize, fontface='bold'),
  show_legend = T
)
alter_fun <- function(x,y,w,h,v){
  n = sum(v)
  h = h*0.9
  w = w*0.9
  grid.rect(x, y, w, h, gp = gpar(fill = "grey95", col = NA))
  if(n) grid.rect(x,y-h*0.5+1:n/n*h, w, 1/n*h, gp=gpar(fill=oncocol[names(which(v))], col=NA), just = "top")
}
draw_panel <- function(meta,all0,panel.genelist){
  panel.genelist <- sort(panel.genelist)
  
  all1 <- all0
  all1 <- all1[,c("ID1","Gene_Name","Annotation")]
  colnames(all1) <- c("ID1","GENE","VAR_TYPE_SX")
  mut <- all1[!grepl("2",all1$ID1),]
  mut <- unique(mut[, c("ID1","GENE","VAR_TYPE_SX")])
  mut <- data.table::dcast(mut, ID1~GENE, fun.aggregate = function(x){paste(x, collapse = ";")})
  row.names(mut) <- mut$ID1
  mut <- base::t(mut[,-1])
  mut <- mut[panel.genelist, ]
  
  annoInfo <- data.frame(meta[,c(1,2,9)])
  rownames(annoInfo) <- annoInfo[,"ID1"]
  annoInfo <- annoInfo[,-1,drop=F]
  annoInfo$type <- paste0(annoInfo$cancer,"_",annoInfo$type)
  annoInfo$type <- gsub("cfDNA1|cfDNA2","cancer",annoInfo$type)
  
  mutN <- t(mut)
  mutN <- merge(annoInfo,mutN,by="row.names")
  colnames(mutN)[1] <- "ID1"
  mutN <- mutN[do.call(order,c(mutN[,-c(1)],list(decreasing=T))),]
  mutN <- mutN[order(mutN[,2],decreasing=F),]
  
  fig2a_sample_order <- c(colnames(fig2a_landscape_cancer),colnames(fig2a_landscape_suspect))
  mutN <- mutN[match(fig2a_sample_order,mutN$ID1),]
  
  mut <- mut[, match(mutN$ID1, colnames(mut))]
  annoInfo <- annoInfo[match(mutN$ID1, rownames(annoInfo)),,drop=F]
  
  colsplit <- annoInfo[,1,drop=F]
  colsplit$type <- gsub(".*_","",colsplit$type)
  
  annotationsize <- 16
  rownamessize <- 20
  topAnno <- ComplexHeatmap::HeatmapAnnotation(
    type= annoInfo[,1],
    col = list(type = c("Gallbladder_cancer" = "#E68FAC", 
                        "Cholangio_cancer" = "#B3446C", 
                        "Gallbladder_suspect" = "#8DB600",
                        "Cholangio_suspect" = "#008755")
    ),
    annotation_name_gp =gpar(fontsize = 0, fontface='bold'),
    show_legend = F
  )
  
  alter_fun <- function(x,y,w,h,v){
    n=sum(v)
    h = h*0.9
    w = w*0.9
    grid.rect(x, y, w, h, gp = gpar(fill = "grey95", col = NA))
    if(n) grid.rect(x,y-h*0.5+1:n/n*h, w, 1/n*h, gp=gpar(fill=oncocol[names(which(v))], col=NA), just = "top")
  }
  
  onco <- oncoPrint(mat = mut, alter_fun = alter_fun, 
                    get_type = function(x) unique(strsplit(x, ";")[[1]]),
                    col = oncocol,
                    heatmap_legend_param = list(title='Alternations', 
                                                title_position=c('topleft'), 
                                                #nrow=length(names(table(mut)))-1,
                                                title_gp=gpar(fontsize=rownamessize+2, fontface='bold'),
                                                labels_gp=gpar(fontsize=rownamessize+2)),
                    column_names_gp =  gpar(fontsize = rownamessize, fontface='bold'),
                    row_names_gp = gpar(fontsize = rownamessize, fontface='bold'),
                    row_title_gp= gpar(fontsize = rownamessize, fontface='bold'),
                    column_title_gp= gpar(fontsize = rownamessize,fontface='bold'),
                    show_pct = F, 
                    pct_gp = gpar(fontsize = annotationsize),
                    row_names_side = "left",
                    show_heatmap_legend=T,
                    gap = unit(0.1, 'mm'),
                    column_gap = unit(2, "mm"),
                    row_gap = unit(0.1, "mm"),
                    remove_empty_columns = F,
                    column_split = colsplit,
                    top_annotation = topAnno ,
                    right_annotation=NULL,
                    column_order = mutN$ID1,
                    row_order = panel.genelist)
  draw(onco)
  return(mut)
}

change_pro_abbre <- function(old_cha){
  ncha <- c()
  for(i in old_cha){
    if(i==""){
      cha <- i
    }else{
      cha0 <- gsub("\\D+","",i)
      cha <- gsub("^p\\.","",i)
      cha <- strsplit(cha,"[0-9]+")[[1]]
      if(cha[1] %in% pro$Let3){
        cha1 <- pro$Let1[cha[1]==pro$Let3]}
      else{
        cha1 <- cha[1]
      }
      if(cha[2] %in% pro$Let3){
        cha2 <- pro$Let1[cha[2]==pro$Let3]
      }else{
        cha2 <- cha[2]
      }
      cha <- paste0("p.",cha1,cha0,cha2)
    }
    ncha <- c(ncha,cha)
  }
  return(ncha)
}
draw_mutation <- function(meta,fig2a_sample_order){
  dt1 <- read.csv("./data/Gallbladder.csv")
  dt2 <- read.csv("./data/Cholangiocarcinoma.csv")
  dt1$cancer <- "Gallbladder"
  dt2$cancer <- "Cholangiocarcinoma"
  dat <- rbind(dt1,dt2)
  dat <- dat[!grepl("2",dat$ID1),]
  panel.genelist <- read.csv("./data/panel.genelist.csv")
  panel.genelist <- panel.genelist$x
  
  all0 <- dat
  all0$HGVS.p <- change_pro_abbre(all0$HGVS.p)
  all0$HGVS.c <- paste0(all0$Feature_ID, all0$HGVS.c)
  all0$HGVS.p1 <-  pro_change$Protein.Change[match(all0$HGVS.c,pro_change$HGVSc)]
  all0$HGVS.p[all0$HGVS.p==""] <- all0$HGVS.p1[all0$HGVS.p==""]
  all0$HGVS.p[is.na(all0$HGVS.p)] <- ""
  all0$Gene_Name <- paste0(all0$Gene_Name," (",all0$HGVS.p,")")
  all0 <- all0[all0$ID1 %in% meta$ID1,]
  
  stat.cho <- select_diff_gene(all0)
  stat.cho <- stat.cho[!grepl("\\(\\)",stat.cho$GENE),]
  stat.cho <- stat.cho[order(-stat.cho$mut_sum),]
  panel.genelist1 <-  as.character(stat.cho$GENE)
  panel.genelist1 <- panel.genelist1[gsub(" \\(.*\\)","",panel.genelist1) %in% panel.genelist]
  #rank panel.genelist1
  panel.genelist2 <- gsub(" \\(.*\\)","",panel.genelist1)
  panel.genelist3 <- c()
  for (i in panel.genelist) {
    panel.index <- grep(i,panel.genelist2)
    panel.genelist3 <- c(panel.genelist3,panel.index)
  }
  panel.genelist1 <- panel.genelist1[panel.genelist3]
  
  all1 <- all0[,c("ID1","Gene_Name","Annotation","cancer")]
  colnames(all1) <- c("ID1","GENE","VAR_TYPE_SX","Cancer")
  all1 <- all1[!is.na(all1$ID1),]
  mut <- all1[!grepl("2",all1$ID1),]
  mut <- unique(mut[, c("ID1","GENE","VAR_TYPE_SX")])
  mut <- data.table::dcast(mut, ID1~GENE, fun.aggregate = function(x){paste(x, collapse = ";")})
  row.names(mut) <- mut$ID1
  mut <- base::t(mut[,-1])
  mut <- mut[panel.genelist1, ]
  
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
  topAnno <- ComplexHeatmap::HeatmapAnnotation(
    type= annoInfo[,1],
    cfDNA = annoInfo[,2],
    CA199= annoInfo[,3],
    CEA=annoInfo[,4],
    col = list(type = c("Gallbladder_cancer" = "#E68FAC", 
                        "Cholangio_cancer" = "#B3446C", 
                        "Gallbladder_suspect" = "#8DB600",
                        "Cholangio_suspect" = "#008755"),
               cfDNA=col_fun_cfDNA,
               CA199=col_fun_CA199,
               CEA=col_fun_CEA),
    annotation_name_gp =gpar(fontsize = annotationsize, fontface='bold'),
    show_legend = T
  )
  
  alter_fun <- function(x,y,w,h,v){
    n=sum(v)
    h = h*0.9
    w = w*0.9
    grid.rect(x, y, w, h, gp = gpar(fill = "grey95", col = NA))
    if(n) grid.rect(x,y-h*0.5+1:n/n*h, w, 1/n*h, gp=gpar(fill=oncocol[names(which(v))], col=NA), just = "top")
  }
  
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
                    column_gap = unit(2, "mm"),
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
                    row_order = panel.genelist1)
  draw(onco)
  print(mut)
}
km_plot <- function(mx.km,dat){
  p <- ggsurvplot(mx.km, data = dat, palette = c("#e41a1c", "#377eb8"), 
                  conf.int = F,
                  pval.size = 11,
                  xlab = 'Time (Days)', 
                  conf.int.style='step', legend.title='', title = NULL,
                  legend =  c(0.8, 0.9),legend.labs = c("High", "Low"),
                  pval = T, pval.method = F,
                  risk.table = T, tables.height = 0.25, 
                  fontsize = 8,
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
  print(p)
}
dist2d <- function(a,b,c) {
  v1 <- b - c
  v2 <- a - b
  m <- cbind(v1,v2)
  d <- abs(det(m))/sqrt(sum(v1*v1))
  return(d)
} 
AFV <- function(x,y){
  dd0 <- c()
  NN <- length(x)
  for(i in 1:NN){
    pp0 <- c(as.numeric(x[i]),as.numeric(y[i]))
    dd1 <- dist2d(pp0,c(0,0),c(1,1)) 
    dd0 <- c(dd0,dd1)
  }
  return(median(dd0,na.rm = T))
}




