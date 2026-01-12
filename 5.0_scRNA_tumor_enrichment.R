
#~/R/my_projects/JAK_STAT_BC/5.0_scRNA_tumor

if(T){
  rm(list = ls())
  setwd("~/R/my_projects/JAK_STAT_BC/5.0_scRNA_tumor")
  folder_path <- "./out_tumor"
  # 检查文件夹是否存在
  if (!dir.exists(folder_path)) {
    # 如果文件夹不存在，则创建它
    dir.create(folder_path)
    print(paste("Folder created at", folder_path))
  } else {
    print(paste("Folder already exists at", folder_path))
  }
  library(ggplot2);library(ggpubr);library(RColorBrewer);library(tibble)
  my36colors <-c("#0072B5","#E18727","#F37C95","#20854E", '#D6E7A3', '#57C3F3', '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 
  library(Seurat);library(SCP);library(qs);library(dplyr)
}

if(T){
  #
  genes<-c(
    "AKT3", "SPRY3", "SPRY1", "SPRY2", "STAM2", "IRF9", "PIAS3", "IL24", 
    "CISH", "IL22RA2", "SOCS4", "CNTF", "CNTFR", "CREBBP", "CSF2", 
    "CSF2RA", "CSF2RB", "CSF3", "CSF3R", "CSH1", "CTF1", "IL23R", 
    "SPRED1", "IFNLR1", "SPRED2", "EP300", "EPO", "EPOR", "AKT1", 
    "AKT2", "CLCF1", "PIK3R5", "CBLC", "GH1", "GH2", "GHR", 
    "IFNL2", "IFNL3", "IFNL1", "GRB2", "IL19", "SOCS7", "IFNE", 
    "IFNA1", "IFNA2", "IFNA4", "IFNA5", "IFNA6", "IFNA7", "IFNA8", 
    "IFNA10", "IFNA13", "IFNA14", "IFNA16", "IFNA17", "IFNA21", 
    "IFNAR1", "IFNAR2", "IFNB1", "IFNG", "IFNGR1", "IFNGR2", "IFNW1", 
    "IL2", "IL2RA", "IL2RB", "IL2RG", "IL3", "IL3RA", "IL4", 
    "IL4R", "IL5", "IL5RA", "IL6", "IL6R", "IL6ST", "IL7", 
    "IL7R", "IL9", "IL9R", "IL10", "IL10RA", "IL10RB", "IL11", 
    "IL11RA", "IL12A", "IL12B", "IL12RB1", "IL12RB2", "IL13", 
    "IL13RA1", "IL13RA2", "IL15", "IL15RA", "JAK1", "JAK2", "JAK3", 
    "LEP", "LEPR", "LIF", "LIFR", "MPL", "MYC", "OSM", 
    "IL20", "IL21R", "IL22", "IL23A", "PIAS4", "PIK3CA", "PIK3CB", 
    "PIM1", "PIK3CD", "PIK3CG", "PIK3R1", "PIK3R2", "IL20RA", 
    "IL20RB", "IL26", "PRL", "PRLR", "IFNK", "PTPN6", "PTPN11", 
    "IL22RA1", "IL21", "CCND1", "BCL2L1", "CRLF2", "SOS1", 
    "SOS2", "STAT1", "STAT2", "STAT3", "STAT4", "STAT5A", 
    "STAT5B", "STAT6", "TPO", "TYK2", "STAM", "SPRY4", "PIK3R3", 
    "TSLP", "PIAS1", "SOCS1", "CBL", "CBLB", "SOCS2", "CCND2", 
    "CCND3", "SOCS3", "PIAS2", "OSMR", "SOCS5"
  )
  genes_1 <- c(
    "JAK1", "JAK2", "JAK3", 
    "STAT1", "STAT2", "STAT3", 
    "STAT4", "STAT5A", "STAT5B", 
    "STAT6"
  )
  interest_gene <-"STAT4"
  
  my_cancer<-c("BRCA")
  
  genelist<- list("Jak-STAT_pathway"=genes,"Core_Jak-STAT_pathway"=genes_1
                  #"Proliferation"= c("BIRC5", "CCNB1", "CDC20", "NUF2", "CEP55", "NDC80", "MKI67", "PTTG1", "RRM2", "TYMS","UBE2C")
                  ) #https://www.nature.com/articles/s41588-021-00911-1
  input_gene= c("Jak-STAT_pathway","Core_Jak-STAT_pathway") #genes#c(genes,"H3K36me","5mC" )
  project = "Jak_STAT_BC"
}
#qsave(sceList,paste0(folder_path,"/sce_Tumor.qs"))
#sce_tumor <-qread(paste0(folder_path,"/sce_Tumor.qs"))#only tumor
#sceList <-qread("sceList_score.qs")#all
sceList <-qread(paste0(folder_path,"/sce_Tumor.qs"))#only tumor

#0.load data----
if(F){
  sceList <-qread("sceList_score.qs") #qread("sce_Epi.qs")#8k
  table(sceList@meta.data[["celltype_major"]])#tumor 24489
  sceList <-subset(sceList,celltype_major=="Cancer Epithelial")
  table(sceList@meta.data[["celltype_minor"]])
  # Cancer Basal SC  Cancer Cycling  Cancer Her2 SC  Cancer LumA SC  Cancer LumB SC 
  # 4312            5359            3708            7742            3368 
  ##score----
  #load genelist http://biocc.hrbmu.edu.cn/CancerSEA/goDownload
  gene_set<-read.csv("./out_tumor/CancerSEA.csv")
  #gene_set$gene <-rownames(gene_set)
  list<- split(as.matrix(gene_set)[,2], gene_set[,1])
  list<- c(list,genelist)
  
  library(BiocParallel)
  register(MulticoreParam(workers = 24, progressbar = TRUE))#32
  #score
  sceList <- AddModuleScore(object = sceList,#sce,
                            features = list,#gene,
                            name = names(list )#"score" #ctrl = 100,
  ) #
  #colnames(sceList@meta.data)[15] <- 'Score' #ncol(sceList@meta.data)
  colnames(sceList@meta.data)[16:32] <- gsub("\\d","",colnames(sceList@meta.data)[16:32] )
  sceList$Score <-sceList$`Jak-STAT_pathway`
  
}
# 1.0 score and group----
##分组High low; 查看比例和分布；通路活性相关性
sceList$group <- ifelse(sceList$Score > mean(sceList$Score),"High","Low")
table(sceList$group)
# High   Low 
# 11646 12843 

##barplot 柱状图 比例
library(scRNAtoolVis)
Idents(sceList)<-sceList$celltype_minor
names(sceList@reductions)[3] <-"UMAP"#"umap"

##柱状图barplot in every sample ----
cellRatioPlot(object = sceList,
              sample.name = "orig.ident",
              celltype.name = "celltype_minor", #col.width  fill.col
              fill.col=my36colors#c("darkgoldenrod2","seagreen","steelblue") 
)+theme(axis.text.x = element_text(angle=90, hjust=1,vjust = 1))
ggsave(paste0(folder_path,"/cellRatioPlot_major.pdf"),width = 7,height = 3,family="serif")#4.5-3 #74 #44  #6.5/9 5

cellRatioPlot(object = subset(sceList,sceList$group=="High"),
              sample.name = "orig.ident",
              celltype.name = "celltype_minor", #col.width  fill.col
              fill.col=my36colors#c("darkgoldenrod2","seagreen","steelblue") 
)+
  theme(axis.text.x = element_text(angle=90, hjust=1,vjust = 1))+
  cellRatioPlot(object = subset(sceList,sceList$group=="Low"),
                sample.name = "orig.ident",
                celltype.name = "celltype_minor", #col.width  fill.col
                fill.col=my36colors#c("darkgoldenrod2","seagreen","steelblue") 
  )+
  theme(axis.text.x = element_text(angle=90, hjust=1,vjust = 1))
ggsave(paste0(folder_path,"/cellRatioPlot_tumor_high-low.pdf"),width = 14,height = 3.5,family="serif")#4.5-3 #74 #44  #6.5/9 5

#mean ratio in high/low
Idents(sceList) <-sceList$celltype_minor
Cellratio <- prop.table(table(Idents(sceList), sceList$group), margin = 2)  
Cellratio <-data.frame(Cellratio)
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
colnames(cellper)[1]<-"id"
colourCount = length(unique(Cellratio$Var1)) #3

ggplot(Cellratio) + 
  geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",width = 0.7,size = 0.5,colour = NA)+ 
  #geom_bar(aes(x =Var2, y= Freq, fill = Var1),stat = "identity",position = "dodge")+ 
  theme_classic2() +
  labs(x='',y = 'Ratio',fill="")+
  scale_fill_manual(values = my36colors )+ # brewer.pal(12,"Paired")
  scale_color_manual(values = my36colors )
#coord_flip()+
#theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"))
ggsave(paste0(folder_path,"/Bar_major_Score_high-low.pdf"),width = 3,height = 2.5,family="serif")#4.5-3 #74 #44  #6.5/9 5


### 查看整体分布Ro/e---
library(Startrac)
#remotes::install_local("sscVis-master")
#remotes::install_local("master.zip")
R_oe <- calTissueDist(sceList@meta.data,
                      byPatient = F,
                      colname.cluster = "celltype_minor",
                      colname.patient = "orig.ident",
                      colname.tissue = "group",
                      method = "chisq", 
                      min.rowSum = 0) 
R_oe
if(T){
  library(Startrac)
  library(ggplot2)
  library(tictoc)
  library(ggpubr)
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(circlize)
  library(tidyverse)
  library(sscVis)
}
col_fun <- colorRamp2(c(min(R_oe, na.rm = TRUE), 1, max(R_oe, na.rm = TRUE)),
                      c("skyblue", "white", "pink"))
#col_fun <- colorRamp2(c(0, 0.5, 1.5, 2), c("blue", "green",  "orange", "red"))

p <-Heatmap(as.matrix(R_oe),
            show_heatmap_legend = TRUE, 
            cluster_rows = TRUE, 
            cluster_columns = TRUE,
            row_names_side = 'right', 
            show_column_names = TRUE,
            show_row_names = TRUE,
            col = col_fun,
            row_names_gp = gpar(fontsize = 12),
            column_names_gp = gpar(fontsize = 12),
            heatmap_legend_param = list(
              title = "Ro/e",  # 自定义图注名称
              at = seq(0.5, 2, by = 0.5), # 例刻度的位置/自己的数据必须修改一下！
              labels = seq(0.5, 2, by = 0.5) # 每个刻度的标签/自己的数据必须修改一下！
            ),
            cell_fun = function(j, i, x, y, width, height, fill) {
              grid.text(sprintf("%.2f", R_oe[i, j]), x, y, gp = gpar(fontsize = 10, col = "black"))
            }
)
pdf(paste0(folder_path,"/cellRatioPlot_major_Score_All_high-low-Roe.pdf"),width = 3.5,height = 3,family="serif")#4.5-3 #74 #44  #6.5/9 5
print(p)
dev.off()

### 显示 + - 号
p <-Heatmap(as.matrix(R_oe),
            show_heatmap_legend = TRUE, 
            cluster_rows = TRUE, 
            cluster_columns = TRUE,
            row_names_side = 'right', 
            show_column_names = TRUE,
            show_row_names = TRUE,
            col = col_fun,
            row_names_gp = gpar(fontsize = 12),  # 坐标字体
            column_names_gp = gpar(fontsize = 12), # 坐标字体
            column_names_rot = 90, # 将列名倾斜 45°
            heatmap_legend_param = list(
              title = "Ro/e",  
              at = seq(0.5, 2, by = 0.5), 
              labels = seq(0.5, 2, by = 0.5) 
            ),
            cell_fun = function(j, i, x, y, width, height, fill) {
              value <- R_oe[i, j]
              if (value > 1) {
                label <- "+++"
              } else if (value > 0.8) {
                label <- "++"
              } else if (value >= 0.2) {
                label <- "+"
              } else if (value > 0) {
                label <- "+/-"
              } else if (value == 0) {
                label <- "-"
              } else {
                label <- NA  # 如果有负值或其他异常值
              }
              grid.text(label, x, y, gp = gpar(fontsize = 16 , col = "black"))  # 调整里面
            }
)
pdf(paste0(folder_path,"/cellRatioPlot_major_Score_All_high-low-Roe-1.pdf"),width = 3.5,height = 3,family="serif")#4.5-3 #74 #44  #6.5/9 5
print(p)
dev.off()
## +++, Ro/e > 1;
## ++, 0.8 < Ro/e ≤ 1;
## +, 0.2 ≤ Ro/e ≤ 0.8;
## +/−, 0 < Ro/e < 0.2;
## −, Ro/e = 0

##exp in samples and subtype----
#样本得分和表达: dotplot and violin
for (j in c("Score",colnames(sceList@meta.data)[17:30] )){ #,interest_gene colnames(sceList@meta.data)[16:30]
  print(j)
  for (i in c("orig.ident","subtype") ){
    print(i)
    ##
    library(MySeuratWrappers)
    Idents(sceList) <-sceList$orig.ident
    VlnPlot(sceList, features = j, group.by  = i,
            stacked=T,pt.size= 0, sort='increasing',
            cols = my36colors, #颜色
            #direction = "horizontal", #水平作图
            x.lab = '', y.lab = '')+NoLegend()+labs(y=j,title="")+theme(axis.text.x = element_text(angle = 90,size = 12))
    width = ifelse(length(unique(sceList@meta.data[,i]))/5 >4,length(unique(sceList@meta.data[,i]))/5,3.5)
    height = ifelse(nchar(max(as.character(unique(sceList@meta.data[,i]))))>10,4,3.5)
    ggsave(paste0(folder_path,"/VlnPlot_Score_",j,"_",i,".pdf"),width = width,height = height,family="serif")
    
    #dotplot
    DotPlot(sceList , features =j,group.by = i,
            #cols = c("skyblue", "pink")
    )+ #Minor_Type orig.ident
      #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
      #scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
      scale_color_viridis_c(option = "A")+
      labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
      theme_bw()+
      theme(#panel.grid = element_blank(), 
        axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
    hei1 <-ifelse(i=="orig.ident",height,height/2)
    ggsave(paste0(folder_path,"/Dotplot_Score_",j,"_",i,".pdf"),width = 3,height = hei1,family="serif")
    
    ##
  }
  
}
##dot plot all score----
#display.brewer.all()
for (i in c("group") ){ #"group","celltype_minor"
  print(i)
  
  for (j in c("Blues","Reds","Greens","PuBuGn","GnBu") ){
    print(j)
    DotPlot(sceList , 
            #features =c( colnames(sceList@meta.data)[17:30] ) , #"Score",
            features =c("Score", colnames(sceList@meta.data)[17:30] ) , #
            group.by = i,
            #cols = c("skyblue", "pink")
    )+coord_flip()+ #theme_bw()+ #去除背景，旋转图片
      #scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
      #scale_color_viridis_c()+
      scale_color_gradientn(colors = brewer.pal(6, j) ) +#"YlOrRd"
      labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
      theme_bw(base_size = 14)+
    theme(axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
    #hei1 <-ifelse(i=="orig.ident",height,height/2)
    #ggsave(paste0(folder_path,"/Dotplot_all_Score_",j,"_",i,".pdf"),width = 4,height = 3.5,family="serif")
    #ggsave(paste0(folder_path,"/Dotplot_all_Score_minor_",j,"_",i,".pdf"),width = 4.5,height = 4.5,family="serif")
    ggsave(paste0(folder_path,"/Dotplot_all_Score_group_",j,"_",i,".pdf"),width = 4,height = 3.5,family="serif")
    
  }
   
}
##Cor 计算相关性----
#x="Score";y=unique(gene_set$Pathway)[1]
for (y in unique(gene_set$Pathway) ){
  print(y)
  x="Score"
  ##
  methods <-c("pearson","spearman")
  for (j in methods){
    print(j)
    df <-sceList@meta.data[c(x,y,"celltype_minor")]
    df <- na.omit(df) #删除没有配对的NA
    colnames(df)<-c("X","Y","Type")
    ggplot(df, aes(x=X, y=Y)) + 
      #xlim(-20,15) + ylim(-15,10) +
      labs(x = x, y = y) +
      geom_point(aes(color=Type),size = 1.0,alpha=0.2) + #,color="skyblue"
      geom_smooth(method ='lm', size=0.5,color="red") + #navy
      #ggrepel::geom_text_repel(aes(label=cancer))+ #geom_label_repel
      stat_cor(method = j,size = 5,color="red") +
      scale_colour_manual(values = my36colors ) +
      theme_bw() + 
      theme(axis.text.x = element_text(size = 14), 
            axis.text.y = element_text(size = 14), 
            axis.title.x = element_text(size = 16), legend.position = "none",
            axis.title.y = element_text(size = 16),panel.grid = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(size = 12))
    ggsave(file=paste0(folder_path,"/Cor_point_",x,"_",y,"_",j,".pdf"),width = 3,height = 3,family="serif")
  }
  ##
}
##Cor heatmap----
library(ggcorrplot)

df_cor <- sceList@meta.data[,c("Score",unique(gene_set$Pathway) ) ]
cormtcars <- round(cor(df_cor), 3) #round()函数自定义小数点后位数

if(T){
  #ggcorrplot(cormtcars,lab=T)
  #ggcorrplot(cormtcars,method = "circle",lab=T)
  
  pmtcars <- cor_pmat(df_cor)
  ggcorrplot(cormtcars,hc.order = T,method = "circle",  #分等级聚类重排矩阵
             ggtheme = ggplot2::theme_void(base_size = 15), #主题修改
             colors = c("CornflowerBlue","white","Salmon"), #自定义颜色，看自己喜欢，或是参考好看的文献Figure用法。
             lab = T,lab_size = 4,    #相关系数文本字体大小
             tl.cex = 12,             #坐标轴字体大小
             p.mat = pmtcars,         #添加显著性信息
             sig.level = 0.01,        #显著性水平
             outline.color = "white",#外框颜色 "gray"
             insig = "blank",pch = 4,                 #不够显著的色块进行标记，pch表示选择不同的标记方法，可以尝试其他数字表示什么标记方法
             pch.cex = 10)            #不显著标记的大小，使用insig = "blank"将不显著的空白处理
  ggsave(file=paste0(folder_path,"/ggcorrplot-p.pdf"),width = 6.5,height = 6.5,family="serif")
  
  ##显示一半##corrplot
  library(corrplot)
  mycol2 <- colorRampPalette(c("navy","white", "pink"), alpha = TRUE)
  
  pdf(file=paste0(folder_path,"/corrplot-p.pdf"),width = 5,height = 5,family="serif")
  corrplot(cormtcars, method = c('square'), type = c('lower'), 
           col = mycol2(100),
           outline = F,#'grey', #是否为图形添加外轮廓线，默认FLASE，可直接TRUE或指定颜色向量
           order = c('AOE'), #排序/聚类方式选择："original", "AOE", "FPC", "hclust", "alphabet"
           diag = FALSE, #是否展示对角线结果，默认TRUE
           tl.cex = 1.2, #文本标签大小
           tl.col = 'black', #文本标签颜色
           addgrid.col= 'grey' #格子轮廓颜色
  )
  dev.off()
  
  pdf(file=paste0(folder_path,"/corrplot-p-1.pdf"),width = 5,height = 5,family="serif")
  corrplot(cormtcars, method = c('square'), type = c('lower'), 
           col = mycol2(100),
           outline = 'grey', #是否为图形添加外轮廓线，默认FLASE，可直接TRUE或指定颜色向量
           order = c('AOE'), #排序/聚类方式选择："original", "AOE", "FPC", "hclust", "alphabet"
           diag = FALSE, #是否展示对角线结果，默认TRUE
           tl.cex = 1.2, #文本标签大小
           tl.col = 'black', #文本标签颜色
           addgrid.col= 'grey' #格子轮廓颜色
  )
  #下三角图添加不显著叉号：
  corrplot(cormtcars, add = TRUE,
           #method = c('number'), 
           type = c('lower'),
           col = mycol2(100),
           order = c('AOE'), 
           diag = FALSE, 
           number.cex = 0.9,#tl.cex=5,
           tl.pos = 'n', 
           cl.pos = 'n',
           p.mat = pmtcars,
           insig = "pch"
  )
  dev.off()
  
  #COR PHEATMAP GGPLOT2
  #
  library(tidyverse)
  library(RColorBrewer)
  library(readr)
  # 拆分数据
  df1 <- as.data.frame(cormtcars)#read_tsv('./out/data_1.tsv',skip = 1) %>% slice(1:8)
  df2 <-as.data.frame(pmtcars)
  df2$`R-value` <- rownames(df2)
  df2 <- df2 %>%  #p matrix
    pivot_longer(-`R-value`) %>% 
    mutate(value=as.numeric(value),
           p_signif=symnum(value, corr = FALSE, na = FALSE,
                           cutpoints = c(0, 0.001, 0.01,0.05,1),
                           symbols = c("***", "**", "*", "")))
  # 定义颜色
  col <- colorRampPalette(brewer.pal(11,"RdBu")[3:9])(100)
  # 热图绘制
  p1 <- df1 %>% pivot_longer(-`R-value`) %>% 
    mutate(name=factor(name,levels = unique(name)),
           `R-value` =factor(
             `R-value`,levels = unique(`R-value`))) %>% 
    mutate(value=as.numeric(value)) %>% 
    ggplot(aes(`R-value`,name)) +
    #geom_tile(fill="white",color="black",linewidth=0.3) +
    geom_point(aes(size=value,color=value)) +
    geom_text(data=df2,aes(label=p_signif),color="black") +
    scale_color_gradientn(colours =rev(col),
                          breaks=c(0,0.4,0.8),
                          labels=c(0,0.4,0.8)) +
    guides(size="none") +
    #scale_y_discrete(expand = c(0,0)) +scale_x_discrete(expand = c(0,0)) +
    labs(x=NULL,y=NULL) +
    theme_test() +
    theme(axis.ticks = element_blank(),
          legend.title = element_blank(),
          axis.text.y=element_text(color="black"),
          axis.text.x=element_text(
            angle=90,color="black",vjust=0.5),
          plot.margin = margin(0.5,0,0.5,0.5, unit = "cm"))
  p1
  ggsave(paste0(folder_path,"/corrplot-p-2.pdf"),width = 6,height = 5,family="serif")
  
  #select only Score
  df1$`R-value` <- rownames(df1)
  df1 <- df1[c("Score","R-value")]
  df1 <- df1[-1,]
  df2 <- df2[df2$`R-value` =="Score",]
  df2 <- df2[-1,]
  p1 <- df1 %>% pivot_longer(-`R-value`) %>% 
    mutate(name=factor(name,levels = unique(name)),
           `R-value` =factor(
             `R-value`,levels = unique(`R-value`))) %>% 
    mutate(value=as.numeric(value)) %>% 
    ggplot(aes(`R-value`,name)) +
    #geom_tile(fill="white",color="black",linewidth=0.3) +
    geom_point(aes(size=value,color=value)) +
    geom_text(data=df2,aes(x=name,y=`R-value`,label=p_signif),color="black") +
    #scale_color_gradientn(colours =rev(col),breaks=c(-0.2,0,0.2,0.4,0.6,0.8,1),labels=c(-0.2,0,0.2,0.4,0.6,0.8,1)) +
    scale_color_gradientn(colours =rev(col) )+
    guides(size="none") +
    #scale_y_discrete(expand = c(0,0)) +scale_x_discrete(expand = c(0,0)) +
    labs(x=NULL,y=NULL) +
    theme_test(base_size = 14) +
    theme(axis.ticks = element_blank(),
          legend.title = element_blank(),
          axis.text.y=element_text(color="black"),
          #axis.text.x=element_text(angle=90,color="black",vjust=0.5,hjust = 1),
          plot.margin = margin(0.5,0,0.5,0.5, unit = "cm"))+
    coord_flip()
  p1
  ggsave(paste0(folder_path,"/corrplot-p-3.pdf"),width = 2.5,height = 3.5,family="serif")
  
}




##boxplot TNBC vs Non-TNBC----
unique(sceList@meta.data$subtype)
df_epi_score<-sceList@meta.data#[sceList@meta.data$celltype_major %in% c("Normal Epithelial","Cancer Epithelial"),]
df_epi_score$Cell_type <-ifelse(sceList@meta.data$subtype=="TNBC","TNBC","Non-TNBC")
table(df_epi_score$Cell_type)
# Non-TNBC     TNBC 
# 13653    10836 

library(ggplot2);library(EnvStats);library(ggpubr)
#table(all_merge$Efficacy);colnames(all_merge)
ggplot(df_epi_score,
       aes(x = Cell_type,#reorder(Cell_type,-Score), 
           y =Score,fill =Cell_type) )+ #,fill =Efficacy y=Srps
  # geom_point(aes(color =Cell_type),alpha=0.5,size=0.2,
  #            position=position_jitterdodge(jitter.width = 0.45,
  #                                          jitter.height = 0,
  #                                          dodge.width = 0.8))+
  geom_boxplot(alpha=0.7,width=0.55,
               position=position_dodge(width=0.8),
               size=0.05,outlier.colour = NA)+
  labs(y="Score",x= NULL)+
  geom_violin(alpha=0.4,width=0.9,
              position=position_dodge(width=0.8),
              size=0.25)+
  stat_n_text()+ #library(EnvStats) # 显示样本数NR 38  N 31
  scale_color_manual(values = c( "#56B4E9","#F37C95") )+
  scale_fill_manual(values = c(  "#56B4E9","#F37C95") )+
  #scale_color_manual(values = cols)+
  #facet_grid(.~Treatment,space = "free",scales = "free")+ #根据treatment分页
  stat_compare_means(#aes(group = Group) ,   #library(ggpubr)
    #comparisons=list(c("Non-responder","Responder")),#my_comparisons,
    label = "p.format",#"p.format p.signif 1.1e-06
    method = "wilcox.test", #wilcox t
    show.legend= F,#删除图例中的"a"
    label.x=1.2,bracket.size=0.1,vjust=0.1,#label.y=0.35,
    #label.y = max(log2(GDSC2$Camptothecin_1003)),
    hide.ns = T,size=4)+
  #ylim(c(0.1,4))+ #y轴范围0.36
  #scale_x_discrete(labels= c("Normal","Cancer"))+
  theme_bw(base_size = 14,base_family = "serif")+#facet_grid(~Treatment)+
  theme(axis.title  = element_text(size=14,family = "serif"), axis.text = element_text(size=12),legend.position ="none",  
        panel.grid = element_blank(),legend.key = element_blank(),axis.title.x =element_blank() )#+ 
ggsave(paste0(folder_path,"/boxplot_score_epi_vs_tumor.pdf"),width = 2.5,height = 2.5,family="serif")


#save qs

# 2.0 Enrichment----
table(sceList$group)
##2.1 subtype enrichment and high/low in subtype----
#for (i in unique(sce$Minor_Type) ){
if(T){
  #print(i)
  sce1 <-sceList#subset(sce,sce$Minor_Type==i)
  #sce1$group <- ifelse(sce1$Score > mean(sce1$Score),"High","Low")
  for (j in c("group") ){
    print(j)
    ##run
    sce1 <- RunDEtest(srt = sce1, group_by = j, fc.threshold = 1, only.pos = FALSE)#slow
    VolcanoPlot(srt = sce1, group_by = j)
    wid=ifelse(length(unique(sce1@meta.data[,j]) )==3,14,9 )
    ggsave(paste0(folder_path,"/VolcanoPlot_subtype_",i,"_",j,".pdf"),width = wid,height = 4.5,family="serif")
    
    #DEGs <- paste0("sce1@tools$DEtest_",j,"$AllMarkers_wilcox")
    DEGs <- sce1@tools$DEtest_group$AllMarkers_wilcox#[[2]]$AllMarkers_wilcox
    write.csv(DEGs,paste0(folder_path,"/DEGs_subtype_",i,"_",j,"_all.csv"))
    DEGs <- DEGs[with(DEGs, avg_log2FC > 1/2 & p_val < 0.05), ]#p_val_adj abs(avg_log2FC)
    write.csv(DEGs,paste0(folder_path,"/DEGs_subtype_",i,"_",j,".csv"))
    
    # Annotate features with transcription factors and surface proteins
    sce1 <- AnnotateFeatures(sce1, species = "Homo_sapiens",#"Mus_musculus", 
                             db = c("GO_BP", "KEGG")#c("TF") #, "CSPA"
    )
    ht <- FeatureHeatmap(
      srt = sce1, group.by = j, features = DEGs$gene, feature_split = DEGs$group1,
      species = "Homo_sapiens",cluster_columns = T,cluster_rows = T,#"Mus_musculus", 
      db = c("GO_BP", "KEGG"), anno_terms = TRUE, #, "WikiPathway"
      #feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
      height = 5, width = 4
    )
    print(ht$plot) 
    ggsave(paste0(folder_path,"/FeatureHeatmap_subtype_",i,"_",j,".pdf"),width = 14,height = 6,family="serif")
    
    #enrich in subtype
    sce1 <- RunEnrichment(
      srt = sce1, group_by = j, db = c("GO_BP","KEGG"), #species = "Mus_musculus",
      DE_threshold = "avg_log2FC > log2(1.5) & p_val_adj < 0.05"
    )
    EnrichmentPlot(
      srt = sce1, group_by = j,db = c("GO_BP"),only_sig = T, #group_use = c("Ductal", "Endocrine"),
      plot_type = "bar", topTerm = 10,palette ="Paired"#palcolor =my36colors
    )
    ggsave(paste0("./out/EnrichmentPlot_subtype_",i,"_",j,"_GO.pdf"),width = wid,height = 4,family="serif")
    EnrichmentPlot(
      srt = sce1, group_by = j,db = c("KEGG"), only_sig = T, #group_use = c("Ductal", "Endocrine"),
      plot_type = "bar", topTerm = 20,palette ="Paired"#palcolor =my36colors
    )
    ggsave(paste0(folder_path,"/EnrichmentPlot_subtype_",i,"_",j,"_KEGG.pdf"),width = wid,height = 8,family="serif")
    #dotplot
    EnrichmentPlot(srt = sce1, group_by = j,  db = c("GO_BP","KEGG"),only_sig = T,
                   #theme_use = ggplot2::theme_classic, theme_args = list(base_size = 10),
                   palette="Blues",plot_type = "comparison") #YlOrRd
    ggsave(paste0(folder_path,"/EnrichmentPlot_subtype_",i,"_",j,"_KEGG_dot.pdf"),width = 14,height = 6,family="serif")
    
    ##irGSEA
    library(irGSEA)
    Idents(sce1)<-sce1@meta.data[,j]
    sce1<- irGSEA.score(object = sce1, assay = "RNA", slot = "data", seeds = 123, 
                        ncores = 16, min.cells = 3, min.feature = 0, custom = F, 
                        geneset = NULL, msigdb = T, species = "Homo sapiens", 
                        category = "H", subcategory = NULL, geneid = "symbol", 
                        method = c(  "ssgsea"), #"AUCell", "UCell", "singscore",#slow
                        aucell.MaxRank = NULL, ucell.MaxRank = NULL, kcdf = 'Gaussian') 
    #整合差异基因集
    result.dge <- irGSEA.integrate(object = sce1, 
                                   group.by = j, 
                                   metadata = NULL, col.name = NULL, method = c(  "ssgsea"))#"AUCell", "UCell", "singscore",
    #可视化
    result.dge[["ssgsea"]][["Name"]] <-gsub("HALLMARK-","",result.dge[["ssgsea"]][["Name"]])
    irGSEA.heatmap.plot <- irGSEA.heatmap( object= result.dge,
                                           cluster.color=my36colors[1:length(unique(result.dge[["ssgsea"]][["cluster"]]))],
                                           direction.color=c("skyblue","pink"),#significance.color = c("grey70","pink"),
                                           heatmap.width = 10,rowname.fointsize = 8,
                                           method= "ssgsea", top= 50,show.geneset= NULL)
    pdf(paste0(folder_path,"/irGSEA_HALLMARK_subtype_",i,"_",j,".pdf"),width = 6,height = 8,family="serif")
    print(irGSEA.heatmap.plot)
    dev.off()
    
    
    ##end run
  }

  
}

##2.2 enrichment of DEG in your interest subtype----
input_subtype <-"Tumor"
#Idents(sce1)<-sce1@meta.data[,j]
#DEG <-read.csv("./out/DEGs_subtype_Basal_group.csv",row.names = 1)
DEG <- FindAllMarkers(sce1, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
#DEG %>% group_by(group) #%>% top_n(n = 3, wt = avg_log2FC)

input_gene_enrich <-DEG[DEG$cluster=="High",]$gene #192
gene_type<-"ALL" #

#run
###A GO----
library(clusterProfiler)
trans = bitr(input_gene_enrich, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")
library(dplyr)
eego <- enrichGO(gene          = trans$ENTREZID,
                 #universe      = names(geneList),
                 OrgDb         = "org.Hs.eg.db",
                 ont           = "ALL",
                 pAdjustMethod = "BH",
                 #pvalueCutoff  = 0.01,
                 qvalueCutoff  = 0.05,
                 readable      = TRUE)
#plot
df_go <- data.frame(eego) %>%
  group_by(ONTOLOGY) %>%
  slice_head(n = 10) %>% #前10个,前5个？
  arrange(desc(pvalue))
ratio <- lapply(df_go$GeneRatio,function(x){as.numeric(eval(parse(text = x)))}) %>% unlist()
df_go$ratio <- ratio
df_go$Description <- factor(df_go$Description,levels = df_go$Description)
#plot
library(ggplot2)
ggplot(df_go) +
  ggforce::geom_link(aes(x = 0,y = Description,
                         xend = -log10(pvalue),yend = Description,
                         alpha = after_stat(index),#stat(index),
                         color = ONTOLOGY,
                         size = 10),#after_stat(index)),#流星拖尾效果
                     n = 500,
                     #color = "#FF0033",
                     show.legend = F) +
  geom_point(aes(x = -log10(pvalue),y = Description),
             color = "black",
             fill = "white",size = 4,shape = 21) +
  geom_line(aes(x = ratio*100,y = Description,group = 1),
            orientation = "y",linewidth = 1,color = "#FFCC00") + #线条-比例 #FFCC00 #E59CC4
  scale_x_continuous(sec.axis = sec_axis(~./100,
                                         labels = scales::label_percent(),
                                         name = "Percent of geneRatio")) +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.text = element_text(face = "bold.italic"),
        axis.text = element_text(color = "black")) +
  ylab("") + xlab("-log10 Pvalue") +
  facet_wrap(~ONTOLOGY,scales = "free",ncol = 1) + #3
  scale_color_brewer(palette = "Set2")#OrRd Set1 Paired
#scale_color_manual(values = c('#CCE0F5', '#CCC9E6', '#625D9E') )
ggsave(paste0(folder_path,"/DEG_",input_subtype,"_",gene_type,"_GO.pdf"),family="serif",width = 6,height = 7)

###B KEGG----
kegg <- enrichKEGG(
  gene          = trans$ENTREZID,
  keyType     = "kegg",
  organism   = 'hsa',
  pvalueCutoff      = 0.05,
  pAdjustMethod     = "BH",
  qvalueCutoff  = 0.05
)
library(DOSE)
kegg<-setReadable(kegg, OrgDb = "org.Hs.eg.db", keyType="ENTREZID")

#2.1 enrichplot
logFC <- DEG[ DEG$gene %in% input_gene_enrich,]$avg_log2FC
names(logFC) <- DEG[ DEG$gene %in% input_gene_enrich,]$gene 

library(enrichplot)
pdf(paste0(folder_path,"/DEG_",input_subtype,"_",gene_type,"_KEGG_cnetplot.pdf"),family="serif",width = 10,height = 6)
enrichplot::cnetplot(kegg,showCategory = 10,#foldChange = logFC,
                     circular=F,colorEdge = TRUE)
enrichplot::cnetplot(kegg,showCategory = 10,foldChange = logFC,
                     circular=F,colorEdge = TRUE)
# cnetplot(kegg, showCategory = 10, #选择top10的pathway ，这里也可以用包含pathway名称的向量           
#          color.params = list(foldChange = logFC, #用上面的logFC值标注基因的颜色
#                                   edge = T)) #显示pathway的标注
heatplot(kegg, foldChange=logFC,
         #showCategory = 8,
         symbol = "dot", #"rect" 矩形
         pvalue = NULL,
         label_format = 30)
dev.off()


ego <- kegg@result
ego <- ego [order(ego$Count,decreasing = T),]
ego <- ego[ego$pvalue <=0.05,] #P > 0.05
library(dplyr)
ego <- ego %>% mutate(row_number = row_number())

ego <- ego[c(2,4,6,10,17,32,45,49,71),]
library(ggplot2);library(ggsci);library("scales")
ego$Description <- factor(ego$Description,levels = c(ego$Description))

#2.2 barplot
mytheme<- theme(axis.title = element_text(size = 13),
                axis.text = element_text(size = 11),
                plot.title = element_text(size = 14,
                                          hjust= 0.5,
                                          face= "bold"),
                legend.title = element_text(size = 13),
                legend.text = element_text(size = 11))

for (i in c("RdPu","PRGn","PiYG","RdYlGn","Blues","Greens") ){
  print(i)
  p<- ggplot(data = ego,
             aes(x = Count, y = Description, fill = -log10(pvalue)) )+
    #scale_fill_viridis_c(alpha = 0.7)+
    scale_fill_distiller(palette = i,direction = -1) +#RCOLORBREWER #RdPu YlOrRd PiYG
    ##"RdPu"紫色 PRGn紫色绿色 #PiYG 紫色绿色  RdYlGn橙绿 
    geom_bar(stat = "identity", width = 0.8) +
    labs(x = "Count", y = "", title = "") + 
    theme_classic(base_size = 14)+# mytheme+
    theme(axis.text.y = element_blank(),axis.ticks.y = element_blank() )+
    geom_text(size=4, aes(x = 0.05, label = Description), hjust = 0, color="black" ) # hjust = 0,左对齐
    
  p
  ggsave(paste0(folder_path,"/DEG_",input_subtype,"_",gene_type,"_",i,"_KEGG_bar-1.pdf"),
         family="serif",width = 5,height = 3)#7 3.5
  
  #2.2 dotplot
  p<- ggplot(data = ego,
             aes(x = Count, y = Description) )+ #, fill = -log10(pvalue)
    scale_color_distiller(palette = i,direction = -1) + #RdPu YlOrRd PiYG
    ##PRGn紫色绿色 #PiYG 紫色绿色  RdYlGn橙绿 
    #geom_bar(stat = "identity", width = 0.8) +
    geom_point(aes(size=Count,color=-log10(pvalue) ))+
    labs(x = "Count", y = "", title = "KEGG") + 
    theme_bw()+ mytheme#+
  #scale_color_viridis_c(alpha = 0.7)
  p
  ggsave(paste0(folder_path,"/DEG_",input_subtype,"_",gene_type,"_",i,"_KEGG_dot.pdf"),family="serif",width = 7,height = 4)
  
}

#2.4 barplot+基因
mytheme1 <- theme(
  axis.title = element_text(size = 13),
  axis.text = element_text(size = 11),
  axis.text.y = element_blank(), # 在自定义主题中去掉 y 轴通路标签:
  axis.ticks.length.y = unit(0,"cm"),
  plot.title = element_text(size = 13, hjust = 0.5, face = "bold"),
  #legend.title = element_text(size = 13),
  #legend.text = element_text(size = 11),
  legend.position = "none",
  plot.margin = margin(t = 5.5, r = 10, l = 5.5, b = 5.5)
)
library(RColorBrewer)
p <- ggplot(data = ego, aes(x = -log10(pvalue), y = Description, fill =Description )) +  #, fill = Cluster
  #scale_fill_manual(values =c('#6bb9d2', '#d55640')) +
  #scale_fill_manual(values =my36colors) +
  scale_fill_manual(values =rev( colorRampPalette(brewer.pal(11, "Paired"))(length(ego$Description)) )  ) +
  geom_bar(stat = "identity", width = 0.5, alpha = 0.8) +
  scale_x_continuous(expand = c(0,0)) + # 调整柱子底部与y轴紧贴
  #scale_y_continuous(expand = c(0,-1)) +
  labs(x = "-log10(pvalue)", y = "Pathway", title = "KEGG") +
  # x = 0.61 用数值向量控制文本标签起始位置
  geom_text(size=3.8, aes(x = 0.05, label = Description),hjust = 0) + # hjust = 0,左对齐
  #geom_text(size=3, aes(x = 0.05, label = geneID), hjust = 0, vjust =2.5, color="skyblue" ) + # hjust = 0,左对齐
  geom_text(size=3, aes(x = 0.05, label = geneID), hjust = 0, vjust =-1, color="skyblue" ) + # hjust = 0,左对齐
  theme_classic() + 
  mytheme1 
#ylim(min(ego$Description),max(ego$Description))
#scale_y_discrete(labels = \(x) sub(pattern = "_", replacement = " ", x, fixed = TRUE)) 
p
ggsave(paste0(folder_path,"/DEG_",gene_type,"_KEGG_bar_genes.pdf"),width = 5,height = 4,family="serif")#7.5

###C HALLMARK ?----

#compare clusterprofiler

#3.0 CellChat----
##prepare data----
##计算样本中肿瘤细胞平均得分
# av <-AverageExpression(sce_tumor , 
#                        assays = "RNA")#样本中基因的平均表达
# av=av[[1]] 
# cg=names(tail(sort(apply(av, 1, sd)),1000)) 
# pheatmap::pheatmap(cor(av[cg,]))
Srps <-data.frame(srps=sce_tumor$Score,sample=sce_tumor$orig.ident )#99301-1 如何计算每个样本比值？
table(Srps$sample)

mean_srps<-data.frame()
for (i in unique(Srps$sample)){
  mean_1 <- data.frame(Srps=mean(Srps[Srps$sample==i,]$srps),id=i)
  mean_srps<-rbind(mean_srps,mean_1)
} #20 sample
mean_srps$Type <-ifelse(mean_srps$Srps > mean(mean_srps$Srps),"High","Low")
table(mean_srps$Type)
# High  Low #Tumor mean Score
# 12    8
#plot
ggplot(mean_srps,aes(x=reorder(id,-Srps),y=Srps,color=Type))+
  geom_point()+labs(x="",y="Mean Score")+
  geom_hline(yintercept = mean(mean_srps$Srps),color="gray",linetype="dashed", size = 0.5 )+
  theme_classic(base_size = 14) +theme(axis.text.x = element_text(angle = 90,hjust = 1,vjust = 0.5) )+
  scale_color_manual(values = c("pink","skyblue") )+
  coord_flip()
ggsave(paste0(folder_path,"/Tumor_sample_mean_score-1.pdf"),family="serif",width = 3,height = 3.5)#5 2.5
write.csv(mean_srps,paste0(folder_path,"/Tumor_sample_mean_score.csv") )

##merge in sce all
sceList$Type <- ifelse(sceList$orig.ident %in% mean_srps[mean_srps$Type=="High",]$id,"High","Low" )
# High   Low 
# 54489 44815 

sceList$Cell_type <-sceList$celltype_major


##run cellchat----
if(T){
  library(qs)
  library(Seurat)
  library(ggplot2)
  library(clustree)
  library(cowplot)
  library(dplyr)
  library(data.table)
  library(tibble)
  library(SCP)
  library(CellChat)
  library(patchwork)
  my36colors <-c( '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 
  color_used = c("#B84D64","#864A68","#E32D32","#5E549A","#8952A0","#384B97","#911310",
                 "#7CA878","#35A132","#6B70B0","#20ACBD","#959897",'skyblue',
                 "#F4A2A3","#B6CCD7","#AF98B5","#E01516","#FFFFFF")
}

for (i in unique(sceList@meta.data$Type) ){
  library(BiocParallel)
  register(MulticoreParam(workers = 32, progressbar = TRUE)) #32
  library()
  print(i)
  if(T){
    data.input =GetAssayData(sceList, assay = "RNA", slot = "data")#sceList@assays$RNA@data #data_humanSkin$data  # normalized data matrix sceList@assays$RNA@data#
    meta =sceList@meta.data #data_humanSkin$meta # a dataframe with rownames containing cell mata data
    
    cell.use = rownames(meta[meta$Type == i,]) #分组！  #rownames(meta)[meta$condition == 'LS']# extract the cell names from disease data
    data.input = data.input[, cell.use]# Subset the input data for CelChat analysis
    meta = meta[cell.use, ]
    meta$labels<-meta$Cell_type
    unique(meta$labels) # check the cell labels
    
    ##
    cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels")
    cellchat <- setIdent(cellchat, ident.use = 'labels') #将label设置为显示的默认顺序
    CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
    #showDatabaseCategory(CellChatDB)
    #选择数据库中特定子集进行细胞通讯分析:
    CellChatDB.use <- subsetDB(CellChatDB, search = 'Secreted Signaling') #可选择Secreted Signaling、ECM-Receptor或Cell-Cell Contact
    #将数据库添加到CellChat对象中(DB)
    cellchat@DB <- CellChatDB.use
    print("------------initial--------------------")
    # 对信号基因的表达数据取子集以节省计算成本
    cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
    future::plan("multisession", workers = 32) #4 do parallel
    # 识别过表达基因
    cellchat <- identifyOverExpressedGenes(cellchat)#slow
    # 识别配体-受体对
    cellchat <- identifyOverExpressedInteractions(cellchat)
    # 将配体、受体投射到PPI网络
    cellchat <- projectData(cellchat, PPI.human)
    cellchat <- computeCommunProb(cellchat, raw.use = TRUE)#slow >20min
    # triMean is used for calculating the average gene expression per cell group. 
    # Error in data.use[RsubunitsV, ] : subscript out of bounds #有负值，改代码GetAssayData(sceList, assay = "RNA", slot = "data")
    ###如果特定细胞群中只有少数细胞，则过滤掉细胞间的通信
    cellchat <- filterCommunication(cellchat, min.cells = 10)
    # df.net <- subsetCommunication(cellchat)
    # df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) 
    # df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb"))
    cellchat <- computeCommunProbPathway(cellchat)
    cellchat <- aggregateNet(cellchat)
    print("------------start plot--------------------")
    groupSize <- as.numeric(table(cellchat@idents))
    
    pdf(paste0(folder_path,"/cellchat_circle_",i,".pdf"),width = 10,height = 5,family = "serif")
    par(mfrow = c(1,2), xpd=TRUE)
    netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
    netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
    dev.off()
    
    pdf(paste0(folder_path,"/cellchat_circle_cell_",i,".pdf"),width = 10,height = 10,family = "serif")
    mat <- cellchat@net$weight
    par(mfrow = c(3,4), xpd=TRUE)
    for (y in 1:nrow(mat)) {
      mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
      mat2[y, ] <- mat[y, ]
      netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
    }
    dev.off()
    
    # Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
    cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") #slow
    
    pdf(paste0(folder_path,"/cellchat_heatmap_pattern_",i,".pdf"),width = 10,height = 6,family = "serif")
    ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing",color.heatmap ="Reds")
    ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming",color.heatmap ="Reds")
    ht1 + ht2
    dev.off()
    
    # show all the significant interactions (L-R pairs) from some cell groups (defined by 'sources.use') to other cell groups (defined by 'targets.use')
    pdf(paste0(folder_path,"/cellchat_bubble_tumor_",i,".pdf"),width = 6,height = 6,family = "serif")
    netVisual_bubble(cellchat,  sources.use = 5, remove.isolate = T)#sources.use = 4, targets.use = c(5:11),
    dev.off()
    #netVisual_bubble(cellchat,  sources.use = 3, remove.isolate = T)#sources.use = 4, targets.use = c(5:11),
    
    # pdf(paste0(folder_path,"/cellchat_scatter_dot_",i,".pdf"),width = 5,height = 5,family = "serif")
    # netAnalysis_signalingRole_scatter(cellchat)
    # dev.off()
    
    saveRDS(cellchat, file = paste0(folder_path,"/cellchat_",i,".rds") ) #save
    print("------------save success--------------------")
  }
  #gc()
}

##compare two groups----
#load 
if(T){
  #rm(list = ls())
  #setwd("~/R/my_projects/RCC/6.0_cellchat")
  library(CellChat)
  library(ComplexHeatmap)
  library(patchwork)
  high <-readRDS(paste0(folder_path,"/cellchat_High.rds") ) #请加载cellchat对应各分组对象
  low <-readRDS(paste0(folder_path,"/cellchat_Low.rds") )
  high<- updateCellChat(high)
  low<- updateCellChat(low)
  #查看两个数据集细胞分群，并检查是否一致：
  identical(levels(high@idents),levels(low@idents))
  #table(cellchat@idents$joint )
  #levels(cellchat@idents$joint)#
  levels(high@idents)
  number_cancer_cluster=3 #选择肿瘤细胞
}

#手动-调整图大小
if(T){
  ####1.0 load data----
  #合并cellchat对象:
  object.list <- list(low = low,
                      high = high) #对照组(NL)在前，比较组(LS)在后，注意顺序
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  cellchat
  
  ###1.总体比较：通讯数目与通讯强度差异---
  print("1.总体比较：通讯数目与通讯强度差异")
  p1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
  p2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
  #pdf("./out/cellchat_bar_high-low.pdf",width = 4,height = 3,family = "serif")
  p1 + p2
  #dev.off()
  ggsave(paste0(folder_path,"/cellchat_bar_high-low.pdf"),width = 4,height = 3,family = "serif")
  print(paste0(folder_path,"cellchat_bar_high-low.pdf"))
  ####2.细胞亚群水平的通讯差异----
  print("2.细胞亚群水平的通讯差异")
  #2.1 细胞通讯差异网络图
  pdf(paste0(folder_path,"/cellchat_net_high-low.pdf"),width = 10,height = 5,family = "serif")
  par(mfrow = c(1,2), xpd = TRUE)
  netVisual_diffInteraction(cellchat, weight.scale = T)
  netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
  dev.off()
  #2.2 细胞通讯差异热图
  pdf(paste0(folder_path,"/cellchat_net_heatmap_high-low.pdf"),width = 6,height = 3.5,family = "serif")
  p3 <- netVisual_heatmap(cellchat,font.size = 10,font.size.title = 12)
  p4 <- netVisual_heatmap(cellchat, measure = "weight",font.size = 10,font.size.title = 12)
  p3 + p4
  dev.off()
  #ggsave("./out/cellchat_net_heatmap_high-low.pdf",width = 8,height = 4,family = "serif")
  
  ####3.信号通路水平的通讯差异 !!!!!!!!!!!!!!!----
  print("3.信号通路水平的通讯差异 !!")
  #3.1 组间富集信号通路差异条形图
  ##基于信息流或互作数对信号通路进行排序
  #pdf("./out/cellchat_rankNet_high-low.pdf",width = 8,height = 4,family = "serif")
  p5 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE,color.use=c("skyblue","pink") ) #堆叠
  p6 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE,color.use=c("skyblue","pink")) #不堆叠
  p5 + p6
  #dev.off()
  ggsave(paste0(folder_path,"/cellchat_rankNet_high-low.pdf"),width = 8,height = 4,family = "serif")
  
  #3.2 传出信号通路水平热图
  i = 1
  pathway.union <- union(object.list[[i]]@netP$pathways,
                         object.list[[i+1]]@netP$pathways)
  pathway.union
  #Reds Blues Purples Greens  GnBu OrRd
  ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]],
                                          pattern = "outgoing", #传出
                                          signaling = pathway.union,
                                          title = names(object.list)[i],
                                          width = 4.5,color.heatmap="Purples",
                                          height = 6.5)
  ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]],
                                          pattern = "outgoing", #传出
                                          signaling = pathway.union,
                                          title = names(object.list)[i+1],
                                          width = 4.5,color.heatmap="Purples",
                                          height = 6.5)
  pdf(paste0(folder_path,"/cellchat_net_outgoing_high-low.pdf"),width = 8,height = 6,family = "serif")
  draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
  dev.off()
  
  #3.3 传入信号通路水平热图
  ht3 = netAnalysis_signalingRole_heatmap(object.list[[i]],
                                          pattern = "incoming", #传入
                                          signaling = pathway.union,
                                          title = names(object.list)[i],
                                          width = 4.5, height = 6.5,
                                          color.heatmap = "Blues")
  ht4 = netAnalysis_signalingRole_heatmap(object.list[[i+1]],
                                          pattern = "incoming", #传入
                                          signaling = pathway.union,
                                          title = names(object.list)[i+1],
                                          width = 4.5, height = 6.5,
                                          color.heatmap = "Blues")
  pdf(paste0(folder_path,"/cellchat_net_incoming_high-low.pdf"),width = 8,height = 6,family = "serif")
  draw(ht3 + ht4, ht_gap = unit(0.5, "cm"))
  dev.off()
  
  #3.4 总体信号通路水平热图
  print("3.4 总体信号通路水平热图")
  ht5 = netAnalysis_signalingRole_heatmap(object.list[[i]],
                                          pattern = "all", #总体
                                          signaling = pathway.union,
                                          title = names(object.list)[i],
                                          width = 4.5, height = 6.5,
                                          color.heatmap = "Reds")
  ht6 = netAnalysis_signalingRole_heatmap(object.list[[i+1]],
                                          pattern = "all", #总体
                                          signaling = pathway.union,
                                          title = names(object.list)[i+1],
                                          width = 4.5, height = 6.5,
                                          color.heatmap = "Reds")
  pdf(paste0(folder_path,"/cellchat_net_all_high-low.pdf"),width = 8,height = 6,family = "serif")
  draw(ht5 + ht6, ht_gap = unit(0.5, "cm"))
  dev.off()
  
  ##scatter_dot
  pdf(paste0(folder_path,"/cellchat_scatter_dot.pdf"),width = 3.5,height = 3,family = "serif")
  netAnalysis_signalingRole_scatter(object.list[[i]])+ggtitle( names(object.list[i]) )
  netAnalysis_signalingRole_scatter(object.list[[i+1]])+ggtitle( names(object.list[i+1]) )
  dev.off()
  
  ####4.配受体对水平通讯差异!!!!!!!!----
  print("4.配受体对水平通讯差异!!!!!!!!")
  #4.1 总配受体对概率差异气泡图
  pdf(paste0(folder_path,"/cellchat_bubble_all_high-low.pdf"),width = 21,height = 10,family = "serif")
  netVisual_bubble(cellchat,angle.x = 45,comparison = c(1, 2) ) #all
  dev.off()
  
  levels(cellchat@idents$joint) #查看细胞亚群
  pdf(paste0(folder_path,"/cellchat_bubble_choose_high-low.pdf"),width = 6,height = 6,family = "serif")
  netVisual_bubble(cellchat,
                   sources.use = number_cancer_cluster, #Malignant cells
                   #targets.use = c(5:12),
                   font.size = 12,font.size.title = 14,
                   comparison = c(1, 2),
                   angle.x = 45)
  dev.off()
  
}

##3.2 cellchat re-plot bar intensity----
if(T){
  #rm(list = ls())
  #setwd("~/R/my_projects/RCC/6.0_cellchat")
  library(CellChat)
  library(ComplexHeatmap)
  library(patchwork)
  high <-readRDS(paste0(folder_path,"/cellchat_High.rds") ) #请加载cellchat对应各分组对象
  low <-readRDS(paste0(folder_path,"/cellchat_Low.rds") )
  high<- updateCellChat(high)
  low<- updateCellChat(low)
  #overall
  #overall <- readRDS(paste0(folder_path,"/cellchat_overall.rds") )
  #overall <- updateCellChat(overall)
  #查看两个数据集细胞分群，并检查是否一致：
  identical(levels(high@idents),levels(low@idents))
  #table(cellchat@idents$joint )
  #levels(cellchat@idents$joint)#
  levels(high@idents)
  #number_cancer_cluster=9 #选择T细胞
}

#prepare data
object.list <- list(high=high,low = low) #overall=overall,
#cellchat <- mergeCellChat(object.list, add.names = names(object.list))

for (i in names(object.list) ){
  print(i)
  #netAnalysis_signalingRole_scatter(high)+netAnalysis_signalingRole_scatter(low)
  p <-netAnalysis_signalingRole_scatter(object.list[[i]])
  df <-data.frame(
    outgoing=p$data$x,incoming=p$data$y,celltype=p$data$labels,count=p$data$Count
  )
  
  ggplot(df,aes(x="Count",y=celltype) )+
    geom_point(aes(size=count,fill=count,color=count), shape = 21) + 
    #scale_size_area(max_size = 5) +
    #geom_point(aes(x = variable, y = ID,fill = ID, size = value), shape = 21) + 
    labs(x = "", y = "") + 
    theme_classic(base_size = 12)+
    scale_color_continuous(low = "skyblue", high = "pink")+
    scale_fill_continuous(low = "skyblue", high = "pink")+
    guides(fill = "none",size = guide_legend(ncol = 1)) 
  ggsave(paste0("./out/Cellchat_out_in_dot_",i,".pdf"),width = 3,height = 2.5,family="serif")
  
  df <-data.frame(
    outgoing=p$data$x,incoming=p$data$y,celltype=p$data$labels#,count=p$data$Count
  )
  
  library(reshape2)
  df_1 <- melt(df, id.vars = "celltype", variable.name = "variable", value.name = "value")
  ggplot(df_1)+
    geom_bar(aes(x = celltype,y = ifelse(variable == 'outgoing',-value,value),fill = factor(variable)),stat = 'identity',position = 'stack')+
    ylab('Strength')+xlab("")+
    labs(fill = 'level')+scale_fill_manual(values = c('#ffc60b','#51c2d5'))+
    #geom_hline(yintercept = 0)+theme(axis.line.x = element_blank())+
    theme_classic(base_size = 12)+
    coord_flip()
  #remove "-" in X axis
  ggsave(paste0("./out/Cellchat_out_in_bar_",i,".pdf"),width = 4,height = 2.5,family="serif")
}

save_pheatmap_pdf <- function(x, filename, width=width, height=height) {
  library(grid)
  #x$gtable$grobs[[1]]$gp <- gpar(lwd = 0.05)#修改聚类树线条宽度
  #x$gtable$grobs[[2]]$vp <- vpar(lwd = 0.05)
  
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height,family = "serif")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}#保存热图的函数

df_1 <-data.frame()
for (i in names(object.list) ){
  print(i)
  #netAnalysis_signalingRole_scatter(high)+netAnalysis_signalingRole_scatter(low)
  p <-netAnalysis_signalingRole_scatter(object.list[[i]])
  df <-data.frame(
    outgoing=p$data$x,incoming=p$data$y,celltype=p$data$labels,count=p$data$Count,group=i
  )
  df_1 <-rbind(df_1,df)
  rownames(df) <- df$celltype
  ##plot
  for (k in c("RdPu","PRGn","PiYG","RdYlGn","Spectral") ){
    print(k)
    p <- pheatmap::pheatmap(df[1:2],cluster_rows = F,cluster_cols = F,scale = "column",
                            border_color = NA,gaps_col =1,cellwidth =  15,cellheight = 15,
                            color=#colorRampPalette( c("skyblue", "white", "pink"))(10)#c("navy", "white", "firebrick3")
                              alpha(rev(RColorBrewer::brewer.pal(9,k) ), 0.7)
    )
    print(p)
    save_pheatmap_pdf(p, paste0("./out/Cellchat_pheatmap_out_in_",i,"_",k,".pdf"),width=4, height=4)
    
  }  
  
}  
##all communication counts-df_1
ggplot(df_1)+
  geom_bar(aes(x = celltype,y = ifelse(group == 'high',-count,count),fill = factor(group)),stat = 'identity',position = 'stack')+
  ylab('Count')+xlab("")+
  labs(fill = 'Group')+scale_fill_manual(values = c('#ffc60b','#51c2d5'))+
  #geom_hline(yintercept = 0)+theme(axis.line.x = element_blank())+
  theme_classic(base_size = 12)+
  coord_flip()
#remove "-" in X axis using AI software
ggsave(paste0("./out/Cellchat_out_in_bar_all.pdf"),width = 4,height = 2.5,family="serif")

library(reshape2)
df_1_1 <- melt(df_1[-c(4,5)], id.vars = "celltype", variable.name = "variable", value.name = "value")
df_1_2 <- melt(df_1[-c(1,2,4)], id.vars = "celltype", variable.name = "variable", value.name = "value")
df_1_1$group <- rep(df_1_2$value,2)

ggplot(df_1_1)+
  geom_bar(aes(x = celltype,y = ifelse(variable == 'outgoing',-value,value),fill = factor(variable)),stat = 'identity',position = 'stack')+
  ylab('Strength')+xlab("")+
  labs(fill = 'level')+scale_fill_manual(values = c('#ffc60b','#51c2d5'))+
  #geom_hline(yintercept = 0)+theme(axis.line.x = element_blank())+
  theme_classic(base_size = 12)+
  coord_flip()+facet_wrap(~group)+
  theme(
    # 去掉分面框的背景边框
    strip.background = element_blank()#strip.text = element_blank()# 去掉分面标题
  )
#remove "-" in X axis
ggsave(paste0("./out/Cellchat_out_in_bar_high_low.pdf"),width = 5,height = 2.5,family="serif")


#all_group pheatmap
df_3 <-data.frame(matrix(ncol=0,nrow=nrow(df)))
rownames(df_3) <- rownames(df)
for (i in names(object.list) ){
  print(i)
  df_2 <-df_1[df_1$group==i,]
  rownames(df_2)<-df_2$celltype
  df_2 <-df_2[,1:2]
  colnames(df_2) <- paste0(colnames(df_2),"_",i)
  df_3 <- cbind(df_2,df_3)
}

for (k in c("RdPu","PRGn","PiYG","RdYlGn","Spectral") ){
  print(k)
  p <- pheatmap::pheatmap(df_3,cluster_rows = F,cluster_cols = F,scale = "column",
                          border_color = NA,gaps_col =c(2,4),cellwidth =  15,cellheight = 15,
                          color=#colorRampPalette( c("skyblue", "white", "pink"))(10)#c("navy", "white", "firebrick3")
                            alpha(rev(RColorBrewer::brewer.pal(9,k) ), 0.7)
  )
  print(p)
  save_pheatmap_pdf(p, paste0("./out/Cellchat_pheatmap_out_in_all_",k,".pdf"),width=4, height=4)
  
}



#using AI to patch figures:https://www.nature.com/articles/s41588-025-02341-9/figures/1


####5.单个/特定信号通路水平差异可视化
print("5.单个/特定信号通路水平差异可视化")
#pathways.show <- c("VEGF") #选择目标信号通路-根据上文
#批量输出所有通路？
# 查看受体配体库，用于pathway_name和interaction_name参数
unique(cellchat@DB$interaction$pathway_name) 
#pathways<-cellchat@DB$interaction
#write.csv(cellchat@DB$interaction,"./out/cellchat_pathways.csv")
pathways.shows <- intersect(object.list[["low"]]@netP[["pathways"]],object.list[["high"]]@netP[["pathways"]] )
pathways.shows
#unique(cellchat@DB$interaction$pathway_name)#c("VEGF","SPP1") #unique(cellchat@DB$interaction$pathway_name)[1] #140

#for (pathways.show in unique(cellchat@DB$interaction$pathway_name) )
if(T){
  
  for (pathways.show in pathways.shows){
    print(pathways.show)
    weight.max <- getMaxWeight(object.list,
                               slot.name = c("netP"),
                               attribute = pathways.show) #控制不同数据集的边权重
    #paste0(folder_path,
    pdf(paste0(folder_path,"/cellchat_net_pathway_",pathways.show,"_high-low.pdf"),width = 8,height = 4,family = "serif")
    par(mfrow = c(1,2), xpd = TRUE)
    for (i in 1:length(object.list)) {
      netVisual_aggregate(object.list[[i]],
                          signaling = pathways.show,
                          layout = "circle",
                          edge.weight.max = weight.max[1],
                          edge.width.max = 10,#pt.title = 1,
                          signaling.name = paste(pathways.show, names(object.list)[i]))
    }
    dev.off()
    
    #5.2 使用热图：
    #pathways.show <- c("CXCL")
    pdf(paste0(folder_path,"/cellchat_heatmap_pathway_",pathways.show,"_high-low.pdf"),width = 6,height = 5,family = "serif")
    
    par(mfrow = c(1,2), xpd = TRUE)
    ht <- list()
    for (i in 1:length(object.list)) {
      ht[[i]] <- netVisual_heatmap(object.list[[i]],font.size = 12,font.size.title = 14,
                                   signaling = pathways.show,
                                   color.heatmap = "Reds",
                                   title.name = paste(pathways.show, "signaling ",names(object.list)[i]))
    }
    ComplexHeatmap::draw(ht[[1]] + ht[[2]], ht_gap = unit(0.5, "cm"))
    dev.off()
    
    
    
  }
  
  ##dotplot
  # netVisual_bubble(cellchat,sources.use = number_cancer_cluster, 
  #                  signaling = pathways.show,
  #                  remove.isolate = TRUE, max.dataset = 2)
  ##tumor --> T/NK
  pdf(paste0(folder_path,"/cellchat_bubble_high-low_Tumor_vs_T.pdf"),width = 3.5,height = 3.5,family = "serif")
  netVisual_bubble(cellchat,
                   sources.use = c(3), #Malignant cells
                   targets.use = c(9),#T
                   color.heatmap = c( "viridis"),
                   signaling = pathways.shows,#[1:5],#c("VEGF","SPP1"),#pathways.shows[1:5],
                   font.size = 12,font.size.title = 14,
                   comparison = c(1, 2),remove.isolate = TRUE,
                   angle.x = 45)
  dev.off()
}

##show unique pathway in high
pathways.shows <- object.list[["high"]]@netP[["pathways"]] 
pathways.shows
for (pathways.show in pathways.shows ){
  print(pathways.show)
  pdf(paste0(folder_path,"/cellchat_heatmap_pathway_",pathways.show,"_high_low.pdf"),width = 4.5,height = 4,family = "serif")
  p<-netVisual_heatmap(object.list[["high"]],font.size = 12,font.size.title = 14,
                    signaling = pathways.show,
                    color.heatmap = "Greens",#"Reds",
                    title.name = paste(pathways.show, "signaling ",names(object.list)["high"]))
  print(p)
  dev.off()
}


#qsave(sceList,paste0(folder_path,"/sce_Tumor.qs"))#SAVE qs!

#4.0 immune molecules----
#https://www.nature.com/articles/s41591-025-03776-7/figures/2
#immune-stimulating and immune checkpoint molecules
# immune_genes <-c("IL6","IL6R","CD27","CD80","CD86","IFNG","PRF1","GZMB","CXCL9","CXCL10",
#                  "IDO1","ICOS","BTLA","TIGIT","LAG3","CTLA4","HAVCR2","PDCD1","CD274","PDCD1LG2"
#                  )
immune_genes_list<-list(
  Immune_stimulators=c("IL6","IL6R","CD27","CD80","CD86","IFNG","PRF1","GZMB","CXCL9","CXCL10"),
  Immune_checkpoint=c("IDO1","ICOS","BTLA","TIGIT","LAG3","CTLA4","HAVCR2","PDCD1","CD274","PDCD1LG2")
)
for (i in c("RdPu","PRGn","PiYG","RdYlGn","Blues","Greens") ){
  print(i)

}
DotPlot(sce_tumor , features =immune_genes_list,#immune_genes,
        group.by = "group",
        #cols = c("skyblue", "pink")
)+ #Minor_Type orig.ident
  #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
  #scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
  #scale_color_viridis_c(option = "A")+
  scale_color_distiller(palette = i,direction = 1)+ #
  labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
  theme_bw(base_size = 14)+
  theme(panel.grid = element_blank(), axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
#hei1 <-ifelse(i=="orig.ident",height,height/2)
ggsave(paste0(folder_path,"/Dotplot_immune_genes_",i,"-1.pdf"),width = 6,height = 4,family="serif")#6 2
