
#~/R/my_projects/JAK_STAT_BC/10.0_anti-PD1_Nat.Med

if(T){
  rm(list = ls())
  setwd("~/R/my_projects/JAK_STAT_BC/10.0_anti-PD1_Nat.Med")
  folder_path <- "./out"
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
  # library(Seurat);library(SCP);library(qs);library(dplyr)
  # library(tidyverse)
  # library(forestploter)
  # library(GSVA)
  # library(SummarizedExperiment)
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
sceList <-qs::qread("anti-PD1-Nat.Med.qs")
#"anti-PD1-Nat.Med-Tcells.qs"
#ref: https://www.nature.com/articles/s41591-021-01323-8
#re-load using fastmnn function in SCP package

#0.load data----
if(F){
  sceList <- readRDS("~/R/sc_bc/anti-PD1-Nat.Med/out/sce_all_fastmnn.rds")
  library(cols4all);library(ggplot2) #
  mycol <- c4a('10',9)
  clinical <- openxlsx::read.xlsx("~/R/sc_bc/anti-PD1-Nat.Med/download_data/Supplementary Table.xlsx",sheet = 3)#42-11
  table(clinical$Type)
  # ER-/PR-\r\n/HER2+  ER+/PR-\r\n/HER2-  ER+/PR+/\r\nHER2-  ER+/PR+/\r\nHER2+  
  #   4                  1                 16                  2 
  # TNBC  
  # 19 
  
  ##SCORE
  sceList <- AddModuleScore(object = sceList,#sce,
                            features = genelist,#gene,
                            name = names(genelist )#"score" #ctrl = 100,
  ) #
  colnames(sceList@meta.data)[15:16] <- gsub("\\d","",colnames(sceList@meta.data)[15:16] )
  sceList$Score <-sceList$`Jak-STAT_pathway`
  sceList$group <- ifelse(sceList$Score > mean(sceList$Score),"High","Low")
  table(sceList$group)
  # High    Low 
  # 103999 122636 
  sceList$orig.ident  <-sceList$patient_id
}

#1.0 expression----
DimPlot(sceList, reduction = "umap",raster=T,label = F,group.by = "cellType",
        pt.size = 1.0,cols =my36colors ) +NoAxes() +ggtitle("")
ggsave("./out/Umap_celltype.pdf", width = 4.5, height = 3.5,family="serif")

CellDimPlot(
  srt = sceList, group.by = c("cellType"),#label = T,#, "Standardpca_SNN_res.0.8"
  #reduction = "StandardUMAP2D", 
  palcolor = my36colors,raster = T,pt.size = 1.8,
  theme_use = "theme_blank")
ggsave("./out/Umap_celltype-SCP.pdf", width = 4.5, height = 3.5,family="serif")

for (i in c("Score",genes_1)){
  print(i)

  FeatureDimPlot(sceList, features = i,raster = T,pt.size = 2,
                 reduction = "umap",palette = "Reds",#"YlOrRd",#"Reds",#
                 theme_use = "theme_blank",label = F,show_stat = F,title = ""
  )
  ggsave(paste0(folder_path,"/FeatureDimPlot_",i,".pdf"),width = 4.5,height = 3.5,family="serif")
}

FeaturePlot(sceList, features = 'Score',raster=T,cols = c("skyblue","pink","red"), pt.size = 1)+ NoAxes()# +
  #theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))#"skyblue","grey","pink"
ggsave(paste0(folder_path,"/FeaturePlot_1.pdf"),width = 4.5,height = 3.5,family="serif")

FeaturePlot(sceList, features = 'Score',#Srps #cols =c("grey","pink","red"),#cols = c("skyblue","grey","pink"), 
            pt.size = 2,order = T)+ 
  NoAxes() +scale_colour_gradientn(colours = RColorBrewer::brewer.pal(n =9, name ="Reds"))+ggtitle("") #+NoLegend() Blues Reds rev
#theme(panel.border = element_rect(fill=NA,color="black", size=1, linetype="solid"))+ggtitle("")#
ggsave(paste0(folder_path,"/FeaturePlot_2.pdf"),width = 4.5,height = 3.5,family="serif")


##dotplot
#"BC_type" "timepoint" "expansion"
#DotPlot(sceList , features = "Score",group.by = "cellType" )#,assay = "mnn.reconstructed" #+NoAxes()

for (j in c("Score",genes_1)){ #,interest_gene ,colnames(sceList@meta.data)[16:21] 
  print(j)
  for (i in c("BC_type" ,"timepoint" ) ){ #"orig.ident","subtype",
    print(i)
    ##
    library(MySeuratWrappers)
    Idents(sceList) <-sceList$orig.ident
    VlnPlot(sceList, features = j, group.by  = i,
            stacked=T,pt.size= 0, sort='increasing',
            cols = my36colors, #颜色
            #direction = "horizontal", #水平作图
            x.lab = '', y.lab = '')+NoLegend()+labs(y=j,title="")+theme(axis.text.x = element_text(angle = 0,size = 12,hjust = 0.5))
    width = ifelse(length(unique(sceList@meta.data[,i]))/5 >4,length(unique(sceList@meta.data[,i]))/5,3.5)
    height = ifelse(nchar(max(as.character(unique(sceList@meta.data[,i]))))>10,4,3.5)
    ggsave(paste0(folder_path,"/VlnPlot_Score_",j,"_",i,".pdf"),width = width,height = height,family="serif")
    
    for (k in c("RdPu","PRGn","PiYG","Reds") ){ #"RdYlGn","Blues","Greens"
      print(k)
      #dotplot
      DotPlot(sceList , features =j,group.by = i,
              #cols = c("skyblue", "pink")
      )+ #Minor_Type orig.ident
        #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
        #scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
        #scale_color_viridis_c(option = "A")+
        scale_color_distiller(palette = k,direction = -1) +
        labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
        theme_bw(base_size = 12)#+
        #theme(#panel.grid = element_blank(), 
        #  axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
      hei1 <-ifelse(i=="orig.ident",height,height/2)
      ggsave(paste0(folder_path,"/Dotplot_Score_",j,"_",i,"_",k,".pdf"),width = 3,height = hei1+0.8,family="serif")
      
    }
    
    ##
  }
  
}

#,"expansion"
for (i in c("expansion") ){ #"orig.ident","subtype",
  print(i)
  for (k in c("RdPu","PRGn","PiYG","Reds","RdYlGn","Blues","Greens") ){ #
    print(k)
    #dotplot
    DotPlot(sceList[,sceList@meta.data$expansion !="n/a"] , 
            features =c("Score",genes_1),group.by = i,
            #cols = c("skyblue", "pink")
    )+ #Minor_Type orig.ident
      coord_flip()+ #theme_bw()+ #去除背景，旋转图片
      #scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
      #scale_color_viridis_c(option = "A")+
      scale_color_distiller(palette = k,direction = -1) +
      labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
      theme_bw(base_size = 12)+
      theme(axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5))
    #theme(#panel.grid = element_blank(), 
    #  axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
    #hei1 <-ifelse(i=="orig.ident",height,height/2)
    width <-ifelse(i=="BC_type",3.5,3)
    ggsave(paste0(folder_path,"/Dotplot_Score_core_",i,"_",k,".pdf"),width = width,height = 3.5,family="serif")
    
  }
}

#all dot in one
for (i in c("BC_type" ,"timepoint" ) ){ #"orig.ident","subtype",
  print(i)
  for (k in c("RdPu","PRGn","PiYG","Reds","RdYlGn","Blues","Greens") ){ #
    print(k)
    #dotplot
    DotPlot(sceList , features =c("Score",genes_1),group.by = i,
            #cols = c("skyblue", "pink")
    )+ #Minor_Type orig.ident
      coord_flip()+ #theme_bw()+ #去除背景，旋转图片
      #scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
      #scale_color_viridis_c(option = "A")+
      scale_color_distiller(palette = k,direction = -1) +
      labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
      theme_bw(base_size = 12)+
      theme(axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5))
    #theme(#panel.grid = element_blank(), 
    #  axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
    #hei1 <-ifelse(i=="orig.ident",height,height/2)
    width <-ifelse(i=="BC_type",3.5,3)
    ggsave(paste0(folder_path,"/Dotplot_Score_core_",i,"_",k,".pdf"),width = width,height = 3.5,family="serif")
    
  }
  
}

### 查看整体分布Ro/e---
library(Startrac)
#remotes::install_local("sscVis-master")
#remotes::install_local("master.zip")
R_oe <- calTissueDist(sceList@meta.data,
                      byPatient = F,
                      colname.cluster = "cellType",
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
pdf("./out/cellRatioPlot_major_Score_All_high-low-Roe.pdf",width = 3.5,height = 3,family="serif")#4.5-3 #74 #44  #6.5/9 5
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
pdf("./out/cellRatioPlot_major_Score_All_high-low-Roe-1.pdf",width = 3.5,height = 3,family="serif")#4.5-3 #74 #44  #6.5/9 5
print(p)
dev.off()
## +++, Ro/e > 1;
## ++, 0.8 < Ro/e ≤ 1;
## +, 0.2 ≤ Ro/e ≤ 0.8;
## +/−, 0 < Ro/e < 0.2;
## −, Ro/e = 0

#2.0 Cor----
##cell ratio in sample
Cellratio <- prop.table(table(sceList$cellType ,#
                              sceList$patient_id), margin = 2)#计算各组样本不同细胞群比例
Cellratio <-data.frame(Cellratio)

library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
colnames(cellper)[1]<-"id"

##mean score in smaple
Srps <-data.frame(srps=sceList$Score,sample=sceList$orig.ident )#99301-1 如何计算每个样本比值？
table(Srps$sample)

mean_srps<-data.frame()
for (i in unique(Srps$sample)){
  mean_1 <- data.frame(Srps=mean(Srps[Srps$sample==i,]$srps),id=i)
  mean_srps<-rbind(mean_srps,mean_1)
} #20 sample
mean_srps$Group <-ifelse(mean_srps$Srps > mean(mean_srps$Srps),"High","Low")
colnames(mean_srps)[1] <- "Score"
table(mean_srps$Group)
# High  Low 
# 19   23 
all_merge <- merge(mean_srps,cellper,by="id")  #

##lollipop
##plot Cor 相关性-棒棒糖图lollipop----
y=colnames(all_merge)[4:ncol(all_merge)]#y #[2:30]
x="Score"

outTab <- data.frame()#首先，新建一个空的数据框，用于保存后续分析结果
for (j in y ) {
  print(j)
  corT=cor.test(all_merge[,x],all_merge[,j],method="pearson") # spearman
  cor=corT$estimate
  pValue=corT$p.value
  outVector <- cbind(j,as.numeric(cor),as.numeric(pValue) )
  #print(c(j,cor,pValue))  #全是正相关
  outTab<- rbind(outTab,outVector)
}
colnames(outTab)<-c("Cell","Cor","pvalue")
outTab[2:3]<- as.data.frame(sapply(outTab[2:3],as.numeric )) #转为数值
outTab$'-log10(pvalue)' <- -log10(outTab$pvalue) #-log10(0.05) = 1.30
outTab$Sign <-ifelse(outTab$pvalue <0.05,"*"," ")
#outTab$Cell <-gsub("T_cells_","",outTab$Cell) gsub("c./*_","",outTab$Cell)

p <-ggdotchart(outTab, x = "Cell",y = "Cor",
               color ="-log10(pvalue)",
               sorting = "descending",#排序“ascending”, “descending”, “none”
               add.params = list(color = "lightgray"), #画棒棒
               add = "segments",#画棒棒
               dot.size = "Cor",#-log10(pvalue), #"cor"
               xlab="",ylab="Correlation Coefficient", 
               #label.select = list(criteria = " `pvalue` < 0.05 "),
               #label = T, #小数点太多了
               ggtheme = theme_bw()) + 
  #coord_flip()+  #翻转坐标轴
  scale_color_gradient(low = "skyblue",high="pink")+  #颜色范围orange
  font("x.text", size = 13, vjust = 0.1,angle = 0)+  #X坐标轴
  #font("y.text",size = 13,face = "bold")+  #坐标加粗
  #geom_hline(yintercept=c(-0.4,0.4), linetype="dashed",size=0.1)+ #添加虚线
  #scale_color_continuous(low = "skyblue",high="orange")+
  scale_y_continuous(limits = c(-0.6,1))+ #c(-0.6,1.0) #!!!!optimal!
  geom_hline(yintercept=c(0),size=0.1)+
  geom_text(size=5, aes( label = Sign),vjust=1) +#hjust = 1,
  geom_text(size=3, aes( label = round(Cor,2) ),vjust=-1) +#hjust = 1,
  theme(axis.text.y = element_text(size = 12 ),panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5) )  #,legend.position = c(0.7,0.7)
p
write.csv(outTab,paste0(folder_path,"/outTab_cor_ratio_celltype.csv") )
ggsave(paste0(folder_path,"/cor_ratio_celltype_lollipop.pdf"),family="serif",width = 4,height = 3)#4-3

##T cell proportion
library(ggplot2) #绘图
library(ggsignif);library(ggpubr) #添加统计检验
library(ggdist) #云雨图
ggplot(all_merge, aes(x =Group, y = T_cell*100,fill=Group)) +
  geom_jitter(mapping = aes(color=Group),width = .1, alpha = 0.5,size=0.9,show.legend = FALSE)+ #绘制散点图
  #geom_jitter(size=1.5,show.legend = FALSE,alpha=0.8)+
  geom_boxplot(position = position_nudge(x = 0.2),width=0.3,alpha=0.7, #0.14  0.1
               size=0.1,outlier.size = 0,outlier.alpha =0)+ #绘制箱线图，并通过position设置偏移
  stat_halfeye(mapping = aes(fill=Group),width = 0.2, .width = 0, justification = -1.2, point_colour = NA,alpha=0.2) + #绘制云雨图，并通过position设置偏移
  scale_fill_manual(values = my36colors)+   #映射云雨图和箱线图的颜色
  scale_color_manual(values = my36colors)+  #映射散点的颜色
  #expand_limits(x = c(0.5, 3.8))+ #扩展画板，若显示不全，请根据你的数据范围手动调整或删除此行
  ylim(0,70)+ #控制y轴显示范围，若显示不全，请根据你的数据范围手动调整或删除此行
  xlab("") + #STAT4 Srps#设置X轴标题
  ylab("T cell fraction (%)")+ #ggtitle("A visual case")+  #设置主标题
  stat_compare_means(method = "wilcox",label.x=1.3,label = "p.signif")+#p.formatgeom_signif(comparisons = comb_list,step_increase = .1,map_signif_level = TRUE,vjust = 0.5,hjust= 0)+
  theme_bw(base_size = 14)+theme(legend.position="none",panel.grid=element_blank() )+theme(aspect.ratio = 2/1)#2/1
ggsave(paste0("./out/boxplot_","Type","_Tcell.pdf"),width = 3,height = 3,family="serif")  

#3.0 ROC to NE/E----
#pre Score predict NE/E after post
all_merge$ID <- gsub(".*_","",all_merge$id)
clinical$ID <-gsub(" ","",clinical$ID)
all_merge <-merge(clinical,all_merge,by="ID")  #42
all_merge$treat <-all_merge$`Pre/post-.menopausal.state`
write.csv(all_merge,"all_merge.csv")



df_pre <- all_merge[all_merge$`E/NE` != "ND ",]#40
df_pre <- df_pre[df_pre$treat=="Pre ",]#14
#ROC
library(pROC)
roc_data <- roc(df_pre$`E/NE`,#subset(df_pre_on,df_pre$features.plot=="STAT4")$timepoint_On,#df_pre_1$timepoint, #on pre E/NE
                df_pre$Score,#subset(df_pre_on,df_pre$features.plot=="STAT4")$avg.exp_On,#df_pre_1$pct.exp,
                ci=TRUE) 
roc_data[["ci"]] 
roc_data[["auc"]]  #ALL  0.7054:0.5211-0.8896; pre 0.8:0.5442-1

#batch ROC
df_1 <-data.frame()
for (i in colnames(df_pre)[c(13,15:22)] ){
  print(i)
  roc_res<-roc(df_pre$`E/NE`,
               df_pre[,i],
               ci=TRUE) 
  #print(unique(df_pre_1$features.plot) );
  print(round(roc_res[["ci"]],2) )
  df <-data.frame(name=i,
                  L95=round(roc_res[["ci"]],2)[1],
                  H95=round(roc_res[["ci"]],2)[3],
                  AUC=round(roc_res[["ci"]],2)[2])
  print(df)
  df_1<-rbind(df_1,df)
}
df_1 

#plot
library(tidyverse)
df_auc_all <- arrange(df_1, AUC)#排序
# df_auc_all$group_col <- c( rep("#e7a40e", 2),#"skyblue",#e7a40e"
#                      rep("#1c6891", 7)#, rep("#a59d70", 5), rep("#4f4a30", 3)
# )
df_auc_all$group_col <-ifelse(df_auc_all$name=="Score","skyblue","#1c6891")
df_auc_all$name <-factor(df_auc_all$name,levels = unique(df_auc_all$name),ordered = T)
#srps stat4
#df_auc_all<- subset(df_auc_all,df_auc_all$name %in% c("STAT4","Srps") )


ggplot(df_auc_all  )+ 
  #geom_hline(yintercept = 0, linewidth = 0.3)+
  # 线条：
  #geom_linerange(aes(name, ymin = L95, ymax = H95, color = name), show.legend = F)+
  geom_linerange(aes(name, ymin = 0, ymax = AUC), show.legend = F,color = df_auc_all$group_col)+ #, color = group_col
  geom_text(aes(x=name,y=AUC+0.15,label = AUC),color=df_auc_all$group_col,size =3 )+
  geom_point(aes(name, AUC, size=AUC ),color = df_auc_all$group_col)+ #color = group_col,
  #facet_grid(~group)+       #AUC c(0.5,0.55,0.5,0.55)
  # annotate("rect",
  #          xmin = c(0.5,7.5),  #7.5 8.5 top ssGSEA
  #          xmax = c(7.5,8.5),
  #          ymin = 0.4, ymax = 1, alpha = 0.2, fill = rev(unique(df_1$group_col))) + #背景填充
  #annotate("text", label = df_auc_all$AUC, x = df_auc_all$name, y = df_auc_all$AUC+0.1,size=3.0 )+ #,  colour =df_1$group_col
  
  # annotate("rect",
  #          xmin = c(0.5,6.5,7.5),  #6.5 7.5 top ssGSEA
  #          xmax = c(6.5,7.5,14.5),
  #          ymin = 0.4, ymax = 1, alpha = 0.1, fill =unique(df_auc_all$group_col)) + #rev(unique(df_1$group_col))
  scale_y_continuous(expand = c(0,0),limits = c(0,1.2))+
  scale_size_continuous(range=c(2,3.5)) + #点大小
  xlab("")+
  ylab("AUC")+
  theme_bw()+
  theme(legend.position = "none",panel.grid.minor.y  = element_blank(),
        panel.grid.major.x = element_blank(),axis.text.x = element_text(color =df_auc_all$group_col) 
  )+ #,axis.text.y = element_text(color =df_1$name )
  #coord_flip()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),strip.background.x = element_rect(fill = "white", colour = "white" )) 

#ggsave(paste0("./out/",Sys.Date(),"_pre_STAT4_SRPS_AUC_ENE.pdf"),width = 2.5,height = 2,family="serif")
ggsave(paste0("./out/",Sys.Date(),"_pre_score_cells_AUC_on_ENE.pdf"),width = 3,height = 2,family="serif")
write.csv(df_auc_all,paste0("./out/",Sys.Date(),"_pre_score_cells_AUC_on_ENE.csv") )

#only show immune cells
df_auc_all <- df_auc_all[c(3,4,5,6,7,9),]
ggplot(  df_auc_all
       )+ 
  #geom_hline(yintercept = 0, linewidth = 0.3)+：
  #geom_linerange(aes(name, ymin = L95, ymax = H95, color = name), show.legend = F)+
  geom_linerange(aes(name, ymin = 0, ymax = AUC), show.legend = F,linewidth=0.2,linetype = 2)+ #, color = group_col
  geom_text(aes(x=name,y=AUC+0.11,label = AUC,color=name),#color=df_auc_all$group_col,
            size =3 )+
  geom_point(aes(name, AUC, size=AUC,color=name )#,color = df_auc_all$group_col
             )+ #color = group_col,
  #facet_grid(~group)+       #AUC c(0.5,0.55,0.5,0.55)
  # annotate("rect",
  #          xmin = c(0.5,7.5),  #7.5 8.5 top ssGSEA
  #          xmax = c(7.5,8.5),
  #          ymin = 0.4, ymax = 1, alpha = 0.2, fill = rev(unique(df_1$group_col))) + #背景填充
  #annotate("text", label = df_auc_all$AUC, x = df_auc_all$name, y = df_auc_all$AUC+0.1,size=3.0 )+ #,  colour =df_1$group_col
  
  # annotate("rect",
  #          xmin = c(0.5,6.5,7.5),  #6.5 7.5 top ssGSEA
  #          xmax = c(6.5,7.5,14.5),
  #          ymin = 0.4, ymax = 1, alpha = 0.1, fill =unique(df_auc_all$group_col)) + #rev(unique(df_1$group_col))
  scale_y_continuous(expand = c(0,0),limits = c(0,1.0))+
  scale_size_continuous(range=c(2,3.5)) + #点大小
  xlab("")+ylab("AUC")+
  theme_bw()+
  theme(legend.position = "none",panel.grid.minor.y  = element_blank(),
        panel.grid.major.x = element_blank(),axis.text.x = element_text(color =my36colors ) #df_auc_all$group_col 
  )+ #,axis.text.y = element_text(color =df_1$name )
  #coord_flip()+
  theme(axis.text.x = element_text(angle = 90,vjust = 0.5,hjust = 1),strip.background.x = element_rect(fill = "white", colour = "white" ))+
  scale_color_manual(values = my36colors )+scale_fill_manual(values = my36colors )
  
ggsave(paste0("./out/",Sys.Date(),"_pre_score_cells_AUC_on_ENE-immune_cells.pdf"),width = 2.2,height = 2,family="serif")

#4.0 pheatmap score and cells----

###plot pheatmap
library(pheatmap);library(RColorBrewer)
#样本注释
immu_cibersort_score<-all_merge
all_merge$ID <- paste0("P",all_merge$ID)
rownames(immu_cibersort_score)<-all_merge$ID
immu_cibersort_score$Type <- gsub("\r\n","",immu_cibersort_score$Type)

immu_cibersort_score<- arrange(immu_cibersort_score, Group)#排序
ann_colors=list(#FDR=c('[0,0.05]'='#FDDCA9','(0.05,0.5]'="#DDC9E3"),
  Group=c("High"="pink","Low"="skyblue" ))

save_pheatmap_pdf <- function(x, filename, width, height) {
  library(grid)
  x$gtable$grobs[[1]]$gp <- gpar(lwd = 0.02)#修改聚类树线条宽度
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height,family = "serif")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}#保存热图的函数

pheat<-pheatmap(na.omit(t(immu_cibersort_score[,c(13,15:22)])) , #t(immu_cibersort_score[,1:22])
                show_colnames = F,cluster_rows = T,cluster_cols = F,
                scale = "column",
                #cellwidth=3,cellheight = 1,
                annotation_col =immu_cibersort_score[c("Group")],annotation_colors = ann_colors,
                gaps_col= table(immu_cibersort_score$Group)[[1]] ,border_color = NA, #gaps_row= c(7),
                annotation_names_row = F, annotation_names_col = F,
                color=alpha(rev(RColorBrewer::brewer.pal(10,"Spectral") ), 0.7) #RdYlGn Spectral   #PiYG "BrBG" 褐浅蓝  "PiYG"紫绿
                #colorRampPalette( c("skyblue", "white", "firebrick3"))(10),#
)
save_pheatmap_pdf(pheat,"./out/pheatmap.pdf",5,2.5)

pdf("./out/cibersort_pheatmap.pdf",6,3.5,family = "serif")
print(pheat)
dev.off()

#计算显著性后重新绘制热图!#T检验显著性表格
#exp_int[,i]
High <- rownames(immu_cibersort_score[immu_cibersort_score$Group =="High",])#==
Low <-rownames(immu_cibersort_score[immu_cibersort_score$Group !="High",])#!=

pvalue_interest <-data.frame() #T.test
exp_int<- immu_cibersort_score[c(13,15:22)]#exprSet[rownames(exprSet) %in% gene_more,] #raw FPKM #gene_more 兴趣基因集
for (i in 1:ncol(exp_int) ){ 
  print(i)
  pwilcox <- t.test(exp_int[rownames(exp_int) %in% High,i],exp_int[rownames(exp_int) %in% Low,i]) #样本数 High vs Low
  fc <- mean(exp_int[rownames(exp_int) %in% High,i] )/mean( exp_int[rownames(exp_int) %in% Low,i] )#drug/control  手动输入？ 对照组在前！
  pvalue_interest<-rbind(pvalue_interest,data.frame(pvalue=pwilcox$p.value,Fold=fc,Symbol=colnames(exp_int)[i] ) )
}
pvalue_interest$Sign. <- ifelse(pvalue_interest$pvalue>= 0.05,"ns",
                                ifelse(pvalue_interest$pvalue < 0.001,"***",
                                       ifelse(pvalue_interest$pvalue < 0.01,"**","*")
                                ) 
)
pvalue_interest$FC <-ifelse(pvalue_interest$Fold <=0.5,"Decrease",
                            ifelse(pvalue_interest$Fold >1,"Increase","Unchange") )
rownames(pvalue_interest) <-pvalue_interest$Symbol
#图注
anno_row_p <- data.frame(
  Sign.=pvalue_interest$Sign.,FoldChange=pvalue_interest$FC) #, size = 20, replace = TRUE  adj.P.Val
rownames(anno_row_p) <- rownames(pvalue_interest)
ann_colors_p=list(FoldChange=c('Unchange'='#CCECFF','Decrease'="skyblue",'Increase'="#4a6fe3"),  #'Unchange'='#CAB2D6',
                  Group=c("High"="#E18727","Low"= "#0072B5"),
                  treat=c("Pre "= '#E59CC4', "Post "='#AB3282'),
                  `E/NE`=c("E "= "#3DB54A","NE "= "#A0D082","ND "="white"),
                  Type=c("ER-/PR-/HER2+ "="#8DD3C7", "ER+/PR+/HER2- "= "#BEBADA", "ER+/PR+/HER2+ "= "#FB8072", "ER+/PR-/HER2- "= "#80B1D3","TNBC "="#a50f15"),
                  Sign.=c("ns"="grey","*"="pink",'**'="#D33F6A","***"="red" ) )

#c( "#A0D082","#3DB54A", "#09713A")
#c("#c00000", "#FF9933", "#FFCC33", "#746CB1", "#BD77B2")

#plot
int_p <-pheatmap(t(exp_int),
                 show_colnames = F,cluster_rows = T,cluster_cols = F,
                 scale = "column",treeheight_row = 0,#hide tree
                 cellwidth=10,cellheight = 10,cutree_rows =2,
                 annotation_col =immu_cibersort_score[c("Type","E/NE","treat","Group")],#annotation_colors = ann_colors,
                 gaps_col= table(immu_cibersort_score$Group)[[1]] ,border_color = NA, #gaps_row= c(7), 
                 annotation_names_row = F, annotation_names_col = F,
                 annotation_row = anno_row_p,annotation_colors = ann_colors_p,
                 color=alpha(rev(RColorBrewer::brewer.pal(10,"Spectral") ), 0.8) #PiYG RdYlGn Spectral   #PiYG "BrBG" 褐浅蓝  "PiYG"紫绿
                 #colorRampPalette( c("skyblue", "white", "firebrick3"))(10),#
)
pdf(paste0("./out/pheatmap_sign.pdf"),10,2.5,family = "serif") 
print(int_p)#names(cohort[y])
dev.off()  
pdf(paste0("./out/pheatmap_sign-1.pdf"),10,10,family = "serif") #save legend
print(int_p)#
dev.off() 


##plot boxplot
library(reshape2)
bar_df <-melt(immu_cibersort_score,id.vars=c("Group") ) #exp_int, id.vars=c("Cell_type")
my36colors <-colors#c( '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3', '#57C3F3', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 

library(ggplot2);library(ggpubr)
p_method <-c("wilcox.test","t.test")
for (j in p_method){
  ggplot(bar_df,#
         aes(variable,value,color=Group))+  #age_grou
    # geom_point(alpha=0.7,size=0.2,
    #            position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
    #                                          jitter.height = 0,
    #                                          dodge.width = 0.7))+
    geom_boxplot(alpha=1,width=0.7,fill=NA,
                 position=position_dodge(width=0.8),
                 size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                 outlier.stroke = 0.5)+
    # geom_violin(alpha=0.2,width=0.9,
    #             position=position_dodge(width=0.8),
    #             size=0.25)+
    labs(x="",y="Cell proportion")+
    scale_color_manual(values =c("pink","skyblue") )+ # #blue c("#0066CC","#66CCFF","#9DC3E6")
    theme_bw() +#labs(color="",title = gene,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
    theme(text = element_text(size=16),axis.title  = element_text(size=16), axis.text = element_text(size=14), 
          axis.text.x = element_text(angle = 90,hjust = 1),
          panel.grid = element_blank(),legend.key = element_blank() ,legend.position="top" ) + #
    # stat_compare_means(#aes(group = Group) ,
    # comparisons=list(c("Her2+","Luminal"),c("Her2+","TNBC"),c("Luminal","TNBC") ),
    # label = "p.signif",#"p.format p.signif 9.4e-12# method = "wilcox",
    # show.legend= F,#删除图例中的"a"# label.x=1.5,bracket.size=0.1,#vjust=0.5,
    # #label.y = max(gene_exp_clin$STAT4)-1,# hide.ns = F,size=4)+
    stat_compare_means(method = j,size=4, #"wilcox.test"
                       show.legend= F,label = "p.signif",#p.signif p.format
                       label.x =0.75,label.y =max(bar_df$value)-0.05)#
  ggsave(paste0("./out/cibersort_boxplot_",j,"_",names(cohort[y]),".pdf"),width = 12,height = 6,family="serif")
  
  ##single boxplot
  for (i in as.character(unique(bar_df$variable)) ){
    for (j in p_method[1] ){ #wilcox.test
      ggplot(bar_df[bar_df$variable %in% i,],#
             aes(Group,value,color=Group))+
        geom_point(alpha=0.1,size=0.5,
                   position=position_jitterdodge(jitter.width = 0.75,#0.45  1.2
                                                 jitter.height = 0,
                                                 dodge.width = 0.8))+
        geom_boxplot(alpha=1,width=0.75,fill=NA,
                     position=position_dodge(width=0.8),
                     size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                     outlier.stroke = 0.5)+
        # geom_violin(alpha=0.2,width=0.9,
        #             position=position_dodge(width=0.8),
        #             size=0.25)+
        labs(title = i,x="",y="Cell proportion")+
        scale_color_manual(values =c("#0066CC","#66CCFF","#9DC3E6") )+ # #blue c("#0066CC","#66CCFF","#9DC3E6")c("#EDA065","#66CCFF","#7EC7A7")
        theme_bw() +#labs(color="",title = gene,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
        theme(text = element_text(size=16),axis.title  = element_text(size=14), axis.text = element_text(size=14), 
              axis.text.x = element_text(angle = 90,hjust = 1),
              panel.grid = element_blank(),legend.key = element_blank() ,legend.position="none" ) + #
        stat_compare_means(method = j,size=8, #"wilcox.test"
                           show.legend= F,label = "p.signif",#p.signif p.format
                           label.x =1.5,label.y = max(bar_df[bar_df$variable %in% i,][3]) -max(bar_df[bar_df$variable %in% i,][3])/20)
      #,label.y =max(bar_df[bar_df$variable %in% i,]$value)-0.01)#label.x =0.75,
      ggsave(paste0("./out/cibersort_boxplot_",names(cohort[y]),"_",i,"_",j,".pdf"),width = 2.5,height = 3,family="serif")
      print(i)
    }
  }
  
}

#6.0 enrichment----
#https://www.nature.com/articles/s41591-021-01323-8

#qsave(sceList,"anti-PD1-Nat.Med-Tcells.qs")
#save 

