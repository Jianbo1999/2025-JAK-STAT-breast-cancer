

#~/R/my_projects/JAK_STAT_BC/6.0_scRNA_normal_Epi
if(T){
  rm(list = ls())
  setwd("~/R/my_projects/JAK_STAT_BC/6.0_scRNA_normal_Epi")
  folder_path <- "./out"
  # 检查文件夹是否存在
  if (!dir.exists(folder_path)) {
    # 如果文件夹不存在，则创建它
    dir.create(folder_path)
    print(paste("Folder created at", folder_path))
  } else {
    print(paste("Folder already exists at", folder_path))
  }
  library(ggplot2);library(ggpubr);library(RColorBrewer);library(dplyr);library(tibble)
  my36colors <-c("#0072B5","#E18727","#F37C95","#20854E", '#D6E7A3', '#57C3F3', '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 
  library(Seurat);library(SCP);library(qs)
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
  
  genelist<- list("Jak-STAT_pathway"=genes,"Core_Jak-STAT_pathway"=genes_1)
  input_gene= c("Jak-STAT_pathway","Core_Jak-STAT_pathway") #genes#c(genes,"H3K36me","5mC" )
  project = "Jak_STAT_BC"
}
sce <-qread("./GSE161529-normal-BC-Epi/GSE161529_Epi.qs")

#1. GSE161529-normal-BC-Epi----
##1.1 load data----
sce <-readRDS("./GSE161529-normal-BC-Epi/SeuratObject_NormEpi.rds")
sce <- UpdateSeuratObject(object = sce)
DefaultAssay(sce)<-"RNA"
# 生成0.1至1.5的连续分辨率（步长0.1）
sce <- FindNeighbors(sce, dims = 1:30)  # 基于前30个PCs
resolutions <- seq(0, 0.5, by = 0.05)#seq(0.1, 1.5, by = 0.1)
names(sce@graphs)<-gsub("integrated_","RNA_",names(sce@graphs))
#Provided graph.name not present in Seurat object
sce <- FindClusters(sce, resolution = resolutions, verbose = FALSE)#slow
library(clustree)
tree_base <- clustree(sce, prefix = "RNA_snn_res.") 
tree_base
ggsave("./out/tree.pdf",width = 5,height = 7,family="serif")

sce$seurat_clusters <-sce$RNA_snn_res.0.05

#umap
sce <- RunUMAP(sce, reduction = "pca", dims = 1:30, verbose = FALSE )

for (i in c( colnames(sce@meta.data)[c(1,4,6)],"Minor_Type" ) ){
  # print(DimPlot(sce,group.by = i))
  DimPlot(sce,group.by = i,raster = T)
  ggsave(paste0(folder_path,"/UMAP_",i,".pdf"),width = 4.2,height = 3,family="serif")
  CellDimPlot(
    srt = sce, group.by =i,#label = T,#, "Standardpca_SNN_res.0.8"
    #reduction = "StandardUMAP2D", 
    palette = "Paired",raster = T,pt.size = 2,
    theme_use = "theme_blank"
  )
  #ggsave(paste0(folder_path,"/UMAP_",i,"-SCP.pdf"),width = 4.2,height = 3,family="serif")
  CellDimPlot(
    srt = sce, group.by =i,#label = T,#, "Standardpca_SNN_res.0.8"
    #reduction = "StandardUMAP2D", 
    palcolor = my36colors,raster = T,pt.size = 2,
    theme_use = "theme_blank"
  )
  #ggsave(paste0(folder_path,"/UMAP_",i,"-SCP-mycolors.pdf"),width = 4.2,height = 3,family="serif")
  
  #tsne
  DimPlot(sce,group.by = i,reduction = "tsne",raster = T)
  #ggsave(paste0(folder_path,"/TSNE_",i,".pdf"),width = 4.2,height = 3,family="serif")
  CellDimPlot(
    srt = sce, group.by =i,#label = T,#, "Standardpca_SNN_res.0.8"
    reduction = "tsne", 
    palette = "Paired",raster = T,pt.size = 2,
    theme_use = "theme_blank"
  )
  #ggsave(paste0(folder_path,"/TSNE_",i,"-SCP.pdf"),width = 4.2,height = 3,family="serif")
  
  }

# basal (e.g., KRT5, ACTA2, MYLK, SNAI2), luminal progenitor (TNFRSF11A (RANK), KIT), and mature luminal cells (ESR1, PGR, FOXA1) 
markers <-c( 
            "TNFRSF11A", "KIT",#luminal progenitor
            "ESR1", "PGR", "FOXA1",#mature luminal cell
            "ACTA2", "MYLK", "SNAI2"#basal "KRT5",
            #"CD44","CD24"
)
markers_list <- list(
  Basal =c("ACTA2", "MYLK", "SNAI2"),
  Luminal_progenitor=c("TNFRSF11A", "KIT"),
  Mature_luminal_cell=c("ESR1", "PGR", "FOXA1")
  
)
Idents(sce)<-sce$seurat_clusters

DotPlot(sce,features = markers)+RotatedAxis()
DotPlot(sce , features =markers,#group.by = i,
        #cols = c("skyblue", "pink")
)+ #Minor_Type orig.ident
  #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
  scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
  labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
  theme_bw()+
  theme(#panel.grid = element_blank(), 
    axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
ggsave(paste0(folder_path,"/Dotplot_markers.pdf"),width = 3.5,height = 2,family="serif")
#split
DotPlot(sce , features =markers_list,#group.by = i,
        #cols = c("skyblue", "pink")
)+ #Minor_Type orig.ident
  #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
  scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
  labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
  #theme_classic()+
  theme(#panel.grid = element_blank(), 
    axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
ggsave(paste0(folder_path,"/Dotplot_markers-split.pdf"),width = 5,height = 2.5,family="serif")

#2，3 basal
#1 mature luminal cell
#0 4  luminal progenitor
Basal =c(2,3)
Luminal_progenitor=c(0,4)
Mature_luminal=c(1)

current.cluster.ids <- c(Basal,Luminal_progenitor,Mature_luminal)
new.cluster.ids <- c(rep("Basal",length(Basal)),
                     rep("Luminal_progenitor",length(Luminal_progenitor)),
                     rep("Mature_luminal",length(Mature_luminal))
)

sce@meta.data$Minor_Type <- 
  plyr::mapvalues(x = sce$seurat_clusters, from = current.cluster.ids, to = new.cluster.ids)

table(sce@meta.data$Minor_Type)
# Luminal_progenitor     Mature_luminal              Basal 
# 27758              17112               8846

#re plot UMAP
Idents(sce)<-sce$Minor_Type
DotPlot(sce , features =markers,#group.by = i,
        #cols = c("skyblue", "pink")
)+ #Minor_Type orig.ident
  #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
  scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
  labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
  #theme_classic()+
  theme(#panel.grid = element_blank(), 
    axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
ggsave(paste0(folder_path,"/Dotplot_markers-split-Minor_Type.pdf"),width = 5,height = 2.5,family="serif")
#aver exp
library(scRNAtoolVis)
pdf(paste0(folder_path,"/exp_markers-Minor_Type.pdf"), width = 3, height = 3,family="serif")
averageHeatmap(object = sce,#border = T,
               group.by = "Minor_Type",row_title = "",#border = "grey",
               #htCol = c("#339933", "#FFCC00", "#FF0033"),#热图颜色 #
               column_split = 1:length(unique(sce$Minor_Type)),#row_km = T,
               annoCol=T,myanCol = my36colors[1:length(unique(sce$Minor_Type))], #c("darkgoldenrod2","seagreen","steelblue"),#my36colors,
               markerGene = markers)
dev.off()

##1.2 score and exp plot----
###add score----
library(BiocParallel)
register(MulticoreParam(workers = 24, progressbar = TRUE))#32
#score
sce <- AddModuleScore(object = sce,#sce,
                          features = genelist,#gene,
                          name = names(genelist)#"score" #ctrl = 100,
) 
colnames(sce@meta.data)[19] <- 'Score' #ncol(sce@meta.data)
###score plot----
###样本表达异质性
library(Seurat);library(ggplot2);library(ggpubr)

for (j in c("Score",interest_gene)){
  print(j)
  for (i in c("Minor_Type","group") ){
    print(i)
    library(MySeuratWrappers)
    Idents(sce) <-sce$orig.ident
    VlnPlot(sce, features = c("Score"), group.by  = i,
            stacked=T,pt.size= 0, sort='increasing',
            cols = my36colors, #颜色
            #direction = "horizontal", #水平作图
            x.lab = '', y.lab = '')+NoLegend()+labs(y=j,title="")+theme(axis.text.x = element_text(angle = 90,size = 12))
    width = ifelse(length(unique(sce@meta.data[,i]))/5 >4,length(unique(sce@meta.data[,i]))/5,3.5)
    height = ifelse(nchar(max(as.character(unique(sce@meta.data[,i]))))>10,4,3.5)
    ggsave(paste0(folder_path,"/VlnPlot_Score_",j,"_",i,".pdf"),width = width,height = height,family="serif")
    
    #dotplot
    DotPlot(sce , features =c("Score"),group.by = i,
            #cols = c("skyblue", "pink")
    )+ #Minor_Type orig.ident
      #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
      scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
      labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
      theme_bw()+
      theme(#panel.grid = element_blank(), 
        axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
    ggsave(paste0(folder_path,"/Dotplot_Score_",j,"_",i,".pdf"),width = 3.5,height = height/2,family="serif")
    
    #山脊图
    RidgePlot(sce,sort='increasing',
              features = "Score", group.by  = i,cols=my36colors)
    #width = ifelse(length(unique(sce@meta.data[,i]))/5 >4,length(unique(sce@meta.data[,i]))/5,3.5)
    height = ifelse(nchar(max(as.character(unique(sce@meta.data[,i]))))>10,
                    nchar(max(as.character(unique(sce@meta.data[,i]))))/2,
                    4)
    ggsave(paste0(folder_path,"/RidgePlot_Score_",j,"_",i,".pdf"),width = width*1.5,height = height,family="serif")
    
    
  }
}

###  JKs高低组，分布柱状图----
#sce$patient <-sce$group
sce$group <- ifelse(sce$Score > mean(sce$Score),"High","Low")
table(sce$group)
# High   Low 
# 24774 28942 
library(scRNAtoolVis)
Idents(sce)<-sce$Minor_Type
#names(sce@reductions)[3] <-"UMAP"#"umap"
#umap+count
# pdf("./out/scatterCellPlot.pdf", width = 5, height = 4,family="serif")
# scatterCellPlot(object = sce,rm.axis = F,point.size = 0.5,#cell.id = "Cell_Type",
#                 color = rev(RColorBrewer::brewer.pal(9, 'Paired'))#ggsci::pal_npg()(9)
# )
# dev.off()
####柱状图barplot in every sample ----
cellRatioPlot(object = sce,
              sample.name = "patient",
              celltype.name = "Minor_Type", #col.width  fill.col
              fill.col=brewer.pal(12,"Paired")#[c(1,3,9)]#my36colors#c("darkgoldenrod2","seagreen","steelblue") 
)+
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave("./out/cellRatioPlot_major.pdf",width = 7,height = 3,family="serif")#4.5-3 #74 #44  #6.5/9 5

cellRatioPlot(object = subset(sce,sce$group=="High"),
              sample.name = "patient",
              celltype.name = "Minor_Type", #col.width  fill.col
              fill.col=brewer.pal(12,"Paired")#my36colors#c("darkgoldenrod2","seagreen","steelblue") 
)+
  theme(axis.text.x = element_text(angle=90, hjust=1),legend.position = "none")+
  cellRatioPlot(object = subset(sce,sce$group=="Low"),
                sample.name = "patient",
                celltype.name = "Minor_Type", #col.width  fill.col
                fill.col=brewer.pal(12,"Paired")#my36colors#c("darkgoldenrod2","seagreen","steelblue") 
  )+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        axis.title.y = element_blank(),axis.ticks.y = element_blank(),axis.text.y = element_blank())
ggsave("./out/cellRatioPlot_major_Score_All_high-low.pdf",width = 8,height = 3.5,family="serif")#4.5-3 #74 #44  #6.5/9 5

#mean ratio in high/low
Cellratio <- prop.table(table(Idents(sce), sce$group), margin = 2)  
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
  scale_fill_manual(values = brewer.pal(12,"Paired") )+scale_color_manual(values = brewer.pal(12,"Paired") )
  #coord_flip()+
  #theme(panel.border = element_rect(fill=NA,color="black", linewidth=0.5, linetype="solid"))
#p test?
ggsave("./out/Bar_major_Score_high-low.pdf",width = 3,height = 2.5,family="serif")#4.5-3 #74 #44  #6.5/9 5

#pie


### 查看整体分布Ro/e---
library(Startrac)
#remotes::install_local("sscVis-master")
#remotes::install_local("master.zip")
R_oe <- calTissueDist(sce@meta.data,
                      byPatient = F,
                      colname.cluster = "Minor_Type",
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
  ##
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
}

##1.3 monocle+cytotrace----
###1.3.1 SCP Slingshot----

sce <- RunPAGA(
  srt = sce, group_by = "Minor_Type",
  linear_reduction = "PCA", nonlinear_reduction = "UMAP"
)
#ModuleNotFoundError: No module named 'llvmlite.llvmpy'  0.38
#numba 0.55.2 requires llvmlite<0.39,>=0.38.0rc1, but you have llvmlite 0.31.0 which is incompatible.
#pip uninstall numba #llvmlite 0.31.0 
#pip install numba==0.48.0
#pip install llvmlite==0.31.0
#module 'numpy' has no attribute 'long' #1.21.6
##pip install numpy==1.18.5#降低版本1.20.0！#1.24.4 #==1.20.3  pip uninstall numpy
PAGAPlot(srt = sce, reduction = "UMAP", label = TRUE, label_insitu = TRUE, label_repel = TRUE)

# ##Velocity analysis
# sce <- RunSCVELO(
#   srt = sce, group_by = "Minor_Type",
#   linear_reduction = "PCA", nonlinear_reduction = "UMAP"
# )#AttributeError: module 'numpy' has no attribute 'long'
# VelocityPlot(srt = sce, reduction = "UMAP", group_by = "Minor_Type")

sce <- RunSlingshot(srt = sce, group.by = "Minor_Type", reduction = "umap",show_plot = F)
#show_plot = F RunSlingshot`grobs` must be a single grob or a list of grobs, not a list matrix. UMAP / umap
CellDimPlot(sce, group.by = "Minor_Type", reduction = "umap", lineages = paste0("Lineage", 1) , lineages_span = 0.1)
#Error in `gtable_add_grob()`
FeatureDimPlot(sce, features = paste0("Lineage", 1),raster = T,pt.size = 2,
               reduction = "umap",palette = "Reds",#"YlOrRd",
               theme_use = "theme_blank",label = F,show_stat = F,title = ""
                )
ggsave("./out/FeatureDimPlot_Slingshot_type.pdf",width = 5,height = 5,family="serif")
FeatureDimPlot(sce, features = "Score",raster = T,pt.size = 2,
               reduction = "umap",palette = "Blues",#palcolor = c("skyblue","grey","pink"),#
               label = F,show_stat = F,title = "",
               theme_use = "theme_blank")
ggsave("./out/FeatureDimPlot_Slingshot_score.pdf",width = 5,height = 5,family="serif")

# FeatureDimPlot(sce, features = "Score",raster = T,pt.size = 2,split.by = "group",
#                reduction = "umap",palette = "Blues",#palcolor = c("skyblue","grey","pink"),#
#                label = F,show_stat = F,title = "",
#                theme_use = "theme_blank")

#FeaturePlot(sce, features = "Score",order = F)

#
DynamicPlot(
  srt = sce, lineages = c("Lineage1"), group.by = "Minor_Type",
  features = genes_1,
  compare_lineages = TRUE, compare_features = FALSE
)
#pheatmap
sce <- RunDynamicFeatures(srt = sce, lineages = c("Lineage1"), n_candidates = 200)
ht <- DynamicHeatmap(
  srt = sce, lineages = c("Lineage1"),#, "Lineage2"
  use_fitted = TRUE, n_split = 6, reverse_ht = "Lineage1",
  #species = "Mus_musculus", 
  db = "GO_BP", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE,
  heatmap_palette = "viridis", 
  cell_annotation = "Minor_Type",
  #separate_annotation = list("SubCellType", c("Nnat", "Irx1")), 
  separate_annotation_palette = c("Paired", "Set1"),
  #feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  pseudotime_label = 25, pseudotime_label_color = "red",
  height = 5, width = 2
)
print(ht$plot)

###1.3.2 cytotrace----
library(CytoTRACE)
table(is.na(sce$Minor_Type))
####提取表型文件
table(sce$Minor_Type)
phe <- sce$Minor_Type
phe = as.character(phe)
names(phe) <- rownames(sce@meta.data)
####提取表达矩阵
mat_exp <- as.matrix(sce@assays$RNA@counts)

results <- CytoTRACE(mat = mat_exp)#core ?
plotCytoGenes(results, numOfGenes = 10,outputDir = "./cytotrace_out")
plotCytoTRACE(results, phenotype = phe,gene = "CCR7",outputDir = "./cytotrace_out")

#CytoTRACE2 
#devtools::install_github("digitalcytometry/cytotrace2", subdir = "cytotrace2_r")
# library(CytoTRACE2)
# remotes::install_local("./cytotrace2-main.zip")
cyto_res <-read.table("./cytotrace_outCytoTRACE_plot_table.txt")
#add to meta
sce <- AddMetaData(sce,cyto_res[1])

#plot umap
for (i in c("Score","CytoTRACE",interest_gene)){
  print(i)
  FeatureDimPlot(sce, features = i,raster = T,pt.size = 2,
                 reduction = "umap",palette = "YlOrRd",#"Reds",#
                 theme_use = "theme_blank",label = F,show_stat = F,title = ""
  )
  ggsave(paste0("./out/FeatureDimPlot_",i,".pdf"),width = 5,height = 5,family="serif")
  
}

##Cor计算相关性
x="Score";y="CytoTRACE"
methods <-c("pearson","spearman")
for (j in methods){
  print(j)
  df <-sce@meta.data[c(x,y,"Minor_Type")]
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
  ggsave(file=paste0('./out/Cor_',x,"_",y,"_",j,'.pdf'),width = 3,height = 3,family="serif")
  
}

##Dynamic features
# sce <- RunDynamicFeatures(srt = sce, lineages = c("CytoTRACE"), n_candidates = 200)#slow
# ht <- DynamicHeatmap(
#   srt = sce, lineages = c("CytoTRACE"),
#   use_fitted = TRUE, n_split = 6,# reverse_ht = "Lineage1",
#   #species = "Mus_musculus", 
#   db = "GO_BP", anno_terms = TRUE, anno_keys = TRUE, anno_features = TRUE,
#   heatmap_palette = "viridis", cell_annotation = "Minor_Type",
#   #separate_annotation = list("SubCellType", c("Nnat", "Irx1")), separate_annotation_palette = c("Paired", "Set1"),
#   #feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
#   pseudotime_label = 25, pseudotime_label_color = "red",
#   height = 5, width = 2
# )
# print(ht$plot)
#dot plot/pseudotime 拟时序
DynamicPlot(
  srt = sce, lineages = "CytoTRACE", group.by = "Minor_Type",
  features = c(genes_1,"Score"),#c("Plk1", "Hes1", "Neurod2", "Ghrl", "Gcg", "Ins2"),
  compare_lineages = TRUE, compare_features = FALSE
)#too large add_point = TRUE,
ggsave("./out/Dotplot_pseudotime_CytoTRACE.pdf",width = 12,height = 8,family="serif")
DynamicPlot(
  srt = sce, lineages = "CytoTRACE", group.by = "Minor_Type",
  features = c(genes_1,"Score"),#c("Plk1", "Hes1", "Neurod2", "Ghrl", "Gcg", "Ins2"),
  compare_lineages = TRUE, compare_features = FALSE,add_point = F
) 
ggsave("./out/Dotplot_pseudotime_CytoTRACE_nopoint.pdf",width = 12,height = 6,family="serif")

#plot boxplot
FeatureStatPlot(
  srt = sce, group.by = "Minor_Type", bg.by = "Minor_Type",palcolor = my36colors,bg_palcolor =my36colors,bg_apha =0.1,
  stat.by = c(genes_1,"Score","CytoTRACE"), add_box = TRUE,legend.position = "none",
  comparisons = list(
    c("Luminal_progenitor", "Mature_luminal"),
    c("Mature_luminal", "Basal"),
    c("Luminal_progenitor", "Basal")
  )
)
ggsave("./out/Boxplot_Score_genes_CytoTRACE.pdf",width = 9,height = 9,family="serif")

###1.3.3 monocle 2----
if(T){
  library(Seurat)
  library(monocle)
  library(viridis)
}
DimPlot(sce, label = T, pt.size = 1)
#monocle pieline function
#pay attention---Seuratobject version 5.0.0 
ks_run_Monocle2 <- function(object, #seurat obj or expression matrix (建议数据格式转为matrix,如果数据量大转化为稀疏矩阵as(as.matrix(data), "sparseMatrix"))
                            layer, #used when object is a seurat obj
                            assay, #used when object is a seurat obj
                            lowerDetectionLimit = 0.1, 
                            VARgenesM=c("dispersionTable","seurat","differentialGeneTest"),
                            cellAnno=NULL, 
                            define_root=F,
                            root_state,
                            reverse=NULL
){
  
  if(class(object)[1] == 'Seurat') {
    data <- GetAssayData(object=object, layer=layer, assay=assay)#get expression matrix data
    pd <- new("AnnotatedDataFrame", data = object@meta.data)
    fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
    fd <- new("AnnotatedDataFrame", data = fData)
    
    if(all(data == floor(data))) {
      expressionFamily <- negbinomial.size()
    } else if(any(data < 0)){
      expressionFamily <- uninormal()
    } else {
      expressionFamily <- tobit()
    }
    
    #Creates a new CellDateSet object.
    monocle_cds <- newCellDataSet(data,
                                  phenoData = pd, 
                                  featureData = fd,
                                  lowerDetectionLimit=0.1,
                                  expressionFamily=expressionFamily)
    
  }else{
    print("This fucntions only apply for a seurat obj")
  }
  

  #Estimate size factors and dispersions
  #数据处理
  monocle_cds <- estimateSizeFactors(monocle_cds)#size factor标准化细胞之间的mRNA差异
  monocle_cds <- estimateDispersions(monocle_cds)
  
  #质量控制-filter cells
  monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
  # print(head(fData(monocle_cds)))
  # print(head(pData(monocle_cds)))
  # expressed_genes <- row.names(subset(fData(mouse_monocle), num_cells_expressed >= 10))
  monocle_cds <- monocle_cds[fData(monocle_cds)$num_cells_expressed >= 10, ]
  
  
  #select methods for VariableFeatures
  if(VARgenesM=="dispersionTable"){
    
    disp_table <- dispersionTable(monocle_cds)
    ordering_genes <- subset(disp_table,
                             mean_expression >= 0.1 &
                               dispersion_empirical >= 1.5* dispersion_fit)$gene_id
    
  }
  
  
  if(VARgenesM=="seurat"){
    
    ordering_genes <- VariableFeatures(FindVariableFeatures(object, assay = "RNA"), assay = "RNA")
    
  }
  
  
  if(VARgenesM=="differentialGeneTest"){
    
    diff_test_res <- differentialGeneTest(monocle_cds,fullModelFormulaStr = paste0("~",cellAnno))##~后面是表示对谁做差异分析的变量
    diff_test_res_sig <- diff_test_res[order(diff_test_res$qval,decreasing=F),]
    
    ordering_sce <- diff_test_res_sig[diff_test_res_sig$qval< 0.01,]
    
    if(nrow(ordering_sce)>3000){
      
      ordering_genes <- ordering_sce$gene_short_name[1:3000]
      
    }else{
      
      ordering_genes <- rdering_sce$gene_short_name
    }
    
  }
  
  
  #Marks genes for clustering
  monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes)
  plot_ordering_genes(monocle_cds)
  
  
  #cluster
  monocle_cds <- reduceDimension(monocle_cds, max_components = 2,reduction_method = 'DDRTree')
  
  #order cells
  monocle_cds <- orderCells(monocle_cds, reverse=reverse)
  
  if(define_root){
    monocle_cds <- monocle_cds <- orderCells(monocle_cds,root_state = root_state)
  }
  
  return(monocle_cds)
  
}
##run monocle2
#超过五万会报错Error in sequence.default(n..1, from = seq.int(s.1, length(df), s.1), 
Idents(sce)<-sce$Minor_Type
sce <- subset(sce, downsample = 8000) # 随机抽取1w个细胞
# Luminal_progenitor     Mature_luminal              Basal 
# 27758              17112               8846
table(Idents(sce)) # 查看抽样后各群体的细胞数量
# Luminal_progenitor     Mature_luminal              Basal 
# 8000               8000               8000 
sce_CDS <- ks_run_Monocle2(object = sce,#
                           layer = 'counts',
                           assay = "RNA",
                           VARgenesM="dispersionTable",
                           cellAnno = "Minor_Type")

##
sce_CDS <- qread("./out/sce_CDS.qs")

##Plot monocle2
#可视化
plot_cell_trajectory(sce_CDS,show_branch_points = F,color_by = "Pseudotime",cell_size = 0.5)+
  scale_color_viridis()+theme_blank()
ggsave(paste0(folder_path,"/Monocle2_pseudotime.pdf"),width = 4.2,height = 3,family="serif")

plot_cell_trajectory(sce_CDS,show_branch_points = F,color_by = "State",cell_size = 0.5)+
  scale_color_brewer(palette="Set1")+theme_blank()#scale_color_manual(values=my36colors)#
ggsave(paste0(folder_path,"/Monocle2_State.pdf"),width = 4.2,height = 3,family="serif")

plot_cell_trajectory(sce_CDS,show_branch_points = F,color_by = "Minor_Type",cell_size = 0.5)+
  scale_color_manual(values=my36colors)+theme_blank()
ggsave(paste0(folder_path,"/Monocle2_Minor_Type.pdf"),width = 4.5,height = 3,family="serif")

plot_cell_trajectory(sce_CDS,show_branch_points = F,color_by = "Score",cell_size = 0.5)+
  scale_color_viridis()+theme_blank()
ggsave(paste0(folder_path,"/Monocle2_Score.pdf"),width = 4.2,height = 3,family="serif")

plot_cell_trajectory(sce_CDS,show_branch_points = F,color_by = "Score",cell_size = 0.5)+
  scale_color_viridis(option = "A")+theme_blank()
ggsave(paste0(folder_path,"/Monocle2_Score-1.pdf"),width = 4.2,height = 3,family="serif")

#结合可视化cytroTRACE2结果
plot_cell_trajectory(sce_CDS,show_branch_points = F,color_by = "CytoTRACE",cell_size = 0.5)+
  scale_colour_gradientn(colours = (c(  "#3B0F70FF", "#8C2981FF", "#DE4968FF", "#FE9F6DFF", "#FCFDBFFF")),
                         na.value = "transparent",#"#000004FF",
                         limits=c(0,1),
                         breaks = c(0,1),
                         labels=c("(More diff.)", "(Less diff.)"),
                         name = "Relative\norder \n" ,
                         guide = guide_colorbar(frame.colour = "black", ticks.colour = "black"))+
  theme_blank()+#theme_dr(arrow = grid::arrow(length = unit(0, "inches")))+#坐标轴主题修改
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())
ggsave(paste0(folder_path,"/Monocle2_Umap_cytotrace.pdf"),width = 4.2,height = 3,family="serif")

#pheatmap
pdf(paste0(folder_path,"/monocle2-plot_pseudotime_heatmap-raw.pdf"),width = 4,height = 9,family = "serif")
p <- plot_pseudotime_heatmap(sce_CDS[intersect(rownames(sce_CDS),genes1),], 
                             num_cluster = 4, 
                             show_rownames = T, 
                             return_heatmap = T,#add_annotation_col=pData(sce_CDS)["Pseudotime"],
                             #add_annotation_col= bin[c(2)],#2,3 can't add color in celltypes 
                             hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
p
dev.off()

#pheatmap-visCluster
gene = intersect(rownames(sce_CDS),genes1)
library(ClusterGVis)
df <- plot_pseudotime_heatmap2(sce_CDS[gene,],
                               num_clusters = 4,                               
                               cores = 1)#
#visCluster(object = df,plot.type = "heatmap")
pdf(file = paste0(folder_path,"/monocle2-visCluster-1.pdf"),height = 6,width = 4,onefile = F)
#gene = sample(df$wide.res$gene,20,replace = F)
visCluster(object = df, pseudotime_col = c("skyblue","red"),         
           plot.type = "heatmap",#plot.type = "both",#
           markGenes = intersect(gene,genes_1)  )
dev.off()

#plot point # 横坐标是按照pseudotime排好顺序的。
plot_genes_in_pseudotime(sce_CDS[intersect(rownames(sce_CDS),genes_1),],
                         color_by = "celltype",cell_size = 0.5,
                         ncol= 2 )
ggsave(paste0(folder_path,"/monocle2_point-1.pdf"),width = 7,height = 9,family="serif") #save
plot_genes_in_pseudotime(sce_CDS[intersect(rownames(sce_CDS),genes_1),],
                         color_by = "celltype",cell_size = 0.01,relative_expr=0.5,
                         ncol= 2 )
ggsave(paste0(folder_path,"/monocle2_point-2.pdf"),width = 7,height = 9,family="serif") #save


##拟时差异基因
##Branch分析及热图可视化
if(F){
branch<-length(sce_CDS@auxOrderingData[[sce_CDS@dim_reduce_type]]$branch_points)
outputDir=folder_path
sce_CDS$celltype <- sce_CDS$Minor_Type

  #limit CPU
  library(BiocParallel)
  register(MulticoreParam(workers = 48, progressbar = TRUE))
  if (branch != 0){ 
    for(i in 1:branch){
      print(paste0("Find branch ",i," related genes  need  some  time ",lubridate::now()))
      BEAM_res <-BEAM(sce_CDS, branch_point = i, cores = 4,progenitor_method = 'duplicate')
      BEAM_res <- BEAM_res[order(BEAM_res$qval),]
      BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
      heatmap_gene<-row.names(BEAM_res)[order(BEAM_res$qval)][1:50]
      heat<-plot_genes_branched_heatmap(sce_CDS[heatmap_gene,],
                                        branch_point = i,
                                        num_clusters = 2,
                                        cores = 4,
                                        use_gene_short_name = T,
                                        show_rownames = T,
                                        return_heatmap=T)
      png(paste0(outputDir,'/',"branch",i,"_pseudo_heatmap.png"),w=8,h=8,res=300,units="in")
      grid::grid.newpage()
      grid::grid.draw(heat$ph_res$gtable)
      dev.off()
      pdf(paste0(outputDir,'/',"branch",i,"_pseudo_heatmap.pdf"))
      grid::grid.newpage()
      grid::grid.draw(heat$ph_res$gtable)
      dev.off()
      write.csv(BEAM_res,paste0(outputDir,'/',"branch",i,"_pseudo_related_gene.xls"))
      ## 筛选每个branch 分支上变化显著的3个基因展示：
      branched_genes <- row.names(BEAM_res)[order(BEAM_res$qval)][1:3]   #top3
      p <-plot_genes_branched_pseudotime(sce_CDS[branched_genes,],
                                         branch_point = i,
                                         color_by = 'celltype',#Minor_Type
                                         ncol = 1) +
        theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
      ggsave(p,file=paste0(outputDir,'/',"branch",i,"_genes_branched_pseudotime_celltype",".pdf",sep=""),width = 9, height = 8,limitsize = FALSE,family="serif")
      #interest gene exp
      p <-plot_genes_branched_pseudotime(sce_CDS[intersect(row.names(sce_CDS),genes_1),],
                                         branch_point = i,cell_size = 0.5,
                                         color_by = 'celltype',#Minor_Type
                                         ncol = 2) +
        theme(plot.title = element_text(hjust = 0.5),legend.position = "right")+
        scale_color_manual(values = my36colors)
      print(p)
      ggsave(file=paste0(outputDir,'/',"branch",i,"_genes_branched_pseudotime_celltype_genes_1",".pdf",sep=""),width = 7, height = 9,limitsize = FALSE,family="serif")
      
      p <-plot_genes_branched_pseudotime(sce_CDS[branched_genes,],
                                         branch_point = i,
                                         color_by = "celltype",#'orig.ident',
                                         ncol = 1) +
        theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
      p
      ggsave(p,file=paste0(outputDir,'/',"branch",i,"_genes_branched_pseudotime_orig",".pdf",sep=""),width = 6, height = 8,limitsize = FALSE,family="serif")
      
      genes_2 <- intersect(genes_1,rownames(sce_CDS))
      p <-plot_genes_branched_pseudotime(sce_CDS[genes_2,],
                                         #branch_point = i,
                                         color_by = "celltype",#'orig.ident',
                                         ncol = 2) +
        theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
      p
      ggsave(p,file=paste0(outputDir,'/',"branch",i,"_genes_branched_pseudotime_celltype",".pdf",sep=""),width = 6, height = length(genes_2)*0.75,limitsize = FALSE,family="serif")
      
    }
  }else {
    print("No branch points, no BEAM execution")
    branch = ""
  }
  
  ##Top50差异基因做热图，并分成了3个clusters
  ##disp_table <- dispersionTable(sce_CDS)
  ##ordering_genes <- subset(disp_table,mean_expression >= 0.1 &dispersion_empirical >= 1 * dispersion_fit)
  sce_CDS_expressed_genes <- rownames(fData(sce_CDS))
  diff_test_res <- differentialGeneTest(sce_CDS[sce_CDS_expressed_genes,],fullModelFormulaStr = "~sm.ns(Pseudotime)",cores = 1)
  
  sig_gene_names <- diff_test_res %>% dplyr::filter(qval < 0.01) %>% arrange(qval) %>% row.names()  
  diff_test_res1 <-diff_test_res %>% rownames_to_column("gene")
  write.table(diff_test_res1, file = paste0(outputDir,'/','diff_test_ordering_genes.xls'), row.names = F,col.names = T,quote = F,sep="\t")
  #write.table(ordering_genes,file=paste0(outputDir,'/',"ordering_genes.csv"),row.names = F,col.names = T,quote = F,sep="\t")
  write.table(sig_gene_names, file = paste0(outputDir,'/','ordering_genes_sig_gene_names.all.xls'), row.names = F,col.names = F,quote = F)
  #筛选top50 展示热图
  tryCatch({
    p <-plot_pseudotime_heatmap(sce_CDS[sig_gene_names[1:50],],
                                num_clusters = 3,
                                cores = 4,
                                show_rownames = T,
                                return_heatmap=T)
    ggsave(p,file=paste0(outputDir,'/',"pseudotime_dependent_gene_heatmap_top50.pdf"),width = 6, height = 5,limitsize = FALSE)
    #拆分拟时序分析轴上每个cluster的基因（拆分top50个基因）
    clusters <- cutree(p$tree_row, k = 2)
    clustering <- as.data.frame(clusters)%>%rownames_to_column(var = "gene")
    colnames(clustering) <- c("gene","cluster")
    write.table(clustering,file=paste0(outputDir,'/',"pseudotime_dependent_gene_heatmap_top50_gene_by_cluster.xls"),sep="\t",quote = F,col.names = T,row.names = F)
  }, error = function(e) {
    print("Fewer than 50 significant genes")
  })
  
} #slow 

##extract data
plotdf=pData(sce_CDS)[25:27]
colnames(plotdf)
library(ggridges)
ggplot(plotdf, aes(x=Pseudotime,y=celltype,fill=celltype))+
  geom_density_ridges(scale=1) +
  geom_vline(xintercept = c(5,10),linetype=2)+
  scale_y_discrete("")+
  theme_minimal()+
  theme(
    panel.grid = element_blank()
  )+scale_fill_manual(values = my36colors)


sce <-AddMetaData(sce,meta=plotdf)
##Cor计算相关性
x="Score";y="Pseudotime"#"CytoTRACE"
methods <-c("pearson","spearman")
for (j in methods){
  print(j)
  df <-sce@meta.data[c(x,y,"Minor_Type")]
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
  ggsave(file=paste0(folder_path,"/Cor_",x,"_",y,"_",j,".pdf"),width = 3,height = 3,family="serif")
  
}
#boxplot
FeatureStatPlot(
  srt = sce, group.by = "Minor_Type", bg.by = "Minor_Type",palcolor = my36colors,bg_palcolor =my36colors,bg_apha =0.1,
  stat.by = c("Pseudotime"), add_box = TRUE,legend.position = "none",
  comparisons = list(
    c("Luminal_progenitor", "Mature_luminal"),
    c("Mature_luminal", "Basal"),
    c("Luminal_progenitor", "Basal")
  )
)
ggsave(paste0(folder_path,"/Boxplot_Pseudotime.pdf"),width = 3,height = 3,family="serif")

#umap Pseudotime
for (palette in c("Blues","YlOrRd","Reds") ){
  print(palette)
  FeatureDimPlot(sce, features = "Pseudotime",raster = T,pt.size = 2,
                 reduction = "umap",palette = palette,#palcolor = c("skyblue","grey","pink"),#
                 label = F,show_stat = F,title = "",
                 theme_use = "theme_blank")
  ggsave(paste0(folder_path,"/FeatureDimPlot_monocle_sce_Pseudotime_umap_",palette,".pdf"),width = 5,height = 5,family="serif")
  
}

##change root 更换起点
sce_CDS <-orderCells(sce_CDS,root_state = 2) #State
#sce_CDS <-orderCells(sce_CDS,reverse=T) #reverse
##re-plot

##monocle+ggplot2热图+point curve----
##ggplot重做基因拟时细胞散点图
if(T){
  library(monocle)
  library(tidyverse)
  library(ggridges)
  library(RColorBrewer)
  library(scales)
  library(dplyr)
  library(tidytree)
  library(viridis)
  library(scales);library(pheatmap);library(reshape2)
  
  #sce_CDS <- qread("sce_Epi_CDS.qs")
  Time_genes<- intersect(row.names(sce_CDS),genes )#show genes?
  
  #col anno
  Binner <- function(cds_object){
    df <-pData(cds_object)#data.frame(pData(cds_object[,cells_subset]))
    df <- df[,c("Pseudotime", "State","celltype")]
    df <- df[order(df$Pseudotime, decreasing = F),]
    len <- length(df$Pseudotime)
    bin<-round(len/100)
    #dat <-data.frame()#
    State <- c()
    value <- c()
    Pseudotime <-c()
    Pseudotime_value<-c()
    celltype <-c()#
    for(i in 0:99){
      if(i < 99){
        start <- 1+(bin*i)
        stop <- bin+(bin*i)
        
        value <- median(as.numeric(as.vector(df$State[c(start:stop)])))
        Pseudotime_value <- median(as.numeric(as.vector(df$Pseudotime[c(start:stop)])))
        celltype <-c(celltype,df$celltype[start])#i
        State <- c(State, value)
        Pseudotime <- c(Pseudotime, Pseudotime_value)
        #dat <-cbind(dat,data.frame(Pseudotime=Pseudotime))
      }
      else{
        State <- c(State, value)
        Pseudotime <- c(Pseudotime, Pseudotime_value)
        celltype <-c(celltype,df$celltype[i])#
      }
    }
    return(data.frame(State=State,Pseudotime=Pseudotime,celltype=celltype ) ) #as.data.frame(State)
  }
  bin <- Binner(sce_CDS)
  
  pdf(paste0(folder_path,"/monocle2-plot_pseudotime_heatmap.pdf"),width = 4,height = 9,family = "serif")
  p <- plot_pseudotime_heatmap(sce_CDS[Time_genes,], 
                               num_cluster = 4, 
                               show_rownames = T, 
                               return_heatmap = T,
                               add_annotation_col= bin[c(2)],#2,3 can't add color in celltypes 
                               hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
  print(p)
  dev.off()
  
  pdf(paste0(folder_path,"/monocle2-plot_pseudotime_heatmap-1.pdf"),width = 4,height = 9,family = "serif")
  p <- plot_pseudotime_heatmap(sce_CDS[Time_genes,], 
                               num_cluster = 4, 
                               show_rownames = T, 
                               return_heatmap = T,add_annotation_col= bin[2],
                               hmcols =colorRampPalette(rev(brewer.pal(9, "PRGn")))(100))
  print(p)
  dev.off()
  
  pdf(paste0(folder_path,"/monocle2-plot_pseudotime_heatmap-2.pdf"),width = 4,height = 8,family = "serif")
  p <- plot_pseudotime_heatmap(sce_CDS[Time_genes,], 
                               num_cluster = 4, 
                               show_rownames = T, 
                               return_heatmap = F,add_annotation_col= bin[2],
                               hmcols = viridis(256))
  print(p)
  dev.off()
  
  #数据提取，pheatmap重现#提取数据
  cds_subset <- sce_CDS
  newdata <- data.frame(Pseudotime = seq(min(pData(cds_subset)$Pseudotime),
                                         max(pData(cds_subset)$Pseudotime),
                                         length.out = 100)) #sample 100
  #newdata <- pData(cds_subset)["Pseudotime"]#data.frame(Pseudotime =pData(cds_subset)$Pseudotime ) 
  m <- genSmoothCurves(cds_subset[Time_genes,], 
                       trend_formula = '~sm.ns(Pseudotime, df=3)',  
                       relative_expr = T, new_data = newdata)
  m=m[!apply(m,1,sum)==0,]
  m <- log10(m+1) #log处理数据
  #数据缩放
  m=m[!apply(m,1,sd)==0,]
  m=Matrix::t(scale(Matrix::t(m),center=TRUE))
  m[m > 3] = 2#热图最大值
  m[m <- 3] = -2#热图最小值，可自行调整
  
  ##排序设置
  callback = function(hc, mat){
    sv = svd(t(mat))$v[,1]
    dend = reorder(as.dendrogram(hc), wts = sv)
    as.hclust(dend)
  }
  
  #先做一个普通热图
  p1 <- pheatmap(m, 
                 useRaster = T,
                 cluster_cols=FALSE, 
                 cluster_rows=T, 
                 show_rownames=T, 
                 show_colnames=F, 
                 clustering_method = "ward.D2",
                 cutree_rows=4,
                 filename=NA,
                 border_color = NA,
                 fontsize_row = 8,
                 color=colorRampPalette(c("navy","white","firebrick3"))(100),
                 clustering_callback = callback)
  
  
  #列注释
  annotation_col = data.frame(
    pseudotime =newdata$Pseudotime# rescale(newdata$Pseudotime, to = c(-1, 1))
  )
  row.names(annotation_col) <- colnames(m)
  #行注释
  annotation_row <- data.frame(Cluster=factor(cutree(p1$tree_row, 4)))
  row.names(annotation_row) <- rownames(m)
  rowcolor <- c("#85B22E","#E29827","#922927",'#57C3F3') 
  names(rowcolor) <- c("1","2","3","4") #类型颜色
  
  #注释颜色修改
  ann_colors <- list(pseudotime=viridis(100),
                     #colorRampPalette(rev(brewer.pal(9, "PRGn")))(100),#
                     Cluster=rowcolor) #颜色设置
  
  
  p2 <- pheatmap(m, 
                 useRaster = T,
                 cluster_cols=FALSE, 
                 cluster_rows=T, 
                 show_rownames=T, 
                 show_colnames=F, 
                 clustering_method = "ward.D2",
                 cutree_rows=4, #cluster
                 filename=NA,
                 border_color = NA,
                 fontsize_row = 8,
                 color=colorRampPalette(c("navy","white","firebrick3"))(100),#viridis(100),
                 annotation_col = annotation_col,
                 annotation_colors=ann_colors,
                 annotation_row = annotation_row,
                 clustering_callback = callback,
                 annotation_names_col = F,
                 annotation_names_row = F,
                 main="Pseudotime")
  
  #数据提取，ggplot重现#丑#数据就是前面提取的数据#但是需要将数据转化为ggplot的长格式
  heat_gg <- m
  heat_gg <- as.data.frame(heat_gg)
  heat_gg <- heat_gg%>% mutate(gene=row.names(.)) %>% melt()#转化为ggplot画图需要的长列表
  
  p3 <- ggplot(heat_gg,aes(x=variable,y=gene,fill=value))+
    geom_raster()+
    scale_fill_gradient2(low="#003366", high="#990033", mid="white",
                         name='Expression')+
    theme(axis.title = element_blank(),
          panel.grid = element_blank(),
          axis.ticks = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_text(colour = 'black',
                                     size = 8))+
    scale_y_discrete(position = "right")
  #层次聚类
  library(ggtree)
  library(tidytree)
  d <- m%>% as.data.frame()
  hc <- hclust(dist(d), method = "ward.D2")
  clus <-cutree(hc, 4)
  d1 = data.frame(label=names(clus),member=factor(clus))
  ph <- ggtree(as.phylo(hc))
  
  #行注释
  cluster <- d1 %>%
    ggplot(aes(x=1, y=label,fill=member))+
    geom_tile() + 
    scale_y_discrete(position="right") +
    theme_minimal()+xlab(NULL) + ylab(NULL) +
    theme(axis.text.y = element_blank(),
          axis.text.x = element_blank())+
    labs(fill = "cluster")+
    scale_fill_manual(values = c("#708090",'#68A180','#F3B1A0', '#D6E7A3'))
  
  
  #列注释
  group <- colnames(m) %>% as.data.frame() %>% 
    mutate(group=newdata$Pseudotime) %>%
    mutate(p="group") %>%
    ggplot(aes(.,y=1,fill=group))+
    geom_tile() + 
    theme_minimal()+xlab(NULL) + ylab(NULL) +
    theme(axis.text = element_blank(),
          panel.grid = element_blank())+
    scale_fill_viridis()+
    #scale_fill_gradientn(colours = c("#85B22E","#E29827",'#57C3F3',"#922927"))+
    scale_x_discrete(expand = c(0,0))+
    # geom_segment(aes(x = 5, y = 1, xend = 95, yend = 1),
    #              arrow = arrow(length = unit(0.1, "inches"),
    #                            type = 'closed'))+
    theme(plot.margin = margin(0,0,0,0))
  
  
  library(aplot)
  p4 <- p3 %>%insert_left(cluster, width = 0.04)%>%
    insert_top(group, height =0.02)%>%
    insert_left(ph,width=0.1)
  p4
  #山脊图#这里拼图各显神通吧，AI修饰最好
  library(ggridges)
  plotdf=pData(sce_CDS)
  # plotdf$celltype <- factor(plotdf$celltype_minor,levels = c("PMN(4)","PMN(1)","PMN(3)","PMN(5)",
  #                                      "PMN(2)","PMN(6)","PMN(0)","PMN(7)"))
  #plotdf$celltype <-plotdf$celltype_minor
  #细胞拟时顺序是按照下面的图来设定的
  # plot_cell_trajectory(sce_CDS, cell_size = 2.2, color_by = "celltype") +
  #   facet_wrap(~celltype, nrow = 2)
  
  p5 <- ggplot(plotdf, aes(x=Pseudotime,y=celltype,fill=celltype,color=celltype))+
    geom_density_ridges(scale=1) +
    scale_y_discrete(position = 'right')+
    theme_minimal(base_size = 16)+
    theme(panel.grid = element_blank(),
          axis.title = element_blank(),
          axis.text = element_blank()#element_text(colour = 'black', size=8)
    )+
    #scale_x_continuous(position = 'top')+
    scale_fill_manual(values = my36colors )+scale_color_manual(values = my36colors )
  p5
  
  # p6 <- ggplot(heat_gg,aes(x=variable,y=gene,fill=value))+
  #   geom_raster()+
  #   scale_fill_gradient2(low="#003366", high="#990033", mid="white",
  #                        name='Pseudotime')+
  #   theme(axis.title = element_blank(),
  #         panel.grid = element_blank(),
  #         axis.ticks = element_blank(),
  #         axis.text.x = element_blank(),
  #         axis.text.y = element_text(colour = 'black',
  #                                    size = 8))+
  #   scale_y_discrete(position = "right")
  
  # library(patchwork)
  # p5+p6+plot_layout(ncol = 1,heights=c(1,5),guides = 'collect')
  #手动拼图p4+p5
  p4 #趋势相反？
  ggsave(paste0(folder_path,"/monocle2_Epi_ggplot2-heatmap_celltype.pdf"),width = 5,height = 10,family="serif")
  
  p5 
  ggsave(paste0(folder_path,"/monocle2_Epi_ridges_celltype.pdf"),width = 5.5,height = 3,family="serif")
  
  ##
  genes_2 <-intersect(exprs(sce_CDS)@Dimnames[[1]],genes_1)
  
  genes_exp <- list()
  for(i in 1:length(genes_2)){
    A <- log2(exprs(sce_CDS)[genes_2[i],]+1)
    A <- as.data.frame(A)
    genes_exp[[i]] <- A
  }
  
  gene_exp <- do.call(cbind, genes_exp)
  colnames(gene_exp) <- genes_2
  #将上述几个基因的拟时表达添加到monocle
  pData(sce_CDS) = cbind(pData(sce_CDS), gene_exp)
  
  #提取作图数据，只需要游基因表达和拟时即可
  data <- pData(sce_CDS)
  #data$celltype <-data[,celltype]
  colnames(data)
  #选择需要的列即可，我这里的origin.ident就是分组
  data<-data %>% select("celltype","Pseudotime","CytoTRACE","Score",genes_2)
  
  features <- c("CytoTRACE","Score",genes_2)#c("Anxa1", "Ncf1","Ltf","Camp","Ngp","Chil3","S100a8","S100a9")
  plist <- list()
  #ggplot作图
  for (i in 1:length(features)){
    df <- data[, colnames(data)%in% c("celltype","Pseudotime",features[i])]
    colnames(df) <- c("celltype","Pseudotime",'gene')
    p <- ggplot(df, aes(x=Pseudotime, 
                        y=gene)) + 
      geom_point(aes(color=celltype), size=0.5)+
      geom_smooth(method = "loess",level = 0.95,
                  formula = y~x, color='red',se=T)+
      labs(y=features[i],x="Pseudotime")+
      scale_fill_manual(values = my36colors)+scale_color_manual(values = my36colors)+
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black",size = 0.5),
            axis.title = element_text(size = 14, colour = 'black'),#element_blank(),
            axis.text  = element_text(size = 13, colour = 'black')#,legend.position = 'none'
      )#+ggtitle(features[i])
    p
    plist[[i]] <- p
  }
  CombinePlots(plist, ncol = 2)
  ggsave(paste0(folder_path,"/monocle2_point.pdf"),width = 10,height = 12,family="serif") #save
  
  #no point
  for (i in 1:length(features)){
    df <- data[, colnames(data)%in% c("celltype","Pseudotime",features[i])]
    colnames(df) <- c("celltype","Pseudotime",'gene')
    p <- ggplot(df, aes(x=Pseudotime, 
                        y=gene)) + 
      #geom_point(aes(color=celltype), size=0.5)+
      geom_smooth(method = "loess",level = 0.95,
                  formula = y~x, color='red',se=T)+
      labs(y=features[i],x="Pseudotime")+
      theme(panel.grid.major=element_blank(),
            panel.grid.minor=element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.line = element_line(colour = "black",size = 0.5),
            axis.title = element_text(size = 14, colour = 'black'),#element_blank(),
            axis.text  = element_text(size = 13, colour = 'black'),
            legend.position = 'none')#+ggtitle(features[i])
    p
    plist[[i]] <- p
  }
  CombinePlots(plist, ncol = 2)
  ggsave(paste0(folder_path,"/monocle2_no_point.pdf"),width = 6,height = 12,family="serif") #save
  
}

##差异基因拟时序热图优化----
if(T){
  genes1 <- genes#deg1 %>% arrange(desc(morans_I)) %>% pull(gene_short_name)
  cds_list <- list(sce_CDS)
  
  library(ComplexHeatmap)
  library(RColorBrewer)
  library(tidyverse)
  library(circlize)
  
  ###对矩阵行列重新排序
  library(Biobase);library(Seurat)
  # matrix_ori1 <- exprs(cds_list[[1]])[match(genes1, rownames(rowData(cds_list[[1]]))), order( )]
  # 
  # rownames(cds_list[[1]]$Pseudotime)
  # pseudotime(cds_list[[1]])
  # 
  # pData(cds_list[[1]])["Pseudotime"]
  #pData(sce_CDS)$celltype <-pData(sce_CDS)[,celltype]
  
  ###定义色板
  cols <- my36colors[1:length(unique(pData(sce_CDS)$celltype))]#c("#928ac7", "#e0fbd0", "#bdb9da", "#EEA2AD", "#fcb454", "#f7d96e", "#fae2b2", "#9ae76e")
  names(cols) <- levels(factor(pData(sce_CDS)$celltype))#levels(factor(obj_list[[1]]$sub_type))
  
  ###增加列注释信息
  meta <- pData(sce_CDS) %>% rownames_to_column(var = "CellID")##obj_list[[1]]@meta.data
  #order
  meta <-arrange(meta,Pseudotime)
  
  #exp
  matrix_ori1 <- exprs(cds_list[[1]])[intersect(genes1, rownames(cds_list[[1]])),
                                      meta$CellID ]
  matrix1 <- t(apply(matrix_ori1,1,function(x){smooth.spline(x, df=3)$y}))
  matrix1 <- t(apply(matrix1,1,function(x){(x-mean(x))/sd(x)}))
  rownames(matrix1) <- rownames(matrix_ori1)
  colnames(matrix1) <- colnames(matrix_ori1)
  
  cell_anno <- data.frame(cell = colnames(matrix1)) %>% 
    left_join(., meta, by=c("cell" = "CellID")) %>% 
    column_to_rownames(var = "cell")
  ha1 <- HeatmapAnnotation(Pseudotime=cell_anno$Pseudotime,
                           #Score=cell_anno$Score,
                           celltype = cell_anno$celltype,
                           col = list(celltype= cols))
  p1 <- Heatmap(
    matrix1,
    name = "z-score",
    col = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
    show_row_names = FALSE,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8),
    clustering_method_rows = "ward.D2",
    clustering_method_columns = "ward.D2",
    row_title_rot  = 0,
    cluster_rows = TRUE,
    top_annotation = ha1,
    cluster_row_slices = FALSE,
    cluster_columns = FALSE,
    use_raster = TRUE,
    border = F)
  print(p1)
  #基因表达-颜色异常 注释 #Score=cell_anno$Score,
  
  set.seed(1234)
  km1 <- kmeans(matrix1, 4)
  
  p2 <- Heatmap(
    matrix1,
    name = "z-score",
    col = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
    show_row_names = F,
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 12),
    split = km1$cluster,
    row_title_rot = 90,
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    cluster_columns = FALSE,
    top_annotation = ha1,
    use_raster = TRUE,
    border = TRUE
  )
  print(p2)
  #可以看到Kmeans根据基因的表达pattern将所有差异基因聚为4类，我们将顺序做一个调整，并标出差异基因中的转录因子。这里我们使用SCENIC包中提供的人的转录因子基因列表，未做筛选而是直接取转录因子列表与差异基因的交集。
  
  ###添加行注释-兴趣基因
  library(data.table)
  tfs <- genes_1#fread("/path/to/scenic_resources/hs_hgnc_tfs.txt", header = FALSE)
  #tfs <- unique(tfs$V1)
  
  km_tfs <- rownames(matrix1)[rownames(matrix1) %in% tfs]
  row_anno1 <- rowAnnotation(foo = anno_mark(at = match(km_tfs, rownames(matrix1)),
                                             labels = rownames(matrix1)[match(km_tfs, rownames(matrix1))]),
                             gp=gpar(fontsize = 4))
  
  p3 <- Heatmap(
    matrix1,
    name = "z-score",
    col = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
    show_row_names = FALSE,right_annotation = row_anno1,
    #show_row_names =T,#
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8),
    split = km1$cluster,#factor(km1$cluster, levels = c(2,4,1,3)),
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    cluster_columns = FALSE,
    top_annotation = ha1,
    use_raster = F,
    border = TRUE
  )
  pdf(paste0(folder_path,"/monocle2_Epi_Heatmap_celltype-0.pdf"),width = 5.5,height = 10,family="serif")
  print(p3)#save
  dev.off()
  
  p3 <- Heatmap(
    matrix1,
    name = "z-score",
    col = colorRamp2(seq(from=-2,to=2,length=11),rev(brewer.pal(11, "Spectral"))),
    #show_row_names = FALSE,right_annotation = row_anno1,
    show_row_names =T,#
    show_column_names = FALSE,
    row_names_gp = gpar(fontsize = 8),
    split = km1$cluster,#factor(km1$cluster, levels = c(2,4,1,3)),
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    cluster_columns = FALSE,
    top_annotation = ha1,
    use_raster = F,
    border = TRUE
  )
  pdf(paste0(folder_path,"/monocle2_Epi_Heatmap_celltype-1.pdf"),width = 5.5,height = 10,family="serif")
  print(p3)#save
  dev.off()
  
  #最后我们分别提取每个cluster的基因进行GO term富集，富集到的主要条目与文章一致。
  library(clusterProfiler);library(org.Hs.eg.db)
  model <- "org.Hs.eg.db" #human
  
  deg1 <- names(km1$cluster[km1$cluster=='1'])
  full <- bitr(deg1, fromType = "SYMBOL" ,toType = c("ENSEMBL","ENTREZID") , OrgDb = model)
  go1 <- enrichGO(gene = full$ENTREZID , OrgDb = model, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
  go1 <- as.data.frame(go1)#head(as.data.frame(go4), 5)
  write.csv(go1,paste0(folder_path,"/GO_cluster_1.csv") )
  
  deg2 <- names(km1$cluster[km1$cluster=='2'])
  full <- bitr(deg2, fromType = "SYMBOL" ,toType = c("ENSEMBL","ENTREZID") , OrgDb = model)
  go2 <- enrichGO(gene = full$ENTREZID , OrgDb = model, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
  go2 <- as.data.frame(go2)#head(as.data.frame(go4), 5)
  write.csv(go2,paste0(folder_path,"/GO_cluster_2.csv") )
  
  deg3 <- names(km1$cluster[km1$cluster=='3'])
  full <- bitr(deg3, fromType = "SYMBOL" ,toType = c("ENSEMBL","ENTREZID") , OrgDb = model)
  go3 <- enrichGO(gene = full$ENTREZID , OrgDb = model, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
  go3 <- as.data.frame(go3)#head(as.data.frame(go4), 5)
  write.csv(go3,paste0(folder_path,"/GO_cluster_3.csv") )
  
  deg4 <- names(km1$cluster[km1$cluster=='4'])
  full <- bitr(deg4, fromType = "SYMBOL" ,toType = c("ENSEMBL","ENTREZID") , OrgDb = model)
  go4 <- enrichGO(gene = full$ENTREZID , OrgDb = model, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
  go4 <- as.data.frame(go4)#head(as.data.frame(go4), 5)
  write.csv(go4,paste0(folder_path,"/GO_cluster_4.csv") )
  
  #https://mp.weixin.qq.com/s?__biz=MzkyOTU0MDgwMA==&mid=2247484443&idx=1&sn=9cc0240ee7a1593f207f1f315797174d&chksm=c3cfe099c924323c490945d1ff933fd883b3560d67bc0957c4ad61c466014e060e48f22af2fd&mpshare=1&scene=1&srcid=0809owc1nwRjujdchUIPK1HE&sharer_shareinfo=44f7c232a4691198802cf4c891a8dffa&sharer_shareinfo_first=44f7c232a4691198802cf4c891a8dffa#rd
  #save and patch in AI software
  
}


#save sce_CDS
#qsave(sce_CDS,"./out/sce_CDS.qs")

##1.4 enrichment in high/low----
###1.4.1 irGSEA----
library(irGSEA)

#Idents(sce) <- sce$seurat_annotations #设置亚群，后续选择注释后的数据
#计算评分#slow！
sce<- irGSEA.score(object = sce, assay = "RNA", slot = "data", seeds = 123, 
                   ncores = 8, min.cells = 3, min.feature = 0, custom = F, 
                   geneset = NULL, msigdb = T, species = "Homo sapiens", 
                   category = "H", subcategory = NULL, geneid = "symbol", 
                   method = c(  "ssgsea"), #"AUCell", "UCell", "singscore",#slow
                   aucell.MaxRank = NULL, ucell.MaxRank = NULL, kcdf = 'Gaussian') 
#整合差异基因集
result.dge <- irGSEA.integrate(object = sce, 
                               group.by = "group", 
                               metadata = NULL, col.name = NULL, method = c(  "ssgsea"))#"AUCell", "UCell", "singscore",
#可视化
result.dge[["ssgsea"]][["Name"]] <-gsub("HALLMARK-","",result.dge[["ssgsea"]][["Name"]])
irGSEA.heatmap.plot <- irGSEA.heatmap( object= result.dge,
                                       cluster.color=my36colors[1:length(unique(result.dge[["ssgsea"]][["cluster"]]))],
                                       direction.color=c("skyblue","pink"),#significance.color = c("grey70","pink"),
                                       heatmap.width = 12,
                                       method= "ssgsea", top= 50,show.geneset= NULL)
pdf("./out/irGSEA_HALLMARK_type.pdf",width = 6,height = 8,family="serif")
irGSEA.heatmap.plot
dev.off()


###1.4.2 clusterProfiler ----
#找一下各组的差异基因,构建一下数据
DefaultAssay(sce) <- "RNA"
#Idents(sce)
df  <- FindAllMarkers(sce, only.pos = TRUE,
                      min.pct = 0.25,
                      logfc.threshold = 0.75)
df_sig  <- df[df$p_val_adj < 0.05, ]
#write.csv(df_sig, file = "./out/df_sig.csv")

group <- data.frame(gene=df_sig$gene,#构建分组，转化genesymbol-gene ID。
                    group=df_sig$cluster)
library(clusterProfiler)
Gene_ID <- bitr(df_sig$gene, fromType="SYMBOL", 
                toType="ENTREZID", 
                OrgDb="org.Hs.eg.db")
data  <- merge(Gene_ID,group,by.x='SYMBOL',by.y='gene')

data_GO <- compareCluster(
  ENTREZID~group, 
  data=data, 
  fun="enrichGO", 
  OrgDb="org.Hs.eg.db",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05
)
# data_GO_sim <- simplify(data_GO, 
#                         cutoff=0.7, 
#                         by="p.adjust", 
#                         select_fun=min)
# dotplot(data_GO_sim, showCategory=5,font.size = 8)
data_GO_p <- data_GO@compareClusterResult
write.csv(data_GO_p,"./out/GO_type.csv")#save and choose your interest pathway !

#read data and plot
data_GO_p <- read.csv("data_GO_p.csv", header = T)
library(forcats)
data_GO_p$Description <- as.factor(data_GO_p$Description)
data_GO_p$Description <- fct_inorder(data_GO_p$Description)

#top 5
# data_GO_p <- data.frame(data_GO_p) %>%
#   group_by(group) %>%
#   slice_head(n = 5) %>% #前10个,前5个？
#   arrange(desc(pvalue))

ggplot(data_GO_p, aes(Cluster, Description)) +
  geom_point(aes(fill=p.adjust, size=Count), shape=21)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
        axis.text = element_text(color = 'black', size = 10))+
  scale_fill_gradient(low="skyblue",high="pink")+
  labs(x=NULL,y=NULL)#+
  #coord_flip()
ggsave("./out/GO_type.pdf",width = 5,height = 5,family="serif")

##KEGG
data_GO <- compareCluster(
  ENTREZID~group, 
  data=data, 
  fun="enrichKEGG", 
  #OrgDb="org.Hs.eg.db",
  #ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.1
)
data_KEGG_p <- data_GO@compareClusterResult
write.csv(data_KEGG_p,"./out/data_KEGG_p.csv")#save and choose your interest pathway !


ggplot(data_KEGG_p, aes(Cluster, Description)) +
  geom_point(aes(fill=p.adjust, size=Count), shape=21)+
  theme_bw()+
  theme(axis.text.x=element_text(angle=90,hjust = 1,vjust=0.5),
        axis.text = element_text(color = 'black', size = 10))+
  scale_fill_gradient(low="skyblue",high="pink")+
  labs(x=NULL,y=NULL)#+
#coord_flip()
ggsave("./out/GO_type.pdf",width = 5,height = 5,family="serif")
#save qs!
#qsave(sce,"./GSE161529-normal-BC-Epi/GSE161529_Epi.qs")

###1.4.3 SCP enrichment----
#Differential expression analysis
sce <- RunDEtest(srt = sce, group_by = "Minor_Type", fc.threshold = 1, only.pos = FALSE)
VolcanoPlot(srt = sce, group_by = "Minor_Type")
ggsave("./out/VolcanoPlot_Minor_Type.pdf",width = 14,height = 4.5,family="serif")

DEGs <- sce@tools$DEtest_Minor_Type$AllMarkers_wilcox
DEGs <- DEGs[with(DEGs, avg_log2FC > 1 & p_val < 0.05), ]#p_val_adj
# Annotate features with transcription factors and surface proteins
sce <- AnnotateFeatures(sce, species = "Homo_sapiens",#"Mus_musculus", 
                        db = c("GO_BP", "KEGG")#c("TF") #, "CSPA"
                        )
ht <- FeatureHeatmap(
  srt = sce, group.by = "Minor_Type", features = DEGs$gene, feature_split = DEGs$group1,
  species = "Homo_sapiens",#"Mus_musculus", 
  db = c("GO_BP", "KEGG"), anno_terms = TRUE, #, "WikiPathway"
  #feature_annotation = c("TF", "CSPA"), feature_annotation_palcolor = list(c("gold", "steelblue"), c("forestgreen")),
  height = 5, width = 4
)
print(ht$plot)
ggsave("./out/FeatureHeatmap_Minor_Type.pdf",width = 14,height = 6,family="serif")
#enrich in subtype
sce <- RunEnrichment(
  srt = sce, group_by = "Minor_Type", db = c("GO_BP","KEGG"), #species = "Mus_musculus",
  DE_threshold = "avg_log2FC > log2(1.5) & p_val_adj < 0.05"
)
EnrichmentPlot(
  srt = sce, group_by = "Minor_Type",db = c("GO_BP"),only_sig = T, #group_use = c("Ductal", "Endocrine"),
  plot_type = "bar", topTerm = 10,palette ="Paired"#palcolor =my36colors
)
ggsave("./out/EnrichmentPlot_Minor_Type_GO.pdf",width = 14,height = 8,family="serif")
EnrichmentPlot(
  srt = sce, group_by = "Minor_Type",db = c("KEGG"), only_sig = T, #group_use = c("Ductal", "Endocrine"),
  plot_type = "bar", topTerm = 10,palette ="Paired"#palcolor =my36colors
)
ggsave("./out/EnrichmentPlot_Minor_Type_KEGG.pdf",width = 14,height = 6,family="serif")
#dotplot
EnrichmentPlot(srt = sce, group_by = "Minor_Type",  db = c("GO_BP","KEGG"),only_sig = T,
               #theme_use = ggplot2::theme_classic, theme_args = list(base_size = 10),
               palette="Blues",plot_type = "comparison") #YlOrRd
ggsave("./out/EnrichmentPlot_Minor_Type_dot.pdf",width = 14,height = 6,family="serif")

#High low in basal!

#2.0 Enrichment in subtype high vs low----
##2.1 subtype enrichment and high/low in subtype----
for (i in unique(sce$Minor_Type) ){
  print(i)
  sce1 <-subset(sce,sce$Minor_Type==i)
  sce1$group <- ifelse(sce1$Score > mean(sce1$Score),"High","Low")
  for (j in c("group") ){
    print(j)
    ##run
    sce1 <- RunDEtest(srt = sce1, group_by = j, fc.threshold = 1, only.pos = FALSE)#slow
    VolcanoPlot(srt = sce1, group_by = j)
    wid=ifelse(length(unique(sce@meta.data[,j]) )==3,14,9 )
    ggsave(paste0("./out/VolcanoPlot_subtype_",i,"_",j,".pdf"),width = wid,height = 4.5,family="serif")
    
    #DEGs <- paste0("sce1@tools$DEtest_",j,"$AllMarkers_wilcox")
    DEGs <- sce1@tools[[2]]$AllMarkers_wilcox
    write.csv(DEGs,paste0("./out/DEGs_subtype_",i,"_",j,"_all.csv"))
    DEGs <- DEGs[with(DEGs, avg_log2FC > 1/2 & p_val < 0.05), ]#p_val_adj abs(avg_log2FC)
    write.csv(DEGs,paste0("./out/DEGs_subtype_",i,"_",j,".csv"))
    
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
    ggsave(paste0("./out/FeatureHeatmap_subtype_",i,"_",j,".pdf"),width = 14,height = 6,family="serif")
    
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
      plot_type = "bar", topTerm = 10,palette ="Paired"#palcolor =my36colors
    )
    ggsave(paste0("./out/EnrichmentPlot_subtype_",i,"_",j,"_KEGG.pdf"),width = wid,height = 4,family="serif")
    #dotplot
    EnrichmentPlot(srt = sce1, group_by = j,  db = c("GO_BP","KEGG"),only_sig = T,
                   #theme_use = ggplot2::theme_classic, theme_args = list(base_size = 10),
                   palette="Blues",plot_type = "comparison") #YlOrRd
    ggsave(paste0("./out/EnrichmentPlot_subtype_",i,"_",j,"_KEGG_dot.pdf"),width = 14,height = 6,family="serif")
    
    ##irGSEA
    library(irGSEA)
    Idents(sce1)<-sce@meta.data[,j]
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
    pdf(paste0("./out/irGSEA_HALLMARK_subtype_",i,"_",j,".pdf"),width = 6,height = 8,family="serif")
    print(irGSEA.heatmap.plot)
    dev.off()
    
    
    ##end run
  }
  if(T){
    ##progeny
    library(progeny);library(tidyr)
    ## We compute the Progeny activity scores and add them to our Seurat object
    ## as a new assay called Progeny. 
    sce1 <- progeny(sce1, scale=FALSE, organism="Human", top=500, perm=1, 
                    return_assay = TRUE)
    ## We can now directly apply Seurat functions in our Progeny scores. 
    ## For instance, we scale the pathway activity scores. 
    sce1 <- Seurat::ScaleData(sce1, assay = "progeny") 
    ## We transform Progeny scores into a data frame to better handling the results
    progeny_scores_df <- 
      as.data.frame(t(GetAssayData(sce1, slot = "scale.data", 
                                   assay = "progeny"))) %>%
      rownames_to_column("Cell") %>%
      gather(Pathway, Activity, -Cell) 
    CellsClusters <- data.frame(Cell = names(Idents(sce1)), 
                                CellType = as.character(sce1$group),
                                stringsAsFactors = FALSE)
    ## We match Progeny scores with the cell clusters.
    progeny_scores_df <- inner_join(progeny_scores_df, CellsClusters)
    
    ## We summarize the Progeny scores by cellpopulation
    summarized_progeny_scores <- progeny_scores_df %>% 
      group_by(Pathway, CellType) %>%
      summarise(avg = mean(Activity), std = sd(Activity))
    
    ## We prepare the data for the plot
    summarized_progeny_scores_df <- summarized_progeny_scores %>%
      dplyr::select(-std) %>%   
      spread(Pathway, avg) %>%
      data.frame(row.names = 1, check.names = FALSE, stringsAsFactors = FALSE) 
    
    paletteLength = 100
    myColor = colorRampPalette(c("skyblue", "white","red"))(paletteLength)
    
    progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                          length.out=ceiling(paletteLength/2) + 1),
                      seq(max(summarized_progeny_scores_df)/paletteLength, 
                          max(summarized_progeny_scores_df), 
                          length.out=floor(paletteLength/2)))
    progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
                            fontsize_row = 12, 
                            color=myColor, breaks = progenyBreaks, 
                            main = "PROGENy", angle_col = 45,#gaps_row = 8,
                            treeheight_col = 0,  border_color = "grey70")#"NA"
    pdf(paste0("./out/Progeny_",i,"_group.pdf"),width = 3,height = 5,family="serif")
    print(progeny_hmap)
    dev.off() #add sign. ?
    
    #progeny boxplot
    library(ggplot2);library(ggsci);library(ggpubr)
    colnames(progeny_scores_df)
    ggplot(data = progeny_scores_df,
           aes(x = CellType, y = Activity, fill = CellType))+
      facet_wrap(Pathway~. , scales='free',ncol = 7) +
      scale_fill_manual(values = c("#D14039","#008ECB" ) ) + #c("#008ECB", "#EA921D", "#D14039")
      scale_color_manual(values = c("#D14039","#008ECB" ) ) +
      geom_violin(alpha=0.4, position = position_dodge(width = .75),
                  size=0.1, color="black") + # 边框线黑色 black
      geom_boxplot(notch = TRUE, outlier.size = -1,
                   width=0.5,
                   color="black", lwd=0.1, alpha = 0.7)+ # 背景色透明化
      # geom_point(shape = 21, size=2, 
      #            position = position_jitterdodge(), 
      #            color="NA", alpha=1)+ # 边框线黑色 black
      ylab(expression("Pathway activity")) +
      xlab("")  +
      #ggtitle("TCGA")+
      stat_compare_means(method = "t.test",
                         label.x = 1.5,#label.y = max(pathways_long$value),
                         #comparisons=list(c("ClusterA","ClusterB"),c("ClusterB","ClusterC"),c("ClusterA","ClusterC")),
                         step.increase = 0.05,
                         label="p.signif",
                         vjust = 0.7,
      )+#mytheme1+
      theme(legend.position='none')+theme_bw()
    ggsave(paste0("./out/Progeny_",i,"_group_boxplot.pdf"),width = 12,height = 4.5,family="serif")
    
  }#end progeny
  
}

##2.2 enrichment of DEG in your interest subtype----
input_subtype <-"Basal"
DEG <-read.csv("./out/DEGs_subtype_Basal_group.csv",row.names = 1)

input_gene_enrich <-DEG$gene #80
gene_type<-"ALL" #


#run
###1.0 GO----
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
ggsave(paste0("./out/DEG_",input_subtype,"_",gene_type,"_GO.pdf"),family="serif",width = 4,height = 7)

###2.0 KEGG----
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
pdf(paste0("./out/DEG_",input_subtype,"_",gene_type,"_KEGG_cnetplot.pdf"),family="serif",width = 10,height = 6)
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

ego <- ego[c(2,3,6,13,14,22,25,26,35,36,38,45,49,67,76,81),]
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

p<- ggplot(data = ego,
           aes(x = Count, y = Description, fill = -log10(pvalue)) )+
  scale_fill_viridis_c(alpha = 0.7)+
  #scale_fill_distiller(palette = "RdPu",direction = -1) + #RdPu YlOrRd PiYG
  ##PRGn紫色绿色 #PiYG 紫色绿色  RdYlGn橙绿 
  geom_bar(stat = "identity", width = 0.8) +
  labs(x = "Count", y = "", title = "KEGG") + 
  theme_bw()+ mytheme
p
ggsave(paste0("./out/DEG_",input_subtype,"_",gene_type,"_KEGG_bar.pdf"),family="serif",width = 7,height = 4)
#2.2 dotplot
p<- ggplot(data = ego,
           aes(x = Count, y = Description) )+ #, fill = -log10(pvalue)
  #scale_fill_distiller(palette = "RdPu",direction = -1) + #RdPu YlOrRd PiYG
  ##PRGn紫色绿色 #PiYG 紫色绿色  RdYlGn橙绿 
  #geom_bar(stat = "identity", width = 0.8) +
  geom_point(aes(size=Count,color=-log10(pvalue) ))+
  labs(x = "Count", y = "", title = "KEGG") + 
  theme_bw()+ mytheme+
  scale_color_viridis_c(alpha = 0.7)
p
ggsave(paste0("./out/DEG_",input_subtype,"_",gene_type,"_KEGG_dot.pdf"),family="serif",width = 7,height = 4)

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
ggsave(paste0("./out/DEG_",gene_type,"_KEGG_bar_genes.pdf"),width = 5,height = 6.5,family="serif")#7.5

#input all DEG and run: GSEA progeny Dorothea