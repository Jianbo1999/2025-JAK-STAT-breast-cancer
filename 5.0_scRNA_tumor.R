
#~/R/my_projects/JAK_STAT_BC/5.0_scRNA_tumor
if(T){
  rm(list = ls())
  setwd("~/R/my_projects/JAK_STAT_BC/5.0_scRNA_tumor")
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
  
  genelist<- list("Jak-STAT_pathway"=genes,"Core_Jak-STAT_pathway"=genes_1)
  input_gene= c("Jak-STAT_pathway","Core_Jak-STAT_pathway") #genes#c(genes,"H3K36me","5mC" )
  project = "Jak_STAT_BC"
}
#
sceList <-qread("sceList_score.qs")
#sce_CDS <- qread("sce_Epi_CDS.qs") 上皮抽样1000 monocle2

#0.load data----
##css integration; data from GES176078 #https://www.nature.com/articles/s41588-021-00911-1
#sceList <-readRDS("~/R/sc_bc/out_1/mydata_CSS_20230925.rds")
##plot umap
if(F){
  DimPlot(sceList,group.by = "celltype_major",label = TRUE,order = T,
          raster = T,
          cols = my36colors)+labs(title = "") #9个 reduction = "umap",
  ggsave("./out/Umap_celltype_major.pdf", width = 5.5, height = 4,family="serif")
  DimPlot(sceList, group.by = "celltype_subset",label = F)#太多
  DimPlot(sceList, group.by = "celltype_minor",label = F,cols = my36colors)
  CellDimPlot(
    srt = sceList, group.by = c("celltype_major"),#label = T,#, "Standardpca_SNN_res.0.8"
    #reduction = "StandardUMAP2D", 
    palette = "Set3",raster = T,pt.size = 1.5,
    theme_use = "theme_blank"
  )
  ggsave("./out/Umap_celltype_major-SCP-Set3.pdf", width = 5, height = 4,family="serif")
  CellDimPlot(
    srt = sceList, group.by = c("celltype_major"),#label = T,#, "Standardpca_SNN_res.0.8"
    #reduction = "StandardUMAP2D", 
    palcolor = my36colors,raster = T,pt.size = 1.5,
    theme_use = "theme_blank"
  )
  ggsave("./out/Umap_celltype_major-SCP.pdf", width = 5, height = 4,family="serif")
  
  CellDimPlot(
    srt = sceList, group.by = c("celltype_major"),label = T,#, "Standardpca_SNN_res.0.8"
    palcolor = my36colors,raster = T,pt.size = 1.5,
    label.size = 3,label.fg ="black",label.bg = NA,
    reduction = "StandardUMAP2D", theme_use = "theme_blank"
  )
  ggsave("./out/Umap_celltype_major-SCP-label.pdf", width = 5, height = 4,family="serif")
  #celltype_minor
  CellDimPlot(
    srt = sceList, group.by = c("celltype_minor"),#label = T,#, "Standardpca_SNN_res.0.8"
    #reduction = "StandardUMAP2D", 
    palcolor = my36colors,raster = T,pt.size = 1.5,
    theme_use = "theme_blank"
  )
  ggsave("./out/Umap_celltype_minor-SCP.pdf", width = 10, height = 4,family="serif")
  CellDimPlot(
    srt = sceList, group.by = c("celltype_minor"),#label = T, error
    palcolor = my36colors,raster = T,pt.size = 1.5,
    label.size = 3,label.fg ="black",label.bg = NA,
    reduction = "StandardUMAP2D", theme_use = "theme_blank"
  )
  #ggsave("./out/Umap_celltype_minor-SCP-label.pdf", width = 10, height = 4,family="serif")
  
  ##plot QC
  for (i in c("nCount_RNA","nFeature_RNA", "percent.mt") ){
    print(i)
    library(MySeuratWrappers) 
    Idents(sceList) <-sceList$orig.ident
    VlnPlot(sceList, features = i, 
            stacked=T,pt.size= 0, 
            #cols = my36colors, #颜色 
            #direction = "horizontal", #水平作图 
            x.lab = '', y.lab = '')+NoLegend()+theme(axis.text.x = element_text(angle = 90,size = 12))
    ggsave(paste0(folder_path,"/VlnPlot_",i,".pdf"),width = 6.5,height = 4,family="serif")
    
  }
  ##plot main marker
  
  markers <- c(
    "CD19","CD79A","MS4A1",  #B_cells
    "COL1A1","C1R","DCN", #"PDGFRA", "PDGFRB",#"ACTA2", #CAFs
    "EPCAM","CDH1","KRT18",#Cancer Epithelial Normal Epithelial
    "PECAM1","CD34","VWF",#endothelial cells
    "LYZ","CD14","CD68",#Myloid
    "IGHG1","MZB1","SDC1",#"CD79A",#Plasmablasts
    "ACTA2", "PDGFRB","MCAM",#PVL
    "CD3D","CD3E" ,"CD3G"#T-cells
  )
  #table(sceList$celltype_major)
  DotPlot(sceList, features = markers, group.by = "celltype_major") + #Main_Type
    #coord_flip() +
    scale_color_binned()+labs(title = NULL,x ="",y ="")+theme(axis.text.x = element_text(angle = 90,hjust = 1) )
  ggsave("./out/DotPlot_markers_celltype_major.pdf", width = 7, height = 3,family="serif")
  DotPlot(sceList,group.by = "celltype_major",
          features =markers,#immune_markers, 
          cols = c("white", "red")) +
    #coord_flip() +scale_color_binned()+
    labs(title = NULL,x ="",y ="")+theme(axis.text.x = element_text(angle = 90,hjust = 1) )
  ggsave("./out/DotPlot_markers_celltype_major-1.pdf", width = 7, height = 3,family="serif")
  
  library(scRNAtoolVis)
  Idents(sceList) <-sceList$celltype_major
  pdf(paste0(folder_path,"/exp_markers-celltype_major.pdf"), width = length(unique(sceList$celltype_major))/2, height = length(markers)/5,family="serif")
  averageHeatmap(object = sceList,#border = T,
                 group.by = "celltype_major",#border = "grey",
                 #htCol = c("#339933", "#FFCC00", "#FF0033"),#热图颜色 #
                 column_split = 1:length(unique(sceList$celltype_major)),#row_km = T,
                 annoCol=T,myanCol = my36colors[1:length(unique(sceList$celltype_major))], #c("darkgoldenrod2","seagreen","steelblue"),#my36colors,
                 markerGene = markers)
  dev.off()
  #1.0 score----
  #~/R/my_projects/RCC_lactylation/6.0_single_cell
  library(BiocParallel)
  register(MulticoreParam(workers = 24, progressbar = TRUE))#32
  #score
  sceList <- AddModuleScore(object = sceList,#sce,
                            features = genelist,#gene,
                            name = names(genelist)#"score" #ctrl = 100,
  ) #
  colnames(sceList@meta.data)[15] <- 'Score' #ncol(sceList@meta.data)
  #save
  #qs::qsave(sceList,"sceList_score.qs")
  
}

##score exp----
###样本表达异质性
library(Seurat);library(ggplot2);library(ggpubr)

for (i in c("orig.ident","subtype", "celltype_major","celltype_minor") ){
  print(i)
  library(MySeuratWrappers)
  Idents(sceList) <-sceList$orig.ident
  VlnPlot(sceList, features = "Score", group.by  = i,
          stacked=T,pt.size= 0, sort='increasing',
          cols = my36colors, #颜色
          #direction = "horizontal", #水平作图
          x.lab = '', y.lab = '')+NoLegend()+labs(y="Score",title="")+theme(axis.text.x = element_text(angle = 90,size = 12))
  width = ifelse(length(unique(sceList@meta.data[,i]))/5 >4,length(unique(sceList@meta.data[,i]))/5,3.5)
  height = ifelse(nchar(max(as.character(unique(sceList@meta.data[,i]))))>10,4,3.5)
  ggsave(paste0(folder_path,"/VlnPlot_Score_",i,".pdf"),width = width,height = height,family="serif")
  
  #dotplot
  DotPlot(sceList , features =c("Score"),group.by = i,
          #cols = c("skyblue", "pink")
  )+ #celltype_major orig.ident
    #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
    scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
    labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
    theme_bw()+
    theme(#panel.grid = element_blank(), 
      axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
  ggsave(paste0(folder_path,"/Dotplot_Score_",i,".pdf"),width = 3.5,height = height,family="serif")
  
  #山脊图
  RidgePlot(sceList,sort='increasing',
            features = "Score", group.by  = i,cols=my36colors)
  #width = ifelse(length(unique(sceList@meta.data[,i]))/5 >4,length(unique(sceList@meta.data[,i]))/5,3.5)
  height = ifelse(nchar(max(as.character(unique(sceList@meta.data[,i]))))>10,
                  nchar(max(as.character(unique(sceList@meta.data[,i]))))/2,
                  4)
  ggsave(paste0(folder_path,"/RidgePlot_Score_",i,".pdf"),width = width*1.5,height = height,family="serif")

  
}

#Score exp in cancer epi --> subtype
for (i in c("orig.ident","subtype") ){
  print(i)
  VlnPlot(subset(sceList,celltype_major=="Cancer Epithelial"), features = "Score", group.by  = i,
          stacked=T,pt.size= 0, sort='increasing',
          cols = my36colors, #颜色 
          #direction = "horizontal", #水平作图 
          x.lab = '', y.lab = '')+NoLegend()+labs(y="Score",title="")+theme(axis.text.x = element_text(angle = 90,size = 12))
  width = ifelse(length(unique(sceList@meta.data[,i]))/5 >4,length(unique(sceList@meta.data[,i]))/5,3)
  ggsave(paste0(folder_path,"/VlnPlot_Score_cancer_cell_",i,".pdf"),width = width,height = 3.5,family="serif")
  #山脊图
  RidgePlot(subset(sceList,celltype_major=="Cancer Epithelial"),sort='increasing',
            features = "Score", group.by  = i,cols=my36colors)
  ggsave(paste0(folder_path,"/RidgePlot_Score_cancer_cell_",i,".pdf"),width = 3.5,height = width,family="serif")
  #dotplot
  DotPlot(subset(sceList,celltype_major=="Cancer Epithelial") , features =c("Score"),group.by = i,
          #cols = c("skyblue", "pink")
  )+ #celltype_major orig.ident
    #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
    scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
    labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
    theme_bw()+
    theme(#panel.grid = element_blank(), 
      axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
  ggsave(paste0(folder_path,"/Dotplot_Score_cancer_cell_",i,".pdf"),width = 3.5,height = height,family="serif")
  
}
#umap

##gene exp 
for (j in c(interest_gene)){ #,interest_gene ,colnames(sceList@meta.data)[16:21] 
  print(j)
  for (i in c("celltype_major") ){ #"celltype_subset", "orig.ident","subtype",
    print(i)
    ##
    library(MySeuratWrappers)
    Idents(sceList) <-sceList$orig.ident
    VlnPlot(sceList, features = j, group.by  = i,
            stacked=T,pt.size= 0, sort='increasing',
            cols = my36colors, #颜色
            #direction = "horizontal", #水平作图
            x.lab = '', y.lab = '')+NoLegend()+labs(y=j,title="")+theme(axis.text.x = element_text(angle = 90,size = 12,vjust = 0.5))
    width = ifelse(length(unique(sceList@meta.data[,i]))/5 >4,length(unique(sceList@meta.data[,i]))/5,3.5)
    height = ifelse(nchar(max(as.character(unique(sceList@meta.data[,i]))))>10,4,3.5)
    #ggsave(paste0(folder_path,"/VlnPlot_",j,"_",i,".pdf"),width = width,height = height,family="serif")
    
    for (k in c("RdPu","PRGn","PiYG","RdYlGn","Blues","Greens") ){
      print(k)
      #dotplot
      DotPlot(sceList , features =j,group.by = i#,dot.scale = 5,scale = TRUE
              #cols = c("skyblue", "pink")
      )+ #Minor_Type orig.ident
        #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
        #scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
        #scale_color_viridis_c(option = "A")+
        scale_color_distiller(palette = k,direction = -1) +
        labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
        theme_bw(base_size = 14)+scale_size_continuous(range = c(1, 8))#3 8调整点比例
      #theme(#panel.grid = element_blank(), 
      #  axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
      hei1 <-ifelse(i=="orig.ident",height,height/2)
      ggsave(paste0(folder_path,"/Dotplot_",j,"_",i,"_",k,"-1.pdf"),width = 3.7,height = hei1+0.8,family="serif")
      
    }
    
    ##
  }
  
}

##score exp in subtype boxplot----
unique(sceList@meta.data$celltype_major)
df_epi_score<-sceList@meta.data[sceList@meta.data$celltype_major %in% c("Normal Epithelial","Cancer Epithelial"),]
df_epi_score$Cell_type <-df_epi_score$celltype_major


library(ggplot2);library(EnvStats);library(ggpubr)
#table(all_merge$Efficacy);colnames(all_merge)
ggplot(df_epi_score,
       aes(x = reorder(Cell_type,-Score), y =Score,fill =Cell_type) )+ #,fill =Efficacy y=Srps
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
  scale_color_manual(values = c("#F37C95", "#56B4E9") )+
  scale_fill_manual(values = c( "#F37C95", "#56B4E9") )+
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
  scale_x_discrete(labels= c("Normal","Cancer"))+
  theme_bw(base_size = 14,base_family = "serif")+#facet_grid(~Treatment)+
  theme(axis.title  = element_text(size=14,family = "serif"), axis.text = element_text(size=12),legend.position ="none",  
        panel.grid = element_blank(),legend.key = element_blank(),axis.title.x =element_blank() )#+ 

ggsave("./out/boxplot_score_epi_vs_tumor.pdf",width = 2.5,height = 2.5,family="serif")

#2.0 group:proportion and enrichment analysis----
#cancer epi Score--> high/low subgroups
##2.1  JKs高低组，分布柱状图----
plot(density(sceList$Score))
plot(density(subset(sceList,sceList$celltype_major=="Cancer Epithelial")$Score))
plot(subset(sceList,sceList$celltype_major=="Cancer Epithelial")$Score)

sceList$group <- ifelse(sceList$Score > mean(sceList$Score),"High","Low")
table(sceList$group)
# High   Low 
# 48324 50980 
library(scRNAtoolVis)
Idents(sceList)<-sceList$celltype_major
names(sceList@reductions)[3] <-"UMAP"#"umap"
#umap+count
# pdf("./out/scatterCellPlot.pdf", width = 5, height = 4,family="serif")
# scatterCellPlot(object = sceList,rm.axis = F,point.size = 0.5,#cell.id = "Cell_Type",
#                 color = rev(RColorBrewer::brewer.pal(9, 'Paired'))#ggsci::pal_npg()(9)
# )
# dev.off()
####柱状图barplot in every sample ----
cellRatioPlot(object = sceList,
              sample.name = "orig.ident",
              celltype.name = "celltype_major", #col.width  fill.col
              fill.col=my36colors#c("darkgoldenrod2","seagreen","steelblue") 
              )+
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave("./out/cellRatioPlot_major.pdf",width = 7,height = 3,family="serif")#4.5-3 #74 #44  #6.5/9 5

cellRatioPlot(object = subset(sceList,sceList$group=="High"),
              sample.name = "orig.ident",
              celltype.name = "celltype_major", #col.width  fill.col
              fill.col=my36colors#c("darkgoldenrod2","seagreen","steelblue") 
)+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  cellRatioPlot(object = subset(sceList,sceList$group=="Low"),
                sample.name = "orig.ident",
                celltype.name = "celltype_major", #col.width  fill.col
                fill.col=my36colors#c("darkgoldenrod2","seagreen","steelblue") 
  )+
  theme(axis.text.x = element_text(angle=90, hjust=1))
ggsave("./out/cellRatioPlot_major_Score_All_high-low.pdf",width = 14,height = 3.5,family="serif")#4.5-3 #74 #44  #6.5/9 5

#mean ratio in high/low
Idents(sceList) <-sceList$celltype_major
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
ggsave("./out/Bar_major_Score_high-low.pdf",width = 3,height = 2.5,family="serif")#4.5-3 #74 #44  #6.5/9 5


### 查看整体分布Ro/e---
library(Startrac)
#remotes::install_local("sscVis-master")
#remotes::install_local("master.zip")
R_oe <- calTissueDist(sceListList@meta.data,
                        byPatient = F,
                        colname.cluster = "celltype_major",
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
###miloR抽样差异分析----
if(T){
  library(miloR)
  library(SingleCellExperiment)
  library(scater)
  library(scran)
  library(dplyr)
  library(patchwork)
  library(Seurat)
}
Idents(sceList)=sceList$group
DimPlot(sceList,label = T)
sceList$group<- factor(sceList$group, levels=unique(sceList$group))#样本分组
#sc.epi$orig.ident=as.character(sc.epi$orig.ident)
#miloR输入对象是SingleCellExperiment，所以我们是常用的seurat对象的话，转化一下
sc.epi <- subset(sceList, downsample = 100)#1000 too long
sc.epi$orig.ident=as.character(sc.epi$orig.ident)
sc.epi <- as.SingleCellExperiment(sc.epi)



sc.epi_milo <- miloR::Milo(sc.epi)
reducedDim(sc.epi_milo,"UMAP") <- reducedDim(sc.epi,"UMAP")
#构建KNN-Graph
sc.epi_milo <- miloR::buildGraph(sc.epi_milo, k = 30, d = 50)#slow!

#在图上定义邻域
sc.epi_milo <- makeNhoods(sc.epi_milo, 
                          prop = 0.2, #定义要随机抽样的图顶点的比例，通常为0.1-0.2
                          k = 30, #建议使用与buildGraph一样的k值
                          d=50, #KNN降维数，建议使用与buildGraph一样的d值
                          refined = TRUE)
#量化每个邻域的细胞数量
sc.epi_milo <- countCells(sc.epi_milo, meta.data = data.frame(colData(sc.epi_milo)), 
                          sample="orig.ident")
#sample是指每个样本
#构建分组矩阵
traj_design <- data.frame(colData(sc.epi_milo))[,c("orig.ident", "group")]#分别是重复样本ID和分组
traj_design$orig.ident <- as.factor(traj_design$orig.ident)
traj_design <- distinct(traj_design)
rownames(traj_design) <- traj_design$orig.ident
#计算细胞距离
sc.epi_milo <- calcNhoodDistance(sc.epi_milo, d=50)
#两个分组之间的差异分析
da_results <- testNhoods(sc.epi_milo, 
                         design = ~ group, 
                         design.df = traj_design)
#为可视化建立一个抽象的邻域图
sc.epi_milo <- buildNhoodGraph(sc.epi_milo)
#Plot neighbourhood graph

plotUMAP(sc.epi_milo, colour_by = "celltype.3")

plotNhoodGraphDA(sc.epi_milo, da_results, alpha=0.1) +
  scale_fill_gradient2(low="#070091",#修改颜色
                       mid="lightgrey",
                       high="#910000", 
                       name="log2FC",
                       limits=c(-5,5),
                       oob=squish)
#
da_results <- annotateNhoods(sc.epi_milo, da_results, coldata_col = "celltype.3")
plotDAbeeswarm(da_results, group.by = "celltype.3") +
  scale_color_gradient2(low="#070091",
                        mid="lightgrey",
                        high="#910000",
                        limits=c(-5,5),
                        oob=squish) +
  labs(x="", y="Log2 Fold Change") +
  theme_bw(base_size=10)+
  theme(axis.text = element_text(colour = 'black'))
plotDAbeeswarm(da_results,group.by = "celltype.3", alpha=1)+geom_boxplot(outlier.shape = NA)+
  scale_color_gradient2(midpoint=0, low="blue", mid="white",high="red", space ="Lab" )+
  geom_hline(yintercept=0, linetype="dashed")
ggsave("./out/miloR-All_high-low-Roe.pdf",width = 5,height = 4,family="serif")#4.5-3 #74 #44  #6.5/9 5

##肿瘤Score分组 High-low enrichment


##计算样本表达比例/肿瘤表达比例

#3.0 monocle+cytotrace----
#cancer epi + normal epi
sce <-subset(sceList,celltype_major %in% c("Cancer Epithelial","Normal Epithelial") )
for (i in names(sce@reductions) ){
  print(i)
  sce@reductions[i] <-NULL
}

table(sce$celltype_major)
# Cancer Epithelial Normal Epithelial 
# 24489              4338 
table(sce$celltype_minor)
# Cancer Basal SC      Cancer Cycling      Cancer Her2 SC 
# 4312                5359                3708 
# Cancer LumA SC      Cancer LumB SC Luminal Progenitors 
# 7742                3368                1990 
# Mature Luminal       Myoepithelial 
# 1260                1088 
sce$Minor_Type <- sce$celltype_minor #cell type
##re-umap
library(BiocParallel)
register(MulticoreParam(workers = 48, progressbar = TRUE))
sce <-Standard_SCP(sce,cluster_resolution = 0.1)
##plot uamp Epi
#"StandardUMAP2D"
CellDimPlot(
  srt = sce, group.by ="Minor_Type",#label = T,#, "Standardpca_SNN_res.0.8"
  reduction = "StandardUMAP2D", 
  palette = "Paired",raster = T,pt.size = 2,
  theme_use = "theme_blank"
)
ggsave(paste0(folder_path,"/UMAP_Epi_reUMAP.pdf"),width = 4.2,height = 3,family="serif")

##3.1 CytoTRACE----
#run
if(T){
  library(CytoTRACE)
  ####提取表型文件
  phe <- sce$Minor_Type
  phe = as.character(phe)
  names(phe) <- rownames(sce@meta.data)
  ####提取表达矩阵
  mat_exp <- as.matrix(sce@assays$RNA@counts)
  results <- CytoTRACE(mat = mat_exp,ncores = 16)
  plotCytoGenes(results, numOfGenes = 10,outputDir = "./cytotrace_out")
  plotCytoTRACE(results, phenotype = phe,gene = interest_gene,outputDir = "./cytotrace_out")
  plotCytoTRACE(results, phenotype = phe,gene = interest_gene,
                emb = sce@reductions[["StandardUMAP2D"]]@cell.embeddings ,
                outputDir = "./cytotrace_out-umap")
  
  ####plot
  #boxplot
  cyto_res <-as.data.frame(results$CytoTRACE)#read.table("./cytotrace_outCytoTRACE_plot_table.txt")
  colnames(cyto_res)<-"CytoTRACE"
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
  
  #plot exp violin
  FeatureStatPlot(
    srt = sce, group.by = "Minor_Type", bg.by = "Minor_Type",palcolor = my36colors,bg_palcolor =my36colors,bg_apha =0.1,
    stat.by = c(genes_1,"Score","CytoTRACE"), add_box = TRUE,legend.position = "none" #,
    # comparisons = list(
    #   c("Luminal_progenitor", "Mature_luminal"),
    #   c("Mature_luminal", "Basal"),
    #   c("Luminal_progenitor", "Basal")
    # )
  )
  ggsave("./out/Boxplot_Score_genes_CytoTRACE.pdf",width = 12,height = 9,family="serif")
  #violin
  library(MySeuratWrappers)
  for (i in c("CytoTRACE",interest_gene) ){
    print(i)
    VlnPlot(sce, features = i, group.by  = "Minor_Type",
            stacked=T,pt.size= 0, #sort='increasing',
            cols = my36colors #颜色
            #direction = "horizontal", #水平作图 x.lab = '', y.lab = 'CytoTRACE'
    )+NoLegend()+labs(y=i,x="",title="")+theme(axis.text.x = element_text(angle = 90,size = 12))
    ggsave(paste0(folder_path,"/VlnPlot_Epi_vs_cancer_",i,".pdf"),width = 4,height = 3.5,family="serif")
    
  }
   
}
#save Epi 
#qsave(sce,"sce_Epi.qs")
##3.2 monocle----
if(T){
  library(Seurat)
  library(monocle)
  library(viridis)
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
    # expressed_genes <- row.names(subset(fData(sce_CDS), num_cells_expressed >= 10))
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
}

##run monocle2
#超过五万会报错Error in sequence.default(n..1, from = seq.int(s.1, length(df), s.1), 
Idents(sce)<-sce$Minor_Type
table(sce$Minor_Type)
sce <- subset(sce, downsample = 1000)
#run
sce_CDS <- ks_run_Monocle2(object = sce,#
                           layer = 'counts',
                           assay = "RNA",
                           VARgenesM="dispersionTable",
                           cellAnno = "Minor_Type")

##Plot monocle2 #可视化
plot_cell_trajectory(sce_CDS,show_branch_points = F,color_by = "Pseudotime",cell_size = 0.5)+
  scale_color_viridis()+theme_blank()
ggsave("./out/Monocle2_pseudotime.pdf",width = 4.2,height = 3,family="serif")

plot_cell_trajectory(sce_CDS,show_branch_points = F,color_by = "State",cell_size = 0.5)+
  scale_color_brewer(palette="Set1")+theme_blank()#scale_color_manual(values=my36colors)#
ggsave("./out/Monocle2_State.pdf",width = 4.2,height = 3,family="serif")

plot_cell_trajectory(sce_CDS,show_branch_points = F,color_by = "Minor_Type",cell_size = 0.5)+
  scale_color_manual(values=my36colors)+theme_blank()
ggsave("./out/Monocle2_Minor_Type.pdf",width = 4.5,height = 3,family="serif")

plot_cell_trajectory(sce_CDS,show_branch_points = F,color_by = "Score",cell_size = 0.5)+
  scale_color_viridis()+theme_blank()
ggsave("./out/Monocle2_Score.pdf",width = 4.2,height = 3,family="serif")

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
ggsave("./out/Monocle2_Umap_cytotrace.pdf",width = 4.2,height = 3,family="serif")

##change root 更换起点
sce_CDS <-orderCells(sce_CDS,root_state = 1) #State
#sce_CDS <-orderCells(sce_CDS,reverse=T) #reverse
##re-plot

##Cor and Ridges----
if(T){
  #extract data - ridges
  plotdf=pData(sce_CDS)[25:27]
  colnames(plotdf)
  library(ggridges)
  ggplot(plotdf, aes(x=Pseudotime,y=celltype,color=celltype,fill=celltype))+
    geom_density_ridges(scale=1) +
    #geom_vline(xintercept = c(5,10),linetype=2)+
    scale_y_discrete("")+
    theme_minimal()+
    theme(
      panel.grid = element_blank()
    )+scale_fill_manual(values = my36colors)+scale_color_manual(values = my36colors)
  ggsave("./out/Ridges_Pseudotime.pdf",width = 5,height = 3,family="serif")
  
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
    ggsave(file=paste0('./out/Cor_',x,"_",y,"_",j,'.pdf'),width = 3,height = 3,family="serif")
    
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
  ggsave("./out/Boxplot_Pseudotime.pdf",width = 3,height = 3,family="serif")
  
  #umap Pseudotime
  for (palette in c("Blues","YlOrRd","Reds") ){
    print(palette)
    FeatureDimPlot(sce, features = "Pseudotime",raster = T,pt.size = 2,
                   reduction = "umap",palette = palette,#palcolor = c("skyblue","grey","pink"),#
                   label = F,show_stat = F,title = "",
                   theme_use = "theme_blank")
    ggsave(paste0("./out/FeatureDimPlot_monocle_sce_Pseudotime_umap_",palette,".pdf"),width = 5,height = 5,family="serif")
    
  }
  for (i in c("Pseudotime") ){
    print(i)
    VlnPlot(sce, features = i, group.by  = "Minor_Type",#log = T,
            stacked=T,pt.size= 0, #sort='increasing',
            cols = my36colors #颜色
            #direction = "horizontal", #水平作图 x.lab = '', y.lab = 'CytoTRACE'
    )+NoLegend()+labs(y=i,x="",title="")+theme(axis.text.x = element_text(angle = 90,size = 12))
    ggsave(paste0(folder_path,"/VlnPlot_Epi_vs_cancer_",i,".pdf"),width = 4,height = 3.5,family="serif")
  } 
  #save monocle2 results
  qsave(sce_CDS,"sce_Epi_CDS.qs")
}

##monocle+ggplot2热图+point curve----
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
  
  pdf("./out/monocle2-plot_pseudotime_heatmap.pdf",width = 4,height = 9,family = "serif")
  p <- plot_pseudotime_heatmap(sce_CDS[Time_genes,], 
                               num_cluster = 4, 
                               show_rownames = T, 
                               return_heatmap = T,
                               add_annotation_col= bin[c(2)],#2,3 can't add color in celltypes 
                               hmcols = colorRampPalette(c("navy","white","firebrick3"))(100))
  print(p)
  dev.off()
  
  pdf("./out/monocle2-plot_pseudotime_heatmap-1.pdf",width = 4,height = 9,family = "serif")
  p <- plot_pseudotime_heatmap(sce_CDS[Time_genes,], 
                               num_cluster = 4, 
                               show_rownames = T, 
                               return_heatmap = T,add_annotation_col= bin[2],
                               hmcols =colorRampPalette(rev(brewer.pal(9, "PRGn")))(100))
  print(p)
  dev.off()
  
  pdf("./out/monocle2-plot_pseudotime_heatmap-2.pdf",width = 4,height = 9,family = "serif")
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
  plotdf$celltype <-plotdf$celltype_minor
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
  p4
  ggsave("./out/monocle2_Epi_ggplot2-heatmap_celltype.pdf",width = 5,height = 10,family="serif")
  
  p5 
  ggsave("./out/monocle2_Epi_ridges_celltype.pdf",width = 5.5,height = 3,family="serif")
  
  ##ggplot重做基因拟时细胞散点图
  genes_exp <- list()
  for(i in 1:length(genes_1)){
    A <- log2(exprs(sce_CDS)[genes_1[i],]+1)
    A <- as.data.frame(A)
    genes_exp[[i]] <- A
  }
  
  gene_exp <- do.call(cbind, genes_exp)
  colnames(gene_exp) <- genes_1
  #将上述几个基因的拟时表达添加到monocle
  pData(sce_CDS) = cbind(pData(sce_CDS), gene_exp)
  
  #提取作图数据，只需要游基因表达和拟时即可
  data <- pData(sce_CDS)
  colnames(data)
  #选择需要的列即可，我这里的origin.ident就是分组
  data<-data %>% select("celltype","Pseudotime","CytoTRACE","Score",genes_1)
  
  features <- c("CytoTRACE","Score",genes_1)#c("Anxa1", "Ncf1","Ltf","Camp","Ngp","Chil3","S100a8","S100a9")
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
  ggsave("./out/monocle2_Epi_point.pdf",width = 10,height = 12,family="serif") #save
  
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
  ggsave("./out/monocle2_Epi_no_point.pdf",width = 6,height = 12,family="serif") #save
  
}

##差异基因拟时序热图优化----
library(ComplexHeatmap)
library(RColorBrewer)
library(tidyverse)
library(circlize)
genes1 <- genes#deg1 %>% arrange(desc(morans_I)) %>% pull(gene_short_name)
cds_list <- list(sce_CDS)

if(T){
  ###对矩阵行列重新排序
  library(Biobase);library(Seurat)
  # matrix_ori1 <- exprs(cds_list[[1]])[match(genes1, rownames(rowData(cds_list[[1]]))), order( )]
  # 
  # rownames(cds_list[[1]]$Pseudotime)
  # pseudotime(cds_list[[1]])
  # 
  # pData(cds_list[[1]])["Pseudotime"]
  
  
  ###定义色板
  cols <- my36colors[1:8]#c("#928ac7", "#e0fbd0", "#bdb9da", "#EEA2AD", "#fcb454", "#f7d96e", "#fae2b2", "#9ae76e")
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
                           Score=cell_anno$Score,
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
    split = factor(km1$cluster, levels = c(2,4,1,3)),
    cluster_rows = FALSE,
    cluster_row_slices = FALSE,
    cluster_columns = FALSE,
    top_annotation = ha1,
    use_raster = TRUE,
    border = TRUE
  )
  pdf("./out/monocle2_Epi_Heatmap_celltype-1.pdf",width = 5.5,height = 10,family="serif")
  print(p3)#save
  dev.off()
  
  #最后我们分别提取每个cluster的基因进行GO term富集，富集到的主要条目与文章一致。
  library(clusterProfiler);library(org.Hs.eg.db)
  model <- "org.Hs.eg.db"
  
  deg1 <- names(km1$cluster[km1$cluster=='4'])
  full <- bitr(deg1, fromType = "SYMBOL" ,toType = c("ENSEMBL","ENTREZID") , OrgDb = model)
  go1 <- enrichGO(gene = full$ENTREZID , OrgDb = model, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
  go1 <- head(as.data.frame(go1), 5)
  
  deg2 <- names(km1$cluster[km1$cluster=='3'])
  full <- bitr(deg2, fromType = "SYMBOL" ,toType = c("ENSEMBL","ENTREZID") , OrgDb = model)
  go2 <- enrichGO(gene = full$ENTREZID , OrgDb = model, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
  go2 <- head(as.data.frame(go2), 5)
  
  deg3 <- names(km1$cluster[km1$cluster=='2'])
  full <- bitr(deg3, fromType = "SYMBOL" ,toType = c("ENSEMBL","ENTREZID") , OrgDb = model)
  go3 <- enrichGO(gene = full$ENTREZID , OrgDb = model, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
  go3 <- head(as.data.frame(go3), 5)
  
  deg4 <- names(km1$cluster[km1$cluster=='1'])
  full <- bitr(deg4, fromType = "SYMBOL" ,toType = c("ENSEMBL","ENTREZID") , OrgDb = model)
  go4 <- enrichGO(gene = full$ENTREZID , OrgDb = model, ont = "BP", pAdjustMethod = "fdr", pvalueCutoff = 0.05, readable = TRUE)
  go4 <- head(as.data.frame(go4), 5)
  
  #https://mp.weixin.qq.com/s?__biz=MzkyOTU0MDgwMA==&mid=2247484443&idx=1&sn=9cc0240ee7a1593f207f1f315797174d&chksm=c3cfe099c924323c490945d1ff933fd883b3560d67bc0957c4ad61c466014e060e48f22af2fd&mpshare=1&scene=1&srcid=0809owc1nwRjujdchUIPK1HE&sharer_shareinfo=44f7c232a4691198802cf4c891a8dffa&sharer_shareinfo_first=44f7c232a4691198802cf4c891a8dffa#rd
  #save and patch in AI software
  
}

##Branch分析及热图可视化----
if(T){
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
                                         color_by = 'orig.ident',
                                         ncol = 1) +
        theme(plot.title = element_text(hjust = 0.5),legend.position = "right")
      p
      ggsave(p,file=paste0(outputDir,'/',"branch",i,"_genes_branched_pseudotime_orig",".pdf",sep=""),width = 9, height = 8,limitsize = FALSE,family="serif")
      
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
    ggsave(p,file=paste0(outputDir,'/',"pseudotime_dependent_gene_heatmap_top50.pdf"),width = 9, height = 5,limitsize = FALSE)
    #拆分拟时序分析轴上每个cluster的基因（拆分top50个基因）
    clusters <- cutree(p$tree_row, k = 2)
    clustering <- as.data.frame(clusters)%>%rownames_to_column(var = "gene")
    colnames(clustering) <- c("gene","cluster")
    write.table(clustering,file=paste0(outputDir,'/',"pseudotime_dependent_gene_heatmap_top50_gene_by_cluster.xls"),sep="\t",quote = F,col.names = T,row.names = F)
  }, error = function(e) {
    print("Fewer than 50 significant genes")
  })
  
} #slow 

#4.0 enrichment ----
##4.1 TF activity progeny in interested pseudotime trajectory ----
#run

if(T){
  i="Epi_all"
  ##progeny
  library(progeny);library(tidyr)
  set.seed(123)
  ## We compute the Progeny activity scores and add them to our Seurat object
  ## as a new assay called Progeny. 
  
  #sce <-subset(sceList,celltype_major %in% c("Cancer Epithelial","Normal Epithelial") )
  #table(sce@meta.data$celltype_minor)
  #Tumor or Epi
  sce$type <- ifelse(sce$celltype_minor %in% c("Luminal Progenitors","Mature Luminal","Myoepithelial"),"Normal","Tumor" )
  # Normal  Tumor 
  # 4338  24489 
  
  sce <- progeny(sce, scale=FALSE, organism="Human", top=500, perm=1, 
                  return_assay = TRUE)
  ## We can now directly apply Seurat functions in our Progeny scores. 
  ## For instance, we scale the pathway activity scores. 
  sce <- Seurat::ScaleData(sce, assay = "progeny") 
  ## We transform Progeny scores into a data frame to better handling the results
  progeny_scores_df <- 
    as.data.frame(t(GetAssayData(sce, slot = "scale.data", 
                                 assay = "progeny"))) %>%
    rownames_to_column("Cell") %>%
    gather(Pathway, Activity, -Cell) 
  CellsClusters <- data.frame(Cell = names(Idents(sce)), 
                              CellType = as.character(sce$celltype_minor), #sce$group
                              type=as.character(sce$type), 
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
  library(pheatmap)
  paletteLength = 100
  myColor = colorRampPalette(c("skyblue", "white","red"))(paletteLength)
  
  progenyBreaks = c(seq(min(summarized_progeny_scores_df), 0, 
                        length.out=ceiling(paletteLength/2) + 1),
                    seq(max(summarized_progeny_scores_df)/paletteLength, 
                        max(summarized_progeny_scores_df), 
                        length.out=floor(paletteLength/2)))
  progeny_hmap = pheatmap(t(summarized_progeny_scores_df[,-1]),fontsize=14, 
                          fontsize_row = 12, cluster_cols = F,gaps_col= 5,#treeheight_col = 0,#不显示聚类树
                          #cluster_rows  = F,
                          color=myColor, breaks = progenyBreaks, 
                          #main = "PROGENy", 
                          angle_col = 90, 
                          border_color = "grey70")#"NA"
  pdf(paste0("./out/Progeny_",i,"_group.pdf"),width = 4.6,height = 5,family="serif")# 3 5
  print(progeny_hmap)
  dev.off() #add sign. ?
  
  ##Epi vs tumor
  ggplot(data = progeny_scores_df,
         aes(x = reorder(CellType,-Activity), # CellType x = type,#
             y = Activity, fill = type))+
    facet_wrap(Pathway~. , scales='free',ncol = 7) +
    scale_fill_manual(values = my36colors ) + #my36colors_1 c("#008ECB", "#EA921D", "#D14039")
    scale_color_manual(values = my36colors) + # c("#D14039","#008ECB" )
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
    # stat_compare_means(#method = "t.test",
    #                    label.x = 1.5,#label.y = max(pathways_long$value),
    #                    # comparisons=list(
    #                    #   c("Cancer Her2 SC","Cancer LumA SC"),
    #                    #   c("Cancer Her2 SC","Luminal Progenitors"),
    #                    #   c("Cancer LumA SC","Luminal Progenitors")  ),
    #                    step.increase = 0.05,
    #                    label="p.signif",
    #                    vjust = 0.7,
    # )+#mytheme1+
    theme_bw(base_size = 12,base_family = "serif")+
    theme(#legend.position='none',
          axis.text.x = element_text(angle = 90,hjust = 1))
  #ggsave(paste0("./out/Progeny_",i,"_group_boxplot-epi-tumor-2group.pdf"),width = 10,height = 4.5,family="serif")#13 6.5
  ggsave(paste0("./out/Progeny_",i,"_group_boxplot-epi-tumor-1.pdf"),width = 13,height = 6.5,family="serif")#13 6.5
  
  
  ##macth color in boxplot
  my36colors_1 <-c('#57C3F3',"#F37C95",'#D6E7A3')#,"#20854E"
  #progeny boxplot
  library(ggplot2);library(ggsci);library(ggpubr)
  #colnames(progeny_scores_df)
  #progeny_scores_df <-subset(progeny_scores_df,progeny_scores_df$CellType %in% c("Luminal Progenitors","Cancer Her2 SC","Cancer LumA SC" ) )#"Cancer LumB SC",
  #progeny_scores_df$CellType <- factor(progeny_scores_df$CellType, levels = c("Luminal Progenitors", "Cancer Her2 SC",  "Cancer LumA SC")) #"Cancer LumB SC",  # 指定顺序
  
  ggplot(data = progeny_scores_df,
         aes(x = CellType, y = Activity, fill = CellType))+
    facet_wrap(Pathway~. , scales='free',ncol = 7) +
    scale_fill_manual(values = my36colors ) + #my36colors_1 c("#008ECB", "#EA921D", "#D14039")
    scale_color_manual(values = my36colors) + # c("#D14039","#008ECB" )
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
    # stat_compare_means(#method = "t.test",
    #                    label.x = 1.5,#label.y = max(pathways_long$value),
    #                    comparisons=list(
    #                      c("Cancer Her2 SC","Cancer LumA SC"),
    #                      c("Cancer Her2 SC","Luminal Progenitors"),
    #                      c("Cancer LumA SC","Luminal Progenitors")  ),
    #                    step.increase = 0.05,
    #                    label="p.signif",
    #                    vjust = 0.7,
    # )+#mytheme1+
    theme_bw(base_size = 12,base_family = "serif")+theme(legend.position='none',axis.text.x = element_text(angle = 90,hjust = 1))
  ggsave(paste0("./out/Progeny_",i,"_group_boxplot-2.pdf"),width = 12,height = 6.5,family="serif")#15 6
  
}#end 
##show plot




##4.2 Enrichment in subtypes


##tumor high/low enrichment and cellchat
#5.0 cellchat


#
#T cell dysfunctional and cytotoxic effector gene signature scores
#The cytotoxic gene list consists of 12 genes that translate to effector cytotoxic proteins
#(GZMA, GZMB, GZMH, GZMK, GZMM, GNLY, PRF1 and FASLG) and well-described cytotoxic T cell activation markers (IFNG, TNF, IL2R and IL2).

#5.0 exp and Cor of MHC molecules----

