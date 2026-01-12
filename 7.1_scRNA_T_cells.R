
#aim:in T cells

if(T){
  rm(list = ls())
  setwd("~/R/my_projects/JAK_STAT_BC/7.0_scRNA_Immune_cells/7.1_scRNA_T_cells")
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
  
  genelist<- list("Jak-STAT_pathway"=genes,"Core_Jak-STAT_pathway"=genes_1
                  #"Proliferation"= c("BIRC5", "CCNB1", "CDC20", "NUF2", "CEP55", "NDC80", "MKI67", "PTTG1", "RRM2", "TYMS","UBE2C")
  ) #https://www.nature.com/articles/s41588-021-00911-1
  input_gene= c("Jak-STAT_pathway","Core_Jak-STAT_pathway") #genes#c(genes,"H3K36me","5mC" )
  project = "Jak_STAT_BC"
}

#sceList <-qread("~/R/my_projects/JAK_STAT_BC/7.0_scRNA_Immune_cells/sce_Tcell.qs")
sceList <-qread("sce_Tcell.qs")
#qsave(sceList,"sce_Tcell.qs")

#0. exp----
##亚群表达----
for (j in c("Score",interest_gene)){ #,interest_gene ,colnames(sceList@meta.data)[16:21] 
  print(j)
  for (i in c("celltype_subset") ){ #"orig.ident","subtype",
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
    
    for (k in c("RdPu","PRGn","PiYG","RdYlGn","Blues","Greens") ){
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
        theme_bw()+
        theme(#panel.grid = element_blank(), 
          axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
      hei1 <-ifelse(i=="orig.ident",height,height/2)
      ggsave(paste0(folder_path,"/Dotplot_Score_",j,"_",i,"_",k,".pdf"),width = 4,height = hei1+0.8,family="serif")
      
      }
  
    ##
  }
  
}
##样本：细胞比例~均值 相关性----
#计算样本T细胞亚群百分百
Cellratio <- prop.table(table(sceList$celltype_subset ,#celltype_minor1206细亚群  #celltype_major 20221204
                              sceList$orig.ident), margin = 2)#计算各组样本不同细胞群比例
Cellratio <-data.frame(Cellratio)
Cellratio$Freq <- Cellratio$Freq*100
library(reshape2)
cellper <- dcast(Cellratio,Var2~Var1, value.var = "Freq")#长数据转为宽数据
colnames(cellper)[1]<-"id"

#计算样本Score均值
Srps <-data.frame(srps=sceList$Score,sample=sceList$orig.ident )#99301-1 如何计算每个样本比值？
table(Srps$sample)

mean_srps<-data.frame()
for (i in unique(Srps$sample)){
  mean_1 <- data.frame(Srps=mean(Srps[Srps$sample==i,]$srps),id=i)
  mean_srps<-rbind(mean_srps,mean_1)
}
colnames(mean_srps)[1] <- "Score"
#合并细胞比例到mean_srps
all_merge <- merge(mean_srps,cellper,by="id")  #26-14

##plot Cor 相关性-棒棒糖图lollipop----
y=colnames(all_merge)[3:14]#y #[2:30]
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
  scale_y_continuous(limits = c(-0.5,0.5))+ #c(-0.5,1.0)
  geom_hline(yintercept=c(0),size=0.1)+
  geom_text(size=5, aes( label = Sign),vjust=1) +#hjust = 1,
  theme(axis.text.y = element_text(size = 12 ),panel.grid = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5) )  #,legend.position = c(0.7,0.7)
p
write.csv(outTab,paste0(folder_path,"/outTab_cor_ratio_minor.csv") )
ggsave(paste0(folder_path,"/cor_ratio_minor_lollipop.pdf"),family="serif",width = 4,height = 4)#4-3

##批量相关性  dot 散点图 ----
exp<-column_to_rownames(all_merge,"id")
x="Score"
ys=colnames(all_merge)[3:14]#y #[2:30]

methods <-c("pearson","spearman")
for (j in methods){
  print(paste0(j," :"))
    print(x)
    for (y in ys){
      print(y)
      #run cor
      df <-exp[c(x,y)]
      df <- na.omit(df) #删除没有配对的NA
      colnames(df)<-c("X","Y")
      ggplot(df, aes(x=X, y=Y)) + 
        #xlim(-20,15) + ylim(-15,10) +
        labs(x = x, y = "% of T cells",title=y) +
        geom_point(color="skyblue",size = 1.5,alpha=0.5) + #aes(color=cancer)
        geom_smooth(method ='lm', size=0.5,color="navy") + #skyblue
        #ggrepel::geom_text_repel(aes(label=cancer))+ #geom_label_repel
        stat_cor(method = j,size = 6,color="red") +
        scale_colour_manual(values = colors ) +
        theme_bw() + 
        theme(axis.text.x = element_text(size = 14), 
              axis.text.y = element_text(size = 14), 
              axis.title.x = element_text(size = 16), 
              axis.title.y = element_text(size = 16),panel.grid = element_blank(),
              legend.title = element_blank(),
              legend.text = element_text(size = 12))
      ggsave(file=paste0(folder_path,'/Cor_',x,"_",y,"_",j,'.pdf'),width = 3,height = 3.1,family="serif")
      
    }
  
}

#1.0 monocle+cytotrace----
for (i in names(sceList@reductions) ){
  print(i)
  sceList@reductions[i] <-NULL
}
celltype="celltype_minor"

##re-umap
library(BiocParallel)
register(MulticoreParam(workers = 48, progressbar = TRUE))
sceList <-Standard_SCP(sceList,cluster_resolution = 0.1) #,cluster_resolution = 0.1
##plot uamp Epi
#"StandardUMAP2D"
CellDimPlot(
  srt = sceList, group.by =celltype,#label = T,#, "Standardpca_SNN_res.0.8"
  reduction = "StandardUMAP2D", 
  palette = "Paired",raster = T,pt.size = 2,
  theme_use = "theme_blank"
)
ggsave(paste0(folder_path,"/UMAP_reUMAP.pdf"),width = 4.2,height = 3,family="serif")

##cytotrace----

if(T){
  library(CytoTRACE)
  ####提取表型文件
  phe <- sceList@meta.data[,celltype]
  phe = as.character(phe)
  names(phe) <- rownames(sceList@meta.data)
  ####提取表达矩阵
  mat_exp <- as.matrix(sceList@assays$RNA@counts)
  results <- CytoTRACE(mat = mat_exp,ncores = 24) #slow
  
  plotCytoGenes(results, numOfGenes = 10,outputDir = paste0(folder_path,"/") )#
  plotCytoTRACE(results, phenotype = phe,gene = interest_gene,outputDir = paste0(folder_path,"/") )
  plotCytoTRACE(results, phenotype = phe,gene = interest_gene,
                emb = sceList@reductions[["StandardUMAP2D"]]@cell.embeddings ,
                outputDir = paste0(folder_path,"/"))
  
  ####plot
  #boxplot
  cyto_res <-as.data.frame(results$CytoTRACE)#read.table("./cytotrace_outCytoTRACE_plot_table.txt")
  colnames(cyto_res)<-"CytoTRACE"
  #add to meta
  sceList <- AddMetaData(sceList,cyto_res[1])
  
  #plot umap #re-cluster?
  for (i in c("Score","CytoTRACE",interest_gene)){
    print(i)
    FeatureDimPlot(sceList, features = i,raster = T,pt.size = 2,
                   reduction = "umap",palette = "YlOrRd",#"Reds",#
                   theme_use = "theme_blank",label = F,show_stat = F,title = ""
    )
    ggsave(paste0(folder_path,"/FeatureDimPlot_",i,".pdf"),width = 5,height = 5,family="serif")
  }
  
  ##Cor计算相关性
  x="Score";y="CytoTRACE"
  methods <-c("pearson","spearman")
  for (j in methods){
    print(j)
    df <-sceList@meta.data[c(x,y,celltype)]
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
  
  #dot plot/pseudotime 拟时序
  DynamicPlot(
    srt = sceList, lineages = "CytoTRACE", group.by = celltype,
    features = c(genes_1,"Score"),#c("Plk1", "Hes1", "Neurod2", "Ghrl", "Gcg", "Ins2"),
    compare_lineages = TRUE, compare_features = FALSE,add_rug =F
  )#too large add_point = TRUE,
  ggsave(paste0(folder_path,"/Dotplot_pseudotime_CytoTRACE.pdf"),width = 12,height = 8,family="serif")
  DynamicPlot(
    srt = sceList, lineages = "CytoTRACE", group.by = celltype,
    features = c(genes_1,"Score"),#c("Plk1", "Hes1", "Neurod2", "Ghrl", "Gcg", "Ins2"),
    compare_lineages = T, compare_features = FALSE,add_point = F,add_rug =F
  ) 
  ggsave(paste0(folder_path,"/Dotplot_pseudotime_CytoTRACE_nopoint.pdf"),width = 12,height = 6,family="serif")
  
  #plot exp violin
  FeatureStatPlot(
    srt = sceList, group.by = celltype, bg.by =celltype,palcolor = my36colors,bg_palcolor =my36colors,bg_apha =0.1,
    stat.by = c(genes_1,"Score","CytoTRACE"), add_box = TRUE,legend.position = "none" #,
    # comparisons = list(
    #   c("Luminal_progenitor", "Mature_luminal"),
    #   c("Mature_luminal", "Basal"),
    #   c("Luminal_progenitor", "Basal")
    # )
  )
  ggsave(paste0(folder_path,"/Boxplot_Score_genes_CytoTRACE.pdf"),width = 12,height = 9,family="serif")
  #violin
  library(MySeuratWrappers)
  for (i in c("CytoTRACE",interest_gene) ){
    print(i)
    VlnPlot(sceList, features = i, group.by  = celltype,
            stacked=T,pt.size= 0, #sort='increasing',
            cols = my36colors #颜色
            #direction = "horizontal", #水平作图 x.lab = '', y.lab = 'CytoTRACE'
    )+NoLegend()+labs(y=i,x="",title="")+theme(axis.text.x = element_text(angle = 90,size = 12))
    ggsave(paste0(folder_path,"/VlnPlot_Tcells_",i,".pdf"),width = 4,height = 3.5,family="serif")
    
  }

}

##monocle 2----
if(T){
  library(Seurat)
  library(monocle)
  library(viridis)
}
#names(sceList@reductions)[4] <-"UMAP"#"umap"
DimPlot(sceList, reduction = "UMAP",label = T, pt.size = 1)
#reduction = "UMAP",Error in reduction %||% Seurat:::DefaultDimReduc(object = object) :    object 'DefaultDimReduc' not found

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
Idents(sceList)<-sceList@meta.data[,celltype]
table(Idents(sceList))
# T cells CD8+    T cells CD4+        NK cells Cycling T-cells       NKT cells 
# 11436           19087            1844            1522            1112 

sce$celltype
sce <- subset(sceList, downsample = 1000) # 随机抽取1k个细胞
sce$celltype <-sce@meta.data[,celltype]

table(Idents(sce)) # 查看抽样后各群体的细胞数量

sce_CDS <- ks_run_Monocle2(object = sce,#
                           layer = 'counts',
                           assay = "RNA",
                           VARgenesM="dispersionTable",
                           cellAnno = celltype)
##Plot monocle2
#可视化
if(T){
  plot_cell_trajectory(sce_CDS,show_branch_points = F,color_by = "Pseudotime",cell_size = 0.5)+
    scale_color_viridis()+theme_blank()
  ggsave(paste0(folder_path,"/Monocle2_pseudotime.pdf"),width = 4.2,height = 3,family="serif")
  
  plot_cell_trajectory(sce_CDS,show_branch_points = F,color_by = "State",cell_size = 0.5)+
    scale_color_brewer(palette="Set1")+theme_blank()#scale_color_manual(values=my36colors)#
  ggsave(paste0(folder_path,"/Monocle2_State.pdf"),width = 4.2,height = 3,family="serif")
  
  plot_cell_trajectory(sce_CDS,show_branch_points = F,color_by = celltype,cell_size = 0.5)+
    scale_color_manual(values=my36colors)+theme_blank()
  ggsave(paste0(folder_path,"/Monocle2_Minor_Type.pdf"),width = 4.5,height = 3,family="serif")
  
  plot_cell_trajectory(sce_CDS,show_branch_points = F,color_by = "Score",cell_size = 0.5)+
    scale_color_viridis()+theme_blank()
  ggsave(paste0(folder_path,"/Monocle2_Score.pdf"),width = 4.2,height = 3,family="serif")
  
  for (i in c("RdPu","PRGn","PiYG","RdYlGn","Blues","Greens") ){
    print(i)
    plot_cell_trajectory(sce_CDS,show_branch_points = F,color_by = "Score",cell_size = 0.5)+
      scale_color_distiller(palette = i,direction = -1) +theme_blank()
    ggsave(paste0(folder_path,"/Monocle2_Score_",i,".pdf"),width = 4.2,height = 3,family="serif")
  }
  
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
  
}

##change root 更换起点
sce_CDS <-orderCells(sce_CDS,root_state = 2) #choose State
#sce_CDS <-orderCells(sce_CDS,reverse=T) #reverse
##re-plot
##Cor and Ridges----
if(T){
  #extract data - ridges
  plotdf=pData(sce_CDS)#[25:27]# Pseudotime Pseudotime
  colnames(plotdf)
  plotdf=plotdf[,c(celltype,"Pseudotime","Score")]
  colnames(plotdf)[1]<-"celltype"
  
  library(ggridges)
  ggplot(plotdf, aes(x=Pseudotime,y=celltype,color=celltype,fill=celltype))+
    geom_density_ridges(scale=1) +
    #geom_vline(xintercept = c(5,10),linetype=2)+
    scale_y_discrete("")+
    theme_minimal()+
    theme(
      panel.grid = element_blank()
    )+scale_fill_manual(values = my36colors)+scale_color_manual(values = my36colors)
  ggsave(paste0(folder_path,"/Ridges_Pseudotime.pdf"),width = 5,height = 3,family="serif")
  
  sce <-AddMetaData(sce,meta=plotdf)
  ##Cor计算相关性
  x="Score";y="Pseudotime"#"CytoTRACE"
  methods <-c("pearson","spearman")
  for (j in methods){
    print(j)
    df <-sce@meta.data[c(x,y,celltype)]
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
    srt = sce, group.by = celltype, bg.by = celltype,palcolor = my36colors,bg_palcolor =my36colors,bg_apha =0.1,
    stat.by = c("Pseudotime"), add_box = TRUE,legend.position = "none"#,
    # comparisons = list(
    #   c("Luminal_progenitor", "Mature_luminal"),
    #   c("Mature_luminal", "Basal"),
    #   c("Luminal_progenitor", "Basal")
    # )
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
  for (i in c("Pseudotime") ){
    print(i)
    VlnPlot(sce, features = i, group.by  = celltype,#log = T,
            stacked=T,pt.size= 0, #sort='increasing',
            cols = my36colors #颜色
            #direction = "horizontal", #水平作图 x.lab = '', y.lab = 'CytoTRACE'
    )+NoLegend()+labs(y=i,x="",title="")+theme(axis.text.x = element_text(angle = 90,size = 12))
    ggsave(paste0(folder_path,"/VlnPlot_Epi_vs_cancer_",i,".pdf"),width = 4,height = 3.5,family="serif")
  } 
  #save monocle2 results
  qsave(sce_CDS,"sce_Tcell_CDS.qs")
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
    df <- df[,c("Pseudotime", "State",celltype)]
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
  
  pdf(paste0(folder_path,"/monocle2-plot_pseudotime_heatmap-2.pdf"),width = 4,height = 9,family = "serif")
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
  p4 #趋势相反？
  ggsave(paste0(folder_path,"/monocle2_Epi_ggplot2-heatmap_celltype.pdf"),width = 5,height = 10,family="serif")
  
  p5 
  ggsave(paste0(folder_path,"/monocle2_Epi_ridges_celltype.pdf"),width = 5.5,height = 3,family="serif")
  
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
  #data$celltype <-data[,celltype]
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

genes1 <- genes#deg1 %>% arrange(desc(morans_I)) %>% pull(gene_short_name)
cds_list <- list(sce_CDS)

if(T){
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


#2.0 CellChat in all T cells----
#see mean group in code 7.0
# High   Low  mean
# 16778 18223 

sceList$Type <- sceList$group#ifelse(sceList$orig.ident %in% as.character(mean_srps[mean_srps$Type=="High",]$id),"High","Low" )
table(sceList$Type)
sceList$Cell_type <-sceList$celltype_minor #sceList$celltype_subset
table(sceList$Cell_type)

##run CellChat in T cell----
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
}

##plot CellChat----
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
  number_cancer_cluster=1 #选择兴趣细胞
}

#手动-调整图大小
if(T){
  ####1.0 load data----
  #合并cellchat对象:
  object.list <- list(low = low,
                      high = high) #对照组(NL)在前，比较组(LS)在后，注意顺序
  cellchat <- mergeCellChat(object.list, add.names = names(object.list))
  #cellchat
  
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
  pdf(paste0(folder_path,"/cellchat_net_all_noweight_vs_weight.pdf"),width = 10,height = 5,family = "serif")
  par(mfrow = c(1,2), xpd = TRUE)
  netVisual_diffInteraction(cellchat, weight.scale = T)
  netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
  dev.off()
  ###high or low group
  # i=1
  # netVisual_diffInteraction(object.list$low)
  
  #2.2 细胞通讯差异热图
  pdf(paste0(folder_path,"/cellchat_net_heatmap_all.pdf"),width = 6,height = 3.5,family = "serif")
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
  ggsave(paste0(folder_path,"/cellchat_rankNet_high-low.pdf"),width = 8,height = 2.5,family = "serif")
  
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
  pdf(paste0(folder_path,"/cellchat_net_outgoing_", names(object.list[i]),"_", names(object.list[i+1]),".pdf"),width = 8,height = 6,family = "serif")
  draw(ht1 + ht2)#, ht_gap = unit(0.5, "cm")
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
  pdf(paste0(folder_path,"/cellchat_bubble_all_high-low.pdf"),width = 8,height = 4,family = "serif")#21 10
  netVisual_bubble(cellchat,angle.x = 45,comparison = c(1, 2) ) #all
  dev.off()
  
  levels(cellchat@idents$joint) #查看细胞亚群
  pdf(paste0(folder_path,"/cellchat_bubble_choose_high-low.pdf"),width = 6,height = 4,family = "serif")
  netVisual_bubble(cellchat,
                   sources.use = number_cancer_cluster, #Malignant cells
                   #targets.use = c(1:8),
                   font.size = 12,font.size.title = 14,
                   comparison = c(1, 2),
                   angle.x = 45)
  dev.off()
  
}


#for (pathways.show in unique(cellchat@DB$interaction$pathway_name) )
if(T){
  
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
    pdf(paste0(folder_path,"/cellchat_heatmap_pathway_",pathways.show,"_high-low.pdf"),width = 6,height = 3.5,family = "serif")
    
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
  
  ##show unique pathway in high
  pathways.shows <- object.list[["high"]]@netP[["pathways"]] 
  pathways.shows
  for (pathways.show in pathways.shows ){
    print(pathways.show)
    pdf(paste0(folder_path,"/cellchat_heatmap_pathway_",pathways.show,"_only_high.pdf"),width = 4.5,height = 4,family = "serif")
    p<-netVisual_heatmap(object.list[["high"]],font.size = 12,font.size.title = 14,
                         signaling = pathways.show,
                         color.heatmap = "Greens",#"Reds",
                         title.name = paste(pathways.show, "signaling ",names(object.list)["high"]))
    print(p)
    dev.off()
  }
  
}
#保存sce T cells cytotrace+monocle2
qsave(sceList,"sce_Tcell.qs")

#3.0 enrichment----
#T细胞 mean score group high/low in 7.0_scRNA_Immune_cells.R

#4.0 co-stimulators, co-inhibitors, and T-function markers----
##计算平均表达
library(ggplot2)
library(Seurat)
library(dplyr)
library(patchwork)

DimPlot(sceList, label = T)
#dat=AverageExpression(sceList, group.by = "celltype_minor")$RNA
dat=AverageExpression(sceList, group.by = "celltype_subset")$RNA
#dat1=AverageExpression(sceList, group.by = "group")$RNA
#dat<-cbind(dat,dat1)

##pheatmap
library(pheatmap)
# pheatmap(dat[VariableFeatures(sceList),], border_color = NA,
#          scale = 'row',cluster_cols = F,
#          clustering_method = "ward.D2",
#          show_rownames = F#,main="PBMC 3k"
#          )

immune_molecules<-read.csv("T_immune_molecules.csv")
#https://advanced.onlinelibrary.wiley.com/doi/10.1002/advs.202500722

rownames(immune_molecules)<-immune_molecules$id
# gene.set=unique(c("PTPRC", "CD3D","CD3E", "CD3G","CD4", "CD8A", 
#                   'CCR7','IL7R', 'S100A4', 'KLRG1',
#                   'GZMK',
#                   "GNLY", "NKG7", 
#                   "CD14", "LYZ", 'CST3',
#                   "CD79A", "CD79B", "MS4A1", 
#                   "FCER1A", "FCGR3A", 
#                   "PPBP", "PF4"))
gene.set=immune_molecules$id
pheatmap(dat[gene.set,], border_color = "white",
         scale = 'row',cluster_cols = F,cluster_rows  = F,
         clustering_method = "ward.D2",#gaps_col = 5,#
         gaps_row = c(28,28+28,28+28+4),
         show_rownames = T,
         annotation_row = immune_molecules["type"]
         )
table(immune_molecules$type)
# Co-inhibitor Co-stimulator    T-function 
# 28            28             4 
save_pheatmap_pdf <- function(x, filename, width=width, height=height) {
  library(grid)
  x$gtable$grobs[[1]]$gp <- gpar(lwd = 0.05)#修改聚类树线条宽度
  #x$gtable$grobs[[2]]$vp <- vpar(lwd = 0.05)
  
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height,family = "serif")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}#保存热图的函数

p<-pheatmap(dat[gene.set,], scale = 'row',
            cluster_col = F,cluster_rows  = F,
            #legend_breaks = c( -8,seq(from=-0.2, to=max(a), by=2)),
            legend_labels = c("NA",seq(from=0, to=8, by=2)),
            #gaps_col = 5,
            gaps_row = c(28,28+28,28+28+4),
            border = T,border_color =NA, #gaps_col=3, "black"
            display_numbers = F,na_col = "#DDDDDD", #
            #color = c(colorRampPalette(colors = c("skyblue","firebrick3"))(100)), #colorRampPalette(c("navy", "white", "firebrick3"))(50)
            #cellwidth = 15,cellheight = 15,
            annotation_row = immune_molecules["type"],
            legend = TRUE,fontsize = 11,#angle_col = 45#,filename = "1.pdf" 
            color=#colorRampPalette( c("navy", "white", "firebrick3"))(10)#
            alpha(rev(RColorBrewer::brewer.pal(11,"Spectral") ), 0.8)#,   #"BrBG" 褐浅蓝  "PiYG"紫绿
            #alpha(colorRampPalette(colors = c("orange", "white", "purple"))(10),0.6) ,#colorRampPalette( c("navy", "white", "firebrick3"))(10) ,#
            
)
#save_pheatmap_pdf(p, paste0(folder_path,"/pheatmap_T_molecules-3.pdf"),width=4, height=8)
save_pheatmap_pdf(p, paste0(folder_path,"/pheatmap_T_subset_molecules.pdf"),width=5, height=9)

#group
# p<-pheatmap(scale(dat1[gene.set,]), #scale = 'column',#scale = 'row',
#             cluster_col = F,cluster_rows  = F,
#             #legend_breaks = c( -8,seq(from=-0.2, to=max(a), by=2)),
#             legend_labels = c("NA",seq(from=0, to=8, by=2)),
#             gaps_col = 1,
#             gaps_row = c(28,28+28,28+28+4),
#             border = T,border_color =NA, #gaps_col=3, "black"
#             display_numbers = F,na_col = "#DDDDDD", #
#             #color = c(colorRampPalette(colors = c("skyblue","firebrick3"))(100)), #colorRampPalette(c("navy", "white", "firebrick3"))(50)
#             #cellwidth = 15,cellheight = 15,
#             annotation_row = immune_molecules["type"],
#             legend = TRUE,fontsize = 11#,angle_col = 45#,filename = "1.pdf" 
# ) #no sence

##complexpheatmap


##high low score boxplot
list <-split(immune_molecules[1],immune_molecules[2])

library(BiocParallel)
register(MulticoreParam(workers = 24, progressbar = TRUE))#32
#score
sceList <- AddModuleScore(object = sceList,#sce,
                          features = list,#gene,
                          name = names(list )#"score" #ctrl = 100,
) #
colnames(sceList@meta.data)[27:29] <- gsub("\\d","",colnames(sceList@meta.data)[27:29] )

##dot violin in high/low----
for (j in unique(immune_molecules$type) ){ #,interest_gene ,colnames(sceList@meta.data)[16:21] 
  print(j)
  for (i in c("group") ){ #"orig.ident","subtype",
    print(i)
    ##
    library(MySeuratWrappers)
    Idents(sceList) <-sceList$orig.ident
    VlnPlot(sceList, features = j, group.by  = i,
            stacked=T,pt.size= 0, sort='increasing',
            cols = my36colors, #颜色
            #direction = "horizontal", #水平作图
            x.lab = '', y.lab = '')+NoLegend()+labs(y=j,title="")+theme(axis.text.x = element_text(angle = 90,size = 12,vjust = 0.5))
    width = ifelse(length(unique(sceList@meta.data[,i]))/5 >4,length(unique(sceList@meta.data[,i]))/5,2)
    height = ifelse(nchar(max(as.character(unique(sceList@meta.data[,i]))))>10,4,3)
    ggsave(paste0(folder_path,"/VlnPlot_Score_",j,"_",i,".pdf"),width = width,height = height,family="serif")
    
    for (k in c("RdPu","PRGn","PiYG","RdYlGn","Blues","Greens") ){
      print(k)
      #dotplot
      DotPlot(sceList , features =j,
              group.by = i,
              #cols = c("skyblue", "pink")
      )+ #Minor_Type orig.ident
        #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
        #scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
        #scale_color_viridis_c(option = "A")+
        scale_color_distiller(palette = k,direction = -1) +
        labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
        theme_bw()+
        theme(#panel.grid = element_blank(), 
          axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
      hei1 <-ifelse(i=="orig.ident",height,height/2)
      ggsave(paste0(folder_path,"/Dotplot_Score_",j,"_",i,"_",k,".pdf"),width = 2.5,height = hei1+0.8,family="serif")
      
    }
    
    ##
  }
}

#three in one dot
for (k in c("RdPu","PRGn","PiYG","RdYlGn","Blues","Greens") ){
  print(k)
  #dotplot
  DotPlot(sceList , features =unique(immune_molecules$type),
          group.by = i,
          #cols = c("skyblue", "pink")
  )+ #Minor_Type orig.ident
    #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
    #scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
    #scale_color_viridis_c(option = "A")+
    scale_color_distiller(palette = k,direction = -1) +
    labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
    theme_bw()+
    theme(#panel.grid = element_blank(), 
      axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
  hei1 <-ifelse(i=="orig.ident",height,height/2)
  ggsave(paste0(folder_path,"/Dotplot_Tscore_",i,"_",k,".pdf"),width = 3.5,height = hei1+0.8,family="serif")
  
}

##dot violin in minor type: celltype_minor----
for (j in unique(immune_molecules$type) ){ #,interest_gene ,colnames(sceList@meta.data)[16:21] 
  print(j)
  for (i in c("celltype_minor") ){ #"orig.ident","subtype",
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
    height = ifelse(nchar(max(as.character(unique(sceList@meta.data[,i]))))>10,4,3)
    ggsave(paste0(folder_path,"/VlnPlot_Score_",j,"_",i,".pdf"),width = width,height = height,family="serif")
    
    for (k in c("RdPu","PRGn","PiYG","RdYlGn","Blues","Greens") ){
      print(k)
      #dotplot
      DotPlot(sceList , features =j,
              group.by = i,
              #cols = c("skyblue", "pink")
      )+ #Minor_Type orig.ident
        #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
        #scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
        #scale_color_viridis_c(option = "A")+
        scale_color_distiller(palette = k,direction = -1) +
        labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
        theme_bw()+
        theme(#panel.grid = element_blank(), 
          axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
      hei1 <-ifelse(i=="orig.ident",height,height/2)
      ggsave(paste0(folder_path,"/Dotplot_Score_",j,"_",i,"_",k,".pdf"),width = 3.5,height = hei1+0.8,family="serif")
      
    }
    
    ##
  }
}

#three in one dot
for (k in c("RdPu","PRGn","PiYG","RdYlGn","Blues","Greens") ){
  print(k)
  #dotplot
  DotPlot(sceList , features =unique(immune_molecules$type),
          group.by = i,
          #cols = c("skyblue", "pink")
  )+ #Minor_Type orig.ident
    #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
    #scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
    #scale_color_viridis_c(option = "A")+
    scale_color_distiller(palette = k,direction = -1) +
    labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
    theme_bw()+
    theme(#panel.grid = element_blank(), 
      axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
  hei1 <-ifelse(i=="orig.ident",height,height/2)
  ggsave(paste0(folder_path,"/Dotplot_Tscore_",i,"_",k,".pdf"),width = 4.2,height = hei1+0.8,family="serif")
  
}

##dot violin in subset type----
for (j in unique(immune_molecules$type) ){ #,interest_gene ,colnames(sceList@meta.data)[16:21] 
  print(j)
  for (i in c("celltype_subset") ){ #"orig.ident","subtype",
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
    ggsave(paste0(folder_path,"/VlnPlot_Score_",j,"_",i,".pdf"),width = width,height = height,family="serif")
    
    for (k in c("RdPu","PRGn","PiYG","RdYlGn","Blues","Greens") ){
      print(k)
      #dotplot
      DotPlot(sceList , features =j,
              group.by = i,
              #cols = c("skyblue", "pink")
      )+ #Minor_Type orig.ident
        #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
        #scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
        #scale_color_viridis_c(option = "A")+
        scale_color_distiller(palette = k,direction = -1) +
        labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
        theme_bw()+
        theme(#panel.grid = element_blank(), 
          axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
      hei1 <-ifelse(i=="orig.ident",height,height/2)
      ggsave(paste0(folder_path,"/Dotplot_Score_",j,"_",i,"_",k,".pdf"),width = 4,height = hei1+0.8,family="serif")
      
    }
    
    ##
  }
}

#three in one dot
for (k in c("RdPu","PRGn","PiYG","RdYlGn","Blues","Greens") ){
  print(k)
  #dotplot
  DotPlot(sceList , features =unique(immune_molecules$type),
          group.by = i,
          #cols = c("skyblue", "pink")
  )+ #Minor_Type orig.ident
    #coord_flip()+ #theme_bw()+ #去除背景，旋转图片
    #scale_colour_gradientn( values= seq(0,1,0.2),colours = c( '#330066', '#336699', '#66CC66', '#FFCC33'))+ #颜色渐变设置 
    #scale_color_viridis_c(option = "A")+
    scale_color_distiller(palette = k,direction = -1) +
    labs( x=NULL, y=NULL)+#guides(size=guide_legend(order= 3))+
    theme_bw()+
    theme(#panel.grid = element_blank(), 
      axis.text.x=element_text(angle= 90,hjust = 1,vjust= 0.5)) 
  hei1 <-ifelse(i=="orig.ident",height,height/2)
  ggsave(paste0(folder_path,"/Dotplot_Tscore_",i,"_",k,".pdf"),width = 4.2,height = hei1+0.8,family="serif")
  
}

##Cor Score 相关性----
for (y in unique(immune_molecules$type) ){ #unique(gene_set$Pathway) colnames(sceList@meta.data)[16:21]
  print(y)
  x="Score"
  ##
  methods <-c("pearson","spearman")
  for (j in methods){
    print(j)
    df <-sceList@meta.data[c(x,y,"celltype_minor")] #celltype_minor
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

#Cor in NKT
table(sceList@meta.data[["celltype_minor"]])
for (y in unique(immune_molecules$type) ){ #unique(gene_set$Pathway) colnames(sceList@meta.data)[16:21]
  print(y)
  x="Score"
  ##
  methods <-c("pearson","spearman")
  for (j in methods){
    print(j)
    sceList1<-subset(sceList,celltype_minor== "NKT cells")
    df <-sceList1@meta.data[c(x,y,"celltype_minor")] #celltype_minor
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
    ggsave(file=paste0(folder_path,"/Cor_point_",x,"_",y,"_",j,"_NKT.pdf"),width = 3,height = 3,family="serif")
  }
  ##
}

#Cor in T_cells_c10_NKT_cells_FCGR3A
table(sceList@meta.data[["celltype_subset"]])

for (y in unique(immune_molecules$type) ){ #unique(gene_set$Pathway) colnames(sceList@meta.data)[16:21]
  print(y)
  x="Score"
  ##
  methods <-c("pearson","spearman")
  for (j in methods){
    print(j)
    sceList1<-subset(sceList,celltype_subset== "T_cells_c10_NKT_cells_FCGR3A")
    df <-sceList1@meta.data[c(x,y,"celltype_minor")] #celltype_minor
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
    ggsave(file=paste0(folder_path,"/Cor_point_",x,"_",y,"_",j,"_c10_NKT_cells_FCGR3A.pdf"),width = 3,height = 3,family="serif")
  }
  ##
}

