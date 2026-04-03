
if(T){
  rm(list = ls())
  setwd("~/R/my_projects/JAK_STAT_BC/13.0_TNBC")
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



#1.0 load and RCTD----
##1.1 load ref sce----
if(F){
  # 加载必要包
  library(Seurat)
  library(RCTD)
  library(Matrix)
  sc <-qs::qread("~/R/my_projects/JAK_STAT_BC/10.0_anti-PD1_Nat.Med/anti-PD1-Nat.Med.qs")
  
  # 1. 提取原始计数矩阵并确保数值为整数（dgCMatrix保留double类型）
  counts <- sc@assays$RNA@counts
  # 处理浮点误差（如1.0000001 → 1.0）
  counts@x <- round(counts@x)
  # 验证counts数值为整数（可选）
  if (sum(counts@x %% 1 != 0, na.rm = TRUE) > 0) {
    warning("计数矩阵中存在非整数数值，请确认是原始UMI计数！")
  }
  
  colnames(sc@meta.data)
  table(sc$cellType)
  # B_cell      Cancer_cell Endothelial_cell       Fibroblast 
  # 13995            65885            12464            37439 
  # Mast_cell     Myeloid_cell              pDC           T_cell 
  # 1145            23146              938            71623 
  
  # 2. 准备细胞类型（与原有代码一致）
  Idents(sc) <- "cellType"
  cluster <- as.factor(sc$cellType)
  names(cluster) <- colnames(sc)
  
  # 3. 重构nUMI：确保为整数类型且数值为整数（核心修复）
  ## 方法1：从counts矩阵重新计算nUMI（最可靠，避免Seurat对象的数值误差）
  nUMI <- colSums(counts)  # 直接计算每个细胞的UMI总数，保证为整数
  names(nUMI) <- colnames(counts)  # 确保列名匹配
  
  ## 方法2：若坚持使用sc$nCount_RNA，强制转换为整数
  # nUMI <- as.integer(round(sc$nCount_RNA))
  # names(nUMI) <- colnames(sc)
  
  # 验证nUMI的整数属性（关键：必须通过）
  stopifnot(all(nUMI %% 1 == 0, na.rm = TRUE))
  stopifnot(is.integer(nUMI) || is.numeric(nUMI) && all(nUMI == as.integer(nUMI)))
  
  # 4. 构建Reference对象
  reference <- Reference(counts, cluster, nUMI)
}
if(T){
  # 加载必要包
  library(Seurat)
  library(RCTD)
  library(Matrix)
  sc <-readRDS("~/R/sc_bc/out_1/mydata_CSS_20230925.rds")
  
  # 1. 提取原始计数矩阵并确保数值为整数（dgCMatrix保留double类型）
  counts <- sc@assays$RNA@counts
  # 处理浮点误差（如1.0000001 → 1.0）
  counts@x <- round(counts@x)
  # 验证counts数值为整数（可选）
  if (sum(counts@x %% 1 != 0, na.rm = TRUE) > 0) {
    warning("计数矩阵中存在非整数数值，请确认是原始UMI计数！")
  }
  
  colnames(sc@meta.data)
  table(sc$celltype_major)
  # B-cells              CAFs Cancer Epithelial       Endothelial 
  # 3181              6511             24489              7546 
  # Myeloid Normal Epithelial      Plasmablasts               PVL 
  # 9553              4338              3390              5295 
  # T-cells 
  # 35001 
  
  # 2. 准备细胞类型（与原有代码一致）
  Idents(sc) <- "celltype_major"
  cluster <- as.factor(sc$celltype_major)
  names(cluster) <- colnames(sc)
  
  # 3. 重构nUMI：确保为整数类型且数值为整数（核心修复）
  ## 方法1：从counts矩阵重新计算nUMI（最可靠，避免Seurat对象的数值误差）
  nUMI <- colSums(counts)  # 直接计算每个细胞的UMI总数，保证为整数
  names(nUMI) <- colnames(counts)  # 确保列名匹配
  
  ## 方法2：若坚持使用sc$nCount_RNA，强制转换为整数
  # nUMI <- as.integer(round(sc$nCount_RNA))
  # names(nUMI) <- colnames(sc)
  
  # 验证nUMI的整数属性（关键：必须通过）
  stopifnot(all(nUMI %% 1 == 0, na.rm = TRUE))
  stopifnot(is.integer(nUMI) || is.numeric(nUMI) && all(nUMI == as.integer(nUMI)))
  
  # 4. 构建Reference对象
  reference <- Reference(counts, cluster, nUMI)
  
  ##run RCTD
  library(BiocParallel)
  register(MulticoreParam(workers = 16, progressbar = TRUE))
  
  ###run RCTD in visium list
  spatial_list_RCTD<-list() #save new visium in list
  for (i in 1:length(names(visium_list))  ){ #测试 #1:length(names(visium_list))
    print(i);print(names(visium_list)[i])
    visium<-visium_list[[i]] #list visium_list[1]
    
    coords<-GetTissueCoordinates(visium ) #, image = names(visium_list)[1]
    colnames(coords) <- c('x','y')
    
    counts <- visium@assays$Spatial@counts
    # Create SpatialRNA object
    query <- SpatialRNA(coords, counts, colSums(counts))
    # Run RCTD
    print("------------start RCTD ---------")
    RCTD <- create.RCTD(query, reference, max_cores = 16)#32 slow  #8
    RCTD.full <- run.RCTD(RCTD, doublet_mode = "full") #"doublet"
    barcodes <- colnames(RCTD.full@spatialRNA@counts)
    weights <- RCTD.full@results$weights
    #norm_weights <- normalize_weights(weights)#could not find function "normalize_weights"
    # Normalize weights
    results <- RCTD.full@results
    norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
    norm_weights<-as.data.frame( norm_weights )
    colnames(norm_weights)[3] <- "Cancer_cell"
    colnames(norm_weights)[6] <-"Normal_Epi"
    # 
    #norm_weights 添加注释列 得分手动注释
    
    ##cancer > mean
    norm_weights_1 <-norm_weights
    norm_weights_1$Type <-ifelse(norm_weights_1$Cancer_cell > mean(norm_weights_1$Cancer_cell),"Cancer_cell","")
    table(norm_weights_1$Type)
    df<-apply(norm_weights[, !colnames(norm_weights) %in% c("Type", "Cancer_cell")] ,#norm_weights[c(1,3:5)],
              1,function(x){which.max(x)} )#最大的行号-除cancer外
    df<-as.data.frame(df)
    norm_weights_1$index <-df[,1]
    norm_weights_1$Type <-ifelse(norm_weights_1$Type=="Cancer_cell","Cancer_cell",
                                 colnames(norm_weights_1)[norm_weights_1$index] )
    table(norm_weights_1$Type)
    ##直接最大得分-细胞类型 cell type
    df1<-as.data.frame(
      apply(norm_weights,1,function(x){which.max(x)} )
    )
    norm_weights_1$index1 <-df1[,1]
    norm_weights_1$Type_1 <-colnames(norm_weights_1)[norm_weights_1$index1]
    #colnames(norm_weights_1)
    table(norm_weights_1$Type_1)
    #colnames(norm_weights_1)[5]<-"T&NK"
    
    cell_type_names <- colnames(norm_weights)#RCTD.full@cell_type_info$info[[2]] #list of cell type names
    spatialRNA <- RCTD.full@spatialRNA
    resultsdir <- paste0(folder_path,"/",Sys.Date(),"_","RCTD_plot_",names(visium_list[i]) ) 
    # 'RCTD_Plots' 输出文件夹## you may change this to a more accessible directory on your computer.
    dir.create(resultsdir) #names(visium_list[i])
    print("------------finish RCTD ---------")
    write.csv(norm_weights_1,paste0(resultsdir,"/",Sys.Date(),"_",names(visium_list)[i],"_norm_weights_1.csv"))
    
    # Add RCTD results to Seurat object
    print(table(norm_weights_1$Type))#print(colnames(results[["results_df"]]))
    print(table(norm_weights_1$Type_1))#
    #sc <- AddMetaData(sc, metadata = RCTD.full@results$results_df)#
    visium <- AddMetaData(visium, metadata = norm_weights_1[,c("Type","Type_1")] )#RCTD.full@results$results_df
      #for 循环 list 保存！
    print("------------plot ---------")
    # make the plots
    # Plots the confident weights for each cell type as in full_mode (saved as
    # 'results/cell_type_weights_unthreshold.pdf')
    plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)
    # Plots all weights for each cell type as in full_mode. (saved as
    # 'results/cell_type_weights.pdf')
    plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights)
    # Plots the weights for each cell type as in doublet_mode. (saved as
    # 'results/cell_type_weights_doublets.pdf')
    plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet,
                         results$results_df)
    # Plots the number of confident pixels of each cell type in 'full_mode'. (saved as
    # 'results/cell_type_occur.pdf')
    plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
    print(paste0("------Finish ",names(visium_list)[i]," ---------"))
    #save
    spatial_list_RCTD[[names(visium_list)[i]]]<-visium
    rm(visium)
    
  }
  qs::qsave(spatial_list_RCTD,paste0(folder_path,"/TNBC_Visium_GSE210616_list_RCTD_full_",Sys.Date(),"-anno-GES176078.qs") ) #保存 full or doublet
}

##1.2 load visium list----
visium_list<-qread("./data/TNBC_Visium_GSE210616_list.qs")

##1.3 run RCTD----
if (T){
  #visium_list<- spatial_list_cluster#
  
  library(BiocParallel)
  register(MulticoreParam(workers = 16, progressbar = TRUE))
  
  ###run RCTD in visium list
  spatial_list_RCTD<-list() #save new visium in list
  for (i in 1:length(names(visium_list))  ){ #测试 #1:length(names(visium_list))
    print(i);print(names(visium_list)[i])
    visium<-visium_list[[i]] #list visium_list[1]
    
    coords<-GetTissueCoordinates(visium ) #, image = names(visium_list)[1]
    colnames(coords) <- c('x','y')
    
    counts <- visium@assays$Spatial@counts
    # Create SpatialRNA object
    query <- SpatialRNA(coords, counts, colSums(counts))
    # Run RCTD
    print("------------start RCTD ---------")
    RCTD <- create.RCTD(query, reference, max_cores = 8)#32 slow  #8
    RCTD.full <- run.RCTD(RCTD, doublet_mode = "full") #"doublet"
    barcodes <- colnames(RCTD.full@spatialRNA@counts)
    weights <- RCTD.full@results$weights
    #norm_weights <- normalize_weights(weights)#could not find function "normalize_weights"
    # Normalize weights
    results <- RCTD.full@results
    norm_weights = sweep(results$weights, 1, rowSums(results$weights), '/') 
    norm_weights<-as.data.frame( norm_weights )

    # 
    #norm_weights 添加注释列 得分手动注释
    
    ##cancer > mean
    norm_weights_1 <-norm_weights
    norm_weights_1$Type <-ifelse(norm_weights_1$Cancer_cell > mean(norm_weights_1$Cancer_cell),"Cancer_cell","")
    table(norm_weights_1$Type)
    df<-apply(norm_weights[, !colnames(norm_weights) %in% c("Type", "Cancer_cell")] ,#norm_weights[c(1,3:5)],
               1,function(x){which.max(x)} )#最大的行号-除cancer外
    df<-as.data.frame(df)
    norm_weights_1$index <-df[,1]
    norm_weights_1$Type <-ifelse(norm_weights_1$Type=="Cancer_cell","Cancer_cell",
                                 colnames(norm_weights_1)[norm_weights_1$index] )
    table(norm_weights_1$Type)
    ##直接最大得分-细胞类型 cell type
    df1<-as.data.frame(
      apply(norm_weights,1,function(x){which.max(x)} )
    )
    norm_weights_1$index1 <-df1[,1]
    norm_weights_1$Type_1 <-colnames(norm_weights_1)[norm_weights_1$index1]
    #colnames(norm_weights_1)
    table(norm_weights_1$Type_1)
    #colnames(norm_weights_1)[5]<-"T&NK"
    
    cell_type_names <- RCTD.full@cell_type_info$info[[2]] #list of cell type names
    spatialRNA <- RCTD.full@spatialRNA
    resultsdir <- paste0(folder_path,"/",Sys.Date(),"_","RCTD_plot_",names(visium_list[i]) ) 
    # 'RCTD_Plots' 输出文件夹## you may change this to a more accessible directory on your computer.
    dir.create(resultsdir) #names(visium_list[i])
    print("------------finish RCTD ---------")
    write.csv(norm_weights_1,paste0(resultsdir,"/",Sys.Date(),"_",names(visium_list)[i],"_norm_weights_1.csv"))
    
    # Add RCTD results to Seurat object
    print(table(norm_weights_1$Type))#print(colnames(results[["results_df"]]))
    print(table(norm_weights_1$Type_1))#
    #sc <- AddMetaData(sc, metadata = RCTD.full@results$results_df)#
    visium <- AddMetaData(visium, metadata = norm_weights_1[,c("Type","Type_1")] )#RCTD.full@results$results_df
    #print(table(visium@meta.data[["first_type"]]))
    #print(table(visium@meta.data[["second_type"]]))
    
    # Add RCTD results as a new assay
    #print("------------add RCTD norm matrix---------")
    #visium[["RCTD"]] <- CreateAssayObject(counts = t(as.matrix(norm_weights)))
    #Error: Cannot add a different number of cells than already present
    
    # if (length(sc@assays$rctd_full@key) == 0) {
    #   sc@assays$rctd_full@key <- "rctd_full_"
    # }
    #qs::qsave(visium,paste0(folder_path,Sys.Date(),"_","RCTD_test.qs") ) #保存
    #for 循环 list 保存！
    print("------------plot ---------")
    # make the plots
    # Plots the confident weights for each cell type as in full_mode (saved as
    # 'results/cell_type_weights_unthreshold.pdf')
    plot_weights(cell_type_names, spatialRNA, resultsdir, norm_weights)
    # Plots all weights for each cell type as in full_mode. (saved as
    # 'results/cell_type_weights.pdf')
    plot_weights_unthreshold(cell_type_names, spatialRNA, resultsdir, norm_weights)
    # Plots the weights for each cell type as in doublet_mode. (saved as
    # 'results/cell_type_weights_doublets.pdf')
    plot_weights_doublet(cell_type_names, spatialRNA, resultsdir, results$weights_doublet,
                         results$results_df)
    # Plots the number of confident pixels of each cell type in 'full_mode'. (saved as
    # 'results/cell_type_occur.pdf')
    plot_cond_occur(cell_type_names, resultsdir, norm_weights, spatialRNA)
    # makes a map of all cell types, (saved as
    # 'results/all_cell_types.pdf')
    # plot_all_cell_types(norm_weights_1[6],#results$results_df, 
    #                     spatialRNA@coords, cell_type_names, resultsdir)
    
    
    
    # doublets
    #obtain a dataframe of only doublets
    #doublets <- results$results_df[results$results_df$spot_class == "doublet_certain",]
    # Plots all doublets in space (saved as
    # 'results/all_doublets.pdf')
    #plot_doublets(spatialRNA, doublets, resultsdir, cell_type_names)
    # Plots all doublets in space for each cell type (saved as
    # 'results/all_doublets_type.pdf')
    #plot_doublets_type(spatialRNA, doublets, resultsdir, cell_type_names)
    # a table of frequency of doublet pairs
    #doub_occur <- table(doublets$second_type, doublets$first_type)
    # Plots a stacked bar plot of doublet ocurrences (saved as
    # 'results/doublet_stacked_bar.pdf')
    #plot_doub_occur_stack(doub_occur, resultsdir, cell_type_names)
    print(paste0("------Finish ",names(visium_list)[i]," ---------"))
    #save
    spatial_list_RCTD[[names(visium_list)[i]]]<-visium
    rm(visium)
    
  }
  
  qs::qsave(spatial_list_RCTD,paste0(folder_path,"/TNBC_Visium_GSE210616_list_RCTD_full_",Sys.Date(),".qs") ) #保存 full or doublet
  
}

#2.0 Cor----

#spatial_list_RCTD <- qread("./out/TNBC_Visium_GSE210616_list_RCTD_full_2026-04-01.qs")
#load GES176078 annot
spatial_list_RCTD <-qread("./out/TNBC_Visium_GSE210616_list_RCTD_full_2026-04-02-anno-GES176078.qs")

genelist <- c("CD74")  # 你的基因列表
mygene <- "MIF"  # 要分析的特定基因
methods <- c("pearson", "spearman")  # 相关性方法

##2.1 all cell Cor----
# 初始化结果框
cor_results <- data.frame(
  Sample = character(),
  Group_col = character(),
  Method = character(),
  Correlation = numeric(),
  P_value = numeric(),
  stringsAsFactors = FALSE
)

# 设定方法
methods <- c("pearson", "spearman")

# 主循环
for (i in names(spatial_list_RCTD)) {
  message("Processing: ", i)
  sce <- spatial_list_RCTD[[i]]
  
  # ===================== 基因匹配（自动取交集，防错）
  all_genes <- rownames(sce)
  gene1 <- intersect(genelist, all_genes)
  gene2 <- intersect(mygene, all_genes)
  
  if(length(gene1)==0 | length(gene2)==0) {
    message("Skip ", i, " (gene not found)")
    next
  }
  
  # ===================== 提取表达到 meta.data
  genes_use <- c(gene1, gene2)
  sce@meta.data[, genes_use] <- t(as.matrix(
    GetAssayData(sce, assay = "SCT", slot = "data")[genes_use, ]
  ))
  
  # ===================== 遍历 Type / Type_1
  for (k in c("Type", "Type_1")) {
    message("  Group: ", k)
    
    # 取数据
    df_cor <- sce@meta.data[, c(k, gene1, gene2)]
    colnames(df_cor) <- c("celltype", "gene", "Score")
    
    # 过滤空数据
    if(nrow(df_cor) < 3) {
      message("  Skip (too few cells)")
      next
    }
    
    # ===================== 计算相关
    cor_p <- cor.test(df_cor$gene, df_cor$Score, method = "pearson")
    cor_s <- cor.test(df_cor$gene, df_cor$Score, method = "spearman")
    
    # ===================== 遍历方法并存储结果
    for (j in methods) {
      message("    Method: ", j)
      
      if(j == "pearson") {
        r <- cor_p$estimate
        pval <- cor_p$p.value
      } else {
        r <- cor_s$estimate
        pval <- cor_s$p.value
      }
      
      # 安全绑定结果
      cor_results <- rbind(cor_results, data.frame(
        Sample = i,
        Group_col = k,
        Method = j,
        Correlation = as.numeric(r),
        P_value = as.numeric(pval)
      ))
      
      # ===================== 绘图
      p <- ggplot(df_cor, aes(x = gene, y = Score)) +
        labs(x = gene1, y = gene2, title = i) +
        geom_point(color = "skyblue", size = 1.5, alpha = 0.5) +
        geom_smooth(method = 'lm', size = 0.5, color = "pink") +
        stat_cor(method = j, size = 6, color = "red") +
        scale_colour_manual(values = "black") +  # 按需调整颜色
        theme_bw() +
        theme(
          axis.text.x = element_text(size = 14),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_text(size = 14),
          axis.title.y = element_text(size = 14),
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(size = 14)
        )
      
      print(p)
      
      # 保存
      ggsave(
        filename = paste0(folder_path, "/Cor_", i, "_", k, "_", gene1, "_", gene2, "_", j, ".pdf"),
        plot = p, width = 3, height = 3,family="serif"
      )
    }
  }
}

# 最后查看结果
head(cor_results)

#all cell Cor
print(cor_results)
write.csv(cor_results, paste0(folder_path, "Visium_Correlation_Results.csv"), row.names = FALSE)

#plot all Cor

cor_results <- cor_results[cor_results$Group_col=="Type",]#
for (j in methods) {
  print(j)
  #outTab <- na.omit(cor_results)
  #outTab <- outTab[outTab$Method==j,]
  
  for (i in unique(cor_results$Group_col) ){
    print(i)
    outTab <- na.omit(cor_results)
    outTab <- outTab[outTab$Group_col ==i & outTab$Method==j,]
    
    colnames(outTab)[4:5]<-c("Cor","pvalue")
    outTab[4:5]<- as.data.frame(sapply(outTab[4:5],as.numeric )) #转为数值
    outTab$'-log10(pvalue)' <- -log10(outTab$pvalue) #-log10(0.05) = 1.30
    outTab$Sign <-ifelse(outTab$pvalue <0.05,"*"," ")
    #outTab$Cell <-gsub("T_cells_","",outTab$Cell) gsub("c./*_","",outTab$Cell)
    
    p <-ggdotchart(outTab, x = "Sample",y = "Cor",
                   color ="-log10(pvalue)",
                   sorting = "descending",#排序“ascending”, “descending”, “none”
                   add.params = list(color = "lightgray"), #画棒棒
                   add = "segments",#画棒棒
                   dot.size = "Cor",#-log10(pvalue), #"cor"
                   xlab="",ylab="Correlation Coefficient", 
                   #label.select = list(criteria = " `pvalue` < 0.05 "),
                   #label = T, #小数点太多了
                   ggtheme = theme_bw(base_size = 12)) + 
      #coord_flip()+  #翻转坐标轴
      scale_color_gradient(low = "skyblue",high="pink")+  #颜色范围orange
      font("x.text", size = 12, vjust = 0.1,angle = 0)+  #X坐标轴
      #font("y.text",size = 13,face = "bold")+  #坐标加粗
      #geom_hline(yintercept=c(-0.4,0.4), linetype="dashed",size=0.1)+ #添加虚线
      #scale_color_continuous(low = "skyblue",high="orange")+
      #scale_y_continuous(limits = c(-0.5,0.5))+ #c(-0.5,1.0)
      geom_hline(yintercept=c(0),size=0.1)+
      geom_text(size=6, aes( label = Sign),vjust=1) +#hjust = 1,
      theme(axis.text.y = element_text(size = 12 ),panel.grid = element_blank(),
            axis.text.x = element_text(angle = 90, vjust = 0.5) )  #,legend.position = c(0.7,0.7)
    
    p
    ggsave( paste0(folder_path, "/Cor_",mygene,"_all_",genelist,"_",mygene,"_",j,".pdf"),
            family="serif",width = 8,height = 3)
    
  }
  
}
#negative Cor

##2.2 cell exp mean----
#MIF-tumor CD74-T 患者平均值 分组？正相关JSP
if(T){
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
  genelist_1<- list("Jak-STAT_pathway"=genes,"Core_Jak-STAT_pathway"=genes_1)
}
pair_expr <- data.frame()

# 循环处理样本
for (i in names(spatial_list_RCTD) ) {
  message("Processing: ", i)
  sceList <- spatial_list_RCTD[[i]]
  ##add score
  library(BiocParallel)
  register(MulticoreParam(workers = 24, progressbar = TRUE))#32
  #score
  sceList <- AddModuleScore(object = sceList,#sce,
                            features = genelist_1,#gene,
                            name = names(genelist_1)#"score" #ctrl = 100,
  ) #
  colnames(sceList@meta.data)[10] <- 'Score' 
  
  
  # ===================== 1. 提取癌症细胞（必须有）
  target_cells_1 <- colnames(sceList)[sceList$Type_1 %in% c("Cancer_cell","T_cell")]
  expr_mat_1 <- GetAssayData(sceList, assay = "SCT", slot = "data")[, target_cells_1, drop=F]
  
  g1 <- gene1
  g2 <- gene2
  
  # 癌症基因 g1（一定有长度）
  g1_expr <- if (g1 %in% rownames(expr_mat_1)) expr_mat_1[g1, ] else rep(0, length(target_cells_1))
  g2_expr <- if (g2 %in% rownames(expr_mat_1)) expr_mat_1[g1, ] else rep(0, length(target_cells_1))
  
  avg_expr <- (g1_expr + g2_expr) / 2
  temp_df <- data.frame(
    Patient = i,
    id = colnames(expr_mat_1),
    Score = sceList$Score[colnames(expr_mat_1)] ,
    g1 = as.numeric(g1_expr),
    g2 = as.numeric(g2_expr),
    Avg_Expr = as.numeric(avg_expr)
  )
  
  pair_expr <- rbind(pair_expr, temp_df)
}

##Score with mean MIF-CD74
library(ggplot2);library(ggpubr);library(ggsci)
ggplot(pair_expr, aes(x=Score, y=Avg_Expr)) + 
  #xlim(-20,15) + ylim(-15,10) +
  labs(x = "JSP score", y = "MIF-CD74 mean expression") +
  geom_point(aes(color=Patient),size = 2,alpha=0.3) +
  geom_smooth(method ='lm', size=0.5,color="skyblue") +
  stat_cor(method = "spearman",size = 6,color="red") +
  #scale_colour_manual(values = colors ) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave(file=paste0(folder_path,'/Cor_score_MIF-CD74.pdf'),width = 6,height = 3.5,family="serif")

#calculate mean by patient
library(dplyr)
patient_avg <- pair_expr %>%
  group_by(Patient) %>%
  summarise(
    across(c(Score, g1, g2, Avg_Expr), mean),
    .groups = "drop"
  )

my43colors <- c(
  "#0072B5", "#E18727", "#F37C95", "#20854E", "#D6E7A3", "#57C3F3", "#E5D2DD",
  "#53A85F", "#F1BB72", "#F3B1A0", "#476D87", "#E95C59", "#E59CC4", "#AB3282",
  "#23452F", "#BD956A", "#8C549C", "#585658", "#9FA3A8", "#E0D4CA", "#5F3D69",
  "#C5DEBA", "#58A4C3", "#E4C755", "#F7F398", "#AA9A59", "#E63863", "#E39A35",
  "#C1E6F3", "#6778AE", "#91D0BE", "#B53E2B", "#712820", "#DCC1DD", "#CCE0F5",
  "#CCC9E6", "#625D9E", "#68A180", "#3A6963", "#968175",
  "#7FBC41", "#FFC107", "#8E24AA"
)
ggplot(patient_avg, aes(x=Score, y=Avg_Expr)) + 
  #xlim(-20,15) + ylim(-15,10) +
  labs(x = "JSP score", y = "MIF-CD74 mean expression") +
  geom_point(aes(color=Patient),size = 2,alpha=0.7) +
  geom_smooth(method ='lm', size=0.5,color="skyblue") +
  stat_cor(method = "spearman",size = 6,color="red") +
  scale_colour_manual(values = my43colors ) +
  theme_bw() + 
  theme(axis.text.x = element_text(size = 16), 
        axis.text.y = element_text(size = 16), 
        axis.title.x = element_text(size = 18), 
        axis.title.y = element_text(size = 18),
        legend.title = element_blank(),
        legend.text = element_text(size = 12))
ggsave(file=paste0(folder_path,'/Cor_score_MIF-CD74-patient.pdf'),width = 6,height = 3.5,family="serif")


for (i in names(spatial_list_RCTD) ) {
  message("Processing: ", i)
  sceList <- spatial_list_RCTD[[i]]
  print(table(sceList$Type_1))#Type no T cell
}

#3.0 distance and co-localization----
##3.1 co-localization----
#MIF+tumor CD74+ T cells 
j="Spectral"
spatial_list_RCTD_1 <- list()
for (i in names(spatial_list_RCTD) ) {
  message("Processing: ", i)
  sce <- spatial_list_RCTD[[i]]
  ##add score
  library(BiocParallel)
  register(MulticoreParam(workers = 24, progressbar = TRUE))#32
  #score
  sce <- AddModuleScore(object = sce,#sce,
                            features = genelist_1,#gene,
                            name = names(genelist_1)#"score" #ctrl = 100,
  ) #
  colnames(sce@meta.data)[10] <- 'Score'
  ##add gene
  genes_use <- c(gene1, gene2)
  sce@meta.data[, genes_use] <- t(as.matrix(
    GetAssayData(sce, assay = "SCT", slot = "data")[genes_use, ]
  ))
  ##save
  spatial_list_RCTD_1[[i]]<-sce
}
#qs::qsave(spatial_list_RCTD_1,paste0(folder_path,"/TNBC_Visium_GSE210616_list_RCTD_full_",Sys.Date(),"-anno-GES176078.qs") ) #add score


##spatial plot
celltype_colors <- c(
  "B-cells"      = "#476D87",    
  "CAFs"         = "#E18727",  
  "Plasmablasts"  = "#F37C95",   
  "Endothelial"  =   "#0072B5" , 
  "Myeloid"      =   "#20854E", 
  "Normal_Epi"   = "#AB3282",   
  "Cancer_cell"= "#57C3F3",
  "PVL"          = "#53A85F",   
  "T-cells"      = "#E95C59"    
)


for (i in names(spatial_list_RCTD_1) ){
  print(i)
  sce <- spatial_list_RCTD_1[[i]]
  #only HE
  SpatialDimPlot(sce,group.by= "Type",label.size=3,repel=T, label = F,# 报错
                 stroke=NA, #cols = my36colors,#描边 label = T报错
                 #image.alpha=0,
                 pt.size.factor = 0) + NoLegend()+
    #scale_fill_brewer(palette = "Paired" )+
    scale_fill_manual(values = celltype_colors)+
    ggtitle(i)
  ggsave(paste0(folder_path,"/SpatialDimPlot_",i,"_HE.pdf"),width = 4,height = 3,family="serif")
  
  
  for (j in c("Type","Type_1") ){
    print(j)
    #-HE
    SpatialDimPlot(sce,group.by= j,label.size=3,repel=T, label = T,# 报错
                   stroke=NA, #cols = my36colors,#描边 label = T报错
                   image.alpha=0,
                   pt.size.factor = 5.5) + #NoLegend()+
      #scale_fill_brewer(palette = "Paired" )+
      scale_fill_manual(values = celltype_colors)+
      ggtitle(i)
    ggsave(paste0(folder_path,"/SpatialDimPlot_",i,"_",j,".pdf"),width = 4,height = 3,family="serif")
    
  }

  
  for (K in c("MIF","CD74") ){
    print(K)
    SpatialFeaturePlot(
      object = sce,
      features = K,#"PCSK9"  # 使用原始基因名确保匹配
      pt.size.factor = 5.5,    # 点的大小（可根据需求调整）
      stroke = NA,             # 去除点的描边
      image.alpha = 0          # 保留空间背景图像（设为0则隐藏）
    ) +labs(title = i)+#ggtitle(paste0(i))+
      scale_fill_gradientn(
        colors = rev(brewer.pal(n = 8, name = "Spectral")),  # c("Spectral") ){ #"Reds","PuBuGn","GnBu","BrBG","PiYG","RdBu",
        na.value = "grey90"
      )
    ggsave(paste0(folder_path,"/SpatialFeaturePlot_",i,"_",K,".pdf"),width = 4,height = 3,family="serif")
    
    
  }
  
}


##3.2 distance----

###run demo----

if(F){
  library(tidyverse)
  library(phenoptr)
  #devtools::install_github("akoyabio/phenoptr")
  #install.packages("~/R/my_projects/JAK_STAT_BC/13.0_TNBC/data/phenoptr-main", repos = NULL, type = "source")
  #https://github.com/akoyabio/phenoptr/blob/main/data/sample_cell_seg_data.rda
  
  

  #load("~/R/my_projects/JAK_STAT_BC/13.0_TNBC/data/sample_cell_seg_data.rda")
  csd <- sample_cell_seg_data
  csd <- csd %>% filter(Phenotype!='other')
  table(csd$Phenotype)
  # CD68+   CD8+    CK+ FoxP3+ 
  #   417    228   2257    228 
  
  csd <- csd[4:7]
  colnames(csd )
  
  distances <- find_nearest_distance(csd)
  glimpse(distances)
  
  csd_with_distance <- bind_cols(csd, distances)
  
  # merged_with_distance <- merged %>%
  #   group_by(`Sample Name`) %>%
  #   do(bind_cols(., find_nearest_distance(.)))
  
  ####Analyzing nearest neighbor distances
  csd_with_distance %>% group_by(Phenotype) %>% 
    select(Phenotype, starts_with('Distance to')) %>% 
    summarize_all(~round(mean(.), 1))
  
  colnames(csd_with_distance)
  ggplot(csd_with_distance, aes(`Distance to CD8+`, color=Phenotype)) +
    geom_density(size=0.5)+
    scale_color_manual(values = my36colors)+theme_bw()
  
  ggplot(csd_with_distance, aes(`Distance to CD68+`, color=Phenotype)) +
    geom_density(size=0.5)+
    scale_color_manual(values = my36colors)+theme_bw()
  
}


###my distance----
celltype_colors_1 <- c(
  "B-cells"      = "#476D87",    
  "CAFs"         = "#E18727",  
  "Plasmablasts"  = "#F37C95",   
  "Endothelial"  =   "#0072B5" , 
  "Myeloid"      =   "#20854E", 
  "Normal_Epi"   =  "#E95C59", 
  "Cancer_cell"= "#57C3F3",
  "PVL"          = "#53A85F",   
  "T-cells"      = "#AB3282", 
  "MIF+ Cancer_cell"="red"
)


for (i in names(spatial_list_RCTD_1) ){
  sce <- spatial_list_RCTD_1[[i]]
  if(!"T-cells" %in% names(table(sce$Type_1)) || table(sce$Type_1)["T-cells"]==0) next
  
  ###
  print(i)
  # [1] "093A"
  # [1] "093B"
  # [1] "093C"
  # [1] "094A"
  # [1] "095A"
  # [1] "095B"
  # [1] "119C"
  # [1] "120C"
  # [1] "120D"
  # [1] "395D"
  # [1] "397C"
  # [1] "397D"
  ##count
  # sce@meta.data[, genes_use] <- t(as.matrix(
  #   GetAssayData(sce, assay = "Spatial", slot = "data")[genes_use, ]
  # ))
  sce$Cell_Annotation <- ifelse( (sce$Type_1 == "Cancer_cell") & (sce$MIF > mean(sce$MIF)),"MIF+ Cancer_cell", 
                                 sce$Type_1)
  table(sce$Cell_Annotation)
  
  SpatialDimPlot(sce,group.by= "Cell_Annotation",label.size=3,repel=T, label = T,# 报错
                 stroke=NA, #cols = my36colors,#描边 label = T报错
                 image.alpha=0,
                 pt.size.factor = 5.5) + #NoLegend()+
    #scale_fill_brewer(palette = "Paired" )+
    scale_fill_manual(values = celltype_colors_1)+
    ggtitle(i)
  ggsave(paste0(folder_path,"/SpatialDimPlot_",i,"_Cell_Annotation.pdf"),width = 4,height = 3,family="serif")
  SpatialDimPlot(sce,group.by= "Cell_Annotation",label.size=3,label = F,# 报错
                 stroke=NA, #cols = my36colors,#描边 label = T报错
                 image.alpha=0,
                 pt.size.factor = 5.5) + #NoLegend()+
    #scale_fill_brewer(palette = "Paired" )+
    scale_fill_manual(values = celltype_colors_1)+
    ggtitle(i)
  ggsave(paste0(folder_path,"/SpatialDimPlot_",i,"_Cell_Annotation-nolabel.pdf"),width = 4,height = 3,family="serif")
  
  
  # 提取目标信息并生成标准数据框
  library(tidyverse)
  library(phenoptr)
  cell_info <- data.frame(
    "Phenotype" = sce@meta.data$Cell_Annotation,      # Type重命名为Phenotype
    "Cell ID" = rownames(sce@meta.data),  # 细胞ID
    "Cell X Position" = sce@images$slice1@coordinates$col,  # X坐标
    "Cell Y Position" = sce@images$slice1@coordinates$row   # Y坐标
  )
  colnames(cell_info) <-c( "Phenotype", "Cell ID","Cell X Position" ,"Cell Y Position")
  
  distances <- find_nearest_distance(cell_info)
  #glimpse(distances)
  csd_with_distance <- bind_cols(cell_info, distances)
  csd_with_distance %>% group_by(Phenotype) %>% 
    select(Phenotype, starts_with('Distance to')) %>% 
    summarize_all(~round(mean(.), 1))
  
  ggplot(csd_with_distance, aes(`Distance to MIF+ Cancer_cell`, color=Phenotype)) +
    geom_density(size=0.3)+
    scale_color_manual(values = celltype_colors)+theme_bw(base_size = 12)
  ggsave(paste0(folder_path,"/Distance_",gene2,"+_all_cells_",i,".pdf"),
         width = 4.5,height = 2.5,family="serif")
  
  # ggplot(csd_with_distance[csd_with_distance$Phenotype=="T-cells",], aes(`Distance to MIF+ Cancer_cell`, color=Phenotype)) +
  #   geom_density(size=0.5)+
  #   scale_color_manual(values = celltype_colors)+theme_bw(base_size = 12)+
  # ggplot(csd_with_distance[csd_with_distance$Phenotype=="T-cells",], aes(`Distance to Cancer_cell`, color=Phenotype)) +
  #   geom_density(size=0.5)+
  #   scale_color_manual(values = celltype_colors)+theme_bw(base_size = 12)
  colnames(csd_with_distance)
  
  if(!"Distance to Cancer_cell" %in% colnames(csd_with_distance)) next
  if(T){
    ggplot(csd_with_distance, aes(`Distance to Cancer_cell`, color=Phenotype)) +
      geom_density(size=0.3)+
      scale_color_manual(values = celltype_colors)+theme_bw(base_size = 12)
    ggsave(paste0(folder_path,"/Distance_",gene2,"-_all_cells_",i,".pdf"),
           width = 4.5,height = 2.5,family="serif")
    
    df1 <-csd_with_distance[csd_with_distance$Phenotype=="T-cells",c("Phenotype","Distance to Cancer_cell")]
    colnames(df1)[2]<- c("Distance to Cancer_cell")
    df1$Group <- "MIF- Cancer_cell"
    df2 <-csd_with_distance[csd_with_distance$Phenotype=="T-cells",c("Phenotype","Distance to MIF+ Cancer_cell")]
    colnames(df2)[2]<- c("Distance to Cancer_cell")
    df2$Group <- "MIF+ Cancer_cell"
    df <- rbind(df1,df2)
    
    ggplot(df, aes(`Distance to Cancer_cell`, color=Group)) +
      geom_density(size=0.5)+labs(y="T cell density",title = i)+
      scale_color_manual(values = c("MIF+ Cancer_cell"="pink","MIF- Cancer_cell"="skyblue") )+
      theme_bw(base_size = 12)+theme(panel.grid.minor = element_blank() )
    ggsave(paste0(folder_path,"/Distance_",gene2,"+_T_cells_",i,".pdf"),
           width = 4.5,height = 2.5,family="serif")
  }
  
  
  ##
}

##annotate MIF+ Cancer cell
#i="093A"

# colnames(sce@meta.data)
# #table(sce$Type)
# table(sce$Type_1)
#plot(sce$MIF)#"Spatial"

