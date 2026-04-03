
#~/R/my_projects/JAK_STAT_BC/13.0_TNBC

if(T){
  rm(list = ls())
  setwd("~/R/my_projects/JAK_STAT_BC/13.0_TNBC")
  #输出文件夹!!!!!!!
  program_name="out"
  
  folder_path <- program_name#"./out" # 检查文件夹是否存在
  if (!dir.exists(folder_path)) {   # 如果文件夹不存在，则创建它
    dir.create(folder_path)
    print(paste("Folder created at", folder_path))
  } else {
    print(paste("Folder already exists at", folder_path))
  }
  
  my36colors <-c("#0072B5","#E18727","#F37C95","#20854E", '#D6E7A3', '#57C3F3', '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 
  
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
  )#geneset
  genes_1 <- c(
    "JAK1", "JAK2", "JAK3", 
    "STAT1", "STAT2", "STAT3", 
    "STAT4", "STAT5A", "STAT5B", 
    "STAT6"
  )#core gene
  interest_gene <-"STAT4"
  
  my_cancer<-c("BRCA")
  
  genelist<- list("Jak-STAT_pathway"=genes,"Core_Jak-STAT_pathway"=genes_1
                  #"Proliferation"= c("BIRC5", "CCNB1", "CDC20", "NUF2", "CEP55", "NDC80", "MKI67", "PTTG1", "RRM2", "TYMS","UBE2C")
  ) #https://www.nature.com/articles/s41588-021-00911-1
  input_gene= c("Jak-STAT_pathway","Core_Jak-STAT_pathway") #genes#c(genes,"H3K36me","5mC" )
  project = "Jak_STAT_BC"
}


#1.0 FUSCC-TNNBC (n=89)----
##load----
#article:https://www.sciencedirect.com/science/article/pii/S1535610819300960?via%3Dihub#app2
#matrix: https://www.biosino.org/node/analysis/detail/OEZ000398
#400+ mrna   67 RFS=1

load("~/R/my_projects/IB_2025/Bioinfo/data/FUSCC-TNBC-protein-mRNA.Rdata")
table(meta$mRNA_Subtype)
# BLIS   IM  LAR  MES #
# 34   20   21   14
table(meta$Lehmann_Subtype)
# BL1 BL2  IM LAR   M MSL UNS 
# 14   5  16  12   7   9  11 
##1.1 ssGSEA score----


exp <- mRNA
tumor_id <- colnames(exp)[grepl("T$", colnames(exp))]
exp <-exp[,tumor_id]#89

sum(is.na(exp))
# 核心：将表达矩阵中所有NA替换为0
exp[is.na(exp)] <- 0


library(GSVA)
##First we should build a parameter object for the desired methodology.R
gsvaPar <- ssgseaParam(exprData = as.matrix(exp), #sample in column 
                       geneSets = genelist,
                       normalize = TRUE)
##Second, we call the gsva() function with the parameter object as first argument. 
gene_score <- gsva(gsvaPar, verbose = FALSE)
gene_score <-as.data.frame(t(gene_score))
gene_score$ID <-rownames(gene_score)

meta$ID <- paste0(meta$sampleID,"T")
gene_score_meta <- merge(gene_score,meta,by="ID")

##1.2 add gene exp----
gene_exp <- exp[interest_gene,]
gene_exp <- as.data.frame(t(gene_exp))
gene_exp$ID <- rownames(gene_exp)
gene_score_meta <- merge(gene_exp,gene_score_meta,by="ID")

##1.3 pathway activity----
library(progeny);library(tidyr)
set.seed(123)
pathways <- progeny(as.matrix(exp), scale=TRUE,
                    organism="Human", #如为小鼠，填 "Mouse"
                    top = 100, perm = 1)
df_activity <- as.data.frame(pathways)
df_activity$ID <- rownames(df_activity)
#add activity
gene_score_meta <- merge(gene_score_meta,df_activity,by="ID")

##1.4 plot----


#boxplot
library(ggplot2);library(ggpubr)
mytheme <- theme(plot.title = element_text(size = 14,color="black",hjust = 0.5),
                 axis.title = element_text(size = 14,color ="black"), 
                 axis.text = element_text(size= 14,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.title.x = element_blank(),
                 #axis.text.x = element_text(angle = 60, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 14),
                 legend.title= element_blank() )#element_text(size= 12) 

for (i in c(interest_gene,"JAK-STAT", "Jak-STAT_pathway","Core_Jak-STAT_pathway") ){
  print(i)
  gene_meta_exp <- gene_score_meta
  #gene_meta_exp$riskscore <- gene_meta_exp$`Jak-STAT_pathway`
  gene_meta_exp$Group <- gene_meta_exp$mRNA_Subtype
  gene_meta_exp$riskscore <- gene_meta_exp[,i]
  ggplot(gene_meta_exp, aes(x = reorder(Group,riskscore), y = riskscore,color=Group ) )+ # STAT4 #  取LOG2不影响显著性
    labs(y=i,x= NULL,title="FUSCC-TNBC")+  
    #geom_boxplot(aes(fill = response),outlier.size = 0.02,size=0.02,outlier.alpha = 0)+ 
    geom_point(alpha=0.7,size=0.5,
               position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                             jitter.height = 0,
                                             dodge.width = 0.7))+
    geom_boxplot(alpha=1,width=0.8,fill=NA,
                 position=position_dodge(width=0.8),
                 size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                 outlier.stroke = 0.5)+
    geom_violin(alpha=0.4,width=0.9,
                position=position_dodge(width=0.8),
                size=0.25)+
    #geom_bar(stat = "identity")+
    EnvStats::stat_n_text()+ #library(EnvStats) # 显示样本数
    scale_fill_manual(values = my36colors )+ #c("skyblue","pink")
    scale_color_manual(values = my36colors)+
    #guides(fill = guide_legend(title = 'STAT4'))+
    theme_bw()+ mytheme + theme(legend.position = "none")#+
  #scale_y_continuous(limits = c(0.5,0.8))
  # stat_compare_means(size=6, #"wilcox.test"
  #                    show.legend= F,label = "p.signif",#p.signif p.format
  #                    label.x =1.5,label.y =max(gene_meta_exp$riskscore)-0.02 )
  ggsave(paste0(folder_path,"/Boxplot_FUSCC-TNBC-",i,".pdf"),width = 3,height = 2.5,family="serif")
  
}

#2.0 TNBC96 (2024_Nat.Comm.) (n=94)----
# article： https://www.nature.com/articles/s41467-024-54145-w#Abs1
# 
# clinical： OS RFS DFS  ;TLS TNM; NO PD-L1
# download：https://doi.org/10.5281/zenodo.8135721、

##load----
count <- readRDS("~/R/my_projects/JAK_STAT_BC/13.0_TNBC/data/bulk_count.RDS")
#clinical <- readRDS("~/R/my_projects/JAK_STAT_BC/13.0_TNBC/data/Clinical.RDS")
clinical <- read.csv("./data/Clinical.csv")
table(clinical$Bareche_molecular_subtype_defined_on_global_pseudobulk)
# BL  IM LAR   M MSL 
# 18  25  11  28  12 
table(clinical$Bareche_molecular_subtype_defined_on_tumor_pseudobulk)
# BL=BLIS  IM LAR   M=MES 
# 41   6  14  33 
table(clinical$Bareche_molecular_subtype_defined_on_stroma_pseudobulk)
# BL  IM LAR   M MSL 
# 5  32   7  14  36

clinical$Group <- clinical$Bareche_molecular_subtype_defined_on_tumor_pseudobulk
#BL basal-like,  IM immunomodulatory, LAR luminal androgen receptor, M mesenchymal, MSL mesenchymal stem-like
table(clinical$Group) #BL=BLIS  IM LAR   M=MES 
library(dplyr)
clinical$Group <- recode(
  clinical$Group,
  "BL" = "BLIS",
  "M" = "MES"
)
table(clinical$Group)
# BLIS   IM  LAR  MES 
# 41    6   14   33 

##2.1 ssGSEA score----


exp <- as.data.frame(log2(count+1)) #count !
#tumor_id <- colnames(exp)[grepl("T$", colnames(exp))]
#exp <-exp[,tumor_id]#89

sum(is.na(exp))
# 核心：将表达矩阵中所有NA替换为0
#exp[is.na(exp)] <- 0


library(GSVA)
##First we should build a parameter object for the desired methodology.R
gsvaPar <- ssgseaParam(exprData = as.matrix(exp), #sample in column 
                       geneSets = genelist,
                       normalize = TRUE)
##Second, we call the gsva() function with the parameter object as first argument. 
gene_score <- gsva(gsvaPar, verbose = FALSE)
gene_score <-as.data.frame(t(gene_score))
gene_score$ID <-paste0(rownames(gene_score),"T")

meta <- clinical
meta$ID <- paste0(meta$ST_TNBC_ID,"T")
gene_score_meta <- merge(gene_score,meta,by="ID")

##1.2 add gene exp----
gene_exp <- exp[interest_gene,]
gene_exp <- as.data.frame(t(gene_exp))
gene_exp$ID <-paste0(rownames(gene_exp),"T")
gene_score_meta <- merge(gene_exp,gene_score_meta,by="ID")

##1.3 pathway activity----
library(progeny);library(tidyr)
set.seed(123)
pathways <- progeny(as.matrix(exp), scale=TRUE,
                    organism="Human", #如为小鼠，填 "Mouse"
                    top = 100, perm = 1)
df_activity <- as.data.frame(pathways)
df_activity$ID <- paste0(rownames(df_activity),"T")
#add activity
gene_score_meta <- merge(gene_score_meta,df_activity,by="ID")

##1.4 plot----

#boxplot
library(ggplot2);library(ggpubr)
mytheme <- theme(plot.title = element_text(size = 14,color="black",hjust = 0.5),
                 axis.title = element_text(size = 14,color ="black"), 
                 axis.text = element_text(size= 14,color = "black"),
                 panel.grid.minor.y = element_blank(),
                 panel.grid.minor.x = element_blank(),
                 axis.title.x = element_blank(),
                 #axis.text.x = element_text(angle = 60, hjust = 1 ),
                 panel.grid=element_blank(),
                 legend.position = "top",
                 legend.text = element_text(size= 14),
                 legend.title= element_blank() )#element_text(size= 12) 

for (i in c(interest_gene,"JAK-STAT", "Jak-STAT_pathway","Core_Jak-STAT_pathway") ){
  print(i)
  gene_meta_exp <- gene_score_meta
  #gene_meta_exp$riskscore <- gene_meta_exp$`Jak-STAT_pathway`
  #gene_meta_exp$Group <- gene_meta_exp$mRNA_Subtype
  gene_meta_exp$riskscore <- gene_meta_exp[,i]
  ggplot(gene_meta_exp, aes(x = reorder(Group,riskscore), y = riskscore,color=Group ) )+ # STAT4 #  取LOG2不影响显著性
    labs(y=i,x= NULL,title="2024_Nat.Comm.")+  
    #geom_boxplot(aes(fill = response),outlier.size = 0.02,size=0.02,outlier.alpha = 0)+ 
    geom_point(alpha=0.7,size=0.5,
               position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                             jitter.height = 0,
                                             dodge.width = 0.7))+
    geom_boxplot(alpha=1,width=0.8,fill=NA,
                 position=position_dodge(width=0.8),
                 size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                 outlier.stroke = 0.5)+
    geom_violin(alpha=0.4,width=0.9,
                position=position_dodge(width=0.8),
                size=0.25)+
    #geom_bar(stat = "identity")+
    EnvStats::stat_n_text()+ #library(EnvStats) # 显示样本数
    scale_fill_manual(values = my36colors )+ #c("skyblue","pink")
    scale_color_manual(values = my36colors)+
    #guides(fill = guide_legend(title = 'STAT4'))+
    theme_bw()+ mytheme + theme(legend.position = "none")#+
  #scale_y_continuous(limits = c(0.5,0.8))
  # stat_compare_means(size=6, #"wilcox.test"
  #                    show.legend= F,label = "p.signif",#p.signif p.format
  #                    label.x =1.5,label.y =max(gene_meta_exp$riskscore)-0.02 )
  ggsave(paste0(folder_path,"/Boxplot_2024_Nat.Comm._",i,".pdf"),width = 3,height = 2.5,family="serif")
  
}

