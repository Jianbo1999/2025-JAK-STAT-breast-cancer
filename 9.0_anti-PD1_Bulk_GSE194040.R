
#~/R/my_projects/JAK_STAT_BC/9.0_anti-PD1_Bulk_GSE194040
if(T){
  rm(list = ls())
  setwd("~/R/my_projects/JAK_STAT_BC/9.0_anti-PD1_Bulk_GSE194040")
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
  library(tidyverse)
  library(forestploter)
  library(GSVA)
  library(SummarizedExperiment)
  
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

#0. load data----
#D:\研究生科研\生物信息学\STAT4\immune_response\GSE194040_Pembrolizumab
load("~/R/my_projects/JAK_STAT_BC/9.0_anti-PD1_Bulk_GSE194040/data/GSE194040_Pembrolizumab_refined.RData")
#1.0 提取合并数据

boxplot(exp_drug[1:100,]) #查看表达4-18
#boxplot(log2(exp_drug[1:100,]) )#2-4

gene_exp <- exp_drug[rownames(exp_drug) %in% genes,] #150-69
gene_exp<- as.data.frame(t(gene_exp))
gene_exp$id <-rownames(gene_exp)#69

#处理临床数据
gene_meta_exp <- meta_drug
# R Responser       NR  Non-responser
gene_meta_exp$response <- ifelse(gene_meta_exp$pcr=="1","R","NR") #31 -38

gene_meta_exp <-merge(gene_meta_exp,gene_exp,by="id")#merge 69-14
table(gene_meta_exp$arm)

#1.0 score ----
library(GSVA)
##First we should build a parameter object for the desired methodology.R
gsvaPar <- ssgseaParam(exprData = as.matrix(exp_drug), #sample in column
                       geneSets = genelist,
                       normalize = TRUE)
##Second, we call the gsva() function with the parameter object as first argument. 
gene_score <- gsva(gsvaPar, verbose = FALSE)
gene_score <-as.data.frame(t(gene_score))
gene_score$id <-rownames(gene_score)
#merge
gene_meta_exp <-merge(gene_meta_exp,gene_score,by="id")
gene_meta_exp$riskscore<-gene_meta_exp$`Jak-STAT_pathway`

gene_meta_exp$riskscore_group <-ifelse(gene_meta_exp$riskscore > median(gene_meta_exp$riskscore),"High","Low")
table(gene_meta_exp$riskscore_group )
# High  Low 
# 34   35

#2.0 绘图----
#2.1 boxplot
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

library(ggplot2);library(ggpubr)
ggplot(gene_meta_exp, aes(x = response, y = riskscore,color=response ) )+ # STAT4 #  取LOG2不影响显著性
  labs(y="Score",x= NULL)+  
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
  scale_fill_manual(values = c("skyblue","pink") )+ #"#1CB4B8", "#EB7369" "skyblue","pink"
  scale_color_manual(values = c("skyblue","pink") )+
  #guides(fill = guide_legend(title = 'STAT4'))+
  theme_bw()+ mytheme + theme(legend.position = "none")+
  stat_compare_means(size=6, #"wilcox.test"
                     show.legend= F,label = "p.signif",#p.signif p.format
                     label.x =1.5,label.y =max(gene_meta_exp$riskscore)-0.02 )
ggsave(paste0(folder_path,"/riskscore_response_1.pdf"),width = 2.5,height = 2.5,family="serif")

#3.0 热图----
library(dplyr)
exp_marker1 <- arrange(gene_meta_exp, response)#排序 69-26
rownames(exp_marker1)<-exp_marker1$id

library(pheatmap)
# pheatmap::pheatmap(exp_marker1[c(genes_1,"riskscore")],cluster_rows = F, 
#                    cellwidth=8,cellheight = 8,
#                    annotation_row=exp_marker1[c(13,24,26)]
# )

p1 <-pheatmap(t(exp_marker1[c(genes_1)]),
              scale = 'row',border = T,border_color =NA,
              cluster_cols  = F, show_colnames = F,
              cellwidth=10,cellheight = 10,gaps_col= 38, #table(exp_marker1$response)
              color=colorRampPalette( c("navy", "white", "firebrick3"))(10),
              annotation_col=exp_marker1[c(13)]
)
p1
#保存热图
save_pheatmap_pdf <- function(x, filename, width, height) {
  library(grid)
  x$gtable$grobs[[1]]$gp <- gpar(lwd = 0.05)#修改聚类树线条宽度
  stopifnot(!missing(x))
  stopifnot(!missing(filename))
  pdf(filename, width=width, height=height,family = "serif")
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}#保存热图的函数
save_pheatmap_pdf(p1,paste0(folder_path,"/pheatmap_3.pdf"),12,3)

#4.0 compare AUC----
#读取文章biomarker
biomarker <-openxlsx::read.xlsx("./data/Table S2 Patient-level biomarker scores.xlsx",rows = c( 2:1000))
colnames(biomarker)[1] <-"id"
#987-54
#匹配出exp_marker3中有的样本 #69-34
exp_marker4 <-merge(exp_marker1,subset(biomarker,biomarker$id %in% exp_marker1$id),by="id")#69-220
rownames(exp_marker4)<-exp_marker4$id
#table s1 Immune pathway mrna 7
#191:198

#绘制带有注释的ssgsea得分表达热图！
exp_marker4 <- arrange(exp_marker4, response)#排序


varis <-c("riskscore",genes_1,colnames(exp_marker4)[191:198] )
#varis <-c("riskscore",genes_1 )
ann_colors = list(
  response=c(NR="skyblue",R="pink")
)
p4 <-pheatmap(t( scale(exp_marker4[varis]) ),  #t(scale(exp_marker4[c(2:7)])) #STATs p3
              cluster_cols  = F, show_colnames = F,cluster_rows  = F,
              cellwidth=10,cellheight = 10,
              gaps_col= 38,gaps_row= c(1,11), border_color = NA,#table(exp_marker1$response)
              color=colorRampPalette( c("navy", "white", "firebrick3"))(10),#alpha(rev(RColorBrewer::brewer.pal(11,"PiYG") ), 0.7),   #"BrBG" 褐浅蓝  "PiYG"紫绿
              #alpha(colorRampPalette(colors = c("orange", "white", "purple"))(10),0.6) ,#colorRampPalette( c("navy", "white", "firebrick3"))(10) ,#
              annotation_col=exp_marker4[c(13,177,218)],annotation_colors = ann_colors,
)
p4
#save_pheatmap_pdf(p4, paste0(folder_path,"/pheatmap_4.pdf"),14,4)
save_pheatmap_pdf(p4, paste0(folder_path,"/pheatmap_4-1.pdf"),14,3)
write.csv(exp_marker4,paste0(folder_path,"/exp_marker4.csv") )


##B. ROC 森林图
#B.1计算响应预测ROC
library(pROC)
# roc_data <- roc(exp_marker2$response,exp_marker2$ssGSEA,ci=TRUE) 
# roc_data[["ci"]];roc_data[["auc"]] 
# round(roc_data[["auc"]],1);round(roc_data[["ci"]],2)[c(1,3)] #round(roc_data[["ci"]],2)#L95 AUC H95

df_1 <-data.frame()
for (i in varis ){
  roc_res<- roc(exp_marker4$response,exp_marker4[,i],ci=TRUE) 
  print(i );print(round(roc_res[["ci"]],2) )
  df <-data.frame(name=i,
                  L95=round(roc_res[["ci"]],2)[1],
                  H95=round(roc_res[["ci"]],2)[3],
                  AUC=round(roc_res[["ci"]],2)[2])
  df_1<-rbind(df_1,df)
}
df_1 #oc in Score


#B.2 绘制ROC森林图
library(ggplot2)
#df_1 <- arrange(df_1, AUC)#排序

df_1$group_col <- c( "skyblue",rep("#e7a40e", 10),#e7a40e"
                     rep("#1c6891", 8)#, rep("#a59d70", 5), rep("#4f4a30", 3)
)
df_1$group_col <-factor(df_1$group_col,levels = unique(df_1$group_col),ordered = T)
df_1$name <-factor(df_1$name,levels = unique(df_1$name)) #,ordered = T

ggplot(df_1)+
  # 0轴竖线：
  geom_hline(yintercept = 0, linewidth = 0.3)+
  # 线条：
  #geom_linerange(aes(name, ymin = L95, ymax = H95, color = name), show.legend = F)+
  geom_linerange(aes(name, ymin = 0.4, ymax = AUC), show.legend = F,color = df_1$group_col)+ #, color = group_col
  #geom_label()+
  
  # 散点：
  geom_point(aes(name, AUC, size=AUC),color = df_1$group_col)+ #color = group_col,
  # annotate("rect",
  #          xmin = c(0.5,7.5),  #7.5 8.5 top ssGSEA
  #          xmax = c(7.5,8.5),
  #          ymin = 0.4, ymax = 1, alpha = 0.2, fill = rev(unique(df_1$group_col))) + #背景填充
  annotate("text", label = df_1$AUC, x = df_1$name, y = df_1$AUC+0.1,size=3.0 )+ #,  colour =df_1$group_col
  # annotate("rect",
  #          xmin = c(0.5,6.5,7.5),  #6.5 7.5 top ssGSEA
  #          xmax = c(6.5,7.5,14.5),
  #          ymin = 0.4, ymax = 1, alpha = 0.1, fill =unique(df_1$group_col)) + #rev(unique(df_1$group_col))
  scale_y_continuous(expand = c(0,0),limits = c(0.4,1))+
  xlab("")+
  ylab("AUC")+
  theme_bw(base_size = 12)+coord_flip()+
  theme(legend.position = "none",axis.text.y = element_text(color =df_1$group_col) ) #,axis.text.y = element_text(color =df_1$name )
ggsave(paste0(folder_path,"/plots_forest_auc_all_2.pdf"), height = 3.5, width = 4,family="serif")
write.csv(df_1,paste0(folder_path,"/ROC_1.csv") )


#ROC----
library(pROC)
roc_data <- roc(exp_marker1$response,exp_marker1$riskscore,ci=TRUE) 
roc_data[["ci"]]
roc_data[["auc"]] #AUC = 0.68
#reportROC::reportROC(gold = exp_marker2$response,predictor = exp_marker2$ssGSEA,important = "se",plot=T)

library(ggplot2)
g <- ggroc(roc_data,legacy.axes = TRUE, alpha = 1, colour = "pink",  size = 1)
g
g + #ggtitle("ROC curve") +
  #geom_ribbon(aes(x=1 - roc_data$specificities, ymin=0, ymax=roc_data$sensitivities),alpha=0.5,fill="skyblue") +#填充
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), color="grey50", linetype="dashed")+
  xlab("1-Specificity (FPR)") + ylab("Sensitivity (TPR)")+
  theme_bw()+
  annotate("text",x=0.6,y=0.2,label=paste(" AUC: ",round(roc_data[["auc"]],1) ,"\n ","95%CI :",round(roc_data[["ci"]],2)[1],"-",round(roc_data[["ci"]],2)[3] ),size=6,colour = "red",family="serif")+
  scale_x_continuous(expand = c(0,0.01)) +scale_y_continuous(expand = c(0,0.01))+ #X轴数据两边c(0,0.02),否则横坐标显示不全
  theme(text=element_text(size=16,family="serif")) +#字体TNR
  theme(plot.margin=unit(rep(1,5),'lines') )
ggsave(paste0(folder_path,"/ROC_riskscore.pdf"),family="serif",width = 4,height=4)
