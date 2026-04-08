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
  genes_2 <- c("CD274","PDCD1LG2","PDCD1")
  interest_gene <-"STAT4"
  
  my_cancer<-c("BRCA")
  
  genelist<- list("Jak-STAT_pathway"=genes,"Core_Jak-STAT_pathway"=genes_1
                  #"Proliferation"= c("BIRC5", "CCNB1", "CDC20", "NUF2", "CEP55", "NDC80", "MKI67", "PTTG1", "RRM2", "TYMS","UBE2C")
  ) #https://www.nature.com/articles/s41588-021-00911-1
  input_gene= c("Jak-STAT_pathway","Core_Jak-STAT_pathway") #genes#c(genes,"H3K36me","5mC" )
  project = "GSE194040 anti-PD-1"
}

#load and score----
#load("~/R/my_projects/JAK_STAT_BC/9.0_anti-PD1_Bulk_GSE194040/data/GSE194040_Pembrolizumab_refined.RData")


boxplot(exp_drug[1:100,]) #查看表达4-18
#boxplot(log2(exp_drug[1:100,]) )#2-4

gene_exp <- exp_drug[rownames(exp_drug) %in% c(genes,genes_2),] #153-69
gene_exp<- as.data.frame(t(gene_exp))
gene_exp$id <-rownames(gene_exp)#69

#处理临床数据
gene_meta_exp <- meta_drug
# R Responser       NR  Non-responser
gene_meta_exp$response <- ifelse(gene_meta_exp$pcr=="1","R","NR") #31 -38

gene_meta_exp <-merge(gene_meta_exp,gene_exp,by="id")#merge 69-14
table(gene_meta_exp$arm)

#score 
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
#multi-cox----

library(survival)
library(survminer)
library(dplyr)
library(broom)
library(ROCR)
library(pROC)

##JSP + PD-L1/2 + PD-1 ----

# ===================== 1. 数据整理：response → N=0, R=1
cox_data <- gene_meta_exp %>%
  select(
    response, PDCD1, CD274, PDCD1LG2, riskscore
  ) %>%
  mutate(
    status = ifelse(response == "NR", 0,
                    ifelse(response == "R", 1, NA))
  ) %>%
  na.omit()

# ===================== 2. 多因素 COX 模型
cox_model <- coxph(
  Surv(status) ~ PDCD1 + CD274 + PDCD1LG2 + riskscore,
  data = cox_data
)

# ===================== 3. 输出模型结果（HR、P值）
summary(cox_model)
tidy(cox_model, exponentiate = T, conf.int = T)

# ===================== 4. 计算模型得分（risk score）
cox_data$model_score <- predict(cox_model, newdata = cox_data, type = "risk")

# roc_obj <- roc(
#   status ~ model_score,
#   data = cox_data,
#   quiet = TRUE
# )
# 
# 
# print(roc_obj$auc)
# 
# # 绘制 ROC 图
# ggroc(roc_obj, legacy.axes = TRUE, size = 1.2, color = "red") +
#   ggtitle(paste0("AUC = ", round(roc_obj$auc, 3))) +
#   geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray") +
#   theme_minimal()

roc_list <- list()
roc_list[["PD-1"]] <- roc(cox_data$status, cox_data$PDCD1, quiet=T)
roc_list[["PD-L1"]] <- roc(cox_data$status, cox_data$CD274, quiet=T)
roc_list[["PD-L2"]] <- roc(cox_data$status, cox_data$PDCD1LG2, quiet=T)
roc_list[["JSP"]] <- roc(cox_data$status, cox_data$riskscore, quiet=T)
roc_list[["model_score"]] <- roc(cox_data$status, cox_data$model_score, quiet=T)

auc_values <- sapply(roc_list, function(x) round(x$auc, 2))
names(roc_list) <- paste0(names(roc_list), " (AUC=", auc_values, ")")
ggroc(roc_list, legacy.axes=TRUE, size=0.5) +
  # scale_color_manual(values=c("PDCD1"="steelblue",
  #                             "CD274"="orange",
  #                             "PDCD1LG2"="purple",
  #                             "JSP"="green3",
  #                             "model_score"="red2")) +
  #scale_color_manual(values = my36colors)+
  scale_color_manual(values =c("#F37C95", "#D6E7A3", "#20854E","#0072B5","#E95C59") )+
  geom_abline(intercept=0,slope=1,linetype="dashed",color="gray50",linewidth=0.2) +
  labs( title=project, color = "") +
  theme_bw(base_size = 12)+#theme_minimal() +
  theme(legend.position="right",panel.grid.major = element_blank() )
ggsave(paste0(folder_path,"/multi-COX-roc-JSP+PDL1-2-PD1.pdf"),width = 4.8,height = 2.7,family="serif")


##JSP + PD-L1 ----
cox_model <- coxph(
  Surv(status) ~ CD274 + PDCD1LG2 + riskscore,
  data = cox_data
)
summary(cox_model)
tidy(cox_model, exponentiate = T, conf.int = T)

cox_data$model_score <- predict(cox_model, newdata = cox_data, type = "risk")

roc_list <- list()
roc_list[["PD-L1"]] <- roc(cox_data$status, cox_data$CD274, quiet=T)
roc_list[["PD-L2"]] <- roc(cox_data$status, cox_data$PDCD1LG2, quiet=T)
roc_list[["JSP"]] <- roc(cox_data$status, cox_data$riskscore, quiet=T)
roc_list[["JSP + PD-L1/2"]] <- roc(cox_data$status, cox_data$model_score, quiet=T)

auc_values <- sapply(roc_list, function(x) round(x$auc, 2))
names(roc_list) <- paste0(names(roc_list), " (AUC=", auc_values, ")")
ggroc(roc_list, legacy.axes=TRUE, size=0.5) +
  scale_color_manual(values =c( "#D6E7A3", "#20854E","#0072B5","#F37C95","#E95C59") )+
  geom_abline(intercept=0,slope=1,linetype="dashed",color="gray50",linewidth=0.2) +
  labs( title=project, color = "") +
  theme_bw(base_size = 12)+#theme_minimal() +
  theme(legend.position="right",panel.grid.major = element_blank() )
#ggsave(paste0(folder_path,"/multi-COX-roc-JSP+PDL1.pdf"),width = 5,height = 2.7,family="serif")
ggsave(paste0(folder_path,"/multi-COX-roc-JSP+PDL1-2.pdf"),width = 4.8,height = 2.7,family="serif")


#Predictive Model Construction
#A multivariable Cox proportional hazards regression model was constructed to integrate the expression values of PDCD1, CD274, PDCD1LG2, and the JAK-STAT pathway (JSP) risk score, aiming to accurately predict clinical treatment response (R: complete response, NR: non-response). A consolidated model risk score was subsequently computed for each sample using the regression coefficients estimated by the coxph() function within the survival package (version 3.7-0).
#To assess and compare the predictive power, receiver operating characteristic (ROC) curve analysis was performed. The area under the curve (AUC) was utilized as the primary metric to evaluate the discriminative ability of single markers and the multi-gene signature. Visualization and quantitative comparison of ROC curves were implemented with the pROC R package (version 1.18.5).