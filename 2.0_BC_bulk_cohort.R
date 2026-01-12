
#~/R/my_projects/JAK_STAT_BC/2.0_BC_bulk_cohort

if(T){
  rm(list = ls())
  setwd("~/R/my_projects/JAK_STAT_BC/2.0_BC_bulk_cohort")
  folder_path <- "./out"
  # 检查文件夹是否存在
  if (!dir.exists(folder_path)) {
    # 如果文件夹不存在，则创建它
    dir.create(folder_path)
    print(paste("Folder created at", folder_path))
  } else {
    print(paste("Folder already exists at", folder_path))
  }
  library(ggplot2);library(ggpubr)
  colors <-c("#0072B5","#E18727","#F37C95","#20854E", '#D6E7A3', '#57C3F3', '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#476D87', '#E95C59', '#E59CC4', '#AB3282', '#23452F', '#BD956A', '#8C549C', '#585658', '#9FA3A8', '#E0D4CA', '#5F3D69', '#C5DEBA', '#58A4C3', '#E4C755', '#F7F398', '#AA9A59', '#E63863', '#E39A35', '#C1E6F3', '#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD', '#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963', '#968175') #颜色设置 
  library(survival)
  library(survminer)
  library(timeROC)
  library(ggsci)#配色
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

load("~/R/my_projects/JAK_STAT_BC/2.0_BC_bulk_cohort/2.0_BC_bulk_cohort.RData")


#清洗好的数据集-only tumor samples
#load("~/R/my_projects/JAK_STAT_BC/2.0_BC_bulk_cohort/BRCA_bulk_TCGA_METABRIC_SCAN_B.RData")
#FUSCC-TNBC ?

#0 prepare TCGA_METABRIC_SCAN_B for exp in subtype----
#G:\bioinfo\breast_cancer_dataset\integrated
meta_all<-readRDS("~/R/my_projects/BC_TF_2025/data/BRCA_bulk_TCGA_METABRIC_SCAN_B_meta.rds")
cohort <-readRDS("~/R/my_projects/BC_TF_2025/data/BRCA_bulk_TCGA_METABRIC_SCAN_B.rds")
#log2+1

##prepare meta data----
##remove OS OS.time
for (y in names(cohort) ) {
  #cohort[[y]] <-cohort[[y]][c("ID","OS.time","OS",colnames(cohort[[y]])[4:ncol(cohort[[y]])])]#<-cohort[y] %>% select(c("ID","OS.time","OS"),everything() ) 
  cohort[[y]] <-cohort[[y]][,-c(1:2) ]
  # cohort[[y]] <-cohort[[y]][c("ID","OS.time","OS",gene) ]#<-cohort[y] %>% select(c("ID","OS.time","OS"),everything() ) 
  # cohort[[y]] <- cohort[[y]][cohort[[y]]$OS.time != 0,]
  # cohort[[y]] <- na.omit(cohort[[y]] )#移除缺失值
} #
#TCGA取tumor #table(type) #N 71 T 531
#tumor_id <-rownames(cohort[["TCGA"]])[(as.numeric(substr(rownames(cohort[["TCGA"]]),14,15)) <10) ==T]
#table(as.numeric(substr(rownames(cohort[["TCGA"]]),14,15)) <10 ) #1-530 5-1 11-71
#cohort[["TCGA"]]<-cohort[["TCGA"]][rownames(cohort[["TCGA"]]) %in% tumor_id,]




##clean meta data----
###TCGA----
meta <-meta_all$TCGA_meta #1097 all tumor

#not run
if(F){
  meta <-tibble::column_to_rownames(meta ,var = "ID" )
  #####pathologic_stage
  meta$pathologic_stage<-meta$ajcc_pathologic_tumor_stage
  meta$pathologic_stage <-gsub("Stage ","",meta$pathologic_stage)
  meta$pathologic_stage<-ifelse(meta$pathologic_stage %in% c("[Discrepancy]","[Not Available]","X"),"Unknown", meta$pathologic_stage)
  meta$pathologic_stage<-ifelse(meta$pathologic_stage %in% c("I","IA","IB"),"I",
                                ifelse(meta$pathologic_stage %in% c("II","IIA","IIB"),"II",
                                       ifelse(meta$pathologic_stage %in% c("III","IIIA","IIIB","IIIC"),"III",meta$pathologic_stage)
                                )
  )
  table(meta$pathologic_stage)
  #T_stage 
  meta$T_stage <- meta$ajcc_tumor_pathologic_pt
  meta$T_stage <-ifelse(meta$T_stage == "TX","Unknown",
                        meta$T_stage
  )
  meta$T_stage <-ifelse(meta$T_stage =="Unknown","Unknown" ,
                        substr(meta$T_stage,2,2)
  )
  table(meta$T_stage) #1       2       3       4 Unknown
  #N_stage 
  meta$N_stage <- meta$ajcc_nodes_pathologic_pn
  meta$N_stage <- substr(meta$N_stage,0,2)
  meta$N_stage <-ifelse(meta$N_stage == "NX","Unknown",
                        meta$N_stage
  )
  table(meta$N_stage) #N0      N1      N2      N3 Unknown 
  #M_stage
  meta$M_stage <- meta$ajcc_metastasis_pathologic_pm
  
  meta$M_stage <-ifelse(meta$M_stage %in% c("MX","cM0 (i+)"),"Unknown",
                        meta$M_stage
  )
  meta$M_stage <- ifelse(meta$M_stage=="Unknown",meta$M_stage,
                         substr(meta$M_stage,0,2)
  )
  table(meta$M_stage)
  ####age
  meta$Age <-as.numeric(meta$age_at_diagnosis)
  ####HER2+ PR TNBC 
  table(meta$ER);table(meta$PR);table(meta$HER2)
  meta$IHC <-ifelse(meta$HER2=="Positive","HER2+", #HER2 +
                    ifelse(meta$HER2=="Negative" & meta$ER=="Positive","Luminal","Unkonwn")   #HER2- ER+
  )
  meta$IHC <-ifelse(meta$HER2=="Negative" & meta$ER=="Negative" & meta$ER=="Negative","TNBC",meta$IHC)   #HER2- ER+
  #meta$IHC <- ifelse(meta$IHC =="Luminal" & meta$PR=="Positive","Luminal A",meta$IHC)
  #meta$IHC <- ifelse(meta$IHC =="Luminal" & meta$PR=="Negative","Luminal B",meta$IHC)
  #all luminal A
  table(meta$IHC)  #Lunimal=HER2- ER+; PR+ --> A,PR- --> B
  # HER2+ Luminal    TNBC Unkonwn 
  # 164     437     126     370 
  #meta_all$TCGA_meta <-meta
}#生存数据丢失


####METABRIC----
meta<-meta_all$METABRIC_meta
meta$T_stage <-meta$Tumor_Stage #NA
meta$T_stage <- ifelse(is.na(meta$T_stage) ==T,"Unknown",meta$T_stage) #720 UNKNOWN
####HER2+ PR TNBC 
table(meta$ER_Status);table(meta$PR_Status);table(meta$HER2_Status)
meta$IHC <-ifelse(meta$HER2_Status=="Positive","HER2+", #HER2 +
                  ifelse(meta$HER2_Status=="Negative" & meta$ER_Status=="Positive","Luminal","Unkonwn")   #HER2- ER+
)
meta$IHC <-ifelse(meta$HER2_Status=="Negative" & meta$ER_Status=="Negative" & meta$ER_Status=="Negative","TNBC",meta$IHC)   #HER2- ER+
meta$IHC <- ifelse(meta$IHC =="Luminal" & meta$PR_Status=="Positive","Luminal A",meta$IHC)
meta$IHC <- ifelse(meta$IHC =="Luminal" & meta$PR_Status=="Negative","Luminal B",meta$IHC)
#luminal A=973 bB=424 #Lunimal=HER2- ER+; PR+ --> A,PR- --> B
meta$IHC <- ifelse(is.na(meta$IHC) ==T,"Unknown",meta$IHC) 
table(meta$IHC)  
# HER2+ Luminal A Luminal B      TNBC   Unknown 
# 247       973       424       333       529 
#meta_all$METABRIC_meta <-meta

####SCAN-B----
meta<-meta_all$SCAN_B_meta  #tumor_size
meta$Age <-as.numeric(meta$age_at_diagnosis)
#meta$er_status== " 0" " 1"
meta$IHC <-ifelse(meta$her2_status==" 1","HER2+", #HER2 +
                  ifelse(meta$her2_status==" 0" & meta$er_status==" 1","Luminal","Unkonwn")   #HER2- ER+
)
meta$IHC <-ifelse(meta$her2_status==" 0" & meta$er_status==" 0" & meta$pgr_status==" 0","TNBC",meta$IHC)   #HER2- ER+
meta$IHC <- ifelse(meta$IHC =="Luminal" & meta$pgr_status==" 1","Luminal A",meta$IHC)
meta$IHC <- ifelse(meta$IHC =="Luminal" & meta$pgr_status==" 0","Luminal B",meta$IHC)
#luminal A=973 bB=424 #Lunimal=HER2- ER+; PR+ --> A,PR- --> B
meta$IHC <- ifelse(is.na(meta$IHC) ==T,"Unknown",meta$IHC) 
table(meta$IHC)  
# HER2+   Luminal Luminal A Luminal B      TNBC   Unkonwn 
# 392       100      2054       123       133       267
meta$N_stage <-ifelse(meta$lymph_node_group==" NodeNegative",0,
                      ifelse(meta$lymph_node_group %in% c(" SubMicroMet"," NA"),"Unknown",meta$lymph_node_group)
                      )
table(meta$N_stage)
# 1to3    4toX       0 Unknown 
# 813     282    1811     163 
#meta_all$SCAN_B_meta <-meta

#save.image("~/R/my_projects/JAK_STAT_BC/2.0_BC_bulk_cohort/BRCA_bulk_TCGA_METABRIC_SCAN_B.RData")
#保存修改好的数据集

#c("pathologic_stage","T_stage","N_stage","M_stage")
#table(is.na(meta_all[["TCGA_meta"]][["OS.time"]])) OS.time 有NA


#1.0 exp in BRCA subtypes----

##1.1 TCGA_METABRIC_SCAN_B ----
#load("~/R/my_projects/JAK_STAT_BC/2.0_BC_bulk_cohort/BRCA_bulk_TCGA_METABRIC_SCAN_B.RData")
#only tumor 
intersect(colnames(meta_all$TCGA_meta),
          intersect(colnames(meta_all$METABRIC_meta),colnames(meta_all$SCAN_B_meta))
          ) #"OS.time" "OS"      "Type" TNBC/Non-TNBC    "Age"     "IHC"
#exp and KM ROC
for (y in names(cohort) ) {
  print(y)
  #cohort[[y]] <-cohort[[y]][c("ID","OS.time","OS",colnames(cohort[[y]])[4:ncol(cohort[[y]])])]#<-cohort[y] %>% select(c("ID","OS.time","OS"),everything() ) 
  #cohort[[y]] <-cohort[[y]][,-c(1:2) ]
  # cohort[[y]] <-cohort[[y]][c("ID","OS.time","OS",gene) ]#<-cohort[y] %>% select(c("ID","OS.time","OS"),everything() ) 
  # cohort[[y]] <- cohort[[y]][cohort[[y]]$OS.time != 0,]
  # cohort[[y]] <- na.omit(cohort[[y]] )#移除缺失值
} 
rownames(meta_all$METABRIC_meta) <-meta_all$METABRIC_meta$Sample_ID
rownames(meta_all$SCAN_B_meta) <-meta_all$SCAN_B_meta$ID
#genes_exp_clinical$sample <-genes_exp_clinical$ID
#save.image("~/R/my_projects/JAK_STAT_BC/2.0_BC_bulk_cohort/BRCA_bulk_TCGA_METABRIC_SCAN_B.RData")


#c("pathologic_stage","T_stage","N_stage","M_stage","IHC","Type")
#TNBC KM !

#km function
plot_km <-function(fit,data,title=NULL){
  print(
    ggsurvplot_list(fit,data = data,pval = T,pval.method = TRUE,##是否添加P值
                    conf.int = F,### 是否添加置信区间
                    legend.title = title,#i, # 设置图例标题title = cancer,
                    legend.labs = c("High", "Low"), # 指定图例分组标签
                    risk.table = F, # 是否添加风险表
                    risk.table.col = "strata", 
                    censor=F,size=0.5,
                    surv.scale="percent",
                    ###linetype = "strata",
                    #surv.median.line = "hv", # 是否添加中位生存线
                    risk.table.y.text.col = F,risk.table.y.text = FALSE,
                    ggtheme = theme_bw(base_family = "serif")+theme(legend.text = element_text(colour = c("red", "blue")))
                    +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
                    +theme(plot.title = element_text(hjust = 0.5,size = 16),axis.title.y.left = element_text(size = 16,vjust = 1),axis.title.x.bottom = element_text(size = 16,vjust = 0))
                    +theme(axis.text.x.bottom = element_text(size = 12,vjust = -0.8,colour = "black")) #,face = "bold"
                    +theme(axis.text.y.left = element_text(size = 12,vjust = 0.5,hjust = 0.5,angle = 90,colour = "black"))
                    +theme(legend.title = element_text(family = "Times",colour = "black",size = 12))
                    +theme(legend.text = element_text(family = "Times",colour = "black",size =12)), # Change ggplot2 theme
                    palette = c("#bc1e5d", "#0176bd"),#"lacent",##如果数据分为两组，就写两种颜色，三组就写三种颜色，具体的颜色搭配参考上一期发布的R语言颜色对比表
                    #xlim=c(10,5000),##x轴的阈值，根据随访数据的天数进行更改
                    xlab = "Years")##随访的时间时天，就写Days，月份就写Months
  ) 
  #print(p1)
}
uni_cox_fun <-function(data2,input_gene){
  outTab=data.frame()
  sigGenes=c("OS.time","OS")#c("futime","fustat")
  for(i in input_gene ){ #trainset[,3:ncol(trainset)]
    cox <- coxph(Surv(OS.time, OS) ~ data2[,i], data = data2)
    coxSummary = summary(cox)
    coxP=coxSummary$coefficients[,"Pr(>|z|)"]
    if(T){ #coxP<0.2
      sigGenes=c(sigGenes,i)
      outTab=rbind(outTab,
                   cbind(id=i,
                         HR= coxSummary$conf.int[,"exp(coef)"] ,
                         HR.95L=coxSummary$conf.int[,"lower .95"],
                         HR.95H=coxSummary$conf.int[,"upper .95"],
                         pvalue=coxSummary$coefficients[,"Pr(>|z|)"])
      )
    }
  }
  colnames(outTab)#;str(outTab)
  outTab[2:5]<-apply(outTab[2:5],2,as.numeric)
  outTab[2:4]<-round(outTab[2:4],2)
  outTab$Sign. <-ifelse(outTab$pvalue<0.001,"***",ifelse(outTab$pvalue<0.01,"**",ifelse(outTab$pvalue<0.05,"*","ns")))
  outTab$pvalue<-round(outTab$pvalue,3)
  # outTab<-dplyr:: arrange(outTab,desc(HR) )#HR排序 ?
  # outTab$id <-as.factor(outTab$id)
  outTab$color <-ifelse(outTab$HR < 1,"skyblue","pink")
  outTab$HR.95H <-ifelse(outTab$HR.95H >5,5,outTab$HR.95H)
  ##plot
  p <-ggplot(outTab, aes(HR, reorder(id,-HR) ))+
    geom_errorbarh(aes(xmax =HR.95H, xmin = HR.95L,color="black"),size= 0.5,height = 0.1) + #color=id
    #geom_point(aes(color=id),size=2.5)+
    geom_point(aes(size=HR,color=color)  )+ #,color="skyblue" id
    #scale_x_continuous(limits= c(-2.5, 2.5))+
    geom_vline(aes(xintercept = 1),color="gray",linetype="dashed", size = 0.5) +
    geom_text(size=4, aes(x = HR+0.1, label = Sign.),color="navy")+#,hjust = 1
    xlab('HR')+ ylab(' ')+
    theme_bw()+
    scale_color_manual(values = c("black","pink","skyblue") )+ #colors
    theme(axis.text.x = element_text(size = 14, color = "black"),axis.text.y = element_text(size = 14, color = "black"),
          panel.grid.minor = element_blank(),panel.grid.major  = element_blank(),
    )+
    theme(title=element_text(size=14),legend.position = "none")
  print(p)
}

if(F){ #load 
  cohort <-readRDS("~/R/my_projects/JAK_STAT_BC/2.0_BC_bulk_cohort/bulk_data/BRCA_bulk_TCGA_METABRIC_SCAN_B.rds")
  meta_all<-readRDS("~/R/my_projects/JAK_STAT_BC/2.0_BC_bulk_cohort/bulk_data/BRCA_bulk_TCGA_METABRIC_SCAN_B_meta.rds")
  y="TCGA"
  #exp<-cohort[[y]][,-c(1:2)]
  #exp<- as.data.frame(t(exp))
  meta_all[[paste0(y,"_meta")]]$clinical_stage  <- meta_all[[paste0(y,"_meta")]]$pathologic_stage
  # meta <-meta_all[[paste0(y,"_meta")]]
  # meta$ID <-rownames(meta)
  # clinical_df<-meta[,c(112:127)]
}
rownames(meta_all[[1]])
rownames(meta_all[[2]]) <- meta_all[[2]]$ID
rownames(meta_all[[3]]) <- meta_all[[3]]$ID

#y=3
for (y in 1:length(cohort) ){ #1:length(cohort)
  print(y);print(names(cohort[y]))
  
  #run!
  genes_exp_clinical <-meta_all[[y]]
  genes_exp_clinical$sample <-rownames(genes_exp_clinical)
  exp <-as.data.frame(t(cohort[[y]][,-c(1:2)] ))  #as.data.frame(t(cohort[[y]]))#panc_22#genes_exp_clinical[,1:152]
  
  ###score
  #exp <- tibble::column_to_rownames(exp, var = "sample")
  library(GSVA)
  #gsva_matrix<- gsva(as.matrix(exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
  ##First we should build a parameter object for the desired methodology.R
  gsvaPar <- ssgseaParam(exprData = as.matrix(exp), 
                         geneSets = genelist,
                         normalize = TRUE)
  ##Second, we call the gsva() function with the parameter object as first argument. 
  gsva_matrix <- gsva(gsvaPar, verbose = FALSE)
  gsva_matrix<-as.data.frame(gsva_matrix);gsva_matrix<-as.data.frame(t(gsva_matrix ))
  gsva_matrix$sample <-rownames(gsva_matrix)
  ###merge in meta
  genes_exp_clinical_1 <-merge(genes_exp_clinical,gsva_matrix,by="sample")
  rownames(genes_exp_clinical_1)<-genes_exp_clinical_1$sample
  #merge gene exp
  genes_exp <-as.data.frame(t(exp[rownames(exp) %in% genes,]))
  genes_exp$sample <-rownames(genes_exp)
  genes_exp_clinical_1 <-merge(genes_exp_clinical_1,genes_exp,by="sample")
  genes_exp_clinical_1$OS <- as.numeric( gsub(" ","",genes_exp_clinical_1$OS ) )
  ###0. KM and exp in all cohort
  #exp only TCGA
  df <- genes_exp_clinical_1
  for (i in c("STAT4")){ #"Jak-STAT_pathway" KM-分组指标
    print(i)
  }
  
  df$Gene <-df[,i]#df$`Jak-STAT_pathway`#df[,i]
  df1 <- df
  #km
  for(km_method in c("median","bestcutoff") ){
    ##OS
    ifelse(km_method == "median",
           df1$Risk <- ifelse(df1$Gene>median(df1$Gene) ,"High","Low"),
           df1$Risk <- ifelse(df1$Gene> surv_cutpoint(df1, time = "OS.time", event = "OS",
                                                      variables = c("Gene"))$cutpoint[,1] ,
                              "High","Low")
    )
    
    fit <- survfit(Surv(OS.time, OS) ~ Risk, data =df1 )
    print(plot_km(fit,df1,title = "Score" ) ) 
    ggsave(filename = paste0("./out/",km_method,"_",names(cohort[y]),"_",i,"_OS.pdf"), height = 3.5, width = 3.5)
    
  }
  
  if(names(cohort[y])=="TCGA"){
    df1$PFI.time <- df1$PFI.time/365
    df1 <- df1[,c("Gene","PFI","PFI.time")]
    df1 <-na.omit(df1)
    for(km_method in c("median","bestcutoff") ){
      #PFI/PFS
      ifelse(km_method == "median",
             df1$Risk <- ifelse(df1$Gene>median(df1$Gene) ,"High","Low"),
             df1$Risk <- ifelse(df1$Gene> surv_cutpoint(df1, time = "PFI.time", event = "PFI",
                                                        variables = c("Gene"))$cutpoint[,1] ,
                                "High","Low")
      )
      
      fit <- survfit(Surv(PFI.time, PFI) ~ Risk, data =df1 )
      print(plot_km(fit,df1,title = "Score" ) ) 
      ggsave(filename = paste0("./out/",km_method,"_",names(cohort[y]),"_",i,"_PFI.pdf"), height = 3.5, width = 3.5)
      
    }
   
  } 
  else{
    next
  }
  
  ###1. COX-all cohort
  #uni_cox_fun(genes_exp_clinical_1,
  #            c(genes_1,"Jak-STAT_pathway","Core_Jak-STAT_pathway") )
  #ggsave(filename = paste0("./out/",names(cohort[y]),"_genes_pathways_COX.pdf"), height = 3.5, width = 4.5,family="serif")

  # ###2. cor with interest gene
#   df_cor <-genes_exp_clinical_1[,c(genes_1,input_gene,"Type")]
#   ##cor  interest_gene with input_gene
#   for (gene in input_gene){
#     print(gene)
#     df_cor1<-df_cor;df_cor1$gene <-df_cor1[,gene];df_cor1$interest_gene <-df_cor1[,interest_gene]
#     library(ggplot2);library(ggpubr);library(ggsci)
#     ggplot(df_cor1, aes(x=interest_gene, y=gene)) + 
#       #xlim(-20,15) + ylim(-15,10) +
#       labs(x = interest_gene, y = gene) +
#       geom_point(aes(color=Type),size = 2,alpha=0.5) +
#       geom_smooth(method ='lm', size=0.5,color="grey30") + #skyblue
#       stat_cor(method = "spearman",size = 6,color="red") +
#       scale_colour_manual(values = colors[c(6,3)] ) +
#       theme_bw() + 
#       theme(axis.text.x = element_text(size = 16), 
#             axis.text.y = element_text(size = 16), 
#             axis.title.x = element_text(size = 18), 
#             axis.title.y = element_text(size = 18),
#             legend.title = element_blank(),
#             legend.text = element_text(size = 12))
#     ggsave(paste0("./out/Cor_dotplot_",names(cohort[y]),"_",interest_gene,"_",gene,".pdf"),family="serif",width = 5,height = 3.5)
#     
# }
#   #ggcor
#   library(ggcorrplot)
#   cormtcars <- round(cor(subset(df_cor,select=-Type)), 2) #round()函数自定义小数点后位数
#   # ggcorrplot(cormtcars,lab=T)
#   # ggsave(paste0("./out/Cor_dotplot_",names(cohort[y]),".pdf"),family="serif",width = 6.5,height = 6.5)
#   ggcorrplot(cormtcars,method = "circle",lab=T)#默认是方形，修改为圆形显示，并标上相关系数
#   ggsave(paste0("./out/Cor_dotplot_",names(cohort[y]),"-1.pdf"),family="serif",width = 6.5,height = 6.5)
#   #使用ggcorrplot包的cor_pmat()函数计算p值：
#   pmtcars <- cor_pmat(subset(df_cor,select=-Type))
#   ggcorrplot(cormtcars,hc.order = T,method = "circle",  #分等级聚类重排矩阵
#              type="full",# upper full lower但是会漏元素
#              ggtheme = theme_void(base_size = 14),#theme_bw(),#ggplot2::theme_void(base_size = 15), #主题修改
#              colors = c("CornflowerBlue","white","Salmon"), #自定义颜色，看自己喜欢，或是参考好看的文献Figure用法。
#              lab = T,lab_size = 3.5,    #相关系数文本字体大小
#              tl.cex = 12,             #坐标轴字体大小
#              p.mat = pmtcars,         #添加显著性信息
#              sig.level = 0.001,        #显著性水平
#              insig = "pch",pch = 4,#insig = "blank"                 #不够显著的色块进行标记，pch表示选择不同的标记方法，可以尝试其他数字表示什么标记方法
#              pch.cex = 5)            #不显著标记的大小，使用insig = "blank"将不显著的空白处理
#   ggsave(paste0("./out/Cor_dotplot_",names(cohort[y]),"-p0.001.pdf"),family="serif",width = 6.5,height = 6.5)
#   ##显示一半
#   library(corrplot)
#   mycol2 <- colorRampPalette(c("skyblue","white", "pink"), alpha = TRUE)
#   pdf(paste0("./out/Cor_dotplot_",names(cohort[y]),"-p0.001-corrplot.pdf"),family="serif",width = 6,height = 6)
#   # corrplot(cormtcars, method = c('circle'), type = c('lower'), #square
#   #          col = mycol2(100),
#   #          outline = 'grey', #是否为图形添加外轮廓线，默认FLASE，可直接TRUE或指定颜色向量
#   #          order = c('AOE'), #排序/聚类方式选择："original", "AOE", "FPC", "hclust", "alphabet"
#   #          diag = FALSE, #是否展示对角线结果，默认TRUE
#   #          tl.cex = 1.2, #文本标签大小
#   #          tl.col = 'black', #文本标签颜色
#   #          addgrid.col= 'grey' #格子轮廓颜色
#   # )
#   corrplot(cormtcars, p.mat =pmtcars, sig.level = 0.001,
#            col = mycol2(100),addgrid.col= 'white',pch.cex = 1,
#            tl.col = 'black',tl.srt = 45)
#   #ggsave(file='./out/corrplot-p-1.pdf',width = 4,height = 4,family="serif")
#   dev.off()
  
if(F){
  ###3.exp and km in cohort
  for (i in c(genes_1,input_gene) ){
    print(i)
    ##临床亚型
    index <- intersect(colnames(genes_exp_clinical_1),c("pathologic_stage","T_stage","N_stage","M_stage","IHC","Type") )
    for (j in index){
      print(j)
      ##exp-boxplot
      Gene = i
      df <- genes_exp_clinical_1
      df$Gene <-df[,i]
      df$tissue <-df[,j]
      df <-df[!(df$tissue=="Unkonwn"),]  #Unknown
      ##
      ggplot(df,aes(x=tissue,y=Gene,color=tissue ))+ #reorder(cancer,-STAT4)
        #geom_boxplot()+
        # geom_boxplot(#width=0.55, #alpha=0.7,
        #   #position=position_dodge(width=0.8),
        #   size=0.05,outlier.colour = NA)+
        # scale_fill_manual(values = c("skyblue","pink") )+ #c( "#56B4E9","#CC79A7")
        # geom_point(alpha=0.5,size=1.5,
        #            position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
        #                                          jitter.height = 0,
        #                                          dodge.width = 0.8))+
        geom_boxplot(alpha=1,width=0.45,fill=NA,position=position_dodge(width=0.8),
                     size=0.4,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                     outlier.stroke = 0.5)+
        geom_violin(alpha=0.2,width=0.9,
                    position=position_dodge(width=0.8),
                    linewidth=0.25)+
        scale_fill_manual(values = colors)+ # c("#7EC7A7","#EDA065"),"#7EC7A7" ,"#66CCFF"
        scale_color_manual(values = colors )+
        #geom_text(aes(y=STAT4+0.3,label=Freq,color="grey10",size = 2) )+
        #scale_y_continuous(expand = c(0,0),limits = c(0,14))+
        labs(x=NULL,y=Gene)+
        EnvStats::stat_n_text()+ #library(EnvStats)
        theme_classic2()+#theme_bw()+
        theme(legend.title = element_blank(),
              legend.position = 'none',
              axis.text = element_text(colour = 'black',size = 12),
              legend.text = element_text(colour = 'black',size = 10),
              axis.title.y= element_text(size = 14),
              #axis.text = element_blank(),
              #axis.title = element_blank(),
              #panel.border = element_blank(),
              panel.grid  = element_blank()
              #plot.title = element_text(hjust = 0.5,size = 20) 
        ) +#rotate_x_text(angle = 90)+
        stat_compare_means(#comparisons = my_comparisons,#aes(group = tissue), 
          label = "p.format",label.x = 1.2,show.legend = F,size=4, #p.format p.signif
          #method = "wilcox.test", 
          hide.ns = F) #, label.y.npc = "top"
      wid=length(unique(df$tissue))
      #    ggsave(paste0("./out/Subtype_",names(cohort[y]),"_",i,"_",j,".pdf"),family="serif",width = wid*1,height = 2.5)
      
      ##km
      #OS overall survival
      for(km_method in c("median","bestcutoff") ){
        ifelse(km_method == "median",
               df$Risk <- ifelse(df$Gene>median(df$Gene) ,"High","Low"),
               df$Risk <- ifelse(df$Gene> surv_cutpoint(df, time = "OS.time", event = "OS",
                                                        variables = c("Gene"))$cutpoint[,1] ,
                                 "High","Low")
        )
        
        fit <- survfit(Surv(OS.time, OS ) ~ Risk, data =df )
        #plot_km(fit,df,title = i) 
        #ggsave(filename = paste0("./out/",km_method,"_",names(cohort[y]),"_",i,"_",j,"_OS.pdf"), height = 3.5, width = 3.5)
      }
      
    }
  }  #genes_1,
  
  ###4.km in IHC shubtypes
  for (i in c(input_gene,interest_gene) ){
    print(i)
    ##临床亚型
    index <- c("Type","IHC")#intersect(colnames(genes_exp_clinical_1),c("pathologic_stage","T_stage","N_stage","M_stage","IHC","Type") )
    for (j in index){
      print(j)
      ##exp-boxplot
      Gene = i
      df <- genes_exp_clinical_1
      df$Gene <-df[,i]
      df$tissue <-df[,j]
      df <-df[!(df$tissue=="Unkonwn"),]  #Unknown
      ##
      
      for (k in unique(df[,j]) ){
        print(k)
        df1 <-df[df$tissue==k,]
        #
        for(km_method in c("median","bestcutoff") ){
          ifelse(km_method == "median",
                 df1$Risk <- ifelse(df1$Gene>median(df1$Gene) ,"High","Low"),
                 df1$Risk <- ifelse(df1$Gene> surv_cutpoint(df1, time = "OS.time", event = "OS",
                                                            variables = c("Gene"))$cutpoint[,1] ,
                                    "High","Low")
          )
          
          #fit <- survfit(Surv(OS.time, OS) ~ Risk, data =df1 )
          #print(plot_km(fit,df1,title = i)) 
          #ggsave(filename = paste0("./out/",km_method,"_",names(cohort[y]),"_",i,"_",j,"_",k,"_OS.pdf"), height = 3.5, width = 3.5)
          
          ##subtype COX OS
          uni_cox_fun(df1,
                      c(genes_1,"Jak-STAT_pathway","Core_Jak-STAT_pathway") )
          #ggsave(filename = paste0("./out/",names(cohort[y]),"_",j,"_",k,"_COX.pdf"), height = 3.5, width = 4.5,family="serif")
          
          
        }
        
        #OS and km in TNBC
        
      }
    }
  }  
  #end run!
}
  
  
}

#gene km
for (y in 1:length(cohort) ){ #1:length(cohort)
  print(y);print(names(cohort[y]))

  #run!
  genes_exp_clinical <-meta_all[[y]]
  genes_exp_clinical$sample <-rownames(genes_exp_clinical)
  #add gene exp
  exp <-as.data.frame(t(cohort[[y]][,-c(1:2)] ))  #as.data.frame(t(cohort[[y]]))#panc_22#genes_exp_clinical[,1:152]
  genes_exp <-as.data.frame(t(exp[rownames(exp) %in% genes,])) #input your genes
  genes_exp$sample <-rownames(genes_exp)
  genes_exp_clinical_1 <-merge(genes_exp_clinical,genes_exp,by="sample")
  #genes_exp_clinical_1$OS <- as.numeric( gsub(" ","",genes_exp_clinical_1$OS ) )
  
  df1 <- genes_exp_clinical_1
  for (i in c("STAT4")){ #KM-分组指标
    print(i)
    df1$Gene <- df1[,i]
    #start km
    for(km_method in c("median","bestcutoff") ){
      ##OS
      ifelse(km_method == "median",
             df1$Risk <- ifelse(df1$Gene>median(df1$Gene) ,"High","Low"),
             df1$Risk <- ifelse(df1$Gene> surv_cutpoint(df1, time = "OS.time", event = "OS",
                                                        variables = c("Gene"))$cutpoint[,1] ,
                                "High","Low")
      )
      
      fit <- survfit(Surv(OS.time, as.numeric(OS)) ~ Risk, data =df1 )
      print(plot_km(fit,df1,title = i ) ) 
      ggsave(filename = paste0("./out/",km_method,"_",names(cohort[y]),"_",i,"_OS.pdf"), height = 3.5, width = 3.5)
      
    }
    
    if(names(cohort[y])=="TCGA"){
      df1$PFI.time <- df1$PFI.time/365
      df1 <- df1[,c("Gene","PFI","PFI.time")]
      df1 <-na.omit(df1)
      for(km_method in c("median","bestcutoff") ){
        #PFI/PFS
        ifelse(km_method == "median",
               df1$Risk <- ifelse(df1$Gene>median(df1$Gene) ,"High","Low"),
               df1$Risk <- ifelse(df1$Gene> surv_cutpoint(df1, time = "PFI.time", event = "PFI",
                                                          variables = c("Gene"))$cutpoint[,1] ,
                                  "High","Low")
        )
        
        fit <- survfit(Surv(PFI.time, PFI) ~ Risk, data =df1 )
        print(plot_km(fit,df1,title = i ) ) 
        ggsave(filename = paste0("./out/",km_method,"_",names(cohort[y]),"_",i,"_PFI.pdf"), height = 3.5, width = 3.5)
        
      }
      
    } 
    else{
      next
    }
    ##end km
  }
}  

#2.0 CIBERSORT----
for (y in 1:length(cohort) ){ #1:length(cohort)
  print(y);print(names(cohort[y]))
  
  #run!
  genes_exp_clinical <-meta_all[[y]]
  genes_exp_clinical$sample <-rownames(genes_exp_clinical)
  exp <- as.data.frame(t(cohort[[y]]))#panc_22#genes_exp_clinical[,1:152]
  
  ###score
  #exp <- tibble::column_to_rownames(exp, var = "sample")
  library(GSVA)
  #gsva_matrix<- gsva(as.matrix(exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
  ##First we should build a parameter object for the desired methodology.R
  gsvaPar <- ssgseaParam(exprData = as.matrix(exp), 
                         geneSets = genelist,
                         normalize = TRUE)
  ##Second, we call the gsva() function with the parameter object as first argument. 
  gsva_matrix <- gsva(gsvaPar, verbose = FALSE)
  gsva_matrix<-as.data.frame(gsva_matrix);gsva_matrix<-as.data.frame(t(gsva_matrix ))
  gsva_matrix$sample <-rownames(gsva_matrix)
  ###merge in meta
  genes_exp_clinical_1 <-merge(genes_exp_clinical,gsva_matrix,by="sample")
  rownames(genes_exp_clinical_1)<-genes_exp_clinical_1$sample
  #merge gene exp
  genes_exp <-as.data.frame(t(exp[rownames(exp) %in% genes,]))
  genes_exp$sample <-rownames(genes_exp)
  #genes_exp$ID <-genes_exp$sample#rownames(genes_exp)
  genes_exp_clinical_1 <-merge(genes_exp_clinical_1,genes_exp,by="sample")
  genes_exp_clinical_1$OS <- as.numeric( gsub(" ","",genes_exp_clinical_1$OS ) )
  
  ##run cibersort
  if(T){
    
    ###prepare data
    exp_input <-exp#as.data.frame(t(cohort$TCGA))#GENE
    exp_input<-2^exp_input-1
    # exp_input$ID <-rownames(exp_input)
    # library(dplyr);library(tidyverse)
    # exp_input <-exp_input %>% select(c("ID"),everything() ) 
    library(IOBR)
    im_cibersort <- deconvo_tme(eset = exp_input,# slow >20min
                                method = "cibersort", arrays = FALSE, perm = 1000)#1000
    colnames(im_cibersort)<-gsub("_CIBERSORT","",colnames(im_cibersort))
    colnames(im_cibersort)[1]<-"sample"
    #merge in metadata
    genes_exp_clinical_1<-merge(genes_exp_clinical_1,im_cibersort[1:23],by="sample")
    #cibersort cell types
    #colnames(im_cibersort)[2:23]
    genes_exp_clinical_1$Group <-ifelse(genes_exp_clinical_1$`Jak-STAT_pathway`>mean(genes_exp_clinical_1$`Jak-STAT_pathway`),"High","Low")
    
    ###plot pheatmap
    library(pheatmap);library(RColorBrewer)
    #样本注释
    immu_cibersort_score<-genes_exp_clinical_1[,c(colnames(im_cibersort)[2:23],"Group")]
    rownames(immu_cibersort_score)<-rownames(genes_exp_clinical_1)
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
    
    # pheat<-pheatmap(na.omit(t(immu_cibersort_score[,1:22])) , #t(immu_cibersort_score[,1:22])
    #                 show_colnames = F,cluster_rows = T,cluster_cols = F,
    #                 scale = "column",
    #                 #cellwidth=3,cellheight = 1,
    #                 annotation_col =immu_cibersort_score[c("Group")],annotation_colors = ann_colors,
    #                 gaps_col= table(immu_cibersort_score$Group)[[1]] ,border_color = NA, #gaps_row= c(7), 
    #                 annotation_names_row = F, annotation_names_col = F,
    #                 color=alpha(rev(RColorBrewer::brewer.pal(10,"Spectral") ), 0.7) #RdYlGn Spectral   #PiYG "BrBG" 褐浅蓝  "PiYG"紫绿
    #                 #colorRampPalette( c("skyblue", "white", "firebrick3"))(10),#
    # )
    # #save_pheatmap_pdf(pheat,"./out/cibersort_pheatmap.pdf",6,3.5) 
    # 
    # pdf("./out/cibersort_pheatmap.pdf",6,3.5,family = "serif") 
    # print(pheat)
    # dev.off()
    
    #计算显著性后重新绘制热图!#T检验显著性表格
    #exp_int[,i]
    High <- rownames(immu_cibersort_score[immu_cibersort_score$Group =="High",])#==
    Low <-rownames(immu_cibersort_score[immu_cibersort_score$Group !="High",])#!=
    
    pvalue_interest <-data.frame() #T.test
    exp_int<- immu_cibersort_score[1:22]#exprSet[rownames(exprSet) %in% gene_more,] #raw FPKM #gene_more 兴趣基因集
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
                      #Group=c("High")
                      Sign.=c("ns"="grey","*"="pink",'**'="#D33F6A","***"="red" ) )
    #plot
    int_p <-pheatmap(t(exp_int),
                     show_colnames = F,cluster_rows = T,cluster_cols = F,
                     scale = "column",
                     #cellwidth=3,cellheight = 1,
                     annotation_col =immu_cibersort_score[c("Group")],#annotation_colors = ann_colors,
                     gaps_col= table(immu_cibersort_score$Group)[[1]] ,border_color = NA, #gaps_row= c(7), 
                     annotation_names_row = F, annotation_names_col = F,
                     annotation_row = anno_row_p,annotation_colors = ann_colors_p,
                     color=alpha(rev(RColorBrewer::brewer.pal(10,"PiYG") ), 0.8) #PiYG RdYlGn Spectral   #PiYG "BrBG" 褐浅蓝  "PiYG"紫绿
                     #colorRampPalette( c("skyblue", "white", "firebrick3"))(10),#
    )
    pdf(paste0("./out/cibersort_pheatmap_sign_",names(cohort[y]),".pdf"),7,4.5,family = "serif") 
    print(int_p)#names(cohort[y])
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
  }
  #end cibersort
}

#2.1 ESTIMATE ----
for (y in 1 ){ #1:length(cohort)
  print(y);print(names(cohort[y]))
  
  #run!
  genes_exp_clinical <-meta_all[[y]]
  genes_exp_clinical$sample <-rownames(genes_exp_clinical)
  exp <- as.data.frame(t(cohort[[y]]))#panc_22#genes_exp_clinical[,1:152]
  
  ###score
  #exp <- tibble::column_to_rownames(exp, var = "sample")
  library(GSVA)
  #gsva_matrix<- gsva(as.matrix(exp), list,method='ssgsea',kcdf='Gaussian',abs.ranking=TRUE)
  ##First we should build a parameter object for the desired methodology.R
  gsvaPar <- ssgseaParam(exprData = as.matrix(exp), 
                         geneSets = genelist,
                         normalize = TRUE)
  ##Second, we call the gsva() function with the parameter object as first argument. 
  gsva_matrix <- gsva(gsvaPar, verbose = FALSE)
  gsva_matrix<-as.data.frame(gsva_matrix);gsva_matrix<-as.data.frame(t(gsva_matrix ))
  gsva_matrix$sample <-rownames(gsva_matrix)
  ###merge in meta
  genes_exp_clinical_1 <-merge(genes_exp_clinical,gsva_matrix,by="sample")
  rownames(genes_exp_clinical_1)<-genes_exp_clinical_1$sample
  #merge gene exp
  genes_exp <-as.data.frame(t(exp[rownames(exp) %in% genes,]))
  genes_exp$sample <-rownames(genes_exp)
  #genes_exp$ID <-genes_exp$sample#rownames(genes_exp)
  genes_exp_clinical_1 <-merge(genes_exp_clinical_1,genes_exp,by="sample")
  genes_exp_clinical_1$OS <- as.numeric( gsub(" ","",genes_exp_clinical_1$OS ) )
  genes_exp_clinical_1$Group <-ifelse(genes_exp_clinical_1$`Jak-STAT_pathway`>mean(genes_exp_clinical_1$`Jak-STAT_pathway`),"High","Low")
  
  ##run ESTIMATE
  library(estimate)
  estimate <- function(dat,pro){
    input.f=paste0(pro,'_estimate_input.txt')
    output.f=paste0(pro,'_estimate_gene.gct')
    output.ds=paste0(pro,'_estimate_score.gct')
    write.table(dat,file = input.f,sep = '\t',quote = F)
    library(estimate)
    filterCommonGenes(input.f=input.f,
                      output.f=output.f ,
                      id="GeneSymbol")
    estimateScore(input.ds = output.f,
                  output.ds=output.ds,
                  platform="illumina")   ## platform
    scores=read.table(output.ds,skip = 2,header = T,check.names = F)
    rownames(scores)=scores[,1]
    scores=t(scores[,3:ncol(scores)])
    library(stringr)
    rownames(scores)=str_replace_all(rownames(scores),'[.]','-') # 这里TCGA样本名里面的-变成.了，进行恢复
    write.csv(scores,file=paste0(pro,"_Stromal_Immune_ESTIMATE.Score.csv")) # 这一步是我增加的
    return(scores)
  }
  exp_input <-exp#as.data.frame(t(cohort$TCGA))#GENE
  exp_input<-2^exp_input-1
  scores=estimate(exp_input,names(cohort[y]) ) #slow
  scores<-as.data.frame(scores);scores$sample <-rownames(scores)
  genes_exp_clinical_1<-merge(genes_exp_clinical_1,scores,by="sample")
  write.csv(genes_exp_clinical_1,paste0("./out/",names(cohort[y]),"_genes_exp_clinical_1.csv") )
  #plot ESTIMATE
  library(reshape2)
  bar_df <-melt(genes_exp_clinical_1[,c("StromalScore","ImmuneScore","ESTIMATEScore","Group")],id.vars=c("Group") )
  library(ggplot2);library(ggpubr)
  p_method <-c("wilcox.test","t.test")
  for (j in p_method){
    ggplot(bar_df,#
           aes(variable,scale(value),color=Group))+  #age_grou
      # geom_point(alpha=0.7,size=0.2,
      #            position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
      #                                          jitter.height = 0,
      #                                          dodge.width = 0.7))+
      geom_boxplot(alpha=1,width=0.7,fill=NA,
                   position=position_dodge(width=0.8),
                   size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                   outlier.stroke = 0.5)+
      geom_violin(alpha=0.2,width=0.9,
                  position=position_dodge(width=0.8),
                  size=0.25)+
      labs(x="",y="Scaled Score")+
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
                         label.x =0.75)#,label.y =max(bar_df$value)-1000)#
    ggsave(paste0("./out/ESTIMATE_boxplot_",j,"_",names(cohort[y]),".pdf"),width = 3,height = 4,family="serif")
    
    ##single boxplot
    for (i in as.character(unique(bar_df$variable)) ){
      for (j in p_method[1] ){ #wilcox.test
        ggplot(bar_df[bar_df$variable %in% i,],#
               aes(Group,scale(value),color=Group))+
          # geom_point(alpha=0.1,size=0.5,
          #            position=position_jitterdodge(jitter.width = 0.75,#0.45  1.2
          #                                          jitter.height = 0,
          #                                          dodge.width = 0.8))+
          geom_boxplot(alpha=1,width=0.75,fill=NA,
                       position=position_dodge(width=0.8),
                       size=0.2,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                       outlier.stroke = 0.5)+
          geom_violin(alpha=0.2,width=0.75,
                      position=position_dodge(width=0.8),
                      size=0.25)+
          labs(title = "",x="",y=paste0("Scaled ",i))+
          scale_color_manual(values =c("pink","skyblue") )+ # #blue c("#0066CC","#66CCFF","#9DC3E6")c("#EDA065","#66CCFF","#7EC7A7")
          theme_bw() +#labs(color="",title = gene,y="Scaled expression",x="")+ #ylab("STAT4 expression" )+xlab("")+
          theme(text = element_text(size=16),axis.title  = element_text(size=14), axis.text = element_text(size=14), 
                axis.text.x = element_text(angle = 90,hjust = 1),
                panel.grid = element_blank(),legend.key = element_blank() ,legend.position="none" ) + #
          stat_compare_means(method = j,size=8, #"wilcox.test"
                             show.legend= F,label = "p.signif",#p.signif p.format
                             label.x =1.5)#,label.y = max(bar_df[bar_df$variable %in% i,][3]) -max(bar_df[bar_df$variable %in% i,][3])/20)
        #,label.y =max(bar_df[bar_df$variable %in% i,]$value)-0.01)#label.x =0.75,
        ggsave(paste0("./out/ESTIMATE_boxplot_",names(cohort[y]),"_",i,"_",j,".pdf"),width = 2.5,height = 3,family="serif")
        print(i)
      }
    }
    
  }
}

#1.2 OS and km in TNBC/non-TNBC

#过滤TCGA-TNBC >30 Day







#3.0 GEO: T vs N ----

##3.1 GEO N vs T and KM----
#GSE37751(N+T=108),GSE38959(TNBC+N=30+19)
###GSE37751----
#load("~/R/my_projects/JAK_STAT_BC/2.0_BC_bulk_cohort/out/GSE37751-BRCA-20250504.Rdata")

#prepare data
meta <-GSE37751_meta
meta$tissue<-ifelse(meta$tissue_type==" Tumor","Tumor","Normal")
table(meta$tissue)#N 47; T 61
exp <-as.data.frame(t(GSE37751[,-c(1:2)])) #colname = sample ID
#exp <-log2(exp+1)
boxplot(exp[1:100,1:10])
cohort <-"GSE37751"
#rm(exp);rm(meta)
#my_ssgsea_plot_NvsT(exp,meta,genelist,cohort)

###GSE38959----
#load("~/R/my_projects/JAK_STAT_BC/2.0_BC_bulk_cohort/out/GSE38959-BRCA-N17-T30.Rdata")
#prepare
table(meta$tissue)#N 17; T 30
exp <-as.data.frame(t(exp)) #colname = sample ID
exp <-log2(exp+1)
boxplot(exp[1:100,1:10])
cohort <-"GSE38959"

#km function
plot_km <-function(fit,data){
  print(
    ggsurvplot_list(fit,data = data,pval = T,pval.method = TRUE,##是否添加P值
                    conf.int = F,### 是否添加置信区间
                    legend.title = "", # 设置图例标题title = cancer,
                    legend.labs = c("High", "Low"), # 指定图例分组标签
                    risk.table = F, # 是否添加风险表
                    risk.table.col = "strata", 
                    censor=F,size=0.5,
                    surv.scale="percent",
                    ###linetype = "strata",
                    #surv.median.line = "hv", # 是否添加中位生存线
                    risk.table.y.text.col = F,risk.table.y.text = FALSE,
                    ggtheme = theme_bw(base_family = "serif")+theme(legend.text = element_text(colour = c("red", "blue")))
                    +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
                    +theme(plot.title = element_text(hjust = 0.5,size = 16),axis.title.y.left = element_text(size = 16,vjust = 1),axis.title.x.bottom = element_text(size = 16,vjust = 0))
                    +theme(axis.text.x.bottom = element_text(size = 12,vjust = -0.8,colour = "black")) #,face = "bold"
                    +theme(axis.text.y.left = element_text(size = 12,vjust = 0.5,hjust = 0.5,angle = 90,colour = "black"))
                    +theme(legend.title = element_text(family = "Times",colour = "black",size = 12))
                    +theme(legend.text = element_text(family = "Times",colour = "black",size =12)), # Change ggplot2 theme
                    palette = "lacent",#c("#bc1e5d", "#0176bd"),##如果数据分为两组，就写两种颜色，三组就写三种颜色，具体的颜色搭配参考上一期发布的R语言颜色对比表
                    #xlim=c(10,5000),##x轴的阈值，根据随访数据的天数进行更改
                    xlab = "Years")##随访的时间时天，就写Days，月份就写Months
  ) 
  #print(p1)
}

##calculate--funtion   #使用前曲线注释ggsave
my_ssgsea_plot_NvsT <-function(exp=exp,meta=meta,genelist=genelist,cohort=cohort){
  #calculate
  library(GSVA)
  gsvaPar <- ssgseaParam(exprData = as.matrix(exp), 
                         geneSets = genelist,
                         normalize = TRUE)
  ##Second, we call the gsva() function with the parameter object as first argument. 
  gsva_matrix <- gsva(gsvaPar, verbose = FALSE)
  gsva_matrix<-as.data.frame(gsva_matrix);gsva_matrix<-as.data.frame(t(gsva_matrix ))
  gsva_matrix$ID <-rownames(gsva_matrix)
  #merge in meta
  genes_exp_clinical <-merge(gsva_matrix,meta,by="ID")
  #merge genes
  genes_exp <-exp[rownames(exp) %in% genes,]#151
  genes_exp <- as.data.frame(t(genes_exp) )
  genes_exp$ID <-rownames(genes_exp) 
  genes_exp_clinical_1 <- merge(genes_exp,genes_exp_clinical,by="ID")
  
  ##plot
  for (i in c(genes_1,input_gene) ){
    print(i)
    Gene = i
    genes_exp_clinical_1$Gene <-genes_exp_clinical_1[,i]
    ##
    ggplot(genes_exp_clinical_1,aes(x=tissue,y=Gene,color=tissue ))+ #reorder(cancer,-STAT4)
      #geom_boxplot()+
      # geom_boxplot(#width=0.55, #alpha=0.7,
      #   #position=position_dodge(width=0.8),
      #   size=0.05,outlier.colour = NA)+
      # scale_fill_manual(values = c("skyblue","pink") )+ #c( "#56B4E9","#CC79A7")
      geom_point(alpha=0.5,size=1.5,
                 position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
                                               jitter.height = 0,
                                               dodge.width = 0.8))+
      geom_boxplot(alpha=1,width=0.45,fill=NA,position=position_dodge(width=0.8),
                   size=0.4,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                   outlier.stroke = 0.5)+
      geom_violin(alpha=0.2,width=0.9,
                  position=position_dodge(width=0.8),
                  linewidth=0.25)+
      scale_fill_manual(values =  c("#7EC7A7","#EDA065"))+ #,"#7EC7A7" ,"#66CCFF"
      scale_color_manual(values =  c("#7EC7A7","#EDA065"))+
      #geom_text(aes(y=STAT4+0.3,label=Freq,color="grey10",size = 2) )+
      #scale_y_continuous(expand = c(0,0),limits = c(0,14))+
      labs(x=NULL,y=Gene)+
      EnvStats::stat_n_text()+ #library(EnvStats)
      theme_classic2()+#theme_bw()+
      theme(legend.title = element_blank(),
            #legend.position = 'none',
            axis.text = element_text(colour = 'black',size = 12),
            legend.text = element_text(colour = 'black',size = 10),
            axis.title.y= element_text(size = 14),
            #axis.text = element_blank(),
            #axis.title = element_blank(),
            #panel.border = element_blank(),
            panel.grid  = element_blank()
            #plot.title = element_text(hjust = 0.5,size = 20) 
      ) +#rotate_x_text(angle = 90)+
      stat_compare_means(#comparisons = my_comparisons,#aes(group = tissue), 
        label = "p.format",label.x = 1.2,show.legend = F,size=4, #p.format p.signif
        method = "wilcox.test", hide.ns = F) #, label.y.npc = "top"
    #ggsave(paste0("./out/",cohort,"_",i,"_TvsN.pdf"),family="serif",width = 3.5,height = 2.5)
    
    if( c("OS") %in% colnames(genes_exp_clinical_1) ){ #判断是否有预后信息
      ##KM
      df1= genes_exp_clinical_1[genes_exp_clinical_1$tissue=="Tumor",] #删除NA？
      df1$Gene <-df1[,i]
      #df1 <-df1[df1$cancer==cancer,] 
      #OS overall survival
      for(km_method in c("median","bestcutoff") ){
        print(km_method)
        ifelse(km_method == "median",
               df1$Risk <- ifelse(df1$Gene>median(df1$Gene) ,"High","Low"),
               df1$Risk <- ifelse(df1$Gene> surv_cutpoint(df1, time = "OS.time", event = "OS",
                                                          variables = c("Gene"))$cutpoint[,1] ,
                                  "High","Low")
        )
        ##
        fit <- survfit(Surv(OS.time, OS) ~ Risk, data =df1 )
        plot_km(fit,df1) 
        ggsave(filename = paste0("./out/",km_method,"_",i,"_",cohort,"_OS.pdf"), height = 3.5, width = 3.5)
      }#  end for km_method
    } #if end
    else{
      print("No survival information in meta !")
      next#break
    }

    
  } #end for i
  
}
my_ssgsea_plot_NvsT(exp,meta,genelist,cohort)



