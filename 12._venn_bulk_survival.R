
#venn JAK-STAT in tumor,epi and Tcell

if(T){
  rm(list = ls())
  setwd("~/R/my_projects/JAK_STAT_BC/12._venn_bulk_survival")
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
  
}


#1.0 load data----

genelist <- list(
  Tumor = c("CISH", "PIAS4", "IL22RA2", "AKT1", "IFNL2", "STAM2", "EPO", 
            "IL19", "SOCS2", "CBLB", "IL13RA1", "STAT5B", "PIK3CB", "IFNGR2", 
            "IL12RB2", "PTPN11", "IFNAR1", "PRLR", "CBLC", "PIK3CG", "PTPN6", 
            "SPRY3", "CSF2RB", "IL20RA", "IFNL3", "TYK2", "BCL2L1", "CCND1", 
            "PIK3R2", "IL2RA", "CTF1", "PIAS3", "PIK3R3", "CNTFR", "CSF3R"),
  Epi = c("CSF3", "CCND2", "IL6R", "STAT3", "TSLP", "IL4R", "PIAS1", "SPRY4", 
          "IL6", "SPRED1", "PIK3CD", "LIFR", "CBL", "PIK3CA", "IL11RA", "JAK1", 
          "IL11"),
  Tcell = c("IL7R", "LEPR", "IL13", "IFNGR1", "OSM", "EPOR", "IL5RA", "CSF2", 
            "LIF", "IL15", "IL2")
)

#load exp and meta
cohort <-readRDS("~/R/my_projects/JAK_STAT_BC/2.0_BC_bulk_cohort/bulk_data/BRCA_bulk_TCGA_METABRIC_SCAN_B.rds")
meta_all<-readRDS("~/R/my_projects/JAK_STAT_BC/2.0_BC_bulk_cohort/bulk_data/BRCA_bulk_TCGA_METABRIC_SCAN_B_meta.rds")

y="TCGA"
exp<-cohort[[y]][,-c(1:2)]
exp<- as.data.frame(t(exp))

meta_all[[paste0(y,"_meta")]]$clinical_stage  <- meta_all[[paste0(y,"_meta")]]$pathologic_stage
#

#2.0 ssGSEA score----
library(GSVA)
##First we should build a parameter object for the desired methodology.R
gsvaPar <- ssgseaParam(exprData = as.matrix(exp), #sample in column
                       geneSets = genelist,
                       normalize = TRUE)
##Second, we call the gsva() function with the parameter object as first argument. 
gene_score <- gsva(gsvaPar, verbose = FALSE)
gene_score <-as.data.frame(t(gene_score))
gene_score$ID <-rownames(gene_score)


#3.0 COX----

##prepare meta
meta <-meta_all[[paste0(y,"_meta")]]
meta$ID <-rownames(meta)
#meta <-meta[,c(112:121)] #OS.time NA?
# meta$overall_survival <- ifelse(meta$vital_status=="Dead", meta$death_days_to, meta$last_contact_days_to) # define the overall survival time
# meta$overall_survival <- as.numeric(meta$overall_survival)/365
clinical_df<-meta[,c(112:127)] #[,c(113:122)]
#整理分期
clinical_df$tumor_stage<- clinical_df$pathologic_stage
clinical_df = clinical_df %>% 
  mutate(tumor_stage = case_when(
    tumor_stage %in% c('I','IA','IB') ~ 1, 
    tumor_stage %in% c('II','IIA','IIB') ~ 2, 
    tumor_stage %in% c('III','IIIA', 'IIIB', 'IIIC') ~ 3, 
    tumor_stage %in% c('IV') ~ 4, 
    TRUE ~ 0
  )) 

#删除 N M Unknown
clinical_df <- clinical_df %>% 
  filter(!(N_stage == "Unknown") | !(M_stage == "Unknown") ) #906
#T N 转为数值
clinical_df$N_stage <- as.numeric(substr(clinical_df$N_stage ,2,2)) #Unknown = NA
clinical_df$M_stage <- as.numeric(substr(clinical_df$M_stage ,2,2))
clinical_df$T_stage <-as.numeric(clinical_df$T_stage)
#gene_score[,4:6] = lapply(gene_score[,4:6], FUN = function(y){as.numeric(y)})


gene_score_1<-merge(clinical_df,gene_score,by="ID")

##single-COX----
uni_cox_fun <-function(data2,input_gene){
  library(survival);library(survminer)
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
    geom_point(aes(size=HR-2,color=color)  )+ #,color="skyblue" id
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
  print(p ) #+scale_x_continuous(trans = 'log2')
  return(outTab)
  
  
}

outTab <-uni_cox_fun(gene_score_1, colnames(gene_score_1)[18:ncol(gene_score_1)] )
write.csv(outTab, paste0(folder_path,"/single-cox_minor_not_factor.csv") )
ggsave(paste0(folder_path,"/single_minor_not_factor.pdf"),width = 5,height = 3.5,family="serif")

##multi-COX----

library(survival)
library(survminer)
cox = coxph(Surv(OS.time, OS) ~ Age+ T_stage+N_stage+ M_stage+tumor_stage+Tumor+Epi+Tcell
              ,
            data = gene_score_1)
ggforest(model = cox, data = gene_score_1, fontsize = 0.6) # .5

#factor; _not_factor
ggsave(paste0(folder_path,"/multicox_minor_not_factor.pdf"),width = 5,height = 4,family="serif")
write.csv(gene_score_1, paste0(folder_path,"/multicox_minor_not_factor.csv") )



#New plot
multicox <-data.frame(id=rownames(summary(cox)$conf.int),# [[1]]
                      HR=summary(cox)$conf.int[,"exp(coef)"],
                      HR.95L=summary(cox)$conf.int[,"lower .95"],
                      HR.95H=summary(cox)$conf.int[,"upper .95"],
                      pvalue=summary(cox)$coefficients[,"Pr(>|z|)"])
str(multicox)
# multicox[,2:5] <- as.data.frame(lapply(multicox[,2:5] ,function(x){round(x,3)} )) #保留2位数字
#multicox$id <-c("Age","Grade","Size","T_stage","N_stage","IHC_score")
multicox$id[3:7] <-colnames(gene_score)[6:10] #c("Grade","Size","IHC_score")
multicox[,2:4] <- as.data.frame(lapply(multicox[,2:4] ,function(x){round(x,2)} ))
multicox[,5] <- round(multicox[,5],4)#format(round(multicox[,5],4), scientific = T, digits = 2)
multicox$id <-rownames(multicox)

multicox$id[6:length(multicox$id)] <-
  colnames(gene_score_1)[as.numeric(gsub("V","",multicox$id[6:length(multicox$id)] ))]

#plot
cox_res <-multicox#outTab_cox#  
forest_table <- cbind(c(" ", cox_res$id),
                      c("HR (95%CI)", paste0(cox_res$HR," (",cox_res$HR.95L,"-",cox_res$HR.95H,")") ),
                      c("Pvalue", cox_res$pvalue)#pvalue
)##设置标签


csize <- data.frame(mean=c(NA, cox_res$HR),
                    lower=c(NA,  cox_res$HR.95L),
                    upper=c(NA, cox_res$HR.95H))##设置图形数据
csize
library("forestplot")  #https://www.jianshu.com/p/b460e3cd3bc5
pdf("./out/multivariate-cox_p0.05_not_factor.pdf",width=8, height=4) 
forestplot(labeltext = forest_table,#outTab_TF, 
           csize,
           graph.pos = 3,
           txt_gp=fpTxtGp(label=gpar(fontfamily = "serif",cex=1.25,fontsize=12),
                          ticks=gpar(cex=1.2),
                          xlab=gpar(cex = 0.2),
                          title=gpar(cex = 1.2)), # 修改label，ticks 字体和大小
           col=fpColors(box="pink", lines="royalblue", zero = "gray50"),##线条颜色设置 ##1c61b6
           zero=1, cex=0.9, lineheight = "auto", boxsize=0.2, colgap=unit(3,"mm"),
           lwd.ci=1, ci.vertices=TRUE, ci.vertices.height = 0.15,  #ci置信区间线条
           #lower = csize$lower,upper = csize$upper,
           #is.summary = c(TRUE,rep(FALSE,14)),#加粗行
           xlog = F)#坐标轴log压缩
dev.off()
write.csv( forest_table, paste0(folder_path,"/multicox_minor-forest_table_not_factor.csv") )

#4.0 KM----
plot_km <-function(fit,data,title=title){
  print(
    ggsurvplot_list(fit,data = data,pval = T,pval.method = TRUE,##是否添加P值
                    conf.int = F,### 是否添加置信区间
                    legend.title = title,#i, # 设置图例标题title = cancer,
                    #legend.labs = c("High", "Low"), # 指定图例分组标签
                    risk.table = F, # 是否添加风险表
                    risk.table.col = "strata", 
                    censor=F,size=0.5,
                    surv.scale="percent",
                    ###linetype = "strata",
                    #surv.median.line = "hv", # 是否添加中位生存线
                    risk.table.y.text.col = F,risk.table.y.text = FALSE,
                    ggtheme = theme_bw(base_family = "serif",base_size = 14)+theme(legend.text = element_text(colour = c("red", "blue")))
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
##HR
uni_cox_fun <-function(data2,input_gene){
  outTab=data.frame()
  sigGenes=c("futime","fustat")
  data2$futime <- data2$OS.time;data2$fustat <- data2$OS
  
  for(i in input_gene ){ #trainset[,3:ncol(trainset)]
    cox <- coxph(Surv(futime, fustat) ~ data2[,i], data = data2)
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

for (i in names(genelist) ){
  print(i)
  Gene = i
  df <- gene_score_1#[gene_score_1$IHC==j,]
  df$Gene <-df[,i]
  ##
  for(km_method in c("median","bestcutoff") ){
    print(km_method)
    ifelse(km_method == "median",
           df$Risk <- ifelse(df$Gene>median(df$Gene) ,"High","Low"),
           df$Risk <- ifelse(df$Gene> surv_cutpoint(df, time = "OS.time", event = "OS",
                                                    variables = c("Gene"))$cutpoint[,1] ,
                             "High","Low")
    )
    
    fit <- survfit(Surv(OS.time, OS ) ~ Risk, data =df )
    plot_km(fit,df,title = i) 
    ggsave(filename = paste0("./out/",km_method,"_",names(cohort[y]),"_",i,"_OS.pdf"), height = 3.5, width = 3)
    
    #PFS/PFI
    df1 <-df[c("Gene","PFI.time","PFI")]
    df1 <- na.omit(df1)
    ifelse(km_method == "median",
           df1$Risk <- ifelse(df1$Gene>median(df1$Gene) ,"High","Low"),
           df1$Risk <- ifelse(df1$Gene> surv_cutpoint(df1, time = "PFI.time", event = "PFI",
                                                      variables = c("Gene"))$cutpoint[,1] ,
                              "High","Low")
    )
    
    fit <- survfit(Surv(PFI.time/365, PFI ) ~ Risk, data =df1 )
    plot_km(fit,df,title = i) 
    ggsave(filename = paste0("./out/",km_method,"_",names(cohort[y]),"_",i,"_PFS.pdf"), height = 3.5, width = 3)
    
  }
  ##
}
