
#~/R/my_projects/JAK_STAT_BC/1.pan-cancer

if(T){
  rm(list = ls())
  setwd("~/R/my_projects/JAK_STAT_BC/1.pan-cancer")
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
#load("~/R/my_projects/JAK_STAT_BC/1.pan-cancer/pan-cancer.R.RData")

#prepare data----
project = "Jak_STAT_BC"
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

my_cancer<-c("BRCA")
genelist<- list("Jak-STAT_pathway"=genes,"Core_Jak-STAT_pathway"=genes_1)
input_gene= c("Jak-STAT_pathway","Core_Jak-STAT_pathway") #genes#c(genes,"H3K36me","5mC" )


#load("~/R/ucsc_pancancer/data/UCSC.RData") #2.2G
#boxplot(panc_22_tumor[1:10]) #肿瘤样本 FPKM
#genes_exp <-panc_22_tumor[rownames(panc_22_tumor) %in% genes,]#
genes_exp <-panc_22[rownames(panc_22) %in% genes,]#151
genes[(genes %in% rownames(genes_exp) ) ==F]   # "IFNLR1" "IFNL2"  "IFNL3"  "IFNL1"

genes_exp <- log2(genes_exp+1)#log2(genes_exp+1)  #11069-151
genes_exp <- as.data.frame(t(genes_exp) )

#处理NA值
genes_exp <-na.omit(genes_exp)#
genes_exp$sample <-rownames(genes_exp) 
genes_exp_clinical <-merge(genes_exp,panc_22_clinical,by="sample") #匹配临床数据9265样本-39
#添加正常样本标记
genes_exp_clinical$tissue <-ifelse(genes_exp_clinical$type=="Tumor",genes_exp_clinical$cancer,paste0(genes_exp_clinical$cancer,"_Normal") ) #加入正常样本信息
table(genes_exp_clinical$tissue)
#删除样本<   30的癌症
table(genes_exp_clinical[genes_exp_clinical$type=="Tumor",]$tissue)
table(genes_exp_clinical$tissue)
# ACC BLCA BRCA CESC CHOL COAD DLBC ESCA  GBM HNSC KICH KIRC KIRP  LGG LIHC LUAD LUSC 
# 79  427 1215  309   45  326   48    3  166  566   91  606  323  529  423  576  552 
# MESO   OV PAAD PCPG PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC  UCS  UVM 
# 87    2  183  187  517  104  265  473    3  139  572  122  190   57   80 
table(genes_exp_clinical$cancer)
#genes_exp_clinical_30 <-genes_exp_clinical[genes_exp_clinical$cancer !=  c("ESCA","OV","STAD"),]
#genes_exp_clinical_30 <- subset(genes_exp_clinical,!(genes_exp_clinical$cancer  %in%   c("ESCA","OV","STAD") ) )#9257-40
#table(genes_exp_clinical_30$cancer)
#删除HNSC异常值
#genes_exp_clinical_30 <-genes_exp_clinical_30[-2670,]#9256 #which.max(genes_exp_clinical_30$STAT4)
#删除正常样本表达？table(genes_exp_clinical$tissue)

##calculate Jak-STAT score
# gene_set<-subset(result1, select = Type)
# gene_set$gene <-rownames(gene_set)
# list<- split(as.matrix(gene_set)[,2], gene_set[,1])
# geneset <- rio::import("ssGSEA-PMID28052254-Signature-28.xlsx",skip = 2)
# geneset <- split(geneset$Metagene,geneset$`Cell type`)

#exp
exp <- panc_22#genes_exp_clinical[,1:152]
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
#merge in meta
genes_exp_clinical_1 <-merge(genes_exp_clinical,gsva_matrix,by="sample")
table(genes_exp_clinical_1$type)
rownames(genes_exp_clinical_1)<-genes_exp_clinical_1$sample
# Normal  Tumor 
# 665   8600 

#1.0 T vs. N---- 
##1.1 all T N----

for (i in input_gene ){
  print(i)
  Gene = i
  genes_exp_clinical_1$Gene <-genes_exp_clinical_1[,i]
  ##
  ggplot(genes_exp_clinical_1,aes(x=type,y=Gene,color=type ))+ #reorder(cancer,-STAT4)
    #geom_boxplot()+
    # geom_boxplot(#width=0.55, #alpha=0.7,
    #   #position=position_dodge(width=0.8),
    #   size=0.05,outlier.colour = NA)+
    # scale_fill_manual(values = c("skyblue","pink") )+ #c( "#56B4E9","#CC79A7")
    # geom_point(aes(fill=type),alpha=0.1,size=1.5,
    #            position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
    #                                          jitter.height = 0,
    #                                          dodge.width = 0.8))+
    geom_boxplot(alpha=1,width=0.45,fill=NA,position=position_dodge(width=0.8),
                 size=0.4,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                 outlier.stroke = 0.5)+
    geom_violin(alpha=0.2,width=0.9,
                position=position_dodge(width=0.8),
                linewidth=0.25)+
    scale_fill_manual(values =  c("#EDA065","#66CCFF"))+ #,"#7EC7A7"
    scale_color_manual(values =  c("#EDA065","#66CCFF"))+
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
      label = "p.signif",label.x = 1.5,show.legend = F,size=4,
      method = "wilcox.test", hide.ns = T) #, label.y.npc = "top"
  ggsave(paste0("./out/",project,"_",i,"_pan_TvsN.pdf"),family="serif",width = 3.5,height = 2.8)
  
  ##
  #1.2 cancer T vs. N----
  
  ggplot(genes_exp_clinical_1,aes(x=cancer,y=Gene,fill=type))+ #reorder(cancer,-STAT4)
    #geom_boxplot()+
    geom_boxplot(#width=0.55, #alpha=0.7,
      #position=position_dodge(width=0.8),
      size=0.05,outlier.colour = NA)+
    scale_fill_manual(values = c("#EDA065","#66CCFF") )+ #c( "#56B4E9","#CC79A7")
    #geom_text(aes(y=STAT4+0.3,label=Freq,color="grey10",size = 2) )+
    #scale_y_continuous(expand = c(0,0),limits = c(0,14))+
    labs(x=NULL,y=Gene)+
    #EnvStats::stat_n_text()+ #library(EnvStats)
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
    ) +rotate_x_text(angle = 90)+
    stat_compare_means(#comparisons = my_comparisons,#aes(group = tissue), 
      label = "p.signif",show.legend = F,size=4, #label.x = 1.5,
      method = "wilcox.test", hide.ns = T) #, label.y.npc = "top"
  ggsave(paste0("./out/",project,"_",i,"_cancers_TvsN.pdf"),family="serif",width = 12,height = 3)
  
}
##1.3 boxplot in single cancer----
data2<-genes_exp_clinical_1#[genes_exp_clinical_1$type=="Tumor",] #tumor  8600
table(data2$cancer)

for (i in c(genes_1,input_gene) ){
  print(i)
  Gene = i
  data2<-genes_exp_clinical_1
  data2$Gene <-data2[,i]
  ##
  for (j in unique(genes_exp_clinical_1$cancer) ){
    print(j)
    data2<-genes_exp_clinical_1[genes_exp_clinical_1$cancer==j,]
    ##
    ggplot(data2,aes(x=type,y=Gene,color=type ))+ #reorder(cancer,-STAT4)
      #geom_boxplot()+
      # geom_boxplot(#width=0.55, #alpha=0.7,
      #   #position=position_dodge(width=0.8),
      #   size=0.05,outlier.colour = NA)+
      # scale_fill_manual(values = c("skyblue","pink") )+ #c( "#56B4E9","#CC79A7")
      # geom_point(aes(fill=type),alpha=0.1,size=1.5,
      #            position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
      #                                          jitter.height = 0,
      #                                          dodge.width = 0.8))+
      geom_boxplot(alpha=1,width=0.45,fill=NA,position=position_dodge(width=0.8),
                   size=0.4,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                   outlier.stroke = 0.5)+
      geom_violin(alpha=0.2,width=0.9,
                  position=position_dodge(width=0.8),
                  linewidth=0.25)+
      scale_fill_manual(values =  c("#EDA065","#66CCFF"))+ #,"#7EC7A7"
      scale_color_manual(values =  c("#EDA065","#66CCFF"))+
      #geom_text(aes(y=STAT4+0.3,label=Freq,color="grey10",size = 2) )+
      #scale_y_continuous(expand = c(0,0),limits = c(0,14))+
      labs(x=NULL,y=Gene,title = j)+
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
        label = "p.signif",label.x = 1.5,show.legend = F,size=4,
        method = "wilcox.test", hide.ns = T) #, label.y.npc = "top"
    ggsave(paste0("./out/",i,"_",j,"_TvsN.pdf"),family="serif",width = 3.5,height = 2.8)
    
  }
  
}


#2.0 expression----
library(pheatmap)
exp<-genes_exp_clinical_1[,c(genes_1,"Jak-STAT_pathway","Core_Jak-STAT_pathway")]
#exp<-t(scale(t(exp)))
exp<-scale(exp)
rownames(exp)<-genes_exp_clinical_1$sample

tumor <-genes_exp_clinical_1[genes_exp_clinical_1$type=="Tumor",]$sample
normal <-genes_exp_clinical_1[genes_exp_clinical_1$type=="Normal",]$sample
#exp<-exp[rownames(exp) %in% c(normal,tumor),]

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

p1 <-pheatmap(exp[c(normal,tumor),],cluster_rows  = F,cluster_cols = T,
         show_colnames = T, show_rownames = F,
         border_color = "white",#cellwidth=10,cellheight = 10,#gaps_col= 38, #table(exp_marker1$response)
         gaps_row = length(normal),
         color=colorRampPalette( c("skyblue", "white", "firebrick3"))(10),
         annotation_row= subset(genes_exp_clinical_1, select = type),annotation_names_row = F#,
         #annotation_colors = ann_colors,
         #display_numbers=TRUE, number_format='%.1f'
)
save_pheatmap_pdf(p1, "./out/pheatmap_exp_1.pdf",4,5)
#jak-stat
p1 <-pheatmap(exp[c(normal,tumor),c("Jak-STAT_pathway","Core_Jak-STAT_pathway")],
              cluster_rows  = T,cluster_cols = T,
              show_colnames = T, show_rownames = F,
              border_color = "white",#cellwidth=10,cellheight = 10,#gaps_col= 38, #table(exp_marker1$response)
              gaps_row = length(normal),
              color=colorRampPalette( c("skyblue", "white", "firebrick3"))(10),
              annotation_row= subset(genes_exp_clinical_1, select = type),annotation_names_row = F#,
              #annotation_colors = ann_colors,
              #display_numbers=TRUE, number_format='%.1f'
)
save_pheatmap_pdf(p1, "./out/pheatmap_exp_3.pdf",3,5)

##cancer type
for (i in unique(genes_exp_clinical_1$cancer) ){
  print(i)
  ##
  exp<-genes_exp_clinical_1[genes_exp_clinical_1$cancer==i,]
  tumor <-exp[exp$type=="Tumor",]$sample
  normal <-exp[exp$type=="Normal",]$sample
  exp<-exp[,c(genes_1,"Jak-STAT_pathway","Core_Jak-STAT_pathway")]
  #exp<-t(scale(t(exp)))
  exp<-scale(exp)
  #rownames(exp)<-genes_exp_clinical_1$sample
  
  ##
  p1 <-pheatmap(exp[c(normal,tumor),],cluster_rows  = F,cluster_cols = T,
                show_colnames = T, show_rownames = F,
                border_color =NA,# "white",#cellwidth=10,cellheight = 10,#gaps_col= 38, #table(exp_marker1$response)
                gaps_row = length(normal),
                color=colorRampPalette( c("skyblue", "white", "firebrick3"))(10),
                annotation_row= subset(genes_exp_clinical_1, select = type),annotation_names_row = F#,
                #annotation_colors = ann_colors,
                #display_numbers=TRUE, number_format='%.1f'
  )
  save_pheatmap_pdf(p1, paste0("./out/pheatmap_exp_",i,"_1.pdf"),3,5)
  
}
#cancer mean score


#3.0 prognosis----

###3.1 COX Univariate Cox hazard analysis in all cancer---
library(survival)
data2<-genes_exp_clinical_1[genes_exp_clinical_1$type=="Tumor",] #tumor  8600
table(data2$cancer)
# ACC BLCA BRCA CESC CHOL COAD DLBC ESCA  GBM HNSC KICH KIRC KIRP  LGG LIHC LUAD LUSC 
# 79  408 1102  306   36  285   48    2  166  522   66  534  291  529  373  517  501 
# MESO   OV PAAD PCPG PRAD READ SARC SKCM STAD TGCT THCA THYM UCEC  UCS  UVM 
# 87    2  179  184  465   94  263  472    3  139  513  120  177   57   80
data2<- data2[!(data2$cancer %in% c("OV","ESCA","STAD")), ]

colnames(data2)[c(178,177)]<-c("futime","fustat")


uni_cox_fun <-function(data2,input_gene){
  outTab=data.frame()
  sigGenes=c("futime","fustat")
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

uni_cox_fun(data2,c("Jak-STAT_pathway","Core_Jak-STAT_pathway") )
ggsave("./out/dot_uniCOX_HR_pathway_all.pdf",width = 4,height = 1.5,family="serif")
uni_cox_fun(data2,genes_1 )
ggsave("./out/dot_uniCOX_HR_core_genes.pdf",width = 3.5,height = 3.5,family="serif")
uni_cox_fun(data2,c(genes_1,"Jak-STAT_pathway","Core_Jak-STAT_pathway") )
ggsave("./out/dot_uniCOX_HR_core_genes_pathway.pdf",width = 4.5,height = 3.5,family="serif")

##genes km----
#all genes in all cancers
library(survminer)
for (j in c(genes_1,input_gene) ){
  print(j)
  data2["riskscore"] <-data2[j] 
  #
  dat <-data2
  cut_value <- surv_cutpoint(dat, #数据集
                             time = "futime", #生存状态
                             event = "fustat", #生存时间
                             variables = c("riskscore") )$cutpoint[,1] #best cutoff
  dat$Risk <- ifelse(dat$riskscore> cut_value,# median(dat$riskscore),#cut_value,
                     "High","Low") # 均值？
  dat$futime<-dat$futime/365
  
  fit <- survfit(Surv(futime, fustat) ~ Risk, data =dat)
  #pdf(paste0("./out/km_",i,"_",j,".pdf"),width=3,height=3.5,family = "serif")
  print(ggsurvplot_list(fit,data = dat,pval = T,pval.method = TRUE,##是否添加P值
                        conf.int = F,### 是否添加置信区间
                        legend.title = j,#title = i,#"Risk", # 设置图例标题
                        legend.labs = c("High", "Low"), #legend = c(0.8, 0.2),# 指定图例分组标签
                        risk.table = F, # 是否添加风险表
                        risk.table.col = "strata", 
                        censor=F,size=0.5,
                        surv.scale="percent",
                        ###linetype = "strata",
                        #surv.median.line = "hv", # 是否添加中位生存线
                        risk.table.y.text.col = F,risk.table.y.text = FALSE,
                        ggtheme = theme_bw(base_family = "serif")#+theme(legend.text = element_text(colour = c("red", "blue")))
                        +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
                        +theme(plot.title = element_text(hjust = 0.5,size = 16),axis.title.y.left = element_text(size = 16,vjust = 1),axis.title.x.bottom = element_text(size = 16,vjust = 0))
                        +theme(axis.text.x.bottom = element_text(size = 12,vjust = -0.8,colour = "black")) #,face = "bold"
                        +theme(axis.text.y.left = element_text(size = 12,vjust = 0.5,hjust = 0.5,angle = 90,colour = "black"))
                        +theme(legend.title = element_text(family = "Times",colour = "black",size = 14))
                        +theme(legend.text = element_text(family = "Times",colour = "black",size =12)), # Change ggplot2 theme
                        palette = c("#bc1e5d", "#0176bd"),#"lacent",#c("#bc1e5d", "#0176bd"),##如果数据分为两组，就写两种颜色，三组就写三种颜色，具体的颜色搭配参考上一期发布的R语言颜色对比表
                        #xlim=c(10,5000),##x轴的阈值，根据随访数据的天数进行更改
                        xlab = "Years")##随访的时间时天，就写Days，月份就写Months
        
        
  )
  #dev.off()
  ggsave(paste0("./out/km_all_cancer_",j,".pdf"),width=3,height=3.5,family = "serif")
}

#all genes in single cancer 
for (j in c(genes_1,input_gene) ){
  print(j)
  data2["riskscore"] <-data2[j] 
  #
  for (i in unique(data2$cancer)){
    print(i)
    dat <-data2[data2$cancer==i,]
    cut_value <- surv_cutpoint(dat, #数据集
                               time = "futime", #生存状态
                               event = "fustat", #生存时间
                               variables = c("riskscore") )$cutpoint[,1] #best cutoff
    dat$Risk <- ifelse(dat$riskscore> cut_value,# median(dat$riskscore),#cut_value,
                       "High","Low") # 均值？
    dat$futime<-dat$futime/365
    fit <- survfit(Surv(futime, fustat) ~ Risk, data =dat)
    #pdf(paste0("./out/km_",i,"_",j,".pdf"),width=3,height=3.5,family = "serif")
    print(ggsurvplot_list(fit,data = dat,pval = T,pval.method = TRUE,##是否添加P值
                          conf.int = F,### 是否添加置信区间
                          legend.title = j,#title = i,#"Risk", # 设置图例标题
                          legend.labs = c("High", "Low"), #legend = c(0.8, 0.2),# 指定图例分组标签
                          risk.table = F, # 是否添加风险表
                          risk.table.col = "strata", 
                          censor=F,size=0.5,
                          surv.scale="percent",
                          ###linetype = "strata",
                          #surv.median.line = "hv", # 是否添加中位生存线
                          risk.table.y.text.col = F,risk.table.y.text = FALSE,
                          ggtheme = theme_bw(base_family = "serif")#+theme(legend.text = element_text(colour = c("red", "blue")))
                          +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
                          +theme(plot.title = element_text(hjust = 0.5,size = 16),axis.title.y.left = element_text(size = 16,vjust = 1),axis.title.x.bottom = element_text(size = 16,vjust = 0))
                          +theme(axis.text.x.bottom = element_text(size = 12,vjust = -0.8,colour = "black")) #,face = "bold"
                          +theme(axis.text.y.left = element_text(size = 12,vjust = 0.5,hjust = 0.5,angle = 90,colour = "black"))
                          +theme(legend.title = element_text(family = "Times",colour = "black",size = 14))
                          +theme(legend.text = element_text(family = "Times",colour = "black",size =12)), # Change ggplot2 theme
                          palette = c("#bc1e5d", "#0176bd"),#"lacent",#c("#bc1e5d", "#0176bd"),##如果数据分为两组，就写两种颜色，三组就写三种颜色，具体的颜色搭配参考上一期发布的R语言颜色对比表
                          #xlim=c(10,5000),##x轴的阈值，根据随访数据的天数进行更改
                          xlab = "Years")##随访的时间时天，就写Days，月份就写Months
          
          
    )
    #dev.off()
    ggsave(paste0("./out/km_all_cancer_",i,"_",j,".pdf"),width=3,height=3.5,family = "serif")
    
  }
}

###3.2 bacth COX Univariate Cox hazard analysis in all single cancer---
#type<- c()
#out_tab<-data.frame()#(HR="",L95="",R95="",pval="",cancer="")
library(survminer)
for (j in c(input_gene) ){
  print(j)
  data2["riskscore"] <-data2[j] 
  ##
  out_tab<-data.frame()
  for (i in  unique(data2$cancer) ) {
    print(i)
    dat <- data2[data2$cancer==i,]
    #dat <- na.omit(dat)
    #plot
    cut_value <- surv_cutpoint(dat, #数据集
                               time = "futime", #生存状态
                               event = "fustat", #生存时间
                               variables = c("riskscore") )$cutpoint[,1] #best cutoff
    dat$Risk <- ifelse(dat$riskscore> cut_value,# median(dat$riskscore),#cut_value,
                       "High","Low") # 均值？
    dat$futime<-dat$futime/365
    
    fit <- survfit(Surv(futime, fustat) ~ Risk, data =dat)
    pval <- survdiff(Surv(futime, fustat) ~ Risk, rho = 0, data =dat)
    pval <- pval$pvalue
    hr <- summary(coxph(Surv(futime, fustat) ~ Risk, data =dat))
    HR <-hr$conf.int[,"exp(coef)"]
    L95 <-hr$conf.int[,"lower .95"]
    R95 <-hr$conf.int[,"upper .95"]
    out_tab<-rbind(out_tab,
                   data.frame(HR=HR,L95=L95,R95=R95,pval=pval,cancer=i)
    )
    #plot km start
#    pdf(paste0("./out/km_",i,"_",j,".pdf"),width=3,height=3.5,family = "serif")
    print(ggsurvplot_list(fit,data = dat,pval = T,pval.method = TRUE,##是否添加P值
                          conf.int = F,### 是否添加置信区间
                          legend.title = j,#title = i,#"Risk", # 设置图例标题
                          legend.labs = c("High", "Low"), #legend = c(0.8, 0.2),# 指定图例分组标签
                          risk.table = F, # 是否添加风险表
                          risk.table.col = "strata", 
                          censor=F,size=0.5,
                          surv.scale="percent",
                          ###linetype = "strata",
                          #surv.median.line = "hv", # 是否添加中位生存线
                          risk.table.y.text.col = F,risk.table.y.text = FALSE,
                          ggtheme = theme_bw(base_family = "serif")#+theme(legend.text = element_text(colour = c("red", "blue")))
                          +theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())
                          +theme(plot.title = element_text(hjust = 0.5,size = 16),axis.title.y.left = element_text(size = 16,vjust = 1),axis.title.x.bottom = element_text(size = 16,vjust = 0))
                          +theme(axis.text.x.bottom = element_text(size = 12,vjust = -0.8,colour = "black")) #,face = "bold"
                          +theme(axis.text.y.left = element_text(size = 12,vjust = 0.5,hjust = 0.5,angle = 90,colour = "black"))
                          +theme(legend.title = element_text(family = "Times",colour = "black",size = 14))
                          +theme(legend.text = element_text(family = "Times",colour = "black",size =12)), # Change ggplot2 theme
                          palette = c("#bc1e5d", "#0176bd"),#"lacent",#c("#bc1e5d", "#0176bd"),##如果数据分为两组，就写两种颜色，三组就写三种颜色，具体的颜色搭配参考上一期发布的R语言颜色对比表
                          #xlim=c(10,5000),##x轴的阈值，根据随访数据的天数进行更改
                          xlab = "Years")##随访的时间时天，就写Days，月份就写Months
          
          
    )
#    dev.off()
    #plot km end 
  }
  if(T){
    #plot HR dot
    #plot https://github.com/yuyang3/pan-B/blob/main/Figure5_survival.R
    plot_df <-out_tab
    plot_df$R95[28]<-plot_df$L95[28]
    plot_df[1:3] <-log2(plot_df[1:3]+1)
    
    plot_df$color[plot_df$pval>= 0.05 & plot_df$HR < 1] <- "Better survival (P > 0.05)"
    plot_df$color[plot_df$pval< 0.05 & plot_df$HR < 1] <- "Better survival"
    plot_df$color[plot_df$pval >= 0.05 & plot_df$HR > 1] <- "Worse survival (P > 0.05)"
    plot_df$color[plot_df$pval < 0.05 & plot_df$HR > 1] <- "Worse survival"
    
    plot_df$color <- factor(plot_df$color, levels = c(
      "Better survival", "Better survival (P > 0.05)",
      "Worse survival (P > 0.05)", "Worse survival"
    ))
    mycolor <- c(
      "Better survival" = "#0F7B9F",
      "Better survival (P > 0.05)" = "#E0F3F8FF",
      "Worse survival (P > 0.05)" = "#FDDBC7FF",
      "Worse survival" = "#C3423F"
    )
    plot_df$significance <- ""
    plot_df$significance[plot_df$pval <= 0.05] <- "*"
    plot_df$significance[plot_df$pval <= 0.01] <- "**"
    plot_df$significance[plot_df$pval <= 0.001] <- "***"
    #plot_df$significance[plot_df$pval <= 0.0001] <- "#"
    
    #plot
    ggplot(data = plot_df, aes(x = reorder(cancer,-HR), y = HR, ymin = L95, ymax = R95, color = color)) +
      geom_pointrange() +
      geom_hline(yintercept = log2(1+1), colour = "grey40", linetype = "dashed", size = 0.2) + # add a dotted line at x=1 after flip
      geom_text(aes(x = cancer, y = max(R95), label = significance), show.legend = FALSE) +
      scale_color_manual(values = mycolor, name = "Survival association") +
      xlab("") +
      ylab("Hazard ratio (95% CI)") +
      theme_classic()+#cowplot::theme_cowplot() +
      theme(
        #text = element_text(size = 7, family = "ArialMT"),
        axis.text.x = element_text(size = 8, angle = 45, hjust = 1),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(hjust = 0.5, size = 8, face = "plain"),
        plot.margin = unit(c(1, 1, 1, 1), "char"),
        axis.line = element_line(linetype = 1, color = "black", size = 0.3),
        axis.ticks = element_line(linetype = 1, color = "black", size = 0.3)
      ) #+ggtitle(paste0(j," score"))
    ggsave(paste0("./out/HR_cancers_",j,"_score.pdf"),width = 9,height = 3,family="serif")
    #pan-cancer !
  }
  ##
}

#4.0 ROC---- 
for (j in c(genes_1,input_gene) ){
  print(j)
  data2["riskscore"] <-data2[j] 
  
  #
  for (i in unique(data2$cancer)){
    print(i)
    df <-data2[data2$cancer==i,]
  ##
    #ROC
    mk=c(1,3,5)
    ROC_rt=timeROC(T=df$futime/365,delta=df$fustat,#df$OS.time,delta=df$OS, 单位对应-年
                   marker=df$riskscore,cause=1,
                   weighting='aalen',
                   times=mk,ROC=TRUE) #slow
    p.dat=data.frame()
    for(x in which(ROC_rt$times>0)){
      lbs=paste0('AUC at ',mk[x],' years: ',round(ROC_rt$AUC[x],2))
      p.dat=rbind.data.frame(p.dat,data.frame(V1=ROC_rt$FP[,x],V2=ROC_rt$TP[,x],Type=lbs))
    }
    
    ggplot(p.dat, aes(x=V1,y=V2, color=Type))+
      #stat_smooth(aes(colour=Type),se = FALSE, size = 1)+
      geom_line(aes(x=V1,y=V2), size = 0.5)+
      scale_color_manual(values = colors)+
      scale_fill_manual(values = colors)+
      #geom_abline(slope = 1, intercept = 0, color = "grey", size = 0.5, linetype = 2)+
      theme_bw() +
      labs(x = "1-Specificity (FPR)", y = "Sensitivity (TPR)",family="serif",color="")+ #size=12,
      theme(
        legend.position=c(0.6,0.20),legend.text = element_text(size = 14),
        legend.background =element_blank(),
        axis.text = element_text(size = 14, color = "black",family="serif"),#坐标轴字体
        axis.title.x = element_text( size = 14, color = "black",family="serif"),
        axis.title.y = element_text( size = 14, color = "black",family="serif"),panel.grid = element_blank()#删除网格
      )+
      scale_x_continuous(expand = c(0,0.05)) +scale_y_continuous(expand = c(0,0.01))
    ggsave(paste0("./out/ROC_",j,"_",i,".pdf"),family="serif",width = 3.5,height = 3.2)
    
    #write.csv(p.dat,"./out/pan-cancer_ROC.csv")
    
  }
}

##all type ROC
rt<-data2
if(T){
  roc_tab<-data.frame(matrix(ncol = 4))
  colnames(roc_tab)<-c("x1","x3","x5","type")
  
  for (i in names(table(rt$cancer))) {
    print(i)
    dat <- rt[rt$cancer==i,]
    dat$futime<-dat$futime/365
    #
    mk=c(1,3,5)
    ROC_rt=timeROC(T=dat$futime,delta=dat$fustat,
                   marker=dat$riskscore,cause=1,
                   weighting='aalen',
                   times=mk,ROC=TRUE)
    
    p.dat=data.frame()
    for(x in which(ROC_rt$times>0)){
      lbs=paste0('AUC ',mk[x],' years: ',round(ROC_rt$AUC[x],2))
      p.dat=rbind.data.frame(p.dat,data.frame(V1=ROC_rt$FP[,x],V2=ROC_rt$TP[,x],Type=lbs))
    }
    roc<-data.frame(
      x1=as.numeric(round(ROC_rt$AUC[1],2)),
      x3=as.numeric(round(ROC_rt$AUC[2],2)),
      x5=as.numeric(round(ROC_rt$AUC[3],2)),
      type=i
    )
    roc_tab<-rbind(roc_tab,roc
    )
  }
  
  
  
  #all roc ggplot2
  roc_tab<-roc_tab[-1,]
  #棒棒糖图
  df <-reshape2::melt(roc_tab, id.vars ="type") #转为长数据
  df$variable<-substr(df$variable,2,2)#as.numeric()
  df %>% ggplot(aes(x=reorder(type,-value),y=value)) +
    #geom_boxplot()+
    geom_segment( aes(x=reorder(type,-value), xend=reorder(type,-value), y=0, yend=value),linewidth=0.5,colour="grey") +
    geom_point(aes(color=variable,size=value))+ #size = type, value
    geom_hline(yintercept = 0.7,color="grey70",linetype="dashed")+
    labs(x="",y="AUC",color="Year",size="AUC")+ #,size="Score"
    theme_bw(base_size = 14)+ylim(0,1)+#theme_classic(base_size = 14)
    theme(panel.grid  = element_blank(),axis.text.x =element_text(angle = 90) )+
    scale_color_manual(values = c("#d7b5d8","#df65b0","#ce1256"))+
    #scale_color_viridis_c()+
    #scale_color_distiller(palette = "Blues",direction = 1)+#scale_color_brewer("Spectral")+
    facet_grid(~variable,scales = "free")+rotate()
  #RColorBrewer::display.brewer.all()
  ggsave("./out/dot_pancancer_roc_1.3.5years.pdf",width = 8,height = 5,family="serif")
  write.csv(df,"./out/pan-cancer_ROC.csv")
  
}

#5.0 my cancer exp and km in subtypes----
#genes_exp_clinical_1
##5.1 histological_type ----
for (j in input_gene){
  print(j)
  df <-genes_exp_clinical_1
  df$Gene <-df[,j]
  df$histological_type <-gsub(";.*","",df$histological_type)
  ##
  for ( cancer in  unique(genes_exp_clinical_1$cancer)  ){ # unique(genes_exp_clinical$cancer) c("LUAD","LUSC","LIHC")
    print(cancer)
    # 分期
    if( length(unique(df[df$cancer==cancer,]$histological_type)) >1 )
    {
      ##
      print(paste0(cancer," has multiple histological_type") )
      #print(table(df[df$cancer==cancer,]$histological_type))
      #plot
      df1 <-df[df$cancer==cancer,] #remove NA
      df1<-df1[is.na(df1$histological_type) ==F,]
      df1<-df1[df1$histological_type !="",]
      print(
        ggplot(df1,
               aes(histological_type,Gene,color= histological_type ) )+
          geom_boxplot(alpha=1,size=0.15,outlier.colour = NA)+
          geom_point(aes(color =histological_type ),alpha=0.5,size=1,
                     position=position_jitterdodge(jitter.width = 0.5,
                                                   jitter.height = 0,
                                                   dodge.width = 0.8)  )+
          ggtitle(cancer)+xlab("")+ylab(j)+#ylab(paste0(genes," expression") )+#xlab(cancer)+
          EnvStats::stat_n_text()+
          scale_color_manual(values = c(RColorBrewer::brewer.pal(8,"Paired"),colors) )+  #colors
          scale_fill_manual(values =c(RColorBrewer::brewer.pal(8,"Paired"),colors) )+ 
          theme_bw()+
          theme(legend.position = "none", 
                legend.text = element_text(size = 12), 
                axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 12), 
                axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12),
                panel.grid  = element_blank()  )+rotate_x_text(angle = 90)+
          stat_compare_means(#comparisons = my_comparisons,#aes(group = tissue), 
            label ="p.format",label.y.npc = "top",# label.x.npc = "center",#label = "p.signif","p.format"
            hide.ns = T) #method = "wilcox.test",
      )
      #save 
      hight <-ifelse( nchar( max(unique(df[df$cancer==cancer,]$histological_type)) ) <30,#nchar(unique(df[df$cancer==cancer,]$histological_type)[1]) < 17,
                      4.5,#4.5,#nchar(unique(df[df$cancer==cancer,]$histological_type)[1])/4,
                      #nchar(unique(df[df$cancer==cancer,]$histological_type)[1])/4  
                      nchar( max(unique(df[df$cancer==cancer,]$histological_type)) )/4
                      )
      ggsave(paste0("./out/","histological_type","_",j,"_",cancer,".pdf"),height = hight,width = 0.8*length(unique(df[df$cancer==cancer,]$histological_type)) )
      ##
    }
    else
    {
      print(paste0(cancer," has no multiple histological_type") )
      next()  
    }
    
    
  } 
}

##5.2 pathologic_stage","clinical_stage","histological_grade分期表达----
table(df$pathologic_stage)
table(df$clinical_stage)
table(df$histological_grade)

for (j in c(genes_1,input_gene) ){
  print(j)
  df <-genes_exp_clinical_1
  df$Gene <-df[,j]
  cancer = my_cancer
  ##
  for (stage in c("pathologic_stage" ) ){ #,"clinical_stage","histological_grade"
    print(stage)
    #plot
    df1 <-df[df$cancer==cancer,] #remove NA
    #df1<-df1[is.na(df1$histological_type) ==F,]
    #df1<-df1[df1[stage] !=c("[Discrepancy]", "[Unknown]",""," "),]
    df1$Stage <-df1[,stage]
    #df1<-df1[df1$Stage != c("[Discrepancy]", "[Unknown]"," ",""),]
    df1<-df1[df1$Stage != "",];df1<-df1[df1$Stage != "[Discrepancy]",];df1<-df1[df1$Stage != "[Unknown]",]
    df1<-df1[df1$Stage != "Stage X",]
    #merge IA+IB
    df1$Stage<- gsub("Stage ","",df1$Stage ) 
    df1$Stage<-ifelse(df1$Stage %in% c("I","IA","IB"),"I",
           ifelse(df1$Stage %in% c("II","IIA","IIB"),"II",
                  ifelse(df1$Stage %in% c("III","IIIA","IIIB","IIIC"),"III","IV")
           )
    )
    #table(df1$Stage)
    df1<-df1[is.na(df1$Stage ) ==F,]
    df1<-df1[c("Stage","Gene")];df1<-na.omit(df1)

    ##
    print(
      ggplot(df1,
             aes(Stage,Gene,color= Stage ) )+
        geom_boxplot(alpha=1,size=0.15,outlier.colour = NA)+
        geom_point(aes(color =Stage ),alpha=0.5,size=1,
                   position=position_jitterdodge(jitter.width = 0.5,
                                                 jitter.height = 0,
                                                 dodge.width = 0.8)  )+
        ggtitle(cancer)+xlab("")+ylab(j )+#xlab(cancer)+
        EnvStats::stat_n_text()+
        scale_color_manual(values = c(RColorBrewer::brewer.pal(8,"Paired"),colors) )+  #colors
        scale_fill_manual(values =c(RColorBrewer::brewer.pal(8,"Paired"),colors) )+ 
        theme_bw()+
        theme(legend.position = "none", 
              legend.text = element_text(size = 12), 
              axis.title.x = element_text(size = 14), axis.text.x = element_text(size = 12), 
              axis.title.y = element_text(size = 14), axis.text.y = element_text(size = 12),
              panel.grid  = element_blank()  )+#rotate_x_text(angle = 90)+
        stat_compare_means(#comparisons = my_comparisons,#aes(group = tissue), 
          label ="p.format",label.y.npc = "top",# label.x.npc = "center",#label = "p.signif","p.format"
          hide.ns = T) #method = "wilcox.test",
    )
    #save 
    # hight <-ifelse( nchar(unique(df[df$cancer==cancer,]$histological_type)[1]) < 17,
    #                 4.5,#nchar(unique(df[df$cancer==cancer,]$histological_type)[1])/4,
    #                 nchar(unique(df[df$cancer==cancer,]$histological_type)[1])/4  )
    # width <-ifelse( length(table(df[df$cancer==cancer,][stage])) < 17,
    #                  1.5*length(unique(df[df$cancer==cancer,][stage])) ,#nchar(unique(df[df$cancer==cancer,]$histological_type)[1])/4,
    #                  length(table(df[df$cancer==cancer,][stage])) /2  )
    #width <-ifelse(stage=="clinical_stage",,)
    ggsave(paste0("./out/",stage,"_",j,"_",cancer,".pdf"),height = 3.5,width =  3.5) #length( unique(df1["Stage"])[[1]])*0.6
    
    
  } #for stage
  
}#for j

#6.0 TCGA_vs_GTEX pan-cancer----
#prepare_GTEX_TCGA_UCSC.R
##6.1 load data----
#load("~/R/my_projects/JAK_STAT_BC/1.pan-cancer/data/GTEX_TCGA_Tagrget.RData")
#load("~/R/my_projects/JAK_STAT_BC/1.pan-cancer/data/GTEX_TCGA_Tagrget_meta.RData")
#remove Target cohort
gtex_meta <-gtex_meta[!(gtex_meta$cohort=="TARGET"),]
table(gtex_meta$cohort)
boxplot(gtex.fpkm.pro[1:100,1:500])#0-5w
##remove batch effect
gtex.fpkm.pro <-gtex.fpkm.pro[,colnames(gtex.fpkm.pro) %in% gtex_meta$Sample]#12933
gtex.fpkm.pro <-log2(gtex.fpkm.pro+1)
boxplot(gtex.fpkm.pro[1:100,1:500])#0-15
##run
## 各选10各样本，箱线图看批次效应
library(dplyr)
box_draw <- gtex.fpkm.pro %>% select(grep("^GTE",colnames(gtex.fpkm.pro))[1:10],grep("^TCGA",colnames(gtex.fpkm.pro))[1:10])
boxplot(box_draw, col = "lightblue",las = 2)
#
library(sva)
batch <- substr(colnames(gtex.fpkm.pro),1,4)#c(rep("TCGA",length(grep("^TCGA",colnames(gtex.fpkm.pro)))),rep("GTEx",length(grep("^GTE",colnames(gtex.fpkm.pro))))) #批次信息
combat_data <- ComBat(gtex.fpkm.pro, batch = batch)
## 看一下去除之后的
combat_data_draw <- combat_data %>% data.frame() %>% select(grep("^GTE",colnames(combat_data))[1:10],grep("^TCGA",colnames(combat_data))[1:10])
boxplot(combat_data_draw, col = "pink",las = 2)

#6.2 score----
exp <- combat_data#gtex.fpkm.pro#
rm(gtex.fpkm.pro)
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
gsva_matrix$Sample <-rownames(gsva_matrix)

#merge in meta
genes_exp_clinical_1 <-merge(gtex_meta,gsva_matrix,by="Sample")
#rownames(genes_exp_clinical_1)<-genes_exp_clinical_1$Sample

#merge genes exp
genes_exp <-exp[rownames(exp) %in% genes,]#155
#genes[(genes %in% rownames(genes_exp) ) ==F]   # "IFNLR1" "IFNL2"  "IFNL3"  "IFNL1"
genes_exp <- log2(genes_exp+1)#log2(genes_exp+1)  #155-19109
genes_exp <- as.data.frame(t(genes_exp) )
genes_exp$Sample <-rownames(genes_exp)
genes_exp_clinical_1 <-merge(genes_exp_clinical_1,genes_exp,by="Sample")
plot(genes_exp_clinical_1$Gene)#删除离群值

genes_exp_clinical_1 <-genes_exp_clinical_1[genes_exp_clinical_1$`Jak-STAT_pathway` >0,]
#write.csv(genes_exp_clinical_1,"./out/genes_exp_clinical_1_GTEX_TCGA.csv")

#6.3 plot----
for (i in input_gene ){
  print(i)
  Gene = i
  genes_exp_clinical_1$Gene <-genes_exp_clinical_1[,i]
  ##
  ggplot(genes_exp_clinical_1,aes(x=group,y=Gene,color=group ))+ #reorder(cancer,-STAT4)
    #geom_boxplot()+
    # geom_boxplot(#width=0.55, #alpha=0.7,
    #   #position=position_dodge(width=0.8),
    #   size=0.05,outlier.colour = NA)+
    # scale_fill_manual(values = c("skyblue","pink") )+ #c( "#56B4E9","#CC79A7")
    # geom_point(aes(fill=group),alpha=0.1,size=1.5,
    #            position=position_jitterdodge(jitter.width = 0.45,#0.45  1.2
    #                                          jitter.height = 0,
    #                                          dodge.width = 0.8))+
    geom_boxplot(alpha=1,width=0.45,fill=NA,position=position_dodge(width=0.8),
                 size=0.4,outlier.shape = NA,#outlier.size = 1.5,outlier.colour = "red",
                 outlier.stroke = 0.5)+
    geom_violin(alpha=0.2,width=0.9,
                position=position_dodge(width=0.8),
                linewidth=0.25)+
    scale_fill_manual(values =  c("skyblue","pink"))+ #"#EDA065","#66CCFF"
    scale_color_manual(values =  c("skyblue","pink"))+
    #geom_text(aes(y=STAT4+0.3,label=Freq,color="grey10",size = 2) )+
    #scale_y_continuous(expand = c(0,0),limits = c(0,14))+
    labs(x=NULL,y=Gene)+#ylim(min(genes_exp_clinical_1$Gene)-0.1,max(genes_exp_clinical_1$Gene)+0.1)+
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
      label = "p.signif",label.x = 1.5,show.legend = F,size=4,
      method = "wilcox.test", hide.ns = T) #, label.y.npc = "top"
  ggsave(paste0("./out/",project,"_",i,"_pan_TvsN_TCGA-GTEX.pdf"),family="serif",width = 3.0,height = 2.8)
  
  ##
  #1.2 cancer T vs. N----
  
  ggplot(genes_exp_clinical_1,aes(x=TCGA,y=Gene,fill=group))+ #reorder(cancer,-STAT4)
    #geom_boxplot()+
    geom_boxplot(#width=0.55, #alpha=0.7,
      #position=position_dodge(width=0.8),
      size=0.05,outlier.colour = NA)+
    scale_fill_manual(values =  c("skyblue","pink"))+ #"#EDA065","#66CCFF"
    scale_color_manual(values =  c("skyblue","pink"))+#geom_text(aes(y=STAT4+0.3,label=Freq,color="grey10",size = 2) )+
    #scale_y_continuous(expand = c(0,0),limits = c(0,14))+
    labs(x=NULL,y=Gene)+#ylim(0.2,max(genes_exp_clinical_1$Gene)+0.1 )+
    #EnvStats::stat_n_text()+ #library(EnvStats)
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
    ) +rotate_x_text(angle = 90)+
    stat_compare_means(#comparisons = my_comparisons,#aes(group = tissue), 
      label = "p.signif",show.legend = F,size=4, #label.x = 1.5,
      method = "wilcox.test", hide.ns = T) #, label.y.npc = "top"
  ggsave(paste0("./out/",project,"_",i,"_cancers_TvsN_TCGA-GTEX.pdf"),family="serif",width = 12,height = 3)
  
}
#save