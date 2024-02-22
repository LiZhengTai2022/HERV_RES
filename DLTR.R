#######package_prepare######
library(DESeq2)
library(pheatmap)
library(EnhancedVolcano)
library(MASS) 
library(tidyr)
library(dplyr)
library(ggplot2)
library(RColorBrewer)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(ReactomePA)
library(ggridges)
library(org.Hs.eg.db)
library(ggpubr)
library(gmodels)
library(ggthemes)

###############文件准备与TPM计算###############
data <- read.table("all_LTR.txt",header = T,skip = 1,check.names = F)
database <- read.csv("HCC_list4_database.txt",header = T,check.names = F)
kb <- data$Length/1000
#countData <- data
#rownames(countData) <- countData[,1]
#countData <- countData[,-1]
countData <- as.matrix(data[,7:ncol(data)])
rownames(countData) <- data[,1]
rpk <- countData / kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)#计算TPM

##############使用DESeq2计算差异#######################
dds <- DESeqDataSetFromMatrix(countData, colData = database, design = ~ condition)

# 计算样本总数的一半
halfOfSamples <- ceiling(ncol(dds) / 2)

# 筛选出至少在一半以上样本中有表达的基因
dds <- dds[rowSums(counts(dds) > 0) >= halfOfSamples, ]

dds <- DESeq(dds)##标准化结果
res <- results(dds,contrast = c("condition", "T", "N"))#三个变量，第一个是需要比较的列名，第二个是分子，第三个是分母
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)),by="row.names",sort=FALSE)

#gene_df <- bitr(resdata$Row.names, fromType = "ENSEMBL", toType = c("SYMBOL","ENTREZID"), OrgDb = org.Hs.eg.db)

#resdata_symbol <- merge(resdata,gene_df,by.x="Row.names",by.y="ENSEMBL")
################绘制火山图#############
#label_mi <- resdata_symbol$SYMBOL
P <- EnhancedVolcano(resdata,lab = resdata$Row.names,x = 'log2FoldChange',y = 'padj',labSize = 4,title = "",subtitle = "",pCutoff = 0.05,FCcutoff = 1,xlim = c(-7, 7),ylim = c(0,15),border = 'full',colAlpha = 1,cutoffLineType = 'twodash',   xlab = bquote(~Log[2]~ 'fold change'),ylab = bquote(~-Log[10]~adjusted~italic(P)))
#P+ggplot2::coord_cartesian(xlim=c(-12, 12),ylim = c(0,350))+ggplot2::scale_x_continuous(breaks=seq(-12,12,2))+ggplot2::scale_y_continuous(breaks=seq(0,350,25))
############热图#########
res_de <- subset(resdata, resdata$padj<0.05, select=c('Row.names','log2FoldChange', 'pvalue','padj'))
res_de_up <- subset(res_de, res_de$log2FoldChange>=1)
res_de_dw <- subset(res_de, res_de$log2FoldChange<=(-1)*1)
# foldchang > +-1########
res_de_up_id <- as.vector(res_de_up$Row.names)
res_up <- resdata[resdata$Row.names%in%res_de_up_id,]
write.csv(res_up,"diff_up.csv")
res_de_dw_id <- as.vector(res_de_dw$Row.names)
res_dw <- resdata[resdata$Row.names%in%res_de_dw_id,]
write.csv(res_up,"diff_dw.csv")
red_de_top <- c(res_de_up_id, res_de_dw_id)
res_diff <- resdata[resdata$Row.names%in%red_de_top,]
write.csv(res_diff,"diff_gens.csv",row.names = F)

tpm_diff <- tpm[rownames(tpm)%in%res_diff$Row.names,]
tpm_up <- tpm[rownames(tpm)%in%res_up$Row.names,]
tpm_dw <- tpm[rownames(tpm)%in%res_dw$Row.names,]
#####PCA
#计算PCA
pca.info <- fast.prcomp(tpm_diff)
#显示PCA计算结果
head(pca.info$rotation)
#加入标签信息
pca.data <- data.frame(sample = rownames(pca.info$rotation),Type = database$condition,pca.info$rotation)
#绘图
ggscatter(pca.data,x = "PC1",y = "PC2",color = "Type") + theme_base()



percentVar <- round(100*attr(pca.data, "percentVar"))
pca <- ggplot(pca.data, aes(PC1, PC2, color=Type)) + 
  geom_point(size=3) +
  scale_color_manual(values=c("forestgreen","red"))+
  ggtitle("PCA") + 
  xlab(paste0("PC1: ", percentVar[1], "% variance")) + 
  ylab(paste0("PC2: ", percentVar[2], "% variance"))+
  theme_bw()
#####计算95%置信
pca+stat_ellipse(aes(color = Type), level = 0.95, show.legend = FALSE)
######
df <- database
rownames(df) <- df[,1]
df[,1] <- NULL
tpm_1 <- t(tpm)
a <- rownames(tpm_diff)
tpm_1 <- tpm_1[,colnames(tpm_1)%in%a]
tpm_out <- merge(tpm_1,df,by="row.names")
tpm_out$condition <- ifelse(tpm_out$condition=="N","0","1")

# 计算矩阵中大于一的元素比例
calc_non_zero_ratio <- function(x) {
  return(sum(x > 1 ) / length(x))
}
non_zero_ratios <- apply(tpm_out, 2, calc_non_zero_ratio)
selected_columns <- which(non_zero_ratios > 1/2)

filtered_mat <- tpm_out[, selected_columns]
b <- tpm_out$condition
filtered_mat$condition <- b
tpm_out <- filtered_mat
tpm_out_heat <- tpm_out
rownames(tpm_out_heat) <- tpm_out_heat[,1]
tpm_out_heat <- tpm_out_heat[,-1]
tpm_out_heat <- tpm_out_heat[,-ncol(tpm_out_heat)]
tpm_out_heat <- t(tpm_out_heat)
a <- rownames(tpm_out_heat)
b <- res_up$SYMBOL
c <- res_dw$SYMBOL
filter_up <- intersect(a,b)
filter_dw <- intersect(a,c)
res_up_path <- res_up[res_up$SYMBOL%in%filter_up,]
res_up_path <- subset(res_up_path,select=c("SYMBOL","ENTREZID","log2FoldChange"))
res_dw_path <- res_dw[res_dw$SYMBOL%in%filter_dw,]
res_dw_path <- subset(res_dw_path,select=c("SYMBOL","ENTREZID","log2FoldChange"))

col <- colorRampPalette(c("white", "#4CAC8E"))(70)
ann_colors <- list(condition= c(N="cyan2",T="darkmagenta"),gender=c(Male="cyan2",Female="darkmagenta"),age=col)#配色分组
coul <- colorRampPalette(brewer.pal(11, "PiYG"))(25)
pheatmap(tpm_diff, cluster_row=T, scale="row",annotation_col=df,annotation_colors = ann_colors,show_rownames = F,show_colnames = F,border_color = NA,cluster_cols = F)

library(mlr3verse)
library(data.table)
tpm_out <- tpm_out[,-1]
col <- gsub('[-]', '', colnames(tpm_out))
colnames(tpm_out) <- col
tpm_out$condition <- as.factor(tpm_out$condition)
colnames(tpm_out) <- make.names(colnames(tpm_out), unique = TRUE)
task = TaskClassif$new(id="HCC", backend=tpm_out, target="condition") 
learner = lrn("classif.ranger", importance = "impurity")#创建时"激活"其变量重要性度量
filter = flt("importance", learner = learner)

set.seed(1234)
max_importance_s <- -Inf
max_importance <- data.table()
for (i in 1:1000) {
  filter$calculate(task)
  # 计算重要性，将其存储在 'importance' 数据表中
  importance <- as.data.table(filter)
  importance_s <- importance$score[1]
  if (importance_s > max_importance_s) {
    max_importance_s <- importance_s
    max_importance <- copy(importance) # 使用 copy() 函数创建一个新的数据表，防止后续修改影响已存储的数据表
  }
}
print(max_importance)
ggdotchart(max_importance[1:50,],x="feature",y="score",color="score",sorting = "descending",add = "segments",ggtheme=theme_pubr(),rotate=T)
write.csv(max_importance,"max_importance.csv")


##########箱线图############
tpm_diff1 <- t(tpm_diff)
#行是样本，列是基因
mydata_1 <- merge(tpm_diff1,database,by.x='row.names',by.y='Run')
#mydata_1 <- mydata_1[,-(ncol(mydata_1))]
rownames(mydata_1) <- mydata_1[,1]
mydata_1 <- mydata_1[,-1]
a <- max_importance$feature[1:15]
a <- gsub('[.]', '|', a)
b <- mydata_1$condition
mydata_1 <- mydata_1[,colnames(mydata_1)[1:ncol(mydata_1)-1]%in%a]
med <- apply(mydata_1,2,median)
med <- as.numeric(med)
mydata_2 <- t(mydata_1)
mydata_2 <- cbind(mydata_2,med)
mydata_2 <- as.data.frame(mydata_2)
mydata_2 <- arrange(mydata_2,desc(med))
mydata_2$med <- NULL
mydata_1 <- t(mydata_2)
mydata_1 <- as.data.frame(mydata_1)
mydata_1$condition <- b
mydata_1
dat2 <- mydata_1
dat3 <- mydata_1[,1:ncol(mydata_1)-1]+1
dat3$condition <- b
c <- colnames(dat3)[1:7]
#dat3 <- merge(dat3,database,by.x='row.names',by.y = 'Run')
#dat3 <- dat3[,-(ncol(dat3))]

#dat4 <-  gather(dat2,key = "gene",value = "expression",-condition)
#ggplot(data = dat4,aes(x = condition,y = expression,color = condition))+geom_boxplot(aes(x = condition,y = expression,color = condition))+theme_bw(base_size = 10)+scale_color_manual(values=c("forestgreen","red"))+facet_wrap(~gene,nrow = 1)+stat_compare_means(method = "wilcox.test",hide.ns = T,label = "p.signif")+geom_jitter(aes(fill=condition),width =0.2,shape = 20,size=2.5)+theme(panel.grid=element_blank())+theme_classic()+scale_y_log10()+facet_wrap(~gene,nrow = 2)
dat4 <-  gather(dat3,key = "gene",value = "expression",-condition)
dat4$gene=factor(dat4$gene,ordered = TRUE,levels = colnames(mydata_1)[1:ncol(mydata_1)-1])
mutate(dat4,gene = factor(gene, levels = colnames(mydata_1)[1:ncol(mydata_1)-1]))
ggplot(data = dat4,aes(x = condition,y = expression,color = condition))+geom_boxplot(aes(x = condition,y = expression,color = condition))+theme_bw(base_size = 10)+scale_color_manual(values=c("forestgreen","red"))+facet_wrap(~gene,nrow = 1)+scale_y_log10()+stat_compare_means(method = "wilcox.test",hide.ns = T,label = "p.signif")+geom_jitter(aes(fill=condition),width =0.2,shape = 20,size=2.5)+theme(plot.subtitle=element_text(size=15),panel.grid=element_blank())+theme_classic()+facet_wrap(~gene,nrow = 5)

print(c)
#####find drugs
library(httr)
library(jsonlite)
library(tidyverse)


get_dgidb_data <- function(gene) {
  url <- paste0("http://www.dgidb.org/api/v2/interactions.json?genes=", gene)
  response <- GET(url)
  data <- fromJSON(content(response, "text"))
  return(data)
}

#######输入目标向量
for (i in c) {
  tryCatch({
    dgidb_data <- get_dgidb_data(i)
    drug_list_i <- dgidb_data$matchedTerms$interactions[[1]]
    drug_listi <- subset(drug_list_i, select = c("drugName", "drugConceptId", "score"))
    out_put <- paste("drug_list_", i, ".csv", sep = "")
    write.csv(drug_listi, out_put)
  }, error = function(e) {
    cat("Error occurred while processing", i, ":\n", conditionMessage(e), "\n")
  })
}



####KEGG######
gene.df <- bitr(res_dw, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
#head(gene.df,2)
options(clusterProfiler.download.method = "wininet")
kk_up <- enrichKEGG(gene = res_up$ENTREZID,organism = 'hsa',pAdjustMethod = "BH",pvalueCutoff = 1)
kk_dw <- enrichKEGG(gene = res_dw$ENTREZID,organism = 'hsa',pAdjustMethod = "BH",pvalueCutoff = 1)
barplot(kk_dw,title="Enrichment KEGG_Dotplot_dw",showCategory=20)
barplot(kk_up,title="Enrichment KEGG_Dotplot_up",showCategory=20)
cnetplot(kk_dw, foldChange=res_dw$log2FoldChange,showCategory = 5, colorEdge = T, node_label="all",circular=T)
cnetplot(kk_dw, foldChange=res_dw$log2FoldChange,showCategory = 5, colorEdge = T, node_label="all",circular=T)
heatplot(kk_dw, foldChange=res_dw$log2FoldChange,showCategory = 20)







