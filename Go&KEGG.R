library(topGO)
library(clusterProfiler)
library(pathview)
library(enrichplot)
library(ReactomePA)
library(ggridges)
library(org.Hs.eg.db)
library(patchwork)
library(ggplot2)

data(geneList, package="DOSE") #富集分析的背景基因集
gene <- names(geneList)[abs(geneList) > 2]

gene.df <- bitr(turquoise_gene$nodeName, fromType = "ENSEMBL", toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
head(gene.df,2)


ggo <- groupGO(gene = gene.df$ENTREZID, OrgDb = org.Hs.eg.db, ont = "CC",level = 3,readable = TRUE)

ego_ALL <- enrichGO(gene = gene.df$ENTREZID, universe = names(geneList), OrgDb = org.Hs.eg.db,ont = "ALL",pAdjustMethod = "BH",pvalueCutoff = 1, qvalueCutoff = 1,readable = TRUE) 

ego_MF <- enrichGO(gene = gene.df$ENTREZID, universe = names(geneList),OrgDb = org.Hs.eg.db,ont = "MF", pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = FALSE)
ego_MF1 <- setReadable(ego_MF, OrgDb = org.Hs.eg.db)

ego_BP <- enrichGO(gene = gene.df$ENTREZID, universe = names(geneList),OrgDb = org.Hs.eg.db,ont = "BP", pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = FALSE)

ego_CC <- enrichGO(gene = gene.df$ENTREZID, universe = names(geneList),OrgDb = org.Hs.eg.db,ont = "CC", pAdjustMethod = "BH",pvalueCutoff = 1,qvalueCutoff = 1,readable = FALSE)

##可视化--点图
dotplot(ego_MF,title="EnrichmentGO_MF_dot")#点图，按富集的数从大到小的
dotplot(ego_CC,title="EnrichmentGO_CC_dot")
dotplot(ego_BP,title="EnrichmentGO_BP_dot")
write.csv(ego_BP,"BP_turquoise.csv")
write.csv(ego_MF,"MF_turquoise.csv")
write.csv(ego_CC,"CC_turquoise.csv")
##可视化--条形图 
barplot(ego_MF, showCategory=35,title="EnrichmentGO_MF")#条状图，按p从小到大排，绘制前20个Term


kk <- enrichKEGG(gene = gene.df$ENTREZID,organism = 'hsa', pvalueCutoff = 1,qvalueCutoff = 1)

dotplot(kk,title="Enrichment KEGG_dot",)
write.csv(kk,"kegg_turquoise.csv")

p3 = cnetplot(kk, foldChange=res_IDs$log2FoldChange,showCategory = 5, colorEdge = T, node_label="all",circular=T)
heatplot(kk, foldChange=res_IDs$log2FoldChange,showCategory = 20)

browseKEGG(kk,"hsa04061")

edox <- setReadable(kk, 'org.Hs.eg.db', 'ENTREZID')
cnetplot(edox, foldChange=res_dw$log2FoldChange, circular = F, colorEdge = TRUE)

#####通路之间 相互关联
x2 <- pairwise_termsim(kk)
emapplot(x2, showCategory = 20)

edox2 <- pairwise_termsim(kk)
p1 <- treeplot(edox2,showCategory = 20)
p2 <- treeplot(edox2, cluster.params = list(method = "average"))
aplot::plot_list(p1, p2, tag_levels='A')


#GSEA输入文件
res_IDs <- subset(res_diff,select=c("log2FoldChange","ENTREZID"))
res_diff_u <- filter(res_IDs,!duplicated(res_IDs$log2FoldChange))
data_sort <- res_diff_u %>%
  arrange(desc(log2FoldChange))
data_sort_1 <- filter(data_sort,!duplicated(data_sort$ENTREZID))
#data_sort <- unique(data_sort$log2FoldChange)
gene_list <- data_sort_1$log2FoldChange
names(gene_list) <- data_sort_1$ENTREZID
#gene_list <- unique(gene_list)
head(gene_list)

res <- gseGO(
  gene_list,    # 根据logFC排序的基因集
  ont = "ALL",    # 可选"BP"、"MF"、"CC"三大类或"ALL"
  OrgDb = org.Hs.eg.db,    # 使用人的OrgDb
  keyType = "ENTREZID",    # 基因id类型
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",    # p值校正方法
)


#gseGO.res <- gseGO(, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont="MF", minGSSize = 5, maxGSSize = 1000, pvalueCutoff=1)

ridgeplot(res, 5) 

gseaplot2(res,1:3,pvalue_table = T,ES_geom = 'line' )


gseKEGG.res <- gseKEGG(geneList = gene_list,keyType = "kegg", minGSSize = 5, maxGSSize = 1000, pvalueCutoff=0.05,verbose = F,organism = "hsa")
gseaplot2(gseKEGG.res ,1:10, pvalue_table = TRUE,ES_geom = 'line')
gseReactome.res <- gsePathway(geneList = gene_list, nPerm = 1000, minGSSize = 10, maxGSSize = 1000, pvalueCutoff=0.05)
gseaplot2(gseReactome.res,1:3, pvalue_table = T)

#####筛选上下调的通路
down_kegg<-gseKEGG.res[gseKEGG.res$pvalue<0.1 & gseKEGG.res$enrichmentScore < -0.5,];down_kegg$group=-1
up_kegg<-gseKEGG.res[gseKEGG.res$pvalue<0.1 & gseKEGG.res$enrichmentScore > 0.5,];up_kegg$group=1

geneList=res_IDs$log2FoldChange
names(geneList)=as.character(res_IDs$ENTREZID)
pro='gsea'
lapply(1:nrow(up_kegg), function(i){ 
  gseaplot2(gseKEGG.res,up_kegg$ID[i],
            title=up_kegg$Description[i],pvalue_table = T)
  ggsave(paste0(pro,'_up_kegg_',
                gsub('/','-',up_kegg$Description[i]),
                '.pdf'))
})


