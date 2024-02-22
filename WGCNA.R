library(WGCNA)
library(dplyr)
library(tidyr)
#####列为基因，行为样本
mydata <- read.csv("LGG-FPKM.csv",check.names = F,header = T)

mydata0 <- mydata %>%
  group_by(GeneSymbol) %>%
  mutate(count = row_number()) %>%
  ungroup()

mydata0 <- mydata0 %>%
  unite("Gene_name", GeneSymbol, count, sep = "", remove = TRUE)

mydata0 <- as.data.frame(mydata0)
rownames(mydata0) <- mydata0[,1]
mydata0$Gene_name <- NULL
print(mydata0)
mydata0 <- t(mydata0)

gsg = goodSamplesGenes(mydata0, verbose = 3)
gsg$allOK

if (!gsg$allOK)
{
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0)
    printFlush(paste("Removing genes:", paste(names(mydata0)[!gsg$goodGenes], collapse = ", ")))
  if (sum(!gsg$goodSamples)>0)
    printFlush(paste("Removing samples:", paste(rownames(mydata0)[!gsg$goodSamples], collapse = ", ")))
  # Remove the offending genes and samples from the data:
  mydata0 = mydata0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(mydata0), method = "average")
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

abline(h = 150000, col = "red")#根据离群点位置划线
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 4000000, minSize = 2)
#####
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==150000)
mydata = mydata0[keepSamples, ]
nGenes = ncol(mydata0)
nSamples = nrow(mydata0)

# 假设 HCC_LTR 已经是一个DataFrame，并且行名（row.names）设置为样本ID
# 首先转换每一列为因子，以确保它们可以被用于模型矩阵

HCC_LTR$LTR88c.448.551 <- factor(HCC_LTR$LTR88c.448.551)
HCC_LTR$LTR13.38.1005 <- factor(HCC_LTR$LTR13.38.1005)
HCC_LTR$LTR37int.2304.2733 <- factor(HCC_LTR$LTR37int.2304.2733)

# 为每个条件创建一个模型矩阵
datTraits_LTR88c <- model.matrix(~0 + HCC_LTR$LTR88c.448.551)
colnames(datTraits_LTR88c) <- levels(HCC_LTR$LTR88c.448.551)

datTraits_LTR13 <- model.matrix(~0 + HCC_LTR$LTR13.38.1005)
colnames(datTraits_LTR13) <- levels(HCC_LTR$LTR13.38.1005)

datTraits_LTR37 <- model.matrix(~0 + HCC_LTR$LTR37int.2304.2733)
colnames(datTraits_LTR37) <- levels(HCC_LTR$LTR37int.2304.2733)

# 合并模型矩阵
datTraits <- cbind(datTraits_LTR88c, datTraits_LTR13, datTraits_LTR37)

# 确保行名与HCC_LTR的行名匹配
rownames(datTraits) <- rownames(HCC_LTR)


sampleTree2 = hclust(dist(mydata0), method = "average")
traitColors = numbers2colors(datTraits, signed = FALSE) #用颜色代表关联度
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits),
                    main = "Sample dendrogram and trait heatmap")

c <- c("SRR10822598","SRR10822554")

b <- merge(mydata0,HCC_LTR,by = "row.names")
b <- b[,c(1,ncol(b))]
rownames(b) <- b[,1]
b$Row.names <- NULL
datTraits <- model.matrix(~0+factor(b$condition))
colnames(datTraits)=levels(factor(b$condition))
rownames(datTraits)=rownames(HCC_LTR)

options(stringsAsFactors = FALSE) 
enableWGCNAThreads() ## 打开多线程

powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(mydata0, powerVector = powers, verbose = 5)
#设置网络构建参数选择范围，计算无尺度分布拓扑矩阵
# Plot the results:
##sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")


#####一步法建立网络#####
cor <- WGCNA::cor

net = blockwiseModules(
  mydata0,
  power = sft$powerEstimate,
  maxBlockSize = nGenes,
  TOMType = "unsigned", minModuleSize = 30,
  reassignThreshold = 0, mergeCutHeight = 0.25,
  numericLabels = TRUE, pamRespectsDendro = FALSE,
  saveTOMs = TRUE,
  saveTOMFileBase = "AS-green-TPM-TOM",
  verbose = 3
)

table(net$colors)

cor<-stats::cor

#####模块可视化
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
table(mergedColors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
## assign all of the gene to their corresponding module 
## hclust for the genes.

######性状与基因模块之间的关系

design=datTraits
moduleColors <- labels2colors(net$colors)
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(mydata0, moduleColors)$eigengenes
MEs = orderMEs(MEs0); ##不同颜色的模块的ME值矩阵(样本vs模块)
moduleTraitCor = cor(MEs, design , use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
# Display the correlation values within a heatmap plot
par(mar = c(5, 10, 3, 2) + 0.1, cex.lab = 1)
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.3,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


#####计算模块与基因的相关性矩阵
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(mydata0, MEs, use = "p"));
## 算出每个模块跟基因的皮尔森相关系数矩阵
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

#####计算性状与基因的相关性矩阵
HP = as.data.frame(design[,1]);
names(LP) = "HP"
geneTraitSignificance = as.data.frame(cor(mydata0, HP, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(HP), sep="");
names(GSPvalue) = paste("p.GS.", names(HP), sep="");

#####把两个相关性矩阵联合起来并指定感兴趣模块进行分析
module = "yellow"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for N24h",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2,abline = TRUE,abline.lty = 2,col = module)

#####针对所有基因绘制分模块热图
nGenes = ncol(mydata0)
nSamples = nrow(mydata0)
geneTree = net$dendrograms[[1]]; 
dissTOM = 1-TOMsimilarityFromExpr(mydata0, power = 3); 
plotTOM = dissTOM^7; 
diag(plotTOM) = NA; 
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

#####模块与性状的关系
# Recalculate module eigengenes
MEs = moduleEigengenes(mydata0, moduleColors)$eigengenes
## 只有连续型性状才能只有计算
## 这里把是否属于 LTR_H 表型这个变量用0,1进行数值化。
HP = as.data.frame(design[,1]);
names(HP) = "HP"
# Add the weight to existing module eigengenes
MET = orderMEs(cbind(MEs, HP))
# Plot the relationships among the eigengenes and the trait
par(cex = 0.9)
plotEigengeneNetworks(MET, "", marDendro = c(0,4,1,2), marHeatmap = c(3,4,1,2), cex.lab = 0.8, xLabelsAngle
                      = 90)
# Plot the dendrogram;
par(cex = 1.0)
## 模块的聚类图
plotEigengeneNetworks(MET, "Eigengene dendrogram", marDendro = c(0,4,2,0),
                      plotHeatmaps = FALSE)
# Plot the heatmap matrix (note: this plot will overwrite the dendrogram plot)
par(cex = 1.0)
## 性状与模块热图
plotEigengeneNetworks(MET, "Eigengene adjacency heatmap", marHeatmap = c(3,4,2,2),
                      plotDendrograms = FALSE, xLabelsAngle = 90)
#####模块包含的基因
# Select module
module = "yellow";
# Select module probes
probes = colnames(mydata0) 
inModule = (moduleColors==module);
modProbes = probes[inModule];


# Recalculate topological overlap
TOM = TOMsimilarityFromExpr(mydata0, power = 6,nThreads = 4); 
# Select module
module = "yellow";
# Select module probes
probes = colnames(mydata0) ## probe就是基因名
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
## 也是提取指定模块的基因名
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
## 模块对应的基因关系矩阵 

#####导出到cytoscape
cyt = exportNetworkToCytoscape(
  modTOM,
  edgeFile = paste("CytoscapeInput-edges-yellow", paste(module, collapse="-"), ".txt", sep=""),
  nodeFile = paste("CytoscapeInput-nodes-yelow", paste(module, collapse="-"), ".txt", sep=""),
  weighted = TRUE,
  threshold = 0.02,
  nodeNames = modProbes, 
  nodeAttr = moduleColors[inModule]
);
#####网络筛选前30个主要基因
nTop = 30;
IMConn = softConnectivity(mydata0[, modProbes]);
top = (rank(-IMConn) <= nTop)
filter <- modTOM[top, top]
#####获取HUbgenes
HubGenes <- chooseTopHubInEachModule(mydata0,moduleColors)
