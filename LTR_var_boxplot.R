data <- read.table("all_LTR_PRJEB50623.txt",header = T,skip = 1,check.names = F)
colnames(data) <- gsub(".sorted.bam", "", colnames(data))
kb <- data$Length/1000
countData <- as.matrix(data[,7:ncol(data)])
rownames(countData) <- data[,1]
rpk <- countData / kb
tpm <- t(t(rpk)/colSums(rpk) * 1000000)#计算TPM
database <- read.csv("PRJEB50623.txt",header = T,check.names = F)
database <- subset(database,select=c("Run","sampling_site"))
database$condition <- ifelse(database$sampling_site=="normal tissue adjacent to neoplasm","N","T")
database$sampling_site <- NULL
database
tpm_LTR <- tpm[rownames(tpm)%in%LTR,]
tpm_LTR <- t(tpm_LTR)
mydata_1 <- merge(tpm_LTR,database,by.x="row.names",by.y="Run")
rownames(mydata_1) <- mydata_1[,1]
mydata_1 <- mydata_1[,-1]
med <- apply(mydata_1,2,median)
med <- as.numeric(med)
mydata_2 <- t(mydata_1)
mydata_2 <- cbind(mydata_2,med)
mydata_2 <- as.data.frame(mydata_2)
mydata_2 <- arrange(mydata_2,desc(med))
mydata_2$med <- NULL
mydata_1 <- t(mydata_2)
mydata_1 <- as.data.frame(mydata_1)
b <- mydata_1$condition
dat3 <- mydata_1[,1:ncol(mydata_1)-1]
dat3 <- dat3 %>%
  mutate(across(everything(), as.numeric))
dat3$condition <- b
#dat3 <- dat3+1
dat4 <-  gather(dat3,key = "gene",value = "expression",-condition)
ggplot(data = dat4, aes(x = condition, y = expression, color = condition)) +
  geom_boxplot(aes(x = condition, y = expression, color = condition)) +
  theme_bw(base_size = 10) +
  scale_color_manual(values = c( "#D2B48C","#EE82EE")) +  # 使用新的十六进制颜色代码
  facet_wrap(~gene, nrow = 1) +
  scale_y_log10() +
  stat_compare_means(method = "wilcox.test", hide.ns = TRUE, label = "p.signif") +
  geom_jitter(aes(fill = condition), width = 0.2, shape = 20, size = 2.5) +
  theme(plot.subtitle = element_text(size = 15), panel.grid = element_blank()) +
  theme_classic() +
  facet_wrap(~gene, nrow = 1)

write.csv(mydata_1,"LTR_var_BC.csv")
