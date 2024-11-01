library("mlr3verse")
library("precrec")
library('GGally')
HCC <- HCC[,-1]
col <- gsub('[-]', '', colnames(HCC))
colnames(HCC) <- col
HCC$condition <- as.factor(HCC$condition)
colnames(HCC) <- make.names(colnames(HCC), unique = TRUE)
task = TaskClassif$new(id="HCC", backend=HCC, target="condition") 
#tpm_out$Sample_Type <- ifelse(tpm_out$Sample_Type=="normal",0,1)
#####数据情况

autoplot(task, type = "pairs")
learners = lrns(c("classif.rpart", "classif.ranger", "classif.svm","classif.naive_bayes","classif.xgboost"),
                predict_type = "prob")

design = benchmark_grid(task, learners, rsmps("cv", folds = 5))
design

bmr = benchmark(design) # 执行基准测试
bmr$aggregate(list(msr("classif.acc"), msr("classif.auc")))

autoplot(bmr, type = "roc")# ROC 曲线
autoplot(bmr, measure = msr("classif.auc")) 
autoplot(bmr, type = "prc")
