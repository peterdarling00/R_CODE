# 1.准备输入数据
#加载临床信息和表达矩阵
rm(list = ls())
# 表达矩阵 exprSet
# 临床信息 meta

#2.构建lasso回归模型，选择需要分析的临床指标
#输入数据是表达矩阵(仅含tumor样本)和每个病人对应的生死（顺序必须一致）
x=t(exprSet)
y=meta$event
library(glmnet)

# 2.1挑选合适的λ值
set.seed(1006)
cv_fit <- cv.glmnet(x=x, y=y)
plot(cv_fit)

fit <- glmnet(x=x, y=y)
plot(fit,xvar = "lambda")

# 2.2 用这两个λ值重新建模
model_lasso_min <- glmnet(x=x, y=y,lambda=cv_fit$lambda.min)
model_lasso_1se <- glmnet(x=x, y=y,lambda=cv_fit$lambda.1se)

# 选中的基因与系数存放于模型的子集beta中，用到的基因有一个s0值，没用的基因只记录了“.”，所以可以用下面代码挑出用到的基因。
head(model_lasso_min$beta,20)
choose_gene_min=rownames(model_lasso_min$beta)[as.numeric(model_lasso_min$beta)!=0]
choose_gene_1se=rownames(model_lasso_1se$beta)[as.numeric(model_lasso_1se$beta)!=0]
length(choose_gene_min)
length(choose_gene_1se)

save(choose_gene_min,file = paste0(proj,"_lasso_choose_gene_min.Rdata"))
save(choose_gene_1se,file = paste0(proj,"_lasso_choose_gene_1se.Rdata"))

# 3.模型预测和评估
# newx参数是预测对象。输出结果lasso.prob是一个矩阵，第一列是min的预测结果，第二列是1se的预测结果，预测结果是概率，或者说百分比，不是绝对的0和1。
# 将每个样本的生死和预测结果放在一起，直接cbind即可。
lasso.prob <- predict(cv_fit, newx=x , s=c(cv_fit$lambda.min,cv_fit$lambda.1se) )
re=cbind(y ,lasso.prob)
head(re)

re=as.data.frame(re)
colnames(re)=c('event','prob_min','prob_1se')
re$event=as.factor(re$event)

# ROC曲线-min
library(pROC)
library(ggplot2)
m <- roc(meta$event, re$prob_min)
g <- ggroc(m,legacy.axes = T,size = 1,color = "#2fa1dd")
auc(m)
#> Area under the curve: 0.9974

g + theme_minimal() +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               colour = "grey", linetype = "dashed")+
  annotate("text",x = .75, y = .25,
           label = paste("AUC of min = ",format(round(as.numeric(auc(m)),2),nsmall = 2)),color = "#2fa1dd")

# ROC曲线-1se
library(pROC)
library(ggplot2)
m <- roc(meta$event, re$prob_1se)
g <- ggroc(m,legacy.axes = T,size = 1,color = "#2fa1dd")
auc(m)
#> Area under the curve: 0.9974

g + theme_minimal() +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               colour = "grey", linetype = "dashed")+
  annotate("text",x = .75, y = .25,
           label = paste("AUC of 1se = ",format(round(as.numeric(auc(m)),2),nsmall = 2)),color = "#2fa1dd")

#两个模型的曲线画在一起
m2 <- roc(meta$event, re$prob_1se)
auc(m2)
#> Area under the curve: 0.9249
g <- ggroc(list(min = m,se = m2),legacy.axes = T,size = 1)

g + theme_minimal() +
  scale_color_manual(values = c("#2fa1dd", "#f87669"))+
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), 
               colour = "grey", linetype = "dashed")+
  annotate("text",x = .75, y = .25,
           label = paste("AUC of min = ",format(round(as.numeric(auc(m),2)),nsmall = 2)),color = "#2fa1dd")+
  annotate("text",x = .75, y = .15,
           label = paste("AUC of 1se = ",format(round(as.numeric(auc(m2)),2),nsmall = 2)),color = "#f87669")