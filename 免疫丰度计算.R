# 免疫细胞丰度计算
# ssGSEA方法
library(tinyarray)
load("Rdata/exp_g.Rdata")
geneset = rio::import("import/mmc3.xlsx",skip = 2)
geneset = split(geneset$Metagene,geneset$`Cell type`)
lapply(geneset[1:3], head)

library(GSVA)
f = "Rdata/immu_cell.Rdata"
if(!file.exists(f)){
  re_immune <- gsva(exp, geneset, method="ssgsea",
             mx.diff=FALSE, verbose=FALSE)
  save(re_immune,file = f) 
}
load(f)
draw_boxplot(re_immune,Group,color = c("#1d4a9d","#e5171b"))

draw_pca(re_immune,Group)








###免疫细胞和核心基因相关性图
library(ggplot2)
dat = data.frame(da)
p1 = ggscatter( dat, 
                x = "CTLA4", y = "CEP55",
                size =1,color = "#1d4a9d",
                add = "reg.line", 
                add.params = list(color = "red"))+
  stat_cor(label.y = 10 )
p2 = ggscatter( dat, 
                x = "CTLA4", y = "KIF4A",
                size = 1,color = "#1d4a9d",
                add = "reg.line", 
                add.params = list(color = "red"))+
  stat_cor(label.y = 10 )
p1+p2