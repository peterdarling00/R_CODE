### 1.3个GEO数据的下载、整理、id转换
library(tinyarray)
library(stringr)

### 1.数据的下载
## 1.1 GEO数据下载和探针名转化为基因名
# a中包含表达矩阵，临床信息，平台信息
# 注意exp的取值范围，通常在15以内，若超出范围需要将exp转化
a = geo_download("GSE10072",by_annopbrobe = F)  #注意exp的取值范围，看exp是否需要转化
group1 = ifelse(str_detect(a$pd$title,"Normal"),"normal","tumor")  #注意修改
table(group1)

b = geo_download("GSE31210",by_annopbrobe = F)
b$exp = log2(b$exp+1)  #转化exp
group2 = ifelse(str_detect(b$pd$`tissue:ch1`,"normal"),"normal","tumor")
table(group2)

x = geo_download("GSE40791",by_annopbrobe = F)
group3 = ifelse(x$pd$source_name_ch1=="normal lung","normal","tumor")
table(group3)

# 探针名id转化为基因名称symbol
find_anno("GPL96")  #注意修改
library(hgu133a.db) #需要根据平台信息进行修改
ids <- toTable(hgu133aSYMBOL)  #需要根据平台信息进行修改
a$exp = trans_array(a$exp,ids) 

find_anno("GPL570")
library(hgu133plus2.db)
ids <- toTable(hgu133plus2SYMBOL)
b$exp = trans_array(b$exp,ids)
x$exp = trans_array(x$exp,ids)

## 1.2 TCGA数据整理
# 从xena下载的fpkm数据
# TCGA下载的表达矩阵，行名为ENSG编号ensembl，需要转化为基因名称symbol
y = read.table("import/TCGA-LUAD.htseq_fpkm.tsv.gz",
               header = T,
               check.names = F,
               row.names=1)
y[1:4,1:4]
group4 = make_tcga_group(y)  #分组信息，normal在前，tumor在后的因子水平
table(group4) 
y = trans_exp(y)   # TCGA行名ENSG编号转化为基因名
y = y[rowSums(y)>0,]  #过滤标准为：在任何一个样本中表达不为0就保留
y[1:4,1:4]


### 2.数据的合并
# 不建议把不同GPL平台的数据合并，更不建议转录组数据和芯片数据合并。这里是按照原文的做法整了一下。
## 2.1 提取出共同的基因
genes = intersect_all(rownames(a$exp),
                      rownames(b$exp),
                      rownames(x$exp),
                      rownames(y))
length(genes)

## 2.2 共同的基因画韦恩图
draw_venn(x = list(GSE10072 = rownames(a$exp),
                   #GSE31210 = rownames(b$exp),
                   GSE40791 = rownames(x$exp),
                   TCGA = rownames(y)),
          "genes")
ggplot2::ggsave("plot/gene_venn.png")

# 都是些什么类型的基因呢
library(AnnoProbe)
anno = annoGene(genes,"SYMBOL")
table(anno$biotypes)

## 2.3 合并表达矩阵
exp = cbind(a$exp[genes,],
            b$exp[genes,],
            x$exp[genes,],
            y[genes,])
exp = as.matrix(exp)

## 2.4 合并分组信息
Group = c(group1,group2,group3,as.character(group4))  #先用字符串类型合并
Group = factor(Group,levels = c("normal","tumor"))    #再统一转化为因子，normal在前，tumor在后
table(Group)


### 3. 批次效应
mod = model.matrix(~as.factor(Group))

# batch里面1:n，n代表为多少个批次 
batch = rep(1:4,times = c(length(group1),length(group2),
                          length(group3),length(group4)))
table(batch)

library(sva)
f = "Rdata/exp.Rdata"
if(!file.exists(f)){
  exp_adj = ComBat(exp,batch = batch,mod = mod, par.prior=TRUE,ref.batch = 1)
  exp = limma::normalizeQuantiles(exp_adj)
  save(exp,file = f)
}
load(f)

# 样本数量太多，抽10个画个箱线图
boxplot(exp[,sample(1:ncol(exp),10)])

### 4.保存表达矩阵和分组信息
save(exp,Group,file = "Rdata/exp_group.Rdata")

### 5.差异分析和画图
fids = data.frame(probe_id = rownames(exp),
                  symbol = rownames(exp))

degs = get_deg_all(exp,Group,ids = fids)
degs$plots

ggplot2::ggsave("plot/deg_plot.png",width = 15,height =6)

### 6.保存表达矩阵、分组信息和差异分析结果
save(exp,Group,degs,file = "Rdata/exp_group_deg.Rdata")

### 7.富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(GOplot)
library(stringr)
library(tinyarray)

load("Rdata/exp_group_deg.Rdata")

# 快速富集
er = quick_enrich(degs$cgs$diff$diffgenes)
er$go.dot 
er$kk.dot

# z-score气泡图和八卦图
ego <- data.frame(er$go) 
ego <- ego[ego$ONTOLOGY=="BP",c(1,2,3,9,7)] 
ego <- ego[1:10,]
ego$geneID <- str_replace_all(ego$geneID,"/",",") 
names(ego)=c("Category","ID","Term","Genes","adj_pval")
genes = data.frame(ID=degs$deg$symbol,
                   logFC=degs$deg$logFC)
head(genes)
circ <- circle_dat(ego,genes)
GOBubble(circ, labels = 15,table.legend = F,ID = F)
GOCircle(circ)


### PPI网络
g = degs$deg$symbol[degs$deg$change!="stable"]
write.table(g,file = "export/g.txt",
            sep = "\t",quote = F,
            row.names = F,col.names = F)

# 把g.txt文件传递给STRING网页工具，导出string_interactions.tsv文件，再把这个文件传入到cytoscape,利用MCODE找到子网络，导出结果mcode_result.txt


