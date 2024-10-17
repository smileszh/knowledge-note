# scRNA-Seq QC

## 目录

- 过滤原因
- 查看数据
- 过滤前可视化
- 识别高表达基因
- 过滤
- 过滤后可视化
- 批量QC脚本







## 过滤原因

单细胞RNA测序（scRNA-seq）数据的质控（Quality Control, QC）是数据分析的关键步骤之一，旨在确保数据的高质量和可靠性。



### 1. 细胞过滤

#### 1.1 低质量细胞

- **过滤标准**：总UMI（Unique Molecular Identifier）数低于某个阈值的细胞。
- **原因**：低UMI数可能表示细胞捕获效率低或细胞死亡，导致不可靠的数据。

#### 1.2 双细胞
- **过滤标准**：检测到的基因数或UMI数显著高于平均水平的细胞。
- **原因**：双细胞（doublets）是两个或多个细胞被误认为一个细胞，可能影响数据的准确性。

#### 1.3 高线粒体基因比例

- **过滤标准**：线粒体基因表达占总表达的比例高于某个阈值的细胞（如20%）。
- **原因**：高线粒体基因比例通常表示细胞处于凋亡或应激状态。

#### 1.4 高核糖体基因比例

- **过滤标准**：核糖体基因表达占总表达的比例高于某个阈值的细胞。
- **原因**：高核糖体基因比例可能表示技术噪音或细胞状态异常。

### 2. 基因过滤

#### 2.1 低表达基因

- **过滤标准**：在很少细胞中表达的基因（如低于3个细胞）。
- **原因**：低表达基因可能是技术噪音，对下游分析贡献有限。


### 3. 批次效应

- **方法**：使用统计方法（如Harmony或MNN）进行批次效应校正。
- **原因**：不同实验批次间的技术差异可能导致数据的系统性偏差。







**当然了，不同数据有不同的特性，也会有不同的过滤方法和参数，只需要先掌握基本流程，再根据自己看到的方法调整即可。我们可以先进行检查，看有没有怪异的样本，再进行调整！**



## 计算线粒体基因比例

### 1. 识别线粒体基因

首先，通过匹配基因名中的前缀来识别线粒体基因。在大多数情况下，线粒体基因的名称以 "MT-" 开头。

```
mito_genes <- rownames(input_sce)[grep("^MT-", rownames(input_sce), ignore.case = TRUE)]
print(mito_genes)  # 这将打印出识别出的线粒体基因，可能是13个基因
```

### 2. 计算线粒体基因百分比

使用 *`PercentageFeatureSet`* 函数来计算每个细胞中线粒体基因的表达百分比，并将结果存储在 *`meta.data`* 的 *`percent_mito`* 列中。

```
input_sce <- PercentageFeatureSet(input_sce, features = mito_genes, col.name = "percent_mito")
```

### 3. 计算和打印五数概括

使用 *`fivenum`* 函数来计算 *`percent_mito`* 列的五数概括，并打印结果。

```
fivenum(input_sce@meta.data$percent_mito)
```



## 计算核糖体基因比例

### 1. 识别核糖体基因

通过匹配基因名中的前缀来识别核糖体基因。通常，核糖体基因的名称以 "Rps" 或 "Rpl" 开头。

```
ribo_genes <- rownames(input_sce)[grep("^Rp[sl]", rownames(input_sce), ignore.case = TRUE)]
print(ribo_genes)  # 这将打印出识别出的核糖体基因

Copyr
```

### 2. 计算核糖体基因百分比

使用 *`PercentageFeatureSet`* 函数来计算每个细胞中核糖体基因的表达百分比，并将结果存储在 *`meta.data`* 的 *`percent_ribo`* 列中。

```
input_sce <- PercentageFeatureSet(input_sce, features = ribo_genes, col.name = "percent_ribo")

Copyr
```

### 3. 计算和打印五数概括

使用 *`fivenum`* 函数来计算 *`percent_ribo`* 列的五数概括，并打印结果。

```
percent_ribo_fivenum <- fivenum(input_sce@meta.data$percent_ribo)
print(percent_ribo_fivenum)
```

## 计算红血细胞基因比例

### 1. 识别血红蛋白基因

通过匹配基因名中的前缀来识别血红蛋白基因。通常，血红蛋白基因的名称以 "Hb" 开头。

```
Hb_genes <- rownames(input_sce)[grep("^Hb[^(p)]", rownames(input_sce), ignore.case = TRUE)]
print(Hb_genes)  # 打印识别出的血红蛋白基因

Copyr
```

### 2. 计算血红蛋白基因百分比

使用 *`PercentageFeatureSet`* 函数来计算每个细胞中血红蛋白基因的表达百分比，并将结果存储在 *`meta.data`* 的 *`percent_hb`* 列中。

```
input_sce <- PercentageFeatureSet(input_sce, features = Hb_genes, col.name = "percent_hb")

Copyr
```

### 3. 计算和打印五数概括

使用 *`fivenum`* 函数来计算 *`percent_hb`* 列的五数概括，并打印结果。

```
percent_hb_fivenum <- fivenum(input_sce@meta.data$percent_hb)
print(percent_hb_fivenum)
```



## 过滤前可视化

```R
# 小提琴图
feats <- c("nFeature_RNA", "nCount_RNA")
p1=VlnPlot(input_sce, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2)   +NoLegend()
p1 
w=length(unique(input_sce$orig.ident))/3+5;w
ggsave(filename="Vlnplot1.pdf",plot=p1,width = w,height = 5)
  

feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2=VlnPlot(input_sce, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims=T) + 
  scale_y_continuous(breaks=seq(0, 100, 5)) +
  NoLegend()
p2	
w=length(unique(input_sce$orig.ident))/2+5;w
ggsave(filename="Vlnplot2.pdf",plot=p2,width = w,height = 5)

# 点图
feats <- c("nFeature_RNA", "nCount_RNA")
p3=FeatureScatter(input_sce, group.by = "orig.ident", features = feats, pt.size = 0.5)
p3
w=length(unique(input_sce$orig.ident))/3+5;w
ggsave(filename="Scatterplot.pdf", plot=p3, width = w,height = 5)
```

## 识别高表达基因

我们也可以画出表达量最高的50个基因，这在单细胞RNA分析中有多种用途和意义。

  1. **识别主要表达基因**

绘制表达量最高的基因的盒图可以帮助识别在数据集中主要表达的基因。这些基因可能在生物学过程中起重要作用，或者在不同细胞类型之间有显著的差异表达。

  2. **数据质量控制**

高表达基因的分布情况可以用于评估数据质量。例如，如果某些基因在所有细胞中都表现出异常高的表达，可能提示数据中存在技术性偏差或污染。

  3. **初步数据探索**

在数据分析的早期阶段，绘制高表达基因的盒图可以帮助研究人员快速了解数据的整体特征。这些信息可以指导后续的分析步骤，如聚类分析、差异表达分析等。

  4. **细胞类型标记**

高表达基因往往包括一些众所周知的细胞类型标记基因。通过观察这些基因的表达分布，可以初步推断数据集中可能存在的细胞类型。

```
C=subset(input_sce.filt,downsample=100)@assays$RNA$counts
dim(C)
C=Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
  
most_expressed <- order(apply(C, 1, median), decreasing = T)[50:1]
  
pdf("TOP50_most_expressed_gene.pdf",width = 15,height = 15)
  
par(cex = 1)
par(mar = c(4,6,4,4))
  
boxplot(as.matrix(Matrix::t(C[most_expressed, ])),
          cex = 0.2, las = 2, 
          ylim=c(0,8),
          xlab = "% total count per cell",
          col = (scales::hue_pal())(50)[50:1], 
          horizontal = TRUE)
  
  
dev.off()
```





## 过滤

### 1. 过滤最少表达基因数的细胞&最少表达细胞数的基因

也就是前面提到的低质量细胞和第表达基因。

在CreateSeuratObject的时候已经是进行了这个过滤操作。所以我们一般不再进行重复过滤。

```
  seurat_object <- CreateSeuratObject(
     counts = ct,
     project = pro,
     min.cells = 5,
     min.features = 300
   )
```

代码创建一个 Seurat 对象，其中：

- *`counts = ct`*: 使用你提供的计数矩阵 *`ct`*。
- *`project = pro`*: 使用你提供的项目名称 *`pro`*。
- *`min.cells = 5`*: 每个基因至少在 5 个细胞中表达才能被保留。
- *`min.features = 300`*: 每个细胞至少表达 300 个基因才能被保留。



但也可以自己调整这个过滤参数再次进行过滤。

```
selected_c <- WhichCells(input_sce, expression = nFeature_RNA > 500)
selected_f <- rownames(input_sce)[Matrix::rowSums(input_sce@assays$RNA$counts > 0 ) > 3]
input_sce.filt <- subset(input_sce, features = selected_f, cells = selected_c)
dim(input_sce) 
dim(input_sce.filt) 
```



### 2. 线粒体/核糖体/红血细胞基因比例(根据上面的violin图)

```
selected_mito <- WhichCells(input_sce.filt, expression = percent_mito < 25)
selected_ribo <- WhichCells(input_sce.filt, expression = percent_ribo > 3)
selected_hb <- WhichCells(input_sce.filt, expression = percent_hb < 1 )
length(selected_hb)
length(selected_ribo)
length(selected_mito)
  
input_sce.filt <- subset(input_sce.filt, cells = selected_mito)
input_sce.filt <- subset(input_sce.filt, cells = selected_ribo)
# input_sce.filt <- subset(input_sce.filt, cells = selected_hb)
dim(input_sce.filt)
  
table(input_sce.filt$orig.ident) 
```

## 过滤后可视化

```
feats <- c("nFeature_RNA", "nCount_RNA")
p1_filtered=VlnPlot(input_sce.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
  NoLegend()
p1
w=length(unique(input_sce.filt$orig.ident))/3+5
ggsave(filename="Vlnplot1_filtered.pdf",plot=p1_filtered,width = w,height = 5)


feats <- c("percent_mito", "percent_ribo", "percent_hb")
p2_filtered=VlnPlot(input_sce.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
  NoLegend()
p2
w=length(unique(input_sce.filt$orig.ident))/2+5
ggsave(filename="Vlnplot2_filtered.pdf",plot=p2_filtered,width = w,height = 5) 
```



## 批量QC脚本

**曾老师**已经写好了一个QC脚本进行单细胞进行指控，这个笔记也是源于此脚本。

```
basic_qc <- function(input_sce){
  #计算线粒体基因比例
  mito_genes=rownames(input_sce)[grep("^MT-", rownames(input_sce),ignore.case = T)] 
  print(mito_genes) #可能是13个线粒体基因
  #input_sce=PercentageFeatureSet(input_sce, "^MT-", col.name = "percent_mito")
  input_sce=PercentageFeatureSet(input_sce, features = mito_genes, col.name = "percent_mito")
  fivenum(input_sce@meta.data$percent_mito)
  
  #计算核糖体基因比例
  ribo_genes=rownames(input_sce)[grep("^Rp[sl]", rownames(input_sce),ignore.case = T)]
  print(ribo_genes)
  input_sce=PercentageFeatureSet(input_sce,  features = ribo_genes, col.name = "percent_ribo")
  fivenum(input_sce@meta.data$percent_ribo)
  
  #计算红血细胞基因比例
  Hb_genes=rownames(input_sce)[grep("^Hb[^(p)]", rownames(input_sce),ignore.case = T)]
  print(Hb_genes)
  input_sce=PercentageFeatureSet(input_sce,  features = Hb_genes,col.name = "percent_hb")
  fivenum(input_sce@meta.data$percent_hb)
  
  #可视化细胞的上述比例情况
  # feats <- c("nFeature_RNA", "nCount_RNA", "percent_mito",
  #            "percent_ribo", "percent_hb")
  feats <- c("nFeature_RNA", "nCount_RNA")
  p1=VlnPlot(input_sce, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
    NoLegend()
  p1 
  w=length(unique(input_sce$orig.ident))/3+5;w
  ggsave(filename="Vlnplot1.pdf",plot=p1,width = w,height = 5)
  
  feats <- c("percent_mito", "percent_ribo", "percent_hb")
  p2=VlnPlot(input_sce, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3, same.y.lims=T) + 
    scale_y_continuous(breaks=seq(0, 100, 5)) +
    NoLegend()
  p2	
  w=length(unique(input_sce$orig.ident))/2+5;w
  ggsave(filename="Vlnplot2.pdf",plot=p2,width = w,height = 5)
  
  
  p3=FeatureScatter(input_sce, "nCount_RNA", "nFeature_RNA", group.by = "orig.ident", pt.size = 0.5)
  p3
  w=length(unique(input_sce$orig.ident))/3+5;w
  ggsave(filename="Scatterplot.pdf", plot=p3, width = w,height = 5)
  
  #根据上述指标，过滤低质量细胞/基因
  #过滤指标1:最少表达基因数的细胞&最少表达细胞数的基因
  # 一般来说，在CreateSeuratObject的时候已经是进行了这个过滤操作
  # 如果后期看到了自己的单细胞降维聚类分群结果很诡异，就可以回过头来看质量控制环节
  # 先走默认流程即可
  if(F){
    selected_c <- WhichCells(input_sce, expression = nFeature_RNA > 500)
    selected_f <- rownames(input_sce)[Matrix::rowSums(input_sce@assays$RNA$counts > 0 ) > 3]
    input_sce.filt <- subset(input_sce, features = selected_f, cells = selected_c)
    dim(input_sce) 
    dim(input_sce.filt) 
  }
  
  input_sce.filt =  input_sce
  
  # par(mar = c(4, 8, 2, 1))
  # 这里的C 这个矩阵，有一点大，可以考虑随抽样 
  C=subset(input_sce.filt,downsample=100)@assays$RNA$counts
  dim(C)
  C=Matrix::t(Matrix::t(C)/Matrix::colSums(C)) * 100
  
  most_expressed <- order(apply(C, 1, median), decreasing = T)[50:1]
  
  pdf("TOP50_most_expressed_gene.pdf",width = 15,height = 15)
  
  par(cex = 1)
  par(mar = c(4,6,4,4))
  
  boxplot(as.matrix(Matrix::t(C[most_expressed, ])),
          cex = 0.2, las = 2, 
          ylim=c(0,8),
          xlab = "% total count per cell",
          col = (scales::hue_pal())(50)[50:1], 
          horizontal = TRUE)
  
  
  dev.off()
  
  rm(C)
  
  #过滤指标2:线粒体/核糖体基因比例(根据上面的violin图)
  selected_mito <- WhichCells(input_sce.filt, expression = percent_mito < 25)
  selected_ribo <- WhichCells(input_sce.filt, expression = percent_ribo > 3)
  selected_hb <- WhichCells(input_sce.filt, expression = percent_hb < 1 )
  length(selected_hb)
  length(selected_ribo)
  length(selected_mito)
  
  input_sce.filt <- subset(input_sce.filt, cells = selected_mito)
  input_sce.filt <- subset(input_sce.filt, cells = selected_ribo)
  input_sce.filt <- subset(input_sce.filt, cells = selected_hb)
  dim(input_sce.filt)
  
  table(input_sce.filt$orig.ident) 
  
  #可视化过滤后的情况
  feats <- c("nFeature_RNA", "nCount_RNA")
  p1_filtered=VlnPlot(input_sce.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 2) + 
    NoLegend()
  w=length(unique(input_sce.filt$orig.ident))/3+5;w 
  ggsave(filename="Vlnplot1_filtered.pdf",plot=p1_filtered,width = w,height = 5)
  
  feats <- c("percent_mito", "percent_ribo", "percent_hb")
  p2_filtered=VlnPlot(input_sce.filt, group.by = "orig.ident", features = feats, pt.size = 0, ncol = 3) + 
    NoLegend()
  w=length(unique(input_sce.filt$orig.ident))/2+5;w 
  ggsave(filename="Vlnplot2_filtered.pdf",plot=p2_filtered,width = w,height = 5) 
  return(input_sce.filt)   
}
```

