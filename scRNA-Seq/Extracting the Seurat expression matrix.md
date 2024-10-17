# Extracting the Seurat expression matrix

在 Seurat 对象中，基因表达矩阵存储在 Assay 对象中。要输出某个 Assay 中的基因表达矩阵，可以按照以下步骤操作。以下示例假设你的 Seurat 对象名为 *`input_sce`*，并且你想要访问和输出默认的 "RNA" Assay 中的基因表达矩阵。



```R
library(Seurat)

# 获取默认 Assay 的名称
default_assay_name <- DefaultAssay(input_sce)

# 获取默认 Assay 对象
default_assay <- input_sce[[default_assay_name]]

# 访问归一化后的基因表达矩阵
normalized_data <- GetAssayData(default_assay, layer = "data")

# 访问原始计数基因表达矩阵
count_data <- GetAssayData(default_assay, layer = "counts")

# 访问缩放后的基因表达矩阵
scaled_data <- GetAssayData(default_assay, layer = "scale.data")

# 输出归一化后的基因表达矩阵
write.csv(as.matrix(normalized_data), file = "normalized_data.csv")

# 输出原始计数基因表达矩阵
write.csv(as.matrix(count_data), file = "count_data.csv")

# 输出缩放后的基因表达矩阵
write.csv(as.matrix(scaled_data), file = "scaled_data.csv")
```

