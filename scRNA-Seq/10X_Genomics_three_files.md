# 10X Genomics三文件

10x Genomics 的细胞基因矩阵通常存储在三个关键文件中，这些文件共同描述了单细胞测序数据。这三个文件分别是：

## barcodes.tsv

一列多行，每一行代表一个细胞。

```txt
AAACCCAAGAGGATCC-1
AAACCCAAGCCGCTTG-1
AAACCCAAGGATGGCT-1
AAACCCAGTCTCGGGT-1
AAACCCAGTGCCCGTA-1
AAACCCATCCCAAGCG-1
AAACCCATCCGCGGAT-1
AAACCCATCGAATCCA-1
AAACGAAAGCCTCTGG-1
AAACGAAAGTCATCCA-1
```



## features/genes.tsv

两列多行，每一行代表一个基因，每行第一个是基因ID，第二个是对应的基因symbol名称。

```txt
ENSMUSG00000086053	Gm15178
ENSMUSG00000100764	Gm29155
ENSMUSG00000102095	C730036E19Rik
ENSMUSG00000100635	Gm29157
ENSMUSG00000100480	Gm29156
ENSMUSG00000051285	Pcmtd1
ENSMUSG00000097797	Gm26901
ENSMUSG00000103067	Gm30414
ENSMUSG00000026312	Cdh7
ENSMUSG00000039748	Exo1
```



## matrix.mtx

三列多行，前两行可以理解为标题。从第三行开始，第一个数C1代表基因（即genes/features.tsv中第C1行对应的基因），第二个数C2代表细胞（即barcodes.tsv中第C2行对应的细胞），第三列是表达量。

```
%%MatrixMarket matrix coordinate integer general
%metadata_json: {"software_version": "cellranger-6.1.0", "format_version": 2}
33904 6739 16639697
18 1 1
36 1 1
63 1 2
111 1 1
118 1 1
135 1 1
167 1 1
```





查看三个文件的行数即对应 **细胞数量**，**基因数量**，以及**有表达量的值的数量**

```
$ wc -l *tsv
    6739 barcodes.tsv # 细胞数量
   33904 genes.tsv # 基因数目
   40643 total # 有表达量的值的数目
```















