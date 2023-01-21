# cpplot
Visualize the results of cell-cell communication analysis based on CellPhoneDB

> 整理了微信公众号“TOP生物信息”上面几篇讲解CellPhoneDB的原创帖子，将其中涉及到的代码稍加整理得到此R包。

[基于CellPhoneDB的细胞通讯分析及可视化 (上篇)——2021-07-24发布](https://blog.csdn.net/qq_38774801/article/details/119060843?spm=1001.2014.3001.5501)

[基于CellPhoneDB的细胞通讯分析及可视化 (下篇)——2021-07-24发布](https://blog.csdn.net/qq_38774801/article/details/119061480?spm=1001.2014.3001.5501)

[【单细胞高级绘图】08.细胞通讯_两组比较_气泡图——2022-08-30发布](https://mp.weixin.qq.com/s?__biz=MzkzMzE5NTM4NA==&mid=2247485827&idx=1&sn=a8429c52caf17526d65f3e7e8f609f0b&chksm=c2517294f526fb8231e4fe4b9c3c53a920a251bf7eacee7d80afea78d39c0e7e412635a0a59d&token=1752071379&lang=zh_CN&scene=21#wechat_redirect)

[【单细胞高级绘图】09.细胞通讯_两组比较_连线图——2022-08-31发布](https://mp.weixin.qq.com/s?__biz=MzkzMzE5NTM4NA==&mid=2247485843&idx=1&sn=741b5436ef7c03fa2ba77335bc2157a9&chksm=c2517284f526fb92512ce955e97372216d9fe29dceff5bdf73eb178eb3cf418e7fec1f14d0e7&token=1752071379&lang=zh_CN&scene=21#wechat_redirect)

> 简而言之，这个包对接的是CellPhoneDB的流程，跑完CellPhoneDB之后，可能需要画几张图，比如：

1. 各种细胞之间互作的数量关系
2. 具体的互作细节（什么细胞之间有什么L-R pair）
3. 如果有两个组都进行了CellPhoneDB的分析，如何比较两组的结果

> 下面代码演示一下

### 0. 下载并加载R包
**最好提前安装几个依赖包：**
RColorBrewer, igraph, reshape2, scales, tidyverse, xlsx

```
devtools::install_github("SiyuanHuang1/cpplot")
library(cpplot)
library(tidyverse)
```
### 1. 各种细胞之间互作的数量关系
这一部分，有三个函数可以实现，分别是：
```
ccc_number_heatmap1(pfile = "test/pvalues.txt") #ggplot对象
ccc_number_heatmap2(pfile = "test/pvalues.txt") #ggplot对象
ccc_number_line(pfile = "test/pvalues.txt",vertex.size = 20) #不是ggplot对象，不能用ggsave保存
```
出图如下：

(图片的解读可以参考我最上面提到的几篇帖子)

![image](https://user-images.githubusercontent.com/37017211/213872028-46679f18-aafa-4dad-86ee-de85db349209.png)
![image](https://user-images.githubusercontent.com/37017211/213872068-44377a1b-3d6b-4f07-8099-6d871d41bd44.png)
![image](https://user-images.githubusercontent.com/37017211/213872101-63653202-78b4-4340-94ec-6edaf8fdfe9d.png)

### 2. 具体的互作细节
```
ccc_bubble(
  pfile="./test/pvalues.txt",
  mfile="./test/means.txt",
  # 下面这些是默认参数，可以不变
  # neg_log10_th = -log10(0.05),
  # means_exp_log2_th = 1,
  # notused.cell = NULL,
  # used.cell = NULL,
  # neg_log10_th2 = 3,
  # means_exp_log2_th2 = c(-4, 6),
  # cell.pair = NULL,
  # gene.pair = NULL,
  # color_palette = c("#313695", "#4575B4", "#ABD9E9", "#FFFFB3", "#FDAE61", "#F46D43","#D73027", "#A50026"),
  # text_size = 12
)
```
![image](https://user-images.githubusercontent.com/37017211/213872479-90b043f3-a6c9-4385-9e73-2346b94e0bc2.png)
```
# 改写参数
ccc_bubble(
  pfile="./test/pvalues.txt",
  mfile="./test/means.txt",
  cell.pair=c("Mcell|Scell","Mcell|NKcell","Mcell|Tcell","Scell|Mcell","NKcell|Mcell","Tcell|Mcell"),
  #这里是自定义的顺序，若是可选细胞对的子集，则只展示子集，若有交集则只展示交集；空值情况下，会根据可选细胞对自动排序
  gene.pair=c("MIF_TNFRSF14","FN1_aVb1 complex","EGFR_MIF")
  #作用同上
)
```
![image](https://user-images.githubusercontent.com/37017211/213872535-652580c8-af8c-4d57-9447-638fe7c8d704.png)

### 3. 两组之间的比较
##### 第1种图
```
### 必要参数
ccc_compare(group1.name = "Old",group2.name = "Young",
            group1.pfile = "cellphonedb/Old/pvalues.txt",group1.mfile="cellphonedb/Old/means.txt",
            group2.pfile="cellphonedb/Young/pvalues.txt",group2.mfile="cellphonedb/Young/means.txt",
            p.threshold = 0.01,thre=1,
            plot.width=105,plot.height=110,filename = "test0121_"
)

### 额外参数
# 比如，这里我想展示EC细胞分别充当cellA和cellB的图
# 也可以指定gene pair
ccc_compare(group1.name = "Old",group2.name = "Young",
            group1.pfile = "cellphonedb/Old/pvalues.txt",group1.mfile="cellphonedb/Old/means.txt",
            group2.pfile="cellphonedb/Young/pvalues.txt",group2.mfile="cellphonedb/Young/means.txt",
            p.threshold = 0.05,thre=1,
            #gene.pair = NULL,
            cell.pair=c(
              paste0("EC|",c("APC","SMC","Mac","DC","Neutrophil")),
              paste0(c("APC","SMC","Mac","DC","Neutrophil"),"|EC")
            ),
            plot.width=18,plot.height=30,filename = "test0121b_"
)
```
![image](https://user-images.githubusercontent.com/37017211/213871689-814b4d9b-3ef8-4b8b-be07-67a4e9c5fa6a.png)

(图片的解读可以参考我最上面提到的几篇帖子)

##### 第2种图
```
ccc_compare2(group1.name = "Old",group2.name = "Young",
             group1.pfile = "cellphonedb/Old/pvalues.txt",group1.mfile="cellphonedb/Old/means.txt",
             group2.pfile="cellphonedb/Young/pvalues.txt",group2.mfile="cellphonedb/Young/means.txt",
             p.threshold = 0.05,thre=0.5,
             cell.pair="EC|APC", #指定ligand产生的细胞|receptor产生的细胞
             plot.width=15,plot.height=30,filename = "test0121_"
)
```
之后会得到一个xlsx表格，画图会用到
```
ccc_line(table.path="test0121_Old2Young.xlsx",ligand.cell="EC",receptor.cell="APC",
         group1.name = "Old",group2.name = "Young",#这五个参数和上一步对应
         ligand.color="#4dbbd6",receptor.color="#90d1c1",
         pt.size=6,
         line.thre1=0.5,line.thre2=6,#line.thre1和上一步的"thre"参数一致，line.thre2可以用来调整线的粗细，值越大，线越细
         file.name="test0121b_",plot.width=25,plot.height=20)
```
然后就能得到这张图：
![image](https://user-images.githubusercontent.com/37017211/213871164-964a1b9a-46c6-4064-b5a3-d1a89963b98b.png)

(图片的解读可以参考我最上面提到的几篇帖子)。

##### 第3种图
```
ccc_compare2(group1.name = "Old",group2.name = "Young",
             group1.pfile = "cellphonedb/Old/pvalues.txt",group1.mfile="cellphonedb/Old/means.txt",
             group2.pfile="cellphonedb/Young/pvalues.txt",group2.mfile="cellphonedb/Young/means.txt",
             p.threshold = 0.05,thre=0.5,
             cell.pair="EC|APC", #指定ligand产生的细胞|receptor产生的细胞
             plot.width=15,plot.height=30,filename = "test0121_"
)
```
这一步跟第2种图一样。后续还要找两组的差异基因
```
library(Seurat)
testseu=readRDS("testseu.rds")
# 此次演示为了加快运行速度，人为减少了数据量，实际分析中找差异基因不建议这么做
selectedCB=sample(testseu@meta.data$CB,1000)
testseu=testseu%>%subset(CB %in% selectedCB)

# 基于分组找差异基因
marker_group=data.frame()
Idents(testseu)="celltype_age"
for ( ci in c("EC","APC") ) {
  tmp.marker <- FindMarkers(
    testseu, logfc.threshold = 0, min.pct = 0.01,
    only.pos = F, test.use = "wilcox",
    ident.1=paste0(ci,"_Old"),ident.2=paste0(ci,"_Young")
  )
  
  tmp.marker$gene=rownames(tmp.marker)
  tmp.marker$cluster_group=ifelse(tmp.marker$avg_log2FC > 0,paste0(ci,"_Old"),paste0(ci,"_Young"))
  tmp.marker$cluster=ci
  tmp.marker=tmp.marker%>%arrange(desc(avg_log2FC))
  
  marker_group=marker_group%>%rbind(tmp.marker)
}
#本次演示的数据集为小鼠数据集，在运行cellphonedb时，进行了基因symbol的转换。
#此处找差异基因得到的symbol为真实基因名，为了让两个分析匹配，DEG表格也应该做基因名转换。
#但是为了简化，此处只是简单地将小鼠基因名转为大写，不是很精确。大家在分析的时候建议严格一点。
marker_group$gene=marker_group$gene %>% toupper()
```
然后借助差异基因，再画图
```
ccc_line2(cpdb.table.path = "test0121_Old2Young.xlsx",marker_group = marker_group,
          ligand.cell = "EC",receptor.cell = "APC",
          group1.name = "Old",group2.name = "Young",
          line.size = 2,file.name = "test0121b_",plot.width = 25,plot.height = 20
)
```
![image](https://user-images.githubusercontent.com/37017211/213871326-ee8baba5-4a17-4b42-9d0c-de1dd9302104.png)

(图片的解读可以参考我最上面提到的几篇帖子)

***
如果对这个R包有什么疑问，可以给我发邮件(huangsiyuan1001@163.com)或者关注我们的公众号(TOP生物信息)到后台反馈。我看GitHub没这两个勤，提issue没这两个高效。

![扫码_搜索联合传播样式-白色版](https://user-images.githubusercontent.com/37017211/213873520-7840d1c6-2ddd-4486-9ce7-aae9cbba333d.png)
