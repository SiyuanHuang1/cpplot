#' bubble plot showing the detailed information of cell-cell communications
#'
#' @param pfile pvalues.txt文件路径
#' @param mfile means.txt文件路径
#' @param neg_log10_th 前面画第一张图的时候用的什么值当作显著性阈值，这里保持一致
#' @param means_exp_log2_th 表达均值也限制一下
#' @param notused.cell 确定不包含的细胞，默认是空值
#' @param used.cell 必须含有的细胞，默认是空值
#' @param neg_log10_th2 限定显著性的最大值
#' @param means_exp_log2_th2 限定表达值的范围
#' @param cell.pair 这里是自定义的顺序，若是可选细胞对的子集，则只展示子集，若有交集则只展示交集；空值情况下，会根据可选细胞对自动排序
#' @param gene.pair 作用同上
#' @param color_palette 颜色编码向量，颜色由浅到深
#' @param text_size 坐标轴文本的大小
#'
#' @return a ggplot2 object
#' @import tidyverse
#' @import RColorBrewer
#' @import scales
#' @import reshape2
#' @export
#'
#' @examples
ccc_bubble <- function(
    pfile, #pvalues.txt文件路径
    mfile, #means.txt文件路径
    neg_log10_th= -log10(0.05), #前面画第一张图的时候用的什么值当作显著性阈值，这里保持一致
    means_exp_log2_th=1, #表达均值也限制一下
    notused.cell=NULL, #确定不包含的细胞，默认是空值
    used.cell=NULL, #必须含有的细胞，默认是空值
    neg_log10_th2=3, #限定显著性的最大值
    means_exp_log2_th2=c(-4,6), #限定表达值的范围
    cell.pair=NULL, #这里是自定义的顺序，若是可选细胞对的子集，则只展示子集，若有交集则只展示交集；空值情况下，会根据可选细胞对自动排序
    gene.pair=NULL, #作用同上
    color_palette = c("#313695", "#4575B4", "#ABD9E9", "#FFFFB3", "#FDAE61", "#F46D43", "#D73027", "#A50026"),
    text_size = 12
){
  #####################################################################################################
  #整个函数的逻辑是先根据一组阈值筛选合适的geneA_geneB_cellA_cellB关系，                              #
  #随后根据候选的geneA_geneB、cellA_cellB提取原始矩阵，进行展示，                                     #
  #举个例子，若geneA_geneB_cellA_cellB关系共有7个，占据了3个geneA_geneB pair，4个cellA_cellB pair。   #
  #最终的图形会展示3乘4=12个点，这其中只有7个是满足p值、表达值要求的。                                #
  #                                                                                                   #
  #                                                                 作者：黄思源                      #
  #                                                                 邮箱：huangsiyuan@pku.edu.cn      #
  #####################################################################################################
  #library(tidyverse)
  #library(RColorBrewer)
  #library(scales)
  #library(reshape2)

  #-----------------------------------------------------------
  pvalues=read.table(pfile,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
  pvalues=pvalues[,c(2,12:dim(pvalues)[2])]
  RMpairs=names(sort(table(pvalues$interacting_pair))[sort(table(pvalues$interacting_pair)) > 1])
  pvalues=pvalues[!(pvalues$interacting_pair %in% RMpairs),]
  pvalues.df1=reshape2::melt(pvalues,id="interacting_pair")
  colnames(pvalues.df1)=c("geneA_geneB","cellA_cellB","pvalue")
  pvalues.df1$neg_log10=-log10(pvalues.df1$pvalue)
  pvalues.df1$geneA_geneB_cellA_cellB=paste(pvalues.df1$geneA_geneB,pvalues.df1$cellA_cellB,sep = ",")

  means=read.table(mfile,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
  means=means[,c(2,12:dim(means)[2])]
  rmpairs=names(sort(table(means$interacting_pair))[sort(table(means$interacting_pair)) > 1])
  means=means[!(means$interacting_pair %in% rmpairs),]
  means.df1=reshape2::melt(means,id="interacting_pair")
  colnames(means.df1)=c("geneA_geneB","cellA_cellB","means_exp")
  means.df1$geneA_geneB_cellA_cellB=paste(means.df1$geneA_geneB,means.df1$cellA_cellB,sep = ",")
  means.df1=means.df1[,c("geneA_geneB_cellA_cellB","means_exp")]

  raw.df=merge(pvalues.df1,means.df1,by="geneA_geneB_cellA_cellB")
  raw.df$means_exp_log2=log2(raw.df$means_exp)

  #-----------------------------------------------------------
  #根据第一组阈值筛选
  final.df=dplyr::filter(raw.df,neg_log10 > neg_log10_th & means_exp_log2 > means_exp_log2_th)

  final.df$geneA=stringr::str_replace(final.df$geneA_geneB,"_.*$","") #此处的geneA geneB有不妥之处，
  final.df$geneB=stringr::str_replace(final.df$geneA_geneB,"^.*_","") #没有考虑到complex的命名规则
  final.df$cellA=stringr::str_replace(final.df$cellA_cellB,"\\|.*$","") #后面尽量不用
  final.df$cellB=stringr::str_replace(final.df$cellA_cellB,"^.*\\|","")

  #-----------------------------------------------------------
  #根据其它规则过滤
  #有明确不呈现在图中的细胞，过滤如下
  if (!is.null(notused.cell)) {
    final.df=final.df[!(final.df$cellA %in% notused.cell),]
    final.df=final.df[!(final.df$cellB %in% notused.cell),]
  }
  #两种细胞相同的话一般不展示
  final.df=final.df[!(final.df$cellA==final.df$cellB),]

  #-----------------------------------------------------------
  #提取原始矩阵
  final.df.gene=unique(final.df$geneA_geneB)
  final.df.cell=unique(final.df$cellA_cellB)
  #pair中必须含有的细胞
  if (!is.null(used.cell)){
    tmp_cell=c()
    for (i in used.cell) {
      tmp_cell=union(tmp_cell,final.df.cell[str_detect(final.df.cell,i)])
    }
    final.df.cell=tmp_cell
  }
  raw.df=raw.df[raw.df$geneA_geneB %in% final.df.gene, ]
  raw.df=raw.df[raw.df$cellA_cellB %in% final.df.cell, ]

  #-----------------------------------------------------------
  #范围修正
  raw.df$neg_log10=ifelse(raw.df$neg_log10 > neg_log10_th2,neg_log10_th2,raw.df$neg_log10)
  raw.df$means_exp_log2=ifelse(
    raw.df$means_exp_log2 > means_exp_log2_th2[2],
    means_exp_log2_th2[2],
    ifelse(
      raw.df$means_exp_log2 < means_exp_log2_th2[1],
      means_exp_log2_th2[1],
      raw.df$means_exp_log2
    )
  )

  #-----------------------------------------------------------
  #cellA-cellB的排列顺序
  raw.df$cellA_cellB=as.character(raw.df$cellA_cellB)
  if (!is.null(cell.pair)) {
    tmp_pair=intersect(cell.pair,unique(raw.df$cellA_cellB))
    raw.df=raw.df[raw.df$cellA_cellB %in% tmp_pair,]
    raw.df$cellA_cellB=factor(raw.df$cellA_cellB,levels = tmp_pair)
  } else {
    tmp_pair=sort(unique(raw.df$cellA_cellB))
    raw.df$cellA_cellB=factor(raw.df$cellA_cellB,levels = tmp_pair)
  }
  #geneA-geneB的排列顺序
  raw.df$geneA_geneB=as.character(raw.df$geneA_geneB)
  if (!is.null(gene.pair)) {
    tmp_pair=intersect(gene.pair,unique(raw.df$geneA_geneB))
    raw.df=raw.df[raw.df$geneA_geneB %in% tmp_pair,]
    raw.df$geneA_geneB=factor(raw.df$geneA_geneB,levels = tmp_pair)
  } else {
    tmp_pair=sort(unique(raw.df$geneA_geneB))
    raw.df$geneA_geneB=factor(raw.df$geneA_geneB,levels = tmp_pair)
  }

  #-----------------------------------------------------------
ggplot(raw.df,aes(cellA_cellB,geneA_geneB))+
    geom_point(aes(size=neg_log10,color=means_exp_log2))+
    scale_color_gradientn("log2 mean\n(molecule 1, \nmolecule 2)",colors = color_palette)+
    scale_size_continuous("-log10(p value)")+
    theme_bw()+
    theme(
      axis.title = element_blank(),
      axis.text.x.bottom = element_text(hjust = 1, angle = 45, size=text_size, color = "black"),
      axis.text.y.left = element_text(size = text_size,color = "black"),
      axis.ticks.length = unit(0.15,"cm")
    )
}
