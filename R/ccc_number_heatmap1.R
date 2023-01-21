#' heatmap showing the number of interactions
#'
#' @param pfile the file path of "pvalues.txt"
#' @param pvalue.threshold
#' @param order.of.celltype
#' @param ccc.number.max 控制图中数值范围
#' @param size.of.text 长度为4的数值向量，依次表示坐标轴title、text、图例title、图例text的文本大小
#' @param color.palette 配色方案，长度为3的颜色编码向量
#'
#' @return a ggplot2 object
#' @import tidyverse
#' @import RColorBrewer
#' @import scales
#' @export
#'
#' @examples
ccc_number_heatmap1=function(
    pfile,
    pvalue.threshold = 0.05,
    order.of.celltype = NULL,
    ccc.number.max = NULL,
    size.of.text = c(18,14,14,12),
    color.palette = c("#4393C3","#ffdbba","#B2182B")
){
  #library(tidyverse)
  #library(RColorBrewer)
  #library(scales)

  pvalues=read.table(pfile,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
  pvalues=pvalues[,12:dim(pvalues)[2]] #此时不关注前11列
  statdf=as.data.frame(colSums(pvalues < pvalue.threshold)) #统计在某一种细胞pair的情况之下，显著的受配体pair的数目；阈值可以自己选
  colnames(statdf)=c("number")

  #排在前面的分子定义为indexa；排在后面的分子定义为indexb
  statdf$indexb=stringr::str_replace(rownames(statdf),"^.*\\|","")
  statdf$indexa=stringr::str_replace(rownames(statdf),"\\|.*$","")
  #设置合适的细胞类型的顺序
  if(is.null(order.of.celltype)){
    rankname=sort(unique(statdf$indexa))
  } else {
    rankname=order.of.celltype
  }

  #转成因子类型，画图时，图形将按照预先设置的顺序排列
  statdf$indexa=factor(statdf$indexa,levels = rankname)
  statdf$indexb=factor(statdf$indexb,levels = rankname)

  #调整图上数值的范围
  if(is.null(ccc.number.max)){
    limits.max = ceiling(max(statdf$number) / 10) * 10
  }else{
    statdf$number[statdf$number > ccc.number.max] = ccc.number.max
    limits.max = ccc.number.max
  }

  ggplot(statdf,aes(x=indexa,y=indexb,fill=number))+
    geom_tile(color="white")+
    scale_fill_gradientn(colours = color.palette,limits=c(0,limits.max))+
    scale_x_discrete("cluster 1 produces molecule 1")+
    scale_y_discrete("cluster 2 produces molecule 2")+
    theme_minimal()+
    theme(
      axis.title = element_text(size = size.of.text[1]),
      axis.text.x.bottom = element_text(hjust = 1, vjust = NULL, angle = 45,size = size.of.text[2],color = "black"),
      axis.text.y.left = element_text(size = size.of.text[2],color = "black"),
      legend.title = element_text(size = size.of.text[3]),
      legend.text = element_text(size = size.of.text[4]),
      panel.grid = element_blank()
    )
}
