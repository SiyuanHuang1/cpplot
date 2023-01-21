#' heatmap showing the total number of interactions
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
ccc_number_heatmap2=function(
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
  pvalues=pvalues[,12:dim(pvalues)[2]]
  statdf=as.data.frame(colSums(pvalues < pvalue.threshold))
  colnames(statdf)=c("number")

  statdf$indexb=stringr::str_replace(rownames(statdf),"^.*\\|","")
  statdf$indexa=stringr::str_replace(rownames(statdf),"\\|.*$","")
  statdf$total_number=0

  for (i in 1:dim(statdf)[1]) {
    tmp_indexb=statdf[i,"indexb"]
    tmp_indexa=statdf[i,"indexa"]
    if (tmp_indexa == tmp_indexb) {
      statdf[i,"total_number"] = statdf[i,"number"]
    } else {
      statdf[i,"total_number"] = statdf[statdf$indexb==tmp_indexb & statdf$indexa==tmp_indexa,"number"]+
        statdf[statdf$indexa==tmp_indexb & statdf$indexb==tmp_indexa,"number"]
    }
  }

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
    limits.max = ceiling(max(statdf$total_number) / 10) * 10
  }else{
    statdf$total_number[statdf$total_number > ccc.number.max] = ccc.number.max
    limits.max = ccc.number.max
  }

  ggplot(statdf,aes(x=indexa,y=indexb,fill=total_number))+
    geom_tile(color="white")+
    scale_fill_gradientn(colours = color.palette,limits=c(0,limits.max))+
    scale_x_discrete("cluster 1")+
    scale_y_discrete("cluster 2")+
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
