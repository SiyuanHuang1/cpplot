#' connecting plot showing the total number of interactions
#'
#' @param pfile the file path of "pvalues.txt"
#' @param pvalue.threshold
#' @param order.of.celltype
#' @param ccc.number.max 控制图中数值范围
#' @param color.palette 节点的颜色，命名的颜色编码向量
#' @param color.line 颜色编码向量，颜色深度需要递进
#' @param width.parameter 调整连线的宽度
#' @param vertex.label.cex 节点标注字体大小
#' @param vertex.size 节点大小
#'
#' @return
#' @import tidyverse
#' @import RColorBrewer
#' @import scales
#' @import igraph
#' @export
#'
#' @examples
ccc_number_line=function(
    pfile,
    pvalue.threshold = 0.05,
    order.of.celltype = NULL,
    ccc.number.max = NULL,
    color.palette = NULL,
    color.line = NULL, #颜色编码向量，颜色深度需要递进
    width.parameter = 0.5, #调整连线的宽度
    vertex.label.cex = 1, #节点标注字体大小
    vertex.size = 30 #节点大小
){
  #library(tidyverse)
  #library(RColorBrewer)
  #library(scales)
  #library(igraph)

  pvalues=read.table(pfile,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
  pvalues=pvalues[,12:dim(pvalues)[2]]
  statdf=as.data.frame(colSums(pvalues < pvalue.threshold))
  colnames(statdf)=c("number")

  statdf$indexb=stringr::str_replace(rownames(statdf),"^.*\\|","")
  statdf$indexa=stringr::str_replace(rownames(statdf),"\\|.*$","")
  #设置合适的细胞类型的顺序
  if(is.null(order.of.celltype)){
    rankname=sort(unique(statdf$indexa))
  } else {
    rankname=order.of.celltype
  }

  A=c()
  B=c()
  C=c()
  remaining=rankname
  for (i in rankname[-length(rankname)]) {
    remaining=setdiff(remaining,i)
    for (j in remaining) {
      count=statdf[statdf$indexa == i & statdf$indexb == j,"number"]+
        statdf[statdf$indexb == i & statdf$indexa == j,"number"]
      A=append(A,i)
      B=append(B,j)
      C=append(C,count)
    }
  }

  statdf2=data.frame(indexa=A,indexb=B,number=C)
  statdf2=statdf2 %>% rbind(statdf[statdf$indexa==statdf$indexb,c("indexa","indexb","number")])
  statdf2=statdf2[statdf2$number > 0,] #过滤掉值为0的观测

  #调整图上数值的范围
  if(is.null(ccc.number.max)){
    limits.max = ceiling(max(statdf2$number) / 10) * 10
  }else{
    statdf2$number[statdf2$number > ccc.number.max] = ccc.number.max
    limits.max = ccc.number.max
  }

  #设置节点和连线的颜色
  if(is.null(color.palette)){
    color1=scales::hue_pal()(length(rankname))
    names(color1)=rankname
  }else{
    color1=color.palette
  }

  if (is.null(color.line)) {
    color2=colorRampPalette(RColorBrewer::brewer.pal(9, "Reds")[3:7])(limits.max) #将颜色分成多少份，取决于互作关系数目的最大值
    names(color2)=1:limits.max #每一份颜色用对应的数字命名
  } else {
    color2=colorRampPalette(color.line)(limits.max)
    names(color2)=1:limits.max
  }


  #做网络图
  ##下面的四行代码相对固定
  net <- igraph::graph_from_data_frame(statdf2[,c("indexa","indexb","number")])
  edge.start <- igraph::ends(net, es=E(net), names=FALSE)
  group <-  igraph::cluster_optimal(net)
  coords <- igraph::layout_in_circle(net, order = order(membership(group)))

  igraph::E(net)$width <- igraph::E(net)$number * width.parameter #将数值映射到连线的宽度，有时还需要微调，这里除以2就是这个目的
  igraph::E(net)$color <- color2[as.character(igraph::E(net)$number)] #用前面设置好的颜色赋给连线，颜色深浅对应数值大小
  igraph::E(net)$label = igraph::E(net)$number #连线的标注
  igraph::E(net)$label.color <- "black" #连线标注的颜色

  igraph::V(net)$label.color <- "black" #节点标注的颜色
    igraph::V(net)$color <- color1[names(igraph::V(net))] #节点的填充颜色，前面已经设置了；V(net)返回节点信息

    #调整节点位置的线条角度
    ##如果没有这两行代码，节点位置的圆圈是向右的
    loop.angle<-ifelse(
      coords[igraph::V(net),1]>0,
      -atan(coords[igraph::V(net),2]/coords[igraph::V(net),1]),
      pi-atan(coords[igraph::V(net),2]/coords[igraph::V(net),1])
    )
    igraph::E(net)$loop.angle[which(edge.start[,2]==edge.start[,1])] <- loop.angle[edge.start[which(edge.start[,2]==edge.start[,1]),1]]

    plot(
      net,
      edge.arrow.size = 0, #连线不带箭头
      edge.curved = 0, #连线不弯曲
      vertex.frame.color = "black", #节点外框颜色
      layout = coords,
      vertex.label.cex = vertex.label.cex, #节点标注字体大小
      vertex.size = vertex.size #节点大小
    )
}
