#' compare the cellphonedb results of two groups
#'
#' @param group1.name 第一组叫什么名字
#' @param group2.name 第二组叫什么名字
#' @param group1.pfile 第一组两个文件的路径
#' @param group1.mfile
#' @param group2.pfile 第二组两个文件的路径
#' @param group2.mfile
#' @param p.threshold
#' @param thre 类似log2FC
#' @param cell.pair 想展示的细胞对
#' @param gene.pair 想展示的gene对
#' @param plot.width
#' @param plot.height
#' @param filename
#'
#' @return a pdf file and a xlsx file
#' @import tidyverse
#' @import RColorBrewer
#' @import scales
#' @import reshape2
#' @import xlsx
#' @export
#'
#' @examples
ccc_compare=function(
    group1.name=NULL,#第一组叫什么名字
    group2.name=NULL,#第二组叫什么名字
    #默认进行【第一组】比【第二组】

    #第一组两个文件的路径
    group1.pfile=NULL,
    group1.mfile=NULL,
    #第二组两个文件的路径
    group2.pfile=NULL,
    group2.mfile=NULL,

    p.threshold = 0.05,
    thre=1, #类似log2FC

    cell.pair=NULL,#想展示的细胞对
    gene.pair=NULL,#想展示的gene对

    plot.width=NULL,plot.height=NULL,#图片的参数

    filename=""
){
  #library(tidyverse)
  #library(RColorBrewer)
  #library(scales)
  #library(reshape2)

  ### 先定义几个重要的函数 #######################################################
  ### 读取函数，整理cellphonedb的输出结果
  getdata=function(
    pfile, #pvalues.txt文件路径
    mfile #means.txt文件路径
  ){
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
    raw.df
  }

  ### 画图函数
  getplot=function(df,geneA_geneB=NULL,cellA_cellB=NULL){
    if(! is.null(cellA_cellB)){
      df=df[df$group1_cellA_cellB %in% cellA_cellB,]
    }
    if(! is.null(geneA_geneB)){
      df=df[df$group1_geneA_geneB %in% geneA_geneB,]
    }

    df=df[,c("group1_geneA_geneB", "group1_cellA_cellB", "group1mean_to_group2mean_log2")]
    #选取有用的列，然后重新命名。
    #这里我没有选择p值相关的列，因为显著性这个东西只要达到一定的标准，具体的大小就已经不重要了
    colnames(df)=c("geneA_geneB","cellA_cellB","group1mean_to_group2mean_log2")

    if (! is.null(cellA_cellB)) {
      df$cellA_cellB=factor(as.character(df$cellA_cellB),levels = cellA_cellB)
    }else{
      df$cellA_cellB=factor(as.character(df$cellA_cellB),levels = sort(unique(as.character(df$cellA_cellB))))
    }
    if(! is.null(geneA_geneB)) {
      df$geneA_geneB=factor(as.character(df$geneA_geneB),levels = geneA_geneB)
    } else {
      df$geneA_geneB=factor(as.character(df$geneA_geneB),levels = sort(unique(as.character(df$geneA_geneB)),decreasing = T))
    }
    #上述代码是对细胞对、基因对进行排序

    p=ggplot(df,aes(cellA_cellB,geneA_geneB))+
      geom_point(aes(size=abs(group1mean_to_group2mean_log2),color=group1mean_to_group2mean_log2))+
      scale_color_gradient2(paste0("log2(",group1.name,"mean_to_",group2.name,"mean)"),high = "#ee3a2c",low = "#5284c1",mid = "white")+
      scale_size_continuous(paste0("abs(log2(",group1.name,"mean_to_",group2.name,"mean))"))+
      scale_x_discrete("")+scale_y_discrete("")+
      theme_bw()+
      theme(axis.text.x.bottom = element_text(hjust = 1, angle = 45, size=12, color = "black"),
            axis.text.y.left = element_text(size=12, color = "black"))
    p
  }

  ### 读取两组cellphonedb的结果 ##################################################
  group1.df=getdata(pfile = group1.pfile,
                    mfile = group1.mfile)
  colnames(group1.df)=paste0("group1_",colnames(group1.df))
  colnames(group1.df)[1]="geneA_geneB_cellA_cellB"

  group2.df=getdata(pfile = group2.pfile,
                    mfile = group2.mfile)
  colnames(group2.df)=paste0("group2_",colnames(group2.df))
  colnames(group2.df)[1]="geneA_geneB_cellA_cellB"

  both.df=dplyr::inner_join(group1.df,group2.df,by="geneA_geneB_cellA_cellB")
  #只比较两组中都出现的互作关系(geneA_geneB_cellA_cellB)，只在一组中出现的互作关系没法【定量】比较

  ### 过滤 #######################################################################
  both.df=dplyr::filter(both.df,group1_pvalue < p.threshold & group2_pvalue < p.threshold)
  #阈值可改
  #筛选出：在两组结果中，都较显著的互作关系

  both.df$p1xp2_log10_neg=both.df$group1_neg_log10 + both.df$group2_neg_log10
  #综合两个p值的指标

  both.df$group1mean_to_group2mean_log2=log2(both.df$group1_means_exp / both.df$group2_means_exp)
  #第一组比第二组
  both.df$group1mean_to_group2mean_log2[both.df$group1mean_to_group2mean_log2 > 3] = 3
  both.df$group1mean_to_group2mean_log2[both.df$group1mean_to_group2mean_log2 < (-3)] = -3
  #调整数值范围

  ### 第一组比第二组高(低)多少，才算高(低)，才被筛选出来 #################################
  both.df.f=both.df[abs(both.df$group1mean_to_group2mean_log2) > thre,]
  #阈值可改

  ### 挑出想要展示的细胞类型 ######################################################################
  if(is.null(cell.pair)){
    selected.df=both.df.f
  }else{
    selected.df=both.df.f[both.df.f$group1_cellA_cellB %in% cell.pair,]
  }

  ### 保存数据 ####################################################################
  #library(xlsx)
  xlsx::write.xlsx(selected.df,file = paste0(filename,group1.name,2,group2.name,".xlsx"),col.names = T,row.names = F)
  #这个表格还没有对基因对进行过滤，所以可能比图上面的基因对多；这有利于自定制绘图

  ### 画图 ########################################################################
  getplot(selected.df,cellA_cellB = cell.pair,geneA_geneB = gene.pair)
  ggplot2::ggsave(paste0(filename,group1.name,2,group2.name,".pdf"),width = plot.width,height = plot.height,units = "cm")

  return("ok!")
}
