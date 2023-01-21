#' compare the cellphonedb results of two groups after considering the L-R direction
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
#' @param gene.pair 想展示的gene对（注意：改变gene.pair会改变PDF的结果，但不会改变Excel；表格始终是gene.pair为null的结果，因为这个结果后续画别的图会用到，所以没让它变）
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
ccc_compare2=function(
    group1.name=NULL,#第一组叫什么名字
    group2.name=NULL,#第二组叫什么名字
    #默认进行第一组比第二组

    #第一组两个文件的路径
    group1.pfile=NULL,
    group1.mfile=NULL,
    #第二组两个文件的路径
    group2.pfile=NULL,
    group2.mfile=NULL,

    p.threshold = 0.05,
    thre=1, #类似log2FC

    cell.pair=NULL,#想展示的细胞对
    gene.pair=NULL,#想展示的gene对（注意：改变gene.pair会改变PDF的结果，但不会改变Excel；表格始终是gene.pair为null的结果，因为这个结果后续画别的图会用到，所以没让它变）

    plot.width=NULL,plot.height=NULL,#图片的参数

    filename=""
){
  # library(tidyverse)
  # library(RColorBrewer)
  # library(scales)
  # library(reshape2)

  ### 先定义几个重要的函数 #######################################################
  ### 读取函数，整理cellphonedb的输出结果
  getdata=function(
    pfile, #pvalues.txt文件路径
    mfile #means.txt文件路径
  ){
    ### 整理p值矩阵
    #对行筛选
    pvalues=read.table(pfile,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
    pvalues=pvalues[!str_detect(pvalues$partner_a,"complex"),]
    pvalues=pvalues[!str_detect(pvalues$partner_b,"complex"),]
    pvalues=pvalues[! (pvalues$receptor_a == "True" & pvalues$receptor_b == "True"),]
    pvalues=pvalues[! (pvalues$receptor_a == "False" & pvalues$receptor_b == "False"),]
    pvalues=pvalues[,c(2,9,12:dim(pvalues)[2])]


    #cellA|cellB 细胞类型相同的列不参与后续分析
    tmpdf=as.data.frame(colnames(pvalues)[-c(1,2)])
    colnames(tmpdf)="cellA_cellB"
    tmpdf$index=rownames(tmpdf)
    tmpdf$cellA=str_replace(tmpdf$cellA_cellB,"\\|.*$","")
    tmpdf$cellB=str_replace(tmpdf$cellA_cellB,"^.*\\|","")
    tmpdf=tmpdf[!tmpdf$cellA == tmpdf$cellB,]
    tmpdf$index=as.numeric(tmpdf$index)
    tmpdf$index=tmpdf$index+2
    pvalues=pvalues[,c(1,2,tmpdf$index)]


    #构建数据框，每行包含一对替换值的位置
    tmpdf$index=1:length(tmpdf$index)
    tmpdf$index=tmpdf$index+2
    tmpdf$compartment=0
    tmpdf$cellAorcellB=""
    for (i in 1:length(tmpdf$index)) {
      rank=which(tmpdf$cellA_cellB %in% paste0(tmpdf[i,"cellB"],"|",tmpdf[i,"cellA"]))
      tmpdf[i,"compartment"] = tmpdf[rank,"index"]

      tmpdf[i,"cellAorcellB"] = paste(sort(c(tmpdf[i,"cellA"],tmpdf[i,"cellB"])),collapse = "or")
    }
    tmpdf=tmpdf[!duplicated(tmpdf$cellAorcellB),]


    #b不是receptor时，人为更换gene pair的顺序，以及调换cell pair值
    for (i in 1:dim(pvalues)[1]) {
      if(pvalues[i,"receptor_b"] == "False"){

        pvalues[i,"interacting_pair"]=paste0(
          str_split(pvalues[i,"interacting_pair"],"_")[[1]][2],"_",
          str_split(pvalues[i,"interacting_pair"],"_")[[1]][1]
        )

        for (j in 1:dim(tmpdf)[1]) {
          tmpvalue=pvalues[i,tmpdf[j,"index"]]
          pvalues[i,tmpdf[j,"index"]]=pvalues[i,tmpdf[j,"compartment"]]
          pvalues[i,tmpdf[j,"compartment"]]=tmpvalue
        }

        pvalues[i,"receptor_b"] = "True"
      }
    }
    pvalues$receptor_b=NULL


    #去除重复的gene pair
    RMpairs=names(sort(table(pvalues$interacting_pair))[sort(table(pvalues$interacting_pair)) > 1])
    pvalues=pvalues[!(pvalues$interacting_pair %in% RMpairs),]
    pvalues.df1=melt(pvalues,id="interacting_pair")
    colnames(pvalues.df1)=c("ligand_receptor","cellA_cellB","pvalue")
    pvalues.df1$neg_log10=-log10(pvalues.df1$pvalue)
    pvalues.df1$ligand_receptor_cellA_cellB=paste(pvalues.df1$ligand_receptor,pvalues.df1$cellA_cellB,sep = ",")


    ### 整理均值矩阵
    #对行筛选
    means=read.table(mfile,header = T,sep = "\t",stringsAsFactors = F,check.names = F)
    means=means[!str_detect(means$partner_a,"complex"),]
    means=means[!str_detect(means$partner_b,"complex"),]
    means=means[! (means$receptor_a == "True" & means$receptor_b == "True"),]
    means=means[! (means$receptor_a == "False" & means$receptor_b == "False"),]
    means=means[,c(2,9,12:dim(means)[2])]


    #细胞类型相同的列不参与后续分析
    tmpdf=as.data.frame(colnames(means)[-c(1,2)])
    colnames(tmpdf)="cellA_cellB"
    tmpdf$index=rownames(tmpdf)
    tmpdf$cellA=str_replace(tmpdf$cellA_cellB,"\\|.*$","")
    tmpdf$cellB=str_replace(tmpdf$cellA_cellB,"^.*\\|","")
    tmpdf=tmpdf[!tmpdf$cellA == tmpdf$cellB,]
    tmpdf$index=as.numeric(tmpdf$index)
    tmpdf$index=tmpdf$index+2
    means=means[,c(1,2,tmpdf$index)]


    #构建数据框，每行包含一对替换值的位置
    tmpdf$index=1:length(tmpdf$index)
    tmpdf$index=tmpdf$index+2
    tmpdf$compartment=0
    tmpdf$cellAorcellB=""
    for (i in 1:length(tmpdf$index)) {
      rank=which(tmpdf$cellA_cellB %in% paste0(tmpdf[i,"cellB"],"|",tmpdf[i,"cellA"]))
      tmpdf[i,"compartment"] = tmpdf[rank,"index"]

      tmpdf[i,"cellAorcellB"] = paste(sort(c(tmpdf[i,"cellA"],tmpdf[i,"cellB"])),collapse = "or")
    }
    tmpdf=tmpdf[!duplicated(tmpdf$cellAorcellB),]


    #b不是receptor时，人为更换gene pair的顺序，以及调换cell pair的两个值
    for (i in 1:dim(means)[1]) {
      if(means[i,"receptor_b"] == "False"){

        means[i,"interacting_pair"]=paste0(
          str_split(means[i,"interacting_pair"],"_")[[1]][2],"_",
          str_split(means[i,"interacting_pair"],"_")[[1]][1]
        )

        for (j in 1:dim(tmpdf)[1]) {
          tmpvalue=means[i,tmpdf[j,"index"]]
          means[i,tmpdf[j,"index"]]=means[i,tmpdf[j,"compartment"]]
          means[i,tmpdf[j,"compartment"]]=tmpvalue
        }

        means[i,"receptor_b"] = "True"
      }
    }
    means$receptor_b=NULL


    #去除重复的gene pair
    rmpairs=names(sort(table(means$interacting_pair))[sort(table(means$interacting_pair)) > 1])
    means=means[!(means$interacting_pair %in% rmpairs),]
    means.df1=melt(means,id="interacting_pair")
    colnames(means.df1)=c("ligand_receptor","cellA_cellB","means_exp")
    means.df1$ligand_receptor_cellA_cellB=paste(means.df1$ligand_receptor,means.df1$cellA_cellB,sep = ",")
    means.df1=means.df1[,c("ligand_receptor_cellA_cellB","means_exp")]


    ###合并两个矩阵
    raw.df=merge(pvalues.df1,means.df1,by="ligand_receptor_cellA_cellB")
    raw.df$means_exp_log2=log2(raw.df$means_exp)
    raw.df
  }

  ### 画图函数
  getplot=function(df,ligand_receptor=NULL,cellA_cellB=NULL){
    if(! is.null(cellA_cellB)){
      df=df[df$group1_cellA_cellB %in% cellA_cellB,]
    }
    if(! is.null(ligand_receptor)){
      df=df[df$group1_ligand_receptor %in% ligand_receptor,]
    }

    df=df[,c("group1_ligand_receptor", "group1_cellA_cellB", "group1mean_to_group2mean_log2")]
    #选取有用的列，然后重新命名。
    #这里我没有选择p值相关的列，因为显著性这个东西只要达到一定的标准，具体的大小就已经不重要了
    colnames(df)=c("ligand_receptor","cellA_cellB","group1mean_to_group2mean_log2")

    if (! is.null(cellA_cellB)) {
      df$cellA_cellB=factor(as.character(df$cellA_cellB),levels = cellA_cellB)
    }else{
      df$cellA_cellB=factor(as.character(df$cellA_cellB),levels = sort(unique(as.character(df$cellA_cellB))))
    }
    if(! is.null(ligand_receptor)) {
      df$ligand_receptor=factor(as.character(df$ligand_receptor),levels = ligand_receptor)
    } else {
      df$ligand_receptor=factor(as.character(df$ligand_receptor),levels = sort(unique(as.character(df$ligand_receptor)),decreasing = T))
    }
    #上述代码是对细胞对、基因对进行排序

    p=df%>%ggplot(aes(cellA_cellB,ligand_receptor))+
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
  colnames(group1.df)[1]="ligand_receptor_cellA_cellB"

  group2.df=getdata(pfile = group2.pfile,
                    mfile = group2.mfile)
  colnames(group2.df)=paste0("group2_",colnames(group2.df))
  colnames(group2.df)[1]="ligand_receptor_cellA_cellB"

  both.df=group1.df%>%inner_join(group2.df,by="ligand_receptor_cellA_cellB")
  #只比较两组中都出现的互作关系(ligand_receptor_cellA_cellB)，只在一组中出现的互作关系没法【定量】比较

  ### 过滤 #######################################################################
  both.df=both.df%>%dplyr::filter(group1_pvalue < p.threshold & group2_pvalue < p.threshold)
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
  # library(xlsx)
  write.xlsx(selected.df,file = paste0(filename,group1.name,2,group2.name,".xlsx"),col.names = T,row.names = F)

  ### 画图 ########################################################################
  getplot(selected.df,cellA_cellB = cell.pair,ligand_receptor = gene.pair)
  ggsave(paste0(filename,group1.name,2,group2.name,".pdf"),width = plot.width,height = plot.height,units = "cm")

  return("ok!")
}
