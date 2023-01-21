#' compare the cellphonedb results of two groups using connecting plot
#'
#' @param cpdb.table.path the xlsx file from last step, ccc_compare2 function
#' @param marker_group 按照我公众号“代码09”教程得到的DEG数据框
#' @param ligand.cell 产生ligand的细胞类型
#' @param receptor.cell 产生receptor的细胞类型
#' @param group1.name
#' @param group2.name
#' @param line.size
#' @param file.name
#' @param plot.width
#' @param plot.height
#'
#' @return a pdf file
#' @import tidyverse
#' @import xlsx
#' @export
#'
#' @examples
ccc_line2=function(
    cpdb.table.path=NULL,
    marker_group=NULL,

    ligand.cell=NULL,
    receptor.cell=NULL,
    group1.name=NULL,
    group2.name=NULL,

    line.size=2,
    file.name="",
    plot.width=25,
    plot.height=20
){
  testdata=read.xlsx(cpdb.table.path,sheetIndex = 1)
  testdata=testdata%>%arrange(group1mean_to_group2mean_log2)
  testdata$ligand=str_replace(testdata$group1_ligand_receptor,"_.*$","")
  testdata$receptor=str_replace(testdata$group1_ligand_receptor,"^.*_","")


  pair.df=testdata[,c("ligand","receptor")]
  ligand=as.character(unique(pair.df$ligand))
  marker_group_ligand=marker_group%>%dplyr::filter(cluster == ligand.cell)
  common.gene=intersect(ligand,marker_group_ligand$gene)
  marker_group_ligand=marker_group_ligand[marker_group_ligand$gene %in% common.gene,]
  rownames(marker_group_ligand)=marker_group_ligand$gene
  marker_group_ligand=marker_group_ligand[common.gene,]
  marker_group_ligand$ligand2num=1:length(marker_group_ligand$gene)
  marker_group_ligand$ligand_x=0


  receptor=as.character(unique(pair.df$receptor))
  marker_group_receptor=marker_group%>%dplyr::filter(cluster == receptor.cell)
  common.gene=intersect(receptor,marker_group_receptor$gene)
  marker_group_receptor=marker_group_receptor[marker_group_receptor$gene %in% common.gene,]
  rownames(marker_group_receptor)=marker_group_receptor$gene
  marker_group_receptor=marker_group_receptor[common.gene,]
  marker_group_receptor$receptor2num=1:length(marker_group_receptor$gene)
  marker_group_receptor$receptor_x=1
  marker_group_ligand$ligand2num=(marker_group_ligand$ligand2num-1) * (max(marker_group_receptor$receptor2num-1) / max(marker_group_ligand$ligand2num-1)) + 1


  testdata2=testdata%>%subset(ligand %in% marker_group_ligand$gene )%>%subset(receptor %in% marker_group_receptor$gene)
  for (i in 1:dim(testdata2)[1]) {
    testdata2$ligand2num[i] = marker_group_ligand[marker_group_ligand$gene == testdata2$ligand[i],"ligand2num"]
    testdata2$receptor2num[i] = marker_group_receptor[marker_group_receptor$gene == testdata2$receptor[i],"receptor2num"]
  }
  testdata2$ligand_x=0
  testdata2$receptor_x=1


  test.df=data.frame(x=c(0,0,1,1),
                     y=c(0,max(marker_group_receptor$receptor2num)+1,0,max(marker_group_receptor$receptor2num)+1),
                     label=c(
                       paste0("Ligand\n",ligand.cell),
                       paste0(ligand.cell,"\nLigand"),
                       paste0("Receptor\n",receptor.cell),
                       paste0(receptor.cell,"\nReceptor")
                     ))


  ggplot()+
    geom_segment(data = testdata2,mapping = aes(x=ligand_x,y=ligand2num,xend=receptor_x,yend=receptor2num,color=group1mean_to_group2mean_log2),size=line.size)+
    scale_color_gradient2(paste0("log2(",group1.name,"mean_to_",group2.name,"mean)"),low = "#595eef",mid = "#dadada",high = "#d962cc")+
    geom_point(data = marker_group_ligand,mapping = aes(x=ligand_x,y=ligand2num,size=-log10(p_val),fill=avg_log2FC),shape=21,color="#969696")+
    geom_text(data = marker_group_ligand,mapping = aes(x=ligand_x-0.05,y=ligand2num,label=gene),hjust=1)+
    geom_point(data = marker_group_receptor,mapping = aes(x=receptor_x,y=receptor2num,size=-log10(p_val),fill=avg_log2FC),shape=21,color="#969696")+
    geom_text(data = marker_group_receptor,mapping = aes(x=receptor_x+0.05,y=receptor2num,label=gene),hjust=0)+
    scale_fill_gradient2(paste0("avg_log2FC_",group1.name,"/",group2.name),high = "#ee3a2c",low = "#5284c1",mid = "white")+
    geom_text(data = test.df,mapping = aes(x=x,y=y,label=label))+
    scale_x_continuous(expand = c(0,0),limits = c(-0.3,1.3))+
    theme_minimal()+
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank()
    )
  ggsave(paste0(file.name,"line2.pdf"),width = plot.width,height = plot.height,units = "cm")
}
