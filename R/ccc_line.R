#' compare the cellphonedb results of two groups using connecting plot
#'
#' @param table.path the xlsx file from last step, ccc_compare2 function
#' @param group1.name
#' @param group2.name
#' @param ligand.cell 产生ligand的细胞类型
#' @param receptor.cell 产生receptor的细胞类型
#' @param ligand.color
#' @param receptor.color
#' @param pt.size
#' @param line.thre1 这两个参数可以控制线的（相对）粗细
#' @param line.thre2
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
ccc_line=function(
    table.path=NULL,
    group1.name = NULL,
    group2.name = NULL,
    ligand.cell=NULL,
    receptor.cell=NULL,

    ligand.color="#4dbbd6",
    receptor.color="#90d1c1",

    pt.size=6,
    line.thre1=0.5,
    line.thre2=6,
    file.name="",plot.width=25,plot.height=20
){
  # library(tidyverse)
  # library(xlsx)


  testdata=read.xlsx(table.path,sheetIndex = 1)
  testdata=testdata%>%arrange(group1mean_to_group2mean_log2)
  testdata$ligand=str_replace(testdata$group1_ligand_receptor,"_.*$","")
  testdata$receptor=str_replace(testdata$group1_ligand_receptor,"^.*_","")


  testdata$ligand=factor(testdata$ligand,levels = unique(testdata$ligand))
  testdata$receptor=factor(testdata$receptor,levels = unique(testdata$receptor))
  testdata$ligand2num=as.numeric(testdata$ligand)
  testdata$receptor2num=as.numeric(testdata$receptor)


  testdata$ligand_x=0
  testdata$receptor_x=1
  testdata$ligand2num=(testdata$ligand2num-1) * (max(testdata$receptor2num-1) / max(testdata$ligand2num-1)) + 1


  receptor.df=unique(testdata[,c("receptor","receptor2num","receptor_x")])
  ligand.df=unique(testdata[,c("ligand","ligand2num","ligand_x")])
  test.df=data.frame(x=c(0,0,1,1),
                     y=c(0,max(testdata$receptor2num)+1,0,max(testdata$receptor2num)+1),
                     label=c(
                       paste0("Ligand\n",ligand.cell),
                       paste0(ligand.cell,"\nLigand"),
                       paste0("Receptor\n",receptor.cell),
                       paste0(receptor.cell,"\nReceptor")
                     )
  )


  num1=min(abs(testdata$group1mean_to_group2mean_log2)) %>% floor()
  num2=max(abs(testdata$group1mean_to_group2mean_log2)) %>% ceiling()


  ggplot()+
    geom_segment(data = testdata,mapping = aes(x=ligand_x,y=ligand2num,xend=receptor_x,yend=receptor2num,color=group1mean_to_group2mean_log2,size=abs(group1mean_to_group2mean_log2)))+
    geom_text(data = receptor.df,mapping = aes(x=receptor_x+0.05,y=receptor2num,label=receptor),hjust=0)+
    geom_point(data = receptor.df,mapping = aes(x=receptor_x,y=receptor2num),shape=21,size=pt.size,color="white",fill=receptor.color)+
    geom_text(data = ligand.df,mapping = aes(x=ligand_x-0.05,y=ligand2num,label=ligand),hjust=1)+
    geom_point(data = ligand.df,mapping = aes(x=ligand_x,y=ligand2num),shape=21,size=pt.size,color="white",fill=ligand.color)+
    geom_text(data = test.df,mapping = aes(x=x,y=y,label=label))+
    scale_color_gradient2(paste0("log2(",group1.name,"mean_to_",group2.name,"mean)"),low = "#595eef",mid = "#dadada",high = "#d962cc")+
    scale_x_continuous(expand = c(0,0),limits = c(-0.3,1.3))+
    scale_size_continuous(paste0("abs(log2(",group1.name,"mean_to_",group2.name,"mean))"),limits = c(line.thre1,line.thre2),breaks = seq(num1,num2,1),labels = seq(num1,num2,1))+
    theme_minimal()+
    theme(
      panel.grid = element_blank(),
      axis.title = element_blank(),
      axis.text = element_blank()
    )
  ggsave(paste0(file.name,"line.pdf"),width = plot.width,height = plot.height,units = "cm")
}
