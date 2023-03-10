% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccc_line2.R
\name{ccc_line2}
\alias{ccc_line2}
\title{compare the cellphonedb results of two groups using connecting plot}
\usage{
ccc_line2(
  cpdb.table.path = NULL,
  marker_group = NULL,
  ligand.cell = NULL,
  receptor.cell = NULL,
  group1.name = NULL,
  group2.name = NULL,
  line.size = 2,
  file.name = "",
  plot.width = 25,
  plot.height = 20
)
}
\arguments{
\item{cpdb.table.path}{the xlsx file from last step, ccc_compare2 function}

\item{marker_group}{按照我公众号“代码09”教程得到的DEG数据框}

\item{ligand.cell}{产生ligand的细胞类型}

\item{receptor.cell}{产生receptor的细胞类型}

\item{plot.height}{}
}
\value{
a pdf file
}
\description{
compare the cellphonedb results of two groups using connecting plot
}
