#' Venn plot for the common CpG island and CpG island shore.
#' 
#' @import VennDiagram
#' 
#' @description This function outputs the venn plot for the common CpG island and CpG island shore regions 
#' that are covered by methylated cytosine sites based on R package VennDiagram.
#' 
#' @param genefeatureall_cpgfeature refers to the input file with Q values of two features.
#' @param title refers to the figure title, with default "Venn plot".
#' @param fillcolor refers to the filled color, with default "cornflowerblue"and "green".
#' 
#' @return Outputs a venn figure in two features.
#'
#' @examples
#' Venn_plot(genefeatureall_cpgfeature)
#' Venn_plot(genefeatureall_cpgfeature_Qvalue)
#' Venn_plot(genefeatureall_cpgfeature_Qvalue, fillcolor = c("red","blue"))
#' 
#' @export


Venn_plot <- function(genefeatureall_cpgfeature_Qvalue, title = "Venn plot", fillcolor = c("cornflowerblue","green")){
  
  # combine chromosome and position as the unique id #
  all_feature1 <- genefeatureall_cpgfeature_Qvalue[genefeatureall_cpgfeature_Qvalue$Methdiff1 != "NaN", ]
  all_feature2 <- genefeatureall_cpgfeature_Qvalue[genefeatureall_cpgfeature_Qvalue$Methdiff2 != "NaN", ]
  
  venn.diagram(x = list(A = paste(all_feature1$start, all_feature1$chr, sep = ""),
                        B = paste(all_feature2$start, all_feature2$chr, sep = "")), 
               filename = "Venn plot.tif", imagetype = "tiff", main = title,
               col = "transparent", fill = fillcolor, alpha = 0.50, cex=1, cat.cex=0.8,
               cat.fontface = "bold", category.names = c("CpG island","Shore"))
  
  # figure directory #
  message(paste("The figure is outputed at ", getwd(), "/Venn plot.tif", sep = ""))
}




 