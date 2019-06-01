#' Heat map plot for chromosomes and features.
#' 
#' @import pheatmap
#' 
#' @description This function outputs the heat map plot for methylation level in different chromosomes of differentially methylated genes with features based on R package pheatmap.
#' 
#' @param regiongeneall_significant refers to the input file of methylation levels with differentially methylated genes or the genes in different gene body features.
#' @param featurename refers to the feature name of the output file from Significant_filter() for genefeatureall_cpgfeature file, that is "CpGisland" or "Shore", with default NULL.
#' @param title refers to the figure title, with the default "Methylation level (%)".
#' @param display_numbers refers to whether to display the methylation value in the figure, with default FALSE.
#' @param number_format refers to the displayed number of the methylation value in round format.
#' @param cluster_rows refers to whether to cluster the row, with the default FALSE.
#' @param cluster_cols refers to whether to cluster the column, with the default TRUE.
#' @param gaps_row refers to whether to divide the row, with the default c(1,2) that divide the rows into three parts by row 1 and row 2.
#' @param gaps_col refers to whether to divide the column, with the NULL.
#' 
#' @return Outputs a heat map figure with methylation levels.
#'
#' @examples
#' Heatmap_plot(regiongeneall_significant)
#' Heatmap_plot(genefeatureall_cpgfeature_significantcpgisland, featurename = "CpGisland", display_numbers = FALSE, title = "Methylation level (%) for genes with CpG island")
#' Heatmap_plot(genefeatureall_cpgfeature_significantshore, featurename = "Shore", title = "Methylation level (%) for genes with shores")
#' Heatmap_plot(genefeatureall_cpgfeature_significantshore, featurename = "Shore", title = "Methylation level (%) for genes with shores", cluster_cols = FALSE)
#' Heatmap_plot(regiongeneall_significant, title = "Methylation level (%) for genes", display_numbers = FALSE)
#' Heatmap_plot(regiongeneall_significant, title = "Methylation level (%) for genes", display_numbers = FALSE, cluster_rows = TRUE, gaps_row = NULL)
#' 
#' @export


Heatmap_plot <- function(regiongeneall_significant, featurename = NULL, title = "Methylation level (%)", 
                         display_numbers = FALSE, number_format = "%.0f", cluster_rows = FALSE, cluster_cols = TRUE, 
                         gaps_row = c(1,2), gaps_col = NULL){
  
  if(is.null(featurename) == TRUE){
    
    # calculate the total group number and get the methylation position #
    groupnum <- length(grep("group", colnames(regiongeneall_significant))) / 2
    grouppos <- grep("group", colnames(regiongeneall_significant))[1:groupnum]
    
    methheat <- array(0,c(groupnum, nrow(regiongeneall_significant)))
    methheat[1:groupnum, ] <- t(regiongeneall_significant[, grouppos])
    colnames(methheat) <- paste(1:ncol(methheat), colnames(regiongeneall_significant$id), sep = "_")
    rownames(methheat) <- paste("Group", 1:groupnum, sep = "")
    
    # factor the chromosome #
    annotation_col = data.frame(Chromosome = regiongeneall_significant$chr)
    
  }else{
    
    # calculate the total group number and get the methylation position #
    groupnum <- length(grep(featurename, colnames(regiongeneall_significant))) / 2
    grouppos <- grep(featurename, colnames(regiongeneall_significant))[1:groupnum]
    
    methheat <- array(0,c(groupnum, nrow(regiongeneall_significant)))
    methheat[1:groupnum, ] <- t(regiongeneall_significant[, grouppos])
    colnames(methheat) <- paste(1:ncol(methheat), colnames(regiongeneall_significant$id), sep = "_")
    rownames(methheat) = paste("Group", 1:groupnum, sep = "")
    
    # factor the chromosome and feature#
    annotation_col = data.frame(Chromosome = regiongeneall_significant$chr, Feature = regiongeneall_significant$feature)
  }
  
  rownames(annotation_col) <- colnames(methheat)
  
  # factor the group #
  annotation_row = data.frame(Group = rownames(methheat))
  rownames(annotation_row) <- rownames(methheat)
  
  # heat map #
  pheatmap(methheat*100, main = title, display_numbers = display_numbers, number_format = number_format, cluster_rows = cluster_rows,
           cluster_cols = cluster_cols, annotation_col = annotation_col, annotation_row = annotation_row, gaps_row = gaps_row,
           gaps_col = gaps_col)
}




 