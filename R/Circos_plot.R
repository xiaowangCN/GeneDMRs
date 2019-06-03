#' Plot the circos.
#' 
#' @import RCircos
#' 
#' @description This function outputs the circos plot for the methylation level and the density of gene, CpG island and CpG island shore on different chromosomes based on R RCircos package. 
#' All the files used in this function should contain chromosome, start position, and end position information that are required for R RCircos package.
#' 
#' @param inputcytofile refers to the output of Cytofile_read() which contains the chromosome information.	
#' @param inputmethfile_QC refers to input file with methylation levels after quality control.
#' @param inputrefseqfile refers to the output of Bedfile_read() which contains the gene information.
#' @param inputcpgifeaturefile refers to the output of Bedfile_read() which contains the CpG island and CpG island shore information.
#' @param labelname refers to the label of gene names which could be the significant genes after Significant_filter(), with default regiongeneall_significant with differentially methylated genes.
#' Sometimes, regiongenealls_significant will have some errors because it has unannotated chromosome name like chrUn_JH584304 or chrUn_NW_018084826v1. Thus, these chromosome names should be removed. 
#' If the labelname is from selfdefinedfile, then the file should contain the headers with chr (chromosome), start (start position), end (end position) and id (gene name).
#' 
#' @param linecolor refers to the colors of the lines plot for different methylation levels, with default NULL (black). 
#' If the linecolor is used, the length of colors should correspond to the length of groups. 
#' 
#' @return Outputs the circus figure with chromosomes, gene labels, the densities of the genes (track 3), CpG islands (track 4) and CpG island shores (track 5)
#' and the methylation levels of different groups from the outermost circle to the innermost circle.
#'
#' @references Hongen Zhang, Paul Meltzer, and Sean Davis. RCircos: an R package for Circos 2D track plots. BMC Bioinformatics, 2013, 14:244.
#' 
#' @examples
#' Circos_plot(inputcytofile, inputmethfile_QC, inputrefseqfile, inputcpgifeaturefile)
#' Circos_plot(inputcytofile, inputmethfile_QC, inputrefseqfile, inputcpgifeaturefile, 
#'             labelname = regiongeneall_significant, linecolor = c("blue1", "orange1", "green1"))
#' Circos_plot(inputcytofile, inputmethfile_QC, inputrefseqfile, inputcpgifeaturefile, 
#'             labelname = selfdefinedfile, linecolor = c("blue", "orange", "green"))
#' 
#' @export


Circos_plot <- function(inputcytofile, inputmethfile_QC, inputrefseqfile, inputcpgifeaturefile, 
                        labelname = regiongeneall_significant, linecolor = NULL){
  
  # divide genome into 1,000,000 windows for better view #
  windowfile <- Window_divide(inputcytofile, windowbp = 1000000)
  
  # get the methylation level for groups #
  windowfilealls <- Methmean_region(inputmethfile_QC, windowfile, chrnum = "all")
  
  # new chromosome + position as Methmean_region delete the unmacthed sites #
  chrpos <- data.frame(chr = windowfilealls$chr, start = windowfilealls$start, end = windowfilealls$end)
  
  groupnum <- length(grep("Methgroup", colnames(windowfilealls)))
  grouppos <- grep("Methgroup", colnames(windowfilealls))
  #for(k in 1:groupnum){
  #  assign(paste("Group", k, sep = "_"), data.frame(chrpos, Meth = windowfilealls[, grouppos[k]]))
  #}
  
  inputcpgi <- filter(inputcpgifeaturefile, cpgfeature %in% "CpGisland")
  inputshore <- filter(inputcpgifeaturefile, cpgfeature %in% "Shores")
  
  # gene, CpG island and shore density #
  out_density <- data.frame(chrpos, density = array(0, c(nrow(chrpos), 3)))
  
  for(i in 1:nrow(chrpos)){
    filter1 <- filter(inputrefseqfile, chr %in% as.vector(unlist(out_density[1]))[i])
    filter2 <- filter(inputcpgi, chr %in% as.vector(unlist(out_density[1]))[i])
    filter3 <- filter(inputshore, chr %in% as.vector(unlist(out_density[1]))[i])
    out_density[i, (ncol(chrpos) + 1)] <- sum(filter1$start > out_density[i,2] & filter1$end <= out_density[i,3])
    out_density[i, (ncol(chrpos) + 2)] <- sum(filter2$start > out_density[i,2] & filter2$end <= out_density[i,3])
    out_density[i, (ncol(chrpos) + 3)] <- sum(filter3$start > out_density[i,2] & filter3$end <= out_density[i,3])
  }
  
  # label for gene names #
  labelnamefile <- data.frame(chr = regiongeneall_significant$chr, start = regiongeneall_significant$start,
                          end = regiongeneall_significant$end, labelname = regiongeneall_significant$id)
  
  
  # full steps in Rcircos #
  # parameter set #
  RCircos.Set.Core.Components(inputcytofile, chr.exclude = NULL, tracks.inside = 7 + groupnum, tracks.outside = 0)
  RC.param <- RCircos.Get.Plot.Parameters()
  RC.param['text.size'] <- 0.4
  RC.param['scatter.color'] <- "purple"
  RC.param['track.background'] <- "white"
  RC.param['grid.line.color'] <- "white"
  RCircos.Reset.Plot.Parameters(RC.param)
  
  # plot area set #
  RCircos.Set.Plot.Area()
  par(mai=c(0, 0.1, 0, 0.1))
  plot.new()
  plot.window(c(-2.5,2.5), c(-2.5, 2.5))
  
  # start plot #
  RCircos.Chromosome.Ideogram.Plot()
  RCircos.Gene.Connector.Plot(labelnamefile, track.num = 1, side = "in", outside.pos = 10)
  RCircos.Gene.Name.Plot(labelnamefile, name.col = 4, track.num = 2,side = "in")
  RCircos.Histogram.Plot(out_density, data.col = 4, track.num = 5, side="in")
  RCircos.Scatter.Plot(out_density, data.col = 5, track.num = 6 , side="in")
  
  # reset scatter parameter #
  RC.param <- RCircos.Get.Plot.Parameters()
  RC.param['scatter.color'] <- "purple1"
  RCircos.Reset.Plot.Parameters(RC.param)
  RCircos.Scatter.Plot(out_density, data.col = 6, track.num = 7, side="in")
  
  # methylation for groups #
  for(j in 1:groupnum){
    RC.param <- RCircos.Get.Plot.Parameters()
    
    # change the line color #
    RC.param['line.color'] <- linecolor[j]
    RCircos.Reset.Plot.Parameters(RC.param)
    RCircos.Line.Plot(windowfilealls[, -1], data.col = grouppos[j] - 1, track.num = 7 + j, side="in")
  }
}




 