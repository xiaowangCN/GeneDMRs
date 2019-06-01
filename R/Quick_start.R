#' Install the dependencies
#' 
#' @import BiocManager
#' 
#' @description This function install the other R packages for the dependencies of GeneDMRs.
#' 
#' @param Dbannotation refers to the annotation dataset for enrichment, with default "org.Mm.eg.db" of mouse.
#' 
#' @return Outputs a list of installed R packages.
#'
#' @examples
#' Prepare_GeneDMRs(Dbannotation = "org.Mm.eg.db")
#' 
#' @export


Prepare_GeneDMRs <- function(Dbannotation = "org.Mm.eg.db"){
  
  # install the dependencies and library them #
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(c("devtools", "clusterProfiler", "corrplot", "dplyr", "genomation", "Hmisc", "KEGG.db", 
                         "pheatmap", "plotrix", "qqman", "RColorBrewer", "RCircos", "VennDiagram", Dbannotation))
  
  require("clusterProfiler")
  require("corrplot")
  require("devtools")
  require("dplyr")
  require("genomation")
  require("Hmisc")
  require("KEGG.db")
  require("pheatmap")
  require("plotrix")
  require("qqman")
  require("RColorBrewer")
  require("RCircos")
  require("VennDiagram")
  require(Dbannotation)
}





#' Quick use the GeneDMRs package for gene based differentially methlated regions
#' 
#' @import dplyr
#' @import qqman
#' @import pheatmap
#' @import clusterProfiler
#' 
#' @description This function outputs a series of results and figures for gene based regions' methylation analysis.
#' 
#' @param paths refers to the path of input file, with default the package path.
#' @param suffixmeth refers to the suffix of methylation file, e.g. ".gz", ".zip" and so on (some files are in text .txt format, then ".txt" or ".txt.gz"), with default ".gz".
#' @param bedfile refers to the file name of bed file for "refseq". This file is downloaded from UCSC website, with default "refseq".
#' @param suffixbed refers to the suffix of bed file, e.g. ".gz", ".zip" and so on (some files are in text .txt format, then ".txt" or ".txt.gz"), with default ".txt".
#' @param Dbannotation refers to the annotation dataset for enrichment, with default "org.Mm.eg.db" of mouse.
#' @param keggorganism refers to the species name for KEGG enrichment, with default "mmu" of mouse.
#' 
#' @return Outputs a series of DMG results.
#'
#' @examples
#' Quick_GeneDMRs()
#' Quick_GeneDMRs(Dbannotation = "org.Mm.eg.db", keggorganism = "mmu")
#' 
#' @export


Quick_GeneDMRs <- function(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""), suffixmeth = ".gz",
                           bedfile = "refseq", suffixbed = ".txt", Dbannotation = "org.Mm.eg.db", keggorganism = "mmu"){
  
  # read the file #
  inputmethfile <- Methfile_read(paths = paths, suffix = suffixmeth)
  inputrefseqfile <- Bedfile_read(paths = paths, bedfile = bedfile, suffix = suffixbed, feature = FALSE)
  
  # quality control #
  inputmethfile_QC <- Methfile_QC(inputmethfile)
  
  # methylation mean #
  regiongeneall <- Methmean_region(inputmethfile_QC, inputrefseqfile, chrnum = "all")
  
  # statistical test #
  regiongeneall_Qvalue <- Logic_regression(regiongeneall)
  
  # sifnificant filter #
  regiongeneall_significant <- Significant_filter(regiongeneall_Qvalue)
  
  # visualization #
  Group_boxplot(regiongeneall)
  Manhattan_plot(regiongeneall_Qvalue, chrlabs = c(1:19,"X","Y"), col = c("red","blue"), genomewideline = -log10(1e-02))
  Chromosome_pieplot(regiongeneall_significant, title = "Distribution of DMG")
  Heatmap_plot(regiongeneall_significant, title = "Methylation level (%) for genes", display_numbers = FALSE)
  
  Enrich_plot(regiongeneall_significant, enrichterm = "GOgroup", Dbannotation = Dbannotation, title = "Biological process for significant gene")
  Enrich_plot(regiongeneall_significant, enrichterm = "GO", Dbannotation = Dbannotation, title = "Go term for significant gene") 
  Enrich_plot(regiongeneall_significant, enrichterm = "pathway", keggorganism = keggorganism, title = "Pathway for significant gene")
  
  return(regiongeneall_significant)
  
  print(paste("The finish time is", date(), sep = " "))
}




 
#' Quick use the GeneDMRs package for differentially methylated cytosine sites
#' 
#' @import dplyr
#' @import qqman
#' 
#' @description This function outputs the differentially methylated cytosine sites (DMCs).
#' 
#' @param paths refers to the path of input file, with default the package path.
#' @param suffixmeth refers to the suffix of methylation file, e.g. ".gz", ".zip" and so on (some files are in text .txt format, then ".txt" or ".txt.gz"), with default ".gz".
#' 
#' @return Outputs DMC results.
#'
#' @examples
#' Quick_DMCs()
#' 
#' @export


Quick_DMCs <- function(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""), suffixmeth = ".gz"){
  
  # read the file #
  inputmethfile <- Methfile_read(paths = paths, suffix = suffixmeth)
  
  # quality control #
  inputmethfile_QC <- Methfile_QC(inputmethfile)
  
  # methylation mean #
  siteall <- Methmean_site(inputmethfile_QC)
  
  # statistical test #
  siteall_Qvalue <- Logic_regression(siteall)
  
  # sifnificant filter #
  siteall_significant <- Significant_filter(siteall_Qvalue, methdiff = 0.1)
  
  # visualization #
  Manhattan_plot(siteall_Qvalue, chrlabs = c(1:19,"X","Y"), col = c("red","blue"), genomewideline = -log10(1e-02))
  Volcano_plot(siteall_Qvalue)
  Chromosome_pieplot(siteall_significant, title = "Distribution of DMC")
  
  return(siteall_significant)
  print(paste("The finish time is", date(), sep = " "))
}




