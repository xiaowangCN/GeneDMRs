#' Gene-based differentially methylated regions analysis (GeneDMRs) and install the dependencies.
#' 
#' @description GeneDMRs is an R package to detect the differentially methylated regions based on genes (DMG), gene body (DMP, DME, DMI), 
#' CpG islands and gene body interacted with CpG island features (e.g., DMG/DMP/DME/DMI_CpG island and DMG/DMP/DME/DMI_CpG island shore). 
#' This function can install the other R packages for the dependencies of GeneDMRs.
#' 
#' @param Dbannotation refers to the annotation dataset for enrichment, with default "org.Mm.eg.db" of mouse.
#' 
#' @return Outputs a list of required R packages.
#'
#' @examples
#' GeneDMRs(Dbannotation = "org.Mm.eg.db")
#' 
#' @export


GeneDMRs <- function(Dbannotation = "org.Mm.eg.db"){
  
  # install the dependencies and library them #
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install(c("devtools", "clusterProfiler", "corrplot", "dplyr", "ffbase", "genomation", "KEGG.db", 
                         "pheatmap", "plotrix", "qqman", "RCircos", "VennDiagram", Dbannotation))
  
  require("clusterProfiler")
  require("corrplot")
  require("devtools")
  require("dplyr")
  require("ffbase")
  require("genomation")
  require("KEGG.db")
  require("pheatmap")
  require("plotrix")
  require("qqman")
  require("RCircos")
  require("VennDiagram")
  require(Dbannotation)
}





#' Quick use the GeneDMRs package for gene based differentially methlated regions.
#' 
#' @import dplyr
#' @import ffbase
#' @import qqman
#' @import pheatmap
#' @import clusterProfiler
#' 
#' @description This function outputs a series of results and figures for gene based regions' methylation analysis.
#' 
#' @param paths refers to the path of input file, with default the package path.
#' @param control_paths refers to the path of control groups, with default NULL.
#' @param case_paths refers to the path of case groups, with default NULL.
#' @param suffixmeth refers to the suffix of methylation file, e.g., ".gz", ".zip" and so on (some files are in text .txt format, then ".txt" or ".txt.gz"), with default ".gz".
#' @param bedfile refers to the file name of bed file for "refseq". This file is downloaded from UCSC website, with default "refseq".
#' @param suffixbed refers to the suffix of bed file, e.g., ".gz", ".zip" and so on (some files are in text .txt format, then ".txt" or ".txt.gz"), with default ".txt".
#' @param Dbannotation refers to the annotation dataset for enrichment, with default "org.Mm.eg.db" of mouse.
#' @param keggorganism refers to the species name for KEGG enrichment, with default "mmu" of mouse.
#' 
#' @return Outputs a series of DMG results.
#'
#' @examples
#' allDMGs <- Quick_GeneDMRs()
#' allDMGs_mouse <- Quick_GeneDMRs(Dbannotation = "org.Mm.eg.db", keggorganism = "mmu")
#' 
#' # if only case and control group (n = 2) paths are provided #
#' controls <- c("C:/Users/GeneDMRs/methdata/1_1.gz", "C:/Users/GeneDMRs/methdata/1_2.gz", "C:/Users/GeneDMRs/methdata/1_3.gz")
#' cases <- c("C:/Users/GeneDMRs/methdata/2_1.gz", "C:/Users/GeneDMRs/methdata/2_1.gz")
#' allDMGs <- Quick_GeneDMRs(paths = "C:/Users/GeneDMRs/methdata", control_paths = controls, case_paths = cases)
#' 
#' @export


Quick_GeneDMRs <- function(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""), suffixmeth = ".gz",
                           control_paths = NULL, case_paths = NULL, bedfile = "refseq", 
                           suffixbed = ".txt", Dbannotation = "org.Mm.eg.db", keggorganism = "mmu"){
  
  # read the file #
  # if case and control group paths are provided #
  if(is.null(control_paths) == F & is.null(case_paths) == F){
    inputmethfile <- Methfile_read(control_paths = controls, case_paths = cases)
    
  }else{
    inputmethfile <- Methfile_read(paths = paths, suffix = suffixmeth)
  }
  
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




 
#' Quick use the GeneDMRs package for differentially methylated cytosine sites.
#' 
#' @import dplyr
#' @import qqman
#' 
#' @description This function outputs the differentially methylated cytosine sites (DMCs).
#' 
#' @param paths refers to the path of input file, with default the package path.
#' @param control_paths refers to the path of control groups, with default NULL.
#' @param case_paths refers to the path of case groups, with default NULL.
#' @param suffixmeth refers to the suffix of methylation file, e.g., ".gz", ".zip" and so on (some files are in text .txt format, then ".txt" or ".txt.gz"), with default ".gz".
#' 
#' @return Outputs DMC results.
#'
#' @examples
#' allDMCs <- Quick_DMCs()
#' 
#' # if only case and control group (n = 2) paths are provided #
#' controls <- c("C:/Users/GeneDMRs/methdata/1_1.gz", "C:/Users/GeneDMRs/methdata/1_2.gz", "C:/Users/GeneDMRs/methdata/1_3.gz")
#' cases <- c("C:/Users/GeneDMRs/methdata/2_1.gz", "C:/Users/GeneDMRs/methdata/2_1.gz")
#' allDMCs <- Quick_DMCs(control_paths = controls, case_paths = cases)
#' 
#' @export


Quick_DMCs <- function(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""), suffixmeth = ".gz",
                       control_paths = NULL, case_paths = NULL){
  
  # read the file #
  # if case and control group paths are provided #
  if(is.null(control_paths) == F & is.null(case_paths) == F){
    inputmethfile <- Methfile_read(control_paths = controls, case_paths = cases)
    
  }else{
    inputmethfile <- Methfile_read(paths = paths, suffix = suffixmeth)
  }
  
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




