# GeneDMRs

Gene-based differentially methylated regions analysis

## Getting Started

### Description

GeneDMRs is an R package to detect the differentially methylated regions based on genes (DMG), gene body (DMP, DME, DMI), CpG islands and gene body interacted with CpG island features (e.g., DMG/DMP/DME/DMI_CpG island and DMG/DMP/DME/DMI_CpG island shore). 

### Dependencies

Use the annotation dataset for enrichment, e.g., "org.Mm.eg.db" of mouse

```
if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")
    
  BiocManager::install(c("devtools", "clusterProfiler", "corrplot", "dplyr", "ffbase", "genomation", 
                         "pheatmap", "plotrix", "qqman", "RCircos", "VennDiagram", "org.Mm.eg.db"))
```

### Installation

```
source("https://install-github.me/xiaowangCN/GeneDMRs")
```

or

```
library("devtools")

install_github("xiaowangCN/GeneDMRs")
```

### User manual

See the GeneDMRs.pdf file

### Sample data

Before starting quickly or starting step by step, the user could download the sample data or the whole folder from "/methdata" for testing (https://github.com/xiaowangCN/GeneDMRs/tree/master/methdata). In the folder "/methdata", "1_1.gz", "1_2.gz" and "1_3.gz" files are the control group, while "2_1.gz" and "2_1.gz" files are case group. The user just needs to give one path for GeneDMRS package, e.g., "paths = paste(system.file(package = "GeneDMRs")" which is the package systme path. 

If the folder is downloaded on the desktop, just run like:
```

allDMGs <- Quick_GeneDMRs(paths = "C:/Users/Desktop/methdata")
```

### Examples

1. If get all differentially methylated genes (DMGs) quickly

```
allDMGs <- Quick_GeneDMRs(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""))
```

Or if it is a case-control design, the user can specify arbitrary file names for case group and control group, separately, where the paths = "C:/Users/GeneDMRs/methdata" here is mainly used for reading bedfile, i.e., inputrefseqfile <- Bedfile_read(paths = "C:/Users/GeneDMRs/methdata"). For example:

```
controls <- c("C:/Users/GeneDMRs/methdata/1_1.gz", "C:/Users/GeneDMRs/methdata/1_2.gz", "C:/Users/GeneDMRs/methdata/1_3.gz")
cases <- c("C:/Users/GeneDMRs/methdata/2_1.gz", "C:/Users/GeneDMRs/methdata/2_1.gz")
allDMGs <- Quick_GeneDMRs(paths = "C:/Users/GeneDMRs/methdata", control_paths = controls, case_paths = cases)
```

2. If get all differentially methylated cytosine sites (DMCs) quickly

```
allDMCs <- Quick_DMCs(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""))
```

Or if it is a case-control design, then the user can specify arbitrary file names for case group and control group, separately, such as:

```
controls <- c("C:/Users/GeneDMRs/methdata/1_1.gz", "C:/Users/GeneDMRs/methdata/1_2.gz", "C:/Users/GeneDMRs/methdata/1_3.gz")
cases <- c("C:/Users/GeneDMRs/methdata/2_1.gz", "C:/Users/GeneDMRs/methdata/2_1.gz")
allDMCs <- Quick_DMCs(control_paths = controls, case_paths = cases)
```

3. If get all DMGs step by step

```
# read the methylation file #
inputmethfile <- Methfile_read(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""), suffix = ".gz")

# or if it is a case-control design #
controls <- c("C:/Users/GeneDMRs/methdata/1_1.gz", "C:/Users/GeneDMRs/methdata/1_2.gz", "C:/Users/GeneDMRs/methdata/1_3.gz")
cases <- c("C:/Users/GeneDMRs/methdata/2_1.gz", "C:/Users/GeneDMRs/methdata/2_1.gz")
inputmethfile <- Methfile_read(control_paths = controls, case_paths = cases)

# read the bedfile #
inputrefseqfile <- Bedfile_read(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""), bedfile = "refseq", suffix = ".txt", feature = FALSE)
  
# quality control #
inputmethfile_QC <- Methfile_QC(inputmethfile)
  
# methylation mean #
regiongeneall <- Methmean_region(inputmethfile_QC, inputrefseqfile, chrnum = "all")
  
# statistical test #
regiongeneall_Qvalue <- Logic_regression(regiongeneall)
  
# sifnificant filter #
regiongeneall_significant <- Significant_filter(regiongeneall_Qvalue)
```

## Author

Xiao Wang, Dan Hao, Haja N. Kadarmideen.

## Maintainer

Xiao Wang. <wangxiao880923@gmail.com>

## Reference

Xiao Wang, Dan Hao and Haja N. Kadarmideen. GeneDMRs: an R package for Gene-based Differentially Methylated Regions analysis. F1000Research 2019, 8(ISCB Comm J):1299 (slides). (https://doi.org/10.7490/f1000research.1117223.1)

Xiao Wang, Dan Hao and Haja N. Kadarmideen. GeneDMRs: an R package for Gene-based Differentially Methylated Regions analysis. bioRxiv 2020, 04.11.037168. (https://doi.org/10.1101/2020.04.11.037168)

Xiao Wang, Dan Hao and Haja N. Kadarmideen. GeneDMRs: an R package for Gene-based Differentially Methylated Regions analysis. Journal of Computational Biology 2021, 28(3), 304-316. (https://doi.org/10.1089/cmb.2020.0081)
