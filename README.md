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
    
  BiocManager::install(c("devtools", "clusterProfiler", "corrplot", "dplyr", "genomation", "KEGG.db", 
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

### Examples

If dependencies need to be installed

```
GeneDMRs(Dbannotation = "org.Mm.eg.db")
```

Before the starting quickly or starting step by step, the user could download the example data or the whole folder from "/methdata" for testing. In the folder "/methdata", "1_1.gz", "1_2.gz" and "1_3.gz" files are the control group, while "2_1.gz" and "2_1.gz" files are case group. The user just needs to give one path for GeneDMRS package, e.g., "paths = paste(system.file(package = "GeneDMRs")" which is the package systme path.

1. If get all differentially methylated genes (DMGs) quickly

```
allDMGs <- Quick_GeneDMRs(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""))
```

2. If get all differentially methylated cytosine sites (DMCs) quickly

```
allDMCs <- Quick_DMCs(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""))
```

3. If get all DMGs step by step

```
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
```

### More about GeneDMRs

```
?GeneDMRs
```

## Author

Xiao Wang, Dan Hao, Haja N. Kadarmideen. Department of Applied Mathematics and Computer Science, Technical University of Denmark

## Maintainer

Xiao Wang. <xiwa@dtu.dk> or <wangxiao880923@gmail.com>

## Reference

Xiao Wang, Dan Hao and Haja N. Kadarmideen. GeneDMRs: an R package for Gene-based Differentially Methylated Regions analysis. F1000Research 2019, 8(ISCB Comm J):1299 (slides) (https://doi.org/10.7490/f1000research.1117223.1)
