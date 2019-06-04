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
    
  BiocManager::install(c("devtools", "clusterProfiler", "corrplot", "dplyr", "genomation", "Hmisc", "KEGG.db", 
                         "pheatmap", "plotrix", "qqman", "RColorBrewer", "RCircos", "VennDiagram", "org.Mm.eg.db"))
```

### Installation

```
source("https://install-github.me/xiaowangCN/GeneDMRs")
```

### User manual

See the GeneDMRs.pdf file

### Examples

If dependencies need to be installed

```
GeneDMRs(Dbannotation = "org.Mm.eg.db")
```

If get all differentially methylated genes (DMGs) quickly

```
allDMGs <- Quick_GeneDMRs()
```

If get all differentially methylated cytosine sites (DMCs) quickly

```
allDMCs <- Quick_DMCs()
```

### More about GeneDMRs

```
?GeneDMRs
```

## Author and contact information

Xiao Wang, Technical University of Denmark, xiwa@dtu.dk or wangxiao880923@gmail.com

## Reference

Xiao Wang and Haja N. Kadarmideen. GeneDMRs: an R package for Gene-based Differentially Methylated Regions analysis. International Society for Computational Biology (ISCB), 2019.
