# GeneDMRs

# Description
GeneDMRs is an R package to detect the differentially methylated regions based on genes (DMG), gene body (DMP, DME, DMI), CpG islands and gene body interacted with CpG island features (e.g. DMG/DMP/DME/DMI_CpG island and DMG/DMP/DME/DMI_CpG island shore). 

# Dependencies
# Use the annotation dataset for enrichment, e.g. "org.Mm.eg.db" of mouse
if (!requireNamespace("BiocManager", quietly = TRUE))

    install.packages("BiocManager")
    
  BiocManager::install(c("devtools", "clusterProfiler", "corrplot", "dplyr", "genomation", "Hmisc", "KEGG.db", 
                         "pheatmap", "plotrix", "qqman", "RColorBrewer", "RCircos", "VennDiagram", "org.Mm.eg.db"))

# Installation
source("https://install-github.me/xiaowangCN/GeneDMRs")

# Use the annotation dataset for enrichment, e.g. "org.Mm.eg.db" of mouse
Prepare_GeneDMRs(Dbannotation = "org.Mm.eg.db")

# User manual
See the GeneDMRs.pdf file

# examples
Quick_GeneDMRs()

Quick_DMCs()

# Contact
Xiao Wang, xiwa@dtu.dk, Technical University of Denmark
