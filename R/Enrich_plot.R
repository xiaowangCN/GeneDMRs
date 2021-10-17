#' Enrich plot for GO terms and pathways.
#'
#' @import clusterProfiler
#' @import org.Mm.eg.db
#'
#' @description This function outputs the groups, GO terms and pathways plot for the enrichment based on R package clusterProfiler and org.Mm.eg.db (Mouse for example).
#'
#' @param regiongenealls_significant refers to the differentially methylated genes.
#' @param adjustpvaluecut refers to the threshold of the adjusted P values for the enrichment, with default 0.1.
#' @param enrichterm refers to the term need to be analyzed, which can be "GOgroup", "GO", "pathway", with default "pathway".
#' @param category refers to whether to divide the enrichments into two categories, i.e., hypo/hyper methylated or down/up regulated, with default TRUE.
#' @param Dbannotation refers to the annotation dataset, with default "org.Mm.eg.db" of mouse.
#' @param keggorganism refers to the species name for KEGG enrichment, with default "mmu" of mouse.
#' @param listnum refers to the list of display number, with default 20.
#' @param title refers to the title of figure, with default "Enrichments for significant gene".
#' @param expressionfile_significant refers to an additional file for differentially expressed genes which includes gene name and Log fold change (LogFC).
#' This is a optional file for GO terms and pathways, with default NULL.
#'
#' @param expressionfile_genetype refers to the gene type of expressionfile_significant file, which can be "REFSEQ", "ENTREZID", "SYMBOL",
#' or other gene types that can be used in clusterProfiler, with default NULL.
#'
#' @return Outputs a dot-plot figure of enrichment.
#'
#' @examples
#' Enrich_plot(regiongenealls_significant, enrichterm = "GOgroup", Dbannotation = "org.Mm.eg.db", title = "Biological process for significant gene")
#' Enrich_plot(regiongenealls_significant, enrichterm = "GO", Dbannotation = "org.Mm.eg.db", title = "Go term for significant gene")
#' Enrich_plot(regiongenealls_significant, enrichterm = "GO", Dbannotation = "org.Hs.eg.db", title = "Go term for significant gene") # for human data #
#' Enrich_plot(regiongenealls_significant, adjustpvaluecut = 0.2, enrichterm = "pathway", Dbannotation = "org.Mm.eg.db", keggorganism = "mmu", title = "Pathway for significant gene")
#' Enrich_plot(regiongenealls_significant, enrichterm = "pathway", keggorganism = "hsa", Dbannotation = "org.Hs.eg.db", category = FALSE, title = "Pathway for significant gene") # for human data #
#'
#' expressionfile_significant <- read.table(paste(system.file(package = "GeneDMRs"), "/methdata/DEgenes.txt", sep=""), header = T) # read DEgene file #
#'
#' Enrich_plot(regiongenealls_significant, adjustpvaluecut = 0.2, enrichterm = "GO", Dbannotation = "org.Mm.eg.db", title = "Go term for significant gene in two categories",
#' expressionfile_significant = expressionfile_significant, expressionfile_genetype = "SYMBOL")
#' Enrich_plot(regiongenealls_significant, enrichterm = "pathway", Dbannotation = "org.Mm.eg.db", keggorganism = "mmu", title = "Pathway for significant gene in two categories",
#' expressionfile_significant = expressionfile_significant, expressionfile_genetype = "SYMBOL")
#'
#' @export


Enrich_plot <- function(regiongenealls_significant, adjustpvaluecut = 0.1, enrichterm = "pathway", category = TRUE, Dbannotation = "org.Mm.eg.db",
                        keggorganism = "mmu", listnum = 20, title = "Enrichment for significant gene",
                        expressionfile_significant = NULL, expressionfile_genetype = NULL){

  # get the gene name #
  genefile <- as.vector(unlist(regiongenealls_significant[1]))
  eg = bitr(genefile, fromType = "REFSEQ", toType= c("ENTREZID", "SYMBOL"), OrgDb = Dbannotation)

  # uniqe ENTREZID of eg #
  eg <- distinct(eg, ENTREZID, .keep_all = TRUE)

  # check enrichterm #
  if(length(grep("GOgroup", enrichterm)) > 0){

    # group GO terms for Biological process (BP) #
    ggo <- groupGO(gene = eg$ENTREZID, OrgDb = Dbannotation, ont = "BP", level = 3, readable = TRUE)
    barplot(ggo, drop = TRUE, showCategory = listnum, title = title)

    # if not "GOgroup" then be "GO" or "pathway" #
  }else{

    # if without gene expression data #
    if(is.null(expressionfile_significant) == TRUE){

      # create new eg for use #
      DMgene_merge <- data.frame(REFSEQ = regiongenealls_significant$id, Methdiff = regiongenealls_significant$Methdiff1)
      eg_merge <- merge(x = eg, y = DMgene_merge, by = "REFSEQ")
      neweg <- data.frame(Entrez = eg_merge$ENTREZID, Methdiff = eg_merge$Methdiff)
      neweg$group <- "Hyper-methylated"
      neweg$group[neweg$Methdiff < 0] <- "Hypo-methylated"

      # check enrichterm #
      if(length(grep("GO", enrichterm)) > 0){
        if(category == TRUE){

          # GO terms between two categories #
          formula_go <- compareCluster(Entrez~group, data = neweg, fun = "enrichGO", ont = "all",
                                       pvalueCutoff = adjustpvaluecut, OrgDb = Dbannotation)
          clusterProfiler::dotplot(formula_go, x = "group", color = "p.adjust", showCategory = listnum, split = NULL, font.size = 14,
                                   title = title)
        }else{

          # GO terms in one category #
          ego <- enrichGO(gene = neweg$Entrez, OrgDb = Dbannotation, ont= "all", pvalueCutoff = adjustpvaluecut)
          clusterProfiler::dotplot(ego, color = "p.adjust", showCategory = listnum, split = NULL, font.size = 12, title = title)
        }

      }else if(length(grep("pathway", enrichterm)) > 0){
        if(category == TRUE){

          # pathways between two categories #
          formula_path <- compareCluster(Entrez~group, data = neweg, fun = "enrichKEGG",
                                         organism = keggorganism, pvalueCutoff = adjustpvaluecut)
          clusterProfiler::dotplot(formula_path, x = "group", color = "p.adjust", showCategory = listnum, split = NULL, font.size = 14,
                                   title = title)
        }else{
          # pathways in one category #
          kk <- enrichKEGG(gene= neweg$Entrez, organism = keggorganism, keyType = "kegg", pvalueCutoff = adjustpvaluecut)
          clusterProfiler::dotplot(kk, color = "p.adjust", showCategory = listnum, split = NULL, font.size = 12, title = title)
        }
      }

    }else{

      # new dataset for DE (expressed) and DM (methylated) genes with +/- #
      names(expressionfile_significant)[1] <- expressionfile_genetype
      DEgene_merge <- merge(x = eg, y = expressionfile_significant, by = expressionfile_genetype)
      DMgene_merge <- data.frame(REFSEQ = regiongenealls_significant$id, Methdiff = regiongenealls_significant$Methdiff1)
      DMDEgene_merge <- merge(x = DEgene_merge, y = DMgene_merge, by = "REFSEQ")

      neweg <- data.frame(Entrez = DMDEgene_merge$ENTREZID, logFC = DMDEgene_merge$LogFC, Methdiff = DMDEgene_merge$Methdiff)
      if(category == TRUE){
        neweg$group <- "Hyper-methylated"
        neweg$group[neweg$Methdiff < 0] <- "Hypo-methylated"
        neweg$othergroup <- "Up-regulated"
        neweg$othergroup[neweg$logFC < 0] <- "Down-regulated"

      }else{
        neweg$group <- "Differentially methylated"
        neweg$othergroup <- "Differentially expressed"
      }

      # check enrichterm #
      if(length(grep("GO", enrichterm)) > 0){

        # GO terms between two categories #
        formula_go <- compareCluster(Entrez~group + othergroup, data = neweg, fun = "enrichGO", ont = "all",
                                     pvalueCutoff = adjustpvaluecut, OrgDb = Dbannotation)
        clusterProfiler::dotplot(formula_go, x = "group", color = "p.adjust", showCategory = listnum, split = NULL, font.size = 14,
                title = title) + ggplot2::facet_grid(~othergroup)

      }else if(length(grep("pathway", enrichterm)) > 0){

        # pathways between two categories #
        formula_path <- compareCluster(Entrez~group + othergroup, data = neweg, fun = "enrichKEGG",
                                       organism = keggorganism, pvalueCutoff = adjustpvaluecut)
        clusterProfiler::dotplot(formula_path, x = "group", color = "p.adjust", showCategory = listnum, split = NULL, font.size = 14,
                title = title) + ggplot2::facet_grid(~othergroup)
      }
    }
  }
}




