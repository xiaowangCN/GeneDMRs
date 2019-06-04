#' Internal Use Function
#' That transforms data for boxplot.
#'
#' @description This function transforms filtered subset of inputmethfile to a dataframe for boxplot.
#' 
#' @param methallchr refers to the filtered subset of inputmethfile.
#'
#' @return Outputs dataframe.
#'  
#' @examples
#' Boxplot_trans(methallchr)
#' 
#' @export


Boxplot_trans <- function(methallchr){
  groupfile <- NULL
  
  #calculate the total group number from input file
  groupnum <- length(grep("_1",colnames(methallchr))) / 2
  
  for(i in 1:groupnum){
    realreplicatenum <- length(grep(paste(i, "_", sep = ""), colnames(methallchr))) / 2
    for(j in 1:realreplicatenum){
      colnum <- grep(paste("Cs", i, "_", j, sep = ""), colnames(methallchr))
      framefile <- data.frame(methallchr[,c(colnum, colnum + 1)], 0, 0, paste(i, "_", j, sep = ""))
      colnames(framefile) <-c( "Cs","Ts", "Read", "Meth", "Sample")
      groupfile <- rbind(groupfile, framefile) 
      colnames(groupfile) <-c( "Cs","Ts", "Read", "Meth", "Sample")
    }
  }
  
  # Calculate the read coverage and methylation level #
  groupfile$Read <- groupfile$Cs + groupfile$Ts
  groupfile$Meth <- groupfile$Cs / groupfile$Read
  
  return(groupfile)
}





#' Boxplot the methylation levels or read numbers in different samples.
#'
#' @import dplyr
#' 
#' @description This function outputs the methylation levels or read numbers of the selected genes or all the genes in the different samples.
#' 
#' @param inputmethfile refers to the input file with methylation levels.
#' @param inputrefseqfile refers to the input of gene regions.
#' @param Meth_plot refers to whether to plot the methylation levels, with default TRUE, otherwise to plot the read numbers.
#' @param ylab refers to the label of y axis, with default "Methylation level".
#' @param refseqname refers to NCBI ID of specific gene, with default NULL.
#' @param col refers to boxplot colors, with default NULL.
#' 
#' @return Outputs a boxplot figure with all the samples.
#'
#' @examples
#' Sample_boxplot(inputmethfile_QC, inputrefseqfile)
#' Sample_boxplot(inputmethfile_QC, inputrefseqfile, refseqname = "NM_001244864")
#' Sample_boxplot(inputmethfile_QC, inputrefseqfile, refseqname = c("NM_001244864", "NM_001244534"))
#' Sample_boxplot(inputmethfile_QC, inputrefseqfile, , ylab = "Methylation level (%)"， refseqname = c("NM_001244864", "NM_001143697", "NM_213902"), col = c("red", "green", "blue"))
#' Sample_boxplot(inputmethfile_QC, inputrefseqfile, Meth_plot = FALSE, ylab = "Read number", col = c("red", "blue"))
#' Sample_boxplot(inputmethfile_QC, inputrefseqfile, Meth_plot = FALSE, ylab = "Read number", refseqname = c("NM_001244864", "NM_001244534"))
#' Sample_boxplot(inputmethfile_QC, inputrefseqfile, Meth_plot = FALSE, ylab = "Read number", 
#' refseqname = c("NM_001244864", "NM_001143697", "NM_213902"), col = c("red", "green", "blue"))
#' 
#' @export


Sample_boxplot <- function(inputmethfile, inputrefseqfile, Meth_plot = TRUE, ylab = "Methylation level", refseqname = NULL, col = NULL){
  
  # total methylation level #
  if(is.null(refseqname) == TRUE){
    
    # transforms a dataframe for boxplot #
    groupfile <- Boxplot_trans(inputmethfile)
    
    if(Meth_plot == TRUE){
      boxplot(Meth*100 ~ Sample, groupfile, ylab = ylab, main = "Overall methylation level", col = col)
    }else{
      boxplot(Read ~ Sample, groupfile, ylab = ylab, main = "Overall read number", col = col)
    }
    
    #check how many gene names
  }else if(length(refseqname)==1){
    
    #check a list of genes from regiongeneall file
    geneall <- filter(inputrefseqfile, refseq %in% refseqname)
    
    methallchr <- filter(inputmethfile, chr %in% geneall$chr, posi >= geneall$start, posi <= geneall$end)
    groupfile <- Boxplot_trans(methallchr)
    
    if(Meth_plot == TRUE){
      boxplot(Meth*100 ~ Sample, groupfile, ylab = ylab, main = refseqname, col = col)
    }else{
      boxplot(Read ~ Sample, groupfile, ylab = ylab, main = refseqname, col = col)
    }
    
  }else{
    
    #check a list of genes from regiongeneall file
    geneall <- filter(inputrefseqfile, refseq %in% refseqname)
    
    # set the subsequent figures # 
    op <- par(mfrow=c(length(refseqname),1))
    
    for(i in 1:length(refseqname)){
      genealli <- geneall[i,]
      methallchr <- filter(inputmethfile, chr %in% genealli$chr, posi >= genealli$start, posi <= genealli$end)
      groupfile <- Boxplot_trans(methallchr)
      
      if(Meth_plot == TRUE){
        boxplot(Meth*100 ~ Sample, groupfile, ylab = ylab, main = refseqname[i], col = col)
      }else{
        boxplot(Read ~ Sample, groupfile, ylab = ylab, main = refseqname[i], col = col)
      }
    }
    
    # At end of plotting, reset to previous settings of par #
    par(op)
    
  }
} 

  
  
  
  
#' Boxplot the methylation levels for groups.
#'
#' @description This function outputs the methylation levels of all the groups in boxplot without considering other features.
#' 
#' @param regiongeneall refers to the input file with group methylation levels. 
#' @param ttest refers to whether to perform the Student t-test, with default TRUE.
#' @param title refers to the figure title, with default "Group boxplot among genes".
#' @param col refers to the boxplot colors, with default NULL.
#' 
#' @return Outputs boxplot figure with groups.
#'
#' @examples
#' Group_boxplot(regiongeneall)
#' Group_boxplot(genebodypromoterall, title = "Three groups among genes in promoter region")
#' Group_boxplot(regiongeneall, ttest = FALSE, title = "Three groups among genes", col = c("red", "green", "blue"))
#' 
#' @export


Group_boxplot <- function(regiongeneall, ttest = TRUE, title = "Group boxplot", col = NULL){
  
  #calculate the total group number from input file
  groupnum <- length(grep("Methgroup",colnames(regiongeneall)))
  
  groupfile <- NULL
    
  for(i in 1:groupnum){
    tmpfile <- regiongeneall[,grep("Methgroup",colnames(regiongeneall))[i]]
    framefile <- data.frame(tmpfile, paste("Group", i, sep = ""))
    groupfile <- rbind(groupfile, framefile)
  }
  
  colnames(groupfile) <-c("Methylation","Group")
  
  boxplot(Methylation*100 ~ Group, groupfile, axes = FALSE, main = title, col = col)
  axis(1, at = 1:groupnum, labels = paste("Group", 1:groupnum, sep = ""),font = 2)
  axis(2, las=1)
  title(ylab = "Methylation (%)", cex.lab = 1.5)
  
  # student t-test #
  if (ttest == TRUE){
    num = 1
    for(m in 1:(groupnum - 1)){
      for(n in (m + 1):groupnum){
        pvalue <- round(t.test(regiongeneall[,grep("Methgroup", colnames(regiongeneall))[m]], regiongeneall[,grep("Methgroup",colnames(regiongeneall))[n]])$p.value, 2)
        
        # when pvalue is right #
        if(!is.null(pvalue)){
          text(1.5, 95, paste("P value (Group ", m, " vs ", n, ") = ", pvalue, sep = ""), pos = 1, offset = num)
          num <- num + 1
        }
      }
    }
  }
}





#' Internal Use Function.
#' That boxplots the methylation levels for one group based on CpG island features.
#'
#' @description This function output the methylation levels in boxplot for one group based on CpG island features,
#' e.g. CpG island and CpG island shore based on groups.
#' 
#' @param genefeatureall_cpgfeature refers to output of Region_match as well. 
#' the output of gene or genebody match with cpgi features, e.g. all genes in promoter or exon or intron regions.
#' 
#' @param groupnum refers to the group number, can be NULL if without groups but the colnames should have "group" name 
#' like CpGislandgroup or Shoregroup, can also be automatically set in the Group_feature_boxplot function.
#' 
#' @param ttest refers to whether perform the Student t-test, with default TRUE.
#' @param cpgfeaturelable refers to CpG island features, with default CpGisland" and "Shore".
#' @param title refers to main text of tilte， with default "Group1".
#' @param col refers to colors, with default NULL.
#' @param yaxislabel refers to the label of y axis, with default TRUE,
#' but the yaxislabel will be set to FALSE when Group_feature_boxplot function is used for only one y axis in multiple groups.
#' 
#' @return Outputs boxplot figure.
#'
#' @examples
#' Cpgfeature_boxplot(genefeatureall_cpgfeature, groupnum = 1, ttest = TRUE, cpgfeaturelable = c("CpGisland", "Shore"), title = "Group1", col = "blue")
#' 
#' @export

Cpgfeature_boxplot <- function(genefeatureall_cpgfeature, groupnum = 1, ttest = TRUE, cpgfeaturelable = c("CpGisland", "Shore"), title = NULL, col = NULL, yaxislabel = TRUE){
  
  # calculate the feature group number #
  cpgfeaturenum <- length(grep("Methgroup1", colnames(genefeatureall_cpgfeature)))
  
  featurefile <- NULL
  
  for(i in 1:cpgfeaturenum){
    tmpfile <- genefeatureall_cpgfeature[, grep(paste("Methgroup", groupnum, sep = ""), colnames(genefeatureall_cpgfeature))[i]]
    framefile <- data.frame(tmpfile, paste("Feature", i, sep = ""))
    featurefile <- rbind(featurefile, framefile)
  }
  
  colnames(featurefile) <-c("Methylation","Feature")
  featurefile <- filter(featurefile, Methylation != "NaN")
  
  boxplot(Methylation*100 ~ Feature, featurefile, axes = FALSE, col = col)
  axis(1, at = 1:cpgfeaturenum, labels = cpgfeaturelable, font = 2)
  if(yaxislabel==TRUE){
    axis(2, las=1)
    title(ylab = "Methylation (%)", cex.lab = 1.5)
  }
  
  # student t-test #
  if (ttest == TRUE){
    
    # if more than two features #
    if(length(cpgfeaturelable) > 2){
      num = 1
      for(m in 1:(cpgfeaturenum - 1)){
        for(n in (m + 1):cpgfeaturenum){
          pvalue <- round(t.test(genefeatureall_cpgfeature[,grep(paste("Methgroup", groupnum, sep = ""), colnames(genefeatureall_cpgfeature))[m]], 
                                 genefeatureall_cpgfeature[,grep(paste("Methgroup", groupnum, sep = ""), colnames(genefeatureall_cpgfeature))[n]])$p.value, 2)
          
          # when pvalue is right #
          if(!is.null(pvalue)){
            plab <- paste("P value ( ", m, " vs ", n, ") = ", pvalue, sep = "")
            text(1.5, y = par()$usr[4]*1.05, xpd=T, label = plab, pos = 1, offset = num, col = col)
            num <- num + 1
          }
        }
      }
      
    }else if(length(cpgfeaturelable) == 2){
      pvalue <- round(t.test(genefeatureall_cpgfeature[,grep(paste("Methgroup", groupnum, sep = ""), colnames(genefeatureall_cpgfeature))[1]], 
                             genefeatureall_cpgfeature[,grep(paste("Methgroup", groupnum, sep = ""), colnames(genefeatureall_cpgfeature))[2]])$p.value, 2)
      
      # when pvalue is right #
      if(!is.null(pvalue)){
        plab <- paste(title, " ( P value = ", pvalue," )", sep = "")
        text(1.5, y = par()$usr[4]*1.05, xpd=T, label = plab, col = col)
      }
    }
  }
}






#' Boxplot the methylation levels for groups based on CpG island features.
#'
#' @description This function outputs the methylation levels in boxplot for one or more groups based on CpG island features,
#' e.g. CpG island and CpG island shore features.
#' 
#' @param genefeatureall_cpgfeature refers to the input file with group methylation levels and CpG island features.
#' @param groupnum refers to the group number, with default "all" for all of the groups.
#' @param ttest refers to whether to perform the Student t-test, with default TRUE.
#' @param cpgfeaturelable refers to CpG island features, with default "CpGisland" and "Shore". Only one CpG island feature can also be available, e.g. "CpGisland".
#' @param title refers to the figure title， with default "Group1", "Group2" and "Group3".
#' @param col refers to the boxplot colors, with default NULL.
#' 
#' @return Outputs a boxplot figure with groups and CpG island features.
#'
#' @examples
#' Group_cpgfeature_boxplot(genefeatureall_cpgfeature, groupnum = 1)
#' Group_cpgfeature_boxplot(genefeatureall_cpgfeature, groupnum = "all", ttest = TRUE, cpgfeaturelable = c("CpGisland", "Shore"), 
#' title = c("Group1", "Group2", "Group3"), col = c("blue", "red", "green"))
#' 
#' @export


Group_cpgfeature_boxplot <- function(genefeatureall_cpgfeature, groupnum = "all", ttest = TRUE, cpgfeaturelable = c("CpGisland", "Shore"), 
                                     title = c("Group1", "Group2", "Group3"), col = NULL){

  if(groupnum == "all"){
      
	  # calculate the total feature and group number from input inputmethfile_QC #
      featurenum <- length(grep("Methgroup1", colnames(genefeatureall_cpgfeature)))
      groupnum <- length(grep("Methgroup",colnames(genefeatureall_cpgfeature))) / featurenum
      yaxislabel = FALSE
	  
    }else{
      yaxislabel = TRUE
    }
  
  # draw the full picture and set the subsequent figures #
  op <- par(mfrow=c(1, groupnum))
  
  for(i in 1:groupnum){
    if(i == 1){
      Cpgfeature_boxplot(genefeatureall_cpgfeature, groupnum = i, ttest = TRUE, cpgfeaturelable, title = title[i], col = col[i], yaxislabel)
      
      # add the y axis #
      axis(2, las=1)
	  title(ylab = "Methylation (%)", cex.lab = 1.5)
      
    }else{
      Cpgfeature_boxplot(genefeatureall_cpgfeature, groupnum = i, ttest = TRUE, cpgfeaturelable, title = title[i], col = col[i], yaxislabel)
    }
  }

  # At end of plotting, reset to previous settings of par #
  par(op)
}
  




#' Boxplot the methylation levels for gene body based on CpG island features.
#'
#' @description This function outputs the methylation levels in boxplot for one or more features of gene body based on CpG island features, 
#' e.g. CpG island and CpG island shore features.
#' 
#' @param genefeatureall_cpgfeature refers to the input file of methylation levels with gene body and CpG island features.
#' @param genebodyname refers to the name of gene body features e.g. promoter, exon, intron and TSSes, with default "promoters", "exons", "introns", "TSSes".
#' @param ttest refers to whether perform the Student t-test, with default TRUE.
#' @param cpgfeaturelable refers to CpG island features, with default "CpGisland" and "Shore".
#' @param title refers to the figure title, with default "Promoter", "Exon", "Intron" and "TSS".
#' @param col refers to the boxplot colors, with default NULL.
#' 
#' @return Outputs a boxplot figure with gene body and CpG island features.
#'
#' @examples
#' Genebody_cpgfeature_boxplot(genefeatureall_cpgfeature)
#' Genebody_cpgfeature_boxplot(genefeatureall_cpgfeature, genebodyname = c("promoters","exons"), 
#' ttest = TRUE, cpgfeaturelable = c("CpGisland", "Shore"), title = c("Promoter", "Exon"), col = c("blue", "red"))
#' 
#' Genebody_cpgfeature_boxplot(genefeatureall_cpgfeature, genebodyname = c("promoters","exons","introns","TSSes"), 
#' ttest = TRUE, cpgfeaturelable = c("CpGisland", "Shore"), title = c("Promoters", "Exons", "Introns", "TSSes"), col = c("blue", "red", "green", "purple"))
#' 
#' @export


Genebody_cpgfeature_boxplot <- function(genefeatureall_cpgfeature, genebodyname = c("promoters","exons","introns","TSSes"), 
                                        ttest = TRUE, cpgfeaturelable = c("CpGisland", "Shore"), title = c("Promoter", "Exon", "Intron", "TSS"), col = NULL){
  
  # calculate the total feature and group number from input inputmethfile_QC #
  featurenum <- length(grep("Methgroup1", colnames(genefeatureall_cpgfeature)))
  groupnum <- length(grep("Methgroup",colnames(genefeatureall_cpgfeature))) / featurenum
  
  # calculate the genebody feature group number #
  genebodyfeaturenum <- length(genebodyname)
  
  # draw the full picture and set the subsequent figures #
  op <- par(mfrow=c(1, genebodyfeaturenum))
  
  for(i in 1:genebodyfeaturenum){
    featurefile <- NULL
    subfile <- filter(genefeatureall_cpgfeature, feature %in% genebodyname[i])
    
    # calculate the feature group number #
    cpgfeaturenum <- length(cpgfeaturelable)
    
    for(j in 1: cpgfeaturenum){
      tmpfile <- subfile[,grep(cpgfeaturelable[j], colnames(subfile))[1:groupnum]]
      
      # unlist the tmfile to make list to one vector #
      framefile <- data.frame(matrix(unlist(tmpfile)), paste("Feature", j, sep = ""))
      featurefile <- rbind(featurefile, framefile)
    }
    colnames(featurefile) <-c("Methylation","Feature")
    featurefile <- filter(featurefile, Methylation != "NaN")
    
    # boxplot #
    boxplot(Methylation*100 ~ Feature, featurefile, axes = FALSE, col = col[i])
    axis(1, at = 1:cpgfeaturenum, labels = cpgfeaturelable, font = 2)
    
    # add one y axis #
    if(i == 1){
      axis(2, las=1)
      title(ylab = "Methylation (%)", cex.lab = 1.5)
    }
    
    # student t-test #
    if (ttest == TRUE){
      
      # if more than two features #
      if(length(cpgfeaturelable) > 2){
        num = 1
        for(m in 1:(cpgfeaturenum - 1)){
          for(n in (m + 1):cpgfeaturenum){
            pvalue <- round(t.test(genefeatureall_cpgfeature[,grep(cpgfeaturelable[m], colnames(genefeatureall_cpgfeature))[1:groupnum]], 
                                   genefeatureall_cpgfeature[,grep(cpgfeaturelable[n], colnames(genefeatureall_cpgfeature))[1:groupnum]])$p.value, 2)
            
            # when pvalue is right #
            if(!is.null(pvalue)){
              plab <- paste("P value ( ", m, " vs ", n, ") = ", pvalue, sep = "")
              text(1.5, 0.95, plab, pos = 1, offset = num, col = col[i])
              num <- num + 1
            }
          }
        }
        
      }else if(length(cpgfeaturelable) == 2){
        pvalue <- round(t.test(genefeatureall_cpgfeature[,grep(cpgfeaturelable[1], colnames(genefeatureall_cpgfeature))[1:groupnum]], 
                               genefeatureall_cpgfeature[,grep(cpgfeaturelable[2], colnames(genefeatureall_cpgfeature))[1:groupnum]])$p.value, 2)
        
        # when pvalue is right #
        if(!is.null(pvalue)){
          plab <- paste(title[i], " ( P value = ", pvalue," )", sep = "")
          text(1.5, y = par()$usr[4]*1.05, xpd=T, label=plab, col = col[i])
        }
      }
    }
  }
  
  # At end of plotting, reset to previous settings of par #
  par(op)
}
 



 