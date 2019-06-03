#' Pie plot based on different chromosomes.
#' 
#' @description This function outputs the pie plot for the percentages of sites or regions in different chromosomes.
#' 
#' @param genefeatureall_cpgfeature_significantcpgisland refers to the input file with chromosomes, which can be files with/without significant filter().
#' @param genefeatureall_cpgfeature_significantshore refers to another input file with chromosomes, e.g. genefeatureall_cpgfeature_significantshore file for comparison, with default NULL.
#' @param methdirection refers to the methylation direction when the input file contains the methylation difference column i.e. Methdiff1 after Logic_regression(), 
#'  which can be "hypo", "hyper" and "both", with the default "both" for both directions.
#' 
#' @param title refers to figure titles, with the default "Pie plot for chromosome".
#' 
#' @return Outputs pie figure in different chromosomes.
#'
#' @examples
#' Chromosome_pieplot(genefeatureall_cpgfeature_significantcpgisland, title = "")
#' Chromosome_pieplot(genefeatureall_cpgfeature_significantcpgisland, title = "CpGisland")
#' Chromosome_pieplot(genefeatureall_cpgfeature_significantcpgisland, genefeatureall_cpgfeature_significantshore = genefeatureall_cpgfeature_significantshore, title = c("CpGisland","Shore"))
#' Chromosome_pieplot(siteall, title = "All cytosine sites") # Only consider the annotated chromosomes and the unannotated chromosomes will be discarded #
#' Chromosome_pieplot(siteall_Qvalue, title = "All cytosine sites")
#' Chromosome_pieplot(siteall_significant, title = "Significant cytosine sites")
#' Chromosome_pieplot(siteall_Qvalue, methdirection = "hyper", title = "Hyper-methylated distribution")
#' Chromosome_pieplot(siteall_significant, methdirection = "hypo", title = "Hypo-methylated pie plot")
#' Chromosome_pieplot(regiongeneall_Qvalue, methdirection = "hyper", title = "Hyper-methylated genes")
#' 
#' @export


Chromosome_pieplot <- function(genefeatureall_cpgfeature_significantcpgisland, genefeatureall_cpgfeature_significantshore = NULL, methdirection = "both", title = "Pie plot for chromosome"){
  
  # select the hypo or hyper methylated data #
  if(methdirection == "hypo"){
    genefeatureall_cpgfeature_significantcpgisland <- genefeatureall_cpgfeature_significantcpgisland[genefeatureall_cpgfeature_significantcpgisland$Methdiff1 < 0, ]
    
  }else if(methdirection == "hyper"){
    genefeatureall_cpgfeature_significantcpgisland <- genefeatureall_cpgfeature_significantcpgisland[genefeatureall_cpgfeature_significantcpgisland$Methdiff1 > 0, ]
    
  }else if(methdirection == "both"){
    genefeatureall_cpgfeature_significantcpgisland <- genefeatureall_cpgfeature_significantcpgisland
  }
  
  # find the unannotated chromosome rows and delete them #
  unqualifiedrow1 <- grep("_", genefeatureall_cpgfeature_significantcpgisland$chr)
  
  if(length(unqualifiedrow1) > 0){
    genefeatureall_cpgfeature_significantcpgisland <- genefeatureall_cpgfeature_significantcpgisland[-unqualifiedrow1,]
  }
  
  # calculate the percentage of chromosome #
  chromtable1 <- data.frame(table(genefeatureall_cpgfeature_significantcpgisland$chr))
  chromtable1 <- chromosmoe_sort(chromtable1)
  percen1 <- chromtable1$freq / nrow(genefeatureall_cpgfeature_significantcpgisland)
  lable1 <- paste(chromtable1$chr, " ", round(percen1*100, 1), "%", sep="")
  
  # if only one dataset genefeatureall_cpgfeature_significantcpgisland #
  if(is.null(genefeatureall_cpgfeature_significantshore) == TRUE){
    
    pie(percen1, labels = lable1, col = rainbow(length(lable1)), main = title[1])
    
  }else{
    
    # set the subsequent figures for two files # 
    op <- par(mfrow=c(1,2))
    
    pie(percen1, labels = lable1, col = rainbow(length(lable1)), main = title[1])
    
    # select the hypo or hyper methylated data #
    if(methdirection == "hypo"){
      genefeatureall_cpgfeature_significantshore <- genefeatureall_cpgfeature_significantshore[genefeatureall_cpgfeature_significantshore$Methdiff2 < 0, ]
      
    }else if(methdirection == "hyper"){
      genefeatureall_cpgfeature_significantshore <- genefeatureall_cpgfeature_significantshore[genefeatureall_cpgfeature_significantshore$Methdiff2 > 0, ]
      
    }else if(methdirection == "both"){
      genefeatureall_cpgfeature_significantshore <- genefeatureall_cpgfeature_significantshore
    }
    
    # find the unannotated chromosome rows and delete them #
    unqualifiedrow2 <- grep("_", genefeatureall_cpgfeature_significantshore$chr)
    
    if(length(unqualifiedrow2) > 0){
      genefeatureall_cpgfeature_significantshore <- genefeatureall_cpgfeature_significantshore[-unqualifiedrow2,]
    }
    
    # calculate the percentage of chromosome #
    chromtable2 <- data.frame(table(genefeatureall_cpgfeature_significantshore$chr))
    chromtable2 <- chromosmoe_sort(chromtable2)
    percen2 <- chromtable2$freq / nrow(genefeatureall_cpgfeature_significantshore)
    lable2 <- paste(chromtable2$chr, " ", round(percen2*100, 1), "%", sep="")
    
    pie(percen2, labels = lable2, col = rainbow(length(lable2)), main = title[2])
    
    # At end of plotting, reset to previous settings of par #
    par(op)
  }
}

  
  
  
  
#' Pie plot based on different features.
#' 
#' @import plotrix
#' 
#' @description This function outputs the pie plot of feature percentages in gene body or CpG island mainly for DMC sites with features.
#' 
#' @param siteall_significant_feature refers to the input file with features, mainly for DMC sites with features.
#' @param methdirection refers to the methylation direction when the input file contains the methylation difference column i.e. Methdiff1 after Logic_regression(), 
#'  which can be "hypo", "hyper" and "both", with the default "both" for both directions.
#' 
#' @param title refers to figure titles, with the default "Pie plot for chromosome".
#' @param threeDplot refers to whether to pie plot in three dimensions based on R pacakge plotrix, with the default TRUE.
#' 
#' @return Outputs a pie figure in different features.
#'
#' @examples
#' Feature_pieplot(siteall_significant_feature)
#' Feature_pieplot(siteall_significant_feature, methdirection = "hypo")
#' Feature_pieplot(siteall_significant_feature, title = c("Gene body", "CpG island"))
#' Feature_pieplot(siteall_significant_feature, title = c("Pie plot for Gene body", "Pie plot for CpG island"), threeDplot = FALSE)
#' Feature_pieplot(siteall_significant_feature, methdirection = "hyper", title = c("Pie plot for Gene body", "Pie plot for CpG island"))
#' 
#' @export


Feature_pieplot <- function(siteall_significant_feature, methdirection = "both", title = "Pie plot for feature", threeDplot = TRUE){
  
  # select the hypo or hyper methylated data #
  if(methdirection == "hypo"){
    siteall_significant_feature <- siteall_significant_feature[siteall_significant_feature$Methdiff1 < 0, ]
    
  }else if(methdirection == "hyper"){
    siteall_significant_feature <- siteall_significant_feature[siteall_significant_feature$Methdiff1 > 0, ]
    
  }else if(methdirection == "both"){
    siteall_significant_feature <- siteall_significant_feature
  }
  
  # count the feature number #
  featurenum <- length(grep("Feature", names(siteall_significant_feature)))
  
  # calculate the percentage of first feature #
  chromtable1 <- table(siteall_significant_feature$Feature1)
  
  # check if two features in one site with "&" #
  duppos <- grep("&", names(chromtable1))
  dupnum <- length(duppos)
  
  if(dupnum > 0){
    for(i in 1:dupnum){
      dupname <- names(chromtable1[duppos[i]])
      dupnamenum <- table(strsplit(dupname, " & ")[[1]])
      
      # add the frequency number to the original feature #
      orifeatnum <- length(names(chromtable1[-duppos]))
      for(j in 1:orifeatnum){
        orifeatnumname <- names(chromtable1[-duppos])[j]
        if(length(dupnamenum[names(dupnamenum)==orifeatnumname]) > 0){
          chromtable1[names(chromtable1)==orifeatnumname] <- chromtable1[names(chromtable1)==orifeatnumname] + 
            (dupnamenum[names(dupnamenum)==orifeatnumname] * chromtable1[duppos[i]])
        }
      }
    }
    
    # recalculate the percentage of features that delete the "&" #
    chromtable1_new <- chromtable1[-duppos]
    percen1 <- chromtable1_new / sum(chromtable1_new)
    lable1 <- paste(rownames(chromtable1_new)," ", round(percen1*100, 1), "%", sep="")
    
  }else{
    percen1 <- chromtable1 / sum(chromtable1)
    lable1 <- paste(rownames(chromtable1)," ", round(percen1*100, 1), "%", sep="")
  }
  
  # if only one feature featurenum == "feature_id + feature1"#
  if(featurenum == 2){
    
    if(threeDplot == TRUE){
      pie3D(percen1, explode = 0.1, radius = 0.7, height = 0.07, labelcex = 1, pty= "m",
            labels = lable1, col = rainbow(length(lable1)), main = title[1])
      
    }else{
      pie(percen1, labels = lable1, col = rainbow(length(lable1)), main = title[1])
    }
    
  }else{
    
    # set the subsequent figures for two features # 
    op <- par(mfrow=c(1,2))
    
    if(threeDplot == TRUE){
      pie3D(percen1, explode = 0.1, radius = 0.7, height = 0.07, labelcex = 1, pty= "m",
            labels = lable1, col = rainbow(length(lable1)), main = title[1])
      
    }else{
      pie(percen1, labels = lable1, col = rainbow(length(lable1)), main = title[1])
    }
    
    # calculate the percentage of second feature #
    chromtable2 <- table(siteall_significant_feature$Feature2)
    
    # check if two features in one site with "&" #
    s_duppos <- grep("&", names(chromtable2))
    s_dupnum <- length(s_duppos)
    
    if(s_dupnum > 0){
      for(i in 1:s_dupnum){
        s_dupname <- names(chromtable2[s_duppos[i]])
        s_dupnamenum <- table(strsplit(s_dupname, " & ")[[1]])
        
        # add the frequency number to the original feature #
        s_orifeatnum <- length(names(chromtable2[-s_duppos]))
        for(j in 1:s_orifeatnum){
          s_orifeatnumname <- names(chromtable2[-s_duppos])[j]
          if(length(s_dupnamenum[names(s_dupnamenum)==s_orifeatnumname]) > 0){
            chromtable2[names(chromtable2)==s_orifeatnumname] <- chromtable2[names(chromtable2)==s_orifeatnumname] + 
              (s_dupnamenum[names(s_dupnamenum)==s_orifeatnumname] * chromtable2[s_duppos[i]])
          }
        }
      }
      
      # recalculate the percentage of features that delete the "&" #
      chromtable2_new <- chromtable2[-s_duppos]
      percen2 <- chromtable2_new / sum(chromtable2_new)
      lable2 <- paste(rownames(chromtable2_new)," ", round(percen2*100, 1), "%", sep="")
      
    }else{
      percen2 <- chromtable2 / sum(chromtable2)
      lable2 <- paste(rownames(chromtable2)," ", round(percen2*100, 1), "%", sep="")
    }
    
    if(threeDplot == TRUE){
      pie3D(percen2, explode = 0.1, radius = 0.7, height = 0.07, labelcex = 1, pty= "m",
            labels = lable2, col = rainbow(length(lable2)), main = title[2])
      
    }else{
      pie(percen2, labels = lable2, col = rainbow(length(lable2)), main = title[2])
    }
    
    # At end of plotting, reset to previous settings of par #
    par(op)
  }
}




 
#' Internal Use Function
#' That sorts the chromosomes from "chr1", "chr2", "chr3"......, "chrX" and "chrY".
#'
#' @description This function sorts the chromosomes.
#' 
#' @param chromtable refers to the input file for sorting the chromosomes.
#'
#' @return Outputs the sorted file by chromosome number.
#'  
#' @examples
#' chromtable <- chromosmoe_sort(chromtable)
#' 
#' @export


chromosmoe_sort <- function(chromtable){
  colnames(chromtable) <- c("chr", "freq")
  chromtable <- chromtable[chromtable$freq != 0,]
  
  # transfer chr column to character #
  chromtable$chr <- as.vector(unlist(chromtable$chr))
  
  # set the chromosome label to number #
  for(i in 1:sum(chromtable$freq != 0)){
    
    chromtable$chr[chromtable$chr == paste("chr", i ,sep="")] <- i
  }
  chromtable$chr[chromtable$chr == "chrX"] <- sum(chromtable$freq != 0) - 1
  chromtable$chr[chromtable$chr == "chrY"] <- sum(chromtable$freq != 0)
  
  # transfer chr column to numeric again for sort #
  chromtable$chr <- as.numeric(chromtable$chr)
  chromtable <- chromtable[order(chromtable$chr), ] 
  
  # paste "chr" to chromosome number #
  chromtable$chr <- paste("chr", chromtable$chr, sep = "")
  chromtable$chr[chromtable$chr == paste("chr", (sum(chromtable$freq != 0) - 1), sep = "")] <- "chrX"
  chromtable$chr[chromtable$chr == paste("chr", sum(chromtable$freq != 0), sep = "")] <- "chrY"
  
  return(chromtable)
}




