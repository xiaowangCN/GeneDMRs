#' Manhattan plot for all cytosines or regions.
#' 
#' @import qqman
#' 
#' @description This function outputs the Manhattan plot for all cytosines or regions in different chromosomes with significant line based on R package qqman.
#' 
#' @param siteall_Qvalue refers to the input file with Q value from DMR_test(), e.g., regiongeneall_Qvalues, genefeatureall_cpgfeature_Qvalue or others with Q values.
#' @param chrlabs refers to the label of chromosomes, with default NULL.
#' @param col refers to the color of plots, with default black and grey.
#' @param suggestiveline refers to the significant line, with default 0.01.
#' @param genomewideline refers to the genome-wide significant line, with default 0.001.
#' 
#' @return Outputs a Manhattan figure with Q values.
#'
#' @examples
#' Manhattan_plot(siteall_Qvalue, ylab = "-log(Q-value)")
#' Manhattan_plot(regiongenealls_Qvalue, chrlabs = c(1:18,"X"), col = c("green","orange"), genomewideline = -log10(1e-02))
#' Manhattan_plot(genefeatureall_cpgfeature_Qvalue, ylab = c("-log(Q value) for CpG island", "-log(Q value) for Shore"), col = c("red","blue"),
#' suggestiveline = -log10(5e-02), genomewideline = -log10(1e-02))                        
#' 
#' @export


Manhattan_plot <- function(siteall_Qvalue, chrlabs = NULL, col = c("black", "grey"), ylab = "-log(Q value)",
                           suggestiveline = -log10(1e-02), genomewideline = -log10(1e-03)){
  
  # transfer chr column to character #
  siteall_Qvalue$chr <- as.vector(unlist(siteall_Qvalue$chr))
  
  # find the unannotated chromosome rows and delete them #
  unqualifiedrow <- grep("_", siteall_Qvalue$chr)
  
  if(length(unqualifiedrow) > 0){
    siteall_Qvalue <- siteall_Qvalue[-unqualifiedrow,]
  }
  
  # set the chromosome label to number as the default of qqman #
  chromtable <- table(siteall_Qvalue$chr)
  
  # set the chromosome label to number as the default of qqman #
  for(i in 1:length(chromtable[chromtable != 0])){
    
    siteall_Qvalue$chr[siteall_Qvalue$chr == paste("chr", i ,sep="")] <- i
  }
  siteall_Qvalue$chr[siteall_Qvalue$chr == "chrX"] <- length(chromtable[chromtable != 0]) - 1
  siteall_Qvalue$chr[siteall_Qvalue$chr == "chrY"] <- length(chromtable[chromtable != 0])
  siteall_Qvalue$chr <- as.numeric(siteall_Qvalue$chr)
  
  # check the site position or regions start or more feature qvalues #
  if(length(grep("posi", colnames(siteall_Qvalue))) > 0){
    
    subtmp <- data.frame(CHR = siteall_Qvalue$chr, BP = siteall_Qvalue$posi, P = siteall_Qvalue$Qvalue1, 
                         SNP = paste("rs", 1: nrow(siteall_Qvalue), sep = ""))
    
    # filter the unvailable Q values #					   
    subtmp <- subtmp[subtmp$P != "NaN", ]	
    subtmp <- subtmp[order(subtmp$CHR), ] 
    manhattan(subtmp, chr = "CHR", bp = "BP", p = "P", snp = "SNP", chrlabs = chrlabs, col = col, ylab = ylab[1],
              suggestiveline = suggestiveline, genomewideline = genomewideline)
    
  }else if(length(grep("Qvalue", colnames(siteall_Qvalue))) > 1){
    
    featurenum <- length(grep("Qvalue", colnames(siteall_Qvalue)))
    featurepos <- grep("Qvalue", colnames(siteall_Qvalue))
    
    # set the subsequent figures # 
    op <- par(mfrow=c(featurenum, 1))
    
    for(i in 1:featurenum){
      
      # for the position of region, the middle point will be used #
      subtmp <- data.frame(CHR = siteall_Qvalue$chr, BP = 0.5*(siteall_Qvalue$start + siteall_Qvalue$end),
                           P = siteall_Qvalue[, featurepos[i]], SNP = siteall_Qvalue$id)
      
      # filter the unvailable Q values #					   
      subtmp <- subtmp[subtmp$P != "NaN", ]
      subtmp <- subtmp[order(subtmp$CHR), ] 
      manhattan(subtmp, chr = "CHR", bp = "BP", p = "P", snp = "SNP", chrlabs = chrlabs, col = col, ylab = ylab[i],
                suggestiveline = suggestiveline, genomewideline = genomewideline)
    }
    
    # At end of plotting, reset to previous settings of par #
    par(op)
    
  }else{
    
    # for the position of region, the middle point will be used #
    subtmp <- data.frame(CHR = siteall_Qvalue$chr, BP = 0.5*(siteall_Qvalue$start + siteall_Qvalue$end),
                         P = siteall_Qvalue$Qvalue1, SNP = siteall_Qvalue$id)
    
    # filter the unvailable Q values #					   
    subtmp <- subtmp[subtmp$P != "NaN", ]	
    subtmp <- subtmp[order(subtmp$CHR), ] 
    manhattan(subtmp, chr = "CHR", bp = "BP", p = "P",snp = "SNP", chrlabs = chrlabs, col = col, ylab = ylab[1],
              suggestiveline = suggestiveline, genomewideline = genomewideline)
  }
}




 
#' Volcano plot for all the cytosines.
#' 
#' @description This function outputs the volcano plot for all the cytosines with Q values and methylation differences.
#' 
#' @param siteall_Qvalue refers to the input file with Q values and methylation differences.
#' @param title refers to the figure title, with default "Volcano for Q value and methylation difference".
#' @param qvalue refers to the threshold of Q values that Q values less than this will be colored, with default 0.01.
#' @param methdiffpercentage refers to the threshold of methylation level (%) differences that 
#' methylation differences larger than this will be colored, with default 5, 10, 15, 20, 25.
#' 
#' @param pointcolor refers to the point plot color, with default "red", "purple", "orange", "yellow", "blue", "green".
#' 
#' @return Outputs a volcano figure.
#'
#' @examples
#' Volcano_plot(siteall_Qvalue)
#' Volcano_plot(siteall_Qvalue, pointcolor = c("red", "blue", "yellow", "purple", "orange", "green"))
#' Volcano_plot(siteall_Qvalue, title = "Volcano plot", qvalue = 0.001, methdiffpercentage = c(10, 15, 20, 30, 40), 
#' pointcolor = c("red", "purple", "orange", "yellow", "blue", "green"))
#' 
#' @export


Volcano_plot <- function(siteall_Qvalue, title = "Volcano for Q value and methylation difference", qvalue = 0.01,
                         methdiffpercentage = c(5, 10, 15, 20, 25), pointcolor = c("red", "purple", "orange", "yellow", "blue", "green")){
  
  siteall_Qvalue$Methdiff1 <- siteall_Qvalue$Methdiff1*100
  with(siteall_Qvalue, plot(Methdiff1, -log10(Qvalue1), pch = 20, main = title, xlim = c(-100, 100), 
                 xlab = "Methylation (%) difference", ylab = "-log10(Qvalue)"))
  
  # Add colored points: red if q < 0.01 #
  with(subset(siteall_Qvalue, Qvalue1 < qvalue ), points(Methdiff1, -log10(Qvalue1), pch = 20, col = pointcolor[1]))
  
  # purple if q < 0.01 & Methdiff1 > 5 #
  with(subset(siteall_Qvalue, Qvalue1 < qvalue & abs(Methdiff1) > methdiffpercentage[1]), points(Methdiff1, -log10(Qvalue1), pch = 20, col = pointcolor[2]))
  
  # orange if q < 0.01 & Methdiff1 > 10 #
  with(subset(siteall_Qvalue, Qvalue1 < qvalue & abs(Methdiff1) > methdiffpercentage[2]), points(Methdiff1, -log10(Qvalue1), pch = 20, col = pointcolor[3]))
  
  # yellow if q < 0.01 & Methdiff1 > 15 #
  with(subset(siteall_Qvalue, Qvalue1 < qvalue & abs(Methdiff1) > methdiffpercentage[3]), points(Methdiff1, -log10(Qvalue1), pch = 20, col = pointcolor[4]))
  
  # blue if q < 0.01 & Methdiff1 > 20 #
  with(subset(siteall_Qvalue, Qvalue1 < qvalue & abs(Methdiff1) > methdiffpercentage[4]), points(Methdiff1, -log10(Qvalue1), pch = 20, col = pointcolor[5]))
  
  # green if q < 0.01 & Methdiff1 > 25 #
  with(subset(siteall_Qvalue, Qvalue1 < qvalue & abs(Methdiff1) > methdiffpercentage[5]), points(Methdiff1, -log10(Qvalue1), pch = 20, col = pointcolor[6]))
}





 