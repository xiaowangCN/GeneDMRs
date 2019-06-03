#' Quality control for the input methylation file.
#' 
#' @import dplyr
#' 
#' @description This function discards the cytosine sites with low read coverage (quantile) or high read coverage (quantile).
#'
#' @param inputmethfile refers to the input of methylation file after Methfile_read().
#' @param low_coveragenum refers to the minimum read coverage to be discarded, with default 10.
#' @param high_coveragenum refers to the maximum read coverage to be discarded, with default NULL.
#' @param low_quantile refers to the minimum quantile of read coverage to be discarded, with default NULL.
#' @param high_quantile refers to the maximum quantile of read coverage to be discarded, with default 99.99.
#' @param coveragewrite refers to whether to write out the read coverage file to the given path, with default TRUE.
#'
#' @return Outputs a data frame contain chromosome, position, and Cs & Ts for different replicates and groups after quality control.
#' 
#' @examples 
#' inputmethfile_QC <- Methfile_QC(inputmethfile)
#' inputmethfile_QC <- Methfile_QC(inputmethfile, low_coveragenum = 20, high_quantile = 99.99)
#' inputmethfile_QC <- Methfile_QC(inputmethfile, low_coveragenum = 10, high_coveragenum = 100, coveragewrite = FALSE)
#' 
#' @export


# So, what are the inputmethfile_QC file looking like? #

 # The output file of this function is with header and combines all the files with the same chromosome and position #
 # It contains choromsome, position, methylated read coverages (Cs) and unmethylated read coverages (Ts) of each sample.
    #chr  posi Cs1_1 Ts1_1 Cs1_2 Ts1_2 Cs1_3 Ts1_3 Cs2_1 Ts2_1 Cs2_2 Ts2_2 Cs2_3 Ts2_3 Cs3_1 Ts3_1 Cs3_2 Cs3_2 Ts3_2 Cs3_3 Ts3_3
 #1  chr1 21897     0    10     1   107     1    43     0    71     0    47     1    69     1    88    0    64     0    19
 #2  chr1 21905     0    10     0   108     0    44     0    71     1    46     0    70     1    88    1    63     0    19
 #3  chr1 21913     0    10     0   106     0    44     0    72     0    47     0    70     0    89    0    64     0    19
 #4  chr1 21917     0    10     1   105     0    44     1    70     0    46     0    70     0    89    0    64     0    19
 #5  chr1 21926     0    10     0   108     0    44     0    72     0    46     0    70     0    89    0    64     0    19
 #6  chr1 21933     0    10     1   107     0    44     1    72     0    46     0    70     0    89    0    64     0    19
 #7  chr1 21935     0    10     0   108     0    44     1    71     0    46     0    69     1    88    0    64     0    19
 #8  chr1 21937     0    10     1   107     0    43     0    73     1    45     0    70     1    88    0    64     0    19
 #9  chr1 21955     0    10     2   105     1    43     3    70     0    47     0    70     2    87    0    64     0    19
 #10 chr1 21957     0    10     2   104     0    44     0    73     1    46     0    70     0    89    1    63     0    19


# QC for inputmethfile #
Methfile_QC <- function(inputmethfile, low_coveragenum = 10, high_coveragenum = NULL, low_quantile = NULL, high_quantile = 99.9, coveragewrite = TRUE){
  
  # sample number #
  samplenum <- (ncol(inputmethfile) - 2) / 2
  
  # Total coverage = methylated + unmethylated #
  coveragetotal <- array(0, c(nrow(inputmethfile), (1 + samplenum)))
  coveragetotal <-  data.frame(inputmethfile[,c(1:2)], coveragetotal)
  
  j = 3
  for(i in 3:(2 + samplenum)){
    coveragetotal[,i] <- inputmethfile[,j] + inputmethfile[,(j+1)]
    j = j + 2
  }
  names(coveragetotal)[3:ncol(coveragetotal)] <- c(paste("coverage", 1:samplenum, sep = ""), "unQCnum")
  
  # set coveragenum by low_quantile or high_quantile for each sample as quantile before coveragenum#
  if(is.null(low_quantile) == FALSE){
  
    # set tmpnum #
    tmpnum = min(coveragetotal[,3:(2 + samplenum)])
    for(p in 3:(2 + samplenum)){
      quantilenum <- quantile(coveragetotal[,p], (low_quantile / 100))
	  if(quantilenum >= tmpnum){
	  tmpnum <- quantilenum
	  }
    }
    low_coveragenum <- tmpnum
  }  
  
  if(is.null(high_quantile) == FALSE){
  
    # set tmpnum #
    tmpnum = max(coveragetotal[,3:(2 + samplenum)])
    for(q in 3:(2 + samplenum)){
      quantilenum <- quantile(coveragetotal[,q], (high_quantile / 100))
	  if(quantilenum <= tmpnum){
	  tmpnum <- quantilenum
	  }
    }
    high_coveragenum <- tmpnum
  }  

  # coverage number less than coveragenum, default 10 and then filtered #
  if(is.null(high_coveragenum) == TRUE){
    unQC <- coveragetotal[,3:(2 + samplenum)] < low_coveragenum
  }else{
    unQC <- coveragetotal[,3:(2 + samplenum)] < low_coveragenum | coveragetotal[,3:(2 + samplenum)] > high_coveragenum
  }
  
  # transform TRUE and FALSE to number 1 and 0 #
  unQC[unQC=="TRUE"] <- 1
  unQC[unQC=="FALSE"] <- 0
  coveragetotal$unQCnum <- apply(unQC, 1, sum)

  # print to the system file #
  if(coveragewrite == TRUE){
    write.table(coveragetotal, paste(system.file(package = "GeneDMRs"), "/methdata/Coverage", sep = ""), col.names = T, row.names = F, quote = F)
    }
  
    # filter the CpG sites with coverage number > 10 #
    tmp <- data.frame(inputmethfile, unQCnum = coveragetotal$unQCnum)
    inputmethfile <- filter(tmp, unQCnum == 0)
    
    return(inputmethfile[,-ncol(inputmethfile)])
}





#' Merge the methylation file after quality control with DMCs 
#' 
#' @description This function merges the methylation file after quality control of all samples with the DMCs after Significant_filter().
#'
#' @param inputmethfile_QC refers to the input methylation file after quality control.
#' @param siteall_significant refers to the input DMCs file.
#'
#' @return Outputs a data frame by merging two input files of inputmethfile_QC and siteall_significant.
#' 
#' @examples 
#' DMC_inputmethfile_QC <- DMC_methfile_QC(inputmethfile_QC, siteall_significant)
#' 
#' @export


DMC_methfile_QC <- function(inputmethfile_QC, siteall_significant){
  
  # unique id by chromosome and position#
  idsig <- data.frame(id = paste(siteall_significant$posi, siteall_significant$chr, sep = ""))
  tmpinputmeth <- data.frame(id = paste(inputmethfile_QC$posi, inputmethfile_QC$chr, sep = ""), inputmethfile_QC)
  output <- merge(x = idsig, y = tmpinputmeth, by = "id")
  output <- output[order(output$chr, output$posi), ]
  rownames(output) <- 1:nrow(output)
  
  return(output[, -1])
}




