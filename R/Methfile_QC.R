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
#' @param samplenum_QC refers to the sample numbers under quanlity control (e.g., samplenum_QC = 3 means that 
#'                     if three of five samples at one cytosine site have unqualified read coverage, then this site 
#'                     will be discarded), with default "all" samples.
#'
#' @return Outputs a data frame contain chromosome, position, and Cs & Ts for different replicates and groups after quality control.
#' 
#' @examples 
#' inputmethfile_QC <- Methfile_QC(inputmethfile)
#' inputmethfile_QC <- Methfile_QC(inputmethfile, low_coveragenum = 20, high_quantile = 99.99)
#' inputmethfile_QC <- Methfile_QC(inputmethfile, low_coveragenum = 10, high_coveragenum = 100, samplenum_QC = 3)
#' 
#' @export


# So, what are the inputmethfile_QC file looking like? e.g. data from mouse #

 # The output file of this function is with header and combines all the files with the same chromosome and position #
 # It contains choromsome, position, methylated read coverages (Cs) and unmethylated read coverages (Ts) of each sample.
 #   chr    posi   Cs1_1 Ts1_1 Cs1_2 Ts1_2 Cs1_3 Ts1_3 Cs2_1 Ts2_1 Cs2_2 Ts2_2
 #1  chr1 3020877    77     2    77     7    49     2    31     4    68     0
 #2  chr1 3020891    73     6    78     6    49     2    33     2    68     0
 #3  chr1 3020946    47     6    96    17    71     9    52     5    71    12
 #4  chr1 3020988    73     1    58     0    57     6    55     2    61     2
 #5  chr1 3021013    74     0    56     2    59     4    49     8    63     0
 #6  chr1 3622136    38     3    52     6    42     5    36     0    68     3
 #7  chr1 3622256    93     1   113     8    84     4    59     0   102     4
 #8  chr1 3622257    93     0    88     5    88     5    52     4    95     1
 #9  chr1 3622311    74     3    90     5   104     3    79     9   102     0
 #10 chr1 3622341    72     5    86     7    96    11    83     5    94     8


# QC for inputmethfile #
Methfile_QC <- function(inputmethfile, low_coveragenum = 10, high_coveragenum = NULL, low_quantile = NULL, 
                        high_quantile = 99.9, samplenum_QC = "all"){
  
  # sample number #
  samplenum <- (ncol(inputmethfile) - 2) / 2
  
  # Total coverage = methylated + unmethylated #
  inputmethfile <-  data.frame(inputmethfile, array(0, c(nrow(inputmethfile), (1 + samplenum))))
  
  j = 3
  for(i in 1:samplenum){
    inputmethfile[,(2+samplenum*2+i)] <- inputmethfile[,j] + inputmethfile[,(j+1)]
    j = j + 2
  }
  
  # coverage columns #
  coveragename_col <- (3+samplenum*2):ncol(inputmethfile)
  names(inputmethfile)[coveragename_col] <- c(paste("coverage", 1:samplenum, sep = ""), "unQCnum")
  
  # set coveragenum by low_quantile or high_quantile for each sample as quantile before coveragenum #
  if(is.null(low_quantile) == FALSE){
  
    # set tmpnum #
    tmpnum1 = min(inputmethfile[,coveragename_col])
    for(p in 1:samplenum){
      quantilenum <- quantile(inputmethfile[,(2+samplenum*2+p)], (low_quantile / 100))
      if(quantilenum >= tmpnum1){
        tmpnum1 <- quantilenum
      }
    }
    
    # choose the larger value between low_quantile and low_coveragenum #
    if(is.null(low_coveragenum) == TRUE){
      low_coveragenum <- tmpnum1
      
    }else if(low_coveragenum < tmpnum1){
      low_coveragenum <- tmpnum1
    }
  }  
  
  if(is.null(high_quantile) == FALSE){
  
    # set tmpnum #
    tmpnum2 = max(inputmethfile[,coveragename_col])
    for(q in 1:samplenum){
      quantilenum <- quantile(inputmethfile[,(2+samplenum*2+q)], (high_quantile / 100))
      if(quantilenum <= tmpnum2){
        tmpnum2 <- quantilenum
      }
    }
    # choose the samller value between high_quantile and high_coveragenum #
    if(is.null(high_coveragenum) == TRUE){
      high_coveragenum <- tmpnum2
      
    }else if(high_coveragenum > tmpnum2){
      high_coveragenum <- tmpnum2
    }
  }  

  ## filter ##
  # coverage number less than coveragenum, default 10 and then filtered #
  if(samplenum_QC == "all"){
    if(is.null(high_coveragenum) == TRUE){
      for(m in 1:samplenum){
        inputmethfile <- inputmethfile[inputmethfile[,2+samplenum*2+m] >= low_coveragenum, ]
      }
      
    }else{
      for(m in 1:samplenum){
        inputmethfile <- inputmethfile[inputmethfile[,2+samplenum*2+m] >= low_coveragenum &
                                         inputmethfile[,2+samplenum*2+m] <= high_coveragenum, ]
      }
    }
    
  }else{
    
    # for samplenum_QC < samplenum #
    if(samplenum_QC > samplenum){
      message("Wrong samplenum_QC value that is larger than total number")
      
    }else{
      
      # coverage number less than coveragenum, default 10 and then filtered #
      if(is.null(high_coveragenum) == TRUE){
        inputmethfile[,coveragename_col] <- inputmethfile[,coveragename_col] < low_coveragenum
      }else{
        inputmethfile[,coveragename_col] <- inputmethfile[,coveragename_col] < low_coveragenum | 
          inputmethfile[,coveragename_col] > high_coveragenum
      }
      
      # transform TRUE and FALSE to number 1 and 0 #
      inputmethfile[,coveragename_col][inputmethfile[,coveragename_col]=="TRUE"] <- 0
      inputmethfile[,coveragename_col][inputmethfile[,coveragename_col]=="FALSE"] <- 1
      inputmethfile$unQCnum <- apply(inputmethfile[,coveragename_col], 1, sum)
      
      # filter the CpG sites with coverage number and samplenum_QC #
      inputmethfile <- filter(inputmethfile, unQCnum < samplenum_QC)
    }
  }
    
    return(inputmethfile[,1:(2+samplenum*2)])
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




