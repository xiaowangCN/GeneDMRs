#' Divide the genome to windows.
#' 
#' @description This function outputs the window regions of the whole genome.
#' 
#' @param inputcytofile refers to the input cyto file with chromosome information.
#' @param windowbp refers to the window length in base pair (bp) to be divided, with default 1,000,000.
#' 
#' @return Outputs a data frame with window regions.
#'
#' @examples
#' windowfile <- Window_divide(inputcytofile)
#' windowfile <- Window_divide(inputcytofile, windowbp = 10000)
#' 
#' @export


Window_divide <- function(inputcytofile, windowbp = 1000000){
  
  # calculate the position for each e.g., 1 million bp window #
  inputcytofile$start <- ceiling(inputcytofile$end / windowbp) #or floor#
  
  output <- array(0, c(sum(inputcytofile$start), 3))
  colnames(output) <- c("chr","start","end")
  
  #start rearrange the position for output file#
  a = 1
  j = 1
  for(i in 1:nrow(inputcytofile)){
    b = 0
    while(j >= a & j < inputcytofile[i,2]){
      output[j,1] <- as.character(inputcytofile$chr)[i]
      output[j,2] <- 1 + b
      b <- b + windowbp
      output[j,3] <- b
      j = j + 1
    }
    #the last position for each chromosome#
    if(j == inputcytofile$start[i]){
      output[j,1] <- as.character(inputcytofile$chr)[i]
      output[j,2] <- 1 + b
      output[j,3] <- inputcytofile$end[i]
      j = j + 1
    }
    a <- inputcytofile$start[i]
    if(i < nrow(inputcytofile)){
      inputcytofile[(i+1),2] <- inputcytofile[(i+1),2] + inputcytofile[i,2]
    }
  }
  
  # make output as list combining chr character with position numbe #
  output <- data.frame(window = paste("Window", 1:nrow(output), sep = ""), chr = output[,1], 
                       start = as.numeric(output[,2]), end = as.numeric(output[,3]))
  
  return(output)
}




 