#' Plot the methylation correlation.
#' 
#' @import corrplot
#' 
#' @description This function outputs the correlation plot for the methylation level of different samples or groups based on R package corrplot.
#' 
#' @param inputmethfile_QC refers to the input file with methylation levels, with default inputmethfile after quality control.
#' @param unmeth_exclude refers to whether to exclude the unmethylated sites or regions, with default TRUE.
#' 
#' @return Outputs a correlation figure.
#'
#' @examples
#' Correlation_plot(inputmethfile_QC)
#' Correlation_plot(siteall)
#' Correlation_plot(regiongenealls)
#' Correlation_plot(genefeatureall_cpgfeature)
#' Correlation_plot(genefeatureall_cpgfeature, unmeth_exclude = FALSE)
#' 
#' @export


Correlation_plot <- function(inputmethfile_QC, unmeth_exclude = TRUE){
  
  # calculate the total feature and group number from input inputmethfile_QC #
  featurenum <- length(grep("Methgroup1", colnames(inputmethfile_QC)))
  
  # judge group or sample #
  if(featurenum > 0){
    methpos <- grep("Methgroup", colnames(inputmethfile_QC))
    samplematrix <- inputmethfile_QC[, methpos]
    
    # NaN is filled by 0 for correlation #
    samplematrix[samplematrix == "NaN"] <- 0
    
  # for sample #  
  }else{
    
    # new matrix for all samples #
    samplepos <- grep("_", colnames(inputmethfile_QC))
    sampledata <- inputmethfile_QC[, samplepos]
    
    # methylation level #
    Cspos <- grep("Cs", colnames(sampledata))
    Tspos <- grep("Ts", colnames(sampledata))
    samplematrix <- sampledata[, Cspos] / (sampledata[, Cspos] + sampledata[, Tspos])
    
    # opitional #
    #samplematrix <- array(0, c(nrow(sampledata), (ncol(sampledata) / 2)))
    #for(i in 1:(ncol(sampledata) / 2)){
    
    # methylation level #
    #samplematrix[, i] <- sampledata[, (2*i - 1)] / (sampledata[, (2*i - 1)] + sampledata[, (2*i)])
    #}
    
    # rename by lapply function for vector strsplit#
    samplename <- colnames(sampledata)[grep("Cs", colnames(sampledata))]
    colnames(samplematrix) <- paste("Methsample", unlist(lapply(samplename, FUN = function(x) {return(strsplit(x, split = "Cs")[[1]][2])})))
  }
    
  # exclude the unmethylated rows #
  if(unmeth_exclude == TRUE){
    samplematrix <- data.frame(samplematrix, unmeth = apply(samplematrix, 1, sum))
    unmeth_samplematrix <- samplematrix[samplematrix$unmeth != 0, ]
    samplematrix <- unmeth_samplematrix[, -grep("unmeth", colnames(unmeth_samplematrix))]
  }
  
  # correlation #
  Maxcor <- cor(as.matrix(samplematrix))
  res <- cor.mtest(as.matrix(samplematrix), conf.level = .95)
  
  corrplot(Maxcor, type="upper", order="hclust", p.mat = res$p, 
           insig = "label_sig",sig.level = c(.001, .01, .05), pch.cex = .9, pch.col = "black")
}




 