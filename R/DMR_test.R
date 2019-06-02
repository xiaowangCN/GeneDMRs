#' Logistical regression analysis for each region or each cytosine site.
#' 
#' @description This function tests each region or each cytosine site by logistical regression model to
#' achieve the P values and then be adjusted to Q values to account for multiple hypothesis testing.
#'
#' @param genefeatureall_cpgfeature refers to the input file with methylation levels to be tested.
#' @param covariates refers to the extra covariates used in the model, with the default NULL.
#' @param adjustedmethod refers to the methods to adjust P values to Q values, with the default "fdr" method.
#' The adjustedmethod could be "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none" methods as well.
#'
#' @param diffgroup refers to methylation difference between two groups, with the default NULL, that is the max group - min group.
#' The two groups can be manually selected e.g. diffgroup = c("group1", "group2").
#'
#' @return Outputs a data frame of region gene or region cpgi or those regions with different features or cytosine sites, 
#' by accompanying with P values, Q values and methylation differences.
#' 
#' @examples 
#' regiongeneall_Qvalue <- Logic_regression(regiongeneall)
#' regiongenealls_Qvalue <- Logic_regression(regiongenealls)
#' regioncpgall_Qvalue <- Logic_regression(regioncpgall , adjustedmethod = "fdr")
#' regiongenebodyall_Qvalue <- Logic_regression(regiongenebodyall, diffgroup = c("group1", "group2"))
#' regiongeneall_cpgfeature_Qvalue <- Logic_regression(regiongeneall_cpgfeature)
#' genefeatureall_cpgfeature_Qvalue <- Logic_regression(genefeatureall_cpgfeature)
#' genefeatureall_Qvalue <- Logic_regression(genefeatureall, adjustedmethod = "bonferroni")
#' 
#' siteall_Qvalue <- Logic_regression(siteall, adjustedmethod = "fdr") # for each cytosine site #
#' siteall_Qvalue <- Logic_regression(siteall, adjustedmethod = "fdr", diffgroup = c("group1", "group2"))
#' 
#' @export


Logic_regression <- function(genefeatureall_cpgfeature, covariates = NULL, adjustedmethod = "fdr", diffgroup = NULL){
  
  # calculate the total feature and group number from input genefeatureall_cpgfeature #
  featurenum <- length(grep("Methgroup1", colnames(genefeatureall_cpgfeature)))
  groupnum <- length(grep("group", colnames(genefeatureall_cpgfeature))) / (2*featurenum)
  
  output <- data.frame(genefeatureall_cpgfeature, array(0,c(nrow(genefeatureall_cpgfeature), featurenum*3)))
  colnames(output) <- c(colnames(genefeatureall_cpgfeature), paste("Pvalue",1:featurenum,sep = ""),
                        paste("Qvalue",1:featurenum,sep = ""), paste("Methdiff",1:featurenum,sep = ""))
  
  # get the model matrix from treatment and (optional) covariates
  grouplevel = 0:(groupnum - 1)
  vars <- as.data.frame(cbind(grouplevel, covariates))
  
  if(ncol(vars) == 1){
    # get formula from model matrix
    formula <- as.formula(paste("~ ", paste(colnames(vars), collapse= "+")))
    
    # full model with all variables
    modelMat <- model.matrix(formula, as.data.frame(vars))
    
  }else{
    
    # if there are more covariates #
    fmlaCov <- as.formula(paste("~ ", paste(colnames(vars)[-1], collapse= "+")))
    modelCov <- model.matrix(as.formula(fmlaCov), as.data.frame(vars[, -1, drop=FALSE]))
  }
  
  # test each region or site #
  for(i in 1:nrow(output)){
    
    # count the reading line #
    if(i==100000){
      message("The calculating line is 100,000 now")
    }
    if(i==500000){
      message("The calculating line is 500,000 now")
    }
    if(i==1000000){
      message("The calculating line is 1,000,000 now")
    }
    
    # with gene or other regions without CpGisland and Shore #
    if(length(grep("CpGisland", colnames(genefeatureall_cpgfeature))) == 0){
      
      methtmp <- t(genefeatureall_cpgfeature[i, grep("Methgroup", colnames(genefeatureall_cpgfeature))])
      readtmp <- t(genefeatureall_cpgfeature[i, grep("Read", colnames(genefeatureall_cpgfeature))])
      
      # with CpGisland and Shore #  
    }else{
      methtmp <- t(genefeatureall_cpgfeature[i, grep("Methgroup", colnames(genefeatureall_cpgfeature))[1:groupnum]])
      readtmp <- t(genefeatureall_cpgfeature[i, grep("Read", colnames(genefeatureall_cpgfeature))[1:groupnum]])
    }
    
    tmpfile <- data.frame(group = row.names(methtmp), methylevel = methtmp[,1], readlevel = readtmp[,1])
    
    if(tmpfile$methylevel[1] != "NaN"){
      
      if(ncol(vars) == 1){
        obj <- glm.fit(modelMat, tmpfile$methylevel, weights = tmpfile$readlevel, family = binomial(link = "logit"))
        
      }else{
        obj <- glm.fit(modelCov, tmpfile$methylevel, weights = tmpfile$readlevel, family = binomial(link = "logit"))
      }
      
      deviance <- obj$null.deviance - obj$deviance
      if(deviance == 0){
        warning("Null model does not deviate from Reduced model.")
      } 
      
      # difference in degrees of freedom #
      ddf=obj$df.null-obj$df.residual 
      
      # assume no dispersion #
      phi = 1
      output[i, (ncol(genefeatureall_cpgfeature) + 1)] <- pchisq(deviance/phi, ddf, lower.tail = FALSE)
      
    }else{
      output[i, (ncol(genefeatureall_cpgfeature) + 1)] <- 1
    }
    
    # more cpg feature #
    if(featurenum > 1){
      tmpshore <- t(genefeatureall_cpgfeature[i, grep("Meth", colnames(genefeatureall_cpgfeature))[(1 + groupnum):(2*groupnum)]])
      readtmpshore <- t(genefeatureall_cpgfeature[i, grep("Read", colnames(genefeatureall_cpgfeature))[(1 + groupnum):(2*groupnum)]])
      tmpshorefile <- data.frame(group = row.names(tmpshore), methylevel = tmpshore[,1], readlevelshore = readtmpshore[1])
      
      if(tmpshorefile$methylevel[1] != "NaN"){
        
        if(ncol(vars) == 1){
          objshore <- glm.fit(modelMat, tmpshorefile$methylevel, weights = tmpshorefile$readlevelshore, family = binomial(link = "logit"))
          
        }else{
          objshore <- glm.fit(modelCov, tmpshorefile$methylevel, weights = tmpshorefile$readlevelshore, family = binomial(link = "logit"))
        }
        
        devianceshore <- objshore$null.deviance - objshore$deviance
        
        if(devianceshore == 0){
          warning("Null model does not deviate from Reduced model.")
        } 
        
        # difference in degrees of freedom #
        ddfshore=objshore$df.null-objshore$df.residual
        
        # assume no dispersion #
        phi = 1
        output[i, (ncol(genefeatureall_cpgfeature) + 2)] <- pchisq(devianceshore/phi, ddfshore, lower.tail = FALSE)
        
      }else{
        output[i, (ncol(genefeatureall_cpgfeature) + 2)] <- 1
      }
    }
  }
  
  # round to 3 digitals for P values #
  # output[, c((ncol(genefeatureall_cpgfeature) + 1), (ncol(genefeatureall_cpgfeature) + 2))] <- 
  # round(output[, c((ncol(genefeatureall_cpgfeature) + 1), (ncol(genefeatureall_cpgfeature) + 2))], digits = 3)
  
  # adjusted pvalue #
  # adjustedmethod = c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none") #
  adjustedoutput <- Pvalue_adjusted(output, adjustedmethod)
  
  # add the methylation difference between two groups #
  adjustedoutput_methdiff <- Meth_difference(adjustedoutput)
  
  return(adjustedoutput_methdiff)
}





#' Internal Use Function
#' Adjust the P values of logistical regression model.
#' 
#' @description This function adjust the P values of each region or each cytosine by logistical regression model.
#'
#' @param output refers to the output of Methmean_match.
#' @param adjustedmethod refers to the P values adjusted methods, with the default fdr.
#' The adjustedmethod could be "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "none" as well.
#'
#' @return Outputs the regiongene or regioncpg or of gene/cpg with different features with Q values.
#' 
#' @examples 
#' output_Qvalue <- Pvalue_adjusted(output, adjustedmethod = "fdr")
#' output_Qvalue <- Pvalue_adjusted(output)
#' output_Qvalue <- Pvalue_adjusted(output, adjustedmethod = "bonferroni")
#' output_Qvalue <- Pvalue_adjusted(output, adjustedmethod = "BH")
#' 
#' @export


Pvalue_adjusted <- function(output, adjustedmethod = "fdr"){

  # calculate the total feature #
  pvaluenum <- length(grep("Pvalue", colnames(output)))
  
  for(i in 1:pvaluenum){
    
    # get the column number of Pvalue and Qvalue #
    colnump <- grep("Pvalue", colnames(output))[i]
    colnumq <- grep("Qvalue", colnames(output))[i]
    
    # adjust P value #
    tmpoutput <- output[order(output[, colnump]), ]
	
	# assume P values of NaN = 1 #
	validrow <- nrow(tmpoutput)
    tmpoutput[, colnumq] <- p.adjust(tmpoutput[, colnump], method = adjustedmethod, n = validrow)
    output <- tmpoutput
  }
  
  # arrange chromosome and postion (posi for cytosine or start for region) #
  if(length(grep("start", colnames(output))) > 0){
    tmpoutput <- arrange(tmpoutput, chr, start)
	
  }else{
    tmpoutput <- arrange(tmpoutput, chr, posi)
  }
  
  return(tmpoutput)
}





#' Internal Use Function
#' Calculate the methylation difference between groups.
#' 
#' @description This function calculate the methylation difference between groups of each region or each cytosine.
#'
#' @param output refers to the output of Methmean_match.
#' @param diffgroup refers to methylation difference between two groups, with the default NULL that is the max group - min group.
#' The two groups can be manually select like diffgroup = c("group1", "group2").
#'
#' @return Outputs methylation difference between groups.
#' 
#' @examples 
#' output_methdiff <- Meth_difference(output)
#' adjustedoutput_methdiff <- Meth_difference(adjustedoutput, diffgroup = c("group1", "group2"))
#' 
#' @export


Meth_difference <- function(output, diffgroup = NULL){
  
  # calculate the total feature and group number from input  #
  featurenum <- length(grep("Methgroup1", colnames(output)))
  groupnum <- length(grep("group", colnames(output))) / (2*featurenum)
  
  for(i in 1: featurenum){
    
    # methylation difference column #
    colnum_methydiff <- grep("Methdiff", colnames(output))[i]
    
    if(is.null(diffgroup)==TRUE){
      
      # get the column number of Pvalue and Qvalue #
      colnum_groupstart <- grep("Methgroup1", colnames(output))[i]
      colnum_groupend <- grep(paste("Methgroup", groupnum, sep = ""), colnames(output))[i]
      
    }else{
      colnum_groupstart <- grep(paste("Meth", diffgroup[1], sep = ""), colnames(output))[i]
      colnum_groupend <- grep(paste("Meth", diffgroup[2], sep = ""), colnames(output))[i]
    }
    
    # calculate the methylation difference #
    output[, colnum_methydiff] <- output[, colnum_groupend] - output[, colnum_groupstart]
  }
  
  return(output)
}





#' Filter the significant regions or cytosine sites.
#' 
#' @description This function filters significant regions or cytosine sites based on Q value and methylation difference.
#'
#' @param genefeatureall_cpgfeature_Qvalue refers to the input file with Q values and methylation differences need to be filtered.
#' @param qvalue refers to the threshold of Q values that Q values larger than this will be discarded, with default 0.01.
#' @param methdiff refers to the threshold of methylation differences that methylation differences less than this will be discarded, with the default 0.
#' @param featureout refers to which feature will be filtered, with default 1. When featureout = 2, it means that the second feature will be filtered and outputted.
#'
#' @return Outputs a data frame of the significant regions or cytosine sites.
#' 
#' @examples 
#' genefeatureall_cpgfeature_significantcpgisland <- Significant_filter(genefeatureall_cpgfeature_Qvalue)
#' genefeatureall_cpgfeature_significantshore <- Significant_filter(genefeatureall_cpgfeature_Qvalue, featureout = 2)
#' genefeatureall_cpgfeature_significantcpgisland <- Significant_filter(genefeatureall_cpgfeature_Qvalue, qvalue = 0.001, methdiff = 0.01, featureout = 1)
#' regiongeneall_cpgfeature_significantcpgisland <- Significant_filter(regiongeneall_cpgfeature_Qvalue, methdiff = 0.05, featureout = 1)
#' regiongeneall_significant <- Significant_filter(regiongeneall_Qvalue)
#' regiongenealls_significant <- Significant_filter(regiongenealls_Qvalue, methdiff = 0.05)
#' siteall_significant <- Significant_filter(siteall_Qvalue)
#' siteall_significant <- Significant_filter(siteall_Qvalue, qvalue = 0.001, methdiff = 0.1)
#' 
#' @export


Significant_filter <- function(genefeatureall_cpgfeature_Qvalue, qvalue = 0.01, methdiff = 0, featureout = 1){
  
  # get the column number of Qvalue and Methylation difference #
  colnumq <- grep("Qvalue", colnames(genefeatureall_cpgfeature_Qvalue))[featureout]
  colnummeth <- grep("Methdiff", colnames(genefeatureall_cpgfeature_Qvalue))[featureout]
  
  # if one feature then output one and if one more features then output one more #
  tmpsigni <- genefeatureall_cpgfeature_Qvalue[genefeatureall_cpgfeature_Qvalue[, colnumq] < qvalue & 
                                                 abs(genefeatureall_cpgfeature_Qvalue[, colnummeth]) > methdiff, ]
  
  # rename the row names #
  rownames(tmpsigni) <- 1:nrow(tmpsigni)
  
  return(tmpsigni)
}





#' Annotate the differentially methylated cytosine (DMC) to features.
#' 
#' @description This function annotates the differentially methylated cytosine (DMC) after statistical test Logic_regression() to gene body or CpG island features.
#'
#' @param siteall_significant refers to the input file with DMC sites.
#' @param featureid refers to whether to include the feature id or not, with the default TRUE. 
#' The feature id will output the id of first file of the featurefile list e.g. the id of inputgenebodyfile.
#' 
#' @param featurefile refers to the input feature files e.g. inputgenebodyfile and inputcpgifeaturefile, 
#' with default two files in a list as featurefile = list(inputgenebodyfile, inputcpgifeaturefile),
#' and it can also be one file without a list e.g. featurefile = inputgenebodyfile.
#'
#' @return Outputs a data frame contains DMC sites with features.
#' 
#' @examples 
#' siteall_significant_feature <- DMC_feature(siteall_significant, featurefile = list(inputgenebodyfile, inputcpgifeaturefile))
#' siteall_significant_feature <- DMC_feature(siteall_significant, featureid = FALSE, featurefile = list(inputgenebodyfile, inputcpgifeaturefile))
#' siteall_significant_feature <- DMC_feature(siteall_significant, featureid = TRUE, featurefile = inputgenebodyfile)
#' 
#' @export


DMC_feature <- function(siteall_significant, featureid = TRUE, featurefile = list(inputgenebodyfile, inputcpgifeaturefile)){
  
  # judge featurefile is one file or more file in list #
  if(is.null(nrow(featurefile)) == TRUE){
    filenum <- length(featurefile)
    
  }else{
    filenum <- 1
  }
  
  # include feature id or not #
  if(featureid == TRUE){
    
    # if with feature id then add one column #
    output <- data.frame(siteall_significant, array(0,c(nrow(siteall_significant), 1 + length(featurefile))))
    
    # rename #
    colnames(output) <- c(colnames(siteall_significant), "Feature1_id", paste("Feature", 1:length(featurefile), sep = ""))
    
    # read the feature file #
    for(i in 1:filenum){
      
      # read the featurefile #
      if(filenum == 1){
        tmpfile <- featurefile
        
      }else{
        tmpfile <- featurefile[[i]]
      }
      
      # count the reading feature file #
      if(i==2){
        message("The second feature file is reading now")
      }
      
      # locate the DMC to the region possibly with features #
      for(j in 1:nrow(output)){
        
        # filter the features #
        tmpfeature <- filter(tmpfile, chr %in% as.vector(unlist(output[1]))[j], start <= output[j,2], end >= output[j,2])
        
        if(nrow(tmpfeature) == 0){
          
          # name the unannotated gene with "" only by reading first feature file #
          if(i == 1){
            output[j, (ncol(siteall_significant) + 1)] <- ""
          }
          
          # name the unannotated region of genebody as intergenic #
          if(names(tmpfile)[5] == "genebody"){
            output[j, (ncol(siteall_significant) + 1 + i)] <- "intergenic"
            
          }else{
            output[j, (ncol(siteall_significant) + 1 + i)] <- "Others"
          }
        }
        
        if(nrow(tmpfeature) == 1){
          
          # name the annotated gene with NCBI id only by reading first feature file #
          if(i == 1){
            output[j, (ncol(siteall_significant) + 1)] <- as.vector(unlist(tmpfeature[1]))
          }
          
          output[j, (ncol(siteall_significant) + 1 + i)] <- as.vector(unlist(tmpfeature[5]))
        }
        
        # sometimes one cytosine located in two features #
        if(nrow(tmpfeature) > 1){
          
		  # how many genes included #
		  lentmp <- length(as.vector(unlist(tmpfeature[1])))
		  
          # name the annotated gene with NCBI id only by reading first feature file #
          if(i == 1){
		  
		    # paste all the included genes #
		    tmpgenepaste <- as.vector(unlist(tmpfeature[1]))[1]
            for(k in 2:lentmp){
              tmpgenepaste <- paste(tmpgenepaste, as.vector(unlist(tmpfeature[1]))[k], sep = " & ")
            }
			
            output[j, (ncol(siteall_significant) + 1)] <- tmpgenepaste
          }
		  
		  # if one cytosine is located in more than four unique features, it is duplicated #
          if(length(unique(as.vector(unlist(tmpfeature[5])))) > 4){
            message(paste("The line", j, "of inputfile in featurefile", i, "is duplicated", sep = " "))
          }
		  
		  # paste all the included gene bodies #
          tmpgenebodypaste <- as.vector(unlist(tmpfeature[5]))[1]
          for(k in 2:lentmp){
            tmpgenebodypaste <- paste(tmpgenebodypaste, as.vector(unlist(tmpfeature[5]))[k], sep = " & ")
          }
		 	  
          output[j, (ncol(siteall_significant) + 1 + i)] <- tmpgenebodypaste
        }
      }
    }
    
  }else{
    
    # if without feature id #
    output <- data.frame(siteall_significant, array(0,c(nrow(siteall_significant), length(featurefile))))
    
    # rename #
    colnames(output) <- c(colnames(siteall_significant), paste("Feature", 1:length(featurefile), sep = ""))
    
    # read the feature file #
    for(i in 1:filenum){
      
      # read the featurefile #
      if(filenum == 1){
        tmpfile <- featurefile
        
      }else{
        tmpfile <- featurefile[[i]]
      }
      
      # count the reading feature file #
      if(i==2){
        message("The second feature file is reading now")
      }
      
      # locate the DMC to the region possibly with features #
      for(j in 1:nrow(output)){
        
        # filter the features #
        tmpfeature <- filter(tmpfile, chr %in% as.vector(unlist(output[1]))[j], start <= output[j,2], end >= output[j,2])
        
        if(nrow(tmpfeature) == 0){
          
          # name the unannotated region of genebody as intergenic #
          if(names(tmpfile)[5] == "genebody"){
            output[j, (ncol(siteall_significant) + i)] <- "intergenic"
            
          }else{
            output[j, (ncol(siteall_significant) + i)] <- "Others"
          }
        }
        
        if(nrow(tmpfeature) == 1){
          output[j, (ncol(siteall_significant) + i)] <- as.vector(unlist(tmpfeature[5]))
        }
        
        # sometimes one cytosine located in two features #
        if(nrow(tmpfeature) > 1){
		
		  # how many genes included #
		  lentmp <- length(as.vector(unlist(tmpfeature[1])))
		  
		  # if one cytosine is located in more than four unique features, it is duplicated #
          if(length(unique(as.vector(unlist(tmpfeature[5])))) > 4){
            message(paste("The line", j, "of inputfile in featurefile", i, "is duplicated", sep = " "))
          }
		  
		  # paste all the included gene bodies #
          tmpgenebodypaste <- as.vector(unlist(tmpfeature[5]))[1]
          for(k in 2:lentmp){
            tmpgenebodypaste <- paste(tmpgenebodypaste, as.vector(unlist(tmpfeature[5]))[k], sep = " & ")
          }
		  
          output[j, (ncol(siteall_significant) + i)] <- tmpgenebodypaste
        }
      }
    }
  }
  
  return(output)
}




