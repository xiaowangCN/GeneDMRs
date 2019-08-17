#' Internal Use Function
#' That calculates the mean methylation difference.
#'
#' @import dplyr
#' 
#' @description This function calculates the mean methylation difference along each read gene.
#' Formula:The mean methylation of a gene of one treatment group was calculated by: 
#' E((E[MR]_ij )/(E[TR]_ij ))*W_ij and W_ij =  (E[TR]_ij )/(EE[TR]_ij ), 
#' where [MR]_ij and [TR]_ij are methylated and total reads number of the involved CpG/DMC j at a given gene of individual i, 
#' n is the total individual number of one treatment group, m is total number of CpG/DMC involved in this gene and 
#' W_ij is the weight of reads of the involved CpG/DMC j of individual i. 
#' Firstly calculate methylation difference through reads of all cytosines along gene 
#' Secondly calculate methylation difference along each group
#' When provide inputcpgifeaturefile to calculate the gene mean methylation difference ,
#' two features will be combined that CpG island and CpG island shore.
#' 
#' @param inputmethfile_QC refers to the inputmethfile after quality control.
#' @param regionchr refers to the filtered subset of inputrefseqfile or inputcpgifile.
#' @param cpgifeaturefile refers to the whether have inputcpgifeaturefile or not.
#'
#' @return Outputs the mean methylation.
#'  
#' @examples
#' Meth_mean(inputmethfile_QC, regionchr)
#' Meth_mean(inputmethfile_QC, regionchr, cpgifeaturefile = inputcpgifeaturefile)
#' 
#' @export


Meth_mean <- function(inputmethfile_QC, regionchr, cpgifeaturefile = NULL){
  
  #calculate the total group number from input inputmethfile_QC
  groupnum <- length(grep("_1",colnames(inputmethfile_QC))) / 2
  
  # if no cpgifeaturefile #
  if(is.null(cpgifeaturefile)==TRUE){
    output <- data.frame(regionchr,array(0,c(nrow(regionchr),2*groupnum)))
    
    for(i in 1:nrow(regionchr)){
      
      # count the reading line #
      if(i==10000){
        message("The calculating line is 10,000 now")
      }
      if(i==30000){
        message("The calculating line is 30,000 now")
      }
      if(i==50000){
        message("The calculating line is 50,000 now")
      }
	  if(i==100000){
        message("The calculating line is 100,000 now")
      }
	  if(i==500000){
        message("The calculating line is 500,000 now")
      }
      
      #read each region from different chromosome and match all cytosines
      filterdat <- filter(inputmethfile_QC, chr %in% regionchr[i,2], posi >= regionchr[i,3], posi <= regionchr[i,4])
      
      for(j in 1:groupnum){
        methread <- summarise_at(filterdat,vars(contains(paste("Cs",j,"_",sep=""))),list(sum))
        unmethread <- summarise_at(filterdat,vars(contains(paste("Ts",j,"_",sep=""))),list(sum))
        output[i,(ncol(regionchr) + j)] <- sum(methread) / (sum(methread) + sum(unmethread))
		output[i,(ncol(regionchr) + groupnum + j)] <- sum(methread) + sum(unmethread)
      }
    }
    colnames(output) <- c(colnames(regionchr), paste("Methgroup",1:groupnum,sep = ""), paste("Readgroup",1:groupnum,sep = ""))
    
    #when the row of regionchr larther than 1, then delete the unmatched NA rows from group 1
    output <- filter(output,!is.na(output$Methgroup1))
    
    #Output the number of deleted rows
    print(paste("The number of unmatched", colnames(regionchr)[1],"is", (nrow(regionchr) - nrow(output)),sep = " "))
    
  }else{
    
    # if with cpgifeaturefile, every cpg site in CpG island will be counted #
    output <- data.frame(regionchr,array(0,c(nrow(regionchr), 2*2*groupnum)))
    
    # filter subset of cpg and shore from cpgifeaturefile #
    subcpg <- filter(cpgifeaturefile, cpgfeature %in% "CpGisland")
    subshore <- filter(cpgifeaturefile, cpgfeature %in% "Shores")
    
    for(i in 1:nrow(regionchr)){
      
      # count the reading line #
      if(i==10000){
        message("The calculating line is 10,000 now")
      }
      if(i==30000){
        message("The calculating line is 30,000 now")
      }
      if(i==50000){
        message("The calculating line is 50,000 now")
      }
	  if(i==100000){
        message("The calculating line is 100,000 now")
      }
	  if(i==500000){
        message("The calculating line is 500,000 now")
      }
      
      #read each region from different chromosome and match all cytosines
      filterdat <- filter(inputmethfile_QC, chr %in% regionchr[i,2], posi >= regionchr[i,3], posi <= regionchr[i,4])
      
      filtercpgdat = filtershoredat <- filterdat
      
      # when filterdat is not empty #
      if(nrow(filterdat) > 0){
        
        # further filter for CpG island or shore feature #
        for(k in 1:nrow(filterdat)){
          tmpcpg <- filter(subcpg, chr %in% filterdat[k,1], start <= filterdat[k,2], end >= filterdat[k,2])
          if(nrow(tmpcpg) == 0){
            filtercpgdat[k,2] <- 0
          }
          tmpshore <- filter(subshore, chr %in% filterdat[k,1], start <= filterdat[k,2], end >= filterdat[k,2])
          if(nrow(tmpshore) == 0){
            filtershoredat[k,2] <- 0
          }
        }
        
        # delete the zero rows of filtercpgdat and filtershoredat #
        filtercpgdat <- filter(filtercpgdat, posi > 0)
        filtershoredat <- filter(filtershoredat, posi > 0)
      }
      
      for(j in 1:groupnum){
        
        # first calculate the subset mean methylation with CpG island feature #
        methreadcpg <- summarise_at(filtercpgdat,vars(contains(paste("Cs",j,"_",sep=""))),list(sum))
        unmethreadcpg <- summarise_at(filtercpgdat,vars(contains(paste("Ts",j,"_",sep=""))),list(sum))
        output[i,(ncol(regionchr) + j)] <- sum(methreadcpg) / (sum(methreadcpg) + sum(unmethreadcpg))
		output[i,(ncol(regionchr) + groupnum + j)] <- sum(methreadcpg) + sum(unmethreadcpg)
        
        # second calculate the subset mean methylation with CpG island shore feature #
        methreadshore <- summarise_at(filtershoredat,vars(contains(paste("Cs",j,"_",sep=""))),list(sum))
        unmethreadshore <- summarise_at(filtershoredat,vars(contains(paste("Ts",j,"_",sep=""))),list(sum))
        output[i,(ncol(regionchr) + 2*groupnum + j)] <- sum(methreadshore) / (sum(methreadshore) + sum(unmethreadshore))
		output[i,(ncol(regionchr) + 3*groupnum + j)] <- sum(methreadshore) + sum(unmethreadshore)
      }
    }
    colnames(output) <- c(colnames(regionchr), paste("Methgroup", 1:groupnum, "_CpGisland", sep = ""), paste("Readgroup", 1:groupnum, "_CpGisland", sep = ""),
                          paste("Methgroup", 1:groupnum, "_Shore", sep = ""), paste("Readgroup", 1:groupnum, "_Shore", sep = ""))
    
    # when the row of regionchr larther than 1, then delete the unmatched NA rows from group 1 #
    output <- filter(output, !is.na(output$Methgroup1_CpGisland) | !is.na(output$Methgroup1_Shore))
    
    # Output the number of deleted rows #
    print(paste("The number of both unmatched", colnames(regionchr)[1],"is", (nrow(regionchr) - nrow(output)),sep = " "))
  }
  
  return(output)
}





#' Calculate the methylation mean for regions.
#'
#' @import dplyr
#' 
#' @description This function outputs the methylation mean for different groups based on gene and CpG island regions by matching with cytosine. 
#' It is also for gene body of promoter, exon, intron and TSSes regions, cgpi feature of CpG island and CpG island shores and their interactive regions e.g., promoter CpG island.
#' 
#' @param inputmethfile_QC refers to the input of methylation file after quality control.
#' @param inputrefseqfile refers to the input file with regions e.g., inputrefseqfile/inputcpgifile with 4 columns or inputgenebodyfile/inputcpgifeaturefile with 5 columns.
#' @param cpgifeaturefile refers to the input of CpG island feature file e.g., inputcpgifeaturefile, with default NULL.
#' If provided, the output file is methylation mean of inputrefseqfile or inputgenebodyfile with CpG island and CpG island shore features.
#' 
#' @param chrnum refers to the chromosome number or all chromosomes (all) or all chromosomes with unannotated sites (alls), with default "all".
#' @param posistart refers to start position if requested, with default NULL.
#' @param posiend refers to end position if requested, with default NULL.
#'
#' @param featureid refers to NCBI ID of specific gene or all the genes, with default NULL.
#' The CpG id can also be used like "cpgi1" or "shore2".
#' @param featurename refers to different gene body features of promoter, exon, intron and TSSes.
#' The CpG island  features can also be used that are "CpGisland" and "Shores".
#' 
#' @return Outputs a data frame of the methylation mean of provided regions with/without different features.
#'
#' @examples
#' Methmean_region(inputmethfile_QC, inputrefseqfile, chrnum = "alls", featureid = c("NM_001244353", "NM_001244864")) # find sepecific gene by NCBI ID #
#' 
#' Methmean_region(inputmethfile_QC, inputrefseqfile, chrnum = "chr1", posistart = 21800, posiend = 21900)
#' regiongenechr <- Methmean_region(inputmethfile_QC, inputrefseqfile, chrnum = c("chr1","chr2"))
#' regiongeneall <- Methmean_region(inputmethfile_QC, inputrefseqfile, chrnum = "all")
#' DMC_regiongeneall <- Methmean_region(DMC_inputmethfile_QC, inputrefseqfile, chrnum = "all") # Calculate DMC first and then recalculate the methylation mean by replacing the RRBS cytosine sites #
#' regiongenealls <- Methmean_region(inputmethfile_QC, inputrefseqfile, chrnum = "alls") # alls include unannotated CpG site like chrUn_NW_018084826v1 #
#' Methmean_region(inputmethfile_QC,inputcpgifile,"chr1", 21800, 21900) # acturally regiongenepart = regioncpgpart #
#' regioncpgchr <- Methmean_region(inputmethfile_QC, inputcpgifile, chrnum = c("chr1","chr2"))
#' regioncpgall <- Methmean_region(inputmethfile_QC, inputcpgifile, chrnum = "all")
#' regioncpgalls <- Methmean_region(inputmethfile_QC, inputcpgifile)
#' 
#' regiongenebodychr <- Methmean_region(inputmethfile_QC, inputgenebodyfile, chrnum = c("chr1","chr2"))
#' regiongenebodyall <- Methmean_region(inputmethfile_QC, inputgenebodyfile, chrnum = "all")
#' regiongenebodyalls <- Methmean_region(inputmethfile_QC, inputgenebodyfile)
#' regioncpgifeaturechr <- Methmean_region(inputmethfile_QC, inputcpgifeaturefile, chrnum = c("chr1","chr2"))
#' regioncpgifeatureall <- Methmean_region(inputmethfile_QC, inputcpgifeaturefile, chrnum = "all")
#' regioncpgifeaturealls <- Methmean_region(inputmethfile_QC, inputcpgifeaturefile)
#' 
#' partgenebody <- Methmean_region(inputmethfile_QC, inputgenebodyfile, featureid = "NM_001244353")
#' partgenebodyexon <- Methmean_region(inputmethfile_QC, inputgenebodyfile, featureid = "NM_001244353", featurename = "exons")
#' partgenebodyall <- Methmean_region(inputmethfile_QC, inputgenebodyfile, featureid = "NM_001244353", featurename = c("promoters","exons","introns","TSSes"))
#' genebodypromoterall <- Methmean_region(inputmethfile_QC, inputgenebodyfile, featureid = "all", featurename = "promoters")
#' genebodyexonall <- Methmean_region(inputmethfile_QC, inputgenebodyfile, featureid = "all", featurename = "exons") 
#' genefeatureall <- Methmean_region(inputmethfile_QC, inputgenebodyfile, featureid = "all", featurename = c("promoters","exons","introns","TSSes")) #long time #
#' partcpgi <- Methmean_region(inputmethfile_QC, inputcpgifeaturefile, featureid = "cpgi1")
#' partshore <- Methmean_region(inputmethfile_QC, inputcpgifeaturefile, featureid = "shore10")
#' cpgislandall <- Methmean_region(inputmethfile_QC, inputcpgifeaturefile, featureid = "all", featurename = "CpGisland")
#' cpgshoreall <- Methmean_region(inputmethfile_QC, inputcpgifeaturefile, featureid = "all", featurename = "Shores") #long time #
#' cpgfeatureall <- Methmean_region(inputmethfile_QC, inputcpgifeaturefile, featureid = "all", featurename = c("CpGisland", "Shores") #long time #
#' 
#' genebodychr_promoter <- Methmean_region(inputmethfile_QC, inputgenebodyfile, chrnum = "chr1", featureid = "all", featurename = "promoters")
#' cpgchr_island <- Methmean_region(inputmethfile_QC, inputcpgifeaturefile, chrnum = "chr1", featureid = "all", featurename = "CpGisland")
#'
#' # when the cpgifeaturefile = inputcpgifeaturefile is provided # 
#' regiongenechr_cpgfeature <- Methmean_region(inputmethfile_QC, inputrefseqfile, cpgifeaturefile = inputcpgifeaturefile, chrnum = c("chr1","chr2"))
#' regiongeneall_cpgfeature <- Methmean_region(inputmethfile_QC, inputrefseqfile, cpgifeaturefile = inputcpgifeaturefile, chrnum = "all")
#' regiongenealls_cpgfeature <- Methmean_region(inputmethfile_QC, inputrefseqfile, cpgifeaturefile = inputcpgifeaturefile, chrnum = "alls")
#' genebodypromoterall_cpgfeature <- Methmean_region(inputmethfile_QC, inputgenebodyfile, cpgifeaturefile = inputcpgifeaturefile, featureid = "all", featurename = "promoters") 
#' genebodyexonall_cpgfeature <- Methmean_region(inputmethfile_QC, inputgenebodyfile, cpgifeaturefile = inputcpgifeaturefile, featureid = "all", featurename = "exons") 
#' genefeatureall_cpgfeature <- Methmean_region(inputmethfile_QC, inputgenebodyfile, cpgifeaturefile = inputcpgifeaturefile, featureid = "all", 
#' featurename = c("promoters","exons","introns","TSSes")) #long time #
#' 
#' # windows # 
#' windowfileall <- Methmean_region(inputmethfile_QC, windowfile, chrnum = "all")
#' windowfilealls <- Methmean_region(inputmethfile_QC, windowfile, chrnum = "alls")
#' 
#' @export


Methmean_region <- function(inputmethfile_QC, inputrefseqfile, cpgifeaturefile = NULL, chrnum = "all", posistart = NULL, posiend = NULL, 
                            featureid = NULL, featurename = NULL){
  
  # rename inputrefseqfile #
  if(length(colnames(inputrefseqfile)) == 4){
    colnames(inputrefseqfile) <- c("id", "chr", "start", "end")
  }else{
    colnames(inputrefseqfile) <- c("id", "chr", "start", "end", "feature")
  }
  
  # if select one chromosome or more chromosomes with position #
  if(chrnum=="alls"){
    
    #when chrnum == "alls that include all chromosomes and undefined chromosomes"
    regionchr <- inputrefseqfile
    
  }else if(chrnum=="all"){   
  
    # find the unannotated chromosome rows and delete them #
    regionchr <- filter(inputrefseqfile, !grepl("_", chr))
    
  }else{
    
    # if given the start and end position #
    if(is.null(posistart)==FALSE & is.null(posiend)==FALSE){
      regionchr <- data.frame(region = paste(chrnum,":",posistart,"-",posiend,sep = ""), chr = chrnum, start = posistart, end = posiend)
    }else{
      
      # if chrnum!= "alls" & "all", then chrnum could be c("chr1","chr2") #
      regionchr <- filter(inputrefseqfile, chr %in% chrnum)
    }
  }
  
  # if chrnum is given like "chr1", "chr2", then the subset of regionchr will be accorded to #
  inputrefseqfile <- regionchr
  
  # if select id or id with features #
  if(is.null(featureid)==FALSE){
    
    # featureid mainly refers to the gene NCBI ID #
    if(featureid=="all"){
      
      # feature refers to promoter, exon, intron and TSSes regions or CpGisland and Shore regions #
      # featureid=="all" means all of genes or cpgis #
      if(is.null(featurename)==TRUE){
        regionchr <- inputrefseqfile
        
      }else{
        
        # when input file is inputgenebodyfile or inputcpgifeaturefile #
		tmpfeature <- unlist(lapply(X = inputrefseqfile$feature, FUN = function(x) {return(strsplit(x, split = "_")[[1]][1])}))
		
		# combine file with no "_" feature name #
		inputrefseqfile <- data.frame(inputrefseqfile, tmpfeature_select = tmpfeature)
        regionchr <- filter(inputrefseqfile, tmpfeature_select %in% featurename)
		regionchr <- regionchr[, -ncol(regionchr)]
      }
      
      # when featureid!="all", and featureid is the specific gene name or cpgi id #
    }else{
      
      if(is.null(featurename)==TRUE){
        regionchr <- filter(inputrefseqfile, id %in% featureid)
        
      }else{
	  
	    # combine file with no "_" feature name #
		tmpfeature <- unlist(lapply(X = inputrefseqfile$feature, FUN = function(x) {return(strsplit(x, split = "_")[[1]][1])}))
		inputrefseqfile <- data.frame(inputrefseqfile, tmpfeature_select = tmpfeature)
        regionchr <- filter(inputrefseqfile, id %in% featureid, tmpfeature_select %in% featurename)
		regionchr <- regionchr[, -ncol(regionchr)]
      }
    }
	
  }else{
    
    # if chrnum is given like "chr1", "chr2"#
    regionchr <- inputrefseqfile
  }
  
  # output file of the mean methylation with chromosome, start, end, group1, group2.... #
  if(is.null(cpgifeaturefile)==TRUE){
    
    output <- Meth_mean(inputmethfile_QC, regionchr)
  }else{
    
    # if have inputcpgifeaturefile #
    output <- Meth_mean(inputmethfile_QC, regionchr, cpgifeaturefile = inputcpgifeaturefile)
  }
  
  return(output)  
}





#' Calculate the methylation mean for cytosine sites.
#' 
#' @import dplyr
#' 
#' @description This function outputs the methylation mean for each cytosine site. It will calculate methylation difference along each group. 
#' 
#' Formula:The mean methylation of a cytosine of one treatment group was calculated by: 
#' E((E[MR]_i )/(E[TR]_i ))*W_i and W_i =  (E[TR]_i )/(EE[TR]_i ), 
#' where [MR]_i and [TR]_i are methylated and total reads number of the CpG, 
#' n is the total individual number of one treatment group, and 
#' W_i is the weight of reads of the CpG of individual i. 
#' 
#' @param inputmethfile_QC refers to the input of methylation file after quality control.
#'
#' @return Outputs a data frame of the methylation mean of provided cytosine sites.
#'  
#' @examples
#' siteall <- Methmean_site(inputmethfile_QC)
#' 
#' @export


Methmean_site <- function(inputmethfile_QC){
  
  # calculate the total group number from input inputmethfile_QC #
  groupnum <- length(grep("_1",colnames(inputmethfile_QC))) / 2
  output <- data.frame(inputmethfile_QC[, c(1,2)],array(0,c(nrow(inputmethfile_QC),2*groupnum)))
  
  for(i in 1:nrow(inputmethfile_QC)){
    
    # count the reading line #
    if(i==10000){
      message("The calculating line is 10,000 now")
    }
    if(i==30000){
      message("The calculating line is 30,000 now")
    }
    if(i==50000){
      message("The calculating line is 50,000 now")
    }
    if(i==100000){
      message("The calculating line is 100,000 now")
    }
	if(i==500000){
      message("The calculating line is 500,000 now")
    }
    
    for(j in 1:groupnum){
      methread <- summarise_at(inputmethfile_QC[i,],vars(contains(paste("Cs",j,"_",sep=""))),list(sum))
      unmethread <- summarise_at(inputmethfile_QC[i,],vars(contains(paste("Ts",j,"_",sep=""))),list(sum))
      output[i,(2 + j)] <- sum(methread) / (sum(methread) + sum(unmethread))
      output[i,(2 + groupnum + j)] <- sum(methread) + sum(unmethread)
    }
  }
  colnames(output) <- c(colnames(inputmethfile_QC)[c(1,2)], paste("Methgroup",1:groupnum,sep = ""), paste("Readgroup",1:groupnum,sep = ""))
  
  # delete the unmatched NA rows from group 1 #
  output <- filter(output, !is.na(output$Methgroup1))
  
  # Output the number of deleted rows #
  print(paste("The number of unmatched sites", "is", (nrow(inputmethfile_QC) - nrow(output)),sep = " "))

  return(output)
}




