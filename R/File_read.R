#' Read the methylation file.
#' 
#' @import dplyr
#' 
#' @description This function reads all of the methylation files and generates one file with all samples including methylated read coverages (Cs) and unmethylated read coverages (Ts).
#' It can automatically test how many samples and how many replicates in each group and the distribute them from 1_1, 1_2 to the final file by headers.
#' The methylation files should be the standard coverage file (i.e., .bismark.cov) outputted from Bismark software.
#' The dataset of the example is the Reduced representation bisulfite sequencing (RRBS) data of DNA methylation for mouse myeloid progenitor tissue from GEO (Accession number: GSE62392)
#' (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE62392).
#' 
#' @param paths refers to the path of methylation file, with default the package path.
#' @param suffix refers to the suffix of methylation file, e.g., ".gz", ".zip" and so on (some files are in text .txt format, then ".txt" or ".txt.gz"), with default ".gz".
#' 
#' @return Outputs a data frame contain chromosome, position, and Cs & Ts for different replicates and groups.
#' 
#' @examples
#' inputmethfile <- Methfile_read() 
#' inputmethfile <- Methfile_read(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""), suffix = ".gz") # test the functions with default parameters #
#' head(inputmethfile)
#' 
#' @export


# So, what are the origanal coverage file and combined methall file looking like? #

 # Coverage file from Bismark without header, so don't add the header. If the cov file is with header, please use "chr", "posi", "Cs", "Ts" for four columns #
 # It contains six columns that are choromsome, start position, end position, methylated percentage, methylated read coverage and unmethylated read coverage #
 # Methylated percentage = Methylated read coverage / (Methylated read coverage and Unmethylated read coverage) #
 #1 chr1	3391	3391	100	1	0
 #2 chr1	3392	3392	0	0	2
 #3 chr1	3396	3396	100	1	0
 #4 chr1	3397	3397	100	2	0
 #5 chr1	3413	3413	100	1	0
 #6 chr1	3414	3414	100	2	0
 #7 chr1	3430	3430	100	1	0
 #8 chr1	3431	3431	100	2	0
 #9 chr1	3437	3437	100	1	0
 #10 chr1	3438	3438	100	2	0

 
 # The output file of this function is with header and combines all the files with the same chromosome and position #
 # It contains choromsome, position, methylated read coverages (Cs) and unmethylated read coverages (Ts) of each sample.
    #chr  posi Cs1_1 Ts1_1 Cs1_2 Ts1_2 Cs1_3 Ts1_3 Cs2_1 Ts2_1 Cs2_2 Ts2_2 Cs2_3 Ts2_3 Cs3_1 Ts3_1 Cs3_2 Cs3_3 Ts3_3
 #1  chr1 21896     0     4     1    49     0    25     4    24     5    28     0    33     1    43     0    0     6
 #2  chr1 21897     0    10     1   107     1    43     0    71     0    47     1    69     1    88     0    0    19
 #3  chr1 21904     0     4     0    50     0    25     0    28     0    33     1    33     0    44     0    0     6
 #4  chr1 21905     0    10     0   108     0    44     0    71     1    46     0    70     1    88     1    0    19
 #5  chr1 21912     0     4     0    50     0    25     0    28     0    33     0    34     0    44     0    0     6
 #6  chr1 21913     0    10     0   106     0    44     0    72     0    47     0    70     0    89     0    0    19
 #7  chr1 21916     0     4     0    50     1    24     0    28     0    33     0    34     0    44     4    0     6
 #8  chr1 21917     0    10     1   105     0    44     1    70     0    46     0    70     0    89     0    0    19
 #9  chr1 21925     0     4     0    50     0    25     0    28     0    33     0    34     0    44     0    0     6
 #10 chr1 21926     0    10     0   108     0    44     0    72     0    46     0    70     0    89     0    0    19


# The function for automatically reading different replicates of different groups #

Methfile_read <- function(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""), suffix = ".gz"){
  # Output the begin time
  print(paste("The start time is ", date(), sep = ""))

  # set the paths #
  setwd(paths)
  
  # count the number of file
  filenum <- 0
  
  # detect the total group number 
  groupnum <- length(grep("_1",dir(paths)))
  
  # from different group to different relicate files
  for(i in 1:groupnum){
    
    # dectect the real replicate number of this group
    # output the real replicate number
    realreplicatenum <- length(grep(paste(i,"_",sep = ""),dir(paths)))
    
	# read the file from each real replicate of each group and combine the same chromosome and position #
    for(j in 1:realreplicatenum){
	
      # read cov file (without header) from Bismark output and save the column of read number of Ts and Cs #
      tmpfile <- read.table(paste(i, "_", j, suffix, sep = ""), header = F)
	  
	  # If the cov file is with header, please use "chr", "posi", "Cs", "Ts" for four columns #
	  if(tmpfile[1,1] == "chr" | tmpfile[1,2] == "posi"){
	    tmpfile <- read.table(paste(i, "_", j, suffix, sep = ""), header = T)
	    tmpfile <- tmpfile[, c(tmpfile$chr, tmpfile$posi, tmpfile$Cs, tmpfile$Ts)]
		
	  }else{
	    tmpfile <- tmpfile[, -c(3,4)]
	  }
	  
	  # rename the header by replicate number of group number #
      colnames(tmpfile) <- c("chr","posi",paste("Cs",i,"_",j,sep = ""),paste("Ts",i,"_",j,sep = ""))
	  
      filenum <- filenum + 1
	  
	  # combine each replicate file to the total dataset #
      ifelse(filenum==1, methallfile <- tmpfile, methallfile <- inner_join(methallfile, tmpfile, by = c('chr'='chr','posi'='posi')))
    }
  }
  
  # Sort by chromosome and position #
  methallfile <- arrange(methallfile,chr,posi)
  
  # Output how many files are read and the running timen #
  print(paste("Total methylation file number is", filenum, "and the reading finish time is", date(), sep = " "))
  
  return(methallfile)
}

# Finally output the methallfile as described before #
# If there are some Warning messages likes this: Column `chr` joining factors with different levels, coercing to character vector. #
# It doesn't matter because the tool dplyr tested that one of the columns was a factor and that factor had different levels in the different datasets. #
# In order not to lose any information, the factors were converted to character values. #





#' Read the standard bedfile of refseq or cpgi downloaded from UCSC.
#' 
#' @import dplyr
#' @import genomation
#' 
#' @description This function reads the bed file of refseq or cpgi and sorts them by chromosome and position. 
#' The dataset of the example are the mouse reference genes and CpG island information that are downloaded from UCSC website 
#' (http://genome.ucsc.edu/cgi-bin/hgTables). 
#' The R package genomation used here can divide the refseq.bed file into several gene body features, e.g., promoter, exon, intron regions and 
#' the cpgi.bed file into CpG island features, e.g., CpG island and CpG island shore.
#' 
#' @param paths refers to the path of bed file, with default the package path.
#' @param bedfile refers to the file name of bed file like "refseq" or "cpgi". This file is downloaded from UCSC website, with default "refseq".
#' @param suffix refers to the suffix of bed file, e.g., ".gz", ".zip" and so on (some files are in text .txt format, then ".txt" or ".txt.gz"), with default ".txt".
#' @param feature refers to whether to read the bed with the features, with default FALSE. 
#' If feature = TRUE, the output of this function will contain the features e.g., promoter, exon, intron or CpG island, CpG island shore based on R package genomation.
#'              
#' @param featurewrite refers to whether to write out the feature file to the given path, with default FALSE.
#' 
#' @return Outputs a data frame contains four columns of chromosome, start position, end position. 
#' If feature = TRUE, the data frame is five columns with the added feature such as genebody or cpgfeature.
#' 
#' @references Akalin A, Franke V, Vlahovicek K, Mason C, Schubeler D (2014). "genomation: a toolkit to summarize, annotate and visualize genomic intervals." 
#' Bioinformatics. doi: 10.1093/bioinformatics/btu775, http://bioinformatics.oxfordjournals.org/content/early/2014/12/04/bioinformatics.btu775.long.
#' 
#' @examples
#' inputrefseqfile <- Bedfile_read() 
#' inputrefseqfile <- Bedfile_read(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""), bedfile = "refseq", suffix = ".txt", feature = FALSE)
#' inputcpgifile <- Bedfile_read(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""), bedfile = "cpgi", suffix = ".txt", feature = FALSE)
#' head(inputrefseqfile)
#' head(inputcpgifile)
#' 
#' inputgenebodyfile <- Bedfile_read(bedfile = "refseq", feature = TRUE, featurewrite = TRUE)
#' inputcpgifeaturefile <- Bedfile_read(bedfile = "cpgi", feature = TRUE, featurewrite = FALSE)
#' head(inputgenebodyfile)
#' head(inputcpgifeaturefile)
#' 
#' @export


# So, what are the origanally standard bedfile and outputed file looking like? #

 # Bedfile is without header, so don't add the header #
 
 # refseq.bed file #
 # It contains 12 columns that are choromsome, chromstart (thickStart) position, chromend (thickEnd) position, #
 # NCBI ID number for mRNA, score, strand, coding start position, coding end position, score, number of exon, length of exon, distance from TSS start position to exon #
 #1 chr1	75465296	75508244	NM_001243664	0	-	75465385	75508212	0	8	201,136,57,148,63,189,96,155,	0,11425,20986,26474,26728,35769,39257,42793,
 #2 chr1	5698507	6731132	NM_001044603	0	+	5698507	6729846	0	11	166,229,122,84,116,137,62,150,84,118,1399,	0,160362,224013,370124,439504,592645,812774,826595,994107,1025825,1031226,
 #3 chr1	73253304	73510244	NM_001167653	0	-	73253304	73510244	0	8	159,33,132,174,75,199,135,278,	0,27538,48205,49420,82314,106754,149871,256662,
 #4 chr1	77585913	77642342	NM_001080206	0	-	77585913	77642342	0	11	209,132,154,77,180,165,150,104,99,97,247,	0,10778,30879,31123,32747,36536,39881,41003,44299,52245,56182,
 #5 chr1	85758623	86017579	NM_001123219	0	-	85758623	86017579	0	10	141,87,111,98,109,156,134,69,78,208,	0,71926,167857,177254,203119,205053,205423,243559,246232,258748,
 #6 chr1	107997304	108131723	NM_001243563	0	+	108003699	108131137	0	12	99,98,222,139,130,49,27,153,46,90,84,667,	0,6374,73642,82787,101201,107787,108541,117080,118263,130344,130541,133752,
 #7 chr1	122633627	122688609	NM_001206341	0	+	122633907	122685089	0	18	304,179,325,136,230,25,116,135,81,168,142,170,133,163,221,155,220,4200,	0,8356,10502,13455,17533,18050,18724,23354,25146,25780,26091,29689,30624,32730,34495,36698,43102,50782,
 #8 chr1	128960122	129013424	NM_214171	0	-	128960728	129013201	0	24	633,59,117,79,69,65,58,78,114,18,37,209,12,170,161,78,86,84,147,169,134,119,70,532,	0,1029,1338,1658,1825,2366,2536,2948,4018,6625,8220,8935,9780,10414,11725,12641,14652,15615,17664,18622,19485,21066,23205,52770,
 #9 chr1	188704259	188925855	NM_001244640	0	-	188705091	188925855	0	9	875,59,47,70,139,208,732,774,241,	0,5651,5795,6353,7825,10764,108217,118752,221355,
 #10 chr1	189690481	189907808	NM_001166316	0	+	189690481	189907808	0	8	89,153,74,104,141,126,122,121,	0,55712,57518,67005,73033,78455,118475,217206,

  # cpgi.bed file #
 # It contains 4 columns that are choromsome, start position, end position, CpG ID #
 #1 chr1	21811	22330	CpG:_48
 #2 chr1	23707	24083	CpG:_41
 #3 chr1	66380	66649	CpG:_22
 #4 chr1	91255	91533	CpG:_19
 #5 chr1	122129	122347	CpG:_20
 #6 chr1	139652	140059	CpG:_32
 #7 chr1	160025	160481	CpG:_33
 #8 chr1	160781	160986	CpG:_20
 #9 chr1	174275	174652	CpG:_36
 #10 chr1	186372	187350	CpG:_119

 
 # The output file of this function is with header and combines with the other features #
 
 # The inputrefseqfile contains NCBI ID of gene, chromosome, start position, end position.
 #        refseq  chr   start     end
 # 1 NM_001244353 chr1   23826   40033
 # 2    NR_128500 chr1 1013706 1013802
 # 3    NR_128506 chr1 1340182 1340263
 # 4    NR_128509 chr1 1665567 1665643
 # 5 NM_001244864 chr1 2541536 2552860
 # 6 NM_001244534 chr1 2576814 2593516
 # 7    NR_128518 chr1 2893741 2893818
 # 8 NM_001007195 chr1 4738530 4883241
 # 9 NM_001044603 chr1 5698507 6731132
 # 10 NM_001143697 chr1 6807761 6883255
 
 # The inputcpgifile contains CpG ID, chromosome, start position, end position.
 #       cpgi  chr  start    end
 # 1   CpG:_48 chr1  21811  22330
 # 2   CpG:_41 chr1  23707  24083
 # 3   CpG:_22 chr1  66380  66649
 # 4   CpG:_19 chr1  91255  91533
 # 5   CpG:_20 chr1 122129 122347
 # 6   CpG:_32 chr1 139652 140059
 # 7   CpG:_33 chr1 160025 160481
 # 8   CpG:_20 chr1 160781 160986
 # 9   CpG:_36 chr1 174275 174652
 # 10 CpG:_119 chr1 186372 187350

 # The inputgenebodyfile contains NCBI ID of gene, chromosome, start position, end position and gene feature.
 #        refseq  chr start   end    genebody
 # 1 NM_001244353 chr1 22826 24826   promoters
 # 2 NM_001244353 chr1 23826 23826       TSSes
 # 3 NM_001244353 chr1 23827 23965       exons
 # 4 NM_001244353 chr1 23966 27359     introns
 # 5 NM_001244353 chr1 27360 27467       exons
 # 6 NM_001244353 chr1 27468 30767     introns
 # 7 NM_001244353 chr1 30768 30849     exons
 # 8 NM_001244353 chr1 30850 33299   introns
 # 9 NM_001244353 chr1 33300 33429     exons
 # 10 NM_001244353 chr1 33430 38650   introns
  
 # The inputcpgifeaturefile contains chromosome, start position, end position and CpGfeature.
 #     cpgi  chr start   end cpgfeature
 # 1  shore1 chr1 19811 21810     Shores
 # 2   cpgi1 chr1 21811 22330  CpGisland
 # 3  shore2 chr1 22331 23706     Shores
 # 4   cpgi2 chr1 23707 24083  CpGisland
 # 5  shore3 chr1 24084 26083     Shores
 # 6  shore4 chr1 64380 66379     Shores
 # 7   cpgi3 chr1 66380 66649  CpGisland
 # 8  shore5 chr1 66650 68649     Shores
 # 9  shore6 chr1 89255 91254     Shores
 # 10  cpgi4 chr1 91255 91533  CpGisland
  
  
Bedfile_read <- function(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""), bedfile = "refseq", suffix = ".txt", feature = FALSE, featurewrite = TRUE){
  
  # set the paths #
  setwd(paths)
  
  # read refseq or cpgi file #
  if(bedfile=="refseq"){
    if(feature == FALSE){
      regionfile <- Regionfile_read(bedfile, suffix)
      
      return(regionfile)
      
    }else{
      
      # region file is the reference gene or cpg bed file (without header) #
      regionfile <- Regionfile_read(bedfile, suffix)
      geneobj <- readTranscriptFeatures(paste(bedfile, ".bed", suffix, sep = ""))
      
      # Output the defined gene features based on refseq.bed file#
      write.table(geneobj, "Genebody", col.names = F, row.names = F, quote = F)
      
      geneobj <- read.table("Genebody")
      
      # check the geneobj file for promoter regions #
      geneobj <- refseqfile_check(geneobj, regionfile)
      
      # if featurewrite == FALSE then delete geneobj file #
      if(featurewrite == FALSE){
        file.remove("Genebody")
      }
      
      return(geneobj)
    }
    
    
  }else if(bedfile=="cpgi"){
    if(feature == FALSE){
      regionfile <- Regionfile_read(bedfile, suffix)
      
      return(regionfile)
      
    }else{
      
      # region file is the reference gene or cpg bed file (without header) #
      cpgiobj <- readFeatureFlank(paste(bedfile, ".bed", suffix, sep = ""), feature.flank.name = c("CpGisland", "Shores"))
      
      # Output the defined gene features based on refseq.bed file #
      write.table(cpgiobj, "Cpgifeature", col.names = F, row.names = F, quote = F)

      cpgiobj <- read.table("Cpgifeature")
      cpgiobj[,6] <- cpgiobj[,2]
      cpgiobj <- cpgiobj[,3:6]
      colnames(cpgiobj) <- c("chr","start","end","cpgfeature")
	  
	  # rename the "cpgi" #
	  cpginum <- sum(cpgiobj$cpgfeature == "CpGisland")
	  cpginame <- c(paste("cpgi", 1:cpginum, sep = ""), paste("shore", 1:(nrow(cpgiobj) - cpginum), sep = ""))
	  cpgiobj <- data.frame(cpgi = cpginame, cpgiobj)
	        
      # sort the file #
      cpgiobj <- arrange(cpgiobj,chr,start)
      
      # if featurewrite == FALSE then delete cpgiobj file #
      if(featurewrite == FALSE){
        file.remove("Cpgifeature")
      }
      
      return(cpgiobj)
    }
    
  }else{
    
    # inform only refseq or cpgi file can be read #
    stop("Wrong bed file names and please use refseq or cpgi name")
  }
}





#' Internal Use Function
#' That reads the bed file.
#'
#' @description This function reads the bed file and sort them.
#' 
#' @param bedfile refers to the bedfile name (e.g., refseq or cpgi) downloaded from UCSC, with the default refseq.
#' @param suffix refers to the compressed file suffix such as ".gz", ".zip" and so on, with the default suffix ".gz".
#'
#' @return Outputs bed file that can be used in the next step.
#'  
#' @examples
#' regionfile <- Regionfile_read(bedfile, suffix)
#' 
#' @export


Regionfile_read <- function(bedfile, suffix){

  # region file is the reference gene or cpg bed file (without header)
  regionfile <- read.table(paste(bedfile, ".bed", suffix, sep = ""), header=F)
  
  regionfile <- cbind(regionfile[, c(4,1,2,3)], 0)
  colnames(regionfile) <- c(bedfile, "chr", "start", "end", "uniid")
  
  # unique the repeated positions #
  regionfile$uniid <- paste(regionfile$chr, regionfile$start, regionfile$end, sep = "_")
  regionfile <- distinct(regionfile, uniid, .keep_all = TRUE)
  
  #sort by chromosome and position
  regionfile <- arrange(regionfile, chr, start)
  
  return(regionfile[,-5])
}





#' Internal Use Function
#' That checks and annotates the promoter to the specific gene. 
#'
#' @description This function sorts and makes sure geneobj file can be used in this package,
#' that mainly works on the promoters without the choromsome and position information.
#' In addition, it can add the genomic number of exons and introns, e.g., first exon or second intron.
#' 
#' @param geneobj refers to the genebody file from readTranscriptFeatures of package genomation.
#' @param regionfile refers to the reference gene or cpg bed file (without header)
#'
#' @return Outputs gene body file.
#'  
#' @examples
#' geneobj <- refseqfile_check(geneobj, regionfile)
#' 
#' @export


refseqfile_check <- function(geneobj, regionfile){
  geneobj <- geneobj[,c(9,3,4,5,2,7)]
  colnames(geneobj) <- c("refseq","chr","start","end","genebody","strand")
  
  ## promoters without the choromsome and position information ##
  tmp_promoter <- filter(geneobj, genebody == "promoters")
  tmp_nopromoter <- filter(geneobj, genebody != "promoters")
  print(paste("The total reading line for promoter is", nrow(tmp_promoter), sep = " "))
  
  for(i in 1:nrow(tmp_promoter)){
    
    # fill up the missing information of promoters #
    if(as.vector(unlist(tmp_promoter$strand))[i] == "+"){
      
      # when the strand is "+", then start position of gene is in the middle of promoter #
      tmp <- filter(regionfile, chr == as.vector(unlist(tmp_promoter$chr))[i], start == (tmp_promoter[i,3] + 1000))
    }else{
      
      # when the strand is "-", then end position of gene is in the middle of promoter #
      tmp <- filter(regionfile, chr == as.vector(unlist(tmp_promoter$chr))[i], end == (tmp_promoter[i,3] + 1000))
    }
    
    # replace '.' by NCBI ID #
    # if the chromosome and position are same but with two NCBI IDs, juse use the first one #
    tmp_promoter[i,1] <- as.vector(unlist(tmp$refseq))[1]
  }
  
  geneobj <- rbind(tmp_promoter, tmp_nopromoter)
  
  # unique the repeated positions #
  geneobj$strand <- paste(geneobj$chr, geneobj$start, geneobj$end, sep = "_")
  geneobj <- distinct(geneobj, strand, .keep_all = TRUE)
  
  
  ## add the number of exon and intron ##
  tmp_exonintron <- filter(geneobj, genebody == "exons" | genebody == "introns")
  tmp_noexonintron <- filter(geneobj, genebody != "exons" & genebody != "introns")
  
  # sort the file #
  tmp_exonintron <- arrange(tmp_exonintron, chr, start, refseq, genebody)
  
  # add one first row and one last row #
  tmp_exonintron <- rbind(tmp_noexonintron[1, ], tmp_exonintron, tmp_noexonintron[1, ])
  numexon <- 1
  numintron <- 1
  for(j in 2:(nrow(tmp_exonintron) - 1)){
    
    # for same gene #
    if(as.vector(unlist(tmp_exonintron$refseq))[j] == as.vector(unlist(tmp_exonintron$refseq))[j + 1]){
      if(as.vector(unlist(tmp_exonintron$genebody))[j] == "exons"){
        tmp_exonintron[j, ncol(tmp_exonintron)] <- paste("exons", numexon, sep = "_")
        numexon <- numexon + 1
        
      }else if(as.vector(unlist(tmp_exonintron$genebody))[j] == "introns"){
        tmp_exonintron[j, ncol(tmp_exonintron)] <- paste("introns", numintron, sep = "_")
        numintron <- numintron + 1
      }
      
    }else{
      
      # last gene body of the same gene #
      if(as.vector(unlist(tmp_exonintron$refseq))[j] == as.vector(unlist(tmp_exonintron$refseq))[j - 1]){
        if(as.vector(unlist(tmp_exonintron$genebody))[j] == "exons"){
          tmp_exonintron[j, ncol(tmp_exonintron)] <- paste("exons", numexon, sep = "_")
          
        }else if(as.vector(unlist(tmp_exonintron$genebody))[j] == "introns"){
          tmp_exonintron[j, ncol(tmp_exonintron)] <- paste("introns", numintron, sep = "_")
        }
        numexon <- 1
        numintron <- 1
      }else{
        
        # for the different gene #
        numexon <- 1
        numintron <- 1
        if(as.vector(unlist(tmp_exonintron$genebody))[j] == "exons"){
          tmp_exonintron[j, ncol(tmp_exonintron)] <- paste("exons", numexon, sep = "_")
          
        }else if(as.vector(unlist(tmp_exonintron$genebody))[j] == "introns"){
          tmp_exonintron[j, ncol(tmp_exonintron)] <- paste("introns", numintron, sep = "_")
        }
      }
    }
  }
  
  # delete the first row, last row and last column #
  tmp_exonintron$genebody <- tmp_exonintron[, ncol(tmp_exonintron)]
  tmp_exonintron <- tmp_exonintron[-c(1, nrow(tmp_exonintron)), -ncol(tmp_exonintron)]
  geneobj <- rbind(tmp_exonintron, tmp_noexonintron[, -ncol(tmp_noexonintron)])
  
  # sort the file #
  geneobj <- arrange(geneobj, chr, start, refseq)
  
  return(geneobj)
}





#' Read the cyto file.
#'
#' @description This function reads the chromosome information from cyto file (cytoBandIdeo.txt) and sort them by chromosome and position. 
#' The dataset of the example is the mouse genome information downloaded from UCSC website 
#' (http://hgdownload.cse.ucsc.edu/goldenPath/mm10/database/cytoBandIdeo.txt.gz). 
#' 
#' @param paths refers to the path of input file, with default the package path.
#' @param cytofile refers to the name of input cyto file that is downloaded from UCSC website, with default "cytoBandIdeo".
#' @param suffix refers to the suffix of input cyto file, e.g., ".gz", ".zip" and so on (some files are in text .txt format, then ".txt" or ".txt.gz"), with default ".txt.gz".
#'
#' @return Outputs a data frame contains chromosome, start position, end position.
#'  
#' @examples
#' inputcytofile <- Cytofile_read()
#' inputcytofile <- Cytofile_read(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""), cytofile = "cytoBandIdeo", suffix = ".txt.gz")
#' 
#' @export


Cytofile_read <- function(paths = paste(system.file(package = "GeneDMRs"), "/methdata", sep=""), cytofile = "cytoBandIdeo", suffix = ".txt.gz"){
  
  # set the paths #
  setwd(paths)
  
  #chromosome info#
  cyto <- read.table(paste(cytofile, suffix, sep = ""), fill = T, header=F)
  
  # five columns or four columns #
  if(ncol(cyto)==4){
    cyto <- cbind(cyto,"gneg")
	cyto[, 4] <- "p"
  }
  
  colnames(cyto) <- c("chr","start","end","band","stain")
  
  
  # find the unannotated chromosome rows and delete them #
  unqualifiedrow1 <- grep("_", cyto$chr)
  
  if(length(unqualifiedrow1) > 0){
    cyto <- cyto[-unqualifiedrow1,]
  }
  
  # find the chromosome M and delete it #
  unqualifiedrow2 <- grep("M", cyto$chr)
  
  if(length(unqualifiedrow2) > 0){
    cyto <- cyto[-unqualifiedrow2,]
  }
  
  # transfer chr column to character #
  cyto$chr <- as.vector(unlist(cyto$chr))
  
  # set the chromosome label to number as the default of qqman #
  chromtable <- table(cyto$chr)
  
  # set the chromosome label to number as the default of qqman #
  for(i in 1:length(chromtable[chromtable != 0])){
    
    cyto$chr[cyto$chr == paste("chr", i ,sep="")] <- i
  }
  cyto$chr[cyto$chr == "chrX"] <- length(chromtable[chromtable != 0]) - 1
  cyto$chr[cyto$chr == "chrY"] <- length(chromtable[chromtable != 0])
  
  # transfer chr column to numeric again for sort #
  cyto$chr <- as.numeric(cyto$chr)
  cyto <- cyto[order(cyto$chr), ] 
  
  # paste "chr" to chromosome number #
  cyto$chr <- paste("chr", cyto$chr, sep = "")
  cyto$chr[cyto$chr == paste("chr", (length(chromtable[chromtable != 0]) - 1), sep = "")] <- "chrX"
  cyto$chr[cyto$chr == paste("chr", length(chromtable[chromtable != 0]), sep = "")] <- "chrY"
  
  return(cyto)
}




