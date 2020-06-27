library(optparse)
library(doParallel)
library(Biostrings)
library(ShortRead)
library(GenomicAlignments)

parser <- OptionParser()
option_list <- list(
    make_option(c("-v", "--viraldir"), action="store", type="character", help="Path to vViral_BAM_directory"),
    make_option(c("-o", "--outfile"), action="store", type="character", help="Path to filtered outfile"),
    make_option(c("-s", "--sample"), action = "store", type="character", help="Sample name"),
    make_option(c("-r", "--rRoot"), action ="store", type="character", help="Path to R directory in pipeline_viral_track repo"),
    make_option(c("-m", "--minReadsMapped"), action ="store", type="character", help="Minimal reads mapped param"),
    make_option(c("-a", "--viralannotation"), action="store", type="character", default="/ifs/projects/adam/VIRAL_TRACK/References/Virusite_annotation_file.txt", help="Path to VirusSite annotation file [default]")
)
    

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE) )
viraldir_path <- opt$viraldir   
outfile <- opt$outfile
sample_name <- opt$sample
rdir <- opt$rRoot
Minimal_read_mapped <- opt$minReadsMapped
Viral_annotation_file <- opt$viralannotation
Virus_database <- read.delim(Viral_annotation_file, header=T, sep="\t")

filter_outfile <- gsub("QC_filtered.pdf", "QC_filtered.txt", outfile)
unfilter_outfile <- gsub("QC_filtered.pdf", "QC_unfiltered.txt", outfile)
log <- paste0(unfilter_outfile, ".log")

pdf_name <- outfile
temp_chromosome_count_path <- paste("STAR.dir/", sample_name, "/", sample_name, "_Count_chromosomes.txt",sep="")
path_to_Log_file <- gsub("QC_filtered.txt", "Sample_Log.final.out", unfilter_outfile)
path_to_Log_file <- paste0("STAR.dir/", sample_name, "/", sample_name, "_Log.final.out")

# Load in auxillary functions
source(paste(rdir, "/Module_viral_track.R", sep = "")) 

location_temp_chromsome_count <- paste(viraldir_path, "/virus_file_names/Virus_chromosomes_count_filtered.csv", sep ="")
temp_chromosome_count <- read.csv(location_temp_chromsome_count, stringsAsFactors=FALSE, header = TRUE, row.names = 1)
if (length(temp_chromosome_count/2) <= 6){
	N_thread <- length(temp_chromosome_count/2)
	} else {
	N_thread == 6
	}
	
##Make Parallelel Environment : Number of threads calculated on the number of unique viruses present 
cl =makeCluster(N_thread)
registerDoParallel(cl)

dir <- paste(viraldir_path,sep = "")
if(length(list.files(dir))==0){
  cat("\t ----------------------------------------------\n", file=log, append=TRUE)
  cat(paste0("\t NO VIRUSES IDENFIFIED WITH MIN READS MAPPED >=", Minimal_read_mapped, ". \n"), file=log, append=TRUE)
  cat("\t NO QC VIRAL METRICS WILL BE CALCULTED. \n", file=log, append=TRUE)
  cat("\t ----------------------------------------------\n", file=log, append=TRUE)
  }  

if(length(list.files(dir))>0){
  cat("\t ----------------------------------------------\n", file=log, append=TRUE)
  cat(paste0("\t ", length(list.files(dir)), " VIRUSES IDENFIFIED WITH MIN READS MAPPED >=", Minimal_read_mapped, ". \n"), file=log, append=TRUE)
  cat("\t QC VIRAL METRICS WILL BE CALCULTED. \n", file=log, append=TRUE)
  cat("\t ----------------------------------------------\n", file=log, append=TRUE)
  QC_result <- NULL
## Generation of the QC repor
  for(i in 1:length(rownames(temp_chromosome_count))) {
    z <- rownames(temp_chromosome_count)[i]
    z <- gsub("\\|", "-", z)
    BAM_file= readGAlignments(paste(viraldir_path,"/",z,".bam",sep = ""),param = ScanBamParam(what =scanBamWhat()))
    #Let's check the diversity of the reads
    Viral_reads = unique(BAM_file@elementMetadata$seq)
    Viral_reads_contents = alphabetFrequency(Viral_reads,as.prob =T )
    Viral_reads_contents = Viral_reads_contents[,c("A","C","G","T")]
    
    #calculate percentage of nucleotides in the viral reads for Read Entropy
    # If only one reads present returns numeric vector rather than dataframe:
    if (class(Viral_reads_contents)=="numeric") {
      Viral_reads_contents = matrix(Viral_reads_contents_mean,ncol = 4)
    }
    
    Viral_reads_contents_mean =colMeans(Viral_reads_contents)
    Read_entropy = sum(-log(Viral_reads_contents_mean)*Viral_reads_contents_mean,na.rm = T)
    #Caclulate the spatial distribution of the mapped reads : how much percent of the genome is mapped ?
    Covered_genome = coverage(BAM_file)[[z]]
    Covered_genome = as.numeric(Covered_genome)
    Spatial_distribution =sum(Covered_genome>0)/length(Covered_genome)
    Covered_genome = rle(sign(Covered_genome))
    Longest_contig = max(Covered_genome$lengths[Covered_genome$values>0])
    
    ##The mean reads quality
    Reads_quality = as.character(BAM_file@elementMetadata$qual)
    Reads_quality = PhredQuality(Reads_quality)
    Reads_quality = as(Reads_quality,"IntegerList")
    Reads_quality = as.numeric(as.matrix(Reads_quality))
    Mean_read_quality = mean(Reads_quality)
    Sd_read_quality = sd(Reads_quality)
    
    ##... the number of mapped reads and unique mapped reads
    N_unique_mapped_reads = sum(BAM_file@elementMetadata$mapq==255) ##Code specific to STAR aligner.... 
    N_mapped_reads = length(BAM_file)
    Percent_uniquely_mapped = N_unique_mapped_reads/N_mapped_reads
    
    ##DUSTy score identifies low-complexity sequences, in a manner inspired by the dust implementation in BLAST
    Mean_dust_score = NA
    Percent_high_quality_reads = NA
    if ("ShortRead"%in%installed.packages()){
      DUST_score = dustyScore(BAM_file@elementMetadata$seq)
      Mean_dust_score = mean(DUST_score)
      Percent_high_quality_reads =  sum(DUST_score<500)/length(DUST_score)
    }
    
    ##Summary Statistics Per Virus 
    QC_temp = c(N_mapped_reads,N_unique_mapped_reads,Percent_uniquely_mapped,
                Mean_read_quality,Sd_read_quality,
                Viral_reads_contents_mean,Read_entropy,Spatial_distribution,Longest_contig,
                Mean_dust_score,Percent_high_quality_reads)
    QC_result = rbind(QC_result, QC_temp)
  }
} else {
  QC_result = data.frame(N_mapped_reads = numeric(),N_unique_mapped_reads=numeric(),Percent_uniquely_mapped=numeric(),
              Mean_read_quality=numeric(),Sd_read_quality=numeric(),
              A = numeric(), C= numeric(), G = numeric(), T = numeric(), Read_entropy=numeric(),Spatial_distribution=numeric(),Longest_contig=numeric(),
              Mean_dust_score=numeric(),Percent_high_quality_reads=numeric())
}

# Load in QC stuff, need to adjust for file and folder names
# And sort out  most likely
QC_result = foreach(i=gsub("\\|", "-", rownames(temp_chromosome_count)),.combine = rbind,.packages = c("GenomicAlignments","ShortRead")) %dopar% {
	
	
	BAM_file= readGAlignments(paste(viraldir_path,"/",i,".bam",sep = ""),param = ScanBamParam(what =scanBamWhat()))
 	a = gsub("-", "|", i)
	#Let's check the diversity of the reads
	Viral_reads = unique(BAM_file@elementMetadata$seq)
	Viral_reads_contents = alphabetFrequency(Viral_reads,as.prob =T )
	Viral_reads_contents = Viral_reads_contents[,c("A","C","G","T")]

	#calculate percentage of nucleotides in the viral reads for Read Entropy
	Viral_reads_contents_mean =colMeans(Viral_reads_contents)
	Read_entropy = sum(-log(Viral_reads_contents_mean)*Viral_reads_contents_mean,na.rm = T)

	#Caclulate the spatial distribution of the mapped reads : how much percent of the genome is mapped ?
	Covered_genome = coverage(BAM_file)[[a]]
	Covered_genome = as.numeric(Covered_genome)
	Spatial_distribution =sum(Covered_genome>0)/length(Covered_genome)
	Covered_genome = rle(sign(Covered_genome))
	Longest_contig = max(Covered_genome$lengths[Covered_genome$values>0])
  
	##The mean reads quality
	Reads_quality = as.character(BAM_file@elementMetadata$qual)
	Reads_quality = PhredQuality(Reads_quality)
	Reads_quality = as(Reads_quality,"IntegerList")
	Reads_quality = as.numeric(as.matrix(Reads_quality))
	Mean_read_quality = mean(Reads_quality)
	Sd_read_quality = sd(Reads_quality)
  
	##... the number of mapped reads and unique mapped reads
	N_unique_mapped_reads = sum(BAM_file@elementMetadata$mapq==255) ##Code specific to STAR aligner.... 
	N_mapped_reads = length(BAM_file)
	Percent_uniquely_mapped = N_unique_mapped_reads/N_mapped_reads
  
  
	##DUSTy score identifies low-complexity sequences, in a manner inspired by the dust implementation in BLAST
	Mean_dust_score = NA
	Percent_high_quality_reads = NA
	if ("ShortRead"%in%installed.packages()){
		DUST_score = dustyScore(BAM_file@elementMetadata$seq)
		Mean_dust_score = mean(DUST_score)
		Percent_high_quality_reads =  sum(DUST_score<500)/length(DUST_score)
	}
  
	##Summary Statistics Per Virus 
	QC_temp = c(N_mapped_reads,N_unique_mapped_reads,Percent_uniquely_mapped,
              Mean_read_quality,Sd_read_quality,
              Viral_reads_contents_mean,Read_entropy,Spatial_distribution,Longest_contig,
              Mean_dust_score,Percent_high_quality_reads)
    }


## ------------------------------------------------------------------------------------
## Now we perform filtering based ont the Calculaed QC statsitics: 
## Editied by LEO -- otherwise fails if only one virus is detected as produces numeric vector rather than dataframe. 
if (class(QC_result)=="numeric"){
  cat('Only one virus detected - output is not a dataframe: converting \n', file="Viral_track_scanning.log", append = TRUE)
  QC_result <- as.data.frame(t(as.data.frame(QC_result)))
}
colnames(QC_result) = c("N_reads","N_unique_reads","Percent_uniquely_mapped",
                        "Mean_read_quality","Sd_read_quality",
                        c("A","C","G","T"),"Sequence_entropy","Spatial_distribution","Longest_contig",
                        "DUST_score","Percent_high_quality_reads")
QC_result <- as.data.frame(QC_result)
rownames(QC_result) = rownames(temp_chromosome_count)
QC_result = QC_result[QC_result$N_unique_reads>0,]

### ------------------------------------------------------------------------------------
## Now we need to extract information on the mapping by itself to get information about the QC
## This Requires R auxillalry Functions: 
Mapping_information = Extraction_Log_final(path_to_Log_file)
Mean_mapping_length = Mapping_information$Length_vector[1]


detected_virus = rownames(QC_result[QC_result$Sequence_entropy>1.2 & QC_result$Longest_contig>3*Mean_mapping_length & QC_result$Spatial_distribution>0.05,])

## Returning QC Metrics for viruses that passed Filtering 
Filtered_QC=QC_result[detected_virus,]



##Exporting the tables of the QC analysis
write.table(file = unfilter_outfile, x = QC_result, quote = F,sep = "\t")
write.table(file = filter_outfile, x = Filtered_QC, quote = F,sep = "\t")



## ------------------------------------------------------------------------------------
# Note removed splice table info as this didnt seem to be used anywhere else?
##------------------------------------------------------------------------------------

## Calculating some additional Statistics needed for Plots: 

## Additional info : % of reads mapped to viral vs host
Read_count_temp = read.table(temp_chromosome_count_path,header = F,row.names = 1)
colnames(Read_count_temp) = c("Chromosome_length","Mapped_reads","Unknown")
Read_count_temp = Read_count_temp[Read_count_temp$Mapped_reads!=0,]
Read_count_temp$Chr <- rownames(Read_count_temp)

## Need to rename human chr to append chr to make it easy to grep:
for (x in 1:22){
  Read_count_temp[["Chr"]] <- with(Read_count_temp, ifelse(Chr == x , paste0("chr", Chr), Chr))
}

Read_count_temp[["Chr"]] <- with(Read_count_temp, ifelse(Chr == "X", paste0("chr", Chr), Chr))
Read_count_temp[["Chr"]] <- with(Read_count_temp, ifelse(Chr == "Y", paste0("chr", Chr), Chr))
Read_count_temp[["Chr"]] <- with(Read_count_temp, ifelse(Chr == "MT", paste0("chr", Chr), Chr))
rownames(Read_count_temp) <- Read_count_temp$Chr
Read_count_temp$Chr <- NULL

## Calculate the percentage of human reads etc 
host_mapping_count = sum(Read_count_temp[grepl(pattern = "chr",rownames(Read_count_temp)),"Mapped_reads"])
viral_mapping_count = sum(Read_count_temp[grepl(pattern = "NC",rownames(Read_count_temp)),"Mapped_reads"])
total_mapping = viral_mapping_count + host_mapping_count
Ratio_host_virus = matrix(data = c(host_mapping_count,viral_mapping_count)/total_mapping,ncol = 1)*100
percent_viral = viral_mapping_count / total_mapping * 100 

## Number of uniquely mapped reads and other reads for the filtered virus 
Mapping_selected_virus = data.frame(Unique_mapping = (Filtered_QC$N_unique_reads),All_mapping = (Filtered_QC$N_reads),row.names = rownames(Filtered_QC))
Mapping_selected_virus = Mapping_selected_virus[order(Mapping_selected_virus$Unique_mapping,decreasing = T),]

## ------------------------------------------------------------------------------------
## Plotting the Statisitics 

Mapping_selected_virus$Name_id <- rownames(Mapping_selected_virus)

if(length(Mapping_selected_virus[, 1])>0) {
  Mapping_selected_virus$Name_sequence<- c()
  for (i in 1:length(Mapping_selected_virus[, 1])){
    z <- Mapping_selected_virus[i, 3]
    print(z)
    z <- unlist(strsplit(z,"|",fixed = T))[2]
    Mapping_selected_virus$Name_sequence[i] <- z
  } 
} else {
  Mapping_selected_virus <- data.frame(Unique_mapping = numeric(), All_mapping = numeric(), Name_id = character(), Name_sequence = character())
}
Mapping_selected_virus <- merge(Mapping_selected_virus, Virus_database, by="Name_sequence")
rownames(Mapping_selected_virus) <- Mapping_selected_virus$Complete_segment_name



## Open PDF file: 
pdf(pdf_name, height = 24, width = 20)
par(las=0,mfrow=c(4,3),mar=c(6,6,6,4))
Color_vector = c("lightskyblue1","orange","grey80")

# Plotting the proportion of uniquely mapped reas, unmapped etc...
barplot(Mapping_information$Mapping_result,ylim=c(0,100),xlim=c(0,5),ylab="Percentage of reads (%)",col=Color_vector,cex.lab=1.5,cex.axis = 1.5, main = "Percentage of Reads Mapped")
legend(x = 1.5,y=50,legend = c("Unmapped","Mapped to multiple loci","Uniquely mapped"),bty="n",fill = Color_vector[length(Color_vector):1],cex = 1.5)
# Size of mapping, insertion and deletion
barplot(Mapping_information$Length_vector,col="black",names.arg = c("Mapping Length","Insertion Length","Deletion length"), main = "Size of Mapped Reads",  horiz = T,xlim=c(0,max(Mapping_information$Length_vector[1])*1.2),xlab="Nucleotide length",cex.lab=1.3,cex.axis = 1.3,cex.names=1.3)
#Rate of mismatch, deletion and insertion
barplot(Mapping_information$Rate_vector,col="black",names.arg = c("Mismatch rate","Insertion rate","Deletion rate"), main = "Rate of Mismatching", horiz = T,xlim=c(0,max(Mapping_information$Rate_vector[1])*1.2),xlab="Rate (%)",cex.lab=1.3,cex.axis = 1.3,cex.names=1.3)
# Ratio 
Color_vector = c("darkred","grey")
barplot(Ratio_host_virus,ylim=c(0,100),xlim=c(0,5),ylab="Mapping events (%)",col=Color_vector,cex.lab=1.5,cex.axis = 1.5, main = "Ratio of Host:Viral Mapping")
legend(x = 1.5,y=50,legend = c("Viral mapping","Host mapping"),bty="n",fill = Color_vector[length(Color_vector):1],cex = 1.5)


# First QC for the viral hits 
if (length(detected_virus) == 0) {
  Color_vector= string.to.colors(factor(rownames(QC_result)%in%detected_virus),colors = c("orange")) ##Viral sequences that passed QC : green
}
if (length(detected_virus) > 0) {
  Color_vector= string.to.colors(factor(rownames(QC_result)%in%detected_virus),colors = c("orange","green")) ##Viral sequences that passed QC : green
}

# Plot Number of Reads > 50 (filtering threshold)
plot(QC_result$N_reads,QC_result$Spatial_distribution*100,pch=21,bg=Color_vector,log="x",cex=1.5,xlab="Number of Mapped Reads",ylab="% Mapped genome",ylim=c(0,100),cex.lab=1.5,main="Viral Summary: Read Count")
abline(h=5,lwd=2,lty=2,col="grey")
abline(v=opt$minReadsMapped,lwd=2,lty=2,col="grey")

# Plot Number of Unique Reads > 50 (filtering threshold)
plot(QC_result$N_unique_reads,QC_result$Spatial_distribution*100,pch=21,bg=Color_vector,log="x",cex=1.5,xlab="Number of Uniquely Mapped Reads",ylab="% Mapped genome",ylim=c(0,100),cex.lab=1.5,main="Viral Summary: Unique Read Count")
abline(h=5,lwd=2,lty=2,col="grey")

#Second QC for the viral hits 
plot(QC_result$Sequence_entropy,QC_result$Spatial_distribution*100,pch=21,bg=Color_vector, cex=1.5,xlab="Sequence Complexity",ylab="% Mapped genome",cex.lab=1.5,ylim=c(0,100),main="Viral Summary: Sequence Complexity")
abline(h=5,lwd=2,lty=2,col="grey")
abline(v=1.2,lwd=2,lty=2,col="grey")

#Third QC for the viral hits 
plot(QC_result$Longest_contig,QC_result$DUST_score,pch=21,bg=Color_vector, cex=1.5,xlab="Longest Contig (nt)",ylab="DUST Score",cex.lab=1.4,main="Viral Summary: DUST Score")
abline(v=3*Mean_mapping_length,lwd=2,lty=2,col="grey")


#Number of reads for each filtered virus
if (length(detected_virus) > 0) {
  par(las=0)
  barplot(Mapping_selected_virus$All_mapping[nrow(Mapping_selected_virus):1],
          col="black",horiz = T,cex.lab=1.3,cex.axis = 1.3,cex.names=1,
          xlim=c(0,max(Mapping_selected_virus$All_mapping)*1.2),xlab="Number of mapped reads",
          names.arg = rownames(Mapping_selected_virus)[nrow(Mapping_selected_virus):1])
  par(las=0)
  barplot(Mapping_selected_virus$Unique_mapping[nrow(Mapping_selected_virus):1],
          col="black",horiz = T,cex.lab=1.3,cex.axis = 1.3,cex.names=1,
          xlim=c(0,max(Mapping_selected_virus$Unique_mapping)*1.2),xlab="Number of uniquely mapped reads",
          names.arg = rownames(Mapping_selected_virus)[nrow(Mapping_selected_virus):1])
	}
dev.off()

## Filtering and PDF plots generated 

