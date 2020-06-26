## ---------------------------
## Script name: VIRAL TRACK: MATRIX MAKER
## Function: Add-On from ViralTrack (@Pierre)
##           Generate the ViralTrack demultiplexed MATRIX FILES from the output of UMI-tools countin viral_filtered_features directory.
##           The output will be: barcodes.tsv, a genomes.tsv and a viral_counts.mtx in sparse matrix form: similar to output of cell ranger.
##           viral_counts.mtx will include Human demultiplexted raw counts and viral demultiplexed raw counts.
##           viral_counts.mtx will also include % human reads, %mitochondrial reads, % total viral reads and % of reads/virus passng QC thresholds.
## 
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

library(optparse)
library(Matrix)
parser <- OptionParser()
option_list <- list( 
  make_option(c("-e", "--expressiontsv"), action="store", type="character", help="Path to Expression_table.tsv from umicounts"),
  make_option(c("-s", "--samplename"), action="store", type="character", help="Sample name"),
  make_option(c("-o", "--outdir"), action="store", type="character", help="Directory of outfile")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE))

## Create the Final Output directy: viral_filtered_features
MTX_dir = opt$outdir
dir.create(MTX_dir)

#Read non-sparse  Expression_table.tsv features from UMI-tools:
file_tsv = opt$expressiontsv
expression_table <- read.table(file = file_tsv, sep = '\t', header = TRUE, row.names = 1)

#Convert to data-frame to allow for Calculation of Summary Statistics
expression_table <- as.data.frame(t(expression_table))
viral_cols <- grep("refseq", colnames(expression_table), value=TRUE)
human_cols <- colnames(expression_table)[!colnames(expression_table) %in% viral_cols]
mito_cols <- grep("MT_1", colnames(expression_table), value=TRUE)

##Calculating meta data percentages:
# Percent human reads per cell:
expression_table$percent_human_reads <- (rowSums(expression_table[, c(human_cols)]))/ (rowSums(expression_table))* 100
# Percent mito reads per cell:
expression_table$percent_mitochondrial_reads <- (expression_table[, c(mito_cols)])/ (rowSums(expression_table))* 100

#Percent viral read per cell - for all viruses:
if(length(viral_cols)>= 1){
  cat("Viruses PRESENT. Calculating group and individual viral percentages")
  if(length(viral_cols)==1){
    cat("one virus present")
    expression_table$percent_viral_reads <- expression_table[, c(viral_cols)]/ (rowSums(expression_table[, c(viral_cols, human_cols)]))* 100
  } else {
    cat("Multiple Viruses Present")
    expression_table$percent_viral_reads <- (rowSums(expression_table[, c(viral_cols)]))/ (rowSums(expression_table[, c(viral_cols, human_cols)]))* 100  
  } 
} else {
  cat("No viral reads present: percent viral = 0 ") 
  expression_table$percent_viral_reads <- "0"
}

# Percent viral reads per cell per virus.
if(length(viral_cols)>= 1){
  cat("Calculating Percentage per Virus")
  for (i in viral_cols){
    y <- strsplit(i, split='|', fixed=TRUE)
    y = unlist(y)
    name=paste0("percent_viral_reads_", y[2], y[4])
    virus <- i
    expression_table[, name] <- (expression_table[, c(virus)])/ (rowSums(expression_table[, c(viral_cols, human_cols)]))* 100
  }
} 

## All percentages calculated
## Converting Table To SPARSE MATRIX FORMAT
expression_table <- (t(expression_table))
expression_table <- Matrix(expression_table, sparse = TRUE)
writeMM(obj = expression_table, file=paste0(MTX_dir, "viral_counts.mtx"))
cellbarcode <- rownames(expression_table)
genome <- colnames(expression_table)
write.table(cellbarcode, file=paste0(MTX_dir, 'barcodes.tsv'), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
write.table(genome, file=paste0(MTX_dir, 'genomes.tsv'), quote=FALSE, sep='\t', col.names = FALSE, row.names = FALSE)
## Output files generated
