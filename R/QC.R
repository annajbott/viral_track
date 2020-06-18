library(optparse)


parser <- OptionParser()
option_list <- list(
    make_option(c("-v", "--viraldir"), action="store", type="character", help="Path to vViral_BAM_directory"),
    make_option(c("-o", "--outfile"), action="store", type="character", help="Filename of QC_filtered outfile"),
    make_option(c("-r", "--rRoot"), action ="store", type="character", help="Path to R directory in pipeline_viral_track repo")
)
    


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE) )
viraldir_path <- opt$viraldir
outfile <- opt$outfile
rdir <- opt$rRoot 

# Load in auxillary functions
source(paste(rdir, "/Module_viral_track.R", sep = "")) 

location_temp_chromsome_count <- paste(viraldir_path, "/virus_file_names/Virus_chromosomes_count_filtered.csv", sep ="")
temp_chromosome_count <- read.csv(location_temp_chromsome_count, stringsAsFactors=FALSE, header = TRUE, row.names = 1)

# Load in QC stuff, need to adjust for file and folder names
# And sort out ^M most likely

QC_result = foreach(i=rownames(temp_chromosome_count),.combine = rbind,.packages = c("GenomicAlignments","ShortRead")) %dopar% {
BAM_file= readGAlignments(paste(viraldir_path,"/",i,".bam",sep = ""),param = ScanBamParam(what =scanBamWhat()))

# Need to copy rest over but cba ...
