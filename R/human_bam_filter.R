library(optparse)

parser <- OptionParser()
option_list <- list(
    make_option(c("-c", "--chromsomeCount"), action="store", type="character", help="Path to <sample>_Count_chromosomes.txt"),
    make_option(c("-o", "--outdir"), action="store", type="character", help="Path to out directory of Viral BAM files"),
    make_option(c("-m", "--minreads", action="store", type="integer", default = 50, help="Minimum number of mapped viral reads, default = [default]")
)
    


opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE) )
chromosome_count_path <- opt$chromsomeCount
outdir <- opt$outdir
minreads <- opt$minreads

temp_chromosome_count_2 = read.table(chromosome_count_path,header = F,row.names = 1)colnames(temp_chromosome_count_2) = c("Chromosome_length","Mapped_reads","Unknown")

Human_chr = c("X","Y","MT",as.character(1:23))
temp_chromosome_count_2 = temp_chromosome_count_2[rownames(temp_chromosome_count_2)%in%Human_chr,]
temp_chromosome_count_2 = temp_chromosome_count_2[temp_chromosome_count_2$Mapped_reads>Minimal_read_mapped,]

write.csv(temp_chromosome_count_2, paste(outdir, "/human_chromosomes_count_filtered.csv", row.names = TRUE)

# Create empty file with virus name as file name 
for(chromosome in rownames(temp_chromosome_count_2)){
    file_name = paste(outdir, "/", chromosome, ".txt", sep = "")
    file.create(file_name)
}
