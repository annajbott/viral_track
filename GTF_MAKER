##Making GTF for viral track from Choromosome Index: 
library(optparse)

parser <- OptionParser()
option_list <- list( 
  make_option(c("-o", "--outputdir"), action="store", default='/well/immune-rep/users/kvi236/VIRUS/TEST', type="character", help="Path to output directory"),
  make_option(c("-i", "--indexgenome"), action="store", type="character", default="/gpfs2/well/immune-rep/users/kvi236/VIRAL_TRACK_REFERENCE_BUILD_273a", help="Path to VIRAL TRACK reference genome [default]")
)
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser, print_help_and_exit = TRUE, args = commandArgs(trailingOnly = TRUE) )

## Reading in chrLength file from the STAR genome Index build: 
backbone <- paste0(opt$indexgenome, 'chrNameLength.txt')
ChrLength <- read.delim(backbone, header=FALSE)

## Making GTF directory in outs 
GTF_dir = paste(opt$outputdir, '/GTF', sep = "")
dir.create(GTF_dir)
Path_GTF = paste(GTF_dir, "/ViralTrack_GTF.gtf", sep="")
fileConn <- file(Path_GTF, open = "w+")

#Loading the viral annotation file 
Viral_annotation = read.delim(Viral_annotation_file)
Viral_annotation = Viral_annotation[Viral_annotation$Name_sequence!=" ",]

for(i in 1:length(ChrLength[,1])){
  z = paste(ChrLength$V1[i],'RefSeq', 'transcript', 1, ChrLength$V2[i],1000,".",".", paste("gene_id ", ChrLength$V1[i] ,"_1;", sep =""), sep='\t')
  cat(z, file=fileConn, sep="\n", append=TRUE)
}

close(fileConn)
