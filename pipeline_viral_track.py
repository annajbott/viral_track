""" =====================
Pipeline viral track
=========================

Overview 
=========

This workflow was designed as a CGAT pipeline representation of R viral_track code by Bost et al.
Viral track can be found at:
- https://www.sciencedirect.com/science/article/pii/S0092867420305687
- https://github.com/PierreBSC/Viral-Track

Additonal functionality has been added. 

Adapted by Lauren Overend and Anna James-Bott.

Any ViralTrack specific questions direct to Lauren.overend@oriel.ox.ac.uk
Any CGAT pipeline specific questions direct to anna.james-bott@st-hildas.ox.ac.uk



Requires:

* Fastq files 
* Reference genome 


Pipeline output
================
<What comes out of pipeline...>


Code
=====

"""

from ruffus import *

import sys
import os

import cgatcore.pipeline as P
import cgatcore.experiment as E


###############################
# Load in parameters and data
###############################
    
# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])


# Find fastq files (unsure of format)
try:
    PARAMS['data']
except NameError:
    DATADIR = "."
else:
    if PARAMS['data'] == 0:
        DATADIR = "."
    elif PARAMS['data'] == 1:
        DATADIR = "data.dir"
    else:
        DATADIR = PARAMS['data']

SEQUENCESUFFIXES = ("*.fastq.gz",
                    "*.fastq.1.gz", "*.fa", "*.fastq")
SEQUENCEFILES = tuple([os.path.join(DATADIR, suffix_name)
                            for suffix_name in SEQUENCESUFFIXES])


#############################
# Make GTF for viral track
#############################

@follows(mkdir("geneset.dir"))
@originate("geneset.dir/ViralTrack_GTF.gtf")
def make_gtf(outfile):
    '''
    Make custom GTF for viraltrack from chromosome index 
    or copy existing gtf into geneset directory
    '''
    R_ROOT = os.path.join(os.path.dirname(__file__), "R")

    index_genome = PARAMS['index_genome_dir']
    viraltrack_loc = PARAMS['viraltrack_gtf']

    # If there is 
    if viraltrack_loc:
        statement = ''' cp %(viraltrack_loc)s %(outfile)s '''

    else:
        statement = ''' Rscript %(R_ROOT)s/GTF_maker.R -i %(index_genome)s -o %(outfile)s '''

    P.run(statement)


#############################
# Build indexes
#############################

# Build star index genome if index_star set to 1 in params.yml
@active_if(PARAMS['index_star_wget'])
@follows(mkdir("index.dir"))
@follows(mkdir("index.dir/HUMAN_GENOME"))
@originate(["index.dir/HUMAN_GENOME/Homo_sapiens.GRCh38.dna.chromosome.Y.fa"])
def wget_human_chromosomes(outfile):
    ''' 
    Download human chromosomes 1-22, X and Y
    Wget ensembl public fastas then gunzip the files.
    '''

    statement = ''' wget --directory-prefix=index.dir/HUMAN_GENOME/ -r -np -nH -nd 
                -A "Homo_sapiens.GRCh38.dna.chromosome.[0-9].fa.gz" 
                ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/ && 
                wget --directory-prefix=index.dir/HUMAN_GENOME/ -r -np -nH -nd 
                -A "Homo_sapiens.GRCh38.dna.chromosome.[12][0-9].fa.gz" 
                ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/ &&
                wget --directory-prefix=index.dir/HUMAN_GENOME/ -r -np -nH -nd 
                -A "Homo_sapiens.GRCh38.dna.chromosome.[XY].fa.gz" 
                ftp://ftp.ensembl.org/pub/release-100/fasta/homo_sapiens/dna/ &&
                gunzip index.dir/HUMAN_GENOME/*.fa.gz && sleep 120
                '''
    
    P.run(statement, to_cluster=False)


@active_if(PARAMS['index_star_wget'])
@follows(mkdir("index.dir/VIRUS_GENOME"))
@originate("index.dir/VIRUS_GENOME/genes.fasta")
def wget_viral_chromosomes(outfile):
    '''
    Download viral reference from VirusSite and unzip
    '''



    statement = ''' wget --directory-prefix=index.dir/VIRUS_GENOME/ 
                http://www.virusite.org/archive/2020.2/genes.fasta.zip && 
                unzip -d index.dir/VIRUS_GENOME index.dir/VIRUS_GENOME/genes.fasta.zip
                '''

    P.run(statement)

@active_if(PARAMS['index_star'])
@follows(mkdir("index.dir/star_index.dir"))
@originate("index.dir/star_index.dir/Genome")
def STAR_index(outfile):

    '''
    Builds star index for human genome and viruses
    '''

    nthreads = PARAMS['star_nThreads']
    # In yml file contain param for any extra fastas (file location), e.g. covid19, "" if none
    extra_fasta = PARAMS['index_extra_fasta']

    statement = ''' STAR --runThreadN %(nthreads)s --limitGenomeGenerateRAM 268822031285 --runMode genomeGenerate 
                --genomeDir index.dir/star_index.dir --genomeFastaFiles 
                index.dir/VIRUS_GENOME/genes.fasta index.dir/HUMAN_GENOME/*.fa 
                %(extra_fasta)s
                '''

    job_memory="100G"

    P.run(statement)



#############################
# STAR mapping and samtools
#############################

@follows(STAR_index)
@follows(mkdir("STAR.dir"))
@transform(SEQUENCEFILES,
           regex("(\S+).fastq.gz"),
           r"STAR.dir/\1/\1_Aligned.sortedByCoord.out.bam")
def STAR_map(infile, outfile):
    '''
    Run STAR mapping with parameters defined in yml
    '''
    
    if PARAMS['index_star']:
        index_genome = "index.dir/star_index.dir"
    else: 
        index_genome = PARAMS['index_genome_dir']
            

    nthreads = PARAMS['star_nThreads']
    nBAMsortingthreads = PARAMS['star_nThreadsort']
    nBins = PARAMS['star_bins']
    min_reads = PARAMS['star_minreads']
    prefix = outfile.replace("Aligned.sortedByCoord.out.bam", "")
    log_file = outfile.replace(".bam", ".log")
    
    if infile.endswith(".gz"):
        gunzip = "--readFilesCommand zcat"
    else:
        gunzip = ""


    if int(PARAMS['star_nThreadsort']) > int(PARAMS['star_nThreads']):
        E.warn("nThreadsort should not be greater than the number of threads for STAR mapping")

    

    statement = ''' STAR --runThreadN %(nthreads)s --outBAMsortingThreadN %(nBAMsortingthreads)s --outBAMsortingBinsN %(nBins)s 
                --genomeDir %(index_genome)s --readFilesIn %(infile)s --outSAMattributes NH HI AS nM NM XS 
                --outSAMtype BAM SortedByCoordinate --twopassMode Basic --outFilterMatchNmin 35 
                --outFilterScoreMinOverLread 0.6 --outFilterMatchNminOverLread 0.6 --outFileNamePrefix %(prefix)s %(gunzip)s > %(log_file)s
                '''

    job_memory = "70G"

    P.run(statement)


### Merge samples here ###
## samtools merge, samtools sort, samtools index


@transform(STAR_map,
           regex("STAR.dir/(\S+)/(\S+)_Aligned.sortedByCoord.out.bam"),
           r"STAR.dir/\1/\1_Aligned.sortedByCoord.out.bam.bai")
def samtools_index(infile, outfile):
    '''
    Index bam file using samtools
    '''

    #log_file = outfile.replace(".bam",".log")

    statement = '''samtools index %(infile)s '''

    job_memory = "20G"

    P.run(statement)


@transform(STAR_map,
           regex("(\S+)/(\S+)_Aligned.sortedByCoord.out.bam"),
           r"\1/\2_Count_chromosomes.txt")
def samtools_chromosome_count(infile, outfile):
    ''' 
    Compute the number of mapped reads for each chromosome/virus
    '''

    statement = '''samtools idxstats %(infile)s > %(outfile)s'''

    job_memory = "20G"

    P.run(statement)

################################
# Filter out human and viruses
################################

@subdivide(samtools_chromosome_count, 
           regex("(\S+)/(\S+)_Count_chromosomes.txt"),
           r"\1/Viral_BAM_files/virus_file_names/*.txt")
def viral_filter(infile, outfile):
    '''
    Uses R script to filter out human chromosomes and viruses
    under threshold number of reads, then create empty txt file for each virus 
    '''
    chromosome_count = infile
    BAM_folder = infile.replace(os.path.basename(infile), "") + "/Viral_BAM_files/"
    outdir = infile.replace(os.path.basename(infile), "") + "/Viral_BAM_files/virus_file_names"
    # for some reason ruffus complains that the output isnt a string - this is a tmp fix
    outfile = str(outfile)
    minreads = PARAMS['star_minreads']

    if os.path.exists(BAM_folder):
        pass
    else:
        os.mkdir(BAM_folder)
        os.mkdir(outdir)

    R_ROOT = os.path.join(os.path.dirname(__file__), "R")

    statement = '''Rscript %(R_ROOT)s/viral_bam_filter.R -c %(chromosome_count)s 
                -m %(minreads)s -o %(outdir)s'''
    
    P.run(statement)


@follows(samtools_index)
@transform(viral_filter,
           regex("STAR.dir/(\S+)/Viral_BAM_files/virus_file_names/(\S+).txt"),
           add_inputs([STAR_map]),
           r"STAR.dir/\1/Viral_BAM_files/\2.bam")
def viral_BAM(infiles, outfile):
    ''' 
    Takes virus names from empty text files and aligned bam 
    and makes bam file for each virus
    '''

    virus_name_file, aligned_bam = infiles
    name_inf = os.path.split(virus_name_file)[0].split("/")[2]
    for bam in aligned_bam:
        name = os.path.split(bam)[0].split("/")[2]
        if name == name_inf:
            correct_bam = bam
        else:
            pass

    virus =  os.path.basename(virus_name_file).replace(".txt", "")
    virus = virus.replace("-","|")

    # The output file contains | which need to be escaped or maybe remove the | in the output files


    statement = """ samtools view -b %(correct_bam)s '%(virus)s' > %(outfile)s """

    P.run(statement)


@subdivide(samtools_chromosome_count, 
           regex("(\S+)/(\S+)_Count_chromosomes.txt"),
           r"\1/Human_BAM_files/human_file_names/*.txt")
def human_filter(infile, outfile):
    '''
    Uses R script to filter out anything not human and any chromosomes
    under threshold number of reads, then create empty txt file for each name 
    '''

    chromosome_count = infile
    BAM_folder = infile.replace(os.path.basename(infile), "") + "/Human_BAM_files/"
    outdir = infile.replace(os.path.basename(infile), "") + "/Human_BAM_files/human_file_names"
    # for some reason ruffus complains that the output isnt a string - this is a tmp fix
    outfile = str(outfile)
    minreads = PARAMS['star_minreads']

    if os.path.exists(BAM_folder):
        pass
    else:
        os.mkdir(BAM_folder)
        os.mkdir(outdir)

    R_ROOT = os.path.join(os.path.dirname(__file__), "R")

    statement = '''Rscript %(R_ROOT)s/human_bam_filter.R -c %(chromosome_count)s 
                -m %(minreads)s -o %(outdir)s '''
    
    P.run(statement)


@transform(human_filter,
           regex("STAR.dir/(\S+)/Human_BAM_files/human_file_names/(\S+).txt"),
           add_inputs([STAR_map]),
           r"STAR.dir/\1/Human_BAM_files/\2.bam")
def human_BAM(infiles, outfile):
    ''' 
    Takes human chromosome names from empty text files and aligned bam 
    and makes bam file for each chromsome
    '''

    human_name_file, aligned_bam = infiles
    human_chrom =  os.path.basename(human_name_file).replace(".txt", "")
    
    name_inf = os.path.split(human_name_file)[0].split("/")[2]
    for bam in aligned_bam:
        name = os.path.split(bam)[0].split("/")[2]
        if name == name_inf:
            correct_bam = bam
        else:
            pass

    statement = """ samtools view -b %(correct_bam)s '%(human_chrom)s' > %(outfile)s """

    P.run(statement)

#######################
# Quality control
#######################

# Need to finish R script ...
@collate(viral_BAM,
        regex("STAR.dir/(\S+)/Viral_BAM_files/\S+.bam"),
        r"QC.dir/\1/QC_filtered.pdf") # Fill in name
def viral_QC(infile, outfile):
    ''' 
    Performs quality control on viral BAM file
    '''
    R_ROOT = os.path.join(os.path.dirname(__file__), "R")

    viral_bam_directory = os.path.dirname(infile[0])
    sample = viral_bam_directory.split("/")[2]
    min_reads_mapped = PARAMS['min_reads_mapped']
    
    statement = '''Rscript %(R_ROOT)s/QC.R --viraldir %(viral_bam_directory)s -o %(outfile)s -r %(R_ROOT)s -s %(sample)s -m %(min_reads_mapped)s'''


    P.run(statement)
    
# Not finished ....
@collate(viral_BAM,
         regex("STAR.dir/(\S+)/Viral_BAM_files/\S+.bam"),
         add_inputs(viral_QC),
         r"STAR.dir/(\S+)/merged_viral_mapping.bam")
def merge_viruses(infiles, outfile):
    '''
    Merging BAM for all QC filtered viruses
    '''

#############################
# De-multiplexing R file
#############################

@collate(human_BAM,
         regex("STAR.dir/./(\S+)/./Human_BAM_files/\S+.bam"),
         r"STAR.dir/\1/merged_human_mapping.bam")
def samtools_merge(infile, outfile):
    '''
    Merge human chromsome BAM files using samtools
    '''

    path_to_BAM = os.path.dirname(infile[0])
    infile_list = ' '.join(infile)
    

    statement = ''' samtools merge %(outfile)s -f %(infile_list)s '''

    P.run(statement)


@follows(mkdir("featureCounts.dir"))
@transform(samtools_merge,
           regex("STAR.dir/(\S+)/merged_human_mapping.bam"),
           add_inputs(make_gtf),
           r"featureCounts.dir/\1/merged_human_mapping.featureCounts.bam")
def feature_counts(infile, outfile):
    '''
    Perform feature counts on the merged BAM file
    '''

    merged_bam, gtf = infile
    sample_name = merged_bam.split("/")[1]
    counts_file = "featureCounts.dir/" + sample_name + "/" + sample_name + "_counts.txt"
    # Move bam file to featurecounts folder and rename 
    bam_out_og = merged_bam
    
    statement = ''' featureCounts -t transcript -M --primary -R BAM -g gene_id 
                    -a  %(gtf)s -o %(counts_file)s %(merged_bam)s && 
                    cp "%(bam_out_og)s" "%(outfile)s" '''


    P.run(statement)


# What directory to put in ??? Have just put in feature counts
@transform(feature_counts,
           regex("(\S+)/(\S+)/(\S+).featureCounts.bam"),
           r"\1/\2/\2_assigned_sorted.bam")
def samtools_sort_index(infile, outfile):
    '''
    Sort and index BAM file produced by feature counts
    ''' 

    statement = ''' samtools sort %(infile)s -o %(outfile)s &&
                    samtools index %(outfile)s '''

    P.run(statement)


@follows(mkdir("UMI_tools.dir"))
@transform(samtools_sort_index,
           regex("featureCounts.dir/(\S+)/(\S+)_assigned_sorted.bam"),
           r"UMI_tools.dir/\1_Expression_table.tsv")
def UMI_tools(infile, outfile):
    '''
    Count using UMI tools
    '''    

    
    statement = ''' umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XT --per-cell 
                    -I %(infile)s -S %(outfile)s --wide-format-cell-counts '''

    P.run(statement)


## Now need one final R script to create matrix of results
@follows(mkdir("viral_filtered_features.dir"))
@transform(UMI_tools,
           regex("\S+/(\S+)_Expression_table.tsv"),
           r"viral_filtered_features.dir/\1/barcodes.tsv")
def matrix_report(infile, outfile):
    ''' 
    Sparse matrix of features and folder containing 
    summarised results, e.g. QC filters and counts etc.
    '''
    R_ROOT = os.path.join(os.path.dirname(__file__), "R")
    sample_name = outfile.split('/')[1]
    outdir = os.path.dirname(outfile)

    statement = ''' Rscript %(R_ROOT)s/matrix_maker.R 
                    -s %(sample_name)s -o %(outdir)s -e %(infile)s '''

    P.run(statement)

    

@follows(STAR_map, samtools_chromosome_count, samtools_index, viral_filter, viral_BAM,
human_filter, human_BAM, viral_QC, merge_viruses, samtools_merge, feature_counts, 
samtools_sort_index, UMI_tools, matrix_report)
def full():
    '''
    Runs everything
    '''

    pass


    
def main(argv=None):
    if argv is None:
        argv = sys.argv
    P.main(argv)

if __name__ == "__main__":
    sys.exit(P.main(sys.argv))


