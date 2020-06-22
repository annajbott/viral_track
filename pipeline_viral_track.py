""" =====================
Pipeline viral track
=========================

Overview 
=========

This workflow was designed as a CGAT pipeline representation of R viral_track code by Bost et al.
Viral track can be found:
- https://www.sciencedirect.com/science/article/pii/S0092867420305687
- https://github.com/PierreBSC/Viral-Track

Adapted by Lauren Overend (Lauren.overend@oriel.ox.ac.uk)

Any CGAT pipeline specific questions direct to anna.james-bott@st-hildas.ox.ac.uk


Requires:

* fastq files 
* Reference genome 
* ... 

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

    
# load options from the config file
PARAMS = P.get_parameters(
    ["%s/pipeline.yml" % os.path.splitext(__file__)[0],
     "../pipeline.yml",
     "pipeline.yml"])

R_ROOT = os.path.join(os.path.dirname(__file__), "R")

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


@follows(mkdir("STAR.dir"))
@transform(SEQUENCEFILES,
           regex("(\S+).fa"),
           r"STAR.dir/\1/\1_Aligned.sortedByCoord.out.bam")
def STAR_map(infile, outfile):
    '''
    Run STAR mapping with parameters defined in yml
    '''
    
    nthreads = PARAMS['star_nThreads']
    nBAMsortingthreads = PARAMS['star_nThreadsort']
    nBins = PARAMS['star_bins']
    index_genome = PARAMS['genome_dir']
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
## samtools merge, samtools sort


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
           r"\1/\2_Count_chromsomes.txt")
def samtools_chromosome_count(infile, outfile):
    ''' 
    Compute the number of mapped reads for each chromosome/virus
    '''

    statement = '''samtools idxstats %(infile)s > %(outfile)s'''

    job.memory = "20G"

    P.run(statement)


@split(samtoools_chromosome_count, 
           regex("(\S+)/(\S+)_Count_chrosomes.txt"),
           r"\1/Viral_BAM_files/virus_file_names/*.txt")
def viral_filter(infiles, outfile):
    '''
    Uses R script to filter out human chromosomes and viruses
    under threshold number of reads, then create empty txt file for each virus 
    '''

    chromosome_count = infile
    BAM_folder = os.path.dirname(os.path.dirname(outfile))
    outdir = os.path.dirname(outfile)
    minreads <- PARAMS['star_minreads']

    os.mkdir(BAM_folder)
    os.mkdir(outdir)

    statement = '''Rscript %(R_ROOT)s/viral_bam_filter.R -c %(chromosome_count)s 
                -m %(minreads)s -o %(outdir)s '''
    
    P.run(statement)

@transform(viral_filter,
           regex("STAR.dir/(\S+)/Viral_BAM_files/virus_file_names/(\S+).txt)
           add_inputs(STAR_map),
           r"STAR.dir/\1/Viral_BAM_files/\2.bam")
def viral_BAM(infiles, outfile):
    ''' 
    Takes virus names from empty text files and aligned bam 
    and makes bam file for each virus
    '''

    virus_name_file, aligned_bam = infiles
    virus =  os.path.basename(virus_name_file).replace(".txt", "")


    statement = """ samtools view -b %(aligned_bam)s '%(virus)s' > %(outfile)s """

    P.run(statement)


@split(samtoools_chromosome_count, 
           regex("(\S+)/(\S+)_Count_chrosomes.txt"),
           r"\1/Human_BAM_files/human_file_names/*.txt")
def human_filter(infiles, outfile):
    '''
    Uses R script to filter out anything not human and any chromosomes
    under threshold number of reads, then create empty txt file for each name 
    '''

    chromosome_count = infile
    BAM_folder = os.path.dirname(os.path.dirname(outfile))
    outdir = os.path.dirname(outfile)
    minreads <- PARAMS['star_minreads']

    os.mkdir(BAM_folder)
    os.mkdir(outdir)

    statement = '''Rscript %(R_ROOT)s/human_bam_filter.R -c %(chromosome_count)s 
                -m %(minreads)s -o %(outdir)s '''
    
    P.run(statement)

@transform(human_filter,
           regex("STAR.dir/(\S+)/Human_BAM_files/human_file_names/(\S+).txt)
           add_inputs(STAR_map),
           r"STAR.dir/\1/Viral_BAM_files/\2.bam")
def human_BAM(infiles, outfile):
    ''' 
    Takes human chromosome names from empty text files and aligned bam 
    and makes bam file for each chromsome
    '''

    human_name_file, aligned_bam = infiles
    human_chrom =  os.path.basename(human_name_file).replace(".txt", "")


    statement = """ samtools view -b %(aligned_bam)s '%(human_chrom)s' > %(outfile)s """

    P.run(statement)


# Need to finish R script ...
@collate(viral_bam,
        regex("STAR.dir/(\S+)/Viral_BAM_files/\S+.bam"),
        r"QC.dir/\1/QC_filtered.txt") # Fill in name
def viral_QC(infile, outfile):
    ''' 
    Performs quality control on viral BAM file
    '''

    viral_bam_directory = os.path.dirname(infile[0])
    
    
    statement = '''Rscript %(R_ROOT)s/QC.R --viraldir %(viral_bam_directory)s -o %(outfile)s -r %(R_ROOT)s'''


    P.run(statement)
    
# Not finished ....
@collate(viral_bam,
         regex("STAR.dir/(\S+)/Viral_BAM_files/\S+.bam"),
         add_inputs(viral_QC),
         r"STAR.dir/(\S+)/merged_viral_mapping.bam")
def merge_viruses(infiles, outfile):
    '''
    Merging BAM for all QC filtered viruses
    '''


@collate(human_BAM,
         regex("STAR.dir/(\S+)/Human_BAM_files/\S+.bam"),
         r"STAR.dir/(\S+)/merged_human_mapping.bam")
def samtools_merge(infiles, outfile):
    '''
    Merge human chromsome BAM files using samtools
    '''

    path_to_BAM = os.path.dirname(infile[0])
    infile_list = ' '.join(infiles)
    

    statement = ''' samtools merge %(outfile)s -f %(infile_list)s '''

    P.run(statement)


# Need task to create pseudo GTF...
@follows(mkdir("featureCounts.dir"))
@transform(samtools_merge,
           regex("STAR.dir/(\S+)/merged_human_mapping.bam"),
           add_inputs(GTF_task)
           r"featureCounts.dir/\1/merged_human_mapping.featureCounts.bam")
def feature_counts(infile, outfile):
    '''
    Perform feature counts on the merged BAM file
    '''

    merged_bam, gtf = infiles
    sample_name = merged_bam.split("/")[1]
    counts_file = "featureCounts.dir/" + sample_name + "/" + sample_name + "_counts.txt"
    # Move bam file to featurecounts folder and rename 
    bam_out_og = merged_bam + ".featureCounts.bam"
    
    statement = ''' featureCounts -t transcript -M --primary -R BAM -g gene_id 
                    -a  %(gtf)s -o %(counts_file)s %(merged_bam)s && 
                    mv %(bam_out_og)s %(outfile)s '''


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
           r"UMI_tools.dir/\1/Expression_table.tsv")
def UMI_tools(infile, outfile):
    '''
    Count using UMI tools
    '''    

    
    statement = ''' umi_tools count --per-gene --gene-tag=XT --assigned-status-tag=XT --per-cell 
                    -I %(infile)s -S %(outfile)s --wide-format-cell-counts '''

    P.run(statement)


## Now need one final R script to create matrix of results


@follows(STAR_map, samtools_index, samtools_chromosome_count, viral_filter, virus_BAM,
human_filter, human_BAM, viral_QC)
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


