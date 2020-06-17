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

@transform(STAR_map,
           regex("STAR.dir/(\S+)/(\S+)_Aligned.sortedByCoord.out.bam"),
           r"STAR.dir/\1/\1_<whatever_file_suffix_is>.bam")
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
           r"\1/\1_Count_chromsomes.txt")
def samtools_chromosome_count(infile, outfile):
    ''' 
    Compute the number of mapped reads for each chromosome/virus
    '''

    statement = '''samtools idxstats %(infile)s > %(outfile)s'''

    job.memory = "20G"

    P.run(statement)


@follows(STAR_map, samtools_index, samtools_chromosome_count)
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


