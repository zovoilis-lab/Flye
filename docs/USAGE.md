ABruijn manual
==============

Quick usage
-----------

    usage: abruijn.py [-h] [--debug] [--resume] [-t THREADS] [-i NUM_ITERS]
                      [-p {pacbio,nano,pacbio_hi_err}] [-k KMER_SIZE]
                      [-o MIN_OVERLAP] [-m MIN_KMER_COUNT] [-x MAX_KMER_COUNT]
                      [--version]
                      reads out_dir coverage
    
    ABruijn: assembly of long and error-prone reads
    
    positional arguments:
      reads                 path to a file with reads in FASTA format
      out_dir               output directory
      coverage              estimated assembly coverage
    
    optional arguments:
      -h, --help            show this help message and exit
      --debug               enable debug output
      --resume              try to resume previous assembly
      -t THREADS, --threads THREADS
                            number of parallel threads (default: 1)
      -i NUM_ITERS, --iterations NUM_ITERS
                            number of polishing iterations (default: 2)
      -p {pacbio,nano,pacbio_hi_err}, --platform {pacbio,nano,pacbio_hi_err}
                            sequencing platform (default: pacbio)
      -k KMER_SIZE, --kmer-size KMER_SIZE
                            kmer size (default: 15)
      -o MIN_OVERLAP, --min-overlap MIN_OVERLAP
                            minimum overlap between reads (default: 5000)
      -m MIN_KMER_COUNT, --min-coverage MIN_KMER_COUNT
                            minimum kmer coverage (default: auto)
      -x MAX_KMER_COUNT, --max-coverage MAX_KMER_COUNT
                            maximum kmer coverage (default: auto)
      --version             show program's version number and exit



Examples
--------

You can try ABruijn assembly on these ready-to-use datasets right away:

### E. Coli Oxfore Nanopore data, released by the Loman lab:

    wget http://nanopore.s3.climb.ac.uk/MAP006-1_2D_pass.fasta
	abruijn.py MAP006-PCR-1_2D_pass.fasta out_abruijn 60 --platform nano --threads 4

Where 60 is the dataset's coverage. Threads argument is optional, 
you may adjust it for your environment. The assembly results will
be placed in 'out_abruijn' directory


Supported Input Data
--------------------

ABruijn was designed for assembly of long reads from both PacBio and 
Oxford Nanopore technologies. Input reads should be in FASTA format,
so you will need to convert raw PacBio / Oxford Nanopore files to FASTA format 
using the corresponding official tools prior running ABruijn.

### PacBio data

ABruijn was tested on the newest P6-C4 chemistry data with error rates 11-15%.
Typically, 20x-30x coverage is enough for bacterial assembly to get a single contig
for each chromosome without structural misassemblies. However, to get the best 
nucleotide-level quaity, you might need deeper coverage. For 55x E. Coli assembly
ABruijn makes about 20 errors per genome (single nucleotide insertions/deletions). 
Below are more detailed error estimates for E. Coli assemblies with lower coverage:

    cov.   errors
    50x    33
	45x    45
	40x    84
	35x    153
	30x    291
	25x    687


If the coverage of your bacterial dataset is significantly higher than 100x, you
might consider to reduce the coverage by filtering out shorter reads - this
might reduce the running time without affecting the quality. However, for some
compliacated genomes (such as enriched with complicated tandem repeats) you might
need all reads for an accurate structural assembly.


### Oxford Nanopore data

We performed our benchmarks with Oxford Nanopore 2D pass reads with error rates 13-19%.
As the reads are usually shorter and less accurate, you might need a deeper coverage 
to get a complete chromosome assembly (60x as in the example above). For lower coverage datasets
(30x) you might need to adjust some parameters (as described below) to get full chromosomes.
Per-nucleotide accuracy is usually lower than with PacBio technology, especially in 
homopolimer regions.

Data Preparation
----------------

ABruijn works directly with raw reads (after base calling) and does not require any 
prior error correction. The software will work on the corrected reads as well, 
but the running time might increase. 

ABruijn automatically detects chimeric reads or reads with low quality ends, 
so you do not need to manually do it before the assembly. However, it is always
worth chechking reads for contamination, since it may affect the selection of solid
kmers.


Algorithm description
---------------------

This is a brief description of ABruijn algorithm. Plase refer to the manuscript
for more detailed information. The assembly pipeline is organized as follows:

* Kmer counting / erronoeus kmer pre-filtering
* Solid kmer selection (kmers with sufficient frequency, that are unlikely to be erroneous)
* Finding read overlaps based on ABruijn graph
* Detection of chimeric sequences
* Contig assembly by read extension

The resulting contig assembly is now simply a concatenation of read parts, 
so it is error-prone. The next pipeline steps are aimed to polish this
draft assembly to a high quality.

* Alignment of all read on draft assembly using BLASR
* Selection of solid regions
* Read alignment is partitioned into mini-alignments (bubbles)
* Error correction of each bubble using the maximum likelihood approach

The polishing part is repeated multiple times (2 by default).


Running time and memory requirements
------------------------------------

Typically, assembly of a bacteria with 55x coverage takes less than an hour on a modern desktop,
while yeast assembly takes about 5 hours. A eukariotyc genome of size 200 Mbp
could be assembled within a day on a computational server with 64 CPUs.


Parameters description
----------------------



Resuming the existing jobs
--------------------------

Type --resume to resume the existing assembly from the last completed step.
