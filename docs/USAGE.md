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

You can try ABruijn assembly on these ready-to-use datasets:

### E. Coli Oxfore Nanopore data, released by the Loman lab:

    wget http://nanopore.s3.climb.ac.uk/MAP006-1_2D_pass.fasta
	abruijn.py MAP006-PCR-1_2D_pass.fasta out_abruijn 60 --platform nano --threads 4

Where 60 is the dataset's coverage. Threads argument is optional, 
you may adjust it for your environment. The assembly results will
be placed in 'out_abruijn' directory.


Supported Input Data
--------------------

ABruijn was designed for assembly of long reads from both PacBio and 
Oxford Nanopore Technologies (ONT). For simplicity, input reads should 
be in FASTA format, you will need to convert raw PacBio / ONT data to FASTA
using the corresponding official tools.

### PacBio data

ABruijn was tested on the newest P6-C4 chemistry data with error rates 11-15%.
Typically, a bacterial WGS project with 20x-30x+ coverage could be assembled 
into a single, structurally concodrant contig per chromosome. However, to get the best 
nucleotide-level quaity, you might need deeper coverage. For 55x E. Coli assembly
ABruijn makes roughly 20 errors (single nucleotide insertions/deletions). 
Below are empirical error estimates for E. Coli assemblies with lower coverage:

    cov.   errors
    50x    33
	45x    45
	40x    84
	35x    153
	30x    291
	25x    687

If the coverage of the bacterial dataset is significantly higher than 100x, you
might consider to reduce the coverage by filtering out shorter reads - this
might reduce the running time and memory footprint without affecting the quality. 
However, for some compliacated genomes (such as enriched with mosaic tandem repeats)
you might need all available reads for an accurate structural assembly.

Assembly of small eukaryotic genomes (such as yeast or drosophila) is also supported.
However, the processing of larger genomes (more than 500 Mb) currently might
be problemmatic due to the memory requirements (see below).


### Oxford Nanopore Technologies data

We performed our benchmarks with ONT 2D pass reads with error rates 13-19%.
Due to the increased error rate, you might need deeper coverage 
to get a complete chromosome assembly (60x as in the E. Coli example above). For lower coverage datasets
(<30x) you might need to adjust some parameters (as described below) to get complete chromosomes.
Due to the biased error pattern, per-nucleotide accuracy is usually lower than with 
PacBio assembly, especially in homopolimer regions.

Data Preparation
----------------

ABruijn works directly with base-called raw reads and does not require any 
prior error correction. The algorithm will work on the corrected reads as well, 
but the running time might increase.

ABruijn automatically detects chimeric reads or reads with low quality ends, 
so you do not need to trim them before the assembly. However, it is always
worth chechking for a possible contamination in reads, since it may affect the 
selection of solid kmers and genome size / coverage estimates.


Running time and memory requirements
------------------------------------

Typically, assembly of a bacteria with 50x coverage takes less than an hour on a modern desktop,
while yeast assembly takes about 5 hours. A eukariotyc genome of size 200 Mbp
could be assembled within a day on a computational server with 64 CPUs.

Memory requirement scales linearly with the genome size and reads coverage.
A rough estimate is 1 Gb of RAM = 1 Mb of genome x coverage / 50.
For example, an assembly of 500 Mb genome with 50x coverage would require
approximately 500 Gb of memory. Below are running times and memory footprints 
for different datasets.

    Genome         Size     Coverage   Wall_clock   CPU_time   RAM
    E. Coli        4. 6Mb   50         44m          2h40m      2 Gb
    X. Oryzae      5.1 Mb   140        2h55m        10h        15 Gb
	S. Cerevisiae  12.2 Mb  120        4h50m        19h20m     28 Gb
	B. Neritina*   200 Mb   30         48h5m        1400h      278 Gb

* B. Neretina assembly also included symbiotic bacteria genomes.

Parameters description
----------------------

### Estimated assembly coverage (required)

ABruijn requires an estimate of genome coverage as input for 
the selection of solid kmers. This estimate could be rough
(e.g. within 0.5x - 2x of a real coverage).

### Kmer size

This parameter controls the size of kmers used to construct the ABruijn graph.
The default kmer size (15) is sutable for most of the genomes under
several hunderds of megabytes in size. You might want to increase it
a bit (16-17) for larger genomes, which will also require more memory
for the processing.


### Minimum / maximum kmer frequency

Defines which kmers to select for ABruijn graph construction.
Kmeres with low frequency are likely to be erroneous, while too
frequent kmers might significantly slow down the coputation.
This parameter depends on the size of the genome, coverage
and sequencing platform. We recommend it to be chosen automatically 
be the assembler.

### Minimum overlap

This sets a minimum overlap size between two reads which is considered significant.
Since the algithm is based on approximate overlaps (without alignment), we require
relatively long overlaps (5000 by default), which is suitable for the most datasets
with coverage 30x+. However, you may decrease this parameter for better contguity
of low-coverage datasets. You may also icrease it to account for a large number of
repetitive elements of a particular size.

### Sequencing platform	

A choice of PacBio / ONT to account for different error patternts.

### Number of polishing iterations

ABruijn first constructs a draft assembly of the genome, which is a 
concatenation of different parts of the raw reads. Then, the draft assembly
is polished into a high quality sequence. Be default, ABruijn runs two polishing 
iterations: the first on corrects the lion's share fo the errors, while
the second iterations may also correct a few (since the alignment of reads
on the assembled sequence may change).

For some ONT datasets you may increase the number of iterations to get a 
better consensus. You may also use only one iteration, if the running time
is more critical than fixing a small portion of errors. 
If you want to use an external polisher, you may set
the number of iterations to 0 and use the draft sequence from the assembly stage.

### Resuming the existing jobs

Use --resume to resume the existing assembly from the last completed step.


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
and is error-prone. The next pipeline steps are aimed to polish this
draft assembly to a high quality.

* Alignment of all read on the draft assembly using BLASR
* Selection of solid regions
* Reads alignment is partitioned into mini-alignments (bubbles)
* Error correction of each bubble using the maximum likelihood approach

The polishing part is repeated multiple times (2 by default).

