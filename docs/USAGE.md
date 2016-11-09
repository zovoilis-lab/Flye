ABruijn manual
==============

Quick usage
-----------

    usage: abruijn [-h] [--debug] [--resume] [-t THREADS] [-i NUM_ITERS]
                   [-p {pacbio,nano,pacbio_hi_err}] [-k KMER_SIZE]
                   [-o MIN_OVERLAP] [-m MIN_KMER_COUNT] [-x MAX_KMER_COUNT]
                   [--version]
                   reads out_dir coverage
    
    ABruijn: assembly of long and error-prone reads
    
    positional arguments:
      reads                 path to reads file (FASTA format)
      out_dir               output directory
      coverage              estimated assembly coverage (integer)
    
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

### E. coli P6-C4 PacBio data

The original dataset is available at the PacBio website (https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly).
We coverted the raw 'bas.h5' file to the FASTA format for the convenience.

    wget https://github.com/fenderglass/datasets/raw/master/pacbio/E.coli_PacBio_40x.fasta
	abruijn E.coli_PacBio_40x.fasta out_pacbio 40 --platform pacbio --threads 4

with '40' being the dataset's coverage, the threads argument being optional 
(you may adjust it for your environment), and 'out_pacbio' being the directory
where the assembly results will be placed.

### E. coli Oxford Nanopore Technologies data

The dataset was originally released by the Loman lab (http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/).

    wget https://github.com/fenderglass/datasets/raw/master/ont/Loman_E.coli_MAP006-1_2D_50x.fasta
	abruijn Loman_E.coli_MAP006-1_2D_50x.fasta out_nano 50 --platform nano --threads 4


Supported Input Data
--------------------

ABruijn was designed for assembly of long reads from both PacBio and 
Oxford Nanopore Technologies (ONT). For simplicity, input reads should 
be in FASTA format - you will need to convert raw PacBio / ONT data to
FASTA format using the corresponding official tools.

### PacBio data

ABruijn was tested on the newest P6-C4 chemistry data with error rates 11-15%.
Typically, a bacterial WGS project with 20x-30x+ coverage can be assembled 
into a single, structurally concordant contig for each chromosome. However, to get the best 
nucleotide-level quaity, you might need deeper coverage. For a 55x E. coli dataset,
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
might consider reducing the coverage by filtering out shorter reads - this
may reduce the running time and memory footprint without affecting the quality. 
However, for some complicated genomes (such as those enriched with mosaic tandem repeats),
incorporating all available reads may be preferred for obtaining a more accurate structural assembly.

Assembly of small eukaryotic genomes (such as yeast or Drosophila) is also supported.
However, the processing of larger genomes (more than 500 Mb) currently might
be problematic due to the memory requirements (see below).


### Oxford Nanopore Technologies data

We performed our benchmarks with ONT 2D pass reads with error rates 13-19%.
Due to the increased error rate, you might need deeper coverage 
to get a complete chromosome assembly (60x as in the E. coli example above). For lower coverage datasets
(<30x), you might need to adjust some parameters (as described below) to get complete chromosomes.
Due to the biased error pattern, per-nucleotide accuracy is usually lower for ONT data than with 
PacBio data, especially in homopolymer regions.

Data Preparation
----------------

ABruijn works directly with base-called raw reads and does not require any 
prior error correction. The algorithm will work on corrected reads as well, 
but the running time might increase.

ABruijn automatically detects chimeric reads or reads with low quality ends, 
so you do not need to trim them before the assembly. However, it is always
worth checking for possible contamination in the reads, since it may affect the 
automatic selection of estimated parameters for solid kmers and genome size / coverage.


Running time and memory requirements
------------------------------------

Typically, assembly of a bacteria with 50x coverage takes less than an hour on a modern desktop,
while yeast assembly takes about 5 hours. A eukaryotic genome of size 200 Mbp
can be assembled within a day on a computational server with 64 CPUs.

The amount of memory required scales linearly with the genome size and read coverage.
A rough formula for estimating memory requirements is:
1 Gb of RAM = 1 Mb of genome x coverage / 50
For example, an assembly of 500 Mb genome with 50x coverage would require
approximately 500 Gb of memory. Below are running times and memory footprints 
for different datasets.

	Genome         Size     Coverage   Wall_clock   CPU_time   RAM
	E. Coli        4. 6Mb   50         44m          2h40m      2 Gb
	X. Oryzae      5.1 Mb   140        2h55m        10h        15 Gb
	S. Cerevisiae  12.2 Mb  120        4h50m        19h20m     28 Gb
	B. Neritina*   200 Mb   30         48h5m        1400h      278 Gb

* B. Neretina assembly also included symbiotic bacteria genomes.

Parameter descriptions
----------------------

### Estimated assembly coverage (required)

ABruijn requires an estimate of genome coverage as input for 
the selection of solid kmers. This can be a very rough estimate
(e.g. within 0.5x - 2x of the real coverage).

### Kmer size

This parameter controls the size of the kmers used to construct the
ABruijn graph. The default kmer size (15) is suitable for most genomes
under several hundreds of megabytes in size. You might want to increase it
a bit (17 or 19) for larger genomes, which will also require more memory
for processing.


### Minimum / maximum kmer frequency

Defines which kmers to select for ABruijn graph construction.
Kmers with low frequency are likely to be erroneous, while extremely
frequent kmers might significantly slow down the computation.
This parameter depends on the size of the genome, coverage
and sequencing platform. We recommend that it be chosen automatically 
by the assembler.

### Minimum overlap length

This sets a minimum overlap length for two reads to be considered truly overlapping.
Since the algithm is based on approximate overlaps (without alignment), we require
relatively long overlaps (5000 by default), which is suitable for most datasets
with coverage 30x+. However, you may decrease this parameter for better contiguity
on low-coverage datasets. You may also increase it to account for the presence of
a large number of long repetitive elements.

### Sequencing platform	

Indicate whether the data was generated from PacBio or ONT sequencing
so ABruijn can account for their different error patterns.

### Number of polishing iterations

ABruijn first constructs a draft assembly of the genome, which is a 
concatenation of a collection of raw read segments. Then, the draft assembly
is polished into a high quality sequence. By default, ABruijn runs two polishing 
iterations: the first one corrects the lion's share of the errors, while
the second and additional iteratiosn may correct a small number of additional
errors (due to improvements on how reads may align to the corrected assembly).

For some ONT datasets you may increase the number of iterations to get a 
better consensus. You may also use only one iteration, if the running time
is more critical than fixing a small portion of errors (since polishing takes
longer for ONT datasets). If you want to use an external polisher, you may set
the number of iterations to 0 and use the draft sequence from the assembly stage.

### Resuming existing jobs

Use --resume to resume a previous run of the assembler that may have terminated
prematurely. The assembly will continue from the last previously completed step.


Algorithm description
---------------------

This is a brief description of the ABruijn algorithm. Plase refer to the manuscript
for more detailed information. The assembly pipeline is organized as follows:

* Kmer counting / erronoeus kmer pre-filtering
* Solid kmer selection (kmers with sufficient frequency, which are unlikely to be erroneous)
* Finding read overlaps based on the ABruijn graph
* Detection of chimeric sequences
* Contig assembly by read extension

The resulting contig assembly is now simply a concatenation of read parts 
and is error-prone. The next pipeline steps are aimed at polishing this
draft assembly into a high-quality one.

* Alignment of all reads to the draft assembly using BLASR
* Selection of solid regions
* Partition the total alignment of all reads into mini-alignments (bubbles)
* Error correction of each bubble using a maximum likelihood approach

The polishing steps are repeated for the output assembly (2 iterations by default).

