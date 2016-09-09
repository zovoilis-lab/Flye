ABruijn manual
==============

Quick usage
-----------

    usage: abruijn.py [-h] [--debug] [--resume] [-t THREADS] [-i NUM_ITERS]
                      [-k KMER_SIZE] [-o MIN_OVERLAP] [--version]
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
      -k KMER_SIZE, --kmer-size KMER_SIZE
                            kmer size (default: 15)
      -o MIN_OVERLAP, --min-overlap MIN_OVERLAP
                            minimum overlap between reads (default: 5000)
      --version             show program version number and exit


Example
-------

A couple of examples of assembly of publicly available datasets.

Data Requirements
-----------------

ABruijn was designed for assembly of reads from both Pacbio and 
Oxford Nanopores technologies.

For both technologies: coverage, minimum quality etc. Give an overview of errors
depending on coverage. Both bacterial / eukaryote

Data Preparation
----------------

Filtering short reads, QC, trimming etc...


Algorithm description
---------------------

Briefly describe different stages of the algorithm


Running time and memory requirements
------------------------------------

Per algorithm stages

Parameters setting
------------------

Coverage (and its importance), kmer size


Resuming the existing jobs
--------------------------
