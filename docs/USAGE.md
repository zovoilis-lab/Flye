    usage: abruijn.py [-h] [--debug] [--resume] [-t THREADS] [-i NUM_ITERS]
                      [-k KMER_SIZE] [-m MIN_COV] [-x MAX_COV] [--version]
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
      -m MIN_COV, --min-cov MIN_COV
                            minimum kmer coverage (default: auto)
      -x MAX_COV, --max-cov MAX_COV
                            maximum kmer coverage (default: auto)
      --version             show program's version number and exit



