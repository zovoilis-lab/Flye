ABruijn assembler
==================

Version: 0.4b

The software is currently under active development.


Description
-----------
ABruijn is a long read assembler based on A-Bruijn graphs. 
It is described at http://biorxiv.org/content/early/2016/04/13/048413.


Install
-------
Just type 'make'. 
ABruijn also requires the Blasr aligner [https://github.com/PacificBiosciences/blasr] 
to be installed and available through PATH environment variable.


Usage
-----
    
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



Citation
--------
Yu Lin, Jeffrey Yuan, Mikhail Kolmogorov, Max W Shen, Pavel Pevzner, 
"Assembly of Long Error-Prone Reads Using de Bruijn Graphs" (submitted)


Third-party
-----------
ABruijn package includes some third-patry software:

* libcuckoo [http://github.com/efficient/libcuckoo]


License
-------
ABruijn is dictributed under BSD license. See *LICENSE* file for details.


Contacts
--------
Please report any problems directly to the github issue tracker.
Also, you can send feedback to fenderglass@gmail.com
