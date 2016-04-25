ABruijn assembler
==================

Version: 0.1b

The software is currently under active development, so not all features
are fully supported yet. Stay tuned for more stable releases!


Description
-----------
ABruijn is a long read assembler based on A-Bruijn graphs. 
It is described at http://biorxiv.org/content/early/2016/04/13/048413.


Install
-------
Just type 'make'. 
ABruijn also requires the Blasr aligner [https://github.com/PacificBiosciences/blasr] 
to be installed.


Usage
-----
    
    usage: abruijn.py [-h] [-o output_dir] [--debug] [-t THREADS] [-k KMER_SIZE]
                      [-m MIN_COV] [--version]
                      reads
    
    ABruijn: assembly of long anderror-prone reads
    
    positional arguments:
      reads                 path to a file with reads in FASTA format
    
    optional arguments:
      -h, --help            show this help message and exit
      -o output_dir, --outdir output_dir
                            output directory (default: abruijn-out)
      --debug               enable debug output (default: False)
      -t THREADS, --threads THREADS
                            number of parallel threads (default: 1)
      -k KMER_SIZE, --kmer-size KMER_SIZE
                            kmer size (default: 15)
      -m MIN_COV, --min-cov MIN_COV
                            minimum kmer coverage (default: 8)
      --version             show program's version number and exit


Citation
--------
Yu Lin, Jeffrey Yuan, Mikhail Kolmogorov, Max W Shen, Pavel Pevzner, 
"Assembly of Long Error-Prone Reads Using de Bruijn Graphs" (submitted)


Third-party
-----------
ABruijn package includes some third-patry software:

* libbf [http://github.com/mavam/libbf]


License
-------
ABruijn is dictributed under BSD license. See *LICENSE* file for details.


Contacts
--------
Please report any problems directly to the github issue tracker.
Also, you can send feedback to fenderglass@gmail.com
