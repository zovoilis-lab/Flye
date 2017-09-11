ABruijn assembler
==================

Version: 2.1b

ABruijn is a de novo assembler for long and noisy reads, such as
those produced by PacBio and Oxford Nanopore Technologies.
The algorithm uses an A-Bruijn graph to find the overlaps between reads
and does not require them to be error-corrected. The package includes a 
polisher module, which produces assembly of high nucleotide-level quality.

Since the version 2.0, ABruijn performs additional repeat analysis
step, which improves the structural accuracy of the resulting sequence. 
The algorithm also produces a graph representation of the final assembly.

ABruijn has moderate memory requirements and is designed to run on a single node.
Typically, assembly of a bacteria with 50x coverage takes less than half an hour 
on a modern desktop and yeast assembly takes about 2 hours. A metagenome of total size 200 Mbp
can be assembled within a day on using 70Gb of memory and 20 CPUs.

Install
-------
See the *docs/INSTALL.md* file.


Usage
-----
See the *docs/USAGE.md* file.


Publications
------------
Yu Lin, Jeffrey Yuan, Mikhail Kolmogorov, Max W Shen, Mark Chaisson and Pavel Pevzner, 
"Assembly of Long Error-Prone Reads Using de Bruijn Graphs", PNAS 2016


Third-party
-----------
ABruijn package includes some third-party software:

* libcuckoo [http://github.com/efficient/libcuckoo]
* BLASR [https://github.com/PacificBiosciences/blasr]


License
-------
ABruijn is distributed under a BSD license. See the *LICENSE* file for details.


Credits
-------

ABruijn was developed in [Pavel Pevzner's lab at UCSD](http://cseweb.ucsd.edu/~ppevzner/)

Code contributions:

* Original assembler code: Yu Lin
* Original polisher code: Jeffrey Yuan
* Current package: Mikhail Kolmogorov


Contacts
--------
Please report any problems directly to the github issue tracker.
Also, you can send feedback to fenderglass@gmail.com
