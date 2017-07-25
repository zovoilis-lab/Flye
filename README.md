ABruijn assembler
==================

Version: 2.0b

ABruijn is a de novo assembler for PacBio and Oxford Nanopore Technologies reads.
The algorithm uses an A-Bruijn graph to find the overlaps between reads
and does not require them to be error-corrected.  First, the algorithm produces
a draft assembly by concatenating different parts of raw reads.
This coarse sequence is then polished into a high quality assembly.

ABruijn works for both bacterial and eukaryotic genomes. Typically, assembly
of a bacteria with 50x coverage takes less than half an hour on a modern desktop,
while 100x yeast assembly takes about two hours. A eukaryotic genome of size 200 Mbp
can be assembled within a day on a computational server.


Install
-------
See the *docs/INSTALL.md* file.


Usage
-----
See the *docs/USAGE.md* file.


Publications
------------
Yu Lin, Jeffrey Yuan, Mikhail Kolmogorov, Max W Shen, Pavel Pevzner, 
"Assembly of Long Error-Prone Reads Using de Bruijn Graphs"
(http://biorxiv.org/content/early/2016/04/13/048413)


Third-party
-----------
ABruijn package includes some third-party software:

* libcuckoo [http://github.com/efficient/libcuckoo]


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
