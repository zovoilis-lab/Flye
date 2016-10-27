ABruijn assembler
==================

Version: 0.4b

ABruijn is a de novo assembler for PacBio and Oxford Nanopore Technologies reads.
The algorithm is using A-Bruijn graph to find the overlaps between reads
and does not require them to be error-corrected.  First, the algorithm produces
a draft assembly by concatenating different parts of raw reads.
This coarse sequence is then polished to a high quality assembly.

ABruijn works for both bacterial and eukaryotic genomes. Typically, assembly
of a bacteria with 50x coverage takes less than an hour on a modern desktop,
while yeast assembly takes about 5 hours. A eukariotyc genome of size 200 Mbp
could be assembled within a day on a computational server.


Install
-------
See *docs/INSTALL.md* file.


Usage
-----
See *docs/USAGE.md* file.


Publications
------------
Yu Lin, Jeffrey Yuan, Mikhail Kolmogorov, Max W Shen, Pavel Pevzner, 
"Assembly of Long Error-Prone Reads Using de Bruijn Graphs"
(http://biorxiv.org/content/early/2016/04/13/048413)


Third-party
-----------
ABruijn package includes some third-patry software:

* libcuckoo [http://github.com/efficient/libcuckoo]


License
-------
ABruijn is dictributed under BSD license. See *LICENSE* file for details.


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
