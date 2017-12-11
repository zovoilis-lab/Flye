Flye assembler (successor of ABruijn)
=====================================

Version: 2.3

Flye is a de novo assembler for long and noisy reads, such as
those produced by PacBio and Oxford Nanopore Technologies.
The algorithm uses an A-Bruijn graph to find the overlaps between reads
and does not require them to be error-corrected. After the initial assembly, 
Flye performs an extra repeat classification and analysis step to improve the structural
accuracy of the resulting sequence. The package also includes a polisher
module, which produces the final assembly of high nucleotide-level quality.

The 2.x software branch has been renamed to Flye because of
the substantial algorithmic differences.

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
Flye package includes some third-party software:

* libcuckoo [http://github.com/efficient/libcuckoo]
* intervaltree [https://github.com/ekg/intervaltree]
* lemon [http://lemon.cs.elte.hu/trac/lemon]
* minimap2 [https://github.com/lh3/minimap2]


License
-------
Flye is distributed under a BSD license. See the *LICENSE* file for details.


Credits
-------

Flye was developed in [Pavel Pevzner's lab at UCSD](http://cseweb.ucsd.edu/~ppevzner/)

Code contributions:

* Original assembler code: Yu Lin
* Original polisher code: Jeffrey Yuan
* Current package: Mikhail Kolmogorov


Contacts
--------
Please report any problems directly to the github issue tracker.
Also, you can send feedback to fenderglass@gmail.com
