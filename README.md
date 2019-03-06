Flye assembler
==============

[![BioConda Install](https://img.shields.io/conda/dn/bioconda/flye.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/flye)

### Version: 2.4.1

Flye is a de novo assembler for single molecule sequencing reads,
such as those produced by PacBio and Oxford Nanopore Technologies.
We support a wide range of datasets, from small bacterial projects
to large mammalian-scale assemblies. The package represets a complete
pipeline: it takes raw PB / ONT reads as input and outputs polished contigs.

Flye is using repeat graphs to identify genomic repeats
and optimally resolve them, which results into accurate assemblies.
In difference to the current state-of-the-art approaches, Flye does
not error-correct input reads, thus the assembly is faster.
Flye also includes a special mode for metagenome assembly.

Latest updates
--------------

### Flye 2.4.1 release (05 Mar 2019)
* Speed and stability improvements for large datasets
* New option `--polish-target` to run Flye polisher on the target sequence

### Flye 2.4 release (14 Jan 2019)
* Metagenome assembly support fully integrated (`--meta` option)
* New Trestle module for resolving simple unbridged repeats
* New `--plasmids` option that recovers short unassmbled plasmids


Manuals
-------

- [Installation instructions](docs/INSTALL.md)
- [Usage](docs/USAGE.md)


Repeat graph
------------

The Flye algorithms are using repeat graph as a core data structure. 
In difference to de Bruijn graphs which require exact k-mer matches,
repeat graphs are built using apporximate sequence matches, thus
can tollerate higher noise of SMS reads.

The edges of repeat graph represent genomic sequence, and nodes define
the junctions. All edges are classified into uniqe and repetitive.
The genome traverses the graph in an unknown way, so as each unique
edge appears exactly once in this traversal. Repeat graphs are useful
for repeat analysis and resolution - which are one of the key 
genome assembly challenges.


<p align="center">
  <img src="docs/graph_example.png" alt="Graph example"/>
</p>

Above is an example of a repeat graph of a bacterial assembly.
Each edge is labeled with its id, length and coverage. Repetitive edges are shown
in color, and unique edges are black. Note that each edge is represented in 
two copies: forward and reverse complement (marked with +/- signs), 
therefore the entire genome is represented in two copies as well. 

In this example, there are two unresolved repeats: (i) a red repeat of 
multiplicity two and length 35k and (ii) a green repeat cluster of multiplicity
three and length 34k - 36k. As the repeats remained unresolved, there are no reads
in the dataset that cover those repeats in full. Five unique edges 
will correspond to five contigs in the final assembly.

Repeat graphs produced by Flye could be visualized using
[AGB](https://github.com/almiheenko/AGB) or [Bandage](https://github.com/rrwick/Bandage).


Third-party
-----------

Flye package includes some third-party software:

* [libcuckoo](http://github.com/efficient/libcuckoo)
* [intervaltree](https://github.com/ekg/intervaltree)
* [lemon](http://lemon.cs.elte.hu/trac/lemon)
* [minimap2](https://github.com/lh3/minimap2)


License
-------

Flye is distributed under a BSD license. See the [LICENSE file](LICENSE) for details.


Credits
-------

Flye is developed in [Pavel Pevzner's lab at UCSD](http://cseweb.ucsd.edu/~ppevzner/)

Code contributions:

* Repeat graph and current package maintaining: Mikhail Kolmogorov
* Trestle module and original polisher code: Jeffrey Yuan
* Original contig extention code: Yu Lin
* Short plasmids recovrey module: Evgeny Polevikov


Publications
------------
Mikhail Kolmogorov, Jeffrey Yuan, Yu Lin and Pavel Pevzner, 
"Assembly of Long Error-Prone Reads Using Repeat Graphs", bioRxiv, 2018
[doi:10.1101/247148](https://doi.org/10.1101/247148)

Yu Lin, Jeffrey Yuan, Mikhail Kolmogorov, Max W Shen, Mark Chaisson and Pavel Pevzner, 
"Assembly of Long Error-Prone Reads Using de Bruijn Graphs", PNAS, 2016
[doi:10.1073/pnas.1604560113](https://www.doi.org/10.1073/pnas.1604560113)


Contacts
--------
Please report any problems to the [issue tracker](https://github.com/fenderglass/Flye).
If possible, please include "flye.log" file from the output directory
for faster feedback. Alternatively, you can write directly to fenderglass@gmail.com.
