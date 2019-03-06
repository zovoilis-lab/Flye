Flye assembler
==============

[![BioConda Install](https://img.shields.io/conda/dn/bioconda/flye.svg?style=flag&label=BioConda%20install)](https://anaconda.org/bioconda/flye)

### Version: 2.4.1

Flye is a de novo assembler for single molecule sequencing reads,
such as those produced by PacBio and Oxford Nanopore Technologies.
It is designed for a wide range of datasets, from small bacterial projects
to large mammalian-scale assemblies. The package represents a complete
pipeline: it takes raw PB / ONT reads as input and outputs polished contigs.
Flye also includes a special mode for metagenome assembly.

Latest updates
--------------

### Flye 2.4.1 release (05 Mar 2019)
* Speed and stability improvements for large datasets
* New option `--polish-target` to run Flye polisher on the target sequence

### Flye 2.4 release (14 Jan 2019)
* Metagenome assembly support fully integrated (`--meta` option)
* New Trestle module for resolving simple unbridged repeats
* New `--plasmids` option that recovers short unassembled plasmids


Manuals
-------

- [Installation instructions](docs/INSTALL.md)
- [Usage](docs/USAGE.md)


Repeat graph
------------

Flye is using repeat graph as a core data structure. 
In difference to de Bruijn graphs (which require exact k-mer matches),
repeat graphs are built using approximate sequence matches, and
can tolerate higher noise of SMS reads.

The edges of repeat graph represent genomic sequence, and nodes define
the junctions. Each edges is classified into unique or repetitive.
The genome traverses the graph (in an unknown way), so as each unique
edge appears exactly once in this traversal. Repeat graphs reveal the
repeat structure of the genome, which helps to reconstruct an optimal assembly.


<p align="center">
  <img src="docs/graph_example.png" alt="Graph example"/>
</p>

Above is an example of the repeat graph of a bacterial assembly.
Each edge is labeled with its id, length and coverage. Repetitive edges are shown
in color, and unique edges are black. Note that each edge is represented in 
two copies: forward and reverse complement (marked with +/- signs), 
therefore the entire genome is represented in two copies. This is necessary
because the orientation of input reads is unknown.

In this example, there are two unresolved repeats: (i) a red repeat of 
multiplicity two and length 35k and (ii) a green repeat cluster of multiplicity
three and length 34k - 36k. As the repeats remained unresolved, there are no reads
in the dataset that cover those repeats in full. Five unique edges 
will correspond to five contigs in the final assembly.

Repeat graphs produced by Flye could be visualized using
[AGB](https://github.com/almiheenko/AGB) or [Bandage](https://github.com/rrwick/Bandage).


Flye benchmarks
---------------

| Genome                   | Data       | Asm.Size  | NG50     | CPU time  | RAM    |
|--------------------------|------------|-----------|----------|-----------|--------|
| [E.coli][ecoli]          | PB 50x     | 4.6 Mb    | 4.6 Mb   | 2 h       | 2 Gb   |
| [C.elegans][ce]          | PB 40x     | 102 Mb    | 2.9 Mb   | 100 h     | 31 Gb  |
| [A.thaliana][at]         | PB 75x     | 120 Mb    | 11.2 Mb  | 240 h     | 46 Gb  |
| [D.melanogaster][dm-ont] | ONT 30x    | 143 Mb    | 13.2 Mb  | 130 h     | 31 Gb  |     
| [D.melanogaster][dm-pb]  | PB 120x    | 142 Mb    | 17.5 Mb  | 190 h     | 75 Gb  |     
| [NA12878][na12878]       | ONT UL 35x | 2.9 Gb    | 20.8 Mb  | 5000 h    | 600 Gb |
| [HG002][hg002]           | PB CCS 30x | 2.9 Gb    | 24.6 Mb  | 900 h     | 300 Gb |
| [HMP mock][hmp]          | PB meta    | 66 Mb     | 2.7 Mb   | 60 h      | 44 Gb  |

[na12878]: https://github.com/nanopore-wgs-consortium/NA12878/blob/master/Genome.md
[ce]: https://github.com/PacificBiosciences/DevNet/wiki/C.-elegans-data-set
[at]: https://downloads.pacbcloud.com/public/SequelData/ArabidopsisDemoData/
[dm-pb]: https://github.com/PacificBiosciences/DevNet/wiki/Drosophila-sequence-and-assembly
[dm-ont]: https://www.ebi.ac.uk/ena/data/view/SRR6702603
[hg002]: https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_15kb/
[ecoli]: https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly
[hmp]: https://github.com/PacificBiosciences/DevNet/wiki/Human_Microbiome_Project_MockB_Shotgun 

The assemblies generated using Flye 2.4 could be downloaded from [Zenodo](https://zenodo.org/record/2586130)

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
* Original contig extension code: Yu Lin
* Short plasmids recovery module: Evgeny Polevikov


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
Please report any problems to the [issue tracker](https://github.com/fenderglass/Flye/issues).
If possible, please include "flye.log" file from the output directory
for faster feedback. Alternatively, you can write directly to fenderglass@gmail.com.
