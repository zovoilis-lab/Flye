Flye manual
===========

Table of Contents
-----------------

- [Quick usage](#quickusage)
- [Examples](#examples)
- [Supported Input Data](#inputdata)
- [Parameter Descriptions](#parameters)
- [Flye output](#output)
- [Repeat graph](#graph)
- [Running Time and Memory Requirements](#performance)
- [Algorithm Description](#algorithm)


## <a name="quickusage"></a> Quick usage

```
usage: flye (--pacbio-raw | --pacbio-corr | --nano-raw |
	     --nano-corr | --subassemblies) file1 [file_2 ...]
	     --genome-size SIZE --out-dir PATH
	     [--threads int] [--iterations int] [--min-overlap int]
	     [--meta] [--plasmids] [--no-trestle] [--debug]
	     [--version] [--help] [--resume]

Assembly of long and error-prone reads

optional arguments:
  -h, --help            show this help message and exit
  --pacbio-raw path [path ...]
                        PacBio raw reads
  --pacbio-corr path [path ...]
                        PacBio corrected reads
  --nano-raw path [path ...]
                        ONT raw reads
  --nano-corr path [path ...]
                        ONT corrected reads
  --subassemblies path [path ...]
                        high-quality contigs input
  -g size, --genome-size size
                        estimated genome size (for example, 5m or 2.6g)
  -o path, --out-dir path
                        Output directory
  -t int, --threads int
                        number of parallel threads [1]
  -i int, --iterations int
                        number of polishing iterations [1]
  -m int, --min-overlap int
                        minimum overlap between reads [auto]
  --asm-coverage int    reduced coverage for initial contig assembly [not set]
  --plasmids            rescue short unassmebled plasmids
  --meta                metagenome / uneven coverage mode
  --no-trestle          skip Trestle stage
  --resume              resume from the last completed stage
  --resume-from stage_name
                        resume from a custom stage
  --debug               enable debug output
  -v, --version         show program's version number and exit

```

Input reads could be in FASTA or FASTQ format, uncompressed
or compressed with gz. Currenlty, raw and corrected reads
from PacBio and ONT are supported. The expected error rates are
<30% for raw and <2% for corrected reads. Additionally,
--subassemblies option performs a consensus assembly of multiple
sets of high-quality contigs. You may specify multiple
files with reads (separated by spaces). Mixing different read
types is not yet supported. --meta option enables the mode
for metagenome/uneven coverage assembly.

You must provide an estimate of the genome size as input,
which is used for solid k-mers selection. Standard size
modificators are supported (e.g. 5m or 2.6g). In case
of metagenome assembly, expected total assembly size
should be provided.

To reduce memory consumption for large genome assemblies,
you can use a subset of the longest reads for initial contig
assembly by specifying --asm-coverage option. Typically,
40x coverage is enough to produce good draft contigs.



## <a name="examples"></a> Examples

You can try Flye assembly on these ready-to-use datasets:

### E. coli P6-C4 PacBio data

The original dataset is available at the 
[PacBio website](https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly).
We coverted the raw ```bas.h5``` file to the FASTA format for the convenience.

    wget https://zenodo.org/record/1172816/files/E.coli_PacBio_40x.fasta
	flye --pacbio-raw E.coli_PacBio_40x.fasta --out-dir out_pacbio --genome-size 5m --threads 4

with ```5m``` being the expected genome size, the threads argument being optional 
(you may adjust it for your environment), and ```out_pacbio``` being the directory
where the assembly results will be placed.

### E. coli Oxford Nanopore Technologies data

The dataset was originally released by the 
[Loman lab](http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/).

    wget https://zenodo.org/record/1172816/files/Loman_E.coli_MAP006-1_2D_50x.fasta
	flye --nano-raw Loman_E.coli_MAP006-1_2D_50x.fasta --out-dir out_nano --genome-size 5m --threads 4


## <a name="inputdata"></a> Supported Input Data

### PacBio data

Flye was tested on raw PacBio reads (P5C3 and P6C4) with error rate ~15%.
Note that Flye assumes that the input files represent PacBio subreads,
e.g. adaptors and noise are trimmed and multiple passes of the same insertion
sequence are separated. This is typically handled by PacBio instruments/toolchains,
however we saw examples of incorrect third-party raw -> fastq conversions, 
which resulted into incorrectly trimmed data. In case Flye is failing to
get reaonable assemblies, make sure that your reads are properly preprocessed.

### Oxford Nanopore data

We performed our benchmarks with raw ONT reads (R7-R9) with error rate ~15%.
Due to the biased error pattern, per-nucleotide accuracy is usually lower for 
ONT data than with PacBio data, especially in homopolymer regions.

### Error-corrected reads input

While Flye was designed for assembly of raw reads (and this is the recommended way),
it also supports error-corrected PacBio/ONT reads as input (use the ```corr``` option).
The parameters are optimized for error rates <2%. If you are getting highly 
fragmented assembly - most likely error rates in your reads are higher. In this case,
consider to assemble using the raw reads instead.

### Consensus of multiple contig sets

```--subassemblies``` input mode generates a consensus of multiple high quality contig assemblies
(such as produced by different short/long read assemblers). The expected error rate
is <1%. You might want to skip the polishing stage with ```--iterations 0``` argument
(however, it might still be helpful to correct small structural errors).


### Input data preparation

Flye works directly with base-called raw reads and does not require any 
prior error correction. Flye automatically detects chimeric reads or reads with low quality ends, 
so you do not need to curate them before the assembly. However, it is always
worth checking for possible contamination in the reads, since it may affect the 
automatic selection of estimated parameters for solid kmers and genome size / coverage.


## <a name="parameters"></a> Parameter descriptions

### Estimated genome size (required)

You must provide an estimate of the genome size as input,
which is used for solid k-mers selection. The estimate could
be rough (e.g. withing 0.5x-2x range) and does not affect
the other assembly stages. Standard size modificators are
supported (e.g. 5m or 2.6g)

### Minimum overlap length

This sets a minimum overlap length for two reads to be considered overlapping.
In the latest Flye versions, this parameter is chosen automatically
based on the read length distribution (reads N90) and does not require manual setting.
Typical value is 3k-5k (and down to 1k for datasets with shorter read length).
Intuitively, we want to set this parameter as high as possible, so the
repeat graph is less tangled. However, higher values might lead to assembly gaps.
In some *rare* cases (for example in case of biased read length distribution)
it makes sense to set this parameter manualy.

### Metagenome mode

Metagenome assembly mode, that is designed for highly non-uniform coverage and 
is sensitive to underrepresented sequence at low coverage (as low as 2x). 
In some examples of simple metagenomes,  we observed that the normal (isolate) 
Flye mode assembled more contigious bacterial
consensus sequence, while the metagenome mode was slightly more fragmented, but
revealed strain mixtures. For relatively complex metagenome `--meta` mode
is the recommended way.

### Contig assembly coverage

Typically, assemblies of large genomes at high coverage require
a lot of RAM. For high coverage assemblies, you can reduce memory usage
by using only a subset of longest reads for initial contig extension
stage (usually, the memory bottleneck). The parameter ```--asm-coverage```
specifies the target coverage of the longest reads. For a typicall assembly, 40x
is enough to produce good initial contigs. Regardless of this parameter,
all reads will be later used for repeat graph analysis.

### Number of polishing iterations

Flye first constructs a draft assembly of the genome, which is a 
concatenation of a collection of raw read segments. In the end, the draft assembly
is polished into a high quality sequence. By default, Flye runs one polishing 
iteration. The number could be increased, which might correct a small number of additional
errors (due to improvements on how reads may align to the corrected assembly; 
especially for ONT datasets). If the parameter is set to 0, the polishing will
not be performed.

### Re-starting from a particular assembly stage

Use ```--resume``` to resume a previous run of the assembler that may have terminated
prematurely (using the same output directory). 
The assembly will continue from the last previously completed step.

You might also resume from a particular stage with ```--resume-from stage_name```,
where ```stage_name``` is a choice of ```assembly, consensus, repeat, polishing```.
For example, you might supply different sets of reads for different stages.

## <a name="output"></a> Flye output

The main output files are:

* `assembly.fasta` - Final assembly. Contains contigs and possibly scaffolds (see below).
* `assembly_graph.{gfa|gv}` - Final repeat graph. Note that the edge sequences might be
different (shorter) than contig sequences, because contigs migh include multiple
graph edges (see below).
* `assembly_info.txt` - Extra information about contigs (such as length or coverage).

Each contig is formed by a single unique graph edge. If possble, unique contigs are
extended with the sequence from flanking unresolved repeats on the graph. Thus,
a contig fully contains the corresponding graph edge (with the same id), but might
be longer then this edge. This is somewhat similar to unitig-contig relation
in OLC assemblers. In a rare case when a repetitive graph edge is not covered by 
the set of "extended" contigs, it will be also output in the assembly file.

Sometimes it is possible to further order contigs into scaffolds based on the 
repeat graph structure. These ordered contigs will be output as a part of scaffold
in the assembly file (with a `scaffold_` prefix).  Since it is hard to give a reliable estimate of the
gap size, those gaps are represented with the default 100 Ns. `assembly_info.txt`
file (below) contains additional information about how scaffolds were formed.

Extra information about contigs/scaffolds is output into the ```assembly_info.txt``` file.
It is a tab-delimeted table with the columns as follows:

* Contig/scaffold id
* Length
* Coverage
* Is circular (representing circular sequence, such as bacterial chromosome or plasmid)
* Is repetitive (represents repeated, rather than unique sequence)
* Multiplicity (inferred multiplicity based on coverage)
* Graph path (repeat grapth path corresponding to this contig/scaffold).
Scaffold gaps are marked with "??" symbols, and '*' symbol denotes a
terminal graph node.

## <a name="graph"></a> Repeat graph

The final repeat (assembly) graph is output into the ```assembly_graph.gv``` file
(also available in .gfa format). It could be visualized using [Graphviz](https://graphviz.gitlab.io/): 
```dot -Tpng -O assembly_graph.gv```. The edges in this graph 
represent genomic sequences, and nodes define the sequence junstions
The genome traverses this graph (in an unknown way) so as each unique edge 
is covered exactly once. The unresolved genomic repeats
are collapsed into the corresponding edges in the graph
(therefore genome structure remain umbigious).

<p align="center">
  <img src="graph_example.png" alt="Graph example"/>
</p>

Above is an example of a final assembly graph of a bacterial genome.
Each edge is labeled with its id, length and coverage. Repetitive edges are shown
in color, while unique edges are black. The clusters of adjacent repeats are shown with the 
same color. Note that each edge is represented in two copies: forward and
reverse complement (marked with +/- signs), therefore the entire genome is
represented in two copies as well. Sometimes (as in this example), forward and reverse-complement
components are clearly separated, but often they form a single connected component
(in case if the genome contain unresolved inverted repeats).

In this example, there are two unresolved repeats: (i) a red repeat of multiplicity two
and length 35k and (ii) a green repeat cluster of multiplicity three and length 34k - 36k.
As the repeats remained unresolved, there are no reads in the dataset that cover
those repeats in full.

Initial assembly graph state (before repeats were resolved) could be found in
the ```20-repeat/graph_before_rr.gv``` file. Additionally, ```.gfa``` versions of
assembly graphs could be found in ```20-repeat``` directory.



## <a name="performance"></a> Running time and memory requirements


Typically, assembly of a bacteria or small eukaryote coverage takes less than half an hour 
on a modern desktop. C. elegans can be assembled within a few hours on a 
computational node. The more detailed benchmarks are below.

|Genome        | Size   | Coverage | CPU time  | RAM    |
|--------------|--------|----------|-----------|--------|
| E.coli       | 5 Mb   | 50       | 170 m     | 2 Gb   |
| S.cerevisiae | 12 Mb  | 30       | 180 m     | 7 Gb   |
| C.elegans    | 100 Mb | 30       | 160 h     | 23 Gb  |
| H.sapiens    | 3 Gb   | 30       | 5000 h    | 500 Gb |


## <a name="algorithm"></a> Algorithm Description

This is a brief description of the Flye algorithm. Please refer to the manuscript
for more detailed information. The draft contig extension is organized as follows:

* Kmer counting / erroneous kmer pre-filtering
* Solid kmer selection (kmers with sufficient frequency, which are unlikely to be erroneous)
* Contig extension. The algorithm starts from a single read and extends it
  with a next overlapping read (overlaps are dynamically detected using the selected
  solid k-mers).

Note that we do not attempt to resolve repeats at this stage, thus
the reconstructed contigs might contain misassemblies. 
Flye then aligns the reads on these draft contigs using minimap2 and
calls a consensus. Afterwards, Flye performs repeat analysis as follows:

* Repeat graph is constructed from the (possibly misassembled) contigs
* In this graph all repeats longer than minimum overlap are collapsed
* The algorithm resolves repeats using the read information and graph structure
* The unbranching paths in the graph are output as contigs

After resolving bridged repeats, Trestle module attemts to reolve simple unbridged
repeats (of multiplicity 2) using the heterogenities between repeat copies.
Finally, Flye performs polishing of the resulting assembly
to correct the remaining errors:

* Alignment of all reads to the current assembly using minimap2
* Partition the alignment into into mini-alignments (bubbles)
* Error correction of each bubble using a maximum likelihood approach

The polishing steps could be repeated, which might slightly increase quality for some datasets.
