Flye manual
===========

Table of Contents
-----------------

- [Quick usage](#quickusage)
- [Examples](#examples)
- [Supported Input Data](#inputdata)
- [Parameter Descriptions](#parameters)
- [Output Files](#output)
- [Running Time and Memory Requirements](#performance)
- [Algorithm Description](#algorithm)


## <a name="quickusage"></a> Quick usage

    usage: flye (--pacbio-raw | --pacbio-corr | --nano-raw |
                 --nano-corr | --subassemblies) file1 [file_2 ...]
                 --genome-size size --out-dir dir_path [--threads int]
                 [--iterations int] [--min-overlap int] [--resume]
                 [--debug] [--version] [--help]
    
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
                            high-quality contig-like input
      -g size, --genome-size size
                            estimated genome size (for example, 5m or 2.6g)
      -o path, --out-dir path
                            Output directory
      -t int, --threads int
                            number of parallel threads (default: 1)
      -i int, --iterations int
                            number of polishing iterations (default: 1)
      -m int, --min-overlap int
                            minimum overlap between reads (default: 5000)
      --resume              resume from the last completed stage
      --resume-from stage\_name
                            resume from a custom stage
      --debug               enable debug output
      -v, --version         show program's version number and exit



## <a name="examples"></a> Examples

You can try Flye assembly on these ready-to-use datasets:

### E. coli P6-C4 PacBio data

The original dataset is available at the PacBio website 
(https://github.com/PacificBiosciences/DevNet/wiki/E.-coli-Bacterial-Assembly).
We coverted the raw 'bas.h5' file to the FASTA format for the convenience.

    wget https://github.com/fenderglass/datasets/raw/master/pacbio/E.coli_PacBio_40x.fasta
	flye --pacbio-raw E.coli_PacBio_40x.fasta --out-dir out_pacbio --genome-size 5m --threads 4

with '5m' being the expected genome size, the threads argument being optional 
(you may adjust it for your environment), and 'out\_pacbio' being the directory
where the assembly results will be placed.

### E. coli Oxford Nanopore Technologies data

The dataset was originally released by the Loman lab 
(http://lab.loman.net/2015/09/24/first-sqk-map-006-experiment/).

    wget https://github.com/fenderglass/datasets/raw/master/ont/Loman_E.coli_MAP006-1_2D_50x.fasta
	abruijn --nano-raw Loman_E.coli_MAP006-1_2D_50x.fasta --out-dir out_nano --genome-size 5m --threads 4


## <a name="inputdata"></a> Supported Input Data

Input reads could be in FASTA or FASTQ format, uncompressed
or compressed with gz. Currenlty, raw and corrected reads
from PacBio and ONT are supported. Additionally, --subassemblies
option does a consensus assembly of high-quality input contigs.
You may specify multiple fles with reads (separated by spaces).
Mixing different read types is not yet supported.


### PacBio data

Flye was tested on the P6-C4 chemistry data with error rates 11-15%.
Typically, a bacterial WGS project with 20x-30x+ coverage can be assembled 
into a single, structurally concordant contig for each chromosome. However, to get the best 
nucleotide-level quaity, you might need deeper coverage. For a 50x E. coli dataset,
Flye makes roughly 30 errors (single nucleotide insertions/deletions). 

If the coverage of the bacterial dataset is significantly higher than 100x, you
might consider decreasing the coverage by filtering out shorter reads - this
should reduce the running time and memory footprint without affecting the quality. 
However, for some complicated genomes (such as those enriched with mosaic tandem repeats),
incorporating all available reads may be preferred for obtaining a more accurate structural assembly.


### Oxford Nanopore Technologies data

We performed our benchmarks with ONT 2D pass reads with error rates 13-19%.
Due to the increased error rate, you might need deeper coverage 
to get a complete chromosome assembly (60x as in the E. coli example above). For low coverage datasets
(<30x) or datasets with shorter read length you might need to adjust some parameters
(as described below) to get complete chromosomes. Due to the biased error pattern, 
per-nucleotide accuracy is usually lower for ONT data than with PacBio data, especially in homopolymer regions.

### Consensus of multiple input contig sets

--subassemblies input mode generates a consensus of multiple high quality contig assemblies
(such as produced by different short read assemblers). The expected error rate
is similar to the one in corrected PacBio or ONT reads. When using this option,
consider decresing the minimum overlap parameter (for example, 1000 instead of 5000).
You might also want to skip the polishing stage with '--iterations 0' argument
(however, it might still be helpful).


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
Since the algorithm is based on approximate overlaps (without computing exact alignment), 
we require relatively long overlaps (5000 by default) - which is suitable for the most datasets
with realtively good read length and coverage. However, you may decrease this parameter for better contiguity
on low-coverage datasets or datasets with shorter mean read length (such as produced with older
PacBio chemistry).

### Number of polishing iterations

Flye first constructs a draft assembly of the genome, which is a 
concatenation of a collection of raw read segments. In the end, the draft assembly
is polished into a high quality sequence. By default, Flye runs one polishing 
iteration. The number could be increased, which might correct a small number of additional
errors (due to improvements on how reads may align to the corrected assembly; 
especially for ONT datasets). If the parameter is set to 0, the polishing will
not be performed (in case you want to use an external polisher).

### Resuming existing jobs

Use --resume to resume a previous run of the assembler that may have terminated
prematurely. The assembly will continue from the last previously completed step.


## <a name="output"></a> Output files

Flye outputs final contigs/scaffolds along with the simplified assembly graph
and extra information about the assembly fragments. All files will be in the root
of the specified output directory.

### TBD

## <a name="performance"></a> Running time and memory requirements


Typically, assembly of a bacteria or small eukaryote coverage takes less than half an hour 
on a modern desktop. C. elegans can be assembled within a few hours on a 
computational node. The more detailed benchmarks are below.

|Genome        | Size   | Coverage | CPU\_time | RAM    |
|--------------|--------|----------|-----------|--------|
| E.coli       | 5 Mb   | 50       | 170 m     | 2 Gb   |
| S.cerevisiae | 12 Mb  | 30       | 180 m     | 7 Gb   |
| C.elegans    | 100 Mb | 30       | 160 h     | 23 Gb  |
| H.sapiens    | 3 Gb   | 30       | 8400 h    | 700 Gb |


## <a name="algorithm"></a> Algorithm Description

This is a brief description of the Flye algorithm. Plase refer to the manuscript
for more detailed information. The assembly pipeline is organized as follows:

* Kmer counting / erronoeus kmer pre-filtering
* Solid kmer selection (kmers with sufficient frequency, which are unlikely to be erroneous)
* Finding read overlaps based on the A-Bruijn graph
* Detection of chimeric sequences
* Contig assembly by read extension

The resulting contig assembly is now simply a concatenation of read parts 
and is error-prone. Flye then aligns the reads on the draft contigs using minimap2 and
calls a rough consensus. Afterwards, the algorithm performs additional repeat analysis
as follows:

* Repeat graph is reconstructed from the assembled sequence
* In this graph all repeats longer than minimum overlap are collapsed
* The algorithm resolves repeats using the read information and graph structure
* The unbranching paths in the graph are output as contigs

Finally, Flye performs polishing of the resulting assembly
to correct the remaining errors:

* Alignment of all reads to the current assembly using minimap2
* Selection of solid regions
* Partition the total alignment of all reads into mini-alignments (bubbles)
* Error correction of each bubble using a maximum likelihood approach

The polishing steps could be repeated, which might slightly increase quality for some datasets.
