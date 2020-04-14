<p align="center"><img src="images/logo.png" alt="Minipolish" width="600"></p>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3752203.svg)](https://doi.org/10.5281/zenodo.3752203)


## Table of contents

* [Introduction](#introduction)
* [Requirements](#requirements)
* [Installation](#installation)
* [Method](#method)
* [Quick usage](#quick-usage)
* [Full usage](#full-usage)
* [Citation](#citation)
* [License](#license)




## Introduction

[Miniasm](https://github.com/lh3/miniasm) is a great long-read assembly tool: straight-forward, effective and very fast. However, it does not include a polishing step, so its assemblies have a high error rate – they are essentially made of stitched-together pieces of long reads.

[Racon](https://github.com/isovic/racon) is a great polishing tool that can be used to clean up assembly errors. It's also very fast and well suited for long-read data. However, it operates on FASTA files, not the [GFA graphs](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) that miniasm makes.

That's where Minipolish comes in. With a single command, it will use Racon to polish up a miniasm assembly, while keeping the assembly in graph form.

It also takes care of some of the other nuances of polishing a miniasm assembly:
* Adding read depth information to contigs
* Fixing sequence truncation that can occur in Racon
* Adding circularising links to circular contigs if not already present (so they display better in [Bandage](https://github.com/rrwick/Bandage))
* 'Rotating' circular contigs between polishing rounds to ensure clean circularisation




## Requirements

Minipolish assumes that you have [minimap2](https://github.com/lh3/minimap2) and [Racon](https://github.com/isovic/racon) installed and available in your PATH. If you can run `minimap2 --version` and `racon --version` on the command line, you should be good to go!

You'll need Python 3.6 or later to run Minipolish (check with `python3 --version`). The only Python package requirement is [Edlib](https://github.com/Martinsos/edlib/tree/master/bindings/python). If you don't already have this package, it will be installed as part of the Minipolish installation process. You'll also need [pytest](https://docs.pytest.org/en/latest/) if you want to run Minipolish's unit tests.




## Installation

### Install from source

You can install Minipolish using [pip](https://pypi.org/project/pip/), either from a local copy:
```bash
git clone https://github.com/rrwick/Minipolish.git
pip3 install ./Minipolish
minipolish --help
```

or directly from GitHub:
```bash
pip3 install git+https://github.com/rrwick/Minipolish.git
minipolish --help
```

If these installation commands aren't working for you (e.g. an error message like `Command 'pip3' not found` or `command 'gcc' failed with exit status 1`) then check out the [installation issues page on the Badread wiki page](https://github.com/rrwick/Badread/wiki/Installation-issues) (a different tool of mine but this wiki page covers the same problems).


### Run without installation

Minipolish can also be run directly from its repository by using the `minipolish-runner.py` script, no installation required:

```bash
git clone https://github.com/rrwick/Minipolish.git
Minipolish/minipolish-runner.py -h
```

If you run Minipolish this way, it's up to you to make sure that [Edlib](https://github.com/Martinsos/edlib/tree/master/bindings/python) is installed for your Python.




## Method

### Step 1: initial Racon polish with constituent reads

Miniasm's assembled contigs are made up of pieces of long reads and therefore has a high error rate – probably around 90% or so, depending on the input reads.

The miniasm GFA file indicates specifically which reads contributed to each contig on the `a` lines. For example:
```
a       utg000001c      0       1834c7d5-151e-d9af-fe1d-6bd9f68d355e:19-126885  +       31415
a       utg000001c      31415   28c30dac-ce92-b0f7-af70-4fd95c448f5b:32-107841  +       4527
a       utg000001c      35942   668e8adb-6b7d-67ba-684f-44d78fc3fe32:85-164862  +       14938
a       utg000001c      50880   0ed28b14-d384-ee7b-d12f-88512a2a829c:43-218068  +       81190
```

Therefore, the first thing Minipolish does is to run Racon on each contig independently, only using the reads which were used to create that contig. This step typically quite fast because it does not involve high read depths, and it can bring the percent identity up to the high 90s.


### Step 2: full Racon polish rounds

Now that the assembly is in better shape, Minipolish does full Racon-polishing rounds – aligning the full read set to the whole assembly and getting a Racon consensus. The default number of polishing rounds is two, but this is configurable with the `--rounds` option.

Minipolish does two things here to ensure that contigs can circularise cleanly. First, it repairs sequence ends as Racon can sometimes truncate them. I.e. if Racon dropped a handful of bases from the start or end of a contig, Minipolish will put them back on. Second, it rotates (i.e. changes the starting position) of circular contigs between polishing rounds. If all goes well, this means that the first base of a circular contig immediately follows the last base – clean circularisation.


### Step 3: contig read depth

Minipolish finishes by doing one more read-to-assembly alignment, this time not to polish but to calculate read depths. These depths are added to the GFA line for each contig (e.g. `dp:f:77.179`) and they will be recognised if the graph is loaded in [Bandage](https://github.com/rrwick/Bandage).


### CIGARs

It is important to note here something that Minipolish does _not_ do: change/fix the CIGAR strings indicating contig overlap. While circular contigs will be connected with an overlap-free link (i.e. a CIGAR of `0M`), links between linear contigs will have overlap.

For example, if miniasm created a graph with this link...
```
L	utg000001l	+	utg000020l	+	77073M	SD:i:86773
```
...then that link will have the same CIGAR in the polished assembly. However, since the sequence was polished, the overlap value (77073) will no longer be quite right.

So take CIGAR overlaps between polished contigs with a grain of salt. They will still indicate the _approximate_ amount of overlap, not the _exact_ amount.



## Quick usage

First use minimap2 and miniasm to make an assembly, then polish it with Minipolish:
```
minimap2 -t 8 -x ava-ont long_reads.fastq.gz long_reads.fastq.gz > overlaps.paf
miniasm -f long_reads.fastq.gz overlaps.paf > assembly.gfa
minipolish -t 8 long_reads.fastq.gz assembly.gfa > polished.gfa
```

This repo contains a small Bash script (`miniasm_and_minipolish.sh`) to do those three steps in a single command. It takes two positional arguments: the long reads file and the number of threads:
```
miniasm_and_minipolish.sh long_reads.fastq.gz 8 > polished.gfa
```



## Full usage

```
usage: minipolish [-t THREADS] [--rounds ROUNDS] [--pacbio] [-h] [--version] reads assembly

Minipolish

Positional arguments:
  reads                          Long reads for polishing (FASTA or FASTQ format)
  assembly                       Miniasm assembly to be polished (GFA format)

Settings:
  -t THREADS, --threads THREADS  Number of threads to use for alignment and polishing (default: 12)
  --rounds ROUNDS                Number of full Racon polishing rounds (default: 2)
  --pacbio                       Use this flag for PacBio reads to make Minipolish use the map-pb
                                 Minimap2 preset (default: assumes Nanopore reads and uses the map-ont
                                 preset)

Other:
  -h, --help                     Show this help message and exit
  --version                      Show program's version number and exit
```



## Citation

If you use Minipolish in your research, you can cite the following paper in which it was introduced:

[Wick RR, Holt KE. Benchmarking of long-read assemblers for prokaryote whole genome sequencing. F1000Research. 2019;8(2138).](https://f1000research.com/articles/8-2138)



## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
