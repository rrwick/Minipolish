# Minipolish

## Introduction

[Miniasm](https://github.com/lh3/miniasm) is a great long-read assembly tool: straight-forward, effective and extremely fast. However, it does not include a polishing step, so its assemblies have a high error rate – they are essentially made up of stiched-together pieces of long reads.

[Racon](https://github.com/isovic/racon) is a great polishing tool that can be used to clean up assembly errors. It's very fast and well suited for long-read data. However, it operates on FASTA files, not the [GFA graphs](https://github.com/GFA-spec/GFA-spec/blob/master/GFA1.md) that miniasm makes.

That's where this script comes in. It uses Racon to polish up a miniasm assembly, while keeping the assembly in graph form, with a single command.

It also takes care of some of the other nuances of polishing a miniasm assembly:
* Adding read depth information to contigs
* Fixing the sequence truncation that can occur in Racon
* Adding circularising links to circular contigs if not already present (so they display better in [Bandage](https://github.com/rrwick/Bandage))
* 'Rotating' circular contigs between polishing rounds to ensure clean circularisation




## Requirements




## Installation

Since Minipolish is just a single script, no installation is required. You can simply clone it and run it:
```
git clone https://github.com/rrwick/Minipolish
Minipolish/minipolish.py --help
```

If you plan on using it often, you can copy it to someplace in your PATH variable for easier access:
```
git clone https://github.com/rrwick/Minipolish
cp Minipolish/minipolish.py ~/.local/bin
minipolish.py --help
```

If you'd like to double-check that everything works as intended, you can run this repo's [automated tests](test).





## Method




## Quick usage

```
minimap2 -x ava-ont long_reads.fastq.gz long_reads.fastq.gz > overlaps.paf
miniasm -f long_reads.fastq.gz overlaps.paf > assembly.gfa
minipolish.py long_reads.fastq.gz assembly.gfa > polished.gfa
```



## Full usage





## License

[GNU General Public License, version 3](https://www.gnu.org/licenses/gpl-3.0.html)
