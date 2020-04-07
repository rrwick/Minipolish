#!/usr/bin/env bash

set -e

# This script combines minimap2, miniasm and minipolish into a single command
# for easier execution.

# It takes two positional arguments:
#  1) a long read file
#  2) the number of threads to use
#  3) clr, ccs, or no for pabcio map-pb, asm20, or Nanopore reads (see README.md)


# Create temporary intermediate files.
overlaps=$(mktemp)".paf"
unpolished_assembly=$(mktemp)".gfa"

if [ -z "$*" ]; then echo -e "\nUsage: ./miniasm_and_minipolish.sh long_reads.fastq.file threads pacbio-setting aligner > output.gfa \nUsage example: ./miniasm_and_minipolish.sh long_reads.fastq.gz 8 no minimap2 > polished.gfa\n" >&2;exit 1; fi

re='^[1-9]+$'
if ! [[ "$2" =~ $re ]]
then
   echo "Threads is not an integer >0 " >&2; exit 1
fi

#if [ "$4" != "winnowmap" ]
if [ "$4" != "winnowmap" ] && [ "$4" != "minimap2" ]
then
    echo "Must enter minimap2 or winnowmap" >&2; exit 1
fi

if [ "$3" == "clr" ] || [ "$3" == "ccs" ]
then
    preset="ava-pb"
elif [ "$3" == "no" ]
then
    preset="ava-ont"
else
    echo "Must enter clr, ccs, or no" >&2; exit 1
fi

# Find read overlaps with minimap2.
minimap2 -x ${preset} -t "$2" "$1" "$1" > "$overlaps"

# Run miniasm to make an unpolished assembly.
miniasm -f "$1" "$overlaps" > "$unpolished_assembly"

# Polish the assembly with minipolish, outputting the result to stdout.
minipolish --aligner "$4" --pacbio "$3" --threads "$2" "$1" "$unpolished_assembly"

# Clean up.
rm "$overlaps" "$unpolished_assembly"
