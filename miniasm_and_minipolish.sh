#!/usr/bin/env bash

# This script combines minimap2, miniasm and minipolish into a single command
# for easier execution.

# It takes two positional arguments:
#  1) a long read file
#  2) the number of threads to use
#  3) clr, ccs, or no for pabcio map-pb, asm20, or Nanopore reads (see README.md)


# Create temporary intermediate files.
overlaps=$(mktemp)".paf"
unpolished_assembly=$(mktemp)".gfa"

if [ "$3" == "clr" ] || [ "$3" == "ccs" ]
then
    preset="ava-pb"
elif [ "$3" == "no" ]
then
    preset="ava-ont"
else
    echo "Must enter clr, ccs, or no"
fi

# Find read overlaps with minimap2.
minimap2 -x ${preset} -t "$2" "$1" "$1" > "$overlaps"

# Run miniasm to make an unpolished assembly.
miniasm -f "$1" "$overlaps" > "$unpolished_assembly"

# Polish the assembly with minipolish, outputting the result to stdout.
minipolish --pacbio "$3" --threads "$2" "$1" "$unpolished_assembly"

# Clean up.
rm "$overlaps" "$unpolished_assembly"
