#!/usr/bin/env python3
"""
Copyright 2019 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Minipolish

This file is part of Minipolish. Minipolish is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Minipolish is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Minipolish.
If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import collections
import pathlib
import random
import tempfile

from .assembly_graph import load_gfa
from .help_formatter import MyParser, MyHelpFormatter
from .misc import print_stderr, iterate_fastq, get_default_thread_count
from .racon import run_racon


__version__ = '0.1.0'


def get_arguments(args):
    parser = MyParser(description='Minipolish', add_help=False, formatter_class=MyHelpFormatter)

    required_args = parser.add_argument_group('Positional arguments')
    required_args.add_argument('reads', type=str,
                               help='Long reads for polishing (FASTA or FASTQ format)')
    required_args.add_argument('assembly', type=str,
                               help='Miniasm assembly to be polished (GFA format)')

    setting_args = parser.add_argument_group('Settings')
    setting_args.add_argument('-t', '--threads', type=int, default=get_default_thread_count(),
                              help='Number of threads to use for alignment and polishing')
    setting_args.add_argument('--rounds', type=int, default=2,
                              help='Number of full Racon polishing rounds')

    other_args = parser.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version',
                            version='Minipolish v' + __version__,
                            help="Show program's version number and exit")

    args = parser.parse_args(args)
    return args


def main(args=None):
    args = get_arguments(args)
    # TODO: check for Racon and minimap2
    random.seed(0)
    graph = load_gfa(args.assembly)
    initial_polish(graph, args.reads, args.threads)
    # TODO: full rounds of polishing where all reads are used to polish the assembly.
    # TODO: assign depths to each segment
    # TODO: redo the link overlaps
    # TODO: output the polished GFA to stdout.


def initial_polish(graph, read_filename, threads):
    # This first round of polishing is done on a per-segment basis and only uses reads which are
    # definitely associated with the segment (because the GFA indicated that they were used to
    # make the segment).
    print_stderr('Initial polishing round')

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = pathlib.Path(tmp_dir)
        tmp_dir = pathlib.Path('/Users/ryan/Desktop/Minipolish_test/temp_test')  # TEMP
        save_per_segment_reads(graph, read_filename, tmp_dir)
        for segment in list(graph.segments.values()):
            seg_read_filename = tmp_dir / (segment.name + '.fastq')
            seg_seq_filename = tmp_dir / (segment.name + '.fasta')
            segment.save_to_fasta(seg_seq_filename)
            fixed_seqs = run_racon(segment.name, seg_read_filename, seg_seq_filename, threads,
                                   tmp_dir)
            fixed_seq = fixed_seqs[segment.name]
            if len(fixed_seq) > 0:
                segment.sequence = fixed_seq
            else:
                del graph.segments[segment.name]

    print_stderr('')


def save_per_segment_reads(graph, read_filename, tmp_dir):
    read_to_segment = collections.defaultdict(list)
    for segment in graph.segments.values():
        for read_name in segment.read_names:
            read_to_segment[read_name].append(segment.name)
    for read_name, seq, qual in iterate_fastq(read_filename):
        if read_name not in read_to_segment:
            continue
        for seg_name in read_to_segment[read_name]:
            seg_read_filename = tmp_dir / (seg_name + '.fastq')
            with open(seg_read_filename, 'at') as seg_read_file:
                seg_read_file.write(f'@{read_name}\n{seq}\n+\n{qual}\n')



if __name__ == '__main__':
    main()
