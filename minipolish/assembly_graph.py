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

import collections

from .misc import print_stderr, get_open_func


class AssemblyGraph(object):
    def __init__(self):
        self.segments = {}  # dictionary of segment name -> segment object
        self.links = {}  # dictionary of segment names -> link object


class Segment(object):
    def __init__(self, gfa_line):
        parts = gfa_line.strip().split('\t')
        assert parts[0] == 'S'
        self.name = parts[1]
        self.sequence = parts[2]
        self.depth = 0.0
        self.read_names = []

    def save_to_fasta(self, filename):
        with open(filename, 'wt') as fasta:
            fasta.write(f'>{self.name}\n{self.sequence}\n')


class Link(object):
    def __init__(self, gfa_line):
        parts = gfa_line.strip().split('\t')
        assert parts[0] == 'L'
        self.name_1 = parts[1]
        self.strand_1 = parts[2]
        self.name_2 = parts[3]
        self.strand_2 = parts[4]
        self.cigar = parts[5]


def load_gfa(filename):
    print_stderr(f'\nLoading {filename}')
    graph = AssemblyGraph()

    # The constituent reads for segments (GFA 'a' lines) will be stored in this dictionary and
    # then added to the segments at the end of this function. This is so we don't have to assume
    # that 'a' lines come after their corresponding 'S' line (though I expect they always do).
    segment_reads = collections.defaultdict(list)

    with get_open_func(filename)(filename) as gfa:
        for line in gfa:
            if line.startswith('S\t'):
                segment = Segment(line)
                graph.segments[segment.name] = segment
            if line.startswith('a\t'):
                segment_name, read_name = parse_a_line(line)
                segment_reads[segment_name].append(read_name)
            if line.startswith('L\t'):
                link = Link(line)
                names = (link.name_1 + link.strand_1, link.name_2 + link.strand_2)
                assert names not in graph.links
                graph.links[names] = link

    for segment_name, read_names in segment_reads.items():
        assert segment_name in graph.segments
        graph.segments[segment_name].read_names = read_names

    seg_count = len(graph.segments)
    base_count = sum(len(s.sequence) for s in graph.segments.values())
    link_count = len(graph.links)
    print_stderr(f'  {seg_count:,} segments ({base_count:,} bp)')
    print_stderr(f'  {link_count:,} links')
    print_stderr('')
    return graph


def parse_a_line(line):
    parts = line.strip().split('\t')
    assert parts[0] == 'a'
    segment_name = parts[1]
    read_name = parts[3].rsplit(':', 1)[0]
    return segment_name, read_name