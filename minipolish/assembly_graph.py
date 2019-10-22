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
import random
import sys

from .log import log, section_header, explanation
from .misc import get_open_func


class AssemblyGraph(object):
    def __init__(self):
        self.segments = {}  # dictionary of segment name -> segment object
        self.links = {}  # dictionary of segment names -> link object

    def add_link(self, link):
        names = (link.name_1 + link.strand_1, link.name_2 + link.strand_2)
        assert names not in self.links
        self.links[names] = link

    def remove_segment(self, seg_name):
        del self.segments[seg_name]
        for link_name in list(self.links.keys()):
            link = self.links[link_name]
            if link.name_1 == seg_name or link.name_2 == seg_name:
                self.links.pop(link_name, None)

    def print_to_stdout(self):
        segment_names = sorted(self.segments.keys())
        for name in segment_names:
            self.segments[name].print_gfa_line_to_stdout()

        link_names = sorted(self.links.keys())
        for name in link_names:
            self.links[name].print_gfa_line_to_stdout()

    def rotate_circular_sequences(self):
        segment_names = sorted(self.segments.keys())
        for name in segment_names:
            if name.endswith('c'):
                positive_link = (name + '+', name + '+')
                negative_link = (name + '-', name + '-')
                assert positive_link in self.links
                assert negative_link in self.links
                assert self.links[positive_link].cigar == '0M'
                assert self.links[negative_link].cigar == '0M'
                segment = self.segments[name]
                if segment.get_length() > 1:
                    rotation = random.randint(1, segment.get_length() - 1)
                    segment.rotate(rotation)
        log()

    def save_to_fasta(self, filename):
        segment_names = sorted(self.segments.keys())
        with open(filename, 'wt') as fasta:
            for name in segment_names:
                segment = self.segments[name]
                fasta.write(f'>{name}\n{segment.sequence}\n')

    def replace_sequences(self, new_seqs):
        for seg_name, new_seq in new_seqs.items():
            self.segments[seg_name].sequence = new_seq

    def set_depths(self, depth_per_contig):
        for seg_name, depth in depth_per_contig.items():
            self.segments[seg_name].depth = depth

    def get_segment_length(self, seg_name):
        return self.segments[seg_name].get_length()

    def build_reverse_links(self):
        """
        Each link in the graph (e.g. utg000001l+ -> utg000002l-) should have a corresponding link
        in the other direction (utg000002l+ -> utg000001l-). This function
        """
        # Group links into their pairs.
        str_to_links = collections.defaultdict(list)
        for link in self.links.values():
            str_to_links[link.get_canonical_link_str()].append(link)

        # For each pair, make the other link if it doesn't already exist.
        for _, links in str_to_links.items():
            assert len(links) == 1 or len(links) == 2
            if len(links) == 1:
                self.add_link(make_reverse_link(links[0]))

        # Sanity check: all links should be paired now.
        str_to_links = collections.defaultdict(list)
        for link in self.links.values():
            str_to_links[link.get_canonical_link_str()].append(link)
        for links in str_to_links.values():
            assert len(links) == 2

    def build_circularising_links(self):
        """
        Any circular contig (i.e. ends with 'c') should have a circularising link with no overlap
        (e.g. utg000001c+ -> utg000001c+ 0M). This function will build those links (in both
        directions if they don't already exist.
        """
        segment_names = sorted(self.segments.keys())
        for name in segment_names:
            if name.endswith('c'):
                link_names_forward = (name + '+', name + '+')
                link_names_reverse = (name + '-', name + '-')
                if link_names_forward not in self.links:
                    # This function assumes that build_reverse_links has already been run, so if
                    # the forward link is missing, the reverse link should be too.
                    assert link_names_reverse not in self.links
                    new_link = Link(f'L\t{name}\t+\t{name}\t+\t0M')
                    self.add_link(new_link)
                    self.add_link(make_reverse_link(new_link))


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

    def print_gfa_line_to_stdout(self):
        print(f'S\t{self.name}\t{self.sequence}\tdp:f:{self.depth:.3f}')

    def get_length(self):
        return len(self.sequence)

    def rotate(self, rotation):
        assert self.name.endswith('c')  # Only circular contigs should be rotated
        log(f'Rotating {self.name} by {rotation:,} bp')
        self.sequence = self.sequence[rotation:] + self.sequence[:rotation]


class Link(object):
    def __init__(self, gfa_line):
        parts = gfa_line.strip().split('\t')
        assert parts[0] == 'L'
        self.name_1 = parts[1]
        self.strand_1 = parts[2]
        self.name_2 = parts[3]
        self.strand_2 = parts[4]
        self.cigar = parts[5]

    def print_gfa_line_to_stdout(self):
        print(f'L\t{self.name_1}\t{self.strand_1}\t{self.name_2}\t{self.strand_2}\t{self.cigar}')

    def get_forward_link_str(self):
        return self.name_1 + self.strand_1 + self.name_2 + self.strand_2

    def get_reverse_link_str(self):
        return self.name_2 + flip_strand(self.strand_2) + self.name_1 + flip_strand(self.strand_1)

    def get_canonical_link_str(self):
        """
        Returns a string describing the link in a consistent direction. I.e. each link in a pair
        should return the same string from this function.
        """
        if self.strand_1 == '+' and self.strand_2 == '+':
            return self.get_forward_link_str()
        elif self.strand_1 == '-' and self.strand_2 == '-':
            return self.get_reverse_link_str()
        elif self.name_1 < self.name_2:
            return self.get_forward_link_str()
        else:
            return self.get_reverse_link_str()


def flip_strand(strand):
    if strand == '+':
        return '-'
    assert strand == '-'
    return '+'


def load_gfa(filename):
    section_header('Loading graph')
    explanation('Loading the miniasm GFA graph into memory.')
    log(filename)
    graph = AssemblyGraph()

    # The constituent reads for segments (GFA 'a' lines) will be stored in this dictionary and
    # then added to the segments at the end of this function. This is so we don't have to assume
    # that 'a' lines come after their corresponding 'S' line (though I expect they always do).
    segment_reads = collections.defaultdict(list)

    with get_open_func(filename)(filename) as gfa:
        for line in gfa:
            if line.startswith('S\t'):
                segment = Segment(line)
                if not (segment.name.endswith('l') or segment.name.endswith('c')):
                    sys.exit(f'Error: contig name ({segment.name}) does not appear to be in a '
                             f'miniasm format')
                graph.segments[segment.name] = segment
            if line.startswith('a\t'):
                segment_name, read_name = parse_a_line(line)
                segment_reads[segment_name].append(read_name)
            if line.startswith('L\t'):
                graph.add_link(Link(line))

    graph.build_reverse_links()
    graph.build_circularising_links()

    for segment_name, read_names in segment_reads.items():
        assert segment_name in graph.segments
        graph.segments[segment_name].read_names = read_names

    seg_count = len(graph.segments)
    base_count = sum(len(s.sequence) for s in graph.segments.values())
    link_count = len(graph.links)
    log(f'  {seg_count:,} segments ({base_count:,} bp)')
    log(f'  {link_count:,} links')
    log()
    return graph


def parse_a_line(line):
    parts = line.strip().split('\t')
    assert parts[0] == 'a'
    segment_name = parts[1]
    read_name = parts[3].rsplit(':', 1)[0]
    return segment_name, read_name


def make_reverse_link(link):
    """
    This function returns a Link object which is the reverse direction of the given link object,
    using the same CIGAR.
    """
    return Link(f'L\t{link.name_2}\t{flip_strand(link.strand_2)}\t'
                f'{link.name_1}\t{flip_strand(link.strand_1)}\t{link.cigar}')
