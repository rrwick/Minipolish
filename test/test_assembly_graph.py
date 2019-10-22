"""
This module contains some tests for Minipolish. To run them, execute `python3 -m pytest` from the
root Minipolish directory.

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

import pathlib
import random
import tempfile

import minipolish.assembly_graph
import minipolish.misc


def test_rotation_1():
    seq = 'TCCAGCGTTG'
    seg = minipolish.assembly_graph.Segment(f'S\tutg000001c\t{seq}')
    assert seg.sequence == 'TCCAGCGTTG'
    seg.rotate(2)
    assert seg.sequence == 'CAGCGTTGTC'
    seg.rotate(2)
    assert seg.sequence == 'GCGTTGTCCA'


def test_rotation_2():
    seq = 'CATATAAGTGTACCCTGCGAATATGGTTCG'
    seg = minipolish.assembly_graph.Segment(f'S\tutg000001c\t{seq}')
    assert seg.sequence == 'CATATAAGTGTACCCTGCGAATATGGTTCG'
    seg.rotate(5)
    assert seg.sequence == 'AAGTGTACCCTGCGAATATGGTTCGCATAT'
    seg.rotate(12)
    assert seg.sequence == 'CGAATATGGTTCGCATATAAGTGTACCCTG'


def test_canonical_link_str():
    link = minipolish.assembly_graph.Link('L\tutg000001c\t+\tutg000001c\t+\t0M')
    assert link.get_canonical_link_str() == 'utg000001c+utg000001c+'

    link = minipolish.assembly_graph.Link('L\tutg000001c\t-\tutg000001c\t-\t0M')
    assert link.get_canonical_link_str() == 'utg000001c+utg000001c+'

    link = minipolish.assembly_graph.Link('L\tutg000001l\t+\tutg000002l\t-\t0M')
    assert link.get_canonical_link_str() == 'utg000001l+utg000002l-'

    link = minipolish.assembly_graph.Link('L\tutg000002l\t+\tutg000001l\t-\t0M')
    assert link.get_canonical_link_str() == 'utg000001l+utg000002l-'

    link = minipolish.assembly_graph.Link('L\tutg000001l\t-\tutg000002l\t+\t0M')
    assert link.get_canonical_link_str() == 'utg000001l-utg000002l+'

    link = minipolish.assembly_graph.Link('L\tutg000002l\t-\tutg000001l\t+\t0M')
    assert link.get_canonical_link_str() == 'utg000001l-utg000002l+'


def test_load_gfa_1():
    """
    Tests a one-linear-contig graph.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        temp_gfa_filename = str(pathlib.Path(tmp_dir) / 'test.gfa')
        with open(temp_gfa_filename, 'wt') as temp_gfa:
            temp_gfa.write('S\tutg000001l\tACGTACGACTACGACTG\n')
        graph = minipolish.assembly_graph.load_gfa(temp_gfa_filename)
    assert len(graph.segments) == 1
    assert len(graph.links) == 0


def test_load_gfa_2():
    """
    Tests a one-circular-contig graph. This graph explicitly includes the circularising link in
    both directions.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        temp_gfa_filename = str(pathlib.Path(tmp_dir) / 'test.gfa')
        with open(temp_gfa_filename, 'wt') as temp_gfa:
            temp_gfa.write('S\tutg000001c\tACGTACGACTACGACTG\n')
            temp_gfa.write('L\tutg000001c\t+\tutg000001c\t+\t0M\n')
            temp_gfa.write('L\tutg000001c\t-\tutg000001c\t-\t0M\n')
        graph = minipolish.assembly_graph.load_gfa(temp_gfa_filename)
    assert len(graph.segments) == 1
    assert len(graph.links) == 2
    assert ('utg000001c+', 'utg000001c+') in graph.links
    assert ('utg000001c-', 'utg000001c-') in graph.links


def test_load_gfa_3():
    """
    Tests a one-circular-contig graph. This graph explicitly includes the circularising link in
    only one direction. The other direction should be made automatically.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        temp_gfa_filename = str(pathlib.Path(tmp_dir) / 'test.gfa')
        with open(temp_gfa_filename, 'wt') as temp_gfa:
            temp_gfa.write('S\tutg000001c\tACGTACGACTACGACTG\n')
            temp_gfa.write('L\tutg000001c\t+\tutg000001c\t+\t0M\n')
        graph = minipolish.assembly_graph.load_gfa(temp_gfa_filename)
    assert len(graph.segments) == 1
    assert len(graph.links) == 2
    assert ('utg000001c+', 'utg000001c+') in graph.links
    assert ('utg000001c-', 'utg000001c-') in graph.links


def test_load_gfa_4():
    """
    Tests a one-circular-contig graph. This graph does not include the circularising link. It
    should be automatically made based on the 'c' in the contig name.
    """
    with tempfile.TemporaryDirectory() as tmp_dir:
        temp_gfa_filename = str(pathlib.Path(tmp_dir) / 'test.gfa')
        with open(temp_gfa_filename, 'wt') as temp_gfa:
            temp_gfa.write('S\tutg000001c\tACGTACGACTACGACTG\n')
        graph = minipolish.assembly_graph.load_gfa(temp_gfa_filename)
    assert len(graph.segments) == 1
    assert len(graph.links) == 2
    assert ('utg000001c+', 'utg000001c+') in graph.links
    assert ('utg000001c-', 'utg000001c-') in graph.links


def test_remove_segment():
    with tempfile.TemporaryDirectory() as tmp_dir:
        temp_gfa_filename = str(pathlib.Path(tmp_dir) / 'test.gfa')
        with open(temp_gfa_filename, 'wt') as temp_gfa:
            temp_gfa.write('S\tutg000001l\tACGTACGACTACGACTG\n')
            temp_gfa.write('S\tutg000002l\tACGTACGACTACGACTG\n')
            temp_gfa.write('S\tutg000003l\tACGTACGACTACGACTG\n')
            temp_gfa.write('L\tutg000001l\t+\tutg000002l\t+\t0M\n')
            temp_gfa.write('L\tutg000001l\t+\tutg000003l\t+\t0M\n')
        graph = minipolish.assembly_graph.load_gfa(temp_gfa_filename)
    assert len(graph.segments) == 3
    assert len(graph.links) == 4
    graph.remove_segment('utg000003l')
    assert len(graph.segments) == 2
    assert len(graph.links) == 2


def test_rotate_circular_sequence():
    random.seed(0)
    with tempfile.TemporaryDirectory() as tmp_dir:
        temp_gfa_filename = str(pathlib.Path(tmp_dir) / 'test.gfa')
        with open(temp_gfa_filename, 'wt') as temp_gfa:
            temp_gfa.write('S\tutg000001c\tACGTACGACTACGACTG\n')
        graph = minipolish.assembly_graph.load_gfa(temp_gfa_filename)
    before_rotate = graph.segments['utg000001c'].sequence
    graph.rotate_circular_sequences()
    after_rotate = graph.segments['utg000001c'].sequence
    assert before_rotate != after_rotate
    assert sorted(before_rotate) == sorted(after_rotate)


def test_parse_a_line():
    a_line = 'a\tutg000001c\t0\t1834c7d5-151e-d9af-fe1d-6bd9f68d355e:19-126885\t+\t31415\n'
    segment_name, read_name = minipolish.assembly_graph.parse_a_line(a_line)
    assert segment_name == 'utg000001c'
    assert read_name == '1834c7d5-151e-d9af-fe1d-6bd9f68d355e'
