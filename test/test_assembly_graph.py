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
