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

import pytest

import minipolish.alignment


def test_alignment_1():
    alignment_line = ['read_name',     # read name
                      '1000',          # read length
                      '100', '900',    # read start/end
                      '+',             # strand
                      'ref_name',      # ref name
                      '10000',         # ref length
                      '2000', '2800',  # ref start/end
                      '700',           # match count
                      '850',           # alignment length
                      '255']           # mapping quality
    alignment_line = '\t'.join(alignment_line) + '\n'
    a = minipolish.alignment.Alignment(alignment_line)
    assert a.read_name == 'read_name'
    assert a.read_length == 1000
    assert a.read_start == 100
    assert a.read_end == 900
    assert a.strand == '+'
    assert a.ref_name == 'ref_name'
    assert a.ref_length == 10000
    assert a.ref_start == 2000
    assert a.ref_end == 2800
    assert a.matching_bases == 700
    assert a.num_bases == 850
    assert a.percent_identity == pytest.approx(100.0 * 700 / 850)
    assert a.get_ref_depth_contribution() == pytest.approx(0.08)


def test_alignment_2():
    with pytest.raises(SystemExit) as e:
        _ = minipolish.alignment.Alignment('this is not PAF format')
    assert e.type == SystemExit
