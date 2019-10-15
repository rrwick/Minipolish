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

import minipolish.racon
import minipolish.misc


def load_seq(name):
    fasta_seqs = minipolish.misc.load_fasta('test/test_racon.fasta')
    fasta_seq = [x for x in fasta_seqs if x[0] == name]
    return fasta_seq[0][1]


def test_fix_ends_1():
    before_seq = load_seq('test_1_before')
    after_seq = load_seq('test_1_after')
    fixed_seq = load_seq('test_1_fixed')
    result = minipolish.racon.fix_sequence_ends_one_pair(before_seq, after_seq)
    assert result == fixed_seq


def test_fix_ends_2():
    before_seq = load_seq('test_2_before')
    after_seq = load_seq('test_2_after')
    fixed_seq = load_seq('test_2_fixed')
    result = minipolish.racon.fix_sequence_ends_one_pair(before_seq, after_seq)
    assert result == fixed_seq


def test_fix_ends_3():
    before_seq = load_seq('test_3_before')
    after_seq = load_seq('test_3_after')
    fixed_seq = load_seq('test_3_fixed')
    result = minipolish.racon.fix_sequence_ends_one_pair(before_seq, after_seq)
    assert result == fixed_seq


def test_fix_ends_4():
    before_seq = load_seq('test_4_before')
    after_seq = load_seq('test_4_after')
    fixed_seq = load_seq('test_4_fixed')
    result = minipolish.racon.fix_sequence_ends_one_pair(before_seq, after_seq)
    assert result == fixed_seq
