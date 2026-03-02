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
import tempfile

import minipolish.racon
import minipolish.misc
import pytest


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


def test_run_racon_one_read_exits_on_no_alignments(monkeypatch):
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = pathlib.Path(tmp_dir)
        read_filename = tmp_dir / 'reads.fastq'
        unpolished_filename = tmp_dir / 'segment.fasta'
        read_filename.write_text('@read_1\nACGT\n+\nIIII\n')
        unpolished_filename.write_text('>segment\nACGTACGT\n')
        monkeypatch.setattr(minipolish.racon.subprocess, 'call', lambda *args, **kwargs: 0)
        with pytest.raises(SystemExit) as e:
            minipolish.racon.run_racon('segment', read_filename, unpolished_filename, 1,
                                       tmp_dir, 'map-ont')
    assert e.type == SystemExit
    assert 'produced no alignments' in str(e.value)


def test_run_racon_no_alignments_exits(monkeypatch):
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = pathlib.Path(tmp_dir)
        read_filename = tmp_dir / 'reads.fastq'
        unpolished_filename = tmp_dir / 'segment.fasta'
        read_filename.write_text('@read_1\nACGT\n+\nIIII\n@read_2\nTGCA\n+\nIIII\n')
        unpolished_filename.write_text('>segment\nACGTACGT\n')
        monkeypatch.setattr(minipolish.racon.subprocess, 'call', lambda *args, **kwargs: 0)
        with pytest.raises(SystemExit) as e:
            minipolish.racon.run_racon('segment', read_filename, unpolished_filename, 1,
                                       tmp_dir, 'map-ont')
    assert e.type == SystemExit
    assert 'produced no alignments' in str(e.value)
