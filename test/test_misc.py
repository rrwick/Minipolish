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

import gzip
import pathlib
import pytest
import tempfile

import minipolish.misc
import minipolish.version


def test_version():
    version_str = minipolish.version.__version__
    assert version_str.count('.') == 2


def test_get_default_thread_count():
    threads = minipolish.misc.get_default_thread_count()
    assert 1 <= threads <= 16


def test_weighted_average_1():
    nums = [80, 90]
    weights = [20, 30]
    average = minipolish.misc.weighted_average(nums, weights)
    assert average == pytest.approx(86.0)


def test_weighted_average_2():
    nums = [80, 90]
    weights = [0, 0]
    average = minipolish.misc.weighted_average(nums, weights)
    assert average == pytest.approx(85.0)


def test_weighted_average_3():
    nums = [80, 90]
    weights = [1, 0]
    average = minipolish.misc.weighted_average(nums, weights)
    assert average == pytest.approx(80.0)


def test_load_fasta():
    with tempfile.TemporaryDirectory() as tmp_dir:
        temp_filename = str(pathlib.Path(tmp_dir) / 'test.fasta')
        with open(temp_filename, 'wt') as temp_fasta:
            temp_fasta.write('>read_1\nACGATCGACT\n')
            temp_fasta.write('>read_2 extra stuff\nGGCGCTCG\n')
        seqs = minipolish.misc.load_fasta(temp_filename)
    assert seqs[0] == ('read_1', 'ACGATCGACT')
    assert seqs[1] == ('read_2', 'GGCGCTCG')


def test_count_lines():
    with tempfile.TemporaryDirectory() as tmp_dir:
        temp_filename = str(pathlib.Path(tmp_dir) / 'test.txt')
        with open(temp_filename, 'wt') as temp_txt:
            temp_txt.write('line_1\nline_2\nline_3\nline_4\n')
        line_count = minipolish.misc.count_lines(temp_filename)
    assert line_count == 4


def test_count_reads_uncompressed():
    with tempfile.TemporaryDirectory() as tmp_dir:
        temp_filename = str(pathlib.Path(tmp_dir) / 'test.fastq')
        with open(temp_filename, 'wt') as temp_fastq:
            temp_fastq.write('@read_1\nACGT\n+\nIIII\n')
            temp_fastq.write('@read_2 extra stuff\nGGCG\n+\nIIII\n')
        read_count = minipolish.misc.count_reads(temp_filename)
    assert read_count == 2


def test_count_reads_gzipped():
    with tempfile.TemporaryDirectory() as tmp_dir:
        temp_filename = str(pathlib.Path(tmp_dir) / 'test.fastq.gz')
        with gzip.open(temp_filename, 'wt') as temp_fastq:
            temp_fastq.write('@read_1\nACGT\n+\nIIII\n')
            temp_fastq.write('@read_2 extra stuff\nGGCG\n+\nIIII\n')
        read_count = minipolish.misc.count_reads(temp_filename)
    assert read_count == 2


def test_get_sequence_file_type_1():
    with tempfile.TemporaryDirectory() as tmp_dir:
        temp_filename = str(pathlib.Path(tmp_dir) / 'test.fasta')
        with open(temp_filename, 'wt') as temp_fasta:
            temp_fasta.write('>read_1\nACGATCGACT\n')
            temp_fasta.write('>read_2 extra stuff\nGGCGCTCG\n')
        file_type = minipolish.misc.get_sequence_file_type(temp_filename)
    assert file_type == 'FASTA'


def test_get_sequence_file_type_2():
    with tempfile.TemporaryDirectory() as tmp_dir:
        temp_filename = str(pathlib.Path(tmp_dir) / 'test.fastq')
        with open(temp_filename, 'wt') as temp_fastq:
            temp_fastq.write('@read_1\nACGT\n+\nIIII\n')
            temp_fastq.write('@read_2 extra stuff\nGGCG\n+\nIIII\n')
        file_type = minipolish.misc.get_sequence_file_type(temp_filename)
    assert file_type == 'FASTQ'


def test_get_sequence_file_type_3():
    with pytest.raises(SystemExit) as e:
        minipolish.misc.get_sequence_file_type('this_file_does_not_exist')
    assert e.type == SystemExit
