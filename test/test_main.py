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

import minipolish.__main__
import minipolish.assembly_graph
import pytest


def test_initial_polish_skips_when_no_a_lines():
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = pathlib.Path(tmp_dir)
        gfa_filename = tmp_dir / 'assembly.gfa'
        reads_filename = tmp_dir / 'reads.fastq'
        gfa_filename.write_text('S\tutg000001l\tACGTACGT\n')
        reads_filename.write_text('@read_1\nACGT\n+\nIIII\n@read_2\nTGCA\n+\nIIII\n')
        graph = minipolish.assembly_graph.load_gfa(gfa_filename)
        minipolish.__main__.initial_polish(graph, reads_filename, 1, tmp_dir, 'map-ont')
    assert sorted(graph.segments.keys()) == ['utg000001l']
    assert graph.segments['utg000001l'].sequence == 'ACGTACGT'


def test_pacbio_and_minimap2_preset_are_mutually_exclusive():
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = pathlib.Path(tmp_dir)
        gfa_filename = tmp_dir / 'assembly.gfa'
        reads_filename = tmp_dir / 'reads.fastq'
        gfa_filename.write_text('S\tutg000001l\tACGTACGT\n')
        reads_filename.write_text('@read_1\nACGT\n+\nIIII\n')
        with pytest.raises(SystemExit) as e:
            minipolish.__main__.get_arguments(['--pacbio', '--minimap2-preset', 'map-hifi',
                                               str(reads_filename), str(gfa_filename)])
    assert e.type == SystemExit
    assert e.value.code != 0
