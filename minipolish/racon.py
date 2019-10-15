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

import edlib
import subprocess

from .misc import print_stderr, count_reads, load_fasta


RACON_PATCH_SIZE = 200


def run_racon(name, read_filename, unpolished_filename, threads, tmp_dir):
    polished_filename = tmp_dir / (name + '_polished.fasta')

    if name is None:
        name = unpolished_filename
    read_count = count_reads(read_filename)
    if read_count <= 1:
        print_stderr(f'  skipping Racon for {name} (not enough reads)')
        return False

    print_stderr(f'  running Racon on {name}:')
    print_stderr(f'    input: {unpolished_filename}')
    print_stderr(f'    reads: {read_filename} ({read_count} reads)')
    print_stderr(f'    output: {polished_filename}')

    # Align with minimap2
    command = ['minimap2', '-t', str(threads), '-x', 'map-ont', unpolished_filename, read_filename]
    alignments = tmp_dir / (name + '.paf')
    minimap2_log = tmp_dir / (name + '_minimap2.log')
    with open(alignments, 'wt') as stdout, open(minimap2_log, 'w') as stderr:
        subprocess.call(command, stdout=stdout, stderr=stderr)

    # Polish with Racon
    command = ['racon', '-t', str(threads), read_filename, str(alignments), unpolished_filename]
    racon_log = tmp_dir / (name + '_racon.log')
    with open(polished_filename, 'wt') as stdout, open(racon_log, 'w') as stderr:
        subprocess.call(command, stdout=stdout, stderr=stderr)

    fixed_seqs = fix_sequence_ends(unpolished_filename, polished_filename)
    return fixed_seqs


def fix_sequence_ends(before_fasta, after_fasta):
    """
    Racon can sometimes drop the ends of sequences when polishing, so this function does some
    alignments and patches this up when it happens.
    """
    before_contigs = load_fasta(before_fasta)
    after_contigs = load_fasta(after_fasta)

    # There should be a one-to-one relationship between the before and after contig names, with the
    # caveat that a contig may be missing in the after group.
    before_names = [x[0] for x in before_contigs]
    after_names = [x[0] for x in after_contigs]
    assert all(a in before_names for a in after_names)

    fixed_seqs = {}
    for before_name, before_seq in before_contigs:
        if before_name not in after_names:
            fixed_seq = ''
        else:
            after_seq = [x for x in after_contigs if x[0] == before_name][0][1]
            fixed_seq = fix_sequence_ends_one_pair(before_seq, after_seq)
        fixed_seqs[before_name] = fixed_seq
    return fixed_seqs


def fix_sequence_ends_one_pair(before_seq, after_seq):
    # We will grab a smaller chunk of the 'after' sequence to semi-globally align into a larger
    # chunk of the 'before' sequence.
    before_size = RACON_PATCH_SIZE * 2
    after_size = RACON_PATCH_SIZE

    # Do the alignment for the beginning of the sequence.
    before_start = before_seq[:before_size]
    after_start = after_seq[:after_size]
    result = edlib.align(after_start, before_start, mode='HW', task='path')
    start_pos = result['locations'][0][0]
    additional_start_seq = before_start[:start_pos]

    # And do the alignment for the end of the sequence.
    before_end = before_seq[-before_size:]
    after_end = after_seq[-after_size:]
    result = edlib.align(after_end, before_end, mode='HW', task='path')
    end_pos = result['locations'][0][1] + 1
    additional_end_seq = before_end[end_pos:]

    return additional_start_seq + after_seq + additional_end_seq
