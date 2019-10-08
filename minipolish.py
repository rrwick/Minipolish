#!/usr/bin/env python3
"""
Copyright 2019 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Minipolish

This program is free software: you can redistribute it and/or modify it under the terms of the GNU
General Public License as published by the Free Software Foundation, either version 3 of the
License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not,
see <https://www.gnu.org/licenses/>.
"""

import argparse
import collections
import gzip
import os
import pathlib
import random
import shutil
import subprocess
import sys
import tempfile

__version__ = '0.1.0'


def get_arguments(args):
    parser = MyParser(description='Minipolish', add_help=False,
                      formatter_class=MyHelpFormatter)

    required_args = parser.add_argument_group('Positional arguments')
    required_args.add_argument('reads', type=str,
                               help='Long reads for polishing (FASTA or FASTQ format)')
    required_args.add_argument('assembly', type=str,
                               help='Miniasm assembly to be polished (GFA format)')

    setting_args = parser.add_argument_group('Settings')
    setting_args.add_argument('--rounds', type=int, default=2,
                              help='Number of full Racon polishing rounds')

    other_args = parser.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='Show this help message and exit')
    other_args.add_argument('--version', action='version',
                            version='Minipolish v' + __version__,
                            help="Show program's version number and exit")

    args = parser.parse_args(args)
    return args


def main(args=None):
    args = get_arguments(args)
    # TODO: check for Racon
    random.seed(0)
    graph = load_gfa(args.assembly)
    initial_polish(graph, args.reads)
    # TODO: full rounds of polishing where all reads are used to polish the assembly.
    # TODO: assign depths to each segment
    # TODO: redo the link overlaps
    # TODO: output the polished GFA to stdout.


def initial_polish(graph, read_filename):
    # This first round of polishing is done on a per-segment basis and only uses reads which are
    # definitely associated with the segment (because the GFA indicated that they were used to
    # make the segment).
    print_stderr('Initial polishing round')

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp_dir = pathlib.Path(tmp_dir)
        # tmp_dir = pathlib.Path('/Users/ryan/Desktop/test')  # TEMP!
        save_per_segment_reads(graph, read_filename, tmp_dir)
        for segment in graph.segments.values():
            seg_read_filename = tmp_dir / (segment.name + '.fastq')
            seg_seq_filename = tmp_dir / (segment.name + '.fasta')
            segment.save_to_fasta(seg_seq_filename)
            polished_seq_filename = tmp_dir / (segment.name + '_polished.fasta')
            if run_racon(segment.name, seg_read_filename, seg_seq_filename, polished_seq_filename):
                # TODO: change the segment sequence to the polished sequence
                pass

    print_stderr('')


def save_per_segment_reads(graph, read_filename, tmp_dir):
    read_to_segment = collections.defaultdict(list)
    for segment in graph.segments.values():
        for read_name in segment.read_names:
            read_to_segment[read_name].append(segment.name)
    for read_name, seq, qual in iterate_fastq(read_filename):
        if read_name not in read_to_segment:
            continue
        for seg_name in read_to_segment[read_name]:
            seg_read_filename = tmp_dir / (seg_name + '.fastq')
            with open(seg_read_filename, 'at') as seg_read_file:
                seg_read_file.write(f'@{read_name}\n{seq}\n+\n{qual}\n')


def run_racon(name, read_filename, seq_filename, polished_filename):
    if name is None:
        name = seq_filename
    read_count = count_reads(read_filename)
    if read_count <= 1:
        print_stderr(f'  skipping Racon for {name} (not enough reads)')
        return False

    print_stderr(f'  running Racon on {name}:')
    print_stderr(f'    input: {seq_filename}')
    print_stderr(f'    reads: {read_filename} ({read_count})')
    print_stderr(f'    output: {polished_filename}')

    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    # TODO
    return True


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


def get_compression_type(filename):
    """
    Attempts to guess the compression (if any) on a file using the first few bytes.
    http://stackoverflow.com/questions/13044562
    """
    magic_dict = {'gz': (b'\x1f', b'\x8b', b'\x08'),
                  'bz2': (b'\x42', b'\x5a', b'\x68'),
                  'zip': (b'\x50', b'\x4b', b'\x03', b'\x04')}
    max_len = max(len(x) for x in magic_dict)

    unknown_file = open(str(filename), 'rb')
    file_start = unknown_file.read(max_len)
    unknown_file.close()
    compression_type = 'plain'
    for file_type, magic_bytes in magic_dict.items():
        if file_start.startswith(magic_bytes):
            compression_type = file_type
    if compression_type == 'bz2':
        sys.exit('Error: cannot use bzip2 format - use gzip instead')
    if compression_type == 'zip':
        sys.exit('Error: cannot use zip format - use gzip instead')
    return compression_type


def get_open_func(filename):
    if get_compression_type(filename) == 'gz':
        return gzip.open
    else:  # plain text
        return open


def get_sequence_file_type(filename):
    """
    Determines whether a file is FASTA or FASTQ.
    """
    if not os.path.isfile(filename):
        sys.exit('Error: could not find {}'.format(filename))
    if get_compression_type(filename) == 'gz':
        open_func = gzip.open
    else:  # plain text
        open_func = open
    with open_func(filename, 'rt') as seq_file:
        try:
            first_char = seq_file.read(1)
        except UnicodeDecodeError:
            first_char = ''
    if first_char == '>':
        return 'FASTA'
    elif first_char == '@':
        return 'FASTQ'
    else:
        raise ValueError('File is neither FASTA or FASTQ')


def iterate_fastq(filename):
    if get_sequence_file_type(filename) != 'FASTQ':
        sys.exit('Error: {} is not FASTQ format'.format(filename))
    with get_open_func(filename)(filename, 'rt') as fastq:
        for line in fastq:
            line = line.strip()
            if len(line) == 0:
                continue
            if not line.startswith('@'):
                continue
            name = line[1:].split()[0]
            sequence = next(fastq).strip()
            _ = next(fastq)
            qualities = next(fastq).strip()
            yield name, sequence, qualities


def count_reads(filename):
    count = 0
    for _ in iterate_fastq(filename):
        count += 1
    return count


END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
DIM = '\033[2m'


class MyParser(argparse.ArgumentParser):
    """
    This subclass of ArgumentParser changes the error messages, such that if the script is run with
    no other arguments, it will display the help text. If there is a different error, it will give
    the normal response (usage and error).
    """
    def error(self, message):
        if len(sys.argv) == 1:  # if no arguments were given.
            self.print_help(file=sys.stderr)
            sys.exit(1)
        else:
            super().error(message)


class MyHelpFormatter(argparse.HelpFormatter):

    def __init__(self, prog):
        terminal_width = shutil.get_terminal_size().columns
        os.environ['COLUMNS'] = str(terminal_width)
        max_help_position = min(max(24, terminal_width // 3), 40)
        self.colours = get_colours_from_tput()
        super().__init__(prog, max_help_position=max_help_position)

    def _get_help_string(self, action):
        """
        Override this function to add default values, but only when 'default' is not already in the
        help text.
        """
        help_text = action.help
        if action.default != argparse.SUPPRESS and action.default is not None:
            if 'default' not in help_text.lower():
                help_text += ' (default: {})'.format(action.default)
            elif 'default: DEFAULT' in help_text:
                help_text = help_text.replace('default: DEFAULT',
                                              'default: {}'.format(action.default))
        return help_text

    def start_section(self, heading):
        """
        Override this method to add bold underlining to section headers.
        """
        if self.colours > 1:
            heading = BOLD + heading + END_FORMATTING
        super().start_section(heading)

    def _format_action(self, action):
        """
        Override this method to make help descriptions dim.
        """
        help_position = min(self._action_max_length + 2, self._max_help_position)
        help_width = self._width - help_position
        action_width = help_position - self._current_indent - 2
        action_header = self._format_action_invocation(action)
        if not action.help:
            tup = self._current_indent, '', action_header
            action_header = '%*s%s\n' % tup
            indent_first = 0
        elif len(action_header) <= action_width:
            tup = self._current_indent, '', action_width, action_header
            action_header = '%*s%-*s  ' % tup
            indent_first = 0
        else:
            tup = self._current_indent, '', action_header
            action_header = '%*s%s\n' % tup
            indent_first = help_position
        parts = [action_header]
        if action.help:
            help_text = self._expand_help(action)
            help_lines = self._split_lines(help_text, help_width)
            first_line = help_lines[0]
            if self.colours > 8:
                first_line = DIM + first_line + END_FORMATTING
            parts.append('%*s%s\n' % (indent_first, '', first_line))
            for line in help_lines[1:]:
                if self.colours > 8:
                    line = DIM + line + END_FORMATTING
                parts.append('%*s%s\n' % (help_position, '', line))
        elif not action_header.endswith('\n'):
            parts.append('\n')
        for subaction in self._iter_indented_subactions(action):
            parts.append(self._format_action(subaction))
        return self._join_parts(parts)


def get_colours_from_tput():
    try:
        return int(subprocess.check_output(['tput', 'colors']).decode().strip())
    except (ValueError, subprocess.CalledProcessError, FileNotFoundError, AttributeError):
        return 1


def get_default_thread_count():
    return min(multiprocessing.cpu_count(), 16)


def print_stderr(message, end='\n'):
    print(message, file=sys.stderr, flush=True, end=end)

if __name__ == '__main__':
    main()
