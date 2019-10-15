"""
This module contains a class for writing output to both stdout and a log file.

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

import textwrap
import shutil
import sys


def log(message='', end='\n'):
    print(message, file=sys.stderr, flush=True, end=end)


def section_header(text):
    log()
    print(bold_yellow_underline(text), file=sys.stderr, flush=True)


END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
YELLOW = '\033[93m'
DIM = '\033[2m'


def bold_yellow_underline(text):
    return YELLOW + BOLD + UNDERLINE + text + END_FORMATTING


def dim(text):
    return DIM + text + END_FORMATTING


def explanation(text, indent_size=4):
    """
    This function writes explanatory text to the screen. It is wrapped to the terminal width for
    stdout but not wrapped for the log file.
    """
    text = ' ' * indent_size + text
    # TODO: make terminal width work for stderr
    terminal_width = shutil.get_terminal_size().columns
    for line in textwrap.wrap(text, width=terminal_width - 1):
        log(dim(line))
    log()
