#!/usr/bin/env python3
"""
This is the Minipolish installation script. Assuming you're in the same directory, it can be run
like this: `python3 setup.py install`, or (probably better) like this: `pip3 install .`

Copyright 2018 Ryan Wick (rrwick@gmail.com)
https://github.com/rrwick/Minipolish

This file is part of Minipolish. Minipolish is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by the Free Software Foundation,
either version 3 of the License, or (at your option) any later version. Minipolish is distributed
in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
details. You should have received a copy of the GNU General Public License along with Minipolish.
If not, see <http://www.gnu.org/licenses/>.
"""

import os
import shutil
import sys

from setuptools import setup
from setuptools.command.install import install


def readme():
    with open('README.md') as f:
        return f.read()


# Get the program version from another file.
__version__ = '0.0.0'
exec(open('minipolish/version.py').read())


setup(name='Minipolish',
      version=__version__,
      description='Minipolish: a tool for polishing miniasm assemblies with Racon',
      long_description=readme(),
      long_description_content_type='text/markdown',
      url='https://github.com/rrwick/Minipolish',
      author='Ryan Wick',
      author_email='rrwick@gmail.com',
      license='GPLv3',
      packages=['minipolish'],
      install_requires=['edlib'],
      entry_points={"console_scripts": ['minipolish = minipolish.__main__:main']},
      include_package_data=True,
      zip_safe=False,
      python_requires='>=3.6')
