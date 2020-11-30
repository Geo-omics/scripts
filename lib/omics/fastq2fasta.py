# Copyright 2019 Regents of The University of Michigan.

# This file is part of geo-omics-scripts.

# Geo-omics-scripts is free software: you can redistribute it and/or
# modify it under the terms of the GNU General Public License as published
# by the Free Software Foundation, either version 3 of the License, or (at
# your option) any later version.

# Geo-omics-scripts is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.

# You should have received a copy of the GNU General Public License along
# with Geo-omics-scripts.  If not, see <https://www.gnu.org/licenses/>.

"""
Convert fastq into fasta
"""

import argparse
import sys

from . import OmicsArgParser


def convert(data, output, check=True):
    """
    Convert data from FASTQ into FASTA format

    :param data: File-like object with input data
    :param output: File-like object for output
    """
    state = 0
    for line in data:
        if check:
            if state == 0 and not line.startswith('@'):
                raise RuntimeError('Input not in FASTQ format? Expected line '
                                   'to start with @: {}'.format(line))
            elif state == 2 and not line == '+\n':
                raise RuntimeError('Input not in FASTQ format? Expected + '
                                   'separator line: {}'.format(line))

        if state == 0:
            output.write('>' + line[1:])
        elif state == 1:
            output.write(line)

        state = (state + 1) % 4

    if state != 0:
        raise RuntimeError('Input not in FASTQ format? Expected total number '
                           'of lines to be multiple of 4, last line: {}'
                           ''.format(line))


def main(argv=None, namespace=None):
    argp = OmicsArgParser(
        prog=__loader__.name.replace('.', ' '),
        description=__doc__,
        project_home=False,
        threads=False,
    )
    argp.add_argument(
        'inputfile',
        metavar='FILE',
        nargs='?',
        type=argparse.FileType('r'),
        default=sys.stdin,
        help='Fastq file to be converted, by default data is read from stdin.'
    )
    argp.add_argument(
        '-o', '--output',
        metavar='FILE',
        nargs='?',
        type=argparse.FileType('w'),
        default=sys.stdout,
        help='Name of output filie.  Write to stdout by default.'
    )
    argp.add_argument(
        '--force', '-f',
        action='store_true',
        help='Overwrite existing files',
    )
    argp.add_argument(
        '--no-check',
        action='store_false',
        dest='check',
        help='Skip sanity check on input data.  By default it is checked that '
             'the input is indeed in fastq format.',
    )
    args = argp.parse_args(args=argv, namespace=namespace)
    convert(args.inputfile, args.output, check=args.check)


if __name__ == '__main__':
    main()
